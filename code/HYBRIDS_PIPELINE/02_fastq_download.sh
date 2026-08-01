#!/usr/bin/env bash
#SBATCH --job-name=fastq_dl
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=fastq_download_%j.out
#SBATCH --error=fastq_download_%j.err

# ==============================================================================
# Script: 02_fastq_download.sh
#
# PURPOSE:
# Automates the high-performance, concurrent downloading of FASTQ files using
# `aria2c` based on parsed ENA metadata manifests. Optimized for HPC SLURM.
# ==============================================================================

set -Eeuo pipefail
IFS=$'\n\t'

# ==============================================================================
# Trap bắt lỗi cực kỳ quan trọng cho production
# ==============================================================================
trap 'log "ERROR at line $LINENO"; exit 1' ERR
trap 'log "Received SIGINT/SIGTERM"; exit 130' INT TERM

# ==============================================================================
# 1. Định nghĩa Helper Functions
# ==============================================================================
ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }
die() { log "FATAL: $*"; exit 1; }

# ==============================================================================
# 2. Xử lý Conda thông minh
# ==============================================================================
if [[ "${CONDA_DEFAULT_ENV:-}" != "fastq_dl" ]]; then
    if command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)"
        conda activate "${CONDA_ENV:-fastq_dl}" || die "Cannot activate conda environment"
    else
        log "WARNING: Conda not found in PATH. Proceeding with system tools..."
    fi
fi

# =========================
# Config; override by env vars
# =========================
META_INPUT="${META_INPUT:-/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID/01_metadata}"

DOWNLOAD_ROOT="${DOWNLOAD_ROOT:-/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID/02_fastq}"
FASTQ_ROOT="${FASTQ_ROOT:-$DOWNLOAD_ROOT/fastq}"
MANIFEST_ROOT="${MANIFEST_ROOT:-$DOWNLOAD_ROOT/_manifests}"
LOG_ROOT="${LOG_ROOT:-$DOWNLOAD_ROOT/_logs}"
FAILED_ROOT="${FAILED_ROOT:-$DOWNLOAD_ROOT/_failed}"
ARIA2_LOG_DIR="${ARIA2_LOG_DIR:-$LOG_ROOT/aria2}"

PLAN_TSV="${PLAN_TSV:-$MANIFEST_ROOT/download_plan.tsv}"
STATUS_TSV="${STATUS_TSV:-$MANIFEST_ROOT/download_status.tsv}"
ERROR_TSV="${ERROR_TSV:-$MANIFEST_ROOT/download_errors.tsv}"
READY_TSV="${READY_TSV:-$MANIFEST_ROOT/ready_fastq.tsv}"
PROBLEM_TSV="${PROBLEM_TSV:-$MANIFEST_ROOT/problem_fastq.tsv}"
SUMMARY_TSV="${SUMMARY_TSV:-$MANIFEST_ROOT/download_summary.tsv}"

RUN_ID="${RUN_ID:-$(date '+%Y%m%d_%H%M%S')}"
PARTS_ROOT="${PARTS_ROOT:-$LOG_ROOT/parts/$RUN_ID}"
STATUS_PARTS_DIR="${STATUS_PARTS_DIR:-$PARTS_ROOT/status}"
ERROR_PARTS_DIR="${ERROR_PARTS_DIR:-$PARTS_ROOT/errors}"

# Tự động đồng bộ số luồng với cấu hình SLURM
JOBS="${SLURM_CPUS_PER_TASK:-8}"
SPLIT="${SPLIT:-8}"
CONNECTIONS_PER_FILE="${CONNECTIONS_PER_FILE:-8}"
MIN_SPLIT_SIZE="${MIN_SPLIT_SIZE:-64M}"
MAX_TRIES="${MAX_TRIES:-10}"
RETRY_WAIT="${RETRY_WAIT:-20}"
CONNECT_TIMEOUT="${CONNECT_TIMEOUT:-30}"
TIMEOUT="${TIMEOUT:-600}"
LOWEST_SPEED_LIMIT="${LOWEST_SPEED_LIMIT:-0}"
MAX_DOWNLOAD_LIMIT="${MAX_DOWNLOAD_LIMIT:-}"
URL_SCHEME="${URL_SCHEME:-http}" 

# Disk protection
REQUIRE_SPACE="${REQUIRE_SPACE:-1}"
RESERVE_BYTES="${RESERVE_BYTES:-21474836480}"   # 20 GiB

# Safety/behavior toggles
DRY_RUN=0
FAIL_ON_ERROR=0
FORCE_REDOWNLOAD=0
KEEP_PARTS="${KEEP_PARTS:-0}"
PROJECTS="${PROJECTS:-}"

SCRIPT_SELF="$(readlink -f "$0" 2>/dev/null || printf '%s' "$0")"

# =========================
# CLI
# =========================
usage(){
  cat <<'USAGE'
Usage:
  ./02_fastq_download.sh [options]

Options:
  -n, --dry-run          Test run, build plan without downloading.
  -j, --jobs N           Number of concurrent downloads (default: auto from SLURM).
  --force-redownload     Force redownload of all files.
  --fail-on-error        Exit with error code if any download fails.
USAGE
}

WORKER_MODE=0
WORKER_LINE=""
parse_args(){
  while (( $# )); do
    case "$1" in
      -n|--dry-run) DRY_RUN=1 ;;
      -j|--jobs) JOBS="${2:?Missing value for --jobs}"; shift ;;
      --force-redownload) FORCE_REDOWNLOAD=1 ;;
      --fail-on-error) FAIL_ON_ERROR=1 ;;
      --worker-line) WORKER_MODE=1; WORKER_LINE="${2:?Missing worker line}"; shift ;;
      -h|--help) usage; exit 0 ;;
      *) printf 'ERROR: unknown argument: %s\n\n' "$1" >&2; usage >&2; exit 2 ;;
    esac
    shift
  done
}

# =========================
# TSV helpers & Validation
# =========================
tsv_escape(){
  local s="${1:-}"
  s="${s//$'\t'/ }"
  s="${s//$'\r'/ }"
  s="${s//$'\n'/ }"
  printf '%s' "$s"
}

write_tsv_line(){
  local first=1 v
  for v in "$@"; do
    if (( first )); then first=0; else printf '\t'; fi
    tsv_escape "$v"
  done
  printf '\n'
}

status_header(){
  write_tsv_line ts run_id file_id project_id run_accession file_index status expected_bytes observed_bytes expected_md5 observed_md5 out_path url message
}

summary_header(){
  write_tsv_line ts run_id total ready problems download_root plan_tsv ready_tsv problem_tsv message
}

is_int(){ [[ "${1:-}" =~ ^[0-9]+$ ]]; }
bytes_of(){ [[ -e "$1" ]] && { stat -c '%s' "$1" 2>/dev/null || wc -c <"$1"; } || echo 0; }
md5_of(){ [[ -f "$1" ]] && md5sum "$1" | awk '{print $1}' || echo ""; }
available_bytes(){ df -PB1 "$1" | awk 'NR==2{print $4}'; }
safe_name(){ tr -c 'A-Za-z0-9_.-' '_' <<<"${1:-x}" | sed 's/_$//'; }

require_tools(){
  local missing=()
  for t in python3 awk sort stat df md5sum aria2c; do
    command -v "$t" >/dev/null 2>&1 || missing+=("$t")
  done
  if (( ${#missing[@]} )); then
    die "Missing required tools: ${missing[*]}"
  fi
}

validate_config(){
  [[ -e "$META_INPUT" ]] || die "META_INPUT not found: $META_INPUT"
  [[ "$JOBS" =~ ^[0-9]+$ ]] && (( JOBS >= 1 )) || die "JOBS must be >=1"
  [[ "$SPLIT" =~ ^[0-9]+$ ]] && (( SPLIT >= 1 )) || die "SPLIT must be >=1"
  [[ "$CONNECTIONS_PER_FILE" =~ ^[0-9]+$ ]] && (( CONNECTIONS_PER_FILE >= 1 )) || die "CONNECTIONS_PER_FILE must be >=1"
  [[ "$URL_SCHEME" == "http" || "$URL_SCHEME" == "https" || "$URL_SCHEME" == "ftp" || "$URL_SCHEME" == "preserve" ]] || die "URL_SCHEME must be http, https, ftp, or preserve"
  is_int "$RESERVE_BYTES" || die "RESERVE_BYTES must be an integer"
}

init_output_tree(){
  mkdir -p "$FASTQ_ROOT" "$MANIFEST_ROOT" "$ARIA2_LOG_DIR" "$FAILED_ROOT" "$STATUS_PARTS_DIR" "$ERROR_PARTS_DIR"
  if [[ ! -f "$STATUS_TSV" ]]; then status_header > "$STATUS_TSV"; fi
  if [[ ! -f "$ERROR_TSV" ]]; then status_header > "$ERROR_TSV"; fi
  if [[ ! -f "$SUMMARY_TSV" ]]; then summary_header > "$SUMMARY_TSV"; fi
}

# =========================
# Build the download plan
# =========================
build_plan(){
  local out_plan="$1"
  python3 - "$META_INPUT" "$FASTQ_ROOT" "$PROJECTS" "$URL_SCHEME" "$out_plan" <<'PY'
import csv, os, re, sys
from pathlib import Path
from urllib.parse import urlparse

meta_input, fastq_root, projects_raw, url_scheme, out_plan = sys.argv[1:6]
meta_input = Path(meta_input)
fastq_root = Path(fastq_root)
project_filter = {x.strip() for x in projects_raw.split(',') if x.strip()}

cols = [
    'file_id', 'project_id', 'study_accession', 'sample_accession', 'experiment_accession',
    'run_accession', 'library_layout', 'file_index', 'fastq_url', 'expected_md5',
    'expected_bytes', 'filename', 'out_dir', 'out_path'
]

def safe_component(x, fallback='unknown'):
    x = (x or '').strip()
    if not x:
        x = fallback
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', x)

def convert_url(url):
    url = (url or '').strip()
    if not url:
        return ''
    if '://' not in url:
        url = 'ftp://' + url
    if url_scheme == 'preserve':
        return url
    parsed = urlparse(url)
    if parsed.scheme in ('ftp', 'http', 'https'):
        return url_scheme + '://' + parsed.netloc + parsed.path
    return url

def filename_from_url(url, run, file_index):
    try:
        base = Path(urlparse(url).path).name
    except Exception:
        base = ''
    if not base:
        base = f'{run}_{file_index}.fastq.gz'
    return safe_component(base, f'{run}_{file_index}.fastq.gz')

manifests = []
if meta_input.is_file():
    manifests.append(meta_input)
elif meta_input.is_dir():
    for p in sorted(meta_input.rglob('*.fastq_files.tsv')):
        if p.name == 'latest.fastq_files.tsv':
            continue
        if not p.is_file():
            continue
        manifests.append(p)

print(f"[DEBUG] Found {len(manifests)} manifest(s) in {meta_input}", file=sys.stderr)

seen = set()
with open(out_plan, 'w', newline='') as oh:
    writer = csv.DictWriter(oh, fieldnames=cols, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for mf in manifests:
        with open(mf, newline='') as fh:
            # Auto-detect delimiter (Tab vs Comma)
            first_line = fh.readline()
            fh.seek(0)
            delim = '\t' if '\t' in first_line else ','
            
            reader = csv.DictReader(fh, delimiter=delim)
            headers = reader.fieldnames or []
            print(f"[DEBUG] -> Parsing {mf.name} | Delimiter: '{repr(delim)}' | Columns: {len(headers)}", file=sys.stderr)
            
            row_count = 0
            skipped_url = 0
            for row in reader:
                row_count += 1
                # Chuẩn hóa headers thành chữ thường, bỏ khoảng trắng dư thừa
                row_lower = {str(k).strip().lower(): str(v).strip() for k, v in row.items() if k}
                
                project = safe_component(row.get('project_id') or mf.parent.name, mf.parent.name)
                if project_filter and project not in project_filter:
                    continue
                
                run = safe_component(row_lower.get('run_accession') or row_lower.get('run'), 'unknown_run')
                file_index = safe_component(str(row_lower.get('file_index') or '1'), '1')
                
                # Quét mọi dạng cột có khả năng chứa URL
                raw_url = (row_lower.get('fastq_url') or 
                           row_lower.get('submitted_ftp') or 
                           row_lower.get('fastq_ftp') or 
                           row_lower.get('sra_ftp') or 
                           row_lower.get('url'))
                           
                url = convert_url(raw_url)
                if not url:
                    skipped_url += 1
                    continue
                
                fname = filename_from_url(url, run, file_index)
                out_dir = fastq_root / project / run
                out_path = out_dir / fname
                file_id = f'{project}____{run}____{file_index}'
                
                dedup_key = str(out_path)
                if dedup_key in seen:
                    continue
                seen.add(dedup_key)
                
                raw_md5 = (row_lower.get('fastq_md5') or row_lower.get('submitted_md5') or '')
                raw_bytes = (row_lower.get('fastq_bytes') or row_lower.get('submitted_bytes') or '')
                
                writer.writerow({
                    'file_id': file_id,
                    'project_id': project,
                    'study_accession': row_lower.get('study_accession', ''),
                    'sample_accession': row_lower.get('sample_accession', ''),
                    'experiment_accession': row_lower.get('experiment_accession', ''),
                    'run_accession': run,
                    'library_layout': row_lower.get('library_layout', ''),
                    'file_index': file_index,
                    'fastq_url': url,
                    'expected_md5': raw_md5,
                    'expected_bytes': raw_bytes,
                    'filename': fname,
                    'out_dir': str(out_dir),
                    'out_path': str(out_path),
                })
            
            print(f"[DEBUG] -> Data lines: {row_count} | Skipped (missing URL): {skipped_url}", file=sys.stderr)
PY
}

summarize_plan(){
  local plan="$1"
  python3 - "$plan" <<'PY'
import csv, sys
from collections import Counter
path = sys.argv[1]
rows = list(csv.DictReader(open(path, newline=''), delimiter='\t'))
projects = Counter(r['project_id'] for r in rows)
bytes_total = sum(int(r['expected_bytes']) for r in rows if (r.get('expected_bytes') or '').isdigit())

def human(n):
    units = ['B','KiB','MiB','GiB','TiB','PiB']
    n = float(n)
    for u in units:
        if n < 1024 or u == units[-1]:
            return f'{n:.2f} {u}'
        n /= 1024
print(f'FASTQ files planned: {len(rows)}')
print(f'Projects: {len(projects)}')
print(f'Expected total size from metadata: {human(bytes_total)}')
print('Files per project:')
for p, c in projects.most_common():
    print(f'  {p}: {c}')
print('\nFirst planned files:')
for r in rows[:10]:
    print(f"  {r['project_id']} / {r['run_accession']} / {r['filename']} -> {r['out_path']}")
PY
}

# =========================
# Per-file worker
# =========================
write_status_part(){
  local file_id="$1" project_id="$2" run_accession="$3" file_index="$4" status="$5" expected_bytes="$6" observed_bytes="$7" expected_md5="$8" observed_md5="$9" out_path="${10}" url="${11}" message="${12}"
  local safe_id part
  safe_id="$(safe_name "$file_id")"
  part="$STATUS_PARTS_DIR/${safe_id}.tsv"

  write_tsv_line "$(ts)" "$RUN_ID" "$file_id" "$project_id" "$run_accession" "$file_index" "$status" "$expected_bytes" "$observed_bytes" "$expected_md5" "$observed_md5" "$out_path" "$url" "$message" > "$part"

  case "$status" in
    success*|ok_existing*) ;;
    *) write_tsv_line "$(ts)" "$RUN_ID" "$file_id" "$project_id" "$run_accession" "$file_index" "$status" "$expected_bytes" "$observed_bytes" "$expected_md5" "$observed_md5" "$out_path" "$url" "$message" > "$ERROR_PARTS_DIR/${safe_id}.tsv" ;;
  esac
}

quarantine_file(){
  local src="$1" reason="$2" project_id="$3" run_accession="$4" filename="$5"
  local dst_dir dst stamp
  stamp="$(date '+%Y%m%d_%H%M%S')"
  dst_dir="$FAILED_ROOT/$reason/$project_id/$run_accession"
  mkdir -p "$dst_dir"
  dst="$dst_dir/${filename}.${stamp}.${reason}"
  mv -f "$src" "$dst"
  [[ -f "$src.aria2" ]] && mv -f "$src.aria2" "$dst.aria2" || true
  printf '%s' "$dst"
}

download_one_line(){
  local line="$1"
  local file_id project_id study_accession sample_accession experiment_accession run_accession library_layout file_index url expected_md5 expected_bytes filename out_dir out_path
  IFS=$'\t' read -r file_id project_id study_accession sample_accession experiment_accession run_accession library_layout file_index url expected_md5 expected_bytes filename out_dir out_path <<< "$line"

  local observed_bytes observed_md5 aria2_log aria2_project_dir
  observed_bytes="$(bytes_of "$out_path")"
  observed_md5=""

  mkdir -p "$out_dir"
  aria2_project_dir="$ARIA2_LOG_DIR/$project_id/$run_accession"
  mkdir -p "$aria2_project_dir"
  aria2_log="$aria2_project_dir/${filename}.aria2.log"

  if [[ -f "$out_path" && "$FORCE_REDOWNLOAD" != "1" ]]; then
    if [[ -n "$expected_md5" && ! -f "$out_path.aria2" ]]; then
      observed_md5="$(md5_of "$out_path")"
      if [[ "$observed_md5" == "$expected_md5" ]]; then
        write_status_part "$file_id" "$project_id" "$run_accession" "$file_index" "ok_existing_md5" "$expected_bytes" "$observed_bytes" "$expected_md5" "$observed_md5" "$out_path" "$url" "Existing file already matches expected MD5"
        return 0
      else
        local qpath
        qpath="$(quarantine_file "$out_path" "bad_md5" "$project_id" "$run_accession" "$filename")"
        log "Pre-existing MD5 mismatch: moved $out_path to $qpath; attempting a fresh download now"
        observed_bytes=0
        observed_md5=""
      fi
    elif is_int "$expected_bytes" && [[ "$observed_bytes" == "$expected_bytes" && ! -f "$out_path.aria2" ]]; then
      write_status_part "$file_id" "$project_id" "$run_accession" "$file_index" "ok_existing_size_only" "$expected_bytes" "$observed_bytes" "$expected_md5" "" "$out_path" "$url" "Existing file size matches expected bytes; no MD5 available"
      return 0
    elif [[ ! -f "$out_path.aria2" ]]; then
      write_status_part "$file_id" "$project_id" "$run_accession" "$file_index" "exists_unverified_skipped" "$expected_bytes" "$observed_bytes" "$expected_md5" "" "$out_path" "$url" "Existing file cannot be verified safely; use --force-redownload if you want to quarantine and redownload"
      return 0
    fi
  fi

  if [[ -f "$out_path" && "$FORCE_REDOWNLOAD" == "1" && ! -f "$out_path.aria2" ]]; then
    local qpath
    qpath="$(quarantine_file "$out_path" "force_redownload" "$project_id" "$run_accession" "$filename")"
    observed_bytes=0
    observed_md5=""
    log "Force-redownload: moved $out_path to $qpath"
  fi

  if [[ "$REQUIRE_SPACE" == "1" && -n "$expected_bytes" ]] && is_int "$expected_bytes"; then
    observed_bytes="$(bytes_of "$out_path")"
    local missing_bytes available need_total
    if (( expected_bytes > observed_bytes )); then missing_bytes=$(( expected_bytes - observed_bytes )); else missing_bytes=0; fi
    available="$(available_bytes "$out_dir")"
    need_total=$(( missing_bytes + RESERVE_BYTES ))
    if (( available < need_total )); then
      write_status_part "$file_id" "$project_id" "$run_accession" "$file_index" "skipped_no_space" "$expected_bytes" "$observed_bytes" "$expected_md5" "" "$out_path" "$url" "Available bytes=$available; need missing bytes + reserve=$need_total. Free space or change DOWNLOAD_ROOT/RESERVE_BYTES, then rerun."
      return 0
    fi
  fi

  local cmd rc
  cmd=(
    aria2c
    --continue=true
    --auto-file-renaming=false
    --allow-overwrite=false
    --file-allocation=none
    --remote-time=true
    --split="$SPLIT"
    --max-connection-per-server="$CONNECTIONS_PER_FILE"
    --min-split-size="$MIN_SPLIT_SIZE"
    --max-tries="$MAX_TRIES"
    --retry-wait="$RETRY_WAIT"
    --connect-timeout="$CONNECT_TIMEOUT"
    --timeout="$TIMEOUT"
    --summary-interval=60
    --console-log-level=warn
    --check-certificate=false
    --dir="$out_dir"
    --out="$filename"
  )
  if [[ -n "$MAX_DOWNLOAD_LIMIT" ]]; then cmd+=(--max-download-limit="$MAX_DOWNLOAD_LIMIT"); fi
  if [[ "$LOWEST_SPEED_LIMIT" != "0" ]]; then cmd+=(--lowest-speed-limit="$LOWEST_SPEED_LIMIT"); fi
  cmd+=("$url")

  printf '[%s] START %s\n' "$(ts)" "${cmd[*]}" >> "$aria2_log"
  if "${cmd[@]}" >> "$aria2_log" 2>&1; then
    rc=0
  else
    rc=$?
  fi
  printf '[%s] END rc=%s\n' "$(ts)" "$rc" >> "$aria2_log"

  observed_bytes="$(bytes_of "$out_path")"
  if (( rc != 0 )); then
    write_status_part "$file_id" "$project_id" "$run_accession" "$file_index" "download_failed" "$expected_bytes" "$observed_bytes" "$expected_md5" "" "$out_path" "$url" "aria2c exit code=$rc; partial file is kept for resume. See $aria2_log"
    return 0
  fi

  if [[ ! -f "$out_path" ]]; then
    write_status_part "$file_id" "$project_id" "$run_accession" "$file_index" "download_missing_output" "$expected_bytes" "0" "$expected_md5" "" "$out_path" "$url" "aria2c returned success but output file was not found. See $aria2_log"
    return 0
  fi

  if [[ -n "$expected_md5" ]]; then
    observed_md5="$(md5_of "$out_path")"
    if [[ "$observed_md5" == "$expected_md5" ]]; then
      write_status_part "$file_id" "$project_id" "$run_accession" "$file_index" "success_md5" "$expected_bytes" "$observed_bytes" "$expected_md5" "$observed_md5" "$out_path" "$url" "Downloaded and MD5 verified"
    else
      local qpath
      qpath="$(quarantine_file "$out_path" "bad_md5" "$project_id" "$run_accession" "$filename")"
      write_status_part "$file_id" "$project_id" "$run_accession" "$file_index" "md5_mismatch_quarantined" "$expected_bytes" "$observed_bytes" "$expected_md5" "$observed_md5" "$out_path" "$url" "Downloaded file failed MD5 and was moved to $qpath. See $aria2_log"
    fi
  elif is_int "$expected_bytes"; then
    if [[ "$observed_bytes" == "$expected_bytes" ]]; then
      write_status_part "$file_id" "$project_id" "$run_accession" "$file_index" "success_size_only" "$expected_bytes" "$observed_bytes" "$expected_md5" "" "$out_path" "$url" "Downloaded and file size matches expected bytes; no MD5 available"
    else
      local qpath
      qpath="$(quarantine_file "$out_path" "bad_size" "$project_id" "$run_accession" "$filename")"
      write_status_part "$file_id" "$project_id" "$run_accession" "$file_index" "size_mismatch_quarantined" "$expected_bytes" "$observed_bytes" "$expected_md5" "" "$out_path" "$url" "Downloaded file failed size check and was moved to $qpath. See $aria2_log"
    fi
  else
    write_status_part "$file_id" "$project_id" "$run_accession" "$file_index" "success_unverified" "$expected_bytes" "$observed_bytes" "$expected_md5" "" "$out_path" "$url" "Downloaded but no MD5/size was available in metadata"
  fi
}

# =========================
# Merge current-run status
# =========================
merge_status_parts(){
  local current_status="$1" current_errors="$2"
  status_header > "$current_status"
  if compgen -G "$STATUS_PARTS_DIR/*.tsv" >/dev/null; then
    cat "$STATUS_PARTS_DIR"/*.tsv >> "$current_status"
  fi
  status_header > "$current_errors"
  if compgen -G "$ERROR_PARTS_DIR/*.tsv" >/dev/null; then
    cat "$ERROR_PARTS_DIR"/*.tsv >> "$current_errors"
  fi
  tail -n +2 "$current_status" >> "$STATUS_TSV"
  tail -n +2 "$current_errors" >> "$ERROR_TSV"
}

build_ready_problem_tables(){
  local plan="$1" current_status="$2" ready_out="$3" problem_out="$4"
  python3 - "$plan" "$current_status" "$ready_out" "$problem_out" <<'PY'
import csv, sys
plan_path, status_path, ready_out, problem_out = sys.argv[1:5]
ready_status = {'success_md5', 'success_size_only', 'success_unverified', 'ok_existing_md5', 'ok_existing_size_only'}
plan = {}
with open(plan_path, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        plan[r['file_id']] = r
status = {}
with open(status_path, newline='') as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        status[r['file_id']] = r
cols = [
    'file_id','project_id','study_accession','sample_accession','experiment_accession','run_accession',
    'library_layout','file_index','fastq_url','expected_md5','expected_bytes','out_path',
    'status','observed_md5','observed_bytes','message'
]
with open(ready_out, 'w', newline='') as ro, open(problem_out, 'w', newline='') as po:
    rw = csv.DictWriter(ro, fieldnames=cols, delimiter='\t', lineterminator='\n')
    pw = csv.DictWriter(po, fieldnames=cols, delimiter='\t', lineterminator='\n')
    rw.writeheader(); pw.writeheader()
    for fid, p in sorted(plan.items(), key=lambda kv: (kv[1]['project_id'], kv[1]['run_accession'], kv[1]['file_index'])):
        s = status.get(fid, {})
        row = {
            'file_id': fid,
            'project_id': p.get('project_id',''),
            'study_accession': p.get('study_accession',''),
            'sample_accession': p.get('sample_accession',''),
            'experiment_accession': p.get('experiment_accession',''),
            'run_accession': p.get('run_accession',''),
            'library_layout': p.get('library_layout',''),
            'file_index': p.get('file_index',''),
            'fastq_url': p.get('fastq_url',''),
            'expected_md5': p.get('expected_md5',''),
            'expected_bytes': p.get('expected_bytes',''),
            'out_path': p.get('out_path',''),
            'status': s.get('status','not_processed'),
            'observed_md5': s.get('observed_md5',''),
            'observed_bytes': s.get('observed_bytes',''),
            'message': s.get('message','No status record produced in current run'),
        }
        (rw if row['status'] in ready_status else pw).writerow(row)
PY
}

count_rows(){ awk 'NR>1{n++} END{print n+0}' "$1"; }
summarize_status_counts(){ awk -F'\t' 'NR>1{c[$7]++} END{for (s in c) print s, c[s]}' "$1" | sort; }

# =========================
# Main
# =========================
main(){
  parse_args "$@"
  require_tools
  validate_config

  if (( WORKER_MODE )); then
    download_one_line "$WORKER_LINE"
    exit 0
  fi

  if (( DRY_RUN )); then
    local tmp_plan="$(mktemp)"
    log "DRY-RUN: building plan only. No output folders/files will be created under DOWNLOAD_ROOT."
    build_plan "$tmp_plan"
    summarize_plan "$tmp_plan"
    rm -f "$tmp_plan"
    exit 0
  fi

  init_output_tree

  # Log SLURM Info
  log "Hostname: $(hostname)"
  log "Job ID: ${SLURM_JOB_ID:-N/A}"
  log "CPUs: ${SLURM_CPUS_PER_TASK:-N/A}"
  log "Memory: ${SLURM_MEM_PER_NODE:-N/A}"

  # Disk Check
  log "Disk usage before run:"
  df -h "$DOWNLOAD_ROOT" || true

  log "Building download plan from: $META_INPUT"
  build_plan "$PLAN_TSV"

  local total_files
  total_files="$(count_rows "$PLAN_TSV")"
  if (( total_files == 0 )); then
    log "No FASTQ files found in metadata manifests. Check META_INPUT=$META_INPUT"
    exit 0
  fi

  # In ra expected size
  log "--- Plan Summary ---"
  while IFS= read -r line; do
      log "$line"
  done < <(summarize_plan "$PLAN_TSV")
  log "--------------------"

  log "Download plan: $PLAN_TSV"
  log "Planned FASTQ files: $total_files"
  log "Output FASTQ root: $FASTQ_ROOT"
  log "Concurrency: JOBS=$JOBS, SPLIT=$SPLIT, CONNECTIONS_PER_FILE=$CONNECTIONS_PER_FILE"
  log "Run ID: $RUN_ID"

  export META_INPUT DOWNLOAD_ROOT FASTQ_ROOT MANIFEST_ROOT LOG_ROOT FAILED_ROOT ARIA2_LOG_DIR
  export RUN_ID PARTS_ROOT STATUS_PARTS_DIR ERROR_PARTS_DIR
  export JOBS SPLIT CONNECTIONS_PER_FILE MIN_SPLIT_SIZE MAX_TRIES RETRY_WAIT CONNECT_TIMEOUT TIMEOUT LOWEST_SPEED_LIMIT MAX_DOWNLOAD_LIMIT
  export REQUIRE_SPACE RESERVE_BYTES FORCE_REDOWNLOAD

  tail -n +2 "$PLAN_TSV" | xargs -r -P "$JOBS" -d '\n' -I{} "$SCRIPT_SELF" --worker-line "{}"

  local current_status current_errors
  current_status="$MANIFEST_ROOT/download_status.${RUN_ID}.tsv"
  current_errors="$MANIFEST_ROOT/download_errors.${RUN_ID}.tsv"
  merge_status_parts "$current_status" "$current_errors"
  build_ready_problem_tables "$PLAN_TSV" "$current_status" "$READY_TSV" "$PROBLEM_TSV"

  local ready_count problem_count
  ready_count="$(count_rows "$READY_TSV")"
  problem_count="$(count_rows "$PROBLEM_TSV")"

  write_tsv_line "$(ts)" "$RUN_ID" "$total_files" "$ready_count" "$problem_count" "$DOWNLOAD_ROOT" "$PLAN_TSV" "$READY_TSV" "$PROBLEM_TSV" "Completed current run" >> "$SUMMARY_TSV"

  log "Current run status counts:"
  summarize_status_counts "$current_status" >&2 || true
  log "Ready FASTQ table: $READY_TSV"
  log "Problem FASTQ table: $PROBLEM_TSV"
  log "Cumulative status log: $STATUS_TSV"
  log "Cumulative error log: $ERROR_TSV"
  log "aria2 logs: $ARIA2_LOG_DIR"
  log "Quarantined files, if any: $FAILED_ROOT"

  if [[ "$KEEP_PARTS" != "1" ]]; then
    rm -rf "$PARTS_ROOT"
  else
    log "Kept worker parts: $PARTS_ROOT"
  fi

  if (( FAIL_ON_ERROR && problem_count > 0 )); then
    exit 1
  fi
}

main "$@"