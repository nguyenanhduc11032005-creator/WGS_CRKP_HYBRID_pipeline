#!/usr/bin/env bash
# ==============================================================================
# Script: 01_fetch_metadata.sh
#
# PURPOSE:
# Automates the bulk retrieval of ENA (European Nucleotide Archive) `read_run` 
# metadata for a list of study/project accessions provided in a CSV file. 
# It handles API requests, validates fields, and transforms the raw ENA 
# responses into flattened, downstream-friendly FASTQ manifests to streamline 
# large-scale genomic download pipelines.
#
# EXPECTED OUTPUTS:
# For each project, a dedicated folder is created under $META_ROOT containing:
#   1. <project_id>.tsv              : The raw, complete metadata TSV from ENA.
#   2. <project_id>.fastq_files.tsv  : A parsed manifest expanding semicolon-separated 
#                                      URLs, MD5s, and sizes into distinct rows.
#   3. _run_log.tsv                  : Granular logging of the download/validation steps.
#
# At the root of $META_ROOT, two global summaries are generated:
#   1. _manifest_projects.tsv        : High-level overview of success/failure per project.
#   2. _manifest_runs.tsv            : A concatenated master TSV of all successfully 
#                                      retrieved run metadata.
#
# TYPICAL USAGE:
# Run with defaults (defined in the script config block):
#   ./01_fetch_metadata.sh
#
# Run with custom input/output paths (via environment variables):
#   CSV_FILE="/path/to/accessions.csv" META_ROOT="/path/to/outdir" ./01_fetch_metadata.sh
#
# ARGUMENTS & FLAGS:
#   -n, --dry-run   : Validates the input CSV and prints planned ENA API requests 
#                     and output paths without executing any network calls or 
#                     writing any files.
#   -h, --help      : Displays the help message.
#
# KEY ENVIRONMENT VARIABLES (Overridable):
#   CSV_FILE        : Input CSV containing accessions.
#   META_ROOT       : Output directory for all metadata and manifests.
#   MAX_RETRIES     : ENA download attempt limit (default: 5).
#   MAX_TIME        : Maximum time allowed for curl to complete (default: 300s).
# ==============================================================================

set -Eeuo pipefail
IFS=$'\n\t'

# =========================
# Config; override by env vars
# =========================
CSV_FILE="${CSV_FILE:-/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID/accessions.csv}"
META_ROOT="${META_ROOT:-/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID/01_metadata}"
GLOBAL_MANIFEST="${GLOBAL_MANIFEST:-$META_ROOT/_manifest_projects.tsv}"
GLOBAL_RUNS="${GLOBAL_RUNS:-$META_ROOT/_manifest_runs.tsv}"
USER_AGENT="${USER_AGENT:-ena-meta-fetch/1.1}"
MAX_RETRIES="${MAX_RETRIES:-5}"
RETRY_DELAY="${RETRY_DELAY:-3}"
CONNECT_TIMEOUT="${CONNECT_TIMEOUT:-20}"
MAX_TIME="${MAX_TIME:-300}"
DRY_RUN=0

# Keep this list conservative. Add more ENA fields only after verifying the field exists.
REQUIRED_FIELDS=(
  study_accession
  sample_accession
  experiment_accession
  run_accession
  scientific_name
  tax_id
  fastq_md5
  fastq_ftp
  fastq_bytes
  sra_ftp
  library_name
  library_strategy
  library_source
  library_selection
  library_layout
  instrument_platform
  instrument_model
  run_alias
  read_count
  base_count
)

# =========================
# CLI helpers
# =========================
usage(){
  cat <<'EOF'
Usage:
  CSV_FILE=/path/to/accessions.csv META_ROOT=/path/to/metadata ./01_fetch_metadata.sh [options]

Options:
  -n, --dry-run   Read and validate the CSV, then print planned ENA requests and output paths.
                  No network requests are sent and no files/directories are created.
  -h, --help      Show this help message.

Environment variables:
  CSV_FILE          Input CSV/TSV/semicolon file containing study/project accessions.
  META_ROOT         Output metadata root directory.
  GLOBAL_MANIFEST   Optional project-level manifest path.
  GLOBAL_RUNS       Optional combined run-level manifest path.
  MAX_RETRIES       Number of explicit download attempts.
  RETRY_DELAY       Delay between attempts, in seconds.
EOF
}

parse_args(){
  while (( $# )); do
    case "$1" in
      -n|--dry-run) DRY_RUN=1 ;;
      -h|--help) usage; exit 0 ;;
      *)
        printf 'ERROR: unknown argument: %s\n\n' "$1" >&2
        usage >&2
        exit 2
        ;;
    esac
    shift
  done
}

# =========================
# Logging helpers
# =========================
ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }

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

init_global_manifest(){
  mkdir -p "$META_ROOT"
  if [[ ! -f "$GLOBAL_MANIFEST" ]]; then
    write_tsv_line ts project_id status bytes rows md5 tsv_path message > "$GLOBAL_MANIFEST"
  fi
}

append_global_manifest(){
  write_tsv_line "$@" >> "$GLOBAL_MANIFEST"
}

init_run_log(){
  local proj_dir="$1"
  local logf="$proj_dir/_run_log.tsv"
  if [[ ! -f "$logf" ]]; then
    write_tsv_line ts stage status field value message > "$logf"
  fi
}

append_run_log(){
  local proj_dir="$1" stage="$2" status="$3" field="$4" value="$5" msg="$6"
  write_tsv_line "$(ts)" "$stage" "$status" "$field" "$value" "$msg" >> "$proj_dir/_run_log.tsv"
}

# =========================
# Utilities
# =========================
bytes_of(){ [[ -f "$1" ]] && { stat -c '%s' "$1" 2>/dev/null || wc -c <"$1"; } || echo 0; }
md5_of(){ [[ -f "$1" ]] && md5sum "$1" | awk '{print $1}' || echo ""; }
fields_csv(){ local IFS=,; printf '%s' "${REQUIRED_FIELDS[*]}"; }

ena_request_url(){
  local project_id="$1" fields="$2"
  python3 - "$project_id" "$fields" <<'PY'
import sys
from urllib.parse import urlencode
project_id, fields = sys.argv[1], sys.argv[2]
params = {
    'accession': project_id,
    'result': 'read_run',
    'fields': fields,
    'format': 'tsv',
    'download': 'true',
    'limit': '0',
}
print('https://www.ebi.ac.uk/ena/portal/api/filereport?' + urlencode(params))
PY
}

validate_project_id(){
  local id="$1"
  # Prevent path traversal and shell/path surprises before using the accession as a directory name.
  [[ "$id" =~ ^[A-Za-z0-9_.-]+$ ]] || return 1
  # Broad ENA/SRA accession prefixes accepted by read_run reports.
  [[ "$id" =~ ^(PRJ|SRP|ERP|DRP|SRX|ERX|DRX|SRR|ERR|DRR|SAM|SRS|ERS|DRS) ]] || return 1
}

tsv_missing_required_fields(){
  local tsv="$1"
  python3 - "$tsv" "$(fields_csv)" <<'PY'
import csv, sys
path, fields_csv = sys.argv[1], sys.argv[2]
required = fields_csv.split(',')
with open(path, newline='') as fh:
    reader = csv.reader(fh, delimiter='\t')
    try:
        header = next(reader)
    except StopIteration:
        print(','.join(required))
        sys.exit(1)
header = [h.strip().lstrip('\ufeff') for h in header]
missing = [f for f in required if f not in header]
print(','.join(missing))
sys.exit(1 if missing else 0)
PY
}

count_data_rows(){
  awk 'END{if(NR>0) print NR-1; else print 0}' "$1"
}

# =========================
# CSV reader: robust to quoted commas, BOM, CRLF
# =========================
read_projects_from_csv(){
  local file="$1"
  [[ -s "$file" ]] || { log "ERROR: CSV not found or empty: $file"; exit 1; }
  python3 - "$file" <<'PY'
import csv, sys
from pathlib import Path
path = Path(sys.argv[1])
candidates = ['study_accession', 'project_accession', 'project_id', 'study_id']
with path.open(newline='', encoding='utf-8-sig') as fh:
    sample = fh.read(4096)
    fh.seek(0)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=',\t;')
    except csv.Error:
        dialect = csv.excel
    reader = csv.DictReader(fh, dialect=dialect)
    if not reader.fieldnames:
        sys.exit('CSV has no header')
    norm = {name.strip().lower(): name for name in reader.fieldnames if name is not None}
    selected = next((norm[c] for c in candidates if c in norm), None)
    if selected is None:
        # fallback to second column, matching the original script behavior
        selected = reader.fieldnames[1] if len(reader.fieldnames) >= 2 else reader.fieldnames[0]
    seen = set()
    for row in reader:
        value = (row.get(selected) or '').strip().strip('"').strip("'")
        if value and value not in seen:
            seen.add(value)
            print(value)
PY
}

# =========================
# Downstream-friendly file manifest
# =========================
build_fastq_file_manifest(){
  local project_id="$1" tsv="$2" out="$3"
  python3 - "$project_id" "$tsv" "$out" <<'PY'
import csv, sys
project_id, tsv, out = sys.argv[1:4]
cols = [
    'project_id', 'study_accession', 'sample_accession', 'experiment_accession',
    'run_accession', 'library_layout', 'file_index', 'fastq_url', 'fastq_md5', 'fastq_bytes'
]

def split_field(x):
    if x is None or x == '':
        return []
    return [v.strip() for v in x.split(';') if v.strip()]

def as_url(x):
    if not x:
        return ''
    if '://' in x:
        return x
    return 'ftp://' + x

with open(tsv, newline='') as fh, open(out, 'w', newline='') as oh:
    reader = csv.DictReader(fh, delimiter='\t')
    writer = csv.DictWriter(oh, fieldnames=cols, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for row in reader:
        urls = split_field(row.get('fastq_ftp'))
        md5s = split_field(row.get('fastq_md5'))
        sizes = split_field(row.get('fastq_bytes'))
        max_n = max(len(urls), len(md5s), len(sizes), 0)
        for i in range(max_n):
            writer.writerow({
                'project_id': project_id,
                'study_accession': row.get('study_accession', ''),
                'sample_accession': row.get('sample_accession', ''),
                'experiment_accession': row.get('experiment_accession', ''),
                'run_accession': row.get('run_accession', ''),
                'library_layout': row.get('library_layout', ''),
                'file_index': i + 1,
                'fastq_url': as_url(urls[i]) if i < len(urls) else '',
                'fastq_md5': md5s[i] if i < len(md5s) else '',
                'fastq_bytes': sizes[i] if i < len(sizes) else '',
            })
PY
}

# =========================
# ENA downloader
# =========================
fetch_project_tsv(){
  local project_id="$1"
  if ! validate_project_id "$project_id"; then
    log "Project $project_id skipped: invalid or suspicious accession format"
    if (( ! DRY_RUN )); then
      append_global_manifest "$(ts)" "$project_id" "skipped" "0" "0" "" "" "Invalid accession format"
    fi
    return 1
  fi

  local proj_dir="$META_ROOT/$project_id"
  local final_tsv="$proj_dir/${project_id}.tsv"
  local tmp_tsv="$proj_dir/.${project_id}.tsv.$$.part"
  local fastq_manifest="$proj_dir/${project_id}.fastq_files.tsv"
  local fields; fields="$(fields_csv)"
  local request_url; request_url="$(ena_request_url "$project_id" "$fields")"

  if (( DRY_RUN )); then
    log "[DRY-RUN] project_id: $project_id"
    log "[DRY-RUN] would create directory: $proj_dir"
    log "[DRY-RUN] would request: $request_url"
    log "[DRY-RUN] would write raw run metadata: $final_tsv"
    log "[DRY-RUN] would write FASTQ file manifest: $fastq_manifest"
    log "[DRY-RUN] would update global manifests: $GLOBAL_MANIFEST ; $GLOBAL_RUNS"
    return 0
  fi

  mkdir -p "$proj_dir"
  init_run_log "$proj_dir"
  rm -f "$tmp_tsv"
  append_run_log "$proj_dir" "prepare" "info" "fields" "$fields" "Built ENA read_run request"

  local attempt=1
  while (( attempt <= MAX_RETRIES )); do
    append_run_log "$proj_dir" "download" "start" "attempt" "$attempt" "Fetching TSV"
    rm -f "$tmp_tsv"
    if curl --get --fail-with-body --location --silent --show-error \
        --user-agent "$USER_AGENT" \
        --connect-timeout "$CONNECT_TIMEOUT" --max-time "$MAX_TIME" \
        --retry 0 \
        --data-urlencode "accession=$project_id" \
        --data-urlencode "result=read_run" \
        --data-urlencode "fields=$fields" \
        --data "format=tsv" \
        --data "download=true" \
        --data "limit=0" \
        "https://www.ebi.ac.uk/ena/portal/api/filereport" \
        -o "$tmp_tsv"; then
      append_run_log "$proj_dir" "download" "ok" "attempt" "$attempt" "Downloaded"
      break
    else
      append_run_log "$proj_dir" "download" "error" "attempt" "$attempt" "curl failed"
      sleep "$RETRY_DELAY"
    fi
    (( attempt++ ))
  done

  if [[ ! -s "$tmp_tsv" ]]; then
    append_run_log "$proj_dir" "download" "failed" "file" "$tmp_tsv" "Empty or missing after retries"
    append_global_manifest "$(ts)" "$project_id" "failed" "0" "0" "" "" "Download failed or empty response"
    rm -f "$tmp_tsv"
    return 1
  fi

  local missing
  if ! missing="$(tsv_missing_required_fields "$tmp_tsv")"; then
    append_run_log "$proj_dir" "validate" "error" "missing_fields" "$missing" "Header missing required fields"
    append_global_manifest "$(ts)" "$project_id" "failed" "$(bytes_of "$tmp_tsv")" "0" "" "$tmp_tsv" "Missing fields: $missing"
    return 1
  fi

  local rows; rows="$(count_data_rows "$tmp_tsv")"
  if (( rows == 0 )); then
    append_run_log "$proj_dir" "validate" "warning" "rows" "$rows" "No run rows returned"
  else
    append_run_log "$proj_dir" "validate" "ok" "rows" "$rows" "TSV header and row count look valid"
  fi

  mv -f "$tmp_tsv" "$final_tsv"
  build_fastq_file_manifest "$project_id" "$final_tsv" "$fastq_manifest"
  ln -sfn "$final_tsv" "$proj_dir/latest.tsv"
  ln -sfn "$fastq_manifest" "$proj_dir/latest.fastq_files.tsv"
  printf '%s\n' "READY $(ts)" > "$proj_dir/READY"

  local md5 bytes
  md5="$(md5_of "$final_tsv")"
  bytes="$(bytes_of "$final_tsv")"
  append_run_log "$proj_dir" "finalize" "ok" "md5" "$md5" "Saved TSV and FASTQ file manifest"
  append_global_manifest "$(ts)" "$project_id" "success" "$bytes" "$rows" "$md5" "$final_tsv" "OK"
  log "Project $project_id: rows=${rows} bytes=${bytes} md5=${md5}"
}

build_global_runs_manifest(){
  local tmp="$GLOBAL_RUNS.part"
  rm -f "$tmp"
  local wrote_header=0
  local tsv project_id
  shopt -s nullglob
  for tsv in "$META_ROOT"/*/*.tsv; do
    [[ "$(basename "$tsv")" == "latest.tsv" ]] && continue
    [[ "$(basename "$tsv")" == *.fastq_files.tsv ]] && continue
    [[ "$(basename "$tsv")" == "_run_log.tsv" ]] && continue
    [[ -s "$tsv" ]] || continue
    project_id="$(basename "$(dirname "$tsv")")"
    if (( ! wrote_header )); then
      { printf 'source_project_id\t'; head -n1 "$tsv" | tr -d '\r'; } > "$tmp"
      wrote_header=1
    fi
    awk -v p="$project_id" 'NR>1{print p "\t" $0}' "$tsv" >> "$tmp"
  done
  shopt -u nullglob
  if (( wrote_header )); then
    mv -f "$tmp" "$GLOBAL_RUNS"
  else
    write_tsv_line source_project_id > "$GLOBAL_RUNS"
    rm -f "$tmp"
  fi
}

# =========================
# Main
# =========================
main(){
  parse_args "$@"
  if (( DRY_RUN )); then
    log "DRY-RUN mode: no directories/files will be created and no ENA requests will be sent."
  else
    init_global_manifest
  fi
  log "Reading projects from CSV: $CSV_FILE"
  mapfile -t projects < <(read_projects_from_csv "$CSV_FILE")
  if (( ${#projects[@]} == 0 )); then
    log "No projects found in CSV."
    exit 0
  fi

  log "Found ${#projects[@]} unique accession/project IDs."
  local pid
  for pid in "${projects[@]}"; do
    log "==> $pid"
    fetch_project_tsv "$pid" || log "Project $pid failed or was skipped."
  done
  if (( DRY_RUN )); then
    log "Dry-run complete. Planned ${#projects[@]} accession/project IDs."
    log "No files were written. No ENA requests were sent."
  else
    build_global_runs_manifest
    log "Done. Project manifest: $GLOBAL_MANIFEST"
    log "Done. Combined run manifest: $GLOBAL_RUNS"
  fi
}

main "$@"