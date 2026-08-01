#!/bin/bash
#SBATCH --job-name=annotation
#SBATCH --output=_annotation_%A_%a.out
#SBATCH --error=_annotation_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-19%5  

# ==============================================================================
# Script: 06_annotation_v2.sh (Module 3 - HIGH-RESOLUTION ANNOTATION)
# Production-Grade Version: Bền bỉ với SLURM array, chống silent errors.
# ==============================================================================

set -euo pipefail
trap 'echo "[FATAL] Lỗi tại dòng $LINENO"; exit 1' ERR

# ---------------------------------------------------------
# 0. XỬ LÝ THAM SỐ VÀ CẤU HÌNH CƠ BẢN
# ---------------------------------------------------------
DRY_RUN=0
RESUME=1
INCLUDE_FAILED_QC=0

for arg in "$@"; do
    case "$arg" in
        --dry-run)            DRY_RUN=1 ;;
        --force)              RESUME=0 ;;
        --include-failed-qc)  INCLUDE_FAILED_QC=1 ;;
        *) echo "[WARNING] Tham số không xác định: $arg (bỏ qua)" ;;
    esac
done

THREADS="${SLURM_CPUS_PER_TASK:-16}"

ROOT_DIR="/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID"
ASM_IN="$ROOT_DIR/05_assembly"
ASM_QC_SUMMARY="$ROOT_DIR/06_assembly_qc/assembly_qc_summary.tsv"
ANNOT_OUT="$ROOT_DIR/07_annotation"
CHECKPOINT_DIR="$ANNOT_OUT/.checkpoints"
SUMMARY="$ANNOT_OUT/annotation_summary.tsv"

ENV_BAKTA="bakta_env"
ENV_INFERNAL="infernal_env"

BAKTA_DB_DIR="${BAKTA_DB_DIR:-$ROOT_DIR/bakta_db/db}"
RFAM_DB_DIR="${RFAM_DB_DIR:-$ROOT_DIR/rfam_db}"
RFAM_CM="$RFAM_DB_DIR/Rfam.cm"
RFAM_CLANIN="$RFAM_DB_DIR/Rfam.clanin"

GENUS="${GENUS:-Klebsiella}"
SPECIES="${SPECIES:-pneumoniae}"
BAKTA_COMPLETE="${BAKTA_COMPLETE:-1}"
FLANK_BP="${FLANK_BP:-3000}"

# ---------------------------------------------------------
# 0.1 HÀM TIỆN ÍCH HỆ THỐNG
# ---------------------------------------------------------
require_cmds() {
    local missing=()
    for cmd in "$@"; do
        if ! command -v "$cmd" >/dev/null 2>&1; then
            missing+=("$cmd")
        fi
    done
    if (( ${#missing[@]} > 0 )); then
        echo "[FATAL] Thiếu lệnh hệ thống: ${missing[*]}"
        exit 1
    fi
}

conda_env_exists() {
    local env_name="$1"
    conda env list | awk -v env="$env_name" '$1 == env { found=1 } END { exit !found }'
}

run_tool() {
    local env_name="$1"
    local log_file="$2"
    shift 2
    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf '    [DRY-RUN] (%s) %q >> %q\n' "$env_name" "$*" "$log_file"
        return 0
    fi
    (
        eval "$(conda shell.bash hook)"
        conda activate "$env_name"
        "$@"
    ) >> "$log_file" 2>&1
}

run_tool_capture() {
    local env_name="$1"
    shift
    (
        eval "$(conda shell.bash hook)"
        conda activate "$env_name"
        "$@"
    )
}

append_summary() {
    local line="$1"
    if [[ "$DRY_RUN" -eq 1 ]]; then return 0; fi
    ( 
        flock -x 200
        printf '%b\n' "$line" >> "$SUMMARY" 
    ) 200>"${SUMMARY}.lock"
}

sanitize_locus_tag() {
    local raw="$1" tag
    tag=$(printf '%s' "$raw" | tr 'a-z' 'A-Z' | tr -cd 'A-Z0-9')
    tag="${tag:0:12}"
    [[ -n "$tag" ]] || tag="SAMPLE"
    printf '%s' "$tag"
}

# ---------------------------------------------------------
# 0.2 HÀM XỬ LÝ DỮ LIỆU CHUYÊN SÂU (BULLETPROOF VERSION)
# ---------------------------------------------------------
compute_cmscan_Z() {
    local fasta="$1"
    local total_nt
    
    # Dùng native awk đếm trực tiếp chiều dài hệ gen từ file FASTA (bỏ qua phụ thuộc esl-seqstat)
    total_nt=$(awk '/^>/ {next} {seqlen += length($0)} END {print seqlen}' "$fasta")
        
    if [[ -z "$total_nt" || ! "$total_nt" =~ ^[0-9]+$ || "$total_nt" -eq 0 ]]; then
        return 1
    fi
    python3 -c "print(f'{($total_nt * 2) / 1e6:.6f}')"
}

filter_cmscan_tblout() {
    local tblout="$1" out_tsv="$2"
    if [[ ! -s "$tblout" ]]; then
        : > "$out_tsv"
        return 1
    fi
    python3 - "$tblout" "$out_tsv" <<'PY'
import csv, sys
tblout, out_tsv = sys.argv[1:3]
cols = ['idx', 'target_name', 'target_acc', 'query_name', 'query_acc', 'clan',
        'mdl', 'mdl_from', 'mdl_to', 'seq_from', 'seq_to', 'strand', 'trunc',
        'pass', 'gc', 'bias', 'score', 'evalue', 'inc', 'olp', 'anyidx',
        'afrct1', 'afrct2', 'winidx', 'wfrct1', 'wfrct2', 'description']
rows = []
with open(tblout) as fh:
    for line in fh:
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.rstrip('\n').split(None, 25)
        if len(parts) != 26:
            continue
        rec = dict(zip(cols, parts))
        if rec.get('olp') == '=':
            continue
        rows.append(rec)

with open(out_tsv, 'w', newline='') as oh:
    w = csv.DictWriter(oh, fieldnames=cols, delimiter='\t', lineterminator='\n')
    w.writeheader()
    for r in rows:
        w.writerow(r)
PY
}

postprocess_annotation() {
    local sample="$1" bakta_tsv="$2" ncrna_tsv="$3" review_out="$4" flank="$5"
    python3 - "$sample" "$bakta_tsv" "$ncrna_tsv" "$review_out" "$flank" <<'PY'
import csv, re, sys

sample, bakta_tsv, ncrna_tsv, review_out, flank = sys.argv[1:6]
flank = int(flank)

counts = {'cds': 0, 'trna': 0, 'rrna': 0}
rep_candidates = []
rep_pattern = re.compile(
    r'(replication\s+initiation|replication\s+protein|^rep[abc]?$|'
    r'incompatibility\s+protein|iteron|ori[ct]|cop[at])',
    re.IGNORECASE,
)

header = None
try:
    with open(bakta_tsv, newline='') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if line.startswith('#'):
                if 'Sequence Id' in line and 'Type' in line:
                    header = line.lstrip('#').split('\t')
                continue
            if header is None or not line:
                continue
            row = line.split('\t')
            rec = dict(zip(header, row))
            ftype = (rec.get('Type') or '').lower()
            if ftype in counts:
                counts[ftype] += 1
            if ftype == 'cds':
                gene = rec.get('Gene', '') or ''
                product = rec.get('Product', '') or ''
                if rep_pattern.search(gene) or rep_pattern.search(product):
                    try:
                        start = int(rec.get('Start', 0))
                        stop = int(rec.get('Stop', 0))
                    except ValueError:
                        continue
                    rep_candidates.append({
                        'contig': rec.get('Sequence Id', ''),
                        'locus_tag': rec.get('Locus Tag', ''),
                        'gene': gene,
                        'product': product,
                        'start': start,
                        'stop': stop,
                        'strand': rec.get('Strand', ''),
                    })
except FileNotFoundError:
    pass

ncrna_hits = []
try:
    with open(ncrna_tsv, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            try:
                sf, st = int(r['seq_from']), int(r['seq_to'])
            except (KeyError, ValueError):
                continue
            ncrna_hits.append({
                'contig': r.get('query_name', ''),
                'rfam_id': r.get('target_acc', ''),
                'rfam_name': r.get('target_name', ''),
                'lo': min(sf, st),
                'hi': max(sf, st),
                'strand': r.get('strand', ''),
            })
except FileNotFoundError:
    pass

with open(review_out, 'w', newline='') as oh:
    out_cols = ['sample', 'contig', 'rep_locus_tag', 'rep_gene', 'rep_product',
                'rep_start', 'rep_stop', 'rep_strand', 'flank_start', 'flank_end',
                'ncrna_hits_in_flank', 'note']
    w = csv.DictWriter(oh, fieldnames=out_cols, delimiter='\t', lineterminator='\n')
    w.writeheader()
    for rc in rep_candidates:
        flo = max(1, rc['start'] - flank)
        fhi = rc['stop'] + flank
        overlaps = []
        for h in ncrna_hits:
            if h['contig'] != rc['contig']:
                continue
            if h['lo'] <= fhi and h['hi'] >= flo:
                overlaps.append(f"{h['rfam_id'] or h['rfam_name']}@{h['lo']}-{h['hi']}({h['strand']})")
        w.writerow({
            'sample': sample,
            'contig': rc['contig'],
            'rep_locus_tag': rc['locus_tag'],
            'rep_gene': rc['gene'],
            'rep_product': rc['product'],
            'rep_start': rc['start'],
            'rep_stop': rc['stop'],
            'rep_strand': rc['strand'],
            'flank_start': flo,
            'flank_end': fhi,
            'ncrna_hits_in_flank': ';'.join(overlaps) if overlaps else '-',
            'note': 'Can xac nhan thu cong (IGV/Artemis) cho copA/copT',
        })

print(f"{counts['cds']}\t{counts['trna']}\t{counts['rrna']}\t{len(rep_candidates)}")
PY
}

prepare_rfam_db() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "  [DRY-RUN] (kiểm tra/cmpress Rfam.cm nếu chưa index)"
        return 0
    fi
    if [[ ! -f "$RFAM_CM" ]]; then
        echo "[FATAL] Không tìm thấy $RFAM_CM."
        exit 1
    fi
    if [[ ! -f "$RFAM_CLANIN" ]]; then
        echo "[FATAL] Không tìm thấy $RFAM_CLANIN. Vui lòng tải file này về thư mục $RFAM_DB_DIR."
        exit 1
    fi
    (
        flock -x 201
        local need_index=0
        for suf in i1f i1i i1m i1p; do
            if [[ ! -f "${RFAM_CM}.${suf}" ]]; then
                need_index=1
                break
            fi
        done
        if [[ "$need_index" -eq 1 ]]; then
            echo "  -> Rfam.cm chưa đủ index, chạy cmpress..."
            run_tool "$ENV_INFERNAL" "$ANNOT_OUT/.rfam_cmpress.log" cmpress "$RFAM_CM"
        fi
    ) 201>"$RFAM_DB_DIR/.cmpress.lock"
}

# ---------------------------------------------------------
# 0.3 KHỞI TẠO MÔI TRƯỜNG & CHỐNG RACE-CONDITION CHUNG
# ---------------------------------------------------------
require_cmds awk sed python3 flock find sort wc mktemp

echo -e "\n=== KIỂM TRA MÔI TRƯỜNG CONDA ==="
if [[ "$DRY_RUN" -eq 0 ]]; then
    for env in "$ENV_BAKTA" "$ENV_INFERNAL"; do
        if ! conda_env_exists "$env"; then
            echo "[FATAL] Conda environment '${env}' không tồn tại!"
            exit 1
        fi
    done
    if [[ ! -d "$BAKTA_DB_DIR" ]]; then
        echo "[FATAL] Không tìm thấy Bakta database tại $BAKTA_DB_DIR."
        exit 1
    fi

    mkdir -p "$ANNOT_OUT" "$CHECKPOINT_DIR"
    (
        flock -x 200
        if [[ ! -f "$SUMMARY" ]]; then
            printf '%s\n' "Sample\tGroup\tStatus\tCDS\ttRNA\trRNA\tncRNA_Hits\tRepGene_Candidates\tTime\tMessage" > "$SUMMARY"
        fi
    ) 200>"${SUMMARY}.lock"
    
    prepare_rfam_db

    echo "Phần mềm sử dụng:"
    ( eval "$(conda shell.bash hook)"; conda activate "$ENV_BAKTA"; bakta --version || true )
    ( eval "$(conda shell.bash hook)"; conda activate "$ENV_INFERNAL"; cmscan -h 2>&1 | sed -n '2p' || true )
else
    mkdir -p "$ANNOT_OUT" 2>/dev/null || true
    prepare_rfam_db
fi

# ---------------------------------------------------------
# 1. QUÉT & KIỂM TRA CHẶT DỮ LIỆU ĐẦU VÀO
# ---------------------------------------------------------
if ! SAMPLE_LIST=$(find "$ASM_IN" -mindepth 2 -maxdepth 2 -type d -exec test -f '{}/assembly.fasta' \; -print | sort); then
    echo "[FATAL] Lỗi khi quét thư mục assembly tại $ASM_IN"
    exit 1
fi

mapfile -t SAMPLE_PATHS <<< "$SAMPLE_LIST"
TOTAL_SAMPLES=${#SAMPLE_PATHS[@]}

if [[ ${TOTAL_SAMPLES} -eq 0 || -z "${SAMPLE_PATHS[0]}" ]]; then
    echo "[ERROR] Không tìm thấy file assembly.fasta nào hợp lệ tại $ASM_IN."
    exit 1
fi

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    IDX=$((SLURM_ARRAY_TASK_ID - 1))
    if [[ $IDX -ge $TOTAL_SAMPLES ]]; then
        echo "[INFO] Task ID $SLURM_ARRAY_TASK_ID vượt quá số lượng mẫu ($TOTAL_SAMPLES)."
        exit 0
    fi
    TARGET_PATHS=("${SAMPLE_PATHS[$IDX]}")
    echo ">> SLURM Array Mode: Chạy task $SLURM_ARRAY_TASK_ID / $TOTAL_SAMPLES"
else
    TARGET_PATHS=("${SAMPLE_PATHS[@]}")
    echo ">> Sequential Mode: Chạy tuần tự $TOTAL_SAMPLES mẫu"
fi

# Fail-fast loop kiểm tra input ngay từ đầu
for spath in "${TARGET_PATHS[@]}"; do
    [[ -d "$spath" ]] || { echo "[ERROR] Không tồn tại thư mục: $spath"; exit 1; }
    [[ -s "$spath/assembly.fasta" ]] || { echo "[ERROR] Thiếu/rỗng assembly.fasta: $spath/assembly.fasta"; exit 1; }
done

# ---------------------------------------------------------
# 2. THỰC THI ANNOTATION (Bakta + Infernal/Rfam)
# ---------------------------------------------------------
for SPATH in "${TARGET_PATHS[@]}"; do
    GROUP=$(basename "$(dirname "$SPATH")")
    SAMPLE=$(basename "$SPATH")
    FASTA_IN="$SPATH/assembly.fasta"

    echo -e "\n=================================================="
    echo "ĐANG XỬ LÝ: $SAMPLE (Nhóm: $GROUP)"

    if [[ "$INCLUDE_FAILED_QC" -eq 0 && -f "$ASM_QC_SUMMARY" ]]; then
        QC_STATUS=$(awk -F'\t' -v s="$SAMPLE" '$1==s {print $3; exit}' "$ASM_QC_SUMMARY")
        if [[ "$QC_STATUS" == "FAILED" ]]; then
            echo "  [SKIP] Mẫu $SAMPLE có Assembly QC = FAILED. Dùng --include-failed-qc để ép chạy."
            continue
        fi
    fi

    SAMPLE_ANNOT_DIR="$ANNOT_OUT/$GROUP/$SAMPLE"
    BAKTA_OUT="$SAMPLE_ANNOT_DIR/Bakta"
    INFERNAL_OUT="$SAMPLE_ANNOT_DIR/Infernal"
    REVIEW_OUT="$SAMPLE_ANNOT_DIR/ManualReview"
    RUN_LOG="$SAMPLE_ANNOT_DIR/annotation.log"

    # Tách checkpoint linh hoạt
    CKPT_BAKTA="$CHECKPOINT_DIR/${GROUP}__${SAMPLE}.bakta.ok"
    CKPT_INFERNAL="$CHECKPOINT_DIR/${GROUP}__${SAMPLE}.infernal.ok"
    CKPT_DONE="$CHECKPOINT_DIR/${GROUP}__${SAMPLE}.annot.ok"

    if [[ "$RESUME" -eq 1 && -f "$CKPT_DONE" ]]; then
        echo "  [INFO] Đã hoàn tất Annotation trước đó (Resume). Bỏ qua."
        continue
    fi

    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$BAKTA_OUT" "$INFERNAL_OUT" "$REVIEW_OUT"
        echo "Bắt đầu Annotation cho $SAMPLE lúc $(date)" >> "$RUN_LOG"
    fi

    START_TIME=$(date +%s)
    ANNOT_FAIL=0
    LOCUS_TAG=$(sanitize_locus_tag "$SAMPLE")

    # --- 2.1 BAKTA ---
    if [[ "$RESUME" -eq 1 && -f "$CKPT_BAKTA" ]]; then
        echo "  [INFO] Bakta đã xong, bỏ qua bước này."
    else
        echo "  -> Chạy Bakta (locus-tag=$LOCUS_TAG)..."
        BAKTA_ARGS=(--db "$BAKTA_DB_DIR" --output "$BAKTA_OUT" --prefix "$SAMPLE"
                    --locus-tag "$LOCUS_TAG" --genus "$GENUS" --species "$SPECIES"
                    --keep-contig-headers --threads "$THREADS" --force)
        if [[ "$BAKTA_COMPLETE" -eq 1 ]]; then
            BAKTA_ARGS+=(--complete)
        fi
        if ! run_tool "$ENV_BAKTA" "$RUN_LOG" bakta "${BAKTA_ARGS[@]}" "$FASTA_IN"; then
            echo "  [ERROR] Bakta thất bại!"
            ANNOT_FAIL=1
        elif [[ "$DRY_RUN" -eq 0 ]]; then
            touch "$CKPT_BAKTA"
        fi
    fi

    # --- 2.2 INFERNAL (cmscan against Rfam) ---
    TBLOUT="$INFERNAL_OUT/${SAMPLE}.tblout"
    if [[ "$RESUME" -eq 1 && -f "$CKPT_INFERNAL" ]]; then
        echo "  [INFO] Infernal đã xong, bỏ qua bước này."
    elif [[ "$ANNOT_FAIL" -eq 0 ]]; then
        echo "  -> Chạy Infernal cmscan (Rfam)..."
        CMSCAN_STDOUT="$INFERNAL_OUT/${SAMPLE}.cmscan.out"
        if [[ "$DRY_RUN" -eq 1 ]]; then
            echo "    [DRY-RUN] esl-seqstat + cmscan --cpu $THREADS -Z <Mb> --cut_ga --rfam --nohmmonly --fmt 2 --clanin $RFAM_CLANIN -o $CMSCAN_STDOUT --tblout $TBLOUT $RFAM_CM $FASTA_IN"
        else
            Z_MB=""
            if ! Z_MB=$(compute_cmscan_Z "$FASTA_IN"); then
                echo "  [ERROR] Không tính được kích thước hệ gen (-Z) cho cmscan."
                ANNOT_FAIL=1
            elif ! run_tool "$ENV_INFERNAL" "$RUN_LOG" cmscan \
                    --cpu "$THREADS" -Z "$Z_MB" --cut_ga --rfam --nohmmonly \
                    --fmt 2 --clanin "$RFAM_CLANIN" \
                    -o "$CMSCAN_STDOUT" --tblout "$TBLOUT" \
                    "$RFAM_CM" "$FASTA_IN"; then
                echo "  [ERROR] Infernal cmscan thất bại!"
                ANNOT_FAIL=1
            else
                touch "$CKPT_INFERNAL"
            fi
        fi
    fi

    # --- 2.3 LỌC KẾT QUẢ CMSCAN & KHOANH VÙNG ỨNG VIÊN ---
    NCRNA_FILTERED="$INFERNAL_OUT/${SAMPLE}.ncRNA_filtered.tsv"
    REVIEW_TSV="$REVIEW_OUT/${SAMPLE}.ori_copAT_candidates.tsv"
    BAKTA_TSV="$BAKTA_OUT/${SAMPLE}.tsv"
    
    CDS_COUNT="-"; TRNA_COUNT="-"; RRNA_COUNT="-"; NCRNA_HITS="-"; REP_CANDIDATES="-"
    
    if [[ "$DRY_RUN" -eq 0 && "$ANNOT_FAIL" -eq 0 ]]; then
        if ! filter_cmscan_tblout "$TBLOUT" "$NCRNA_FILTERED"; then
            echo "  [ERROR] Lọc kết quả cmscan thất bại (tblout rỗng/thiếu)!"
            ANNOT_FAIL=1
        else
            NCRNA_HITS=$(( $(wc -l < "$NCRNA_FILTERED") - 1 ))
            (( NCRNA_HITS < 0 )) && NCRNA_HITS=0
            
            if COUNTS_LINE=$(postprocess_annotation "$SAMPLE" "$BAKTA_TSV" "$NCRNA_FILTERED" "$REVIEW_TSV" "$FLANK_BP"); then
                IFS=$'\t' read -r CDS_COUNT TRNA_COUNT RRNA_COUNT REP_CANDIDATES <<< "$COUNTS_LINE"
            else
                echo "  [ERROR] Bước tổng hợp/khoanh vùng ứng viên thất bại!"
                ANNOT_FAIL=1
            fi
        fi
    fi

    # --- SUMMARY GHI FILE ---
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))

    if [[ "$DRY_RUN" -eq 0 ]]; then
        if [[ "$ANNOT_FAIL" -eq 1 ]]; then
            append_summary "${SAMPLE}\t${GROUP}\tFAILED\t${CDS_COUNT}\t${TRNA_COUNT}\t${RRNA_COUNT}\t${NCRNA_HITS}\t${REP_CANDIDATES}\t${ELAPSED}s\tXem $RUN_LOG"
            continue
        fi
        
        if [[ -s "$BAKTA_TSV" && -s "$NCRNA_FILTERED" && -f "$REVIEW_TSV" ]]; then
            append_summary "${SAMPLE}\t${GROUP}\tSUCCESS\t${CDS_COUNT}\t${TRNA_COUNT}\t${RRNA_COUNT}\t${NCRNA_HITS}\t${REP_CANDIDATES}\t${ELAPSED}s\tOK"
            touch "$CKPT_DONE"
            echo "  [SUCCESS] Annotation hoàn tất cho $SAMPLE (${ELAPSED}s)"
            echo "      CDS: $CDS_COUNT | ncRNA hits: $NCRNA_HITS | Rep-gene candidates: $REP_CANDIDATES"
        else
            echo "  [ERROR] Thiếu output thật dù commands chạy thành công."
            append_summary "${SAMPLE}\t${GROUP}\tFAILED\t${CDS_COUNT}\t${TRNA_COUNT}\t${RRNA_COUNT}\t${NCRNA_HITS}\t${REP_CANDIDATES}\t${ELAPSED}s\tMissing expected outputs"
        fi
    fi
done

echo -e "\n=================================================="
echo " ĐÃ CHẠY XONG ANNOTATION (MODULE 3) CHO TASK NÀY"
echo "=================================================="
date