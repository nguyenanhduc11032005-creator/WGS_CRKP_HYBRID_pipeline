#!/bin/bash
#SBATCH --job-name=deep_dive
#SBATCH --output=_deepdive_%A_%a.out
#SBATCH --error=_deepdive_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-19%5

# ==============================================================================
# Script: 08_deep_dive_mechanisms.sh (Module 5 - Final Production Version)
# PURPOSE: Phân tích chuyên sâu PCN, IS Interruption và Porin (ompK35/36).
# ==============================================================================

set -Eeuo pipefail

trap 'echo "[FATAL] Script bị ngắt đột ngột tại dòng $LINENO!"; exit 1' ERR INT TERM

# ---------------------------------------------------------
# 0. CLI & THAM SỐ
# ---------------------------------------------------------
DRY_RUN=0
RESUME=1

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        --force) RESUME=0 ;;
    esac
done

THREADS="${SLURM_CPUS_PER_TASK:-16}"

ROOT_DIR="${ROOT_DIR:-/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID}"
CLEAN_IN="${CLEAN_IN:-$ROOT_DIR/04_clean_data}"
ASM_IN="${ASM_IN:-$ROOT_DIR/05_assembly}"
ANNOT_IN="${ANNOT_IN:-$ROOT_DIR/07_annotation}"

DEEP_OUT="${DEEP_OUT:-$ROOT_DIR/09_deep_dive_mechanisms}"
CHECKPOINT_DIR="$DEEP_OUT/.checkpoints"
SUMMARY="$DEEP_OUT/deep_dive_summary.tsv"

ENV_MAPPING="${ENV_MAPPING:-mapping_env}"
ENV_ISESCAN="${ENV_ISESCAN:-isescan_env}"

# ---------------------------------------------------------
# 1. HELPERS
# ---------------------------------------------------------
ts() { date '+%Y-%m-%d %H:%M:%S%z'; }
log() { printf '[%s] %s\n' "$(ts)" "$*" >&2; }
die() { log "FATAL: $*"; exit 1; }

check_program() {
    local env_name="$1"
    local prog="$2"
    if ! conda run -n "$env_name" command -v "$prog" >/dev/null 2>&1; then
        die "Thiếu dependency: $prog trong môi trường $env_name"
    fi
}

append_summary() {
    local line="$1"
    [[ "$DRY_RUN" -eq 1 ]] && return 0
    (
        flock -x 200
        [[ -f "$SUMMARY" ]] || echo -e "Sample\tGroup\tStatus\tChromosome_Depth\tPlasmids_PCN\tompK35_Status\tompK36_Status\tIS_Interrupting_Porin\tTime" > "$SUMMARY"
        echo -e "$line" >> "$SUMMARY"
    ) 200>"${SUMMARY}.lock"
}

if [[ "$DRY_RUN" -eq 0 ]]; then
    check_program "$ENV_MAPPING" "bwa"
    check_program "$ENV_MAPPING" "samtools"
    check_program "$ENV_ISESCAN" "isescan.py"
    check_program "$ENV_MAPPING" "python3"
fi

# ---------------------------------------------------------
# 2. MAIN EXECUTION PER SAMPLE
# ---------------------------------------------------------
run_sample() {
    local SPATH="$1"
    local GROUP="$(basename "$(dirname "$SPATH")")"
    local SAMPLE="$(basename "$SPATH")"
    
    local R1_IN="$CLEAN_IN/$GROUP/$SAMPLE/${SAMPLE}_R1.fastq.gz"
    local R2_IN="$CLEAN_IN/$GROUP/$SAMPLE/${SAMPLE}_R2.fastq.gz"
    local FASTA_IN="$SPATH/assembly.fasta"
    local GFF_IN="$ANNOT_IN/$GROUP/$SAMPLE/Bakta/${SAMPLE}.gff3"
    
    log "=================================================="
    log "ĐANG XỬ LÝ DEEP-DIVE: $SAMPLE (Nhóm: $GROUP)"

    if [[ ! -f "$R1_IN" || ! -f "$R2_IN" || ! -f "$FASTA_IN" || ! -f "$GFF_IN" ]]; then
        log "[ERROR] Thiếu dữ liệu đầu vào (R1/R2/Fasta/GFF3) cho $SAMPLE."
        append_summary "${SAMPLE}\t${GROUP}\tFAILED\t-\t-\t-\t-\t-\tMissing Inputs"
        return 1
    fi

    local SAMPLE_OUT="$DEEP_OUT/$GROUP/$SAMPLE"
    local PCN_DIR="$SAMPLE_OUT/PCN"
    local IS_DIR="$SAMPLE_OUT/ISEScan"
    local PORIN_DIR="$SAMPLE_OUT/Porin"
    local RUN_LOG="$SAMPLE_OUT/deep_dive.log"
    local CKPT_DONE="$CHECKPOINT_DIR/${GROUP}__${SAMPLE}.deepdive.ok"

    if [[ "$RESUME" -eq 1 && "$DRY_RUN" -eq 0 && -f "$CKPT_DONE" ]]; then
        log "[RESUME] Đã hoàn tất cho $SAMPLE -> Bỏ qua."
        return 0
    fi

    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$PCN_DIR" "$IS_DIR" "$PORIN_DIR"
        
        echo "Bắt đầu Deep-Dive cho $SAMPLE lúc $(date)" > "$RUN_LOG"
        echo "--- Tool Versions ---" >> "$RUN_LOG"
        conda run -n "$ENV_MAPPING" bwa 2>&1 | grep "Version:" >> "$RUN_LOG" || true
        conda run -n "$ENV_MAPPING" samtools 2>&1 | grep "Version:" | head -n 1 >> "$RUN_LOG" || true
        conda run -n "$ENV_ISESCAN" isescan.py --version >> "$RUN_LOG" 2>&1 || true
        echo "---------------------" >> "$RUN_LOG"
    fi

    local START_TIME="$(date +%s)"
    local FAIL=0

    # ---------------------------------------------------------
    # BƯỚC 2.1: BWA Index an toàn
    # ---------------------------------------------------------
    if [[ "$DRY_RUN" -eq 0 ]]; then
        (
            flock -x 201
            if [[ ! -f "${FASTA_IN}.bwt" || ! -f "${FASTA_IN}.sa" || ! -f "${FASTA_IN}.pac" || ! -f "${FASTA_IN}.ann" || ! -f "${FASTA_IN}.amb" ]]; then
                log "  -> Tạo BWA Index cho $SAMPLE..."
                conda run -n "$ENV_MAPPING" bwa index "$FASTA_IN" >> "$RUN_LOG" 2>&1
            fi
        ) 201>"${FASTA_IN}.lock"
    fi

    # Khởi chạy Song song PCN Mapping và ISEScan
    local BAM_OUT="$PCN_DIR/${SAMPLE}_mapped.bam"
    local DEPTH_OUT="$PCN_DIR/${SAMPLE}_depth.tsv"
    local PID_PCN=""
    local PID_IS=""

    if [[ "$DRY_RUN" -eq 0 ]]; then
        log "  -> [1/2] Chạy song song PCN Mapping và ISEScan..."
        
        # Process 1: Mapping
        (
            conda run -n "$ENV_MAPPING" bash -c "bwa mem -t $THREADS '$FASTA_IN' '$R1_IN' '$R2_IN' | samtools sort -@ $THREADS -o '$BAM_OUT' -"
            conda run -n "$ENV_MAPPING" samtools index "$BAM_OUT"
            conda run -n "$ENV_MAPPING" samtools coverage "$BAM_OUT" > "$DEPTH_OUT"
        ) >> "$RUN_LOG" 2>&1 || { FAIL=1; } &
        PID_PCN=$!

        # Process 2: ISEScan
        (
            conda run -n "$ENV_ISESCAN" isescan.py --seqfile "$FASTA_IN" --output "$IS_DIR" --nthread "$THREADS"
        ) >> "$RUN_LOG" 2>&1 || { FAIL=1; } &
        PID_IS=$!

        wait $PID_PCN || FAIL=1
        wait $PID_IS || FAIL=1
    fi

    # ---------------------------------------------------------
    # BƯỚC 2.2: Porin Extraction & IS Interruption
    # ---------------------------------------------------------
    local PORIN_REPORT="$PORIN_DIR/${SAMPLE}_porin_report.tsv"
    local TMP_PY="$PORIN_DIR/porin_analyzer.py"
    
    local CHR_DP="-"
    local PLASMIDS_PCN="-"
    local K35="-"
    local K36="-"
    local IS_RES="-"
    
    if [[ "$DRY_RUN" -eq 0 && "$FAIL" -eq 0 ]]; then
        log "  -> [2/2] Phân tích vùng gen ompK35/36, Promoter và giao cắt IS..."
        
        # Tạo file Python trực tiếp để tránh lỗi stdin của Conda run
        cat <<'PY' > "$TMP_PY"
import sys, csv, re
from pathlib import Path

gff_file, depth_file, is_dir, out_report, sample_name = map(Path, sys.argv[1:6])
sample_name = str(sample_name)

chrom_contig = None
porins = {'ompK35': None, 'ompK36': None}
largest_contig = None
max_len = 0

try:
    with gff_file.open("r") as f:
        for line in f:
            if line.startswith("##sequence-region"):
                parts = line.strip().split()
                c_id = parts[1]
                c_len = int(parts[3])
                if c_len > max_len:
                    max_len = c_len
                    largest_contig = c_id
                continue
            if line.startswith("#"): continue
            
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS": continue
            
            attr = parts[8]
            if re.search(r"gene=dnaA\b", attr, re.IGNORECASE):
                chrom_contig = parts[0]
            
            if re.search(r"gene=ompk35\b", attr, re.IGNORECASE) or re.search(r"Name=ompK35\b", attr, re.IGNORECASE):
                porins['ompK35'] = {'contig': parts[0], 'start': int(parts[3]), 'end': int(parts[4]), 'strand': parts[6]}
            elif re.search(r"gene=ompk36\b", attr, re.IGNORECASE) or re.search(r"Name=ompK36\b", attr, re.IGNORECASE):
                porins['ompK36'] = {'contig': parts[0], 'start': int(parts[3]), 'end': int(parts[4]), 'strand': parts[6]}
except Exception as e:
    print(f"[Python Error] Lỗi parse GFF3: {e}", file=sys.stderr)

if not chrom_contig:
    chrom_contig = largest_contig

chrom_depth = -1.0
plasmid_pcns = []

try:
    with depth_file.open("r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        contigs_cov = list(reader)
        
        for c in contigs_cov:
            if c['#rname'] == chrom_contig:
                chrom_depth = float(c['meandepth'])
                break
                
        if chrom_depth > 0:
            for c in contigs_cov:
                if c['#rname'] != chrom_contig and int(c['endpos']) > 1000:
                    pcn = float(c['meandepth']) / chrom_depth
                    plasmid_pcns.append(f"{c['#rname']}:{pcn:.1f}X")
except Exception as e:
    print(f"[Python Error] Lỗi tính PCN: {e}", file=sys.stderr)

pcn_str = ";".join(plasmid_pcns) if plasmid_pcns else "No_Plasmids"
chrom_depth_str = f"{chrom_depth:.1f}" if chrom_depth > 0 else "Error"

is_interruptions = []
# Tự động tìm file csv sinh ra bởi ISEScan
is_tsv_files = list(is_dir.glob("prediction/*.csv"))

if is_tsv_files:
    try:
        with is_tsv_files[0].open("r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                is_contig = row.get('seqID', '')
                is_start = int(row.get('isBegin', 0))
                is_end = int(row.get('isEnd', 0))
                is_fam = row.get('family', 'Unknown_IS')
                
                for p_name, p_data in porins.items():
                    if p_data and is_contig == p_data['contig']:
                        
                        if p_data['strand'] == '+':
                            prom_start = max(1, p_data['start'] - 500)
                            prom_end = p_data['start'] - 1
                        else:
                            prom_start = p_data['end'] + 1
                            prom_end = p_data['end'] + 500
                            
                        if (is_start <= p_data['end'] and is_end >= p_data['start']):
                            is_interruptions.append(f"{is_fam}_in_CDS_{p_name}")
                            
                        elif (is_start <= prom_end and is_end >= prom_start):
                            is_interruptions.append(f"{is_fam}_in_Promoter_{p_name}")
    except Exception as e:
        print(f"[Python Error] Lỗi phân tích ISEScan: {e}", file=sys.stderr)

k35_stat = "Intact" if porins['ompK35'] else "Missing_or_Truncated"
k36_stat = "Intact" if porins['ompK36'] else "Missing_or_Truncated"
is_res = ";".join(set(is_interruptions)) if is_interruptions else "None"

try:
    with out_report.open("w") as f:
        f.write(f"{chrom_depth_str}\t{pcn_str}\t{k35_stat}\t{k36_stat}\t{is_res}\n")
except Exception as e:
    print(f"[Python Error] Lỗi ghi report: {e}", file=sys.stderr)
PY

        # Chạy file Python đã tạo thông qua Conda
        conda run -n "$ENV_MAPPING" python3 "$TMP_PY" "$GFF_IN" "$DEPTH_OUT" "$IS_DIR" "$PORIN_REPORT" "$SAMPLE" >> "$RUN_LOG" 2>&1 || FAIL=1
        
        if [[ -s "$PORIN_REPORT" ]]; then
            IFS=$'\t' read -r CHR_DP PLASMIDS_PCN K35 K36 IS_RES < "$PORIN_REPORT"
        else
            FAIL=1
        fi
        
        # Dọn dẹp file tạm
        rm -f "$TMP_PY"
    fi

    local END_TIME="$(date +%s)"
    local ELAPSED=$((END_TIME - START_TIME))

    if [[ "$DRY_RUN" -eq 0 ]]; then
        if [[ "$FAIL" -eq 0 ]]; then
            append_summary "${SAMPLE}\t${GROUP}\tSUCCESS\t${CHR_DP}\t${PLASMIDS_PCN}\t${K35}\t${K36}\t${IS_RES}\t${ELAPSED}s"
            touch "$CKPT_DONE"
            log "  [SUCCESS] Deep-Dive hoàn tất cho $SAMPLE."
        else
            append_summary "${SAMPLE}\t${GROUP}\tFAILED\t${CHR_DP}\t${PLASMIDS_PCN}\t${K35}\t${K36}\t${IS_RES}\t${ELAPSED}s"
            log "  [ERROR] Lỗi xử lý, xem log: $RUN_LOG"
        fi
    fi
}

# ---------------------------------------------------------
# 3. CHUẨN BỊ MẢNG SAMPLES & THỰC THI
# ---------------------------------------------------------
main() {
    [[ -d "$ASM_IN" ]] || die "Không tìm thấy ASM_IN: $ASM_IN"
    
    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$DEEP_OUT" "$CHECKPOINT_DIR"
        (
            flock -x 200
            [[ -f "$SUMMARY" ]] || echo -e "Sample\tGroup\tStatus\tChromosome_Depth\tPlasmids_PCN\tompK35_Status\tompK36_Status\tIS_Interrupting_Porin\tTime" > "$SUMMARY"
        ) 200>"${SUMMARY}.lock"
    fi

    mapfile -t SAMPLE_PATHS < <(find "$ASM_IN" -mindepth 2 -maxdepth 2 -type d -exec test -f '{}/assembly.fasta' \; -print | sort)
    [[ ${#SAMPLE_PATHS[@]} -eq 0 ]] && die "Không tìm thấy assembly hợp lệ trong $ASM_IN"

    local TOTAL_SAMPLES="${#SAMPLE_PATHS[@]}"
    local TARGET_PATHS=()

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        local IDX=$((SLURM_ARRAY_TASK_ID - 1))
        (( IDX >= TOTAL_SAMPLES )) && exit 0
        TARGET_PATHS=("${SAMPLE_PATHS[$IDX]}")
        log "SLURM Array mode: Task $SLURM_ARRAY_TASK_ID / $TOTAL_SAMPLES"
    else
        TARGET_PATHS=("${SAMPLE_PATHS[@]}")
        log "Sequential mode: Chạy $TOTAL_SAMPLES mẫu"
    fi

    for spath in "${TARGET_PATHS[@]}"; do
        run_sample "$spath"
    done
}

main "$@"