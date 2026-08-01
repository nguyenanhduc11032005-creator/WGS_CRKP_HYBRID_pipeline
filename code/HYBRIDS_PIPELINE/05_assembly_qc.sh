#!/bin/bash
#SBATCH --job-name=asm_qc
#SBATCH --output=_asm_qc_%A_%a.out
#SBATCH --error=_asm_qc_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-19%5

# ==============================================================================
# Script: 05_assembly_qc.sh (Final Production - Offline BUSCO Support)
#
# PURPOSE:
# Đánh giá chất lượng bộ gen lai (Hybrid Assembly) bằng QUAST, BUSCO, và CheckM2.
# ==============================================================================

set -euo pipefail

# ---------------------------------------------------------
# 0. XỬ LÝ THAM SỐ VÀ CẤU HÌNH CƠ BẢN
# ---------------------------------------------------------
DRY_RUN=0
RESUME=1

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        --force)   RESUME=0  ;;
        *) echo "[WARNING] Tham số không xác định: $arg (bỏ qua)" ;;
    esac
done

THREADS="${SLURM_CPUS_PER_TASK:-16}"

ROOT_DIR="/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID"
ASM_IN="$ROOT_DIR/05_assembly"
QC_OUT="$ROOT_DIR/06_assembly_qc"
CHECKPOINT_DIR="$QC_OUT/.checkpoints"
SUMMARY="$QC_OUT/assembly_qc_summary.tsv"

# Khai báo tên môi trường Conda
ENV_QUAST="quast"
ENV_BUSCO="busco_env"
ENV_CHECKM2="checkm2"

BUSCO_LINEAGE="enterobacterales_odb10"

# ---------------------------------------------------------
# 0.1 HÀM TIỆN ÍCH CHẠY LỆNH VÀ GHI LOG TRONG SUBSHELL
# ---------------------------------------------------------
run_tool() {
    local env_name="$1"
    local log_file="$2"
    shift 2

    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf '    [DRY-RUN] (%s) ' "$env_name"
        printf '%q ' "$@"
        printf ' >> %q\n' "$log_file"
        return 0
    fi

    (
        eval "$(conda shell.bash hook)"
        conda activate "$env_name"
        "$@"
    ) >> "$log_file" 2>&1
}

append_summary() {
    local line="$1"
    if [[ "$DRY_RUN" -eq 1 ]]; then return 0; fi
    ( flock -x 200; echo -e "$line" >> "$SUMMARY" ) 200>"${SUMMARY}.lock"
}

# ---------------------------------------------------------
# 0.2 KHỞI TẠO MÔI TRƯỜNG & CHỐNG RACE-CONDITION CHUNG
# ---------------------------------------------------------
echo -e "\n=== KIỂM TRA MÔI TRƯỜNG CONDA ==="
eval "$(conda shell.bash hook)"

if [[ "$DRY_RUN" -eq 0 ]]; then
    for env in "$ENV_QUAST" "$ENV_BUSCO" "$ENV_CHECKM2"; do
        if ! conda env list | grep -E -q "(^|[[:space:]])${env}([[:space:]]|$)"; then
            echo "[FATAL] Conda environment '${env}' không tồn tại!"
            exit 1
        fi
    done

    mkdir -p "$CHECKPOINT_DIR"

    (
        flock -x 200
        if [[ ! -f "$SUMMARY" ]]; then
            echo -e "Sample\tGroup\tStatus\tQUAST_N50\tQUAST_Contigs\tBUSCO_Completeness\tCheckM2_Completeness\tCheckM2_Contamination\tMessage" > "$SUMMARY"
        fi
    ) 200>"${SUMMARY}.lock"
    
    echo "Phần mềm sử dụng:"
    ( eval "$(conda shell.bash hook)"; conda activate "$ENV_QUAST"; quast.py --version || true )
    ( eval "$(conda shell.bash hook)"; conda activate "$ENV_BUSCO"; busco --version || true )
    ( eval "$(conda shell.bash hook)"; conda activate "$ENV_CHECKM2"; checkm2 --version || true )
fi

# ---------------------------------------------------------
# 1. QUÉT DỮ LIỆU ĐẦU VÀO TỪ UNICYCLER
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

# ---------------------------------------------------------
# 2. THỰC THI ASSEMBLY QC ĐA MÔI TRƯỜNG
# ---------------------------------------------------------
for SPATH in "${TARGET_PATHS[@]}"; do
    GROUP=$(basename "$(dirname "$SPATH")")
    SAMPLE=$(basename "$SPATH")
    FASTA_IN="$SPATH/assembly.fasta"
    
    echo -e "\n=================================================="
    echo "ĐANG XỬ LÝ: $SAMPLE (Nhóm: $GROUP)"
    
    SAMPLE_QC_DIR="$QC_OUT/$GROUP/$SAMPLE"
    QUAST_OUT="$SAMPLE_QC_DIR/QUAST"
    BUSCO_OUT="$SAMPLE_QC_DIR/BUSCO"
    CHECKM2_OUT="$SAMPLE_QC_DIR/CheckM2"
    RUN_LOG="$SAMPLE_QC_DIR/assembly_qc.log"
    
    CKPT_FILE="$CHECKPOINT_DIR/${GROUP}__${SAMPLE}.asm_qc.ok"

    if [[ "$RESUME" -eq 1 && -f "$CKPT_FILE" ]]; then
        echo "  [INFO] Đã hoàn tất Assembly QC trước đó (Resume). Bỏ qua."
        continue
    fi

    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$QUAST_OUT" "$BUSCO_OUT" "$CHECKM2_OUT"
        echo "Bắt đầu QC cho $SAMPLE lúc $(date)" > "$RUN_LOG"
    fi

    QC_FAIL=0

    # 2.1 QUAST
    echo "  -> Chạy QUAST..."
    if ! run_tool "$ENV_QUAST" "$RUN_LOG" quast.py "$FASTA_IN" -o "$QUAST_OUT" -t "$THREADS" --silent; then
        echo "  [ERROR] QUAST thất bại!"
        QC_FAIL=1
    fi

    # 2.2 BUSCO (Đã thêm cơ chế --offline và trỏ đúng đường dẫn busco_downloads)
    echo "  -> Chạy BUSCO (Lineage: $BUSCO_LINEAGE)..."
    if ! run_tool "$ENV_BUSCO" "$RUN_LOG" busco -i "$FASTA_IN" -o busco_result --out_path "$BUSCO_OUT" -l "$BUSCO_LINEAGE" -m genome -c "$THREADS" -f -q --offline --download_path "$ROOT_DIR/busco_downloads"; then
        echo "  [ERROR] BUSCO thất bại!"
        QC_FAIL=1
    fi

    # 2.3 CheckM2
    echo "  -> Chạy CheckM2..."
    if ! run_tool "$ENV_CHECKM2" "$RUN_LOG" checkm2 predict --threads "$THREADS" -i "$FASTA_IN" -o "$CHECKM2_OUT" -x fasta --force; then
        echo "  [ERROR] CheckM2 thất bại!"
        QC_FAIL=1
    fi

    # ---------------------------------------------------------
    # 3. TRÍCH XUẤT THÔNG SỐ VÀ GHI LOG
    # ---------------------------------------------------------
    if [[ "$DRY_RUN" -eq 0 ]]; then
        if [[ $QC_FAIL -eq 1 ]]; then
            append_summary "${SAMPLE}\t${GROUP}\tFAILED\t-\t-\t-\t-\t-\tXem file log"
            continue
        fi

        VAL_N50="-"
        VAL_CONTIGS="-"
        if [[ -f "$QUAST_OUT/report.tsv" ]]; then
            VAL_N50=$(awk -F'\t' '$1 == "N50" {print $2}' "$QUAST_OUT/report.tsv")
            VAL_CONTIGS=$(awk -F'\t' '$1 == "# contigs" {print $2}' "$QUAST_OUT/report.tsv")
        fi

        VAL_BUSCO="-"
        BUSCO_TXT=$(find "$BUSCO_OUT" -name "short_summary*.txt" -print -quit 2>/dev/null || true)
        if [[ -n "$BUSCO_TXT" && -f "$BUSCO_TXT" ]]; then
            VAL_BUSCO=$(grep -oE 'C:[0-9.]+%' "$BUSCO_TXT" | head -n 1 | sed 's/C://')
        fi

        VAL_CHK_COMP="-"
        VAL_CHK_CONT="-"
        if [[ -f "$CHECKM2_OUT/quality_report.tsv" ]]; then
            VAL_CHK_COMP=$(awk -F'\t' 'NR==2 {print $2}' "$CHECKM2_OUT/quality_report.tsv")
            VAL_CHK_CONT=$(awk -F'\t' 'NR==2 {print $3}' "$CHECKM2_OUT/quality_report.tsv")
        fi

        append_summary "${SAMPLE}\t${GROUP}\tSUCCESS\t${VAL_N50}\t${VAL_CONTIGS}\t${VAL_BUSCO}\t${VAL_CHK_COMP}\t${VAL_CHK_CONT}\tOK"
        touch "$CKPT_FILE"
        echo "  [SUCCESS] QC hoàn tất cho $SAMPLE"
        echo "      QUAST N50: $VAL_N50 | BUSCO: $VAL_BUSCO | CheckM2 C/C: $VAL_CHK_COMP/$VAL_CHK_CONT"
    fi
done

echo -e "\n=================================================="
echo " ĐÃ CHẠY XONG ASSEMBLY QC CHO TASK NÀY"
echo "=================================================="
date