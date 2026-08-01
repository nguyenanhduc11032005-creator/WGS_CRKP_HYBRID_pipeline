#!/bin/bash
#SBATCH --job-name=unicycler_hybrid
#SBATCH --output=_unicycler_%A_%a.out
#SBATCH --error=_unicycler_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --array=1-19%5  # Cập nhật theo 19 mẫu thực tế (15 Clinical + 4 Mutant)

# ==============================================================================
# Script: 04_unicycler.sh (Production Version - Patched)
#
# PURPOSE:
# Thực hiện Hybrid Assembly bằng Unicycler (Illumina + Nanopore).
# Script tích hợp kiểm tra môi trường chặt chẽ, an toàn với SLURM Array, 
# chống race-condition tuyệt đối và fix lỗi Redirection khi chạy Dry-Run.
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
CLEAN_IN="$ROOT_DIR/04_clean_data"
ASM_OUT="$ROOT_DIR/05_assembly"
CHECKPOINT_DIR="$ASM_OUT/.checkpoints"
SUMMARY="$ASM_OUT/unicycler_summary.tsv"
CONDA_ENV_NAME="wgs_crkp"

# ---------------------------------------------------------
# 0.1 KHỞI TẠO MÔI TRƯỜNG & CHỐNG RACE-CONDITION FILE SUMMARY
# ---------------------------------------------------------
echo -e "\n=== KIỂM TRA MÔI TRƯỜNG & CONDA ==="
eval "$(conda shell.bash hook)"

if [[ "$DRY_RUN" -eq 0 ]]; then
    # Kiểm tra rõ ràng sự tồn tại của Conda Env trước khi activate
    if ! conda env list | grep -E -q "(^|[[:space:]])${CONDA_ENV_NAME}([[:space:]]|$)"; then
        echo "[FATAL] Conda environment '${CONDA_ENV_NAME}' không tồn tại trong hệ thống!"
        exit 1
    fi

    if ! conda activate "$CONDA_ENV_NAME"; then
        echo "[FATAL] Không thể activate conda env '${CONDA_ENV_NAME}'."
        exit 1
    fi
    
    mkdir -p "$CHECKPOINT_DIR"

    # Tạo Header Atomic với flock chống ghi chồng giữa các Task Array
    (
        flock -x 200
        if [[ ! -f "$SUMMARY" ]]; then
            echo -e "Sample\tGroup\tStatus\tFasta_Size\tExecution_Time\tMessage" > "$SUMMARY"
        fi
    ) 200>"${SUMMARY}.lock"

    # In và ghi log phiên bản phần mềm
    echo "Phần mềm sử dụng:"
    echo " - Unicycler version: $(unicycler --version 2>&1 || echo 'N/A')"
    echo " - SPAdes version:    $(spades.py --version 2>&1 || echo 'N/A')"
fi

execute() {
    if [[ "$DRY_RUN" -eq 1 ]]; then echo "    [DRY-RUN] $*"; return 0; fi
    "$@"
}

append_summary() {
    local line="$1"
    if [[ "$DRY_RUN" -eq 1 ]]; then return 0; fi
    (
        flock -x 200
        echo -e "$line" >> "$SUMMARY"
    ) 200>"${SUMMARY}.lock"
}

# ---------------------------------------------------------
# 1. QUÉT DỮ LIỆU ĐẦU VÀO (An toàn tuyệt đối với pipefail)
# ---------------------------------------------------------
# Tránh lỗi ngầm từ process substitution <(find ...)
if ! SAMPLE_LIST=$(find "$CLEAN_IN" -mindepth 2 -maxdepth 2 -type d | sort); then
    echo "[FATAL] Lỗi khi thực hiện lệnh 'find' quét dữ liệu tại $CLEAN_IN"
    exit 1
fi

mapfile -t SAMPLE_PATHS <<< "$SAMPLE_LIST"
TOTAL_SAMPLES=${#SAMPLE_PATHS[@]}

if [[ ${TOTAL_SAMPLES} -eq 0 || -z "${SAMPLE_PATHS[0]}" ]]; then
    echo "[ERROR] Không tìm thấy thư mục mẫu hợp lệ nào tại $CLEAN_IN."
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
# 2. THỰC THI UNICYCLER
# ---------------------------------------------------------
for SPATH in "${TARGET_PATHS[@]}"; do
    GROUP=$(basename "$(dirname "$SPATH")")
    SAMPLE=$(basename "$SPATH")
    
    echo -e "\n=================================================="
    echo "ĐANG XỬ LÝ: $SAMPLE (Nhóm: $GROUP)"
    
    # Định vị Input
    R1="$SPATH/${SAMPLE}_R1.fastq.gz"
    R2="$SPATH/${SAMPLE}_R2.fastq.gz"
    NANO="$SPATH/${SAMPLE}_nanopore.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" || ! -f "$NANO" ]]; then
        echo "  [ERROR] Thiếu file input (.fastq.gz) cho $SAMPLE."
        append_summary "${SAMPLE}\t${GROUP}\tFAILED\t0\t-\tMissing Input Files"
        continue
    fi

    # Định vị Output & Log
    SAMPLE_ASM_DIR="$ASM_OUT/$GROUP/$SAMPLE"
    ASM_FASTA="$SAMPLE_ASM_DIR/assembly.fasta"
    RUN_LOG="$SAMPLE_ASM_DIR/unicycler_command.log"
    CKPT_FILE="$CHECKPOINT_DIR/${SAMPLE}.unicycler.ok"

    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$SAMPLE_ASM_DIR"
    fi

    if [[ "$RESUME" -eq 1 && -f "$CKPT_FILE" ]]; then
        echo "  [INFO] Đã hoàn tất trước đó (Resume). Bỏ qua mẫu này."
        continue
    fi

    # Cấu hình lệnh chạy
    CMD="unicycler -1 $R1 -2 $R2 -l $NANO -o $SAMPLE_ASM_DIR -t $THREADS --keep 1"
    
    if [[ "$DRY_RUN" -eq 0 ]]; then
        {
            echo "=========================================="
            echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
            echo "Directory: $PWD"
            echo "Command: $CMD"
            echo "Unicycler Version: $(unicycler --version 2>&1)"
            echo "=========================================="
        } > "$RUN_LOG"
    fi

    echo "  -> Bắt đầu Hybrid Assembly..."
    START_TIME=$(date +%s)
    
    # Thực thi Unicycler (Đã tách logic Dry Run để tránh lỗi Redirection)
    if [[ "$DRY_RUN" -eq 1 ]]; then
        execute $CMD
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
    else
        if ! execute $CMD >> "$RUN_LOG" 2>&1; then
            END_TIME=$(date +%s)
            ELAPSED=$((END_TIME - START_TIME))
            
            # Bẫy lỗi OOM (Out Of Memory) hoặc crash ngầm từ Unicycler/SPAdes
            LOG_TAIL=$(tail -n 5 "$RUN_LOG" | tr '\n' ' ' | tr '\t' ' ')
            echo "  [ERROR] Unicycler crash cho $SAMPLE! Chi tiết lỗi dòng cuối: $LOG_TAIL"
            append_summary "${SAMPLE}\t${GROUP}\tFAILED\t0\t${ELAPSED}s\tError: $LOG_TAIL"
            continue
        fi
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
    fi

    if [[ "$DRY_RUN" -eq 0 ]]; then
        if [[ -s "$ASM_FASTA" ]]; then
            FASTA_SIZE=$(du -h "$ASM_FASTA" | cut -f1)
            append_summary "${SAMPLE}\t${GROUP}\tSUCCESS\t${FASTA_SIZE}\t${ELAPSED}s\tOK"
            touch "$CKPT_FILE"
            echo "  [SUCCESS] Assembly hoàn tất cho $SAMPLE ($FASTA_SIZE) trong ${ELAPSED}s"
        else
            echo "  [ERROR] Chạy xong nhưng không tìm thấy assembly.fasta cho $SAMPLE."
            append_summary "${SAMPLE}\t${GROUP}\tFAILED\t0\t${ELAPSED}s\tMissing assembly.fasta"
        fi
    fi
done

echo -e "\n=================================================="
echo " ĐÃ CHẠY XONG UNICYCLER CHO TASK NÀY"
echo "=================================================="
date