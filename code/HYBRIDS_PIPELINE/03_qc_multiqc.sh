#!/bin/bash
#SBATCH --job-name=hybrid_qc_prod
#SBATCH --output=_qc_%j.out
#SBATCH --error=_qc_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

# ===========================================================
#  QC PIPELINE CHO DỮ LIỆU HYBRID ILLUMINA+NANOPORE (PRODUCTION)
#  * Tính năng 1: Phân nhóm tự động lưu trực tiếp (Không Symlink)
#  * Tính năng 2: Đánh giá QC toàn diện Trước và Sau khi chạy Filtlong
# ===========================================================

set -euo pipefail

# ---------------------------------------------------------
# 0. XỬ LÝ THAM SỐ VÀ CẤU HÌNH CHUNG
# ---------------------------------------------------------
DRY_RUN=0
RESUME=0
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        --resume)  RESUME=1  ;;
        *) echo "[WARNING] Tham số không xác định: $arg (bỏ qua)" ;;
    esac
done

if [[ "$DRY_RUN" -eq 1 ]]; then
    echo -e "\n=================================================="
    echo "⚠️  DRY RUN MODE – không thực thi lệnh thật"
    echo "=================================================="
fi
if [[ "$RESUME" -eq 1 ]]; then
    echo "🔁  RESUME MODE – bỏ qua mẫu/nhánh đã hoàn tất"
fi

THREADS="${SLURM_CPUS_PER_TASK:-8}"
# Điều chỉnh tối ưu cho K. pneumoniae (~100x coverage)
TARGET_BASES_NANO="${TARGET_BASES_NANO:-550000000}"   
MIN_LEN_NANO="${MIN_LEN_NANO:-1000}"

# ---------------------------------------------------------
# CẤU HÌNH ĐƯỜNG DẪN TẬP TRUNG
# ---------------------------------------------------------
ROOT_DIR="/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID"
BASE_DIR="$ROOT_DIR/02_fastq/fastq/SRP573064"
QC_OUT="$ROOT_DIR/03_qc_results"
CLEAN_OUT="$ROOT_DIR/04_clean_data"
CHECKPOINT_DIR="$QC_OUT/.checkpoints"
LOCKFILE="$QC_OUT/.qc.lock"
METADATA_FILE="$ROOT_DIR/sample_metadata.tsv"

# ---------------------------------------------------------
# 0.1 FILE KHOÁ & KIỂM TRA METADATA
# ---------------------------------------------------------
if [[ ! -f "$METADATA_FILE" ]]; then
    echo "[FATAL] Missing metadata file at $METADATA_FILE."
    echo "Metadata is the Source of Truth. Vui lòng kiểm tra lại đường dẫn."
    exit 1
fi

if [[ "$DRY_RUN" -eq 0 ]]; then
    mkdir -p "$QC_OUT" 
    exec 200>"$LOCKFILE"
    if ! flock -n 200; then
        echo "[FATAL] Một tiến trình QC khác đang chạy. Dừng script."
        exit 1
    fi
    trap 'rm -f "$LOCKFILE"' EXIT
fi

# ---------------------------------------------------------
# 0.2 KÍCH HOẠT CONDA & HÀM TIỆN ÍCH
# ---------------------------------------------------------
echo -e "\n=== KIỂM TRA MÔI TRƯỜNG ==="
eval "$(conda shell.bash hook)"
if [[ "$DRY_RUN" -eq 0 ]]; then
    if ! conda activate wgs_crkp; then
        echo "[FATAL] Không thể activate conda env 'wgs_crkp'"
        exit 1
    fi
fi

execute() {
    if [[ "$DRY_RUN" -eq 1 ]]; then echo "    [DRY-RUN] $*"; return 0; fi
    "$@"
}

execute_quiet() {
    if [[ "$DRY_RUN" -eq 1 ]]; then echo "    [DRY-RUN] $*"; return 0; fi
    "$@" > /dev/null 2>&1
}

execute_log() {
    local logfile="$1"; shift
    if [[ "$DRY_RUN" -eq 1 ]]; then echo "    [DRY-RUN] $* > $logfile 2>&1"; return 0; fi
    "$@" > "$logfile" 2>&1
}

# Đọc an toàn mọi định dạng khoảng trắng và ký tự ẩn \r
get_sample_name() {
    awk -v s="$1" '
    { sub(/\r$/, "") }
    NR>1 && ($1==s || $2==s || $6==s || $7==s) {
        print $2; exit
    }' "$METADATA_FILE"
}

check_fastq_integrity() {
    local f="$1"
    if [[ "$f" == *.gz ]]; then gzip -t "$f" 2>/dev/null || return 1; fi
    seqkit head -n 100 "$f" >/dev/null 2>&1 || return 1
    return 0
}

require_single_file() {
    local desc="$1"        
    local -n arr=$2        
    if [[ ${#arr[@]} -ne 1 ]]; then
        echo "  [ERROR] Phát hiện ${#arr[@]} file $desc. Bỏ qua mẫu này."
        return 1
    fi
    return 0
}

# ---------------------------------------------------------
# 1. TÓM TẮT TRẠNG THÁI QC
# ---------------------------------------------------------
SUMMARY="$QC_OUT/qc_summary.tsv"
HEADER="Sample\tIllumina_QC\tNanopore_QC\tIllumina_reads_before\tIllumina_reads_after\tNanopore_bases_after"

if [[ "$DRY_RUN" -eq 0 ]]; then
    mkdir -p "$CHECKPOINT_DIR"
    if [[ "$RESUME" -eq 0 ]]; then
        rm -rf "$CHECKPOINT_DIR"/* 2>/dev/null || true
        echo -e "$HEADER" > "$SUMMARY"
    else
        [[ -f "$SUMMARY" ]] || echo -e "$HEADER" > "$SUMMARY"
    fi
fi

# ---------------------------------------------------------
# 2. VÒNG LẶP CHÍNH – XỬ LÝ TỪNG MẪU
# ---------------------------------------------------------
cd "$BASE_DIR" || exit 1

for DIR_ENTRY in */; do
    DIR_ENTRY="${DIR_ENTRY%/}"
    echo -e "\n=================================================="
    echo "ĐANG XỬ LÝ: $DIR_ENTRY"
    
    SAMPLE_NAME=$(get_sample_name "$DIR_ENTRY")
    if [[ -z "$SAMPLE_NAME" ]]; then
        echo "  [WARNING] Không tìm thấy metadata cho $DIR_ENTRY. Bỏ qua."
        continue
    fi
    echo "  -> MẪU: $SAMPLE_NAME"

    # PHÂN LOẠI NHÓM ĐỂ LƯU THƯ MỤC TRỰC TIẾP TỪ ĐẦU
    GROUP=$(awk -v s="$SAMPLE_NAME" 'NR>1 { sub(/\r$/, ""); if ($2==s) { print $3; exit } }' "$METADATA_FILE")
    if [[ "$GROUP" == "WildType" ]]; then
        GROUP_FOLDER="Clinical"
    elif [[ "$GROUP" == "Mutant" ]]; then
        GROUP_FOLDER="Mutant"
    else
        GROUP_FOLDER="Unassigned"
    fi
    echo "  -> NHÓM: $GROUP_FOLDER"

    ILLUMINA_CKPT="$CHECKPOINT_DIR/${SAMPLE_NAME}.illumina.ok"
    NANOPORE_CKPT="$CHECKPOINT_DIR/${SAMPLE_NAME}.nanopore.ok"

    # KHAI BÁO CẤU TRÚC ĐƯỜNG DẪN 
    SAMPLE_QC="$QC_OUT/$GROUP_FOLDER/$SAMPLE_NAME"
    SAMPLE_CLEAN="$CLEAN_OUT/$GROUP_FOLDER/$SAMPLE_NAME"

    ILL_RAW_QC="$SAMPLE_QC/Illumina/Raw/FastQC"
    ILL_FASTP_QC="$SAMPLE_QC/Illumina/Fastp"
    ILL_CLEAN_QC="$SAMPLE_QC/Illumina/Clean/FastQC"
    ILL_CLEAN_STATS="$SAMPLE_QC/Illumina/Clean/seqkit_stats.tsv"

    NANO_RAW_QC="$SAMPLE_QC/Nanopore/Raw/NanoPlot"
    NANO_PORECHOP_DIR="$SAMPLE_QC/Nanopore/Porechop"
    NANO_TRIMMED_QC="$SAMPLE_QC/Nanopore/Trimmed/NanoPlot"   # Thư mục cho QC trước Filtlong
    NANO_CLEAN_QC="$SAMPLE_QC/Nanopore/Clean/NanoPlot"
    NANO_CLEAN_STAT="$SAMPLE_QC/Nanopore/Clean/NanoStat.txt"
    NANO_CLEAN_STATS="$SAMPLE_QC/Nanopore/Clean/seqkit_stats.tsv"

    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$SAMPLE_CLEAN"
        mkdir -p "$ILL_RAW_QC" "$ILL_FASTP_QC" "$ILL_CLEAN_QC"
        mkdir -p "$NANO_RAW_QC" "$NANO_PORECHOP_DIR" "$NANO_TRIMMED_QC" "$NANO_CLEAN_QC"
    fi

    ILLUMINA_STATUS="PASS"
    NANOPORE_STATUS="PASS"
    ILL_B="-"
    ILL_A="-"
    NANO_A="-"

    # --- NHÁNH ILLUMINA ---
    if [[ "$RESUME" -eq 1 && "$DRY_RUN" -eq 0 && -f "$ILLUMINA_CKPT" ]]; then
        echo "  [Illumina] Đã hoàn tất (resume) -> bỏ qua"
    else
        ILLUMINA_PATH="${BASE_DIR}/${DIR_ENTRY}/illumina"
        if [[ -d "$ILLUMINA_PATH" ]]; then
            mapfile -t R1_FILES < <(find "$ILLUMINA_PATH" -maxdepth 1 -type f \( -name "*_1.fastq" -o -name "*_1.fastq.gz" \))
            mapfile -t R2_FILES < <(find "$ILLUMINA_PATH" -maxdepth 1 -type f \( -name "*_2.fastq" -o -name "*_2.fastq.gz" \))

            if ! require_single_file "R1" R1_FILES || ! require_single_file "R2" R2_FILES; then
                ILLUMINA_STATUS="FAIL"
            else
                R1="${R1_FILES[0]}"
                R2="${R2_FILES[0]}"

                if ! check_fastq_integrity "$R1" || ! check_fastq_integrity "$R2"; then
                    echo "  [ERROR] File Illumina không toàn vẹn."
                    ILLUMINA_STATUS="FAIL"
                else
                    CLEAN_R1="${SAMPLE_CLEAN}/${SAMPLE_NAME}_R1.fastq.gz"
                    CLEAN_R2="${SAMPLE_CLEAN}/${SAMPLE_NAME}_R2.fastq.gz"

                    if ! execute fastqc -t "$THREADS" -q -o "$ILL_RAW_QC" "$R1" "$R2"; then
                        ILLUMINA_STATUS="FAIL"
                    elif ! execute fastp -i "$R1" -I "$R2" \
                         -o "$CLEAN_R1" -O "$CLEAN_R2" \
                         -q 20 -u 20 -l 50 -3 -W 4 -M 20 \
                         --html "$ILL_FASTP_QC/${SAMPLE_NAME}_fastp.html" \
                         --json "$ILL_FASTP_QC/${SAMPLE_NAME}_fastp.json" \
                         --thread "$THREADS" 2>/dev/null; then
                        ILLUMINA_STATUS="FAIL"
                    elif ! execute fastqc -t "$THREADS" -q -o "$ILL_CLEAN_QC" "$CLEAN_R1" "$CLEAN_R2"; then
                        ILLUMINA_STATUS="FAIL"
                    elif ! execute seqkit stats -a "$CLEAN_R1" "$CLEAN_R2" -T -o "$ILL_CLEAN_STATS"; then
                        ILLUMINA_STATUS="FAIL"
                    fi
                fi
            fi
        else
            ILLUMINA_STATUS="MISSING"
        fi

        # Checkpoint xác thực
        if [[ "$DRY_RUN" -eq 0 && "$ILLUMINA_STATUS" == "PASS" ]]; then
            if [[ -s "$CLEAN_R1" && -s "$CLEAN_R2" && -s "$ILL_CLEAN_STATS" && -s "$ILL_FASTP_QC/${SAMPLE_NAME}_fastp.json" ]]; then
                touch "$ILLUMINA_CKPT"
            else
                echo "  [ERROR] Illumina QC báo PASS nhưng thiếu file output thật."
                ILLUMINA_STATUS="FAIL"
            fi
        fi
    fi

    # --- NHÁNH NANOPORE ---
    if [[ "$RESUME" -eq 1 && "$DRY_RUN" -eq 0 && -f "$NANOPORE_CKPT" ]]; then
        echo "  [Nanopore] Đã hoàn tất (resume) -> bỏ qua"
    else
        NANOPORE_PATH="${BASE_DIR}/${DIR_ENTRY}/nanopore"
        if [[ -d "$NANOPORE_PATH" ]]; then
            mapfile -t NANO_FILES < <(find "$NANOPORE_PATH" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \))

            if ! require_single_file "Nanopore" NANO_FILES; then
                NANOPORE_STATUS="FAIL"
            else
                NANO_FQ="${NANO_FILES[0]}"
                if ! check_fastq_integrity "$NANO_FQ"; then
                    echo "  [ERROR] File Nanopore không toàn vẹn."
                    NANOPORE_STATUS="FAIL"
                else
                    TRIMMED_FASTQ="${NANO_PORECHOP_DIR}/${SAMPLE_NAME}_porechop_trimmed.fastq"
                    PORECHOP_LOG="${NANO_PORECHOP_DIR}/${SAMPLE_NAME}_porechop.log"
                    FILTLONG_LOG="${NANO_PORECHOP_DIR}/${SAMPLE_NAME}_filtlong.log"
                    CLEAN_NANO="${SAMPLE_CLEAN}/${SAMPLE_NAME}_nanopore.fastq.gz"

                    if ! execute_quiet NanoPlot -t "$THREADS" --fastq "$NANO_FQ" -o "$NANO_RAW_QC"; then
                        NANOPORE_STATUS="FAIL"
                    elif ! execute_log "$PORECHOP_LOG" porechop -i "$NANO_FQ" -o "$TRIMMED_FASTQ" --threads "$THREADS"; then
                        NANOPORE_STATUS="FAIL"
                    else
                        # QC Dữ liệu Trung gian (Trước Filtlong)
                        echo "  [Nanopore] QC dữ liệu sau Porechop (Trước Filtlong)..."
                        if ! execute_quiet NanoPlot -t "$THREADS" --fastq "$TRIMMED_FASTQ" -o "$NANO_TRIMMED_QC"; then
                            echo "  [WARNING] Lỗi khi chạy NanoPlot cho dữ liệu Trimmed."
                        fi

                        # Chạy Filtlong
                        echo "  [Nanopore] Chạy Filtlong (min_len=${MIN_LEN_NANO}, target_bases=${TARGET_BASES_NANO})..."
                        if [[ "$DRY_RUN" -eq 1 ]]; then
                            echo "    [DRY-RUN] filtlong ... 2> $FILTLONG_LOG | gzip > $CLEAN_NANO"
                        else
                            if ! { filtlong --min_length "$MIN_LEN_NANO" --target_bases "$TARGET_BASES_NANO" "$TRIMMED_FASTQ" 2> "$FILTLONG_LOG" | gzip > "$CLEAN_NANO"; }; then
                                NANOPORE_STATUS="FAIL"
                            elif [[ ! -s "$CLEAN_NANO" ]]; then
                                NANOPORE_STATUS="FAIL"
                            fi
                        fi

                        # QC Dữ liệu Cuối cùng (Sau Filtlong)
                        if [[ "$NANOPORE_STATUS" == "PASS" ]]; then
                            echo "  [Nanopore] QC dữ liệu sau Filtlong..."
                            if ! execute_quiet NanoPlot -t "$THREADS" --fastq "$CLEAN_NANO" -o "$NANO_CLEAN_QC"; then
                                NANOPORE_STATUS="FAIL"
                            elif [[ "$DRY_RUN" -eq 0 ]] && ! NanoStat --fastq "$CLEAN_NANO" -t "$THREADS" > "$NANO_CLEAN_STAT"; then
                                NANOPORE_STATUS="FAIL"
                            elif [[ "$DRY_RUN" -eq 0 ]] && ! seqkit stats -a "$CLEAN_NANO" -T -o "$NANO_CLEAN_STATS"; then
                                NANOPORE_STATUS="FAIL"
                            fi
                        fi
                    fi
                    # Dọn rác
                    [[ "$DRY_RUN" -eq 0 && -f "$TRIMMED_FASTQ" ]] && rm -f "$TRIMMED_FASTQ"
                fi
            fi
        else
            NANOPORE_STATUS="MISSING"
        fi

        # Checkpoint xác thực
        if [[ "$DRY_RUN" -eq 0 && "$NANOPORE_STATUS" == "PASS" ]]; then
            if [[ -s "$CLEAN_NANO" && -s "$NANO_CLEAN_STAT" && -s "$NANO_CLEAN_STATS" ]]; then
                touch "$NANOPORE_CKPT"
            else
                echo "  [ERROR] Nanopore QC báo PASS nhưng thiếu file output thật."
                NANOPORE_STATUS="FAIL"
            fi
        fi
    fi

    # --- TRÍCH XUẤT STATS ---
    if [[ "$DRY_RUN" -eq 0 ]]; then
        # Lấy thông số reads
        if [[ -f "$ILL_FASTP_QC/${SAMPLE_NAME}_fastp.json" ]]; then
            ILL_B=$(grep '"total_reads"' "$ILL_FASTP_QC/${SAMPLE_NAME}_fastp.json" | head -n 1 | tr -d -c '0-9')
            ILL_A=$(grep '"total_reads"' "$ILL_FASTP_QC/${SAMPLE_NAME}_fastp.json" | head -n 2 | tail -n 1 | tr -d -c '0-9')
        fi
        if [[ -f "$NANO_CLEAN_STATS" ]]; then
            NANO_A=$(awk -F'\t' 'NR==2 {print $5}' "$NANO_CLEAN_STATS")
        fi

        # Cập nhật vào QC Summary
        TMP_SUMMARY=$(mktemp)
        awk -F'\t' -v s="$SAMPLE_NAME" '$1 != s' "$SUMMARY" > "$TMP_SUMMARY"
        mv "$TMP_SUMMARY" "$SUMMARY"
        echo -e "${SAMPLE_NAME}\t${ILLUMINA_STATUS}\t${NANOPORE_STATUS}\t${ILL_B}\t${ILL_A}\t${NANO_A}" >> "$SUMMARY"
    fi
done

# ---------------------------------------------------------
# 3. CHẠY MULTIQC SAU KHI HOÀN TẤT VÒNG LẶP
# ---------------------------------------------------------
if [[ "$DRY_RUN" -eq 0 ]]; then
    echo -e "\n=== CHẠY MULTIQC ==="
    rm -rf "$QC_OUT/MultiQC"
    mkdir -p "$QC_OUT/MultiQC"
    
    # Không ignore nhánh Clinical và Mutant nữa
    multiqc "$QC_OUT" -o "$QC_OUT/MultiQC" -n CRKP_Hybrid_QC_Report \
        --quiet --force --ignore .checkpoints \
        --ignore MultiQC --ignore software_versions.txt --ignore qc_summary.tsv

    echo -e "\n=================================================="
    echo " HOÀN TẤT QC PIPELINE"
    echo "=================================================="
    cat "$SUMMARY" | column -t
fi
date