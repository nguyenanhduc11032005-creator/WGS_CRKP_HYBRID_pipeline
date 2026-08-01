#!/bin/bash
#SBATCH --job-name=var_call_pair
#SBATCH --output=_varcall_pair_%A_%a.out
#SBATCH --error=_varcall_pair_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-4

# ==============================================================================
# Script: 09_variant_calling_pairwise.sh (Module 6)
# Chiến lược: Pairwise Comparison (Native vs Mutant)
# CẬP NHẬT CUỐI: Tăng RAM lên 64GB, cấp 50GB cho SnpEff, chặn Race Condition (-noStats), định tuyến thư mục tmp.
# ==============================================================================

set -Eeuo pipefail
trap 'echo "[FATAL] Script bị ngắt đột ngột tại dòng $LINENO! Xem file slurm-*.out/err để biết chi tiết." >&2; exit 1' ERR INT TERM

# ---------------------------------------------------------
# KIỂM TRA CHẾ ĐỘ DRY-RUN
# ---------------------------------------------------------
DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    echo "[DRY-RUN MODE] Đang ở chế độ giả lập. Không có lệnh thực tế nào được thực thi hay ghi file."
fi

THREADS="${SLURM_CPUS_PER_TASK:-16}"
PAIRS_FILE="pairs.tsv"

ROOT_DIR="/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID"
CLEAN_IN="$ROOT_DIR/04_clean_data"
ANNOT_IN="$ROOT_DIR/07_annotation"
VAR_OUT="$ROOT_DIR/10_variant_calling_pairwise"
CHECKPOINT_DIR="$VAR_OUT/.checkpoints"
SUMMARY="$VAR_OUT/pairwise_summary.tsv"

ENV_MAPPING="mapping_env"       
ENV_CALLING="freebayes_env"     
ENV_SNPEFF="snpeff_env"         

# ---------------------------------------------------------
# HÀM TIỆN ÍCH
# ---------------------------------------------------------
ts() { date '+%Y-%m-%d %H:%M:%S%z'; }
log() { printf '[%s] %s\n' "$(ts)" "$*" >&2; }
die() { log "FATAL: $*"; exit 1; }

run_tool() {
    local env_name="$1"
    local log_file="$2"
    shift 2
    
    if [[ "$DRY_RUN" == "true" ]]; then
        log "[DRY-RUN] [conda activate $env_name] $*"
    else
        (
            eval "$(conda shell.bash hook)"
            conda activate "$env_name"
            "$@"
        ) >> "$log_file" 2>&1
    fi
}

append_summary() {
    local line="$1"
    if [[ "$DRY_RUN" == "true" ]]; then
        log "[DRY-RUN] [SUMMARY MOCK] $line"
        return
    fi
    (
        flock -x 200
        [[ -f "$SUMMARY" ]] || echo -e "Native\tMutant\tStatus\tIll_Mapped\tONT_Mapped\tTotal_SNPs\tTotal_Indels\tcopA_Mutations\tTime" > "$SUMMARY"
        echo -e "$line" >> "$SUMMARY"
    ) 200>"${SUMMARY}.lock"
}

# ---------------------------------------------------------
# PARSE TỆP ÁNH XẠ THEO SLURM ARRAY
# ---------------------------------------------------------
[[ -f "$PAIRS_FILE" ]] || die "Không tìm thấy $PAIRS_FILE"
TASK_IDX=${SLURM_ARRAY_TASK_ID:-1}

PAIR_INFO=$(sed -n "${TASK_IDX}p" "$PAIRS_FILE")
[[ -z "$PAIR_INFO" ]] && die "Dòng $TASK_IDX trong $PAIRS_FILE trống."

read -r REF_GROUP REF_SAMPLE MUT_GROUP MUT_SAMPLE <<< "$PAIR_INFO"

log "=================================================="
log "TASK $TASK_IDX: So sánh $MUT_SAMPLE (Mutant) với $REF_SAMPLE (Native)"

# ---------------------------------------------------------
# KIỂM TRA DỮ LIỆU ĐẦU VÀO
# ---------------------------------------------------------
REF_FASTA="$ANNOT_IN/$REF_GROUP/$REF_SAMPLE/Bakta/${REF_SAMPLE}.fna"
REF_GFF="$ANNOT_IN/$REF_GROUP/$REF_SAMPLE/Bakta/${REF_SAMPLE}.gff3"
REF_INFERNAL="$ANNOT_IN/$REF_GROUP/$REF_SAMPLE/Infernal/${REF_SAMPLE}.ncRNA_filtered.tsv"
REF_FAA="$ANNOT_IN/$REF_GROUP/$REF_SAMPLE/Bakta/${REF_SAMPLE}.faa"
REF_FFN="$ANNOT_IN/$REF_GROUP/$REF_SAMPLE/Bakta/${REF_SAMPLE}.ffn"

MUT_R1="$CLEAN_IN/$MUT_GROUP/$MUT_SAMPLE/${MUT_SAMPLE}_R1.fastq.gz"
MUT_R2="$CLEAN_IN/$MUT_GROUP/$MUT_SAMPLE/${MUT_SAMPLE}_R2.fastq.gz"
MUT_NANO="$CLEAN_IN/$MUT_GROUP/$MUT_SAMPLE/${MUT_SAMPLE}_nanopore.fastq.gz"

for file in "$REF_FASTA" "$REF_GFF" "$MUT_R1" "$MUT_R2" "$MUT_NANO"; do
    [[ -f "$file" ]] || die "Thiếu file đầu vào: $file"
done

# ---------------------------------------------------------
# KHỞI TẠO THƯ MỤC & CHECKPOINT
# ---------------------------------------------------------
PAIR_OUT="$VAR_OUT/${REF_SAMPLE}_vs_${MUT_SAMPLE}"
CUSTOM_DB_DIR="$PAIR_OUT/snpeff_db_${REF_SAMPLE}"
MAP_DIR="$PAIR_OUT/mapping"
VCF_DIR="$PAIR_OUT/vcf"
RUN_LOG="$PAIR_OUT/pairwise_varcall.log"
CKPT_DONE="$CHECKPOINT_DIR/${REF_SAMPLE}_vs_${MUT_SAMPLE}.ok"

if [[ "$DRY_RUN" == "true" ]]; then
    log "[DRY-RUN] mkdir -p $PAIR_OUT $CUSTOM_DB_DIR $MAP_DIR $VCF_DIR $CHECKPOINT_DIR"
else
    mkdir -p "$PAIR_OUT" "$CUSTOM_DB_DIR" "$MAP_DIR" "$VCF_DIR" "$CHECKPOINT_DIR"
    if [[ -f "$CKPT_DONE" ]]; then
        log "[RESUME] Đã hoàn tất phân tích cặp này. Bỏ qua."
        exit 0
    fi
    echo "Bắt đầu Pairwise Variant Calling lúc $(date)" > "$RUN_LOG"
fi

START_TIME="$(date +%s)"
FAIL=0

# ---------------------------------------------------------
# BƯỚC 1: XÂY DỰNG SNPEFF DB & BWA INDEX CHO MẪU NATIVE
# ---------------------------------------------------------
if [[ ! -f "$CUSTOM_DB_DIR/snpEff.config" || ! -f "$CUSTOM_DB_DIR/data/${REF_SAMPLE}/snpEffectPredictor.bin" ]]; then
    log "  -> [REFERENCE] Build Custom SnpEff DB tích hợp ncRNA cho $REF_SAMPLE..."
    if [[ "$DRY_RUN" == "true" ]]; then
        log "[DRY-RUN] [Python Script] Xử lý GFF, tích hợp ncRNA và copy sequences/cds/protein."
        log "[DRY-RUN] [conda activate $ENV_SNPEFF] snpEff -Xmx50g -Djava.io.tmpdir=$PAIR_OUT build -c $CUSTOM_DB_DIR/snpEff.config -gff3 -v ${REF_SAMPLE}"
    else
        mkdir -p "$CUSTOM_DB_DIR/data/${REF_SAMPLE}"
        
        python3 - "$REF_GFF" "$REF_INFERNAL" "$CUSTOM_DB_DIR/data/${REF_SAMPLE}/genes.gff" <<'PY'
import sys, csv
gff_in, infernal_in, gff_out = sys.argv[1:4]
with open(gff_in, 'r') as fin, open(gff_out, 'w') as fout:
    fout.write(fin.read())
    try:
        with open(infernal_in, 'r') as finf:
            for i, row in enumerate(csv.DictReader(finf, delimiter='\t')):
                contig = row.get('query_name')
                if not contig: continue
                try: sf, st = int(row['seq_from']), int(row['seq_to'])
                except: continue
                start, end = min(sf, st), max(sf, st)
                strand = row.get('strand', '.')
                name = row.get('target_name', f'ncRNA_{i}')
                gid, tid, eid = f"gene_ncRNA_{i}", f"transcript_ncRNA_{i}", f"exon_ncRNA_{i}"
                
                fout.write(f"{contig}\tInfernal\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid};Name={name};biotype=ncRNA\n")
                fout.write(f"{contig}\tInfernal\ttranscript\t{start}\t{end}\t.\t{strand}\t.\tID={tid};Parent={gid};Name={name};biotype=ncRNA\n")
                fout.write(f"{contig}\tInfernal\texon\t{start}\t{end}\t.\t{strand}\t.\tID={eid};Parent={tid};biotype=ncRNA\n")
    except FileNotFoundError: pass
PY

        cp "$REF_FASTA" "$CUSTOM_DB_DIR/data/${REF_SAMPLE}/sequences.fa"

        # --- FIX CDS HEADERS ---
        if [[ -f "$REF_FFN" ]]; then
            sed 's/^>/>TRANSCRIPT_/' "$REF_FFN" > "$CUSTOM_DB_DIR/data/${REF_SAMPLE}/cds.fa"
            [[ -s "$CUSTOM_DB_DIR/data/${REF_SAMPLE}/cds.fa" ]] || die "cds.fa rỗng sau xử lý"
        fi

        # --- FIX PROTEIN HEADERS ---
        if [[ -f "$REF_FAA" ]]; then
            sed 's/^>/>TRANSCRIPT_/' "$REF_FAA" > "$CUSTOM_DB_DIR/data/${REF_SAMPLE}/protein.fa"
            [[ -s "$CUSTOM_DB_DIR/data/${REF_SAMPLE}/protein.fa" ]] || die "protein.fa rỗng sau xử lý"
        fi

        echo "${REF_SAMPLE}.genome : ${REF_SAMPLE}" > "$CUSTOM_DB_DIR/snpEff.config"
        # Cấp 50GB RAM & Định tuyến thư mục tmp để chống OOM/Full Disk
        run_tool "$ENV_SNPEFF" "$RUN_LOG" bash -c "snpEff -Xmx50g -Djava.io.tmpdir='$PAIR_OUT' build -c '$CUSTOM_DB_DIR/snpEff.config' -gff3 -v '${REF_SAMPLE}'" || FAIL=1
    fi
fi

if [[ ! -f "${REF_FASTA}.bwt" ]]; then
    log "  -> [REFERENCE] Tạo BWA Index..."
    run_tool "$ENV_MAPPING" "$RUN_LOG" bwa index "$REF_FASTA" || FAIL=1
fi

# ---------------------------------------------------------
# BƯỚC 2: HYBRID MAPPING & JOINT CALLING CHO MẪU MUTANT
# ---------------------------------------------------------
BAM_ILL="$MAP_DIR/${MUT_SAMPLE}_illumina.bam"
BAM_ONT="$MAP_DIR/${MUT_SAMPLE}_nanopore.bam"
RAW_VCF="$VCF_DIR/${MUT_SAMPLE}_raw.vcf"
NORM_VCF="$VCF_DIR/${MUT_SAMPLE}_norm.vcf"
ANN_VCF="$VCF_DIR/${MUT_SAMPLE}_annotated.vcf"

if [[ "$FAIL" -eq 0 ]]; then
    log "  -> Mapping Illumina (Mutant -> Native)..."
    run_tool "$ENV_MAPPING" "$RUN_LOG" bash -c "set -o pipefail; bwa mem -t $THREADS '$REF_FASTA' '$MUT_R1' '$MUT_R2' | samtools sort -@ $THREADS -o '$BAM_ILL' -" || FAIL=1

    log "  -> Mapping Nanopore (Mutant -> Native)..."
    run_tool "$ENV_MAPPING" "$RUN_LOG" bash -c "set -o pipefail; minimap2 -ax map-ont -t $THREADS '$REF_FASTA' '$MUT_NANO' | samtools sort -@ $THREADS -o '$BAM_ONT' -" || FAIL=1
fi

if [[ "$FAIL" -eq 0 ]]; then
    log "  -> Indexing BAM files..."
    run_tool "$ENV_MAPPING" "$RUN_LOG" samtools index -@ "$THREADS" "$BAM_ILL" || FAIL=1
    run_tool "$ENV_MAPPING" "$RUN_LOG" samtools index -@ "$THREADS" "$BAM_ONT" || FAIL=1
    
    log "  -> Gọi biến thể chung (FreeBayes)..."
    run_tool "$ENV_CALLING" "$RUN_LOG" bash -c "freebayes -f '$REF_FASTA' -p 1 -F 0.3 -C 5 -b '$BAM_ILL' -b '$BAM_ONT' > '$RAW_VCF'" || FAIL=1
    
    log "  -> Chuẩn hóa VCF (bcftools norm)..."
    run_tool "$ENV_CALLING" "$RUN_LOG" bash -c "set -o pipefail; bcftools filter -i 'QUAL>20' '$RAW_VCF' | bcftools norm -f '$REF_FASTA' -m -any -O v -o '$NORM_VCF'" || FAIL=1

    log "  -> Annotate bằng SnpEff..."
    # Thêm cờ -noStats chặn tạo file rác HTML/TXT. Định tuyến tmpdir chống sập RAM/Disk
    run_tool "$ENV_SNPEFF" "$RUN_LOG" bash -c "snpEff -Xmx50g -Djava.io.tmpdir='$PAIR_OUT' ann -noStats -c '$CUSTOM_DB_DIR/snpEff.config' '${REF_SAMPLE}' '$NORM_VCF' > '$ANN_VCF'" || FAIL=1
    
    if [[ "$FAIL" -eq 0 ]]; then
        log "  -> Nén và tạo index (bgzip/tabix)..."
        run_tool "$ENV_CALLING" "$RUN_LOG" bash -c "bgzip -c '$ANN_VCF' > '${ANN_VCF}.gz'" || FAIL=1
        run_tool "$ENV_CALLING" "$RUN_LOG" tabix -p vcf "${ANN_VCF}.gz" || FAIL=1
        
        log "  -> Thống kê VCF (bcftools stats)..."
        run_tool "$ENV_CALLING" "$RUN_LOG" bash -c "bcftools stats '${ANN_VCF}.gz' > '${VCF_DIR}/${MUT_SAMPLE}_bcftools.stats'" || FAIL=1
    fi
fi

# ---------------------------------------------------------
# BƯỚC 3: TRÍCH XUẤT THỐNG KÊ
# ---------------------------------------------------------
ILL_MAP="-"
ONT_MAP="-"
TOTAL_SNPS="-"
TOTAL_INDELS="-"
COPA_MUTS="-"

if [[ "$FAIL" -eq 0 ]]; then
    log "  -> Trích xuất thống kê..."
    if [[ "$DRY_RUN" == "true" ]]; then
        ILL_MAP="99.0"
        ONT_MAP="95.0"
        TOTAL_SNPS="100"
        TOTAL_INDELS="10"
        COPA_MUTS="0"
    else
        ILL_MAP=$(conda run -n "$ENV_MAPPING" samtools flagstat "$BAM_ILL" 2>/dev/null | awk -F'[(%]' '/mapped \(/ {print $2}' || echo "-")
        ONT_MAP=$(conda run -n "$ENV_MAPPING" samtools flagstat "$BAM_ONT" 2>/dev/null | awk -F'[(%]' '/mapped \(/ {print $2}' || echo "-")
        
        if [[ -f "${VCF_DIR}/${MUT_SAMPLE}_bcftools.stats" ]]; then
            TOTAL_SNPS=$(grep '^SN' "${VCF_DIR}/${MUT_SAMPLE}_bcftools.stats" | grep 'number of SNPs:' | awk '{print $NF}' || echo "0")
            TOTAL_INDELS=$(grep '^SN' "${VCF_DIR}/${MUT_SAMPLE}_bcftools.stats" | grep 'number of indels:' | awk '{print $NF}' || echo "0")
        else
            log "  [WARNING] Không tìm thấy file bcftools.stats"
        fi
        
        COPA_MUTS=$(conda run -n "$ENV_CALLING" bcftools query -f '%INFO/ANN\n' "${ANN_VCF}.gz" 2>/dev/null | awk -F'|' '{ for(i=1; i<=NF; i+=15) if(toupper($i) ~ /COPA/) {print; next} }' | wc -l || echo "0")
    fi
fi

END_TIME="$(date +%s)"
ELAPSED=$((END_TIME - START_TIME))

if [[ "$FAIL" -eq 0 ]]; then
    append_summary "${REF_SAMPLE}\t${MUT_SAMPLE}\tSUCCESS\t${ILL_MAP}%\t${ONT_MAP}%\t${TOTAL_SNPS}\t${TOTAL_INDELS}\t${COPA_MUTS}\t${ELAPSED}s"
    if [[ "$DRY_RUN" == "false" ]]; then
        touch "$CKPT_DONE"
    fi
    log "[SUCCESS] So sánh $REF_SAMPLE và $MUT_SAMPLE hoàn tất (${ELAPSED}s)."
else
    append_summary "${REF_SAMPLE}\t${MUT_SAMPLE}\tFAILED\t-\t-\t-\t-\t-\t${ELAPSED}s"
    log "[ERROR] Phân tích thất bại. Xem chi tiết tại $RUN_LOG và _varcall_pair_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
fi