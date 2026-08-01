#!/bin/bash
#SBATCH --job-name=phylogenomics
#SBATCH --output=_phylo_%A_%a.out
#SBATCH --error=_phylo_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --array=1-1%3
# ^ CẬP NHẬT lại "1-N%<song_song>" cho khớp số CLADE thật sau khi chạy --prepare-clades (xem bên dưới).

# ==============================================================================
# Script: 10_phylogenomics_pangenome.sh (Module 7 - PHYLOGENOMICS & PANGENOME)
#
# PURPOSE:
#   Nhóm các mẫu theo Sequence Type (ST, lấy từ Kleborate ở Module 4), sau đó với
#   từng CLADE-ST (VD: ST11, ST15...) chạy pipeline:
#     Panaroo (core pangenome, alignment core gene)
#       -> Gubbins (lọc bỏ vùng nhiễu do tái tổ hợp / recombination)
#       -> IQ-TREE (dựng cây phát sinh loài ML, ultrafast bootstrap)
#       -> Comparative Genomics: tích hợp dữ liệu PCN / mất porin ompK35-36 / IS
#          chèn ép từ Module 5 (09_deep_dive_mechanisms) thành các file chú giải
#          iTOL (DATASET_BINARY, DATASET_SIMPLEBAR) gắn lên các nhánh của cây.
#
# INPUT (từ các module trước):
#   - $ROOT_DIR/07_annotation/<GROUP>/<SAMPLE>/Bakta/<SAMPLE>.gff3   (Module 3, có kèm ##FASTA)
#   - $ROOT_DIR/08_functional_screening/_combined/st_for_phylogenomics.tsv (Module 4, Kleborate ST)
#   - $ROOT_DIR/09_deep_dive_mechanisms/deep_dive_summary.tsv        (Module 5, PCN/Porin/IS)
#
# CÁCH DÙNG (2 bước):
#   Bước 1 - Chạy trên login node, KHÔNG cần sbatch (chỉ đọc TSV, không tốn tài nguyên):
#       bash 10_phylogenomics_pangenome.sh --prepare-clades
#     -> Nhóm mẫu theo Kleborate_ST, ghi "$ROOT_DIR/11_phylogenomics/clade_manifest.tsv"
#        và in ra số lượng CLADE hợp lệ N (đủ --min-clade-size mẫu, mặc định 3).
#
#   Bước 2 - Cập nhật dòng "#SBATCH --array=1-N%<song_song>" ở trên cho khớp N, rồi:
#       sbatch 10_phylogenomics_pangenome.sh
#     Hoặc chạy tuần tự không qua SLURM (VD debug 1 clade nhỏ):
#       bash 10_phylogenomics_pangenome.sh
#
# THAM SỐ:
#   --prepare-clades        Chỉ build clade_manifest.tsv rồi thoát, không chạy pipeline.
#   --dry-run               Giả lập, in lệnh sẽ chạy nhưng không thực thi/không ghi file.
#   --force                 Bỏ qua checkpoint, chạy lại từ đầu cho các clade đã xong.
#   --min-clade-size N      Số mẫu tối thiểu/ST để được coi là 1 clade hợp lệ (default: 3).
#
# LƯU Ý PHƯƠNG PHÁP:
#   Gubbins vốn được thiết kế cho whole-genome alignment (tọa độ genome liên tục).
#   Áp dụng lên core-gene alignment nối lại từ Panaroo (theo đúng sơ đồ pipeline) là
#   cách tiếp cận phổ biến để lọc SNP nhiễu do tái tổ hợp trước khi dựng cây, nhưng
#   tọa độ trong "*.recombination_predictions.gff" tương ứng vị trí trên alignment
#   đã nối, KHÔNG phải tọa độ trên genome gốc - cần diễn giải cẩn thận nếu muốn map
#   ngược vùng tái tổ hợp lên từng gen cụ thể.
# ==============================================================================

set -Eeuo pipefail
trap 'echo "[FATAL] Script bị ngắt đột ngột tại dòng $LINENO!" >&2; exit 1' ERR INT TERM

# ---------------------------------------------------------
# 0. CLI
# ---------------------------------------------------------
DRY_RUN=0
RESUME=1
PREPARE_ONLY=0
MIN_CLADE_SIZE=3

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)          DRY_RUN=1; shift ;;
        --force)            RESUME=0; shift ;;
        --prepare-clades)   PREPARE_ONLY=1; shift ;;
        --min-clade-size)   MIN_CLADE_SIZE="${2:?Thiếu giá trị cho --min-clade-size}"; shift 2 ;;
        *) echo "[WARNING] Tham số không xác định: $1 (bỏ qua)"; shift ;;
    esac
done

# ---------------------------------------------------------
# 1. CONFIG & AUTO-DETECT
# ---------------------------------------------------------
THREADS="${SLURM_CPUS_PER_TASK:-16}"

ROOT_DIR="${ROOT_DIR:-/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID}"
ANNOT_IN="${ANNOT_IN:-$ROOT_DIR/07_annotation}"
FUNC_OUT="${FUNC_OUT:-$ROOT_DIR/08_functional_screening}"
ST_FILE="${ST_FILE:-$FUNC_OUT/_combined/st_for_phylogenomics.tsv}"
DEEPDIVE_SUMMARY="${DEEPDIVE_SUMMARY:-$ROOT_DIR/09_deep_dive_mechanisms/deep_dive_summary.tsv}"

PHYLO_OUT="${PHYLO_OUT:-$ROOT_DIR/11_phylogenomics}"
CHECKPOINT_DIR="$PHYLO_OUT/.checkpoints"
CLADE_MANIFEST="$PHYLO_OUT/clade_manifest.tsv"
SKIPPED_FILE="$PHYLO_OUT/unclustered_samples.tsv"
SUMMARY="$PHYLO_OUT/phylogenomics_summary.tsv"

ENV_PANAROO="${ENV_PANAROO:-phylo_env}"
ENV_GUBBINS="${ENV_GUBBINS:-phylo_env}"
ENV_IQTREE="${ENV_IQTREE:-phylo_env}"

GUBBINS_BIN="${GUBBINS_BIN:-run_gubbins.py}"

# Tự động nhận diện lệnh iqtree hoặc iqtree2 (Tránh lỗi Command not found)
if [[ -z "${IQTREE_BIN:-}" ]]; then
    if conda run -n "$ENV_IQTREE" command -v iqtree2 >/dev/null 2>&1; then
        IQTREE_BIN="iqtree2"
    elif conda run -n "$ENV_IQTREE" command -v iqtree >/dev/null 2>&1; then
        IQTREE_BIN="iqtree"
    else
        IQTREE_BIN="iqtree2" # Giá trị mặc định để hàm check_program báo lỗi chuẩn
    fi
fi

PANAROO_CORE_THRESHOLD="${PANAROO_CORE_THRESHOLD:-0.95}"
IQTREE_MODEL="${IQTREE_MODEL:-MFP+ASC}"
IQTREE_BB="${IQTREE_BB:-1000}"
IQTREE_ALRT="${IQTREE_ALRT:-1000}"

# ---------------------------------------------------------
# 2. HELPERS
# ---------------------------------------------------------
ts() { date '+%Y-%m-%d %H:%M:%S%z'; }
log() { printf '[%s] %s\n' "$(ts)" "$*" >&2; }
die() { log "FATAL: $*"; exit 1; }

check_program() {
    local env_name="$1" prog="$2"
    if ! conda run -n "$env_name" command -v "$prog" >/dev/null 2>&1; then
        die "Thiếu dependency: $prog trong môi trường $env_name"
    fi
}

run_tool() {
    local env_name="$1" log_file="$2"
    shift 2
    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf '    [DRY-RUN] (%s) ' "$env_name"
        printf '%q ' "$@"
        printf ' >> %q\n' "$log_file"
        return 0
    fi
    (
        set -e
        echo "[$(ts)] COMMAND: $*"
        eval "$(conda shell.bash hook)"
        conda activate "$env_name" || { echo "[ERROR] Lỗi activate env: $env_name" >&2; exit 1; }
        "$@"
    ) >> "$log_file" 2>&1
}

append_summary() {
    local line="$1"
    [[ "$DRY_RUN" -eq 1 ]] && return 0
    (
        flock -x 200
        [[ -f "$SUMMARY" ]] || echo -e "CladeID\tST\tN_Samples\tStatus\tCore_Genes\tCore_Aln_Length_bp\tSNPs_After_Gubbins\tRecomb_Blocks\tTreefile\tTime\tMessage" > "$SUMMARY"
        echo -e "$line" >> "$SUMMARY"
    ) 200>"${SUMMARY}.lock"
}

# ---------------------------------------------------------
# 3. CHUẨN BỊ CLADE MANIFEST
# ---------------------------------------------------------
if [[ "$PREPARE_ONLY" -eq 1 ]]; then
    echo -e "\n=== CHUẨN BỊ CLADE MANIFEST (nhóm mẫu theo Kleborate ST - Module 4) ==="
    [[ -f "$ST_FILE" ]] || die "Không tìm thấy $ST_FILE. Cần chạy xong Module 4 trước: bash 07_functional_screening.sh --aggregate-only"
    mkdir -p "$PHYLO_OUT"

    N_CLADES=$(python3 - "$ST_FILE" "$ANNOT_IN" "$CLADE_MANIFEST" "$SKIPPED_FILE" "$MIN_CLADE_SIZE" <<'PY'
import csv, sys
from pathlib import Path
from collections import defaultdict

st_file, annot_in, out_manifest, out_skipped, min_size = sys.argv[1:6]
min_size = int(min_size)
annot_in = Path(annot_in)

groups = defaultdict(list)   # ST -> [(GROUP, SAMPLE), ...]
skipped = []                 # (Sample, Group, ST, Reason)

with open(st_file, encoding="utf-8-sig", newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        sample = (row.get("Sample") or "").strip()
        group = (row.get("Group") or "").strip()
        st = (row.get("Kleborate_ST") or "").strip()
        if not sample or not group:
            continue
        if not st or st in ("-", "NA", "nan", "0"):
            skipped.append((sample, group, st or "-", "missing_ST"))
            continue
        gff = annot_in / group / sample / "Bakta" / f"{sample}.gff3"
        if not gff.exists() or gff.stat().st_size == 0:
            skipped.append((sample, group, st, "missing_or_empty_gff3"))
            continue
        groups[st].append((group, sample))

clades = []
for st in sorted(groups.keys(), key=lambda x: (len(x), x)):
    members = sorted(set(groups[st]))
    if len(members) < min_size:
        for group, sample in members:
            skipped.append((sample, group, st, f"clade_too_small(<{min_size})"))
        continue
    clade_id = f"ST{st}" if st.isdigit() else f"Clade_{st}"
    clades.append((clade_id, st, members))

with open(out_manifest, "w", newline="") as oh:
    oh.write("CladeID\tST\tN_Samples\tSamples\n")
    for clade_id, st, members in clades:
        sample_str = ",".join(f"{g}:{s}" for g, s in members)
        oh.write(f"{clade_id}\t{st}\t{len(members)}\t{sample_str}\n")

with open(out_skipped, "w", newline="") as oh:
    oh.write("Sample\tGroup\tST\tReason\n")
    for sample, group, st, reason in skipped:
        oh.write(f"{sample}\t{group}\t{st}\t{reason}\n")

print(len(clades))
PY
)
    echo "  -> Đã ghi clade manifest : $CLADE_MANIFEST"
    echo "  -> Mẫu bị loại/lẻ ST      : $SKIPPED_FILE"
    echo "  -> Tổng số CLADE hợp lệ (>= $MIN_CLADE_SIZE mẫu/ST): $N_CLADES"
    echo ""
    if [[ "$N_CLADES" -gt 0 ]]; then
        echo "  ==> Cập nhật '#SBATCH --array=1-${N_CLADES}%<song_song>' ở đầu script rồi chạy: sbatch $0"
    else
        echo "  [WARNING] Không có clade nào đủ điều kiện. Kiểm tra $SKIPPED_FILE hoặc giảm --min-clade-size."
    fi
    exit 0
fi

# ---------------------------------------------------------
# 4. KHỞI TẠO MÔI TRƯỜNG & CHỐNG RACE-CONDITION
# ---------------------------------------------------------
echo -e "\n=== KIỂM TRA MÔI TRƯỜNG CONDA ==="
eval "$(conda shell.bash hook)"

if [[ "$DRY_RUN" -eq 0 ]]; then
    for env in "$ENV_PANAROO" "$ENV_GUBBINS" "$ENV_IQTREE"; do
        if ! conda env list | awk -v e="$env" '$1==e{f=1} END{exit !f}'; then
            die "Conda environment '${env}' không tồn tại!"
        fi
    done
    check_program "$ENV_PANAROO" panaroo
    check_program "$ENV_GUBBINS" "$GUBBINS_BIN"
    check_program "$ENV_IQTREE" "$IQTREE_BIN"

    mkdir -p "$PHYLO_OUT" "$CHECKPOINT_DIR"
    (
        flock -x 200
        [[ -f "$SUMMARY" ]] || echo -e "CladeID\tST\tN_Samples\tStatus\tCore_Genes\tCore_Aln_Length_bp\tSNPs_After_Gubbins\tRecomb_Blocks\tTreefile\tTime\tMessage" > "$SUMMARY"
    ) 200>"${SUMMARY}.lock"

    echo "Phần mềm sử dụng:"
    ( conda activate "$ENV_PANAROO"; panaroo --version 2>&1 || true )
    ( conda activate "$ENV_GUBBINS"; "$GUBBINS_BIN" --version 2>&1 || true )
    ( conda activate "$ENV_IQTREE"; "$IQTREE_BIN" --version 2>&1 | head -n 1 || true )
else
    mkdir -p "$PHYLO_OUT" 2>/dev/null || true
fi

[[ -f "$CLADE_MANIFEST" ]] || die "Không tìm thấy $CLADE_MANIFEST. Hãy chạy trước: bash $0 --prepare-clades"
[[ -f "$DEEPDIVE_SUMMARY" ]] || log "[WARNING] Không thấy $DEEPDIVE_SUMMARY - bước Comparative Genomics (Module 5) sẽ bị bỏ qua/để trống."

# ---------------------------------------------------------
# 5. XỬ LÝ 1 CLADE
# ---------------------------------------------------------
run_clade() {
    local clade_line="$1"
    local CLADE_ID ST N_SAMPLES SAMPLES_STR
    IFS=$'\t' read -r CLADE_ID ST N_SAMPLES SAMPLES_STR <<< "$clade_line"

    log "=================================================="
    log "ĐANG XỬ LÝ CLADE: $CLADE_ID (ST=$ST, N=$N_SAMPLES mẫu)"

    local CLADE_DIR="$PHYLO_OUT/$CLADE_ID"
    local PANAROO_OUT="$CLADE_DIR/Panaroo"
    local GUBBINS_DIR="$CLADE_DIR/Gubbins"
    local IQTREE_DIR="$CLADE_DIR/IQTREE"
    local COMPARATIVE_DIR="$CLADE_DIR/Comparative"
    local RUN_LOG="$CLADE_DIR/phylogenomics.log"

    local CKPT_PANAROO="$CHECKPOINT_DIR/${CLADE_ID}.panaroo.ok"
    local CKPT_GUBBINS="$CHECKPOINT_DIR/${CLADE_ID}.gubbins.ok"
    local CKPT_IQTREE="$CHECKPOINT_DIR/${CLADE_ID}.iqtree.ok"
    local CKPT_DONE="$CHECKPOINT_DIR/${CLADE_ID}.phylo.ok"

    if [[ "$RESUME" -eq 1 && "$DRY_RUN" -eq 0 && -f "$CKPT_DONE" ]]; then
        log "[RESUME] Clade $CLADE_ID đã hoàn tất trước đó -> Bỏ qua."
        return 0
    fi

    local -a GFF_ARR=()
    local -a _pairs=()
    local tok g s gff
    IFS=',' read -r -a _pairs <<< "$SAMPLES_STR"
    for tok in "${_pairs[@]}"; do
        g="${tok%%:*}"; s="${tok#*:}"
        gff="$ANNOT_IN/$g/$s/Bakta/${s}.gff3"
        if [[ ! -s "$gff" ]]; then
            log "  [ERROR] Thiếu GFF3 cho $s: $gff"
            append_summary "${CLADE_ID}\t${ST}\t${N_SAMPLES}\tFAILED\t-\t-\t-\t-\t-\t-\tMissing GFF3: $s"
            return 1
        fi
        GFF_ARR+=("$gff")
    done

    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$PANAROO_OUT" "$GUBBINS_DIR" "$IQTREE_DIR" "$COMPARATIVE_DIR"
        echo "Bắt đầu Phylogenomics cho $CLADE_ID (ST=$ST) lúc $(date)" > "$RUN_LOG"
        printf 'GFF3 input (%d files):\n' "${#GFF_ARR[@]}" >> "$RUN_LOG"
        printf '  %s\n' "${GFF_ARR[@]}" >> "$RUN_LOG"
    fi

    local START_TIME; START_TIME=$(date +%s)
    local FAIL=0

    # --- 5.1 PANAROO ---
    if [[ "$RESUME" -eq 1 && "$DRY_RUN" -eq 0 && -f "$CKPT_PANAROO" ]]; then
        log "  [INFO] Panaroo đã xong, bỏ qua."
    else
        log "  -> [1/3] Chạy Panaroo (core_threshold=$PANAROO_CORE_THRESHOLD)..."
        if ! run_tool "$ENV_PANAROO" "$RUN_LOG" panaroo -i "${GFF_ARR[@]}" -o "$PANAROO_OUT" \
                --clean-mode strict --remove-invalid-genes -a core --core_threshold "$PANAROO_CORE_THRESHOLD" -t "$THREADS"; then
            log "  [ERROR] Panaroo thất bại cho $CLADE_ID!"
            FAIL=1
        elif [[ "$DRY_RUN" -eq 0 ]]; then
            touch "$CKPT_PANAROO"
        fi
    fi

    local CORE_ALN=""
    if [[ "$DRY_RUN" -eq 0 && "$FAIL" -eq 0 ]]; then
        for cand in "$PANAROO_OUT/core_gene_alignment_filtered.aln" "$PANAROO_OUT/core_gene_alignment.aln"; do
            if [[ -s "$cand" ]]; then CORE_ALN="$cand"; break; fi
        done
        if [[ -z "$CORE_ALN" ]]; then
            log "  [ERROR] Không tìm thấy core_gene_alignment(.filtered).aln sau Panaroo."
            FAIL=1
        else
            local SEQ_COUNT
            SEQ_COUNT=$(grep -c "^>" "$CORE_ALN" || echo 0)
            if (( SEQ_COUNT < MIN_CLADE_SIZE )); then
                log "  [ERROR] Core alignment chỉ có $SEQ_COUNT sequences (Cần >= $MIN_CLADE_SIZE). Bỏ qua."
                FAIL=1
            fi
        fi
    fi

    # --- 5.2 GUBBINS ---
    local GUBBINS_PREFIX="$GUBBINS_DIR/${CLADE_ID}"
    if [[ "$RESUME" -eq 1 && "$DRY_RUN" -eq 0 && -f "$CKPT_GUBBINS" ]]; then
        log "  [INFO] Gubbins đã xong, bỏ qua."
    elif [[ "$FAIL" -eq 0 ]]; then
        log "  -> [2/3] Chạy Gubbins (loại tái tổ hợp)..."
        local GUBBINS_INPUT="${CORE_ALN:-$PANAROO_OUT/core_gene_alignment.aln}"
        if ! run_tool "$ENV_GUBBINS" "$RUN_LOG" "$GUBBINS_BIN" --threads "$THREADS" \
                --tree-builder fasttree --prefix "$GUBBINS_PREFIX" "$GUBBINS_INPUT"; then
            log "  [ERROR] Gubbins thất bại cho $CLADE_ID!"
            FAIL=1
        elif [[ "$DRY_RUN" -eq 0 ]]; then
            touch "$CKPT_GUBBINS"
        fi
    fi

    local SNP_ALN="${GUBBINS_PREFIX}.filtered_polymorphic_sites.fasta"
    if [[ "$DRY_RUN" -eq 0 && "$FAIL" -eq 0 && ! -s "$SNP_ALN" ]]; then
        log "  [ERROR] Không tìm thấy $SNP_ALN sau Gubbins."
        FAIL=1
    fi

    # --- 5.3 IQ-TREE ---
    local IQTREE_PREFIX="$IQTREE_DIR/${CLADE_ID}"
    if [[ "$RESUME" -eq 1 && "$DRY_RUN" -eq 0 && -f "$CKPT_IQTREE" ]]; then
        log "  [INFO] IQ-TREE đã xong, bỏ qua."
    elif [[ "$FAIL" -eq 0 ]]; then
        log "  -> [3/3] Chạy IQ-TREE (model=$IQTREE_MODEL, UFBoot=$IQTREE_BB, SH-aLRT=$IQTREE_ALRT)..."
        local -a IQTREE_ARGS=(-s "$SNP_ALN" -m "$IQTREE_MODEL" -bb "$IQTREE_BB" -alrt "$IQTREE_ALRT" \
                               -nt AUTO -ntmax "$THREADS" -pre "$IQTREE_PREFIX")
        [[ "$RESUME" -eq 0 ]] && IQTREE_ARGS+=(-redo)
        if ! run_tool "$ENV_IQTREE" "$RUN_LOG" "$IQTREE_BIN" "${IQTREE_ARGS[@]}"; then
            log "  [ERROR] IQ-TREE thất bại cho $CLADE_ID!"
            FAIL=1
        elif [[ "$DRY_RUN" -eq 0 ]]; then
            touch "$CKPT_IQTREE"
        fi
    fi

    local TREEFILE="${IQTREE_PREFIX}.treefile"
    if [[ "$DRY_RUN" -eq 0 && "$FAIL" -eq 0 && ! -s "$TREEFILE" ]]; then
        log "  [ERROR] Không tìm thấy $TREEFILE sau IQ-TREE."
        FAIL=1
    fi

    # --- 5.4 COMPARATIVE GENOMICS ---
    if [[ "$DRY_RUN" -eq 0 && "$FAIL" -eq 0 ]]; then
        log "  -> Tích hợp dữ liệu Module 5 (PCN/Porin/IS) lên cây ($CLADE_ID)..."
        if ! python3 - "$DEEPDIVE_SUMMARY" "$ST_FILE" "$CLADE_ID" "$SAMPLES_STR" "$COMPARATIVE_DIR" <<'PY' >> "$RUN_LOG" 2>&1
import csv, sys
from pathlib import Path

deepdive_tsv, st_tsv, clade_id, samples_str, out_dir = sys.argv[1:6]
out_dir = Path(out_dir)
out_dir.mkdir(parents=True, exist_ok=True)

clade_samples = []
for tok in samples_str.split(","):
    tok = tok.strip()
    if not tok or ":" not in tok:
        continue
    _, s = tok.split(":", 1)
    clade_samples.append(s)

def read_tsv_dict(path, key_col):
    out = {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return out
    with p.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            k = (row.get(key_col) or "").strip()
            if k:
                out[k] = row
    return out

deepdive = read_tsv_dict(deepdive_tsv, "Sample")
st_info = read_tsv_dict(st_tsv, "Sample")

def avg_pcn(pcn_field):
    if not pcn_field or pcn_field in ("-", "No_Plasmids"):
        return ""
    vals = []
    for tok in pcn_field.split(";"):
        tok = tok.strip()
        if not tok or ":" not in tok:
            continue
        v = tok.split(":")[-1].rstrip("Xx")
        try:
            vals.append(float(v))
        except ValueError:
            continue
    return f"{sum(vals)/len(vals):.2f}" if vals else ""

meta_rows = []
for sample in clade_samples:
    dd = deepdive.get(sample, {})
    sti = st_info.get(sample, {})
    meta_rows.append({
        "Sample": sample,
        "ST": sti.get("Kleborate_ST", "-"),
        "K_locus": sti.get("K_locus", "-"),
        "O_locus": sti.get("O_locus", "-"),
        "Chromosome_Depth": dd.get("Chromosome_Depth", "-"),
        "Avg_Plasmid_PCN": avg_pcn(dd.get("Plasmids_PCN", "")),
        "Plasmids_PCN_raw": dd.get("Plasmids_PCN", "-"),
        "ompK35_Status": dd.get("ompK35_Status", "-"),
        "ompK36_Status": dd.get("ompK36_Status", "-"),
        "IS_Interrupting_Porin": dd.get("IS_Interrupting_Porin", "-"),
    })

meta_cols = ["Sample", "ST", "K_locus", "O_locus", "Chromosome_Depth", "Avg_Plasmid_PCN",
             "Plasmids_PCN_raw", "ompK35_Status", "ompK36_Status", "IS_Interrupting_Porin"]
with (out_dir / f"{clade_id}_comparative_metadata.tsv").open("w", newline="") as oh:
    w = csv.DictWriter(oh, fieldnames=meta_cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in meta_rows:
        w.writerow(r)

with (out_dir / f"{clade_id}_itol_porin_loss.txt").open("w") as oh:
    oh.write("DATASET_BINARY\nSEPARATOR TAB\n")
    oh.write(f"DATASET_LABEL\t{clade_id}_Porin_Loss\n")
    oh.write("COLOR\t#e41a1c\n")
    oh.write("FIELD_SHAPES\t2\t2\n")
    oh.write("FIELD_LABELS\tompK35_loss\tompK36_loss\n")
    oh.write("FIELD_COLORS\t#e41a1c\t#377eb8\n")
    oh.write("DATA\n")
    for r in meta_rows:
        k35 = 1 if r["ompK35_Status"] not in ("Intact", "-") else 0
        k36 = 1 if r["ompK36_Status"] not in ("Intact", "-") else 0
        oh.write(f"{r['Sample']}\t{k35}\t{k36}\n")

with (out_dir / f"{clade_id}_itol_is_interruption.txt").open("w") as oh:
    oh.write("DATASET_BINARY\nSEPARATOR TAB\n")
    oh.write(f"DATASET_LABEL\t{clade_id}_IS_Interruption\n")
    oh.write("COLOR\t#984ea3\n")
    oh.write("FIELD_SHAPES\t2\n")
    oh.write("FIELD_LABELS\tIS_in_porin\n")
    oh.write("FIELD_COLORS\t#984ea3\n")
    oh.write("DATA\n")
    for r in meta_rows:
        v = r["IS_Interrupting_Porin"]
        present = 1 if v and v not in ("-", "None") else 0
        oh.write(f"{r['Sample']}\t{present}\n")

with (out_dir / f"{clade_id}_itol_plasmid_pcn.txt").open("w") as oh:
    oh.write("DATASET_SIMPLEBAR\nSEPARATOR TAB\n")
    oh.write(f"DATASET_LABEL\t{clade_id}_Avg_Plasmid_PCN\n")
    oh.write("COLOR\t#4daf4a\n")
    oh.write("DATA\n")
    for r in meta_rows:
        oh.write(f"{r['Sample']}\t{r['Avg_Plasmid_PCN'] or '0'}\n")

print("OK")
PY
        then
            log "  [WARNING] Bước Comparative Genomics (Module 5 -> iTOL) lỗi, xem $RUN_LOG. Cây/alignment vẫn giữ nguyên kết quả."
        fi
    fi

    # --- 5.5 TRÍCH XUẤT THỐNG KÊ & GHI SUMMARY ---
    local END_TIME; END_TIME=$(date +%s)
    local ELAPSED=$((END_TIME - START_TIME))

    local CORE_GENES="-" ALN_LEN="-" SNPS_GUBBINS="-" RECOMB_BLOCKS="-"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        if [[ -f "$PANAROO_OUT/summary_statistics.txt" ]]; then
            CORE_GENES=$(awk -F'\t' 'tolower($1) ~ /^core genes/ {print $2; exit}' "$PANAROO_OUT/summary_statistics.txt")
            [[ -z "$CORE_GENES" ]] && CORE_GENES="-"
        fi
        if [[ -n "$CORE_ALN" && -s "$CORE_ALN" ]]; then
            ALN_LEN=$(awk '/^>/{if(seq){print length(seq); exit} next} {seq=seq $0} END{if(seq) print length(seq)}' "$CORE_ALN")
            [[ -z "$ALN_LEN" ]] && ALN_LEN="-"
        fi
        if [[ -s "$SNP_ALN" ]]; then
            SNPS_GUBBINS=$(awk '/^>/{if(seq){print length(seq); exit} next} {seq=seq $0} END{if(seq) print length(seq)}' "$SNP_ALN")
            [[ -z "$SNPS_GUBBINS" ]] && SNPS_GUBBINS="-"
        fi
        if [[ -f "${GUBBINS_PREFIX}.recombination_predictions.gff" ]]; then
            RECOMB_BLOCKS=$(grep -vc '^#' "${GUBBINS_PREFIX}.recombination_predictions.gff" || echo 0)
        fi
    fi

    if [[ "$FAIL" -eq 0 ]]; then
        append_summary "${CLADE_ID}\t${ST}\t${N_SAMPLES}\tSUCCESS\t${CORE_GENES}\t${ALN_LEN}\t${SNPS_GUBBINS}\t${RECOMB_BLOCKS}\t${TREEFILE}\t${ELAPSED}s\tOK"
        [[ "$DRY_RUN" -eq 0 ]] && touch "$CKPT_DONE"
        log "  [SUCCESS] Clade $CLADE_ID hoàn tất (${ELAPSED}s). Core genes=$CORE_GENES | SNPs sau Gubbins=$SNPS_GUBBINS | Recomb blocks=$RECOMB_BLOCKS"
    else
        append_summary "${CLADE_ID}\t${ST}\t${N_SAMPLES}\tFAILED\t${CORE_GENES}\t${ALN_LEN}\t${SNPS_GUBBINS}\t${RECOMB_BLOCKS}\t-\t${ELAPSED}s\tXem $RUN_LOG"
        log "  [ERROR] Clade $CLADE_ID thất bại. Xem log: $RUN_LOG"
    fi
}

# ---------------------------------------------------------
# 6. CHỌN CLADE THEO SLURM ARRAY / CHẠY TUẦN TỰ
# ---------------------------------------------------------
TOTAL_CLADES=$(( $(wc -l < "$CLADE_MANIFEST") - 1 ))
(( TOTAL_CLADES > 0 )) || die "$CLADE_MANIFEST không có clade hợp lệ nào (đã chạy --prepare-clades chưa?)."

declare -a TARGET_LINE_NOS=()
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    if (( SLURM_ARRAY_TASK_ID > TOTAL_CLADES )); then
        log "Task ID $SLURM_ARRAY_TASK_ID vượt quá số lượng clade ($TOTAL_CLADES)."
        exit 0
    fi
    TARGET_LINE_NOS=("$((SLURM_ARRAY_TASK_ID + 1))")
    log "SLURM Array mode: task $SLURM_ARRAY_TASK_ID / $TOTAL_CLADES"
else
    for ((i = 2; i <= TOTAL_CLADES + 1; i++)); do TARGET_LINE_NOS+=("$i"); done
    log "Sequential mode: chạy tuần tự $TOTAL_CLADES clade"
fi

for line_no in "${TARGET_LINE_NOS[@]}"; do
    clade_line=$(sed -n "${line_no}p" "$CLADE_MANIFEST")
    [[ -z "$clade_line" ]] && continue
    run_clade "$clade_line"
done

echo -e "\n=================================================="
echo " ĐÃ CHẠY XONG PHYLOGENOMICS (MODULE 7) CHO TASK NÀY"
echo "=================================================="
date