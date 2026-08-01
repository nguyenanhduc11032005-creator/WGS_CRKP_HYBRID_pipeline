#!/bin/bash
#SBATCH --job-name=func_screen
#SBATCH --output=_func_%A_%a.out
#SBATCH --error=_func_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=1-19%5

# ==============================================================================
# Script: 07_functional_screening.sh (ULTIMATE PARALLEL VERSION - FINAL)
# Đã tương thích 100% với custom Kleborate wrapper (tạo folder thay vì file).
# ==============================================================================

set -Eeuo pipefail

trap 'kill $(jobs -p) 2>/dev/null || true; echo "[FATAL] Script bị ngắt đột ngột tại dòng $LINENO!"; exit 1' ERR INT TERM

# ---------------------------------------------------------
# 0. CLI
# ---------------------------------------------------------
DRY_RUN=0
RESUME=1
INCLUDE_FAILED_QC=0
AGGREGATE_ONLY=0

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        --force) RESUME=0 ;;
        --include-failed-qc) INCLUDE_FAILED_QC=1 ;;
        --aggregate-only) AGGREGATE_ONLY=1 ;;
    esac
done

# ---------------------------------------------------------
# 1. CONFIG
# ---------------------------------------------------------
THREADS="${SLURM_CPUS_PER_TASK:-32}"

ROOT_DIR="${ROOT_DIR:-/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID}"
ASM_IN="${ASM_IN:-$ROOT_DIR/05_assembly}"
ANNOT_IN="${ANNOT_IN:-$ROOT_DIR/07_annotation}"
ASM_QC_SUMMARY="${ASM_QC_SUMMARY:-$ROOT_DIR/06_assembly_qc/assembly_qc_summary.tsv}"

FUNC_OUT="${FUNC_OUT:-$ROOT_DIR/08_functional_screening}"
CHECKPOINT_DIR="$FUNC_OUT/.checkpoints"
SUMMARY="$FUNC_OUT/functional_screening_summary.tsv"
COMBINED_DIR="$FUNC_OUT/_combined"

ENV_KLEBORATE="${ENV_KLEBORATE:-kleborate_env}"
ENV_AMRFINDER="${ENV_AMRFINDER:-amrfinder_env}"
ENV_STARAMR="${ENV_STARAMR:-staramr_env}"
ENV_ABRICATE="${ENV_ABRICATE:-abricate_env}"
ENV_MOBSUITE="${ENV_MOBSUITE:-mobsuite_env}"

AMRFINDER_ORGANISM="${AMRFINDER_ORGANISM:-Klebsiella_pneumoniae}"
STARAMR_POINTFINDER_ORGANISM="${STARAMR_POINTFINDER_ORGANISM:-}"

VFDB_MINID="${VFDB_MINID:-80}"
VFDB_MINCOV="${VFDB_MINCOV:-60}"
PLASMIDFINDER_MINID="${PLASMIDFINDER_MINID:-80}"
PLASMIDFINDER_MINCOV="${PLASMIDFINDER_MINCOV:-60}"

# ---------------------------------------------------------
# 2. HELPERS (BASH & PYTHON)
# ---------------------------------------------------------
ts() { date '+%Y-%m-%d %H:%M:%S%z'; }
log() { printf '[%s] %s\n' "$(ts)" "$*" >&2; }
die() { log "FATAL: $*"; exit 1; }

run_tool_bg() {
    local env_name="$1"
    local tool_log="$2"
    shift 2

    mkdir -p "$(dirname "$tool_log")"

    if [[ "$DRY_RUN" -eq 1 ]]; then
        printf '    [DRY-RUN] (%s) %q > %q\n' "$env_name" "$*" "$tool_log"
        return 0
    fi

    (
        eval "$(conda shell.bash hook)"
        conda activate "$env_name"
        echo "[$(ts)] BẮT ĐẦU: $*" > "$tool_log"
        if "$@" >> "$tool_log" 2>&1; then
            echo "[$(ts)] KẾT THÚC (SUCCESS)" >> "$tool_log"
            exit 0
        else
            local rc=$?
            echo "[$(ts)] KẾT THÚC (FAILED - Exit Code: $rc)" >> "$tool_log"
            exit $rc
        fi
    ) &
}

# Python Metric Collector (Hỗ trợ chui vào thư mục Kleborate)
collect_metrics() {
    local kleb_tsv="$1"
    local amrf_tsv="$2"
    local staramr_dir="$3"
    local vfdb_tsv="$4"
    local plasmid_tsv="$5"
    local mob_dir="$6"

    python3 - "$kleb_tsv" "$amrf_tsv" "$staramr_dir" "$vfdb_tsv" "$plasmid_tsv" "$mob_dir" <<'PY'
import csv, re, sys
from pathlib import Path

kleb_tsv, amrf_tsv, staramr_dir, vfdb_tsv, plasmid_tsv, mob_dir = map(Path, sys.argv[1:7])

def clean(x):
    x = "" if x is None else str(x)
    return x.replace("\t", " ").replace("\r", " ").replace("\n", " ").strip() or "-"

def norm_key(x):
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())

def read_tsv(path):
    # Nếu path là folder (do Kleborate wrapper tạo ra), chui vào tìm file thật
    if path.is_dir():
        target = path / "klebsiella_pneumo_complex_output.txt"
        if target.exists() and target.stat().st_size > 0:
            path = target
        else:
            # Fallback nếu tên file đổi khác
            for f in path.glob("*output.txt"):
                if "hAMRonization" not in f.name:
                    path = f
                    break
            else:
                return []
                
    if not path.exists() or path.is_dir() or path.stat().st_size == 0:
        return []
        
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return []
        return [r for r in reader if any((v or "").strip() for v in r.values())]

def get(row, *aliases):
    if not row: return ""
    m = {norm_key(k): v for k, v in row.items() if k is not None}
    for a in aliases:
        an = norm_key(a)
        if an in m and str(m[an]).strip() not in ("", "-"):
            return str(m[an]).strip()
    for a in aliases:
        an = norm_key(a)
        for k, v in m.items():
            if an in k and str(v).strip() not in ("", "-"):
                return str(v).strip()
    return ""

def uniq_join(vals, max_items=60, max_chars=900):
    out, seen = [], set()
    for v in vals:
        v = clean(v)
        if v != "-" and v not in seen:
            seen.add(v)
            out.append(v)
    text = ";".join(out[:max_items])
    if len(out) > max_items: text += f";...(+{len(out)-max_items})"
    return text[:max_chars] + "..." if len(text) > max_chars else (text or "-")

# Kleborate
kleb_rows = read_tsv(kleb_tsv)
kr = kleb_rows[0] if kleb_rows else {}
species = get(kr, "species", "speciesmatch", "Klebsiellaspecies")
st = get(kr, "ST", "mlstst", "sequencetype")
k_locus = get(kr, "K_locus", "KL", "Klocus")
o_locus = get(kr, "O_locus", "OL", "Olocus")
vir_score = get(kr, "virulence_score", "virulencescore")
res_score = get(kr, "resistance_score", "resistancescore")

# AMRFinderPlus
amr_rows = read_tsv(amrf_tsv)
amr_hits = len(amr_rows)
amr_genes = [get(r, "Genesymbol", "Elementsymbol", "Gene", "Name", "Sequencename", "Allele") for r in amr_rows]

# StarAMR
staramr_hits = sum(len(read_tsv(staramr_dir / fn)) for fn in ["resfinder.tsv", "pointfinder.tsv"])

# VFDB via ABRicate
vfdb_rows = read_tsv(vfdb_tsv)
vfdb_hits = len(vfdb_rows)
vfdb_genes = [get(r, "GENE", "ACCESSION", "PRODUCT") for r in vfdb_rows]

# PlasmidFinder via ABRicate
plasmid_rows = read_tsv(plasmid_tsv)
plasmid_hits = len(plasmid_rows)
plasmid_reps = [get(r, "GENE", "ACCESSION", "PRODUCT") for r in plasmid_rows]

# MOB-suite
mob_path = mob_dir / "mobtyper_results.txt"
if not mob_path.exists():
    candidates = sorted(mob_dir.glob("*mobtyper*results*.txt"))
    mob_path = candidates[0] if candidates else mob_path

mob_rows = read_tsv(mob_path)
mob_plasmids = len(mob_rows)
mob_reps = []
for r in mob_rows:
    rep = get(r, "reptypes", "reptype", "replicontype", "replicon")
    if rep:
        mob_reps.extend(part.strip() for part in re.split(r"[;,]+", rep) if part.strip())

fields = [
    clean(species), clean(st), clean(k_locus), clean(o_locus), clean(vir_score), clean(res_score),
    str(amr_hits), uniq_join(amr_genes), str(staramr_hits), str(vfdb_hits), uniq_join(vfdb_genes),
    str(plasmid_hits), uniq_join(plasmid_reps), str(mob_plasmids), uniq_join(mob_reps),
]
print("\t".join(fields))
PY
}

update_summary() {
    local sample_key="$1"; shift
    [[ "$DRY_RUN" -eq 1 ]] && return 0
    (
        flock -x 200
        [[ -f "$SUMMARY" ]] || echo -e "Sample\tGroup\tStatus\tKleborate_species\tKleborate_ST\tK_locus\tO_locus\tVirulence_score\tResistance_score\tAMRFinder_hits\tAMRFinder_genes\tStarAMR_hits\tVFDB_hits\tVFDB_genes\tPlasmidFinder_hits\tPlasmidFinder_replicons\tMOB_plasmids\tMOB_replicons\tTime\tMessage" > "$SUMMARY"
        
        local tmp="$(mktemp "${SUMMARY}.tmp.XXXXXX")"
        awk -F'\t' -v s="$sample_key" 'NR == 1 || $1 != s' "$SUMMARY" > "$tmp"
        printf '%s\n' "$(IFS=$'\t'; echo "$*")" >> "$tmp"
        mv "$tmp" "$SUMMARY"
    ) 200>"${SUMMARY}.lock"
}

sample_has_qc_failed() {
    local sample="$1"
    [[ -f "$ASM_QC_SUMMARY" ]] || return 1
    local qc_status
    qc_status="$(awk -F'\t' -v s="$sample" 'NR > 1 && $1 == s {print $3; exit}' "$ASM_QC_SUMMARY" 2>/dev/null || true)"
    [[ "$qc_status" == "FAILED" ]]
}

# ---------------------------------------------------------
# 3. RUN PER SAMPLE
# ---------------------------------------------------------
run_sample() {
    local SPATH="$1"
    local GROUP="$(basename "$(dirname "$SPATH")")"
    local SAMPLE="$(basename "$SPATH")"
    local FASTA_IN="$SPATH/assembly.fasta"

    echo
    log "=================================================="
    log "ĐANG XỬ LÝ: $SAMPLE ; GROUP=$GROUP"

    if [[ "$INCLUDE_FAILED_QC" -eq 0 ]] && sample_has_qc_failed "$SAMPLE"; then
        log "[SKIP] $SAMPLE có Assembly QC = FAILED."
        update_summary "$SAMPLE" "$SAMPLE" "$GROUP" "SKIPPED" "-" "-" "-" "-" "-" "-" "0" "-" "0" "0" "-" "0" "-" "0" "-" "0s" "Assembly QC failed"
        return 0
    fi

    local SAMPLE_OUT="$FUNC_OUT/$GROUP/$SAMPLE"
    local KLEB_DIR="$SAMPLE_OUT/Kleborate"
    local AMRF_DIR="$SAMPLE_OUT/AMRFinderPlus"
    local STARAMR_DIR="$SAMPLE_OUT/StarAMR"
    local VFDB_DIR="$SAMPLE_OUT/VFDB"
    local PLASMID_DIR="$SAMPLE_OUT/PlasmidFinder"
    local MOB_DIR="$SAMPLE_OUT/MOBSuite"
    
    local MAIN_LOG="$SAMPLE_OUT/functional_screening.log"

    local KLEB_LOG="$SAMPLE_OUT/kleborate.log"
    local AMRF_LOG="$SAMPLE_OUT/amrfinder.log"
    local STAR_LOG="$SAMPLE_OUT/staramr.log"
    local VFDB_LOG="$SAMPLE_OUT/vfdb.log"
    local PLAS_LOG="$SAMPLE_OUT/plasmidfinder.log"
    local MOB_LOG="$SAMPLE_OUT/mobsuite.log"

    if [[ "$RESUME" -eq 0 && "$DRY_RUN" -eq 0 ]]; then
        rm -rf "$SAMPLE_OUT"/*
    fi

    if [[ "$DRY_RUN" -eq 0 ]]; then
        mkdir -p "$KLEB_DIR" "$AMRF_DIR" "$VFDB_DIR" "$PLASMID_DIR" "$MOB_DIR"
    fi

    local CKPT_DONE="$CHECKPOINT_DIR/${GROUP}__${SAMPLE}.functional.ok"
    if [[ "$RESUME" -eq 1 && "$DRY_RUN" -eq 0 && -f "$CKPT_DONE" ]]; then
        log "[RESUME] Đã hoàn tất cho $SAMPLE -> Bỏ qua."
        return 0
    fi

    local START_TIME="$(date +%s)"
    local FAIL=0
    local FAILED_TOOLS=()

    log "  -> Khởi chạy đồng thời các công cụ trên 32 cores (Background jobs)..."

    # 1. Kleborate (ĐÃ FIX: Dùng -o bình thường và để wrapper tự tạo thư mục)
    local KLEB_TSV="$KLEB_DIR/${SAMPLE}_kleborate_out"
    run_tool_bg "$ENV_KLEBORATE" "$KLEB_LOG" \
        kleborate -a "$FASTA_IN" -o "$KLEB_TSV" --preset kpsc
    PID_KLEB=$!

    # 2. AMRFinderPlus
    local AMRF_TSV="$AMRF_DIR/${SAMPLE}.amrfinder.tsv"
    local BAKTA_FAA="$ANNOT_IN/$GROUP/$SAMPLE/Bakta/${SAMPLE}.faa"
    local BAKTA_GFF="$ANNOT_IN/$GROUP/$SAMPLE/Bakta/${SAMPLE}.gff3"
    
    local AMRF_ARGS=(--plus --threads 8 -n "$FASTA_IN" -o "$AMRF_TSV")
    [[ -n "$AMRFINDER_ORGANISM" ]] && AMRF_ARGS+=(-O "$AMRFINDER_ORGANISM")
    [[ -s "$BAKTA_FAA" ]] && AMRF_ARGS+=(-p "$BAKTA_FAA")
    [[ -s "$BAKTA_GFF" ]] && AMRF_ARGS+=(-g "$BAKTA_GFF" --annotation_format bakta)
    
    run_tool_bg "$ENV_AMRFINDER" "$AMRF_LOG" amrfinder "${AMRF_ARGS[@]}"
    PID_AMRF=$!

    # 3. StarAMR
    [[ "$DRY_RUN" -eq 0 ]] && rm -rf "$STARAMR_DIR" 
    local STAR_ARGS=(search -o "$STARAMR_DIR")
    [[ -n "$STARAMR_POINTFINDER_ORGANISM" ]] && STAR_ARGS+=(--pointfinder-organism "$STARAMR_POINTFINDER_ORGANISM")
    STAR_ARGS+=("$FASTA_IN")
    
    run_tool_bg "$ENV_STARAMR" "$STAR_LOG" staramr "${STAR_ARGS[@]}"
    PID_STAR=$!

    # 4. VFDB via ABRicate
    local VFDB_TSV="$VFDB_DIR/${SAMPLE}.vfdb.tsv"
    run_tool_bg "$ENV_ABRICATE" "$VFDB_LOG" \
        bash -c "abricate --db vfdb --threads 4 --minid $VFDB_MINID --mincov $VFDB_MINCOV '$FASTA_IN' > '$VFDB_TSV'"
    PID_VFDB=$!

    # 5. PlasmidFinder via ABRicate
    local PLASMID_TSV="$PLASMID_DIR/${SAMPLE}.plasmidfinder.tsv"
    run_tool_bg "$ENV_ABRICATE" "$PLAS_LOG" \
        bash -c "abricate --db plasmidfinder --threads 4 --minid $PLASMIDFINDER_MINID --mincov $PLASMIDFINDER_MINCOV '$FASTA_IN' > '$PLASMID_TSV'"
    PID_PLAS=$!

    # 6. MOB-suite
    run_tool_bg "$ENV_MOBSUITE" "$MOB_LOG" \
        mob_recon --infile "$FASTA_IN" --outdir "$MOB_DIR" --num_threads 4 --force
    PID_MOB=$!

    # ==========================================
    # ĐỢI TOÀN BỘ BACKGROUND JOBS HOÀN TẤT
    # ==========================================
    wait $PID_KLEB || { FAILED_TOOLS+=("Kleborate"); FAIL=1; }
    wait $PID_AMRF || { FAILED_TOOLS+=("AMRFinderPlus"); FAIL=1; }
    wait $PID_STAR || { FAILED_TOOLS+=("StarAMR"); FAIL=1; }
    wait $PID_VFDB || { FAILED_TOOLS+=("VFDB"); FAIL=1; }
    wait $PID_PLAS || { FAILED_TOOLS+=("PlasmidFinder"); FAIL=1; }
    wait $PID_MOB  || { FAILED_TOOLS+=("MOB-suite"); FAIL=1; }

    [[ -e "$KLEB_TSV" ]] || { log "  [ERROR] Kleborate output rỗng/không tồn tại"; FAILED_TOOLS+=("Kleborate_empty"); FAIL=1; }

    if [[ "$DRY_RUN" -eq 1 ]]; then
        log "[DRY-RUN] Hoàn tất kế hoạch cho $SAMPLE"
        return 0
    fi

    {
        echo "=================================================="
        echo "Sample: $SAMPLE | Start: $(date -d @"$START_TIME")"
        echo "=================================================="
        cat "$KLEB_LOG" "$AMRF_LOG" "$STAR_LOG" \
            "$VFDB_LOG" "$PLAS_LOG" "$MOB_LOG" 2>/dev/null || true
    } > "$MAIN_LOG"

    local END_TIME="$(date +%s)"
    local ELAPSED=$((END_TIME - START_TIME))

    local STATUS MSG METRICS
    if [[ "$FAIL" -eq 0 ]]; then
        STATUS="SUCCESS"
        MSG="OK"
        touch "$CKPT_DONE"
    else
        STATUS="PARTIAL"
        MSG="Failed tools: ${FAILED_TOOLS[*]}"
    fi

    METRICS="$(collect_metrics "$KLEB_TSV" "$AMRF_TSV" "$STARAMR_DIR" "$VFDB_TSV" "$PLASMID_TSV" "$MOB_DIR")"
    IFS=$'\t' read -r -a METRIC_ARR <<< "$METRICS"

    update_summary "$SAMPLE" "$SAMPLE" "$GROUP" "$STATUS" "${METRIC_ARR[@]}" "${ELAPSED}s" "$MSG"

    log "  [$STATUS] $SAMPLE hoàn tất trong ${ELAPSED}s"
}

# ---------------------------------------------------------
# 4. AGGREGATE COMBINED OUTPUTS
# ---------------------------------------------------------
aggregate_results() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        log "[DRY-RUN] would aggregate results into $COMBINED_DIR"
        return 0
    fi

    mkdir -p "$COMBINED_DIR"
    python3 - "$FUNC_OUT" "$COMBINED_DIR" "$SUMMARY" <<'PY'
import csv, sys
from pathlib import Path

root = Path(sys.argv[1])
combined = Path(sys.argv[2])
summary = Path(sys.argv[3])
combined.mkdir(parents=True, exist_ok=True)

def iter_sample_dirs():
    for gdir in sorted(root.iterdir()):
        if not gdir.is_dir() or gdir.name.startswith(".") or gdir.name.startswith("_"): continue
        for sdir in sorted(gdir.iterdir()):
            if sdir.is_dir(): yield gdir.name, sdir.name, sdir

def combine(out_name, resolver):
    out_path = combined / out_name
    wrote_header = False
    with out_path.open("w", newline="") as oh:
        writer = None
        for group, sample, sdir in iter_sample_dirs():
            path = resolver(sdir, sample)
            
            # Resolve directory to actual file for Kleborate
            if path.is_dir():
                target = path / "klebsiella_pneumo_complex_output.txt"
                if target.exists() and target.stat().st_size > 0:
                    path = target
                else:
                    for f in path.glob("*output.txt"):
                        if "hAMRonization" not in f.name:
                            path = f
                            break
            
            if not path.exists() or path.is_dir() or path.stat().st_size == 0: continue
            with path.open(encoding="utf-8-sig", errors="replace", newline="") as fh:
                reader = csv.reader(fh, delimiter="\t")
                try: header = next(reader)
                except StopIteration: continue
                if not wrote_header:
                    writer = csv.writer(oh, delimiter="\t", lineterminator="\n")
                    writer.writerow(["sample", "group", "source_file"] + header)
                    wrote_header = True
                for row in reader:
                    if not any(cell.strip() for cell in row): continue
                    if len(row) < len(header): row += [""] * (len(header) - len(row))
                    writer.writerow([sample, group, str(path)] + row)
        if not wrote_header: oh.write("sample\tgroup\tsource_file\n")

combine("kleborate_all.tsv", lambda sdir, sample: sdir / "Kleborate" / f"{sample}_kleborate_out")
combine("amrfinder_all.tsv", lambda sdir, sample: sdir / "AMRFinderPlus" / f"{sample}.amrfinder.tsv")
combine("staramr_summary_all.tsv", lambda sdir, sample: sdir / "StarAMR" / "summary.tsv")
combine("staramr_resfinder_all.tsv", lambda sdir, sample: sdir / "StarAMR" / "resfinder.tsv")
combine("staramr_pointfinder_all.tsv", lambda sdir, sample: sdir / "StarAMR" / "pointfinder.tsv")
combine("vfdb_all.tsv", lambda sdir, sample: sdir / "VFDB" / f"{sample}.vfdb.tsv")
combine("plasmidfinder_all.tsv", lambda sdir, sample: sdir / "PlasmidFinder" / f"{sample}.plasmidfinder.tsv")
combine("mobsuite_mobtyper_all.tsv", lambda sdir, sample: sdir / "MOBSuite" / "mobtyper_results.txt")
combine("mobsuite_contig_report_all.tsv", lambda sdir, sample: sdir / "MOBSuite" / "contig_report.txt")

st_out = combined / "st_for_phylogenomics.tsv"
with st_out.open("w", newline="") as oh:
    cols = ["Sample", "Group", "Status", "Kleborate_species", "Kleborate_ST", "K_locus", "O_locus"]
    writer = csv.DictWriter(oh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    if summary.exists() and summary.stat().st_size > 0:
        with summary.open(encoding="utf-8-sig", errors="replace", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for r in reader: writer.writerow({c: r.get(c, "-") for c in cols})

print(f"[OK] Combined outputs written to: {combined}", file=sys.stderr)
PY
}

# ---------------------------------------------------------
# 5. MAIN
# ---------------------------------------------------------
main() {
    if (( AGGREGATE_ONLY )); then
        aggregate_results
        exit 0
    fi

    [[ -d "$ASM_IN" ]] || die "Không thấy ASM_IN: $ASM_IN"
    mkdir -p "$FUNC_OUT" "$CHECKPOINT_DIR"

    (
        flock -x 200
        [[ -f "$SUMMARY" ]] || echo -e "Sample\tGroup\tStatus\tKleborate_species\tKleborate_ST\tK_locus\tO_locus\tVirulence_score\tResistance_score\tAMRFinder_hits\tAMRFinder_genes\tStarAMR_hits\tVFDB_hits\tVFDB_genes\tPlasmidFinder_hits\tPlasmidFinder_replicons\tMOB_plasmids\tMOB_replicons\tTime\tMessage" > "$SUMMARY"
    ) 200>"${SUMMARY}.lock"

    mapfile -t SAMPLE_PATHS < <(find "$ASM_IN" -mindepth 2 -maxdepth 2 -type d -exec test -f '{}/assembly.fasta' \; -print | sort)
    [[ ${#SAMPLE_PATHS[@]} -eq 0 ]] && die "Không tìm thấy assembly.fasta trong $ASM_IN"

    local TOTAL_SAMPLES="${#SAMPLE_PATHS[@]}"
    local TARGET_PATHS=()

    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        local IDX=$((SLURM_ARRAY_TASK_ID - 1))
        (( IDX >= TOTAL_SAMPLES )) && exit 0
        TARGET_PATHS=("${SAMPLE_PATHS[$IDX]}")
        log "SLURM Array mode: task $SLURM_ARRAY_TASK_ID / $TOTAL_SAMPLES"
    else
        TARGET_PATHS=("${SAMPLE_PATHS[@]}")
        log "Sequential mode: chạy $TOTAL_SAMPLES mẫu"
    fi

    for spath in "${TARGET_PATHS[@]}"; do
        run_sample "$spath"
    done

    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" && "$DRY_RUN" -eq 0 ]]; then
        aggregate_results
    else
        log "Nếu chạy bằng SLURM array, sau khi tất cả task xong hãy chạy: bash $0 --aggregate-only"
    fi
}

main "$@"