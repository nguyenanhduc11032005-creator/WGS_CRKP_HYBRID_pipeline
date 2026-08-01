#!/bin/bash
#SBATCH --job-name=download_dbs
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=_dl_dbs_%j.out
#SBATCH --error=_dl_dbs_%j.err

set -euo pipefail

# ==============================================================================
# 1. Trap bắt lỗi và Helper Functions
# ==============================================================================
trap 'log "ERROR at line $LINENO"; exit 1' ERR
trap 'log "Received SIGINT/SIGTERM"; exit 130' INT TERM

ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }
die() { log "FATAL: $*"; exit 1; }

# ==============================================================================
# 2. Xử lý Conda thông minh
# ==============================================================================
TARGET_ENV="fastq_dl"
if [[ "${CONDA_DEFAULT_ENV:-}" != "$TARGET_ENV" ]]; then
    if command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)"
        conda activate "${CONDA_ENV:-$TARGET_ENV}" || die "Cannot activate conda environment: $TARGET_ENV"
    else
        log "WARNING: Conda not found in PATH. Proceeding with system tools..."
    fi
fi

# ==============================================================================
# 3. Cấu hình biến & Đường dẫn
# ==============================================================================
ROOT_DIR="/mnt/d18t/ANHDUC/WGS_CRKP_HYBRID"

BAKTA_DB_DIR="${ROOT_DIR}/bakta_db"
BAKTA_URL="https://zenodo.org/records/14916843/files/db.tar.xz"
BAKTA_ARCHIVE="${BAKTA_DB_DIR}/bakta_db_v6.tar.xz"

RFAM_DB_DIR="${ROOT_DIR}/rfam_db"
RFAM_CM_URL="https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
RFAM_CLAN_URL="https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"

# ==============================================================================
# 4. Kiểm tra Dependencies
# ==============================================================================
log "Checking system dependencies..."
for cmd in aria2c tar df gunzip awk; do
    command -v "$cmd" >/dev/null 2>&1 || die "$cmd not found in PATH"
done

# ==============================================================================
# 5. Tải & Giải nén Bakta Database
# ==============================================================================
mkdir -p "${BAKTA_DB_DIR}"
log "Working directory for Bakta: ${BAKTA_DB_DIR}"

AVAILABLE_GB=$(df -BG "${BAKTA_DB_DIR}" | awk 'NR==2 {gsub("G","",$4); print $4}')
log "Available disk space: ${AVAILABLE_GB} GB"
if (( AVAILABLE_GB < 70 )); then
    die "Need at least ~70 GB free space for Bakta DB."
fi

if find "${BAKTA_DB_DIR}" -maxdepth 1 -type d -name "db*" | grep -q .; then
    log "Existing Bakta database detected. Skipping Bakta download."
else
    log "Downloading Bakta Full Database v6.0..."
    aria2c \
        --continue=true \
        --max-connection-per-server=16 \
        --split=16 \
        --min-split-size=10M \
        --file-allocation=none \
        --dir="${BAKTA_DB_DIR}" \
        --out="$(basename "${BAKTA_ARCHIVE}")" \
        "${BAKTA_URL}"
    
    if [[ ! -s "${BAKTA_ARCHIVE}" ]]; then
        die "Downloaded archive is empty."
    fi
    log "Archive size: $(du -sh "${BAKTA_ARCHIVE}" | cut -f1)"
    
    log "Extracting Bakta database..."
    tar \
        --extract \
        --xz \
        --file="${BAKTA_ARCHIVE}" \
        --directory="${BAKTA_DB_DIR}" \
        --checkpoint=10000 \
        --checkpoint-action=echo="%u files extracted"
    
    log "Extraction completed. Removing archive..."
    rm -f "${BAKTA_ARCHIVE}"
fi

# ==============================================================================
# 6. Tải & Giải nén Rfam Database
# ==============================================================================
mkdir -p "${RFAM_DB_DIR}"
log "Working directory for Rfam: ${RFAM_DB_DIR}"

if [[ -f "${RFAM_DB_DIR}/Rfam.cm" ]]; then
    log "Existing Rfam database detected. Skipping Rfam download."
else
    log "Downloading Rfam Database..."
    aria2c \
        --continue=true \
        --max-connection-per-server=16 \
        --split=16 \
        --file-allocation=none \
        --dir="${RFAM_DB_DIR}" \
        "${RFAM_CM_URL}" \
        "${RFAM_CLAN_URL}"

    log "Extracting Rfam.cm.gz..."
    gunzip -f "${RFAM_DB_DIR}/Rfam.cm.gz"
fi

# ==============================================================================
# 7. Tổng kết (Summary)
# ==============================================================================
log "Database installation completed successfully."
log "Bakta total size: $(du -sh "${BAKTA_DB_DIR}" | cut -f1)"
log "Rfam total size: $(du -sh "${RFAM_DB_DIR}" | cut -f1)"