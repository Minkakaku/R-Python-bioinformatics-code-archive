#!/bin/bash

# ==============================================================================
# Script Name: run_xenium_pipeline.sh
# Description: Automates the post-processing pipeline for Xenium spatial data.
#              1. Sanitizes ROI CSV headers.
#              2. Python: Formats data -> mat.h5, meta.csv
#              3. R: Seurat processing -> raw_count.qs
#              4. Cleanup: Removes intermediate files (mat.h5, meta.csv)
# ==============================================================================

# --- Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
SIF_IMAGE="/mnt/sif/my_sandbox/"
PYTHON_SCRIPT="$SCRIPT_DIR/deal_with_data.py"

# --- Safety Settings ---
set -euo pipefail

# --- Helper Functions ---
log_info() { echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $1"; }
log_error() { echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] [ERROR] $1" >&2; }

usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -r <roi_dir> -c <chip_id> [-g <group_csv>]"
    echo ""
    echo "Required Arguments:"
    echo "  -i  Path to raw Xenium output directory (Input)"
    echo "  -o  Path to output directory (Output)"
    echo "  -r  Path to directory containing ROI CSV files"
    echo "  -c  Chip ID string"
    echo ""
    echo "Optional Arguments:"
    echo "  -g  Path to Group/Sample mapping CSV file (Optional)"
    echo "  -h  Show this help message"
    echo ""
    exit 1
}

# --- Argument Parsing ---
GROUP_FILE=""

while getopts "i:o:r:c:g:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        r) ROI_DIR="$OPTARG" ;;
        c) CHIP_ID="$OPTARG" ;;
        g) GROUP_FILE="$OPTARG" ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Validate required arguments
if [[ -z "${INPUT_DIR:-}" ]] || [[ -z "${OUTPUT_DIR:-}" ]] || [[ -z "${ROI_DIR:-}" ]] || [[ -z "${CHIP_ID:-}" ]]; then
    log_error "Missing required arguments."
    usage
fi

# --- Main Pipeline ---

if [[ ! -f "$PYTHON_SCRIPT" ]]; then
    log_error "Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

if [[ ! -d "$SIF_IMAGE" ]] && [[ ! -f "$SIF_IMAGE" ]]; then
    log_error "Singularity image path not found: $SIF_IMAGE"
    exit 1
fi

log_info "Starting Xenium Post-Processing Pipeline"
log_info "----------------------------------------"
log_info "Project Root     : $PROJECT_ROOT"
log_info "Input Directory  : $INPUT_DIR"
log_info "Output Directory : $OUTPUT_DIR"
log_info "ROI Directory    : $ROI_DIR"
log_info "Chip ID          : $CHIP_ID"
if [[ -n "$GROUP_FILE" ]]; then
    log_info "Group File       : $GROUP_FILE"
else
    log_info "Group File       : None (Skipped)"
fi
log_info "----------------------------------------"

# 1. Sanitize CSV Headers
log_info "Step 1: Sanitizing ROI CSV headers..."
if [[ -d "$ROI_DIR" ]]; then
    find "$ROI_DIR" -name "*.csv" -exec sed -i '/^\s*#/d' {} \;
    log_info "Headers cleaned."
else
    log_error "ROI Directory not found: $ROI_DIR"
    exit 1
fi

# 2. Python Processing
log_info "Step 2: Running Python formatting..."

# Construct Python command
PY_CMD="python3.9 $PYTHON_SCRIPT \
    --input_dir '$INPUT_DIR' \
    --output_dir '$OUTPUT_DIR' \
    --roi_dir '$ROI_DIR' \
    --chip_id '$CHIP_ID'"

if [[ -n "$GROUP_FILE" ]]; then
    PY_CMD="$PY_CMD --group_file '$GROUP_FILE'"
fi

singularity exec --bind /cloud2,/mnt,/public4.2,/public4.4,/public5.4,/public6.1,/public6.2,/public6.4,/public7.1,/public7.2,/public7.3,/public7.4,/public7.5,/public9,/public12,/public15,/public20,/public21,/public23,/public28,/public31.2,/public31.3,/public31.4,/public38,/public40,/public42,/publicsc01,/publicsc02,/publicsc03,/publicsc04,/publicsc05,/publicsc06,/publicsc08,/publicsc09,/publicsc10,/publicsc11,/publicsc12,/publicsc13,/publicsc14,/publicsc15,/publicsc16,/publicsc17,/public80 "$SIF_IMAGE" bash -c "
    source /etc/profile && \
    cd '$PROJECT_ROOT' && \
    $PY_CMD
"
log_info "Python processing done."

log_info "----------------------------------------"
log_info "Pipeline finished successfully."
