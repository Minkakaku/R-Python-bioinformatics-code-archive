#!/bin/bash

set -euo pipefail

# --- Default Values ---
# Set the default value for --ncores here
NCORES=12
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_SCRIPT_DIR="$(pwd)"
RCTD_SCRIPT_PATH="${SCRIPT_DIR}/RCTD.R"
LOG_DIR="${PWD}"
RESULT_DIR="${PWD}"
STATUS_DIR="${PWD}"
CONDA_PREFIX_PATH="/mnt/NexData_file/home/lyh/Software/conda_env/lyh"

# --- Help Function ---
print_usage() {
    echo "Usage: $0 -d <target_dir> -r <ref_path> -a <anno_col> -t <type> [-n <ncores>] [-y <spatial_type>] [--script-dir <dir>] [--rctd-script <path>]"
    echo ""
    echo "Generates an RCTD run script for each sample sub-directory within the target directory."
    echo ""
    echo "Required Arguments:"
    echo "  -d, --target-dir  Directory containing the sample sub-directories (e.g., sample1, sample2...)"
    echo "  -r, --ref           Reference file path (passed to RCTD.R as --ref)"
    echo "  -a, --anno-col    Name of the annotation column (passed to RCTD.R as --anno_column)"
    echo "  -t, --type          Data type of the reference file (passed to RCTD.R as --type)"
    echo "  -y, --spatial-type  Data type of per-sample spatial Seurat objects (default: rds)"
    echo ""
    echo "Optional Arguments:"
    echo "  -n, --ncores      Number of cores to use (Default: ${NCORES})"
    echo "      --script-dir  Directory used to save generated shell scripts (default: current directory)"
    echo "      --rctd-script Path to RCTD.R (default: same directory as this script)"
    echo "      --log-dir     Directory used to save per-sample logs"
    echo "      --result-dir  Directory used to save successful per-sample RCTD.txt copies"
    echo "      --status-dir  Directory used to save per-sample status files"
    echo "      --conda-prefix Conda prefix used by generated scripts"
    echo "  -h, --help          Show this help message"
    exit 1
}

# --- Parse Command-Line Arguments ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -d|--target-dir) TARGET_DIR="$2"; shift ;;
        -r|--ref) REF_PATH="$2"; shift ;;
        -a|--anno-col) ANNO_COLUMN="$2"; shift ;;
        -t|--type) DATA_TYPE="$2"; shift ;;
        -y|--spatial-type) SPATIAL_TYPE="$2"; shift ;;
        -n|--ncores) NCORES="$2"; shift ;;
        --script-dir) JOB_SCRIPT_DIR="$2"; shift ;;
        --rctd-script) RCTD_SCRIPT_PATH="$2"; shift ;;
        --log-dir) LOG_DIR="$2"; shift ;;
        --result-dir) RESULT_DIR="$2"; shift ;;
        --status-dir) STATUS_DIR="$2"; shift ;;
        --conda-prefix) CONDA_PREFIX_PATH="$2"; shift ;;
        -h|--help) print_usage ;;
        *) echo "Unknown parameter: $1"; print_usage ;;
    esac
    shift
done

SPATIAL_TYPE="${SPATIAL_TYPE:-rds}"

# --- Validate Required Arguments ---
if [[ -z "$TARGET_DIR" || -z "$REF_PATH" || -z "$ANNO_COLUMN" || -z "$DATA_TYPE" ]]; then
    echo "Error: Missing one or more required arguments."
    echo "(-d, -r, -a, and -t are all required)"
    echo ""
    print_usage
fi

# --- Validate Target Directory ---
if [[ ! -d "$TARGET_DIR" ]]; then
    echo "Error: Target directory not found: $TARGET_DIR"
    exit 1
fi

if [[ ! -f "$RCTD_SCRIPT_PATH" ]]; then
    echo "Error: RCTD.R not found: $RCTD_SCRIPT_PATH"
    exit 1
fi

mkdir -p "$JOB_SCRIPT_DIR"
JOB_SCRIPT_DIR="$(cd "$JOB_SCRIPT_DIR" && pwd)"
mkdir -p "$LOG_DIR" "$RESULT_DIR" "$STATUS_DIR"
LOG_DIR="$(cd "$LOG_DIR" && pwd)"
RESULT_DIR="$(cd "$RESULT_DIR" && pwd)"
STATUS_DIR="$(cd "$STATUS_DIR" && pwd)"

# --- Print Settings ---
echo "Starting script generation..."
echo "=============================="
echo "Target Directory: $TARGET_DIR"
echo "Reference File:   $REF_PATH"
echo "Annotation Column: $ANNO_COLUMN"
echo "Data Type:        $DATA_TYPE"
echo "Cores:            $NCORES"
echo "Output Script Dir: $JOB_SCRIPT_DIR"
echo "RCTD Script Path:  $RCTD_SCRIPT_PATH"
echo "Log Directory:     $LOG_DIR"
echo "Result Directory:  $RESULT_DIR"
echo "Status Directory:  $STATUS_DIR"
echo "Conda Prefix:      $CONDA_PREFIX_PATH"
echo "=============================="

# --- Main Loop ---
# Loop through all sub-directories in $TARGET_DIR
# "${TARGET_DIR}"/*/ ensures we only match directories
for dir_path in "${TARGET_DIR}"/*/; do
    
    # Check if it is actually a directory
    if [[ ! -d "$dir_path" ]]; then
        continue
    fi
    
    # Get the plain name of the sub-directory (e.g., "sample1")
    subdir_name=$(basename "$dir_path")
    
    # Resolve to an absolute path, which is safer for Singularity/R
    # `readlink -f` converts any path (relative or absolute) to a full absolute path
    abs_dir_path=$(readlink -f "$dir_path")
    
    # Construct paths for required files
    spatial_obj_path="${abs_dir_path}/spatial_seurat.rds"
    output_path="${abs_dir_path}" # Output is the directory itself

    # Check if the per-sample Seurat object exists
    if [[ -f "$spatial_obj_path" ]]; then
        
        # Create the script file in the current working directory
        script_name="${subdir_name}_runRCTD.sh"
        script_path="${JOB_SCRIPT_DIR}/${script_name}"
        job_log="${LOG_DIR}/${subdir_name}_runRCTD.log"
        status_file="${STATUS_DIR}/${subdir_name}.status"
        error_file="${STATUS_DIR}/${subdir_name}.error"
        result_copy_path="${RESULT_DIR}/${subdir_name}_RCTD.txt"
        
        echo "Creating script: $script_path (for $subdir_name)"

        # Write the shebang and command to the script file
        # Note: All variables are quoted to handle spaces in paths, etc.
        cat > "$script_path" <<EOF
#!/bin/bash

set -euo pipefail

# Generated by generate_rctd_jobs_en.sh on $(date)
SAMPLE_NAME="${subdir_name}"
SAMPLE_DIR="${abs_dir_path}"
LOG_FILE="${job_log}"
STATUS_FILE="${status_file}"
ERROR_FILE="${error_file}"
RESULT_COPY_PATH="${result_copy_path}"

mkdir -p "\$(dirname "\$LOG_FILE")" "\$(dirname "\$STATUS_FILE")" "\$(dirname "\$RESULT_COPY_PATH")"
rm -f "\$STATUS_FILE" "\$ERROR_FILE"

{
echo "Starting RCTD run for sample [\${SAMPLE_NAME}]..."
echo "Log file: \${LOG_FILE}"
echo "Sample directory: \${SAMPLE_DIR}"

set +e

conda run -p "${CONDA_PREFIX_PATH}" Rscript "${RCTD_SCRIPT_PATH}" \\
    --ref "${REF_PATH}" \\
    --anno_column "${ANNO_COLUMN}" \\
    --type "${DATA_TYPE}" \\
    --spatial_obj "${spatial_obj_path}" \\
    --spatial_type "${SPATIAL_TYPE}" \\
    --ncores ${NCORES} \\
    --output "${output_path}"

cmd_exit=\$?
set -e

if [[ \$cmd_exit -eq 0 && -f "${output_path}/RCTD.txt" ]]; then
    cp "${output_path}/RCTD.txt" "\$RESULT_COPY_PATH"
    echo "SUCCESS" > "\$STATUS_FILE"
    echo "Copied final result to: \$RESULT_COPY_PATH"
    rm -rf "\$SAMPLE_DIR"
    echo "Deleted intermediate sample directory: \$SAMPLE_DIR"
    echo "RCTD run for sample [\${SAMPLE_NAME}] finished successfully."
    exit 0
fi

error_reason=\$(grep -E "Error|Execution halted|cannot|failed|Failure" "\$LOG_FILE" 2>/dev/null | tail -n 20 || true)
if [[ -z "\$error_reason" ]]; then
    error_reason="RCTD command exited with code \$cmd_exit. Check the full log for details."
fi

echo "FAILED" > "\$STATUS_FILE"
printf '%s\n' "\$error_reason" > "\$ERROR_FILE"
echo "RCTD run for sample [\${SAMPLE_NAME}] failed."
echo "Error summary saved to: \$ERROR_FILE"
echo "\$error_reason"
exit 1
} >> "\$LOG_FILE" 2>&1
EOF
        
        # Add execute permission to the generated script
        chmod +x "$script_path"
        
    else
        echo "Warning: Skipping $subdir_name. Missing spatial_seurat.rds file."
    fi
done

echo "=============================="
echo "Script generation complete."
