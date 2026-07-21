#!/bin/bash

# ================= Script Paths =================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RUN_CNS_SCRIPT="${SCRIPT_DIR}/run_CNs_pipeline.py"
PLOT_CNS_SCRIPT="${SCRIPT_DIR}/plot_CNs_pipeline.R"
HEATMAP_SCRIPT="${PROJECT_ROOT}/CNs_V2/heatmap.R"

# ================= Runtime Config =================
# Update these commands first when moving to a new server.
PYTHON_CMD="python3"
RSCRIPT_CMD="Rscript"

# ================= Default Parameters =================
META_PATH=""
RDATA_PATH=""
OUTPUT_DIR=""
DATA_TYPE="qs"
MODE="radius"
GROUP_COL="Group"
REGION_COL="sample"
SAMPLE_COL="sample"
CLUSTER_COL="predicted.celltype"
X_COL="X"
Y_COL="Y"
RADIUS=80
KNN_K=40
N_CLUSTERS="auto"
ELBOW_MAX=20
METRIC="distortion"
INDEX_COL="0"
OBJECT_NAME=""
OUTPUT_OBJECT="add_CNS.qs"

usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required:"
    echo "  --meta_path        Path to the metadata CSV used for CN analysis"
    echo "  --output_dir       Base output directory"
    echo ""
    echo "Optional plotting stage:"
    echo "  --rdata_path       Path to Seurat object file (.qs/.rds/.RData). If omitted, only the Python CN step runs"
    echo "  --data_type        qs | rds | rdata (Default: qs)"
    echo "  --object_name      Object name to use when data_type=rdata"
    echo "  --output_object    Output qs filename after merging CNs (Default: add_CNS.qs)"
    echo ""
    echo "CN parameters:"
    echo "  --mode             knn or radius (Default: radius)"
    echo "  --group_col        Metadata column used for stacked group plots (Default: Group)"
    echo "  --region_col       Region column used for CN analysis (Default: sample)"
    echo "  --sample_col       Metadata column used for spatial plots (Default: sample)"
    echo "  --cluster_col      Cluster column (Default: predicted.celltype)"
    echo "  --x_col            X coordinate column (Default: X)"
    echo "  --y_col            Y coordinate column (Default: Y)"
    echo "  --index_col        CSV index column for pandas, or none (Default: 0)"
    echo "  --radius           Radius for radius mode (Default: 80)"
    echo "  --knn_k            k for knn mode (Default: 40)"
    echo "  --n_clusters       Integer cluster count or auto (Default: auto)"
    echo "  --elbow_max        Max clusters tested in auto mode (Default: 20)"
    echo "  --metric           Yellowbrick elbow metric (Default: distortion)"
    echo "  -h, --help         Show this message"
    exit 1
}

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --meta_path) META_PATH="$2"; shift 2 ;;
        --rdata_path) RDATA_PATH="$2"; shift 2 ;;
        --output_dir) OUTPUT_DIR="$2"; shift 2 ;;
        --data_type) DATA_TYPE="$2"; shift 2 ;;
        --mode) MODE="$2"; shift 2 ;;
        --group_col) GROUP_COL="$2"; shift 2 ;;
        --region_col) REGION_COL="$2"; shift 2 ;;
        --sample_col) SAMPLE_COL="$2"; shift 2 ;;
        --cluster_col) CLUSTER_COL="$2"; shift 2 ;;
        --x_col) X_COL="$2"; shift 2 ;;
        --y_col) Y_COL="$2"; shift 2 ;;
        --radius) RADIUS="$2"; shift 2 ;;
        --knn_k) KNN_K="$2"; shift 2 ;;
        --n_clusters) N_CLUSTERS="$2"; shift 2 ;;
        --elbow_max) ELBOW_MAX="$2"; shift 2 ;;
        --metric) METRIC="$2"; shift 2 ;;
        --index_col) INDEX_COL="$2"; shift 2 ;;
        --object_name) OBJECT_NAME="$2"; shift 2 ;;
        --output_object) OUTPUT_OBJECT="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
done

if [[ -z "$META_PATH" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: --meta_path and --output_dir are required."
    usage
fi

if [[ ! -f "$RUN_CNS_SCRIPT" ]]; then
    echo "Error: Python pipeline script not found: $RUN_CNS_SCRIPT"
    exit 1
fi

if [[ ! -f "$PLOT_CNS_SCRIPT" ]]; then
    echo "Error: R pipeline script not found: $PLOT_CNS_SCRIPT"
    exit 1
fi

if [[ ! -f "$HEATMAP_SCRIPT" ]]; then
    echo "Error: Heatmap script not found: $HEATMAP_SCRIPT"
    exit 1
fi

CN_OUTPUT_DIR="${OUTPUT_DIR}/cn_result"
PLOT_OUTPUT_DIR="${OUTPUT_DIR}/plot_result"
LOG_DIR="${OUTPUT_DIR}/logs"

mkdir -p "${OUTPUT_DIR}" "${CN_OUTPUT_DIR}" "${PLOT_OUTPUT_DIR}" "${LOG_DIR}"

echo "----------------------------------------"
echo "CNs_ana Pipeline"
echo "----------------------------------------"
echo "Project Root      : ${PROJECT_ROOT}"
echo "Mode              : ${MODE}"
echo "Meta Path         : ${META_PATH}"
echo "Seurat Object     : ${RDATA_PATH:-<skip plotting>}"
echo "Output Dir        : ${OUTPUT_DIR}"
echo "Python Cmd        : ${PYTHON_CMD}"
echo "Rscript Cmd       : ${RSCRIPT_CMD}"
echo "----------------------------------------"

echo "Step 1: Running CN analysis..."
PY_CMD="${PYTHON_CMD} \"${RUN_CNS_SCRIPT}\" \
    --meta_path \"${META_PATH}\" \
    --output_dir \"${CN_OUTPUT_DIR}\" \
    --index_col \"${INDEX_COL}\" \
    --region_col \"${REGION_COL}\" \
    --cluster_col \"${CLUSTER_COL}\" \
    --x_col \"${X_COL}\" \
    --y_col \"${Y_COL}\" \
    --mode \"${MODE}\" \
    --radius \"${RADIUS}\" \
    --knn_k \"${KNN_K}\" \
    --n_clusters \"${N_CLUSTERS}\" \
    --elbow_max \"${ELBOW_MAX}\" \
    --metric \"${METRIC}\""

eval "${PY_CMD}" | tee "${LOG_DIR}/01_run_cns.log"

CN_CSV="${CN_OUTPUT_DIR}/add_CNs.csv"
CONFIG_JSON="${CN_OUTPUT_DIR}/run_config.json"
if [[ ! -f "${CN_CSV}" ]]; then
    echo "Error: CN result CSV not found after Python stage: ${CN_CSV}"
    exit 1
fi
if [[ ! -f "${CONFIG_JSON}" ]]; then
    echo "Error: Config JSON not found after Python stage: ${CONFIG_JSON}"
    exit 1
fi

if [[ -n "${RDATA_PATH}" ]]; then
    echo "Step 2: Running R heatmap and plotting..."
    R_CMD="${RSCRIPT_CMD} \"${PLOT_CNS_SCRIPT}\" \
        --work_dir \"${PLOT_OUTPUT_DIR}\" \
        --rdata_path \"${RDATA_PATH}\" \
        --data_type \"${DATA_TYPE}\" \
        --cns_df \"${CN_CSV}\" \
        --config_path \"${CONFIG_JSON}\" \
        --heatmap_script \"${HEATMAP_SCRIPT}\" \
        --group_col \"${GROUP_COL}\" \
        --sample_col \"${SAMPLE_COL}\" \
        --celltype_col \"${CLUSTER_COL}\" \
        --output_object \"${OUTPUT_OBJECT}\""

    if [[ -n "${OBJECT_NAME}" ]]; then
        R_CMD="${R_CMD} --object_name \"${OBJECT_NAME}\""
    fi

    eval "${R_CMD}" | tee "${LOG_DIR}/02_plot_cns.log"
else
    echo "Step 2: Skip plotting stage because --rdata_path was not provided."
fi

echo "----------------------------------------"
echo "Pipeline finished."
echo "CN outputs   : ${CN_OUTPUT_DIR}"
echo "Plot outputs : ${PLOT_OUTPUT_DIR}"
echo "Logs         : ${LOG_DIR}"
echo "----------------------------------------"
