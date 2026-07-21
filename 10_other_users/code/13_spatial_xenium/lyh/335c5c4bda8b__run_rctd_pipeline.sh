#!/usr/bin/env bash
set -euo pipefail

WORKDIR=${WORKDIR:-/mnt/NexData_file/home/lyh/Test/HD_project/RCTD_py}
SCRIPT=${SCRIPT:-/mnt/NexData_file/home/lyh/Script/Xenium/RCTD_py/rctd_pipeline.py}
CONTAINER=${CONTAINER:-jupyterlab_8972}
PYTHON=${PYTHON:-/usr/local/sc_project/.venv/bin/python3}
MODE=${MODE:-multi}
MAX_MULTI_TYPES=${MAX_MULTI_TYPES:-1}
BATCH_SIZE=${BATCH_SIZE:-2000}
UMI_MIN=${UMI_MIN:-100}
COUNTS_MIN=${COUNTS_MIN:-10}
DEVICE=${DEVICE:-cuda}
SKIP_COUNT_CHECK=0
SPATIAL=""
REFERENCE=""
CELL_TYPE_COL=""
LOG=""

usage() {
  cat <<'EOF'
Usage:
  ./run_rctd_pipeline.sh --spatial spatial.h5ad --reference ref.h5ad --cell-type-col cluster \
    --workdir /path/to/workdir [options]

Required:
  --spatial PATH                 Processed spatial h5ad. .X should be raw counts.
  --reference PATH               Reference h5ad. .X should be raw counts.
  --cell-type-col COL            Cell type column in reference.obs.

Options:
  --workdir PATH                 Working directory for logs and outputs. Created if missing.
                                 Default: /mnt/NexData_file/home/lyh/Test/HD_project/RCTD_py
  --mode full|doublet|multi      Default: multi
  --max-multi-types N            Default: 1; only affects multi mode
  --batch-size N                 Default: 2000
  --umi-min N                    RCTD UMI_min. Default: 100
  --counts-min N                 RCTD counts_MIN. Default: 10
  --device cuda|cpu|auto         Default: cuda
  --skip-count-check             Do not check whether .X looks like raw counts
  --container NAME               Default: jupyterlab_8972
  --python PATH                  Python inside container
  --script PATH                  Python pipeline path inside container
                                 Default: /mnt/NexData_file/home/lyh/Script/Xenium/RCTD_py/rctd_pipeline.py
  --log PATH                     Log path
  -h, --help                     Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --spatial) SPATIAL="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --cell-type-col) CELL_TYPE_COL="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --max-multi-types) MAX_MULTI_TYPES="$2"; shift 2 ;;
    --batch-size) BATCH_SIZE="$2"; shift 2 ;;
    --umi-min) UMI_MIN="$2"; shift 2 ;;
    --counts-min) COUNTS_MIN="$2"; shift 2 ;;
    --device) DEVICE="$2"; shift 2 ;;
    --skip-count-check) SKIP_COUNT_CHECK=1; shift ;;
    --workdir) WORKDIR="$2"; shift 2 ;;
    --container) CONTAINER="$2"; shift 2 ;;
    --python) PYTHON="$2"; shift 2 ;;
    --script) SCRIPT="$2"; shift 2 ;;
    --log) LOG="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

case "$MODE" in full|doublet|multi) ;; *) echo "Invalid --mode: $MODE" >&2; exit 2 ;; esac
case "$DEVICE" in cuda|cpu|auto) ;; *) echo "Invalid --device: $DEVICE" >&2; exit 2 ;; esac

for item in SPATIAL REFERENCE CELL_TYPE_COL; do
  if [[ -z "${!item}" ]]; then
    echo "Missing required --$(echo "$item" | tr 'A-Z_' 'a-z-')" >&2
    usage >&2
    exit 2
  fi
done

mkdir -p "$WORKDIR"

STAMP=$(date +%Y%m%d_%H%M%S)
if [[ -z "$LOG" ]]; then
  LOG="$WORKDIR/rctd_${MODE}_pipeline_${STAMP}.log"
fi

cd "$WORKDIR"

echo "WORKDIR=$WORKDIR"
echo "SCRIPT=$SCRIPT"
echo "CONTAINER=$CONTAINER"
echo "PYTHON=$PYTHON"
echo "SPATIAL=$SPATIAL"
echo "REFERENCE=$REFERENCE"
echo "CELL_TYPE_COL=$CELL_TYPE_COL"
echo "MODE=$MODE"
echo "MAX_MULTI_TYPES=$MAX_MULTI_TYPES"
echo "BATCH_SIZE=$BATCH_SIZE"
echo "UMI_MIN=$UMI_MIN"
echo "COUNTS_MIN=$COUNTS_MIN"
echo "DEVICE=$DEVICE"
echo "LOG=$LOG"

COUNT_CHECK_ARG=""
if [[ "$SKIP_COUNT_CHECK" == "1" ]]; then
  COUNT_CHECK_ARG="--skip-count-check"
fi

docker exec "$CONTAINER" bash -lc \
  "mkdir -p '$WORKDIR' && cd '$WORKDIR' && '$PYTHON' '$SCRIPT' \
    --spatial '$SPATIAL' \
    --reference '$REFERENCE' \
    --cell-type-col '$CELL_TYPE_COL' \
    --workdir '$WORKDIR' \
    --mode '$MODE' \
    --max-multi-types '$MAX_MULTI_TYPES' \
    --batch-size '$BATCH_SIZE' \
    --umi-min '$UMI_MIN' \
    --counts-min '$COUNTS_MIN' \
    --device '$DEVICE' \
    $COUNT_CHECK_ARG" 2>&1 | tee "$LOG"
