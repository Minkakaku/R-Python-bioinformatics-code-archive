#!/usr/bin/env bash
set -euo pipefail

SIF="/mnt/NexData_file/home/lyh/Software/sif/lyh.sif"
BIND="/mnt"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${SCRIPT_DIR}/monocle2_pipeline.py"

usage() {
    cat <<'EOF'
Usage:
  run_monocle2.sh -i INPUT.h5ad -o OUTDIR [options]

Required:
  -i, --input PATH              Input h5ad file
  -o, --outdir PATH             Output directory

Container:
  --sif PATH                    Singularity image
                                default: /mnt/NexData_file/home/lyh/Software/sif/lyh.sif
  --bind PATH                   Singularity bind path
                                default: /mnt

Common options:
  --threads INT                 CPU threads, default: 24
  --mean-expr-thresh FLOAT      Ordering gene mean expression threshold, default: 0.1
  --max-genes INT               Max ordering genes, default: 1000
  --ncenter INT                 DDRTree center number; omit to use auto
  --param-gamma FLOAT           DDRTree gamma, default: 10
  --sigma FLOAT                 DDRTree sigma; omit to use default
  --lambda-param FLOAT          DDRTree lambda; omit to use default
  --subtype-key STR             Cell type column, default: subtype
  --group-key STR               Group column, default: group
  --trajectory-colors A,B,C     Obs columns for trajectory plots, default: State,subtype,group
  --root-state STATE            Root state for order_cells
  --reverse                     Reverse pseudotime direction
  --skip-de                     Skip differential test and dynamic heatmap
  --skip-h5ad                   Do not save final h5ad
  --verbose                     Print verbose Monocle messages

Example:
  bash run_monocle2.sh \
    -i /mnt/NexData_file/home/lyh/project/adata.h5ad \
    -o /mnt/NexData_file/home/lyh/project/monocle2_out \
    --ncenter 290 \
    --param-gamma 10 \
    --subtype-key subtype \
    --group-key group
EOF
}

INPUT=""
OUTDIR=""
ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT="$2"
            shift 2
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --sif)
            SIF="$2"
            shift 2
            ;;
        --bind)
            BIND="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --reverse|--skip-de|--skip-h5ad|--verbose)
            ARGS+=("$1")
            shift
            ;;
        --threads|--mean-expr-thresh|--max-genes|--max-components|--method|--ncenter|--param-gamma|--sigma|--lambda-param|--root-state|--state-key|--subtype-key|--group-key|--trajectory-colors|--cell-size|--split-size|--de-qval|--heatmap-top-features|--heatmap-splits)
            ARGS+=("$1" "$2")
            shift 2
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "${INPUT}" || -z "${OUTDIR}" ]]; then
    echo "ERROR: --input and --outdir are required." >&2
    usage >&2
    exit 1
fi

mkdir -p "${OUTDIR}"

singularity exec --bind "${BIND}" "${SIF}" \
    python "${PIPELINE}" \
    --input "${INPUT}" \
    --outdir "${OUTDIR}" \
    "${ARGS[@]}"
