#!/usr/bin/env bash

set -euo pipefail

CONDA_ENV_ROOT="/mnt/NexData_file/home/lyh/Software/conda_env/lyh"
CONDA_BIN="/mnt/NexData_file/home/lyh/Software/miniconda/bin/conda"

SAMPLE_LIST="sample_list.csv"
COMPARE_LIST="compare_list.csv"
SINGLE_SCRIPT="/mnt/NexData_file/home/lyh/Script/scRNA/Cellchat/Cellchat_single.R"
COMPARE_SCRIPT="/mnt/NexData_file/home/lyh/Script/scRNA/Cellchat/Cellchat_compare.R"
REPORT_SCRIPT="/mnt/NexData_file/home/lyh/Script/scRNA/Cellchat/render_cellchat_report.R"
REPORT_RMD="/mnt/NexData_file/home/lyh/Script/scRNA/Cellchat/cellchat_report.Rmd"
OUTDIR="pipeline_output"
SPECIES="human"
CELLTYPE_COL="celltype"
INPUT_OBJECT=""
SAMPLE_COL=""
CORES="1"
RENDER_ONLY="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample_list) SAMPLE_LIST="$2"; shift 2 ;;
    --compare_list) COMPARE_LIST="$2"; shift 2 ;;
    --single_script) SINGLE_SCRIPT="$2"; shift 2 ;;
    --compare_script) COMPARE_SCRIPT="$2"; shift 2 ;;
    --report_script) REPORT_SCRIPT="$2"; shift 2 ;;
    --report_rmd) REPORT_RMD="$2"; shift 2 ;;
    -o|--outdir) OUTDIR="$2"; shift 2 ;;
    -s|--species) SPECIES="$2"; shift 2 ;;
    -c|--celltype_col) CELLTYPE_COL="$2"; shift 2 ;;
    --input_object) INPUT_OBJECT="$2"; shift 2 ;;
    --sample_col) SAMPLE_COL="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --conda_env) CONDA_ENV_ROOT="$2"; shift 2 ;;
    --conda_bin) CONDA_BIN="$2"; shift 2 ;;
    --render_only) RENDER_ONLY="true"; shift 1 ;;
    -h|--help)
      cat <<'EOF'
Usage: bash run_cellchat_pipeline.sh [options]

Required:
  --input_object PATH
      Full Seurat object file (.rds/.qs/.RData). The script will subset this
      object sample-by-sample before running CellChat.

  --sample_col STRING
      Metadata column used to subset samples from the full Seurat object,
      e.g. orig.ident.

  --sample_list PATH
      CSV describing which samples should be analyzed.
      Required column:
        sample
      Logic:
        sample -> used as output sample name
        sample -> also used as subset value matched in --sample_col
      This file should only contain one column.

Options:
  --compare_list PATH
      CSV describing sample-level comparisons.
      Required columns when not empty:
        sample1, sample2
      If the file is empty or only contains the header, compare analysis is skipped.
      Each pair will run one compare task: sample1 vs sample2.

  --single_script PATH
      Single-sample CellChat R script. Default: Cellchat_single.R

  --compare_script PATH
      Compare CellChat R script. Default: Cellchat_compare.R

  --report_script PATH
      Report rendering R script. Default: render_cellchat_report.R

  --report_rmd PATH
      Report template Rmd. Default: cellchat_report.Rmd

  -o, --outdir PATH
      Output root directory. Default: pipeline_output

  -s, --species STRING
      Species passed to CellChat. Default: human

  -c, --celltype_col STRING
      Metadata column containing cell type annotation. Default: celltype

  --cores INT
      Worker count passed to CellChat single analysis. Default: 4

  --conda_env PATH
      Conda env path used with 'conda run -p'.
      Default: /mnt/NexData_file/home/lyh/Software/conda_env/lyh

  --conda_bin PATH
      Path to conda executable. Default: \$HOME/miniconda3/bin/conda

  --render_only
      Do not run analysis. Only regenerate command script, manifests and report step.
      Failed tasks are recorded in logs/failed_tasks.tsv.

Examples:
  bash run_cellchat_pipeline.sh \
    --input_object all_samples.qs \
    --sample_col orig.ident \
    --sample_list sample_list.csv \
    --compare_list compare_list.csv \
    -o pipeline_output

  sample_list.csv example:
    sample
    AA
    BB
    AB
    BA

  compare_list.csv example:
    sample1,sample2
    AA,AB
    BB,BA
EOF
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

if [[ ! -x "${CONDA_BIN}" ]]; then
  echo "conda not found or not executable: ${CONDA_BIN}" >&2
  exit 1
fi

if [[ ! -d "${CONDA_ENV_ROOT}" ]]; then
  echo "conda env path not found: ${CONDA_ENV_ROOT}" >&2
  exit 1
fi

if [[ -z "${INPUT_OBJECT}" ]]; then
  echo "--input_object is required. It should point to the full Seurat object (.rds/.qs/.RData)." >&2
  exit 1
fi

if [[ -z "${SAMPLE_COL}" ]]; then
  echo "--sample_col is required. It should be the metadata column used for subset, e.g. orig.ident." >&2
  exit 1
fi

if [[ ! -f "${INPUT_OBJECT}" ]]; then
  echo "input object not found: ${INPUT_OBJECT}" >&2
  exit 1
fi

WORKDIR="$(pwd)"
mkdir -p "${OUTDIR}/logs" "${OUTDIR}/metadata" "${OUTDIR}/single" "${OUTDIR}/compare"

PIPELINE_SCRIPT="${OUTDIR}/run_pipeline.sh"
LOG_DIR="${OUTDIR}/logs"
SAMPLE_MANIFEST="${OUTDIR}/metadata/sample_manifest.csv"
COMPARE_MANIFEST="${OUTDIR}/metadata/compare_manifest.csv"

{
  echo "#!/usr/bin/env bash"
  echo "set -euo pipefail"
  echo
} > "${PIPELINE_SCRIPT}"

FAILED_TASKS_FILE="${OUTDIR}/logs/failed_tasks.tsv"
echo -e "task_type\tname\tlog_file" > "${FAILED_TASKS_FILE}"

run_and_log() {
  local log_file="$1"
  shift
  "$@" >"${log_file}" 2>&1
}

conda_exec() {
  "${CONDA_BIN}" run -p "${CONDA_ENV_ROOT}" --no-capture-output "$@"
}

append_cmd() {
  printf '%q ' "$@" >> "${PIPELINE_SCRIPT}"
  printf '\n' >> "${PIPELINE_SCRIPT}"
}

record_failure() {
  local task_type="$1"
  local name="$2"
  local log_file="$3"
  printf '%s\t%s\t%s\n' "${task_type}" "${name}" "${log_file}" >> "${FAILED_TASKS_FILE}"
}

resolve_abs() {
  conda_exec python - "$1" <<'PY'
import os
import sys
print(os.path.abspath(sys.argv[1]))
PY
}

parse_samples() {
  conda_exec python - "${SAMPLE_LIST}" <<'PY'
import csv
import os
import sys

sample_list = sys.argv[1]

if not os.path.exists(sample_list):
    raise SystemExit(f"sample list not found: {sample_list}")

with open(sample_list, newline="", encoding="utf-8-sig") as fh:
    reader = csv.DictReader(fh)
    if reader.fieldnames is None:
        raise SystemExit("sample_list.csv is empty")
    field_map = {name.strip().lower(): name for name in reader.fieldnames}

    def pick(*names):
        for name in names:
            if name in field_map:
                return field_map[name]
        return None

    sample_col = pick("sample", "sample_id", "name")

    if sample_col is None:
        raise SystemExit("sample_list.csv must contain sample/sample_id/name column")

    rows = list(reader)
    if not rows:
        raise SystemExit("sample_list.csv has no data rows")

    for row in rows:
        sample = (row.get(sample_col) or "").strip()
        if not sample:
            continue
        print(sample)
PY
}

parse_compares() {
  conda_exec python - "${COMPARE_LIST}" <<'PY'
import csv
import os
import sys

compare_list = sys.argv[1]
if not os.path.exists(compare_list):
    raise SystemExit(0)

with open(compare_list, newline="", encoding="utf-8-sig") as fh:
    content = fh.read().strip()
    if not content or content == "sample1,sample2":
        raise SystemExit(0)

with open(compare_list, newline="", encoding="utf-8-sig") as fh:
    reader = csv.DictReader(fh)
    if reader.fieldnames is None:
        raise SystemExit(0)
    field_map = {name.strip().lower(): name for name in reader.fieldnames}

    def pick(*names):
        for name in names:
            if name in field_map:
                return field_map[name]
        return None

    group1_col = pick("sample1", "left", "sample_left", "group1", "control", "ctrl")
    group2_col = pick("sample2", "right", "sample_right", "group2", "treat", "case")

    if group1_col is None or group2_col is None:
        raise SystemExit("compare_list.csv must contain sample1 and sample2 columns")

    for row in reader:
        group1 = (row.get(group1_col) or "").strip()
        group2 = (row.get(group2_col) or "").strip()
        if group1 and group2:
            print("\t".join([group1, group2]))
PY
}

declare -a SAMPLE_ROWS=()
while IFS= read -r line; do
  [[ -z "${line}" ]] && continue
  SAMPLE_ROWS+=("${line}")
done < <(parse_samples)

if [[ ${#SAMPLE_ROWS[@]} -eq 0 ]]; then
  echo "No valid sample rows found." >&2
  exit 1
fi

echo "sample,input,subset_col,subset_group" > "${SAMPLE_MANIFEST}"

declare -A SAMPLE_EXISTS=()

for row in "${SAMPLE_ROWS[@]}"; do
  sample="${row}"
  subset_group="${row}"
  sample_outdir="${OUTDIR}/single/${sample}"
  mkdir -p "${sample_outdir}"
  log_file="${LOG_DIR}/single_${sample}.log"

  SAMPLE_EXISTS["${sample}"]="1"

  printf '"%s","%s","%s","%s"\n' "${sample}" "$(resolve_abs "${INPUT_OBJECT}")" "${SAMPLE_COL}" "${subset_group}" >> "${SAMPLE_MANIFEST}"

  cmd=(
    "${CONDA_BIN}" run -p "${CONDA_ENV_ROOT}" --no-capture-output Rscript "${SINGLE_SCRIPT}"
    -i "${INPUT_OBJECT}"
    -o "${sample_outdir}"
    -s "${SPECIES}"
    -c "${CELLTYPE_COL}"
    --python_path "${CONDA_ENV_ROOT}/bin/python"
    --cores "${CORES}"
    --subset_col "${SAMPLE_COL}"
    --subset_group "${subset_group}"
  )
  append_cmd "${cmd[@]}"

  if [[ "${RENDER_ONLY}" != "true" ]]; then
    echo ">>> Running single analysis: ${sample}"
    if ! run_and_log "${log_file}" "${cmd[@]}"; then
      echo ">>> Single analysis failed: ${sample}. See log: ${log_file}" >&2
      record_failure "single" "${sample}" "${log_file}"
      continue
    fi
  fi
done

echo "sample1,sample2,outdir" > "${COMPARE_MANIFEST}"

declare -a COMPARE_ROWS=()
while IFS= read -r line; do
  [[ -z "${line}" ]] && continue
  COMPARE_ROWS+=("${line}")
done < <(parse_compares || true)

for row in "${COMPARE_ROWS[@]}"; do
  IFS=$'\t' read -r s1 s2 <<< "${row}"

  if [[ -z "${SAMPLE_EXISTS[${s1}]-}" || -z "${SAMPLE_EXISTS[${s2}]-}" ]]; then
    echo ">>> Skipping compare ${s1} vs ${s2}: sample not found in sample_list.csv"
    continue
  fi

  if [[ ! -f "${OUTDIR}/single/${s1}/cellchat_object.rds" || ! -f "${OUTDIR}/single/${s2}/cellchat_object.rds" ]]; then
    echo ">>> Skipping compare ${s1} vs ${s2}: required single analysis result is missing"
    continue
  fi

  compare_name="${s1}__vs__${s2}"
  compare_outdir="${OUTDIR}/compare/${compare_name}"
  mkdir -p "${compare_outdir}"
  log_file="${LOG_DIR}/compare_${compare_name}.log"
  input_ctrl="${OUTDIR}/single/${s1}/cellchat_object.rds"
  input_treat="${OUTDIR}/single/${s2}/cellchat_object.rds"

  printf '"%s","%s","%s"\n' \
    "${s1}" "${s2}" "$(resolve_abs "${compare_outdir}")" >> "${COMPARE_MANIFEST}"

  cmd=(
    "${CONDA_BIN}" run -p "${CONDA_ENV_ROOT}" --no-capture-output Rscript "${COMPARE_SCRIPT}"
    --input_ctrl "${input_ctrl}"
    --input_treat "${input_treat}"
    --name_ctrl "${s1}"
    --name_treat "${s2}"
    -o "${compare_outdir}"
    --python_path "${CONDA_ENV_ROOT}/bin/python"
  )
  append_cmd "${cmd[@]}"

  if [[ "${RENDER_ONLY}" != "true" ]]; then
    echo ">>> Running compare analysis: ${s1} vs ${s2}"
    if ! run_and_log "${log_file}" "${cmd[@]}"; then
      echo ">>> Compare analysis failed: ${s1} vs ${s2}. See log: ${log_file}" >&2
      record_failure "compare" "${compare_name}" "${log_file}"
      continue
    fi
  fi
done

report_cmd=(
  "${CONDA_BIN}" run -p "${CONDA_ENV_ROOT}" --no-capture-output Rscript "${REPORT_SCRIPT}"
  --input_dir "${OUTDIR}"
  --report_rmd "${REPORT_RMD}"
)
append_cmd "${report_cmd[@]}"

echo ">>> Rendering final report..."
if ! run_and_log "${LOG_DIR}/report_render.log" "${report_cmd[@]}"; then
  echo ">>> Report rendering failed. See log: ${LOG_DIR}/report_render.log" >&2
  record_failure "report" "CellChat_Report" "${LOG_DIR}/report_render.log"
fi

chmod +x "${PIPELINE_SCRIPT}"
echo ">>> Command script written to: $(resolve_abs "${PIPELINE_SCRIPT}")"
if [[ $(wc -l < "${FAILED_TASKS_FILE}") -gt 1 ]]; then
  echo ">>> Some tasks failed. See: $(resolve_abs "${FAILED_TASKS_FILE}")" >&2
fi
echo ">>> Pipeline complete."
