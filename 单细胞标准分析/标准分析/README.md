# 单细胞标准分析

这里是 `tbd` 中可直接复用的正式单细胞分析代码。核心 CPU 流程、命令行入口、样本发现、GPU 扩展和 R/Seurat 公共函数已经集中到同一目录；C4 专用版本与 GPU 亚群分析仍单独保留。

## 文件入口

| 文件 | 用途 |
|---|---|
| `scanpy_pipeline.py` | Scanpy CPU 主流程与可复用函数 |
| `run_scanpy_pipeline.py` | CPU 命令行入口，可运行单步或完整流程 |
| `sample_discovery.py` | 自动发现 10x 样本、推断分组并生成 `sample_sheet` |
| `scanpy_pipeline_gpu.py` | 标准流程的 GPU 加速实现及 CPU 回退 |
| `run_gpu_pipeline_example.py` | GPU 聚类与 marker 分析入口 |
| `scanpy_fun_subcluster_gpu.py` | GPU 亚群聚类函数 |
| `gpu_utils.py` | GPU 环境检测与 CPU/GPU 转换工具 |
| `C4_scanpy_pipeline.py` | 保留的 C4 专用流程 |
| `R single cell/00_config.R` | R/Seurat 参数与输出目录配置 |
| `R single cell/00_functions.R` | R/Seurat 可复用函数，加载时会自动读取配置 |
| `R single cell/sample_discovery.R` | R 版 10x 样本发现与分组 |
| `R single cell/main/*.R` | QC、聚类、marker、注释和比例的简短正式入口 |
| `R single cell/run_main.sh` | R 标准流程分步运行入口 |
| `R single cell/subcluster/*.R` | 亚群聚类、marker、注释和比例的简短正式入口 |
| `R single cell/run_subcluster.sh` | R 亚群流程分步运行入口 |

## Python/Scanpy

先生成样本表，再运行正式流程：

```bash
python sample_discovery.py discover --data-root data --output-manifest config/sample_manifest.tsv
python sample_discovery.py assign --manifest config/sample_manifest.tsv --group-count 2 --output-sample-sheet config/sample_sheet.tsv --output-report config/group_suggestions.tsv
python run_scanpy_pipeline.py all --sample-sheet config/sample_sheet.tsv --out-root results
```

`run_scanpy_pipeline.py` 也支持 `load`、`qc`、`filter`、`cluster`、`marker` 和 `annotate` 单步运行。项目正式运行时，建议在项目 `workflow/` 中新建简短的 `stepXX_name.py` 调用这些函数，不要复制长函数体。

GPU 聚类和 marker 分析：

```bash
python run_gpu_pipeline_example.py --check-env
python run_gpu_pipeline_example.py --adata results/03-Cell_filter/adata_counts.h5ad --out-root results
```

## R/Seurat

`00_functions.R` 会自动定位并加载同目录的 `00_config.R`。新项目只需在 `workflow/stepXX.分析点.R` 中加载函数并调用需要的分析点：

```r
source("R single cell/00_functions.R")
```

项目路径、分组、阈值和聚类参数通过 `00_config.R` 或环境变量设置；已运行过的 workflow 脚本不覆盖，参数调整时新增下一个编号步骤。

完整 R 流程也可以分阶段运行：

```bash
bash "R single cell/run_main.sh" preflight
bash "R single cell/run_main.sh" qc1
bash "R single cell/run_main.sh" after_qc
```

`qc1` 会先生成样本 QC 阈值表；人工确认阈值和样本分组后，再运行 `after_qc`。其他可用模式包括 `qc2`、`cluster`、`markers`、`annotation` 和 `proportion`。

亚群流程使用 `00_config.R` 中的 `SUBCLUSTER` 参数和待分析细胞类型设置：

```bash
bash "R single cell/run_subcluster.sh" all
```

也可使用 `cluster`、`markers`、`annotation` 或 `proportion` 单独运行一个亚群步骤。正式项目中参数调整或重跑时，仍应在项目 `workflow/` 新增下一个编号脚本。

原 `03_findmarkers.R` 中的胰腺癌小鼠 marker、固定对象路径和单 PNG 绘图属于具体项目，原稿保存在 `legacy/original/01_qc_clustering_annotation/03_findmarkers.R`。带服务器内部路径的 Streamlit 比例分析包装器也只保存在同一 legacy 目录，不进入正式流程。

## 输出目录

Python 正式流程在指定的 `results/` 根目录下输出：

```text
results/
├── 01-CellRanger/
├── 02-Data_QC/
├── 03-Cell_filter/
├── 04-Dimension_reduction_cluster/
├── 05-Marker_gene_analysis/
└── 06-Cell_type_annotation/
```

R/Seurat 默认使用 `01_QC/`、`02_GlobalClustering/`、`03_FindMarkers/`、`04_CellAnnotation/`、`05_CellProportion/` 和 `06_Subcluster/`，中间 `.qs2` 对象保存在配置的对象目录。正式项目中建议把这些编号目录统一放在 `results/` 下，把中间对象放到 `data/processed/`。
