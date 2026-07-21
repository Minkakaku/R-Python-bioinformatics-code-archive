# 单细胞分析代码整理库（tbd）

本目录汇总服务器上现有的单细胞分析脚本，目标是作为新项目的代码素材库使用。所有文件均为复制件，原始文件没有移动或修改。

## 推荐使用方式

不要直接把整个目录当作一个固定流水线运行。新项目开始时，从对应专题目录挑选需要的脚本，复制到新项目的 `workflow/`，再按实际执行顺序命名，例如：

```text
step1.qc.R
step2.global_cluster.R
step3.annotation.R
step4.cellchat.R
```

R/Seurat 项目可先参考：

- `00_core_scrna/00_config.R`
- `00_core_scrna/00_functions.R`
- `01_qc_clustering_annotation/`
- `02_subcluster/`

Python/Scanpy 项目可先参考：

- `00_core_scrna/scanpy_pipeline.py`
- `00_core_scrna/run_scanpy_pipeline.py`
- `00_core_scrna/sample_discovery.py`

专题分析入口：

- CellChat：`03_cellchat/core/pip_cellchat.R`
- Monocle2：`04_trajectory/monocle2/pseudotime.R` 或 `platform_monocle2.analysis.R`
- CellPhoneDB：`05_cellphonedb/`
- inferCNV / CopyKAT：`06_infercnv_copykat/`
- SCENIC / pySCENIC：`07_scenic_regulon/`
- 打分、富集和可视化：`08_scores_enrichment/`

## 目录说明

```text
00_core_scrna/                    通用配置、函数、Scanpy 流程和样本发现
01_qc_clustering_annotation/      QC、整合/聚类、marker、注释、比例统计
02_subcluster/                    亚群聚类、marker、注释和比例统计
03_cellchat/                      CellChat 主流程、报告模板及项目示例
04_trajectory/                    Monocle2 和 PAGA 轨迹分析
05_cellphonedb/                   CellPhoneDB 及可视化辅助脚本
06_infercnv_copykat/              inferCNV 和 CopyKAT
07_scenic_regulon/                SCENIC / pySCENIC
08_scores_enrichment/             score、dotplot、富集及相关报告模板
09_spatial_related/               与 Xenium/空间转录组直接相关的可复用脚本
legacy_versions/                  内容不同但不建议作为首选的旧版本
```

## 去重与溯源

- `manifest.tsv` 记录分类、目标文件、原始路径、SHA-256、大小、修改时间、状态和备注。
- `status=copied` 表示实际复制的文件。
- `status=duplicate_source` 表示该来源与已经保留的目标文件哈希完全一致，因此没有再次复制。
- `status=excluded_empty` 表示发现了脚本名，但文件为零字节，没有纳入代码库。
- `CHECKSUMS.sha256` 可用于核验整理后的代码是否被改动。

## 收集范围

本次主要定向检索了：

- `/home/hf/DataFile`
- `/mnt/NexData_file/share/scRNA_platform/scVision5/moudle/tools`
- `/mnt/NexData_file/share/scRNA-tools` 中 CellChat 主流程及报告模板
- `/data/lvdata/sc-data/wyt/cellchat`
- 已定位的 inferCNV、CellChat、Monocle、SCENIC 等专题项目目录

服务器包含约 33 TB 的本地项目盘和约 PB 级挂载存储，因此没有对原始数据盘、Cell Ranger 输出和所有归档结果做无边界遍历。这样可以避免扫描海量 `.rds/.h5ad/.csv`、图片、网页报告和环境文件，也不会把第三方包源码混进本目录。

## 使用前注意

- 这些代码保留原始实现，部分脚本包含项目路径、样本名、分组、细胞类型列名或数据库路径，使用前必须按新项目修改。
- `platform_wrappers/` 多为 Streamlit/平台启动器；真正的分析代码通常在同专题的 R/Python 主脚本中。
- CellPhoneDB 的历史脚本含删除既有输出目录的逻辑，运行前请先确认输出路径。
- CellChat、Monocle2、Seurat、CellPhoneDB、inferCNV、SCENIC 的接口会随版本变化，需在当前环境中小规模试运行。
- 本次只做代码收集、分类、哈希去重和静态整理，没有使用真实数据完整复跑。
- 没有在已检索范围找到 `cluster_gpu.py`；因此本目录未包含该文件。
- 没有找到可确认的独立 Monocle3 正式脚本；当前轨迹分析以 Monocle2 和 PAGA 为主。

## 明确排除

- `.rds`、`.h5ad`、表达矩阵和项目结果表
- PNG/PDF/HTML 等结果文件
- Conda/虚拟环境、数据库和第三方包源码
- Jupyter Notebook（按脚本优先原则不收集）
- 零字节的 `cellchat.R`、`sc_boxplot.py`、`sc_volcano.py`


## 本次整理统计

- 实际复制代码：87 个
- 哈希重复来源：12 个
- 排除零字节脚本：3 个
- 代码总大小：623596 bytes
- 整理日期：2026-07-20
