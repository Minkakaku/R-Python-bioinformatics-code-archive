# R-Python-bioinformatics-code-archive

本目录汇总服务器上的单细胞及其他生物信息分析代码，作为新项目的代码素材和来源归档使用。所有文件均为复制件，原始文件没有移动或修改。

## 推荐使用方式

不要直接把整个目录当作一个固定流水线运行。新项目开始时，从对应专题目录挑选需要的脚本，复制到新项目的 `workflow/`，再按实际执行顺序命名，例如：

```text
step1.qc.R
step2.global_cluster.R
step3.annotation.R
```

正式单细胞分析入口：

- R/Seurat：`单细胞标准分析/标准分析/R single cell/`
- Python/Scanpy：`单细胞标准分析/标准分析/scanpy_pipeline.py`
- Python 命令行：`单细胞标准分析/标准分析/run_scanpy_pipeline.py`
- 样本发现：`单细胞标准分析/标准分析/sample_discovery.py`
- GPU 扩展、C4 专用版和完整使用方式：`单细胞标准分析/标准分析/README.md`

轨迹及 05–09 专题已提取为 `savebio` 下的同级模板目录。原 QC、全局分析和亚群分析入口均已并入正式 R/Python 单细胞模板；CellChat 已独立为单独分析点。

专题分析入口：

- CellChat：独立分析点目录 `../R-cellchat-workflow-template/`，一个完整 R 文件覆盖单样本、双样本和多样本；原始来源保存在 `legacy/original/03_cellchat/`
- 轨迹/拟时序/RNA velocity：独立目录 `../R-Python-trajectory-workflow-template/`，含 Monocle2、Monocle3、Slingshot、tradeSeq、PAGA、scVelo 及 343 份来源代码索引
- 轨迹原始平台和项目脚本：`legacy/original/trajectory_source_archive/`，仅用于归档仓库内溯源
- CellPhoneDB：`../Python-cellphonedb-workflow-template/`
- inferCNV / CopyKAT：`../R-Python-infercnv-copykat-workflow-template/`
- specific infercnvpy notebook：已转换为正式 `.py`；原 notebook 保存在 `legacy/original/specific/`
- SCENIC / pySCENIC：`../R-Python-scenic-regulon-workflow-template/`
- 打分、富集和可视化：`../R-Python-scores-enrichment-workflow-template/`
- Visium HD：`../R-visium-hd-workflow-template/`，正式 R/Seurat 入口，默认读取 Space Ranger 8 µm bins；混合空间来源代码保留在其 `source_catalog/`

## 目录说明

```text
单细胞标准分析/标准分析/          正式 R/Seurat、Python/Scanpy、GPU 与 C4 分析入口
../R-Python-trajectory-workflow-template/     独立轨迹模板：正式入口、全量来源代码与溯源索引
../Python-cellphonedb-workflow-template/    CellPhoneDB 及可视化代码素材
../R-Python-infercnv-copykat-workflow-template/ inferCNV 和 CopyKAT 代码素材
../R-Python-scenic-regulon-workflow-template/ SCENIC / pySCENIC 代码素材
../R-Python-scores-enrichment-workflow-template/ 打分、富集与可视化代码素材
../R-visium-hd-workflow-template/ R/Seurat Visium HD 正式模板及混合空间来源代码
10_other_users/                   其他人员原始哈希归档及独立模板分配说明
11_semantic_catalog/              按用途区分同名代码的语义整理层
legacy/original/                  原样保留的项目特异或平台代码
legacy_versions/                  内容不同但不建议作为首选的旧版本
```

## 同名代码怎么选

不要按文件名直接去重。历史项目里常有多份 `CellChat.R`、`monocle.R` 或 `report.rmd`，名称相同但分别用于单数据集、双组比较、多组比较、绘图或批量运行。

- `11_semantic_catalog/code_catalog.tsv` 提供全库场景索引。
- `11_semantic_catalog/duplicate_original_names.tsv` 集中列出同名的独立版本。
- `../R-Python-trajectory-workflow-template/SOURCE_CATALOG.tsv` 单独索引 343 份轨迹相关唯一代码，并区分 Monocle2/3、BEAM、Slingshot、tradeSeq、PAGA、RNA velocity 和 CytoTRACE。
- `11_semantic_catalog/03_cellchat/` 已把 CellChat 代码按单数据集、双组、多组、组数待确认、比较后处理、绘图和批处理拆开，并使用来源人员与 SHA-256 生成唯一文件名。
- `10_other_users/code/` 始终保留原始归档；语义整理层只建立逐字节副本和索引。
- `10_other_users/README.md` 记录 11 个细分类别到各独立模板的完整分配关系。

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
- `../R-Python-trajectory-workflow-template/` 已提供 Monocle2、Monocle3、Slingshot、tradeSeq、PAGA 和 scVelo 正式入口；历史来源脚本仍需按当前包版本小规模试运行。

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
