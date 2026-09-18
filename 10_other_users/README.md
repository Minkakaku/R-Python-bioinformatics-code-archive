# 其他分析人员的单细胞代码归档

本目录只整理 `/home/hf` 以外的代码，来源主要为：

- `/mnt/NexData_file/home/<分析人员>/`
- `/data/lvdata/sc-data/<分析人员或项目>/`

源文件只读，未移动、未修改。归档文件按 SHA-256 去重；同一代码在多个项目或人员目录中出现时只复制一份，但所有来源均记录在 `manifest.tsv`。

## 目录结构

```text
10_other_users/
├── README.md
├── manifest.tsv          # 所有保留来源、哈希、分类、目标文件和重复状态
├── recommended.tsv       # 排除纯项目特异脚本后的优先浏览清单
├── excluded.tsv          # 第三方、非单细胞和无权限文件
├── summary.tsv           # 按专题、人员和复用等级统计
├── CHECKSUMS.sha256      # 归档代码完整性校验
└── code/
    ├── 03_cellchat/
    ├── 04_cellphonedb/
    ├── 05_trajectory_monocle/
    ├── 06_trajectory_other/
    ├── 07_infercnv_copykat/
    ├── 08_scenic_regulon/
    ├── 09_qc_doublet/
    ├── 10_scanpy/
    ├── 11_clustering_annotation/
    ├── 12_scores_enrichment/
    ├── 13_spatial_xenium/
    └── 15_core_scrna/
```

归档文件名采用：

```text
<SHA256前12位>__<原文件名>
```

这能避免不同人员使用相同文件名时互相覆盖。原始项目路径、人员、完整哈希和修改时间均可从 `manifest.tsv` 追溯。

## 检索方法

扫描同时使用两类条件：

1. 路径和文件名关键词，例如 CellChat、CellPhoneDB、Monocle、SCENIC、inferCNV、CopyKAT、Scanpy、subcluster、Xenium、Visium、spatial。
2. 代码内容特征，例如 `CreateSeuratObject`、`FindClusters`、`createCellChat`、`newCellDataSet`、`learn_graph`、`DoubletFinder`、`RunHarmony`、`sc.pp`、`sc.tl`、`AddModuleScore`。

仅收集非空且小于 2 MB 的 `.R/.Rmd/.py/.sh`。没有收集 notebook、数据、图片、HTML、Cell Ranger 输出、原始测序文件或分析环境。

## 复用等级

- `reusable_core`：名称和用途较明确的主流程、函数或专题代码。
- `runner`：运行命令和 shell 启动脚本。
- `report_template`：R Markdown 报告模板。
- `supporting_code`：绘图、转换、拆分、输出等辅助代码。
- `project_specific`：强项目依赖的分析脚本，适合参考而不适合直接运行。

`recommended.tsv` 收录除 `project_specific` 外的归档版本，可先从这里筛选。

## 重要发现

- 找到了 `cluster_gpu.py`，来源为 ztt 的 `spatialRNA_gradio` 聚类函数目录。
- 找到了多位分析人员修改过的 Monocle3 脚本，而不是只有 Monocle2。
- CellChat、Scanpy 亚群、SCENIC、inferCNV、CopyKAT、CellPhoneDB 和空间分析存在大量项目版本；重复脚本已合并，独立修改版本均保留。

## 排除规则

以下内容会登记到 `excluded.tsv`，但不复制：

- Conda、R library、Docker R 环境、软件安装目录和第三方源码。
- Spaceranger、Xenium Ranger、RAPIDS-singlecell、LIGER、CytoTRACE2 等软件包自身源码。
- 路径明确属于 bulk RNA-seq、蛋白、代谢、基因组等项目，且没有单细胞特征的脚本。
- 无读取权限的文件。

## 使用注意

- 代码来自历史项目，可能含硬编码路径、样本名、细胞类型、数据库位置和软件版本。
- 文件按原内容复制，没有自动重写，也没有使用真实数据完整复跑。
- 同一专题的多个独立哈希代表真实的不同版本，不等同于都适合作为模板；优先查看 `recommended.tsv`、修改时间和原始来源。
- 若准备上传 GitHub，建议从 `recommended.tsv` 再挑选一轮，避免公开项目路径或内部样本信息。

## 常用查询

```bash
cd R-Python-bioinformatics-code-archive/10_other_users

# 查找 GPU 聚类脚本
grep -i 'cluster_gpu.py' manifest.tsv

# 查看所有 Monocle3 来源
grep -i 'monocle3' manifest.tsv

# 查看某位分析人员的独立归档代码
awk -F '\t' '$2 == "sh" && $12 == "copied"' manifest.tsv

# 校验归档代码
sha256sum -c CHECKSUMS.sha256
```


## 整理结果

- 独立归档代码：3657 个
- 来源记录：5531 条
- 推荐浏览记录：2659 条
- 排除记录：2847 条
- 收集期间变更：6 条
- 校验通过：3657/3657
- 目录总大小：52957725 bytes
- 整理日期：2026-07-21

## 独立模板分配

原始 owner/hash 归档继续保留在本目录；以下模板接收的是 `11_semantic_catalog/` 中逐字节相同的语义副本：

| 类别 | 原始文件数 | 模板去向 |
|---|---:|---|
| `03_cellchat` | 381 | 保留在归档仓库；CellChat 正式仓库按既定要求只放一份覆盖单/双/多样本的完整 R 代码 |
| `04_cellphonedb` | 28 | `../../Python-cellphonedb-workflow-template/source_catalog/` |
| `05_trajectory_monocle` | 281 | `../../R-Python-trajectory-workflow-template/source_catalog/` |
| `06_trajectory_other` | 21 | `../../R-Python-trajectory-workflow-template/source_catalog/` |
| `07_infercnv_copykat` | 78 | `../../R-Python-infercnv-copykat-workflow-template/source_catalog/` |
| `08_scenic_regulon` | 62 | `../../R-Python-scenic-regulon-workflow-template/source_catalog/` |
| `09_qc_doublet` | 4 | `../../R-scrna-workflow-template/source_catalog/` |
| `10_scanpy` | 246 | `../../Python-scrna-workflow-template/source_catalog/` |
| `11_clustering_annotation` | 754 | 按语言分到 R/Python 单细胞模板；Shell 在两边共享 |
| `12_scores_enrichment` | 294 | `../../R-Python-scores-enrichment-workflow-template/source_catalog/` |
| `13_spatial_xenium` | 883 | Python/Shell 到 Python Xenium；R/Rmd/Shell 作为来源素材进入 R Visium HD 模板的 `source_catalog/`；Shell 两边共享 |
| `15_core_scrna` | 625 | 按语言分到 R/Python 单细胞模板；Shell 在两边共享 |

每个模板根目录的 `SOURCE_CATALOG.tsv` 记录 method、scenario、workflow role、owner、原始路径、语义副本、SHA-256 和目标文件。Shell 脚本在跨语言流程中有意复制到两个模板，因此八个模板共有 3068 个副本，对应 2974 个唯一 SHA-256。
