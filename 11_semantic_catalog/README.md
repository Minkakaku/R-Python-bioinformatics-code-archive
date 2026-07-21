# 语义整理层

这里是 `10_other_users/code/` 的可检索整理层。原始归档文件不会被移动、改名或修改；因此来源记录和 SHA-256 校验仍然有效。

## 文件

- `code_catalog.tsv`：3,657 份独立代码的统一索引，增加分析场景、工作流角色、数据类型和语义文件名。
- `duplicate_original_names.tsv`：只列出原始文件名重复的版本。相同名字但不同 SHA-256 的代码会分别保留，不按文件名删除。
- `by_analysis/`：全库按专题和分析场景整理的代码视图，每份文件都有唯一语义名称。
- `03_cellchat/`：CellChat 的精细场景视图，额外区分比较组数和单细胞/空间来源。

## CellChat 分类

```text
01_single_dataset/             单个数据集独立推断细胞通讯
02_two_group_comparison/       两组端到端分析并合并比较
03_multi_group_comparison/     三组及以上端到端分析并合并比较
04_comparison_postprocessing/  读取已有 CellChat 对象后做双组或多组比较
05_visualization/              已有 CellChat 结果的绘图
06_batch_runner/               Shell/Python 批量启动脚本
07_report_template/            R Markdown 报告模板
08_preprocessing_support/      数据准备、转换和辅助代码
09_group_comparison_unspecified/ 已确认是组间比较，但静态代码未明确组数
99_unclear/                    静态规则无法可靠判定的版本
```

每个场景再分为 `spatial_source` 和 `scrna_or_unspecified`。文件名采用：

```text
<数据类型>__<来源人员>__<原始文件名>__<SHA256前12位>.<扩展名>
```

例如，原来都叫 `CellChat.R` 的代码会成为不同的语义文件名，并且由所在目录明确区分单数据集、双组比较和多组比较。

## 使用建议

1. 先看 `03_cellchat/summary.tsv`，按场景缩小范围。
2. 再看 `03_cellchat/catalog.tsv` 的 `reuse_level`、`source_path` 和功能标记。
3. 选择代码后复制到新项目的 `workflow/`，按运行顺序改成 `stepN.<用途>.R`。
4. `99_unclear` 和 `project_specific` 只作为历史参考，不建议直接运行。

## 重新生成

```bash
python tools/build_semantic_catalog.py
```

该命令只重建 `11_semantic_catalog/03_cellchat/` 及 TSV 索引，不会修改 `10_other_users/` 的原始归档。
