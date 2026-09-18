# CellChat 场景化代码视图

本目录用于解决“脚本同名、实际用途不同”的问题。代码按行为分类，而不是按原文件名分类。

- `catalog.tsv`：每个版本的分类、原始名称、来源、SHA-256、复用等级和功能标记。
- `summary.tsv`：各场景和数据类型的文件数。
- 场景目录：原始代码的逐字节副本，文件名包含来源人员和哈希，不会相互覆盖。

静态分类主要检查 `createCellChat()`、`computeCommunProb()`、`mergeCellChat()`、`compareInteractions()`、`netVisual_diffInteraction()` 和 `rankNet(mode = "comparison")` 等行为。只有明确检测到两个对象/分组时才标记为双组，明确检测到三组及以上时才标记为多组；无法从静态代码确定组数的比较流程放在 `09_group_comparison_unspecified/`。

这层没有自动修改历史代码中的项目路径、样本名或参数。
