# 1. 加载必要的依赖包
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

# 假设你的 Seurat 对象名为 seurat_obj
# 假设亚群储存在 active.ident 或 metadata 的 "cell_type" 中
# 假设分组（R和N）储存在 metadata 的 "group" 中
# 请根据实际情况修改下方字段名：
Idents(seuset)
cell_type_col <- "leiden" # 如果是当前 identity，如果是列名则写如 "cell_type"
group_col     <- "group"        # 分组字段名（包含 R 和 N）
score_col     <- "cnv_score_cell" # 得分字段名

---
  
  # 2. 提取并准备数据
  # 提取 metadata 转化为 Data Frame
  if (cell_type_col == "active.ident") {
    plot_data <- data.frame(
      cell_type = Idents(seuset),
      group     = seuset[[group_col]][,1],
      cnv_score = seuset[[score_col]][,1]
    )
  } else {
    plot_data <- data.frame(
      cell_type = seuset[[cell_type_col]][,1],
      group     = seuset[[group_col]][,1],
      cnv_score = seuset[[score_col]][,1]
    )
  }

# 过滤掉既不是 R 结构也不是 N 结构的数据（防止有脏数据）
plot_data <- plot_data %>% filter(group %in% c("R", "N"))

  # 3. 统计计算与 P 值导出 (Wilcoxon 秩和检验)
  # 对 10 个亚群分别进行 R vs N 的差异检验
stat_results <- plot_data %>%
group_by(cell_type) %>%
summarise(
  Mean_R = mean(cnv_score[group == "R"], na.rm = TRUE),
  Mean_N = mean(cnv_score[group == "N"], na.rm = TRUE),
  Median_R = median(cnv_score[group == "R"], na.rm = TRUE),
  Median_N = median(cnv_score[group == "N"], na.rm = TRUE),
  # 使用 Wilcoxon 检验计算 P 值
  P_Value = wilcox.test(cnv_score ~ group)$p.value,
  .groups = 'drop'
) %>%
# 针对多重假设检验进行 P 值校正（FDR/BH法）
mutate(FDR = p.adjust(P_Value, method = "BH"))

# 导出统计结果为 CSV 文件
write.csv(stat_results, file = "../inferCNV_score_statistics_R_vs_N.csv", row.names = FALSE)
message("统计报表已成功导出至：inferCNV_score_statistics_R_vs_N.csv")

  
# 4. 绘制小提琴图并自动添加 P 值
# 使用 ggpubr 可以非常方便地在分面图上添加组间统计 P 值
p <- ggplot(plot_data, aes(x = group, y = cnv_score, fill = group)) +
  # 绘制原生小提琴图
  geom_violin(alpha = 0.8, trim = FALSE) + 
  # 调整颜色：R红 N蓝
  scale_fill_manual(values = c("R" = "#E41A1C", "N" = "#377EB8")) +
  # 分面展示 10 个亚群
  facet_wrap(~ cell_type, scales = "free_y", ncol = 5) + 
  # 添加统计学检验
  stat_compare_means(comparisons = list(c("R", "N")), 
                     method = "wilcox.test", 
                     label = "p.signif", 
                     label.y.npc = 0.85) + 
  # 皮肤与标签
  theme_bw() +
  labs(title = "inferCNV Score Comparison (R vs N)",
       x = "Group",
       y = "inferCNV Score") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(), # 移除网格线让图更美观
        panel.grid.minor = element_blank())

# 预览图片
print(p)

# 保存图片到本地
ggsave("../inferCNV_score_violin_leiden.pdf", plot = p, width = 12, height = 8)
print(p)
save.image("../monocleDT_2.Rdata")
