# ================================
# KEGG enrichment bar-line plot
# Input : all.xlsx
# Output: KEGG_enrichment_barline.pdf / png
# ================================

library(readxl)
library(dplyr)
library(ggplot2)

setwd("/mnt/NexData_file/home/cjy/work/move/Trance/售后/转录组售后/T260435")

# 1. 读取数据
kegg.enrich <- read_excel("all.xlsx")

# 2. 数据整理：不筛选，保留 Excel 原始顺序
kegg.enrich <- kegg.enrich %>%
  mutate(
    Description = as.character(Description),
    pvalue = as.numeric(gsub("<", "", trimws(as.character(pvalue)))),
    Count = as.numeric(trimws(as.character(Count))),
    logP = -log10(pvalue)
  )

# 3. 固定通路顺序：Excel 第一行显示在图最上方
kegg.enrich$Description <- factor(
  kegg.enrich$Description,
  levels = rev(kegg.enrich$Description)
)

# 4. 计算 Count 与 logP 的缩放比例
ratio <- max(kegg.enrich$Count, na.rm = TRUE) /
  max(kegg.enrich$logP, na.rm = TRUE)

kegg.enrich <- kegg.enrich %>%
  mutate(
    Count_scaled = Count / ratio
  )

# 5. 为右侧图例预留空间
x_max <- max(kegg.enrich$logP, kegg.enrich$Count_scaled, na.rm = TRUE)
x_limit <- x_max * 1.45

legend_x <- x_max * 1.08
legend_y1 <- nrow(kegg.enrich) - 0.5
legend_y2 <- nrow(kegg.enrich) - 1.4

# 6. 作图
p <- ggplot(kegg.enrich, aes(y = Description)) +
  geom_col(
    aes(x = logP),
    fill = "#1F78B4",
    width = 0.7
  ) +
  geom_point(
    aes(x = Count_scaled),
    color = "#F8761D",
    size = 2
  ) +
  geom_path(
    aes(x = Count_scaled, group = 1),
    color = "#F8761D",
    linewidth = 0.7
  ) +
  annotate(
    "rect",
    xmin = legend_x,
    xmax = legend_x + x_max * 0.035,
    ymin = legend_y1 - 0.12,
    ymax = legend_y1 + 0.12,
    fill = "#1F78B4"
  ) +
  annotate(
    "text",
    x = legend_x + x_max * 0.05,
    y = legend_y1,
    label = "-log10(pvalue)",
    hjust = 0,
    size = 4
  ) +
  annotate(
    "point",
    x = legend_x + x_max * 0.018,
    y = legend_y2,
    color = "#F8761D",
    size = 2
  ) +
  annotate(
    "text",
    x = legend_x + x_max * 0.05,
    y = legend_y2,
    label = "Term candidate\ngene num",
    hjust = 0,
    size = 4
  ) +
  scale_x_continuous(
    limits = c(0, x_limit),
    expand = expansion(mult = c(0, 0)),
    name = "-log10(pvalue)",
    sec.axis = sec_axis(
      ~ . * ratio,
      name = NULL
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 14) +
  labs(
    y = "",
    title = "KEGG Pathway Enrichment"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(10, 40, 10, 10)
  )

# 7. 显示图片
p

# 8. 保存图片
ggsave("KEGG_enrichment_barline.pdf", p, width = 9, height = 6)
ggsave("KEGG_enrichment_barline.png", p, width = 9, height = 6, dpi = 300)