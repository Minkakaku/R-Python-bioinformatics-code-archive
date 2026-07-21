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
kegg.enrich <- read_excel("up.xlsx")

# 2. 数据整理：不筛选，表格里所有通路都画
kegg.enrich <- kegg.enrich %>%
  mutate(
    pvalue = as.numeric(gsub("<", "", trimws(as.character(pvalue)))),
    Count = as.numeric(trimws(as.character(Count))),
    logP = -log10(pvalue)
  ) %>%
  arrange(pvalue)

# 3. 固定通路顺序
kegg.enrich$Description <- factor(
  kegg.enrich$Description,
  levels = rev(kegg.enrich$Description)
)

# 4. 计算双坐标轴比例
ratio <- max(kegg.enrich$Count, na.rm = TRUE) /
  max(kegg.enrich$logP, na.rm = TRUE)

kegg.enrich <- kegg.enrich %>%
  mutate(
    Count_scaled = Count / ratio
  ) %>%
  arrange(Description)

# 5. 作图
p <- ggplot(kegg.enrich, aes(y = Description)) +
  geom_col(
    aes(x = logP),
    fill = "#1f73a7",
    width = 0.7
  ) +
  geom_point(
    aes(x = Count_scaled),
    color = "#f07f30",
    size = 2
  ) +
  geom_path(
    aes(x = Count_scaled, group = 1),
    color = "#f07f30",
    linewidth = 1
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.12)),
    name = "-log10(pvalue)",
    sec.axis = sec_axis(
      ~ . * ratio,
      name = NULL
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    y = "",
    title = "KEGG Pathway Enrichment"
  )

# 6. 显示图片
p

# 7. 保存图片
ggsave("KEGG Pathway Enrichment.pdf", p, width = 9, height = 6)
ggsave("KEGG Pathway Enrichment.png", p, width = 9, height = 6, dpi = 300)