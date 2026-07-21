# ==================== 拟时序分析 - 方法3：Seurat 差异表达基因 ====================
library(Seurat)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(dplyr)
library(hdf5r)
library(patchwork)
library(RColorBrewer)
library(monocle)
library(tibble)

# 创建输出目录
monocle_outdir <- "Monocle_Epi_method3"
dir.create(monocle_outdir, showWarnings = FALSE)
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead",
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c",
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1",
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E",
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d",
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9",
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6",
          "#005995", "#fa625f", "#f17d80")
# 读取 Seurat 对象
cds <- readRDS("cds_deg-3.rds")
cds <- orderCells(cds,root_state=6)

# 绘制轨迹图
p <- plot_cell_trajectory(cds, color_by = "annotation")
pdf(file.path(monocle_outdir, "monocle_cluster-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "annotation") + 
  facet_wrap(~annotation, nrow = 2) +
  scale_color_manual(values = plan)
pdf(file.path(monocle_outdir, "monocle_facet-annotation-3.pdf"), height = 7, width = 12)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

# 保存 cds 对象
saveRDS(cds, file.path(monocle_outdir, "cds_deg-3_new.rds"))

pdf(file.path(monocle_outdir, "cell trajectory.pdf"), width = 8, height = 6)
plot_cell_trajectory(cds,color_by = "Pseudotime", size = 1, show_backbone = TRUE) + theme(legend.position = "right")
plot_cell_trajectory(cds,color_by = "annotation", size = 1, show_backbone = TRUE) + 
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
dev.off()


pdf(file.path(monocle_outdir, "cell trajectory - - split.pdf"), width = 18, height = 9)
plot_cell_trajectory(cds, color_by = "annotation") + facet_wrap("~annotation", nrow = 2) + theme(legend.position = "right") + scale_color_manual(values = plan) 
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 2) + scale_color_manual(values = plan)
dev.off()


pdf(file.path(monocle_outdir, "cell trajectory - annotation - split by group.pdf"), width = 16, height = 5)
local({
  pData(cds)$group <- pData(cds)$group %>% factor(levels = c("CON", "IRE_4", "IRE_7"))
  plot_cell_trajectory(cds, color_by = "annotation", size=1, show_backbone=TRUE) +
    scale_color_manual(values = plan) + 
    facet_wrap(~ group) +
    theme(legend.position = "right")
})
dev.off()



pdf(file.path(monocle_outdir, "cell trajectory - tree.pdf"), width = 15, height = 12)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation") +
  scale_color_manual(values = plan) +
  facet_wrap(~ annotation)
dev.off()

# cell barplot 
p_col_cell_state <-
  cds %>% 
  pData() %>%
  dplyr::select(annotation, State) %>%
  group_by(annotation, State) %>%
  mutate(annotation = as.character(annotation)) %>%
  # mutate(annotation = factor(annotation, levels = c(
  #   'M0',
  #   'M1',
  #   'M2'
  # ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation, y = freq)) +
  geom_col(aes(fill = State)) +
  scale_fill_manual(values = plan) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(monocle_outdir, "cell stat barplot.pdf"), height = 6, width = 8)


cat("方法3（差异表达基因）运行完成，结果保存在：", monocle_outdir, "\n")
