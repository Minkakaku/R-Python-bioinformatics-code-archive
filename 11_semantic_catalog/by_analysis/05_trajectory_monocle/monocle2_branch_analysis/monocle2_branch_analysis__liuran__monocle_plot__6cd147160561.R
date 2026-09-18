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


monocle_outdir <- "/mnt/NexData_file/home/liuran/projects/3_Xenium/20260113_FWSP20252695_Thyroid_Xenium/monocle_results_260407"


plan <- c(
  "#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
  "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead",
  "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c",
  "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1",
  "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E",
  "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d",
  "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9",
  "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6",
  "#005995", "#fa625f", "#f17d80"
)


cds <- readRDS(file.path(monocle_outdir, "cds_deg-3.rds"))
pData(cds)$RCTD_subtype2 <- factor(pData(cds)$RCTD_subtype2, levels = c(paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)))


# monocle plot -----------
### plot_cell_trajectory -------
pdf(file.path(monocle_outdir, "cell trajectory.pdf"), width = 8, height = 6)
plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE) + theme(legend.position = "right")
plot_cell_trajectory(cds, color_by = "RCTD_subtype2", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
dev.off()


pdf(file.path(monocle_outdir, "cell trajectory - - split by RCTD_subtype2.pdf"), width = 25, height = 25)
plot_cell_trajectory(cds, color_by = "RCTD_subtype2") + facet_wrap("~RCTD_subtype2", nrow = 5) + theme(legend.position = "right") + scale_color_manual(values = plan)
dev.off()

pdf(file.path(monocle_outdir, "cell trajectory - - split by state.pdf"), width = 15, height = 10)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 2) + scale_color_manual(values = plan)
dev.off()

pdf(file.path(monocle_outdir, "cell trajectory - annotation - split by group.pdf"), width = 12, height = 5)
local({
  pData(cds)$group <- pData(cds)$group %>% factor(levels = c("PA", "PC"))
  plot_cell_trajectory(cds, color_by = "RCTD_subtype2", size = 1, show_backbone = TRUE) +
    scale_color_manual(values = plan) +
    facet_wrap(~group) +
    theme(legend.position = "right")
})
dev.off()


### tree --------
pdf(file.path(monocle_outdir, "cell trajectory - tree.pdf"), width = 30, height = 30)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "RCTD_subtype2") + theme(legend.position = "right") +
  scale_color_manual(values = plan) +
  facet_wrap(~RCTD_subtype2)
dev.off()


### cell barplot --------
p_col_cell_state <-
  cds %>%
  pData() %>%
  dplyr::select(RCTD_subtype2, State) %>%
  group_by(RCTD_subtype2, State) %>%
  mutate(RCTD_subtype2 = as.character(RCTD_subtype2)) %>%
  mutate(RCTD_subtype2 = factor(RCTD_subtype2, levels = c(
    paste0("Normal_C", 0:6), paste0("Normal_C", 8:11), paste0("Tumor_C", 0:9)
  ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x = RCTD_subtype2, y = freq)) +
  geom_col(aes(fill = State)) +
  scale_fill_manual(values = plan) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(monocle_outdir, "cell stat barplot.pdf"), height = 12, width = 8)


### density plot -----------
colour <- c(
  "#C6DADF", "#E26EB8", "#918FBE", "#A2C5A6", "#EED784", "#D4B796",
  "#F5BDC8", "#F0BDD8", "#77A1EC", "#F8D39C", "#F17591", "#CCC9E6"
)

data <- as.data.frame(pData(cds))

p_density <- ggplot(data, aes(Pseudotime, fill = RCTD_subtype2)) +
  geom_density(alpha = 0.5, color = "black", size = 1) +
  theme_bw() +
  theme( # 进一步修改主题
    panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  ) +
  theme(legend.title = element_text(size = 20, color = "black")) + # 图例标题
  theme(legend.text = element_text(size = 12, color = "black")) + # 图例文字
  theme(
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  theme(
    axis.text = element_text(size = 20, color = "black"),
    axis.ticks = element_line(size = 1), # 轴刻度
    axis.line = element_line(size = 1)
  ) + # 轴线粗细
  theme(
    legend.position = "top",
    legend.direction = "horizontal"
  ) +
  theme(legend.key = element_blank()) +
  guides(fill = guide_legend(nrow = 3, byrow = T)) +
  scale_fill_manual(values = plan) +
  labs(fill = "") +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid")) # 边框粗细
ggsave(p_density, file = file.path(monocle_outdir, "celltype density plot.pdf"), height = 6, width = 12)
ggsave(p_density, file = file.path(monocle_outdir, "celltype density plot.png"), height = 6, width = 12, dpi = 300)


p_density_celltype <- ggplot(data, aes(Pseudotime, fill = RCTD_subtype2)) +
  geom_density(alpha = 0.5, color = "black", size = 1) +
  facet_wrap(~RCTD_subtype2, scales = "free_y", nrow = 2) +
  theme_bw() +
  RotatedAxis() +
  theme(
    strip.text = element_blank(),
    strip.background = element_rect(color = "white", fill = "white"),
    panel.grid = element_blank(),
    legend.position = "top"
  ) +
  scale_fill_manual(values = plan) +
  guides(fill = guide_legend(nrow = 2, byrow = T))
ggsave(p_density_celltype, filename = file.path(monocle_outdir, "celltype density plot split.pdf"), height = 5, width = 20)
ggsave(p_density_celltype, filename = file.path(monocle_outdir, "celltype density plot split.png"), height = 5, width = 20, dpi = 300)






# BEAM 寻找拟时相关基因 ----------
# 会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
# 这一步很慢
### branch_point = 1  ----------
BEAM_res <- BEAM(cds[express_genes, ], branch_point = 1, progenitor_method = "duplicate")
saveRDS(BEAM_res, "BEAM_res.rds")
BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
sum(BEAM_res$pval < 0.05)

top_genes <- BEAM_res %>%
  arrange(qval) %>%
  slice_head(n = 50) %>%
  pull(gene_short_name) %>%
  as.character()

# tmp1 <- plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, pval < 0.05)), ],
tmp1 <- plot_genes_branched_heatmap(
  cds[top_genes, ],
  branch_point        = 1,
  num_clusters        = 3, # 这些基因被分成几个group
  cores               = 1,
  # branch_labels     = c("Cell fate 1", "Cell fate 2"),
  # hmcols            = NULL, #默认值
  # hmcols            = colorRampPalette(c("navy","white","firebrick3"))(50),
  # hmcols            = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
  hmcols              = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(60),
  branch_colors       = c("#979797", "#F05662", "#7990C8"), # pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
  use_gene_short_name = T,
  show_rownames       = T,
  return_heatmap      = T
)
write.csv(BEAM_res, file.path(monocle_outdir, "BEAM_group_branch1.csv"))

pdf(file.path(monocle_outdir, "BEAM_group_branch1_heatmap.pdf"), width = 10, height = 10)
tmp1$ph_res
dev.off()




### branch_point = 2 ----------
BEAM_res <- BEAM(cds[express_genes, ], branch_point = 2, progenitor_method = "duplicate")
BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
sum(BEAM_res$pval < 0.05)


top_genes <- BEAM_res %>%
  arrange(qval) %>%
  slice_head(n = 50) %>%
  pull(gene_short_name) %>%
  as.character()


# tmp1 <- plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, pval < 0.05)), ],
tmp1 <- plot_genes_branched_heatmap(
  cds[top_genes, ],
  branch_point        = 2,
  num_clusters        = 3, # 这些基因被分成几个group
  cores               = 1,
  # branch_labels     = c("Cell fate 1", "Cell fate 2"),
  # hmcols            = NULL, #默认值
  # hmcols            = colorRampPalette(c("navy","white","firebrick3"))(50),
  # hmcols            = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
  hmcols              = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(60),
  branch_colors       = c("#979797", "#F05662", "#7990C8"), # pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
  use_gene_short_name = T,
  show_rownames       = T,
  return_heatmap      = T
)
write.csv(BEAM_res, file.path(monocle_outdir, "BEAM_group_branch2.csv"))

pdf(file.path(monocle_outdir, "BEAM_group_branch2_heatmap.pdf"), width = 10, height = 10)
tmp1$ph_res
dev.off()






plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
plot_cell_trajectory(cds, color_by = "annotation_sub", size = 1, show_backbone = TRUE) + scale_color_manual(values = plan)
plot_cell_trajectory(cds, color_by = "annotation_sub") + facet_wrap("~annotation_sub", nrow = 1)
plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = TRUE) + scale_color_manual(values = plan)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_sub") + scale_color_manual(values = plan)

# p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
# p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
# p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
#
# p1=plot_cell_trajectory(cds, color_by = "CCL5")  + scale_color_gsea()




# de analysis --------
# 这里是把排序基因（express_genes）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
# 如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(
  cds[express_genes, ],
  cores = 1,
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)


# Time_diff <- Time_diff[, c(5, 2, 3, 4, 1, 6, 7)] # 把gene放前面，也可以不改
write.csv(Time_diff, file.path(monocle_outdir, "Time_diff_all.csv"), row.names = F)

# Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# p <- plot_pseudotime_heatmap(cds[Time_genes, ], num_clusters = 4, show_rownames = T, return_heatmap = T)
# ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

# Time_genes <- top_n(Time_diff, n = 50, desc(qval)) %>%
#   pull(gene_short_name) %>%
#   as.character()

Time_genes <- Time_diff %>% 
  slice_max(order_by = qval, n = 100, with_ties = FALSE) %>% 
  pull(gene_short_name) %>% 
  as.character()

p <- plot_pseudotime_heatmap(cds[Time_genes, ], num_clusters = 4, show_rownames = T, return_heatmap = T)
ggsave(file.path(monocle_outdir, "Time_heatmap_Top100.pdf"), p, width = 6, height = 10)

undebug(plot_pseudotime_heatmap)



# readme --------
readme <- "BEAM analysis: 
轨迹分析结果通常包含多个分支，主要是因为细胞的表达模式受各种因素影响而发生变化；本次拟时序分析有两个分支点，结果分别为BEAM_group_branch1_heatmap和BEAM_group_branch2_heatmap

BEAM_group_branch1_heatmap.pdf和Time_heatmap_Top100.pdf
这两个图的横轴是拟时间，从左至右时间依次增加；
纵轴是拟时间相关基因；
颜色越靠近蓝色标识基因表达量越低，越靠近红色标识基因表达量越高；
对于每个基因，从左至右的一条线（也就是在热图中的一行），代表了这个基因在拟时间过程中的表达量变化。每个格的颜色反应了该基因在该拟时间点的表达量。
"
writeLines(readme, file.path(monocle_outdir, "readme.txt"))







# some TF gene
gene_description <- read_tsv("data/gene_description_biomart.txt")
tf_gene <- gene_description$`Gene name`[str_detect(gene_description$`Gene description`, "transcription factor")] %>% .[!is.na(.)]

Time_genes %>% intersect(tf_gene)

p_mono_pseu_tf_gene <-
  plot_genes_in_pseudotime(
    cds[c("SOX6", "DACH2", "EBF1", "DACH1", "RUNX2"), ],
    markers = "State"
  ) +
  scale_color_manual(values = plan)
ggsave(file.path(monocle_outdir, "TF gene pseudotime.pdf"), p_mono_pseu_tf_gene, width = 5, height = 10)

# p_mono_traj_tf_gene <-
#   c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2") %>%
#   map(~{
#     plot_cell_trajectory(cds, color_by = .x)
#   }) %>%
#   patchwork::wrap_plots(nrow = 1)









# BEAM analysis
BEAM_res <- BEAM(cds[express_genes, ], branch_point = 1, progenitor_method = "duplicate")
# BEAM_res <- BEAM(cds[express_genes_monocle,], branch_point = 2, cores = 2)

BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]
head(BEAM_res)
#           gene_short_name         pval         qval
# CD79A               CD79A 2.757676e-73 7.782161e-70
# TCL1A               TCL1A 1.574889e-65 2.222168e-62
# IGLL5               IGLL5 2.356778e-64 2.216942e-61
# S100A9             S100A9 1.504319e-58 1.061297e-55
# S100A8             S100A8 6.028175e-57 3.402302e-54
# LINC00926       LINC00926 3.180527e-55 1.495908e-52
write.csv(BEAM_res, file.path(monocle_out_dir, "BEAM_res1.csv"), row.names = F)
plot_genes_branched_heatmap(
  cds[row.names(subset(BEAM_res, qval < 1e-20)), ],
  branch_point        = 1, # 绘制的是哪个分支
  num_clusters        = 3, # 分成几个cluster，根据需要调整
  cores               = 1,
  use_gene_short_name = T,
  show_rownames       = T
) # 有632个gene，太多了


# 选前100个基因可视化
BEAM_genes <- BEAM_res %>%
  arrange(qval) %>%
  slice_head(n = 100) %>%
  pull(gene_short_name) %>%
  as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes, ],
  branch_point   = 1,
  num_clusters   = 3, 
  show_rownames  = T, 
  return_heatmap = T
)
ggsave(file.path(monocle_out_dir, "BEAM_heatmap1.pdf"), p$ph_res, width = 10, height = 10)
