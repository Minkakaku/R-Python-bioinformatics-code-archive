# ==================== Script 1: Monocle dispersion-based ordering genes ====================
# 拟时序分析 - 使用Monocle内置离散度筛选高变基因

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
library(readr)
library(ggpubr)

# 设置输出目录
monocle_outdir <- "Monocle_Tcell_method1"
dir.create(monocle_outdir, showWarnings = FALSE)

# 读取Seurat对象
scTcell <- readRDS("Epi_anno.rds")
#scTcell$subtype <- gsub("CD4_Tem", "CD8_Tem", scTcell$subtype)
Idents(scTcell) <- "annotation"

# 提取原始计数矩阵
count <- scTcell[["RNA"]]@counts
fdata <- data.frame(gene_short_name = row.names(count), row.names = row.names(count))
pdata <- scTcell@meta.data

# 构建monocle对象
cds <- newCellDataSet(count,
                      featureData = new("AnnotatedDataFrame", data = fdata),
                      phenoData = new("AnnotatedDataFrame", data = pdata),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

# 标准化和估计离散度
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
pData(cds)$UMI <- Matrix::colSums(exprs(cds))

# 颜色方案
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

### 方法1：基于离散度筛选高变基因
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, ordering_genes)

# 绘制排序基因的离散度图
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes.pdf"), height = 10, width = 10)
print(p)
dev.off()

# 降维和排序
cds <- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)

# 保存轨迹对象
saveRDS(cds, file.path(monocle_outdir, "cds_dispersion.rds"))

# 绘制基本轨迹图
p1 <- plot_cell_trajectory(cds, color_by = "subtype")
p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime")
p3 <- plot_cell_trajectory(cds, color_by = "State")

pdf(file.path(monocle_outdir, "trajectory_basic.pdf"), width = 8, height = 6)
print(p1 + scale_color_manual(values = plan) + theme(legend.position = "right"))
print(p2 + theme(legend.position = "right"))
print(p3 + scale_color_manual(values = plan) + theme(legend.position = "right"))
dev.off()

# 分面轨迹
p_facet_subtype <- plot_cell_trajectory(cds, color_by = "subtype") + 
  facet_wrap(~subtype, nrow = 2) + scale_color_manual(values = plan)
p_facet_state <- plot_cell_trajectory(cds, color_by = "State") + 
  facet_wrap(~State, nrow = 2) + scale_color_manual(values = plan)

pdf(file.path(monocle_outdir, "trajectory_facet.pdf"), width = 18, height = 9)
print(p_facet_subtype)
print(p_facet_state)
dev.off()

# 复杂轨迹图
pdf(file.path(monocle_outdir, "trajectory_tree.pdf"), width = 15, height = 12)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "subtype") +
  scale_color_manual(values = plan) + facet_wrap(~ subtype)
dev.off()

# 细胞状态分布柱状图
p_bar <- pData(cds) %>%
  dplyr::select(subtype, State) %>%
  group_by(subtype, State) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(subtype) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x = subtype, y = freq, fill = State)) +
  geom_col() +
  scale_fill_manual(values = plan) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme_pubr() +
  theme(legend.position = "right")
ggsave(p_bar, filename = file.path(monocle_outdir, "cell_state_barplot.pdf"), height = 6, width = 8)

### BEAM分析（分支依赖基因）
BEAM_res <- BEAM(cds[ordering_genes,], branch_point = 1, progenitor_method = 'duplicate')
BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
write.csv(BEAM_res, file.path(monocle_outdir, "BEAM_results.csv"), row.names = FALSE)

# 绘制分支热图（显著基因）
sig_genes <- row.names(subset(BEAM_res, pval < 0.05))
if (length(sig_genes) > 0) {
  hm <- plot_genes_branched_heatmap(cds[sig_genes,],
                                    branch_point = 1,
                                    num_clusters = 3,
                                    cores = 1,
                                    hmcols = colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(60),
                                    branch_colors = c("#979797", "#F05662", "#7990C8"),
                                    use_gene_short_name = TRUE,
                                    show_rownames = TRUE,
                                    return_heatmap = TRUE)
  pdf(file.path(monocle_outdir, "BEAM_heatmap.pdf"), width = 10, height = 10)
  print(hm$ph_res)
  dev.off()
}

# 前100个基因的热图
top100_genes <- BEAM_res %>% arrange(qval) %>% slice_head(n = 100) %>% pull(gene_short_name)
hm_top <- plot_genes_branched_heatmap(cds[top100_genes,],
                                      branch_point = 1,
                                      num_clusters = 3,
                                      cores = 1,
                                      hmcols = colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(60),
                                      branch_colors = c("#979797", "#F05662", "#7990C8"),
                                      use_gene_short_name = TRUE,
                                      show_rownames = TRUE,
                                      return_heatmap = TRUE)
pdf(file.path(monocle_outdir, "BEAM_heatmap_top100.pdf"), width = 10, height = 10)
print(hm_top$ph_res)
dev.off()

### 差异基因与拟时间关系分析
Time_diff <- differentialGeneTest(cds[ordering_genes,],
                                  cores = 1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)]  # 调整列顺序
write.csv(Time_diff, file.path(monocle_outdir, "Time_diff_all.csv"), row.names = FALSE)

# 绘制前100个拟时间相关基因的热图
Time_genes <- Time_diff %>% top_n(100, desc(qval)) %>% pull(gene_short_name)
p_heat <- plot_pseudotime_heatmap(cds[Time_genes,],
                                  num_clusters = 4,
                                  show_rownames = TRUE,
                                  return_heatmap = TRUE)
ggsave(file.path(monocle_outdir, "Time_heatmap_top100.pdf"), p_heat$ph_res, width = 5, height = 10)

### 转录因子基因分析（示例）
if (file.exists("data/gene_description_biomart.txt")) {
  gene_desc <- read_tsv("data/gene_description_biomart.txt")
  tf_genes <- gene_desc$`Gene name`[grepl("transcription factor", gene_desc$`Gene description`, ignore.case = TRUE)]
  tf_genes <- tf_genes[!is.na(tf_genes)]
  
  # 检查Time_genes中的TF
  intersect(Time_genes, tf_genes)
  
  # 绘制选定TF的拟时间表达
  tf_to_plot <- intersect(c("SOX6", "DACH2", "EBF1", "DACH1", "RUNX2"), rownames(cds))
  if (length(tf_to_plot) > 0) {
    p_tf <- plot_genes_in_pseudotime(cds[tf_to_plot,], color_by = "State") +
      scale_color_manual(values = plan)
    ggsave(file.path(monocle_outdir, "TF_pseudotime.pdf"), p_tf, width = 5, height = 10)
  }
}

cat("Script 1 completed. Results saved in", monocle_outdir, "\n")
