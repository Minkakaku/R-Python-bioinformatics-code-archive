# 拟时序 ------------
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


all_data <- readRDS("all_data_monocle_260113.rds")

cells_to_keep <- grepl("^(Tumor_|Normal_)", all_data@meta.data$annotation_sub)
all_data_sub <- all_data[, cells_to_keep]
Idents(all_data_sub) <- "annotation_sub"

table(all_data_sub@meta.data$annotation_sub)
allCells <- names(Idents(all_data_sub))
allType <- levels(Idents(all_data_sub))

choose_Cells <- unlist(lapply(allType, function(x) {
  cgCells <- allCells[Idents(all_data_sub) == x]
  num <- round(length(cgCells) * 0.2)
  cg <- sample(cgCells, num)
  cg
}))

all_data_sub <- subset(all_data_sub, cells = choose_Cells)
table(all_data_sub@meta.data$annotation_sub)


saveRDS(all_data_sub, "all_data_sub.rds")

monocle_outdir <- "Monocle_epithelial"
dir.create(monocle_outdir)
# expr_dat_monocle <- as(LayerData(JoinLayers(sc_fib), "counts"), "sparseMatrix")
# all_data_subll<-subset(all_data_sub,subtype == 'M2'|subtype == 'M1'|subtype == 'Cycling M2')

count <- GetAssayData(object = JoinLayers(all_data_sub)[["RNA"]], layer = "counts")
fdata <- data.frame(gene_short_name = row.names(count), row.names = row.names(count))
pdata <- all_data_sub@meta.data
# 构建monocle对象
cds <- newCellDataSet(count,
  featureData         = new("AnnotatedDataFrame", data = fdata),
  phenoData           = new("AnnotatedDataFrame", data = pdata),
  lowerDetectionLimit = 0.5,
  expressionFamily    = negbinomial.size()
)
# 标准化sacle及筛选，即使原先的Seurat对象已经归一化过，官方推荐在monocle中重新进行归一化处理
# 标准化细胞之间的mRNA的差异
cds <- estimateSizeFactors(cds)
# 离散度值可以帮助我们进行后续的差异分析
cds <- estimateDispersions(cds)

# 在fData(cds)中添加“num_cell_expressed”列，可用于过滤细胞
cds <- detectGenes(cds, min_expr = 0.1)

head(pData(cds))
head(fData(cds))

expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
length(expressed_genes)

# 向Epidata(cds)中添加“UMI”列
pdata$UMI <- Matrix::colSums(exprs(cds))

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

### 筛选高变基因(离散度高的gene) ----------
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-1.pdf"), height = 10, width = 10)
print(p)
dev.off()
cds <- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "annotation_sub")
pdf(file.path(monocle_outdir, "monocle_cluster-1.pdf"), height = 10, width = 10)
print(p)
dev.off()

# cds=orderCells(cds,root_state =3)
p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-1.pdf"))
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "annotation_sub") + facet_wrap(~annotation_sub, nrow = 4)
pdf(file.path(monocle_outdir, "monocle_facet-annotation_sub-1.pdf"), height = 25, width = 25)
print(p)
dev.off()

p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-1.pdf"), height = 10, width = 10)
print(p)
dev.off()

saveRDS(cds, file.path(monocle_outdir, "cds_disp-1.rds"))


### 高变基因 -----------
# group1_data<- group1_data[['RNA']] %>% as('Assay5')
# group1_data <- group1_data%>%NormalizeData(verbose=FALSE)%>%FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
# group1_data <- FindVariableFeatures(group1_data, selection.method = "vst", nfeatures = 2000)
ordergene <- VariableFeatures(all_data_sub)
cds <- setOrderingFilter(cds, ordergene)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-2.pdf"), height = 10, width = 10)
print(p)
dev.off()
cds <- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "annotation_sub")
pdf(file.path(monocle_outdir, "monocle_cluster-2.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-2.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "annotation_sub") + facet_wrap(~annotation_sub, nrow = 4)
pdf(file.path(monocle_outdir, "monocle_facet-annotation_sub-2.pdf"), height = 25, width = 20)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-2.pdf"), height = 10, width = 10)
print(p)
dev.off()

saveRDS(cds, file.path(monocle_outdir, "cds_VariableFeatures-2.rds"))


### 差异基因 ----------
deg.cluster <- FindAllMarkers(JoinLayers(all_data_sub))
express_genes <- subset(deg.cluster, p_val_adj < 0.05)$gene
express_genes <- unique(express_genes)
cds <- setOrderingFilter(cds, express_genes)
p <- plot_ordering_genes(cds)
pdf(file.path(monocle_outdir, "monocle_ordering_genes-3.pdf"), height = 10, width = 10)
print(p)
dev.off()
cds <- reduceDimension(cds, reduction_method = "DDRTree")
cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "annotation_sub")
pdf(file.path(monocle_outdir, "monocle_annotation_sub-3.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
pdf(file.path(monocle_outdir, "monocle_Pseudotime-3.pdf"), height = 10, width = 10)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "annotation_sub") + facet_wrap(~annotation_sub, nrow = 4)
pdf(file.path(monocle_outdir, "monocle_facet-annotation_sub-3.pdf"), height = 25, width = 20)
print(p)
dev.off()
p <- plot_cell_trajectory(cds, color_by = "State")
pdf(file.path(monocle_outdir, "monocle_state-3.pdf"), height = 10, width = 10)
print(p)
dev.off()

saveRDS(cds, file.path(monocle_outdir, "cds_deg-3.rds"))


## 寻找拟时相关基因
# 会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
# 这一步很慢
BEAM_res <- BEAM(cds[ordergene, ], branch_point = 1, progenitor_method = "duplicate")
BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
sum(BEAM_res$pval < 0.05)
tmp1 <- plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, pval < 0.05)), ],
  branch_point = 2,
  num_clusters = 3, # 这些基因被分成几个group
  cores = 1,
  # branch_labels = c("Cell fate 1", "Cell fate 2"),
  # hmcols = NULL, #默认值
  # hmcols = colorRampPalette(c("navy","white","firebrick3"))(50),
  # hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
  hmcols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(60),
  branch_colors = c("#979797", "#F05662", "#7990C8"), # pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
  use_gene_short_name = T,
  show_rownames = T,
  return_heatmap = T
)
write.csv(BEAM_res, "./monocle/group2_branch2.csv")
pdf("./monocle/group2_branch2_heatmap.pdf", width = 10, height = 10)
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

cds <- readRDS(file.path(monocle_outdir, "cds_VariableFeatures-2.rds"))
pData(cds)$annotation_sub <- factor(pData(cds)$annotation_sub, levels = c(paste0("Normal_C", 0:13), paste0("Tumor_C", 0:9)))

pdf(file.path(monocle_outdir, "cell trajectory.pdf"), width = 8, height = 6)
plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE) + theme(legend.position = "right")
plot_cell_trajectory(cds, color_by = "annotation_sub", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = TRUE) +
  scale_color_manual(values = plan) +
  theme(legend.position = "right")
dev.off()


pdf(file.path(monocle_outdir, "cell trajectory - - split by annotation_sub.pdf"), width = 30, height = 25)
plot_cell_trajectory(cds, color_by = "annotation_sub") + facet_wrap("~annotation_sub", nrow = 4) + theme(legend.position = "right") + scale_color_manual(values = plan)
dev.off()

pdf(file.path(monocle_outdir, "cell trajectory - - split by state.pdf"), width = 15, height = 6)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1) + scale_color_manual(values = plan)
dev.off()

pdf(file.path(monocle_outdir, "cell trajectory - annotation - split by group.pdf"), width = 12, height = 5)
local({
  pData(cds)$group <- pData(cds)$group %>% factor(levels = c("PA", "PC"))
  plot_cell_trajectory(cds, color_by = "annotation_sub", size = 1, show_backbone = TRUE) +
    scale_color_manual(values = plan) +
    facet_wrap(~group) +
    theme(legend.position = "right")
})
dev.off()



pdf(file.path(monocle_outdir, "cell trajectory - tree.pdf"), width = 30, height = 20)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_sub") + theme(legend.position = "right") +
  scale_color_manual(values = plan) +
  facet_wrap(~annotation_sub)
dev.off()

# cell barplot
p_col_cell_state <-
  cds %>%
  pData() %>%
  dplyr::select(annotation_sub, State) %>%
  group_by(annotation_sub, State) %>%
  mutate(annotation_sub = as.character(annotation_sub)) %>%
  mutate(annotation_sub = factor(annotation_sub, levels = c(
    paste0("Normal_C", 0:13), paste0("Tumor_C", 0:9)
  ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) %>%
  ggplot(aes(x = annotation_sub, y = freq)) +
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

# de analysis
# 这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
# 如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(
  cds[express_genes_monocle, ],
  cores = 1,
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)

Time_diff <- Time_diff[, c(5, 2, 3, 4, 1, 6, 7)] # 把gene放前面，也可以不改
write.csv(Time_diff, file.path(monocle_out_dir, "Time_diff_all.csv"), row.names = F)

# Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
# ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()
p <- plot_pseudotime_heatmap(cds[Time_genes, ], num_clusters = 4, show_rownames = T, return_heatmap = T)
ggsave(file.path(BEAM_res, "Time_heatmapTop100.pdf"), p, width = 5, height = 10)

undebug(plot_pseudotime_heatmap)


# some TF gene
gene_description <- read_tsv("data/gene_description_biomart.txt")
tf_gene <- gene_description$`Gene name`[str_detect(gene_description$`Gene description`, "transcription factor")] %>% .[!is.na(.)]

Time_genes %>% intersect(tf_gene)

p_mono_pseu_tf_gene <-
  plot_genes_in_pseudotime(
    cds[c("SOX6", "DACH2", "EBF1", "DACH1", "RUNX2"), ],
    markers = "State"
  ) +
  scale_color_manual(values = col_plan_2)

ggsave(file.path(monocle_out_dir, "TF gene pseudotime.pdf"), p_mono_pseu_tf_gene, width = 5, height = 10)

# p_mono_traj_tf_gene <-
#   c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2") %>%
#   map(~{
#     plot_cell_trajectory(cds, color_by = .x)
#   }) %>%
#   patchwork::wrap_plots(nrow = 1)


# BEAM analysis
BEAM_res <- BEAM(cds[ordergene, ], branch_point = 1, progenitor_method = "duplicate")
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
  cds[row.names(subset(
    BEAM_res,
    qval < 1e-20
  )), ],
  branch_point = 1, # 绘制的是哪个分支
  num_clusters = 3, # 分成几个cluster，根据需要调整
  cores = 1,
  use_gene_short_name = T,
  show_rownames = T
) # 有632个gene，太多了
# 选前100个基因可视化
BEAM_genes <- BEAM_res %>%
  arrange(qval) %>%
  slice_head(n = 100) %>%
  pull(gene_short_name) %>%
  as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes, ],
  branch_point = 1,
  num_clusters = 3, show_rownames = T, return_heatmap = T
)
ggsave(file.path(monocle_out_dir, "BEAM_heatmap1.pdf"), p$ph_res, width = 10, height = 10)
