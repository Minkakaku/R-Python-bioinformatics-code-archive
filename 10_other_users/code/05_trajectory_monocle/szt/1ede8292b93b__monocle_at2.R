library(Seurat)
library(tidyverse)
library(monocle)
source("monocle fix bugs.R")

# 2025.2.10
sc_monolce <- readRDS("subcluster/sc_at2_2.rds") 
sc_monolce <- sc_monolce[, sc_monolce$annotation != "doublet"]


DefaultAssay(sc_monolce) <- "RNA" 
# sc_monolce[['sketch']]   <- NULL 

# set default ------- 
monocle_out_dir <- "monocle_project_at2" 
dir.create(monocle_out_dir) 

col_plan_1 <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
                "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", "#de4307", 
                "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", "#fa6e57", "#f69e53", 
                "#248888", "#999999", "#E7475E", "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", 
                "#FF6D31", "#f1ac9d", "#f06966", "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", 
                "#FF7A5A", "#8ED2C9", "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", 
                "#ffc300", "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
                "#f17d80", "#737495")

col_plan_2 <- c("#1F77B4FF","#FF7F0EFF","#2CA02CFF","#D62728FF","#9467BDFF",
                "#8C564BFF","#E377C2FF","#7F7F7FFF","#BCBD22FF","#17BECFFF",
                "#AEC7E8FF","#FFBB78FF","#98DF8AFF","#FF9896FF","#C5B0D5FF",
                "#C49C94FF","#F7B6D2FF","#C7C7C7FF","#DBDB8DFF","#9EDAE5FF")

# monocle out -------
DimPlot(sc_monolce, reduction = 'umap') + scale_color_manual(values = col_plan_1)

VariableFeatures(sc_monolce)

# subset
set.seed(1234)
meta_sc_monolce_subset <-
  sc_monolce[['annotation']] %>%
  rownames_to_column("cell_id") %>%
  group_by(annotation) %>%
  # slice_sample(n = 5000)
  slice_sample(prop = 1)

# cds ---------
expr_dat_monocle <- as(LayerData(JoinLayers(sc_monolce), "counts"), "sparseMatrix")
p_data_monocle <- sc_monolce[[]]
f_data_monocle <- data.frame(
  gene_short_name = rownames(sc_monolce),
  row.names       = rownames(sc_monolce)
)

cds <- 
  newCellDataSet(
    expr_dat_monocle[, meta_sc_monolce_subset$cell_id],
    phenoData   = new("AnnotatedDataFrame", p_data_monocle[meta_sc_monolce_subset$cell_id, ]),
    featureData = new("AnnotatedDataFrame", f_data_monocle),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


# cds <- clusterCells(cds)


##使用seurat选择的高变基因⚠️ ------
# express_genes_monocle <- VariableFeatures(sc_monolce)
# cds <- setOrderingFilter(cds, express_genes_monocle)
# plot_ordering_genes(cds)

# ##使用clusters差异表达基因
# deg.cluster <- FindAllMarkers(JoinLayers(sc_monolce))
# express_genes_monocle <- subset(deg.cluster, p_val_adj<0.05)$gene
# cds <- setOrderingFilter(cds, express_genes_monocle)
# plot_ordering_genes(cds)
# 
#使用monocle选择的高变基因⚠️
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
express_genes_monocle <- disp.genes
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

# #这一步输入的expressed_genes来自于步骤4。
# fData(cds)$num_cells_expressed <- rowSums(exprs(cds) > 0)
# expressed_genes2 <- row.names(subset(fData(cds), num_cells_expressed >= 10))
# diff <- differentialGeneTest(cds[expressed_genes2, ],fullModelFormulaStr="~annotation",cores=1) 
# 
# head(diff)
# 
# ##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
# deg <- subset(diff, qval < 0.01) #选出2829个基因
# deg <- deg[order(deg$qval,decreasing=F),]
# head(deg)
# 
# ##差异基因的结果文件保存
# write.table(deg,file="train.monocle.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
# 
# ## 轨迹构建基因可视化
# ordergene <- rownames(deg) 
# cds <- setOrderingFilter(cds, ordergene)  
# plot_ordering_genes(cds)


cds <- reduceDimension(cds)
undebug(orderCells)
cds <- orderCells(cds)

# Error:
#   ! `nei()` was deprecated in igraph 2.1.0 and is now defunct.
# ℹ Please use `.nei()` instead.
# install.packages('https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz', type = "source", repos = NULL)
# 

# plot ----------
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds,color_by="annotation", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds,color_by="project", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)

cds <- orderCells(cds, root_state = 3)

saveRDS(cds, file = file.path(monocle_out_dir, "cds.rds"))

# some ploting
# cds <- orderCells(cds, root_state = 4)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1,show_backbone=TRUE) 
plot_cell_trajectory(cds, color_by="annotation", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by="group", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation")
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "group", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)



# p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
# p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
# p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
# 
# p1=plot_cell_trajectory(cds, color_by = "CCL5")  + scale_color_gsea()

pdf(file.path(monocle_out_dir, "cell trajectory.pdf"), width = 8, height = 4)

plot_cell_trajectory(cds,color_by="Pseudotime", cell_size =1, show_backbone=TRUE) + theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="project", cell_size =1,show_backbone=TRUE) +
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right") + 
  guides(color = guide_legend(ncol = 2,override.aes = list(size = 3)))

plot_cell_trajectory(cds,color_by="annotation", cell_size =1, show_backbone=TRUE) + 
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right") + 
  guides(color = guide_legend(ncol = 2,override.aes = list(size = 3)))

plot_cell_trajectory(cds, color_by = "State",cell_size =1,show_backbone=TRUE) +
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

dev.off()


pdf(file.path(monocle_out_dir, "cell trajectory - annotation - split by group.pdf"), width = 12, height = 4)
local({
  Biobase::pData(cds)$project <- all_data[['project']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$group)]
  
  plot_cell_trajectory(cds, color_by = "annotation",cell_size = 1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_1) +
    facet_wrap(~ project) +
    theme(legend.position = "right")+ 
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()



pdf(file.path(monocle_out_dir, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 1, y = 2, color_by = "annotation") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation)
dev.off()


pdf(file.path(monocle_out_dir, "barplot tof cell props group by State - group.pdf"), width = 8, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, project) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = project))
dev.off()

pdf(file.path(monocle_out_dir, "barplot of cell props group by State - annotation.pdf"), width = 8, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, annotation ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation ))
dev.off()


pdf(file.path(monocle_out_dir, "Dimplot - group and cell annotation.pdf"), width = 8, height = 5)
DimPlot(
  sc_monolce,
  group.by = "annotation"
)

DimPlot(
  sc_monolce,
  group.by = "project"
)
dev.off()


# cell barplot 
p_col_cell_state <-
  cds %>% 
  Biobase::pData() %>%
  dplyr::select(annotation, State) %>%
  group_by(annotation, State) %>%
  mutate(annotation = as.character(annotation)) %>%
  # mutate(annotation = factor(annotation, levels = c(
  #   'S1 fibroblasts', 
  #   'Peri-cryptal S2 fibroblasts',
  #   'Immediate early-response (IER) fibroblasts',
  #   'S3 fibroblasts',
  #   'Myofibroblasts'
  # ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation, y = freq)) +
  geom_col(aes(fill = State)) +
  scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
ggsave(p_col_cell_state, filename = file.path(monocle_out_dir, "cell stat barplot.pdf"), height = 4, width = 8)


# de analysis -----------
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(
  cds[express_genes_monocle,], 
  cores = 10, 
  fullModelFormulaStr = "~sm.ns(Pseudotime)"
)

# Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
write.csv(Time_diff, file.path(monocle_out_dir, "Time_diff_all.csv"), row.names = F)

# Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
# p=plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
# ggsave("Time_heatmapAll.pdf", p, width = 5, height = 10)

Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave(file.path(monocle_out_dir, "Time_heatmapTop100.pdf"), p, width = 5, height = 10)

undebug(plot_pseudotime_heatmap)


# some TF gene
gene_description <- read_tsv("data/gene_description_biomart.txt")
tf_gene <- gene_description$`Gene name`[str_detect(gene_description$`Gene description`, 'transcription factor')] %>% .[!is.na(.)]

Time_genes %>% intersect(tf_gene)

p_mono_pseu_tf_gene <- 
  plot_genes_in_pseudotime(
    cds[c("SOX6",  "DACH2", "EBF1",  "DACH1", "RUNX2"), ], 
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


# BEAM analysis ------------------
undebug(BEAM)
# gene_for_beam <- VariableFeatures(sc_epi2)
gene_for_beam <- cds %>% fData() %>% {rownames(.)[.$use_for_ordering]}
# BEAM_res_b1 <- BEAM(cds[gene_for_beam,], branch_point = 1, cores = 10) 
BEAM_res_b2 <- BEAM(cds[gene_for_beam,], branch_point = 2, cores = 10) 
# BEAM_res_b3 <- BEAM(cds[gene_for_beam,], branch_point = 3, cores = 10) 
# BEAM_res_b4 <- BEAM(cds[gene_for_beam,], branch_point = 4, cores = 10) 

# saveRDS(BEAM_res_b1, file.path(monocle_out_dir, "BEAM_res_b1.rds"))
saveRDS(BEAM_res_b2, file.path(monocle_out_dir, "BEAM_res_b2.rds"))
# saveRDS(BEAM_res_b3, file.path(monocle_out_dir, "BEAM_res_b3.rds"))
# saveRDS(BEAM_res_b4, file.path(monocle_out_dir, "BEAM_res_b4.rds"))

# BEAM_res_b1 <- BEAM_res_b1[order(BEAM_res_b1$qval),]
# BEAM_res_b1 <- BEAM_res_b1[,c("gene_short_name", "pval", "qval")]

BEAM_res_b2 <- BEAM_res_b2[order(BEAM_res_b2$qval),]
BEAM_res_b2 <- BEAM_res_b2[,c("gene_short_name", "pval", "qval")]

# BEAM_res_b3 <- BEAM_res_b3[order(BEAM_res_b3$qval),]
# BEAM_res_b3 <- BEAM_res_b3[,c("gene_short_name", "pval", "qval")]
# 
# BEAM_res_b4 <- BEAM_res_b4[order(BEAM_res_b4$qval),]
# BEAM_res_b4 <- BEAM_res_b4[,c("gene_short_name", "pval", "qval")]

# write.csv(BEAM_res_b1, file.path(monocle_out_dir, "BEAM_res_b1.csv"), row.names = F)
write.csv(BEAM_res_b2, file.path(monocle_out_dir, "BEAM_res_b2.csv"), row.names = F)
# write.csv(BEAM_res_b3, file.path(monocle_out_dir, "BEAM_res_b3.csv"), row.names = F)
# write.csv(BEAM_res_b4, file.path(monocle_out_dir, "BEAM_res_b4.csv"), row.names = F)

# plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 1e-4)),],
#                             branch_point = 1, #绘制的是哪个分支
#                             num_clusters = 4, #分成几个cluster，根据需要调整
#                             cores = 5,
#                             use_gene_short_name = T,
#                             show_rownames = T) #有632个gene，太多了

#选前100个基因可视化
BEAM_genes_b1 <- BEAM_res_b1 %>%
  dplyr::filter(qval < 0.05) %>%
  dplyr::filter(!str_detect(gene_short_name , "^MT")) %>% 
  # top_n(n = 100, dplyr::desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()

p_b1 <- 
  cds[BEAM_genes_b1,] %>%
  # cds[c('EVX2','CLEC2A','CNGA1', 'TCP10L2','RFX6'),] %>%
  plot_genes_branched_heatmap(
    branch_point = 1, 
    num_clusters = 5, # 3, 
    show_rownames = F, 
    return_heatmap = T
  )

BEAM_genes_b2 <- BEAM_res_b2 %>%
  dplyr::filter(qval < 0.05) %>%
  dplyr::filter(!str_detect(gene_short_name , "^MT")) %>% 
  # top_n(n = 100, dplyr::desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()

p_b2 <- 
  cds[BEAM_genes_b2,] %>%
  # cds[c('EVX2','CLEC2A','CNGA1', 'TCP10L2','RFX6'),] %>%
  plot_genes_branched_heatmap(
    branch_point = 1, 
    num_clusters = 5, # 3, 
    show_rownames = F, 
    return_heatmap = T
  )

BEAM_genes_b3 <- BEAM_res_b3 %>%
  dplyr::filter(qval < 0.05) %>%
  dplyr::filter(!str_detect(gene_short_name , "^MT")) %>% 
  # top_n(n = 100, dplyr::desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()

p_b3 <- 
  cds[BEAM_genes_b3,] %>%
  # cds[c('EVX2','CLEC2A','CNGA1', 'TCP10L2','RFX6'),] %>%
  plot_genes_branched_heatmap(
    branch_point = 1, 
    num_clusters = 5, # 3, 
    show_rownames = F, 
    return_heatmap = T
  )

BEAM_genes_b4 <- BEAM_res_b4 %>%
  dplyr::filter(qval < 0.05) %>%
  dplyr::filter(!str_detect(gene_short_name , "^MT")) %>% 
  # top_n(n = 100, dplyr::desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()

p_b4 <- 
  cds[BEAM_genes_b4,] %>%
  # cds[c('EVX2','CLEC2A','CNGA1', 'TCP10L2','RFX6'),] %>%
  plot_genes_branched_heatmap(
    branch_point = 1, 
    num_clusters = 5, # 3, 
    show_rownames = F, 
    return_heatmap = T
  )

# ggsave(file.path(monocle_out_dir, "BEAM_heatmap_full_b1.pdf"), p_b1$ph_res, width = 6.5, height = 10)
ggsave(file.path(monocle_out_dir, "BEAM_heatmap_full_b2.pdf"), p_b2$ph_res, width = 6.5, height = 10)
# ggsave(file.path(monocle_out_dir, "BEAM_heatmap_full_b3.pdf"), p_b3$ph_res, width = 6.5, height = 10)
# ggsave(file.path(monocle_out_dir, "BEAM_heatmap_full_b4.pdf"), p_b4$ph_res, width = 6.5, height = 10)


# plot_genes_branched_pseudotime(cds[c("CLCA1", "TFF3", "SLC26A3")], color_by = "State")
# plot_genes_branched_pseudotime(cds[c("AREG", "EREG")], color_by = "State") + scale_y_log10()
# plot_genes_branched_pseudotime(cds[c("MLPH", "PLA2G2A", "TNFRSF11A", "MUC4", "SLC28A2")], color_by = "State")
# 
# plot_genes_jitter(cds[c("AREG", "EREG")], color_by = "State")


# branch gene enrichment  
beam_genes_ls <- 
  p_b2$ph_res$tree_row %>% cutree(k = 4) %>%
  split(names(.), .) 


beam_gene_db <-
  beam_genes_ls %>%
  unlist(use.names = FALSE) %>% 
  unique() %>%
  clusterProfiler::bitr(
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = "org.Mm.eg.db"
  ) 

beam_gene_id_ls <-
  beam_genes_ls %>%
  map(function(x){
    xx <- beam_gene_db$ENTREZID[match(x, beam_gene_db$SYMBOL)]
    xx[!is.na(xx)]
  })


all_enrich_beam_ls <- 
  beam_gene_id_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Mm.eg.db",
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.2
    )
  })

# enrichment 
all_enrich_beam_ls <-
  all_enrich_beam_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    xx2 <- xx %>%
      map_chr(function(g){
        beam_gene_db$SYMBOL[match(g, beam_gene_db$ENTREZID)] %>% str_c(collapse = "/")
      })
    # print(xx2)
    # xx2 <-
    #   xx %>%
    #   map_chr(function(g){
    #     clusterProfiler::bitr(
    #       g,
    #       fromType = "ENTREZID",
    #       toType   = "SYMBOL",
    #       OrgDb    = "org.Mm.eg.db"
    #     ) %>%
    #       .$SYMBOL %>%
    #       str_c(collapse = "/")
    #   })
    x@result$geneSymbol <- xx2
    x
  })



# enrichment plot
p_all_enrich_beam_ls <-
  all_enrich_beam_ls %>%
  imap(function(enrich, name){
    if(nrow(as.data.frame(enrich)) < 1) return(NULL)
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(name)
  })

# export
all_enrich_beam_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file.path(monocle_out_dir, "monocle BEAM gene enrichment analysis table - full.xlsx"),
    asTable = TRUE
  )

pdf(
  file.path(monocle_out_dir, "monocle BEAM gene enrichment analysis dotplot - full.pdf"), 
  height = 5, 
  width = 8
)
p_all_enrich_beam_ls %>%  walk(print)
dev.off()