library(Seurat)
library(tidyverse)
library(monocle)
source("monocle/monocle fix bugs.R")
source("monocle/color_pallete.R")

# set default ------- 
monocle_out_dir <- "monocle_macrophage" 
dir.create(monocle_out_dir) 


sp_mac <- readRDS("sp_subcluster_LUSC/sp_mac.rds")


# monocle out -------
DimPlot(sp_mac, reduction = 'harmonyumap') # + scale_color_manual(values = col_plan_1)
DimPlot(sp_mac, reduction = 'harmonyumap', group.by = "group") # + scale_color_manual(values = col_plan_1)

VariableFeatures(sp_mac)

# subset
set.seed(1234)
meta_sp_mac_subset <-
  sp_mac[['annotation_mac']] %>%
  rownames_to_column("cell_id") %>%
  dplyr::filter(
    # !str_detect(annotation_mac, "like|low")
    str_detect(annotation_mac, "Macrophage")
  ) %>% 
  # group_by(annotation_mac) %>%
  slice_sample(n = 10000) %>%
  mutate(annotation_mac = fct_drop(annotation_mac))
# meta_sp_mac_subset <- sp_mac[['annotation_mac']] %>% rownames_to_column("cell_id")

# cds ---------
expr_dat_monocle <- as(LayerData(JoinLayers(sp_mac), "counts"), "sparseMatrix")
p_data_monocle <- sp_mac[[]]
f_data_monocle <- data.frame(
  gene_short_name = rownames(sp_mac),
  row.names       = rownames(sp_mac)
)

cds <- 
  newCellDataSet(
    expr_dat_monocle[, meta_sp_mac_subset$cell_id],
    phenoData   = new("AnnotatedDataFrame", p_data_monocle[meta_sp_mac_subset$cell_id, ]),
    featureData = new("AnnotatedDataFrame", f_data_monocle),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


# cds <- clusterCells(cds)


# #使用seurat选择的高变基因⚠️ ------
# express_genes_monocle <- VariableFeatures(sp_mac)
# cds <- setOrderingFilter(cds, express_genes_monocle)
# plot_ordering_genes(cds)

# ##使用clusters差异表达基因
# deg.cluster <- FindAllMarkers(JoinLayers(sp_mac))
# express_genes_monocle <- subset(deg.cluster, p_val_adj<0.05)$gene
# cds <- setOrderingFilter(cds, express_genes_monocle)
# plot_ordering_genes(cds)
# 
#使用monocle选择的高变基因⚠️
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 0.8 * dispersion_fit)$gene_id
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
plot_cell_trajectory(cds,color_by="annotation_mac", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds,color_by="group", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)

cds <- orderCells(cds, root_state = 5)

saveRDS(cds, file = file.path(monocle_out_dir, "cds.rds")) 

# some ploting
# cds <- orderCells(cds, root_state = 4)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1,show_backbone=TRUE) 
plot_cell_trajectory(cds, color_by="annotation_mac", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by="group", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_1)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_mac")
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "group", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)



# p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
# p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
# p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
# 
# p1=plot_cell_trajectory(cds, color_by = "CCL5")  + scale_color_gsea()

pdf(file.path(monocle_out_dir, "cell trajectory.pdf"), width = 8, height = 4)

plot_cell_trajectory(cds,color_by="Pseudotime", cell_size =0.1, show_backbone=TRUE) + theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="group", cell_size =0.1,show_backbone=TRUE) +
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

plot_cell_trajectory(cds,color_by="annotation_mac", cell_size =0.1, show_backbone=TRUE) + 
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right") + 
  guides(color = guide_legend(ncol = 2,override.aes = list(size = 3)))

plot_cell_trajectory(cds, color_by = "State",cell_size =1,show_backbone=TRUE) +
  scale_color_manual(values = col_plan_2) +
  theme(legend.position = "right")

dev.off()


pdf(file.path(monocle_out_dir, "cell trajectory - annotation - split by group.pdf"), width = 12, height = 4)
local({
  Biobase::pData(cds)$group <- all_data[['group']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$group)]
  
  plot_cell_trajectory(cds, color_by = "annotation_mac",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_1) +
    facet_wrap(~ group) +
    theme(legend.position = "right")+ 
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()



pdf(file.path(monocle_out_dir, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_mac") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation_mac)
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
  group_by(State, annotation_mac ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>% 
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_mac ))
dev.off()


pdf(file.path(monocle_out_dir, "Dimplot of sc_epi - group and cell annotation.pdf"), width = 8, height = 5)
DimPlot(
  sp_mac,
  group.by = "annotation_mac"
)

DimPlot(
  sp_mac,
  group.by = "group"
)
dev.off()


# cell barplot 
p_col_cell_state <-
  cds %>% 
  Biobase::pData() %>%
  dplyr::select(annotation_mac, State) %>%
  group_by(annotation_mac, State) %>%
  mutate(annotation_mac = as.character(annotation_mac)) %>%
  # mutate(annotation_mac = factor(annotation_mac, levels = c(
  #   'S1 fibroblasts', 
  #   'Peri-cryptal S2 fibroblasts',
  #   'Immediate early-response (IER) fibroblasts',
  #   'S3 fibroblasts',
  #   'Myofibroblasts'
  # ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation_mac, y = freq)) +
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
# gene_for_beam <- VariableFeatures(sp_mac)
gene_for_beam <- rownames(cds)
BEAM_res <- BEAM(cds[gene_for_beam,], branch_point = 3, cores = 10) 

saveRDS(BEAM_res, file.path(monocle_out_dir, "BEAM_res.rds"))

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
#           gene_short_name         pval         qval
# CD79A               CD79A 2.757676e-73 7.782161e-70
# TCL1A               TCL1A 1.574889e-65 2.222168e-62
# IGLL5               IGLL5 2.356778e-64 2.216942e-61
# S100A9             S100A9 1.504319e-58 1.061297e-55
# S100A8             S100A8 6.028175e-57 3.402302e-54
# LINC00926       LINC00926 3.180527e-55 1.495908e-52
write.csv(BEAM_res, file.path(monocle_out_dir, "BEAM_res.csv"), row.names = F)

# plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 1e-4)),],
#                             branch_point = 1, #绘制的是哪个分支
#                             num_clusters = 4, #分成几个cluster，根据需要调整
#                             cores = 5,
#                             use_gene_short_name = T,
#                             show_rownames = T) #有632个gene，太多了

#选前100个基因可视化
BEAM_genes <- BEAM_res %>%
  dplyr::filter(qval < 0.05) %>%
  dplyr::filter(!str_detect(gene_short_name , "^MT")) %>% 
  top_n(n = 100, dplyr::desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()

undebug(plot_genes_branched_heatmap)
p <- 
  cds[BEAM_genes,] %>%
  # cds[c('EVX2','CLEC2A','CNGA1', 'TCP10L2','RFX6'),] %>%
  plot_genes_branched_heatmap(
    branch_point = 1, 
    num_clusters = 4, # 3, 
    show_rownames = T, 
    return_heatmap = T
  )
ggsave(file.path(monocle_out_dir, "BEAM_heatmap_full.pdf"), p$ph_res, width = 6.5, height = 10)


plot_genes_branched_pseudotime(cds[c("WNT5A", "ITGB5"), ], color_by = "annotation_mac")
plot_genes_branched_pseudotime(cds[c("AREG", "EREG")], color_by = "State") + scale_y_log10()
plot_genes_branched_pseudotime(cds[c("MLPH", "PLA2G2A", "TNFRSF11A", "MUC4", "SLC28A2")], color_by = "State")

plot_genes_jitter(cds[c("WNT5A", "PLAUR", "SDC4")], color_by = "State")

cds2 <- cds
gene_to_plot <- exprs(cds)[c("WNT5A", "PLAUR", "SDC4"), ] %>% log1p %>% t %>% as.matrix()
pData(cds2) <- cbind(pData(cds2), gene_to_plot)

plot_cell_trajectory(cds2, color_by = "WNT5A") + ggsci::scale_color_gsea()


FeaturePlot(
  sp_mac,
  features = "WNT5A"
) 
DimPlot(
  sp_mac,
  group.by = "annotation_mac"
)
DimPlot(
  sp_mac,
  split.by = "group",
  group.by = "annotation_mac"
)

DimPlot(
  sp_mac,
  split.by = "project",
  group.by = "annotation_mac"
)


VlnPlot(
  sp_mac,
  features = "ABCC1",
  split.by = "group",
  group.by = "annotation_mac",
  pt.size = 0
)

DotPlot(
  sp_mac[, sp_mac$annotation_mac == "ALG11 Epithelial"],
  group.by = "group",
  features = str_split_1("ABCC1,TAP1,ABCA4,ABCC3,ABCB1,CFTR,TAP2,ABCG2,ABCB4,ABCC6,ABCC8", pattern = ",")
)
DotPlot(
  sp_mac[, sp_mac$annotation_mac == "BIRC5 Epithelial"],
  group.by = "group",
  features = str_split_1("ABCC1,TAP1,ABCA4,ABCC3,ABCB1,CFTR,TAP2,ABCG2,ABCB4,ABCC6,ABCC8", pattern = ",")
)
DotPlot(
  sp_mac[, sp_mac$annotation_mac == "CDKN1A Epithelial"],
  group.by = "group",
  features = str_split_1("ABCC1,TAP1,ABCA4,ABCC3,ABCB1,CFTR,TAP2,ABCG2,ABCB4,ABCC6,ABCC8", pattern = ",")
)

DotPlot(
  sp_mac,
  group.by = "annotation_mac",
  features = str_split_1("ABCC1,TAP1,ABCA4,ABCC3,ABCB1,CFTR,TAP2,ABCG2,ABCB4,ABCC6,ABCC8", pattern = ",")
)


DotPlot(
  sp_mac[, sp_mac$annotation_mac == "ALG11 Epithelial"],
  group.by = "group",
  features = str_split_1("CCND2,TP53,CDK4,MDM2,BCL2L1,CDKN1A,CD82,RRM2,PMAIP1,CCNB1,CCNB2,CCNE1", pattern = ",")
)
DotPlot(
  sp_mac[, sp_mac$annotation_mac == "BIRC5 Epithelial"],
  group.by = "group",
  features = str_split_1("CCND2,TP53,CDK4,MDM2,BCL2L1,CDKN1A,CD82,RRM2,PMAIP1,CCNB1,CCNB2,CCNE1", pattern = ",")
)
DotPlot(
  sp_mac[, sp_mac$annotation_mac == "CDKN1A Epithelial"],
  features = str_split_1("CCND2,TP53,CDK4,MDM2,BCL2L1,CDKN1A,CD82,RRM2,PMAIP1,CCNB1,CCNB2,CCNE1", pattern = ",")
)
DotPlot(
  sp_mac,
  group.by = "annotation_mac",
  features = str_split_1("CCND2,TP53,CDK4,MDM2,BCL2L1,CDKN1A,CD82,RRM2,PMAIP1,CCNB1,CCNB2,CCNE1", pattern = ",")
)


FeaturePlot(
  sp_mac,
  features = c('CCND2','TP53','CDK4','CDKN1A','CD82')
)

# branch gene enrichment 
beam_genes_ls <- 
  p$ph_res$tree_row %>% cutree(k = 4) %>%
  split(names(.), .) 

all_gene_id <- 
  rownames(sp_mac2) %>%
  clusterProfiler::bitr(
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = "org.Hs.eg.db"
  ) %>%
  .$ENTREZID

beam_gene_id_ls <-
  beam_genes_ls %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = "org.Hs.eg.db"
    ) %>%
      .$ENTREZID
  })


all_enrich_beam_ls <- 
  beam_gene_id_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      universe = all_gene_id,
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
    xx2 <-
      xx %>%
      map_chr(function(g){
        clusterProfiler::bitr(
          g,
          fromType = "ENTREZID",
          toType   = "SYMBOL",
          OrgDb    = "org.Hs.eg.db"
        ) %>%
          .$SYMBOL %>%
          str_c(collapse = "/")
      })
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


