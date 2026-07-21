library(Seurat)
library(tidyverse)
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(DOSE)
source("monocle fix bugs.R")
source('color_pallete.R')

# set default ------- 
monocle_out_dir <- "monocle_result" 
# dir.create(monocle_out_dir) 


# Tumor ------------------------------------------------------------------
monocle_out_dir <- "monocle_result/Tumor_5_3" 
# dir.create(monocle_out_dir)
Luminal4 <- readRDS('result_8um/Luminal_reannotation3/Luminal4_final.rds')
Luminal4$clusters %>% table()

Luminal4$annotation_monocle <- 
  Luminal4$clusters %>% 
  str_replace_all(c('0' = 'LE_0', 
                    '1' = 'LE_1', 
                    '2' = 'LE_2', 
                    '3' = 'LE_3', 
                    '4' = 'LE_4', 
                    '5' = 'LE_5'))

Luminal4$annotation_monocle %>% table()


NE <- readRDS('result_8um/NE_subtype/filtered/NE_filterd.rds')
NE$Spatial.008um_snn_res.0.25 %>% table()
NE$clusters <- as.character(NE$Spatial.008um_snn_res.0.25)
NE$clusters <- factor(NE$clusters, levels = c(0:4))
NE$clusters %>% table()

# NE$annotation_monocle <- 
#   NE$clusters %>% 
#   str_replace_all(c('0' = 'NE_0', 
#                     '1' = 'NE_1', 
#                     '2' = 'NE_2', 
#                     '3' = 'NE_3', 
#                     '4' = 'NE_4'))

NE$annotation_monocle <- 
  NE$clusters %>% 
  str_replace_all(c('0' = 'NE', 
                    '1' = 'NE', 
                    '2' = 'NE', 
                    '3' = 'NE', 
                    '4' = 'NE'))
NE$annotation_monocle %>% table()


Tumor <- merge(Luminal4, NE)
Tumor$annotation_monocle %>% table()


#Tumor <- subset(Tumor, annotation_monocle %in% c('LE_0', 'LE_3', 'NE_0', 'NE_1', 'NE_2', 'NE_3', 'NE_4'))
Tumor <- subset(Tumor, annotation_monocle %in% c('LE_5', 'LE_3', 'NE'))


ncount_name <- 'nCount_Spatial.008um'
Tumor <- JoinLayers(Tumor)
Tumor <- NormalizeData(Tumor) # <----------------------
Tumor <- FindVariableFeatures(Tumor, nfeatures = 3000)




## subset ---------
set.seed(1234)
meta_Tumor_subset <-
  Tumor[['annotation_monocle']] %>%
  rownames_to_column("cell_id") %>%
  group_by(annotation_monocle) %>%
  slice_sample(n = 1500)


## cds ---------
expr_dat_monocle <- as(LayerData(JoinLayers(Tumor), "data"), "sparseMatrix") # <----------------------
p_data_monocle <- Tumor[[]]
f_data_monocle <- data.frame(
  gene_short_name = rownames(Tumor),
  row.names       = rownames(Tumor)
)


cds <- 
  newCellDataSet(
    expr_dat_monocle[, meta_Tumor_subset$cell_id],
    phenoData   = new("AnnotatedDataFrame", p_data_monocle[meta_Tumor_subset$cell_id, ]),
    featureData = new("AnnotatedDataFrame", f_data_monocle),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)




##使用seurat选择的高变基因⚠️ ------
express_genes_monocle <- VariableFeatures(Tumor)
cds <- setOrderingFilter(cds, express_genes_monocle)
plot_ordering_genes(cds)

cds <- reduceDimension(cds)
cds <- orderCells(cds)

monocle_out_dir2 <- file.path(monocle_out_dir, 'seuratHVGs')
dir.create(monocle_out_dir2)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_State_seuratHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="annotation_monocle", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_annotation_seuratHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_Pseudotime_seuratHVGs.pdf")), height = 6, width = 8)
pdf(file.path(monocle_out_dir2, "barplot of cell props group by State - annotation_seuratHVGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, annotation_monocle ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%

  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_monocle ))
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - State - split by annotation_seuratHVGs.pdf"), width = 12, height = 6)
local({
  Biobase::pData(cds)$annotation_monocle <- Tumor[['annotation_monocle']][rownames(Biobase:::pData(cds)), ]

  # cds2 <- cds[, !is.na(pData(cds)$annotation_monocle)]

  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ annotation_monocle) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_monocle") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation_monocle)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_monocle", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)
dev.off()

saveRDS(cds, file = file.path(monocle_out_dir, "cds_seuratHVGs.rds"))




##使用clusters差异表达基因 ---------------
Idents(Tumor) <- "annotation_monocle"
deg.cluster <- FindAllMarkers(JoinLayers(Tumor))
express_genes_monocle <- subset(deg.cluster, p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds, express_genes_monocle)
plot_ordering_genes(cds)

cds <- reduceDimension(cds)
cds <- orderCells(cds)

monocle_out_dir2 <- file.path(monocle_out_dir, 'clusterDEGs')
dir.create(monocle_out_dir2)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_State_clusterDEGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="annotation_monocle", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_annotation_clusterDEGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_Pseudotime_clusterDEGs.pdf")), height = 6, width = 8)
pdf(file.path(monocle_out_dir2, "barplot of cell props group by State - annotation_clusterDEGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, annotation_monocle ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%

  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_monocle ))
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - State - split by annotation_clusterDEGs.pdf"), width = 12, height = 6)
local({
  Biobase::pData(cds)$annotation_monocle <- Tumor[['annotation_monocle']][rownames(Biobase:::pData(cds)), ]

  # cds2 <- cds[, !is.na(pData(cds)$annotation_monocle)]

  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ annotation_monocle) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_monocle") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation_monocle)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_monocle", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)
dev.off()

saveRDS(cds, file = file.path(monocle_out_dir, "cds_clusterDEGs.rds"))




##使用clusters差异表达基因 ---------------
Idents(Tumor) <- "annotation_monocle"
deg.cluster <- FindAllMarkers(JoinLayers(Tumor))
express_genes_monocle <- subset(deg.cluster, p_val_adj < 0.05 & avg_log2FC > 0.25 & (pct.1 - pct.2) > 0.1)$gene
cds <- setOrderingFilter(cds, express_genes_monocle)
plot_ordering_genes(cds)

cds <- reduceDimension(cds)
cds <- orderCells(cds)

monocle_out_dir2 <- file.path(monocle_out_dir, 'clusterDEGs_filtered')
dir.create(monocle_out_dir2)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_State_clusterDEGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="annotation_monocle", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_annotation_clusterDEGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_Pseudotime_clusterDEGs.pdf")), height = 6, width = 8)
pdf(file.path(monocle_out_dir2, "barplot of cell props group by State - annotation_clusterDEGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, annotation_monocle ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_monocle ))
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - State - split by annotation_clusterDEGs.pdf"), width = 12, height = 6)
local({
  Biobase::pData(cds)$annotation_monocle <- Tumor[['annotation_monocle']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$annotation_monocle)]
  
  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ annotation_monocle) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_monocle") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation_monocle)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_monocle", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)
dev.off()

saveRDS(cds, file = file.path(monocle_out_dir, "cds_clusterDEGs.rds"))




##使用monocle选择的高变基因⚠️ --------------
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
express_genes_monocle <- disp.genes
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds)
cds <- orderCells(cds)

monocle_out_dir2 <- file.path(monocle_out_dir, 'monocleHVGs_NE')
dir.create(monocle_out_dir2)

# cds <- orderCells(cds, root_state = 2)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_State_monocleHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="annotation_monocle", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_annotation_monocleHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_Pseudotime_monocleHVGs.pdf")), height = 6, width = 8)
pdf(file.path(monocle_out_dir2, "barplot of cell props group by State - annotation_monocleHVGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, annotation_monocle ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_monocle ))
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - State - split by annotation_monocleHVGs.pdf"), width = 12, height = 3)
local({
  Biobase::pData(cds)$annotation_monocle <- Tumor[['annotation_monocle']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$annotation_monocle)]
  
  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ annotation_monocle) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_monocle") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation_monocle)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_monocle", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)
dev.off()

saveRDS(cds, file = file.path(monocle_out_dir, "cds_monocleHVGs_1dispersion_NE.rds"))




##使用monocle选择的高变基因⚠️ --------------
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 1.5 * dispersion_fit)$gene_id
express_genes_monocle <- disp.genes
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds)
cds <- orderCells(cds)

monocle_out_dir2 <- file.path(monocle_out_dir, 'monocleHVGs_1.5dispersion_NE')
dir.create(monocle_out_dir2)

# cds <- orderCells(cds, root_state = 3)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_State_monocleHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="annotation_monocle", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_annotation_monocleHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
ggsave(file.path(monocle_out_dir2, str_glue("trajactory_Pseudotime_monocleHVGs.pdf")), height = 6, width = 8)
pdf(file.path(monocle_out_dir2, "barplot of cell props group by State - annotation_monocleHVGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, annotation_monocle ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = annotation_monocle ))
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - State - split by annotation_monocleHVGs.pdf"), width = 12, height = 4)
local({
  Biobase::pData(cds)$annotation_monocle <- Tumor[['annotation_monocle']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$annotation_monocle)]
  
  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ annotation_monocle) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()
pdf(file.path(monocle_out_dir2, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_monocle") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ annotation_monocle)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "annotation_monocle", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)
dev.off()

saveRDS(cds, file = file.path(monocle_out_dir, "cds_monocleHVGs_1.5dispersion_NE.rds"))




# BEAM分析 -------------
gene_for_beam <- VariableFeatures(Tumor)
BEAM_res <- BEAM(cds[gene_for_beam,], branch_point = 1, cores = 10) 

saveRDS(BEAM_res, file.path(monocle_out_dir2, "BEAM_variable.rds"))

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
write.csv(BEAM_res, file.path(monocle_out_dir2, "BEAM_variable.csv"), row.names = F)

# plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 1e-4)),],
#                             branch_point = 1, #绘制的是哪个分支
#                             num_clusters = 4, #分成几个cluster，根据需要调整
#                             cores = 5,
#                             use_gene_short_name = T,
#                             show_rownames = T) #有632个gene，太多了


BEAM_genes <- BEAM_res %>%
  dplyr::filter(qval < 0.05) %>%
  dplyr::filter(!str_detect(gene_short_name , "^MT")) %>% 
  # top_n(n = 100, dplyr::desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()

# undebug(plot_genes_branched_heatmap)
p <- 
  cds[BEAM_genes,] %>%
  # cds[c('EVX2','CLEC2A','CNGA1', 'TCP10L2','RFX6'),] %>%
  plot_genes_branched_heatmap(
    branch_point = 1, 
    num_clusters = 3, # 3, 
    show_rownames = F, 
    return_heatmap = T
  )
ggsave(file.path(monocle_out_dir2, "BEAM_heatmap_full.pdf"), p$ph_res, width = 6.5, height = 10)


write.csv(p$annotation_row, file.path(monocle_out_dir2, 'BEAM_heatmap_genes.csv'))

# cluster1 LE Tumor_0_3
plot_genes_branched_pseudotime(cds[c("PLA2G2A", "KLK2", "KLK3", "NPY")], color_by = "State")
# cluster2 NE Tumor_0_3
plot_genes_branched_pseudotime(cds[c("CALCA", "GRP", "CGA", "CHGB")], color_by = "State") + scale_y_log10()
# cluster3 
plot_genes_branched_pseudotime(cds[c("COL3A1", "MYH11", "SLC9A3R1", "DKK3", "FAM162B")], color_by = "State")

plot_genes_jitter(cds[c("AREG", "EREG")], color_by = "State")

# cluster1 LE Tumor_5_3
plot_genes_branched_pseudotime(cds[c("PLA2G2A", "KLK2", "KLK3")], color_by = "State")
# cluster2 NE Tumor_5_3
plot_genes_branched_pseudotime(cds[c("CALCA", "GRP", "CHGB")], color_by = "State") + scale_y_log10()
# cluster3 
plot_genes_branched_pseudotime(cds[c("OR51E2", "REG4", "SPOCK1", "STEAP1", "GRHL2")], color_by = "State")

## branch gene enrichment ------------------
# http://172.30.10.203:6601/graphics/plot_zoom_png?width=443&height=708
beam_genes_ls <- 
  p$ph_res$tree_row %>% cutree(k = 3) %>%
  split(names(.), .) 


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


### GO ----------------
all_enrich_beam_ls <- 
  beam_gene_id_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
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
    file.path(monocle_out_dir2, "monocle BEAM gene enrichment analysis table - full.xlsx"),
    asTable = TRUE
  )

pdf(
  file.path(monocle_out_dir2, "monocle BEAM gene enrichment analysis dotplot - full.pdf"), 
  height = 5, 
  width = 8
)
p_all_enrich_beam_ls %>%  walk(print)
dev.off()




### ROA --------------------------
# 获取Homo sapiens的Hallmark基因集
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
# 查看Hallmark基因集的前几行
head(hallmark_gene_sets)
# 提取Term2Gene和Term2Name数据框，这是enricher和GSEA函数需要的格式
hallmark_terms <- hallmark_gene_sets %>% 
  dplyr::select(gs_name, entrez_gene) # gs_name是通路名，entrez_gene是Entrez ID
# 如果你用Gene Symbol作基因ID，可以选gs_name和gene_symbol
# hallmark_terms <- hallmark_gene_sets %>% dplyr::select(gs_name, gene_symbol)


# 使用enricher函数进行ORA分析
# TERM2GENE指定基因集与基因的对应关系，TERM2NAME可选，提供基因集的完整名称
ora_result_ls <- 
  beam_gene_id_ls %>% 
  map(function(beam_gene_id){
    clusterProfiler::enricher(
      gene = beam_gene_id, # 差异基因的Entrez ID向量
      pvalueCutoff = 0.05, # p值阈值
      pAdjustMethod = "BH", # p值校正方法，常用"BH"（即FDR）
      qvalueCutoff = 0.2, # q值阈值
      minGSSize = 10, # 基因集最小大小
      maxGSSize = 500, # 基因集最大大小
      TERM2GENE = hallmark_terms[, c("gs_name", "entrez_gene")]
    ) # 重要！指定基因集
  })
  

# 查看富集结果概要
ora_result_ls$`1`

# 如果希望结果中的基因ID显示为更易读的Gene Symbol，可以转换
ora_result_readable_ls <- 
  ora_result_ls %>% 
  map(function(ora_result){
    setReadable(ora_result, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
  })
  
# 查看转换后的结果
head(ora_result_readable_ls$`1`@result)

ora_result_readable_ls %>% 
  iwalk(function(ora_result_readable, cluster){
    openxlsx::write.xlsx(ora_result_readable, file.path(monocle_out_dir2, str_glue('ora_result_cluster{cluster}.xlsx')))
  })


Hallmark_dotplot_ls <- 
  ora_result_readable_ls %>% 
  imap(function(ora_result_readable, name){
    enrichplot::dotplot(
      ora_result_readable,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(str_glue('BEAM Cluster {name}'))
  })

Hallmark_barplot_ls <- 
  ora_result_readable_ls %>% 
  imap(function(ora_result_readable, name){
    barplot(
      ora_result_readable,
      # x = "GeneRatio",
      # color = "p.adjust",
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(str_glue('BEAM Cluster {name}'))
  })


pdf(file.path(monocle_out_dir2, 'BEAM genes enrichment-Hallmark-dotplot.pdf'), height = 5, width = 8)
Hallmark_dotplot_ls %>% 
  walk(function(Hallmark_dotplot){
    print(Hallmark_dotplot)
  })
dev.off()

pdf(file.path(monocle_out_dir2, 'BEAM genes enrichment-Hallmark-barplot.pdf'), height = 5, width = 8)
Hallmark_barplot_ls %>% 
  walk(function(Hallmark_barplot){
    print(Hallmark_barplot)
  })
dev.off()


#cds <- orderCells(cds, root_state = 4)



