library(Seurat)
library(tidyverse)
library(monocle)
source("monocle fix bugs.R")
source('color_pallete.R')

all_CAFs <- readRDS('Fibroblast_subtype/all_CAFs3.rds')
all_CAFs$subtype_reanno2 %>% table()
CAF <- subset(all_CAFs, subtype_reanno2 %in% c("MMP14_eCAFs", "iCAFs", "LOXL2_eCAFs", "vCAFs", "rCAFs"))
CAF$clinical_type <- str_replace_all(CAF$sample, c(
  'P1$' = 'T1', 
  'P2' = 'T1', 
  'P3' = 'T1', 
  'P4' = 'T2', 
  'P5' = 'T2', 
  'P6' = 'T3', 
  'P7' = 'T3', 
  'P8' = 'T4', 
  'P9' = 'T4', 
  'P10' = 'T4'
))
table(CAF$clinical_type, CAF$sample)
CAF$tumor_type <- str_replace_all(CAF$sample, c(
  'P1$' = 'TA', 
  'P2' = 'TA', 
  'P3' = 'SRCC', 
  'P4' = 'TA', 
  'P5' = 'SRCC', 
  'P6' = 'TA', 
  'P7' = 'TA', 
  'P8' = 'SRCC', 
  'P9' = 'SRCC', 
  'P10' = 'TA'
))

CAF <- CAF[, CAF$nCount_RNA > 100]


# DefaultAssay(CAF) <- "Spatial.008um" 
# CAF[['sketch']]   <- NULL 

# set default ------- 

monocle_out_dir <- "monocle_result_caf" 
dir.create(monocle_out_dir)


# dr --------
# monocle_out_dir <- "monocle_result_macrophage/SRCC" 
# dir.create(monocle_out_dir) 
# SRCC_macro2 <- FindVariableFeatures(SRCC_macro2, nfeatures = 3000)
# SRCC_macro2 <- ScaleData(
#   SRCC_macro2,
#   vars.to.regress = c(
#     # "orig.ident",
#     'nCount_RNA'
#     # 'percent.mt'
#   )
# )
# 
# 
# SRCC_macro2 <- RunPCA(SRCC_macro2, features = VariableFeatures(SRCC_macro2))
# 
# pdf(file = file.path(monocle_out_dir, 'pca plot.pdf'), height = 6, width = 8)
# ElbowPlot(SRCC_macro2, ndims = 30)
# dev.off()
# 
# 
# epi_dims <- 15
# SRCC_macro2 <- RunUMAP(SRCC_macro2, dims = 1:epi_dims, return.model = TRUE, assay  = 'RNA')
# 
# 
# pdf(file.path(monocle_out_dir, str_glue("Dimplot dim{epi_dims}.pdf")), height = 6, width = 8)
# DimPlot(SRCC_macro2, group.by = "clinical_type", reduction = "umap")
# DimPlot(SRCC_macro2, group.by = "sample", reduction = "umap")
# # DimPlot(SRCC_macro2, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(SRCC_macro2, group.by = "annotation_full_re", label = TRUE, reduction = "umap") + NoLegend()
# dev.off()
# 
# 
# SRCC_macro2 <- SRCC_macro2 %>%
#   harmony::RunHarmony(
#     dims.use = 1:epi_dims,
#     group.by.vars = "sample",
#     reduction.save = "harmony"
#   )
# 
# SRCC_macro2 <- RunUMAP(
#   SRCC_macro2,
#   dims = 1:epi_dims,
#   return.model = TRUE,
#   reduction = "harmony",
#   reduction.name = "harmonyumap",
#   reduction.key  = "harmonyUMAP_"
# )
# 
# pdf(file.path(monocle_out_dir, str_glue("Dimplot dim{epi_dims} - harmony.pdf")), height = 6, width = 8)
# DimPlot(SRCC_macro2, group.by = "clinical_type", reduction = "harmonyumap")
# DimPlot(SRCC_macro2, group.by = "sample", reduction = "harmonyumap")
# # DimPlot(SRCC_macro2, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(SRCC_macro2, group.by = "annotation_full_re", label = TRUE, reduction = "harmonyumap")
# dev.off()


# ## dotplot ---------
# marker_for_epi <- c(
#   
#   Endothelial = "PECAM1",
#   Endothelial = "VWF",
#   Endothelial = "AQP1",
#   
#   Epithelial  = "EPCAM", 
#   Epithelial  = "CEACAM5", 
#   Epithelial  = "CEACAM6", 
#   Epithelial  = "KRT5",
#   
#   Enterocyte  = "SLC26A3",
#   Enterocyte  = "CA1",
#   Enterocyte  = "GUCA2A",
#   Enterocyte  = "GUCA2B",
#   
#   # Goblet = "PIGR",
#   Goblet = "MUC2",
#   Goblet = "TFF3",
#   
#   Enteroendocrine = "CHGA",
#   
#   # "Fibroblasts" = "SFRP4", 
#   "Fibroblasts" = "PDGFRA", 
#   "Fibroblasts" = "CCDC80",
#   "Fibroblasts" = "COL1A1",
#   "Fibroblasts" = "DCN",
#   "Fibroblasts" = "ACTA2",
# 
#   "Pericytes" = "PDGFRB",
#   "Pericytes" = "RGS5",
#   "Pericytes" = "MCAM",
#   
#   
# 
#   "Smooth muscle cells" = "PLN",
#   "Smooth muscle cells" = "MYLK",
#   "Smooth muscle cells" = "MYH11",
#   "Smooth muscle cells" = "TAGLN",
# 
# 
#   
#   Immune = 'PTPRC',
#   
#   "B cells" = "CD19",
#   "B cells" = "CD79A",
#   "B cells" = "MS4A1",
#   "Plasma" = "TNFRSF17",
#   "Plasma" = "JCHAIN",
#   "Plasma" = "MZB1",
#   "Plasma" = "XBP1",
# 
#   
#   "T cells" = "CD3E",
#   "T cells" = "CD3G",
#   "T cells" = "CD4",
#   "T cells" = "CD8A",
#   
#   "NK" = "GNLY",
#   "NK" = "NKG7",
#   "NK" = "NCAM1",
#   
#   "Proliferating" = "MKI67",
#   "Proliferating" = "TOP2A",
#   
#   "Myeloid" = "ITGAM",
#   "Myeloid" = "LYZ",
#   "Myeloid" = "CD74",
#   "Macrophage" = "CD68",
#   "Macrophage" = "CD163",
#   "Macrophage" = "C1QA",
#   "Macrophage" = "C1QC",
#   
#   "Mast cell" = "CPA3",
#   "Mast cell" = "CTSG",
#   "Mast cell" = "TPSAB1",
#   "Mast cell" = "MS4A2",
#   "Mast cell" = "KIT",
#   
#   
#   "Monocytes" = "CD14", 
#   "Monocytes" = "FCGR3A",
#   
#   "Dendritic cell" = "CD83",
#   'Dendritic cell' = "LILRA4",  
#   
#   'Dendritic cell' = "CLEC9A",
#   'Dendritic cell' = "IDO1",
#   
#   'Dendritic cell' = "CD1C",
#   'Dendritic cell' = "FCER1A",
#   
#   'Neutrophi' = "S100A8",
#   'Neutrophi' = "S100A9",
#   'Neutrophi' = "CSF3R",
#   'Neutrophi' = "FCGR3B",
#   
#   Gial = "MPZ"
# )
# 
# 
# 
# 
# theme_90angle <-
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   )
# 
# Idents(SRCC_macro2) <- "annotation_full_re"
# pdf(file.path(monocle_out_dir, str_glue("dotplot dims{epi_dims} main_marker.pdf")), width = 25, height = 6)
# DotPlot(SRCC_macro2, group.by = "annotation_full_re", features = marker_for_epi, dot.scale = 10) + theme_90angle
# dev.off()
# 
# # saveRDS(SRCC_macro2, file.path(monocle_out_dir, "SRCC_macro2_dr.rds"))
# 
# # SRCC_macro2 <- readRDS("epi_subcluster/SRCC_macro2.rds")
# # SRCC_macro2 <- readRDS("epi_subcluster_cycle2/SRCC_macro2_15.rds")


# monocle out -------

# DimPlot(SRCC_macro2, reduction = 'harmonyumap') + scale_color_manual(values = col_plan_1)

VariableFeatures(CAF)

# subset
set.seed(1234)
meta_CAF_subset <-
  CAF[['clinical_type']] %>%
  rownames_to_column("cell_id") %>%
  group_by(clinical_type) %>%
  slice_sample(n = 5000)
# meta_CAF_subset <- CAF[['clinical_type']] %>% rownames_to_column("cell_id")

# cds ---------
expr_dat_monocle <- as(LayerData(JoinLayers(CAF), "counts"), "sparseMatrix")
p_data_monocle <- CAF[[]]
f_data_monocle <- data.frame(
  gene_short_name = rownames(CAF),
  row.names       = rownames(CAF)
)


cds <- 
  newCellDataSet(
    expr_dat_monocle[, meta_CAF_subset$cell_id],
    phenoData   = new("AnnotatedDataFrame", p_data_monocle[meta_CAF_subset$cell_id, ]),
    featureData = new("AnnotatedDataFrame", f_data_monocle),
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


# cds <- clusterCells(cds)


## 使用seurat选择的高变基因⚠️ ------
CAF <- FindVariableFeatures(CAF, nfeatures = 2000)
express_genes_monocle <- VariableFeatures(CAF)
cds <- setOrderingFilter(cds, express_genes_monocle)
plot_ordering_genes(cds)

## 使用clusters差异表达基因
# Idents(CAF) <- "clinical_type"
# deg.cluster <- FindAllMarkers(JoinLayers(CAF))
# express_genes_monocle <- subset(deg.cluster, p_val_adj<0.05)$gene
# cds <- setOrderingFilter(cds, express_genes_monocle)
# plot_ordering_genes(cds)

## 使用monocle选择的高变基因⚠️
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.001 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
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
# undebug(orderCells)
cds <- orderCells(cds)

# Error:
#   ! `nei()` was deprecated in igraph 2.1.0 and is now defunct.
# ℹ Please use `.nei()` instead.
# install.packages('https://cran.r-project.org/src/contrib/Archive/igraph/igraph_2.0.3.tar.gz', type = "source", repos = NULL)
 

plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir, str_glue("trajactory_State_seuratHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="subtype_reanno2", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir, str_glue("trajactory_annotation_seuratHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds,color_by="clinical_type", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
ggsave(file.path(monocle_out_dir, str_glue("trajactory_clinical_type_seuratHVGs.pdf")), height = 6, width = 8)
plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
ggsave(file.path(monocle_out_dir, str_glue("trajactory_Pseudotime_seuratHVGs.pdf")), height = 6, width = 8)
pdf(file.path(monocle_out_dir, "barplot of cell props group by State - annotation_seuratHVGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, subtype_reanno2 ) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%

  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = subtype_reanno2 ))
dev.off()
pdf(file.path(monocle_out_dir, "barplot of cell props group by State - clinical_type_seuratHVGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, clinical_type) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%

  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = clinical_type))
dev.off()
pdf(file.path(monocle_out_dir, "barplot of cell props group by State - tumor_type_seuratHVGs.pdf"), width = 6, height = 5)
cds %>%
  Biobase::pData() %>%
  group_by(State, tumor_type) %>%
  summarise(count = n()) %>%
  group_by(State) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = tumor_type))
dev.off()
pdf(file.path(monocle_out_dir, "cell trajectory - State - split by annotation_seuratHVGs.pdf"), width = 12, height = 4)
local({
  Biobase::pData(cds)$subtype_reanno2 <- CAF[['subtype_reanno2']][rownames(Biobase:::pData(cds)), ]

  # cds2 <- cds[, !is.na(pData(cds)$subtype_reanno2)]

  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ subtype_reanno2) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()
pdf(file.path(monocle_out_dir, "cell trajectory - State - split by clinical_type_seuratHVGs.pdf"), width = 12, height = 4)
local({
  Biobase::pData(cds)$clinical_type <- CAF[['clinical_type']][rownames(Biobase:::pData(cds)), ]

  # cds2 <- cds[, !is.na(pData(cds)$clinical_type)]

  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ clinical_type) +
    theme(legend.position = "right")+
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()




# plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_State_clusterDEGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds,color_by="subtype_reanno2", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_annotation_clusterDEGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds,color_by="clinical_type", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_clinical_type_clusterDEGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_Pseudotime_clusterDEGs.pdf")), height = 6, width = 8)
# pdf(file.path(monocle_out_dir, "barplot of cell props group by State - annotation_clusterDEGs.pdf"), width = 6, height = 5)
# cds %>%
#   Biobase::pData() %>%
#   group_by(State, subtype_reanno2 ) %>%
#   summarise(count = n()) %>%
#   group_by(State) %>%
#   mutate(freq = count / sum(count)) %>%
# 
#   ggplot(aes(x = State, y = freq)) +
#   geom_col(aes(fill = subtype_reanno2 ))
# dev.off()
# pdf(file.path(monocle_out_dir, "barplot of cell props group by State - clinical_type_clusterDEGs.pdf"), width = 6, height = 5)
# cds %>%
#   Biobase::pData() %>%
#   group_by(State, clinical_type) %>%
#   summarise(count = n()) %>%
#   group_by(State) %>%
#   mutate(freq = count / sum(count)) %>%
# 
#   ggplot(aes(x = State, y = freq)) +
#   geom_col(aes(fill = clinical_type))
# dev.off()
# pdf(file.path(monocle_out_dir, "barplot of cell props group by State - tumor_type_clusterDEGs.pdf"), width = 6, height = 5)
# cds %>%
#   Biobase::pData() %>%
#   group_by(State, tumor_type) %>%
#   summarise(count = n()) %>%
#   group_by(State) %>%
#   mutate(freq = count / sum(count)) %>%
#   
#   ggplot(aes(x = State, y = freq)) +
#   geom_col(aes(fill = tumor_type))
# dev.off()
# pdf(file.path(monocle_out_dir, "cell trajectory - State - split by annotation_clusterDEGs.pdf"), width = 12, height = 4)
# local({
#   Biobase::pData(cds)$subtype_reanno2 <- CAF[['subtype_reanno2']][rownames(Biobase:::pData(cds)), ]
# 
#   # cds2 <- cds[, !is.na(pData(cds)$subtype_reanno2)]
# 
#   plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
#     scale_color_manual(values = col_plan_2) +
#     facet_wrap(~ subtype_reanno2) +
#     theme(legend.position = "right")+
#     guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
# })
# dev.off()
# pdf(file.path(monocle_out_dir, "cell trajectory - State - split by clinical_type_clusterDEGs.pdf"), width = 12, height = 4)
# local({
#   Biobase::pData(cds)$clinical_type <- CAF[['clinical_type']][rownames(Biobase:::pData(cds)), ]
# 
#   # cds2 <- cds[, !is.na(pData(cds)$clinical_type)]
# 
#   plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
#     scale_color_manual(values = col_plan_2) +
#     facet_wrap(~ clinical_type) +
#     theme(legend.position = "right")+
#     guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
# })
# dev.off()




# plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_State_monocleHVGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds,color_by="subtype_reanno2", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_annotation_monocleHVGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds,color_by="clinical_type", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_clinical_type_monocleHVGs.pdf")), height = 6, width = 8)
# plot_cell_trajectory(cds, color_by="Pseudotime", size=1, show_backbone=TRUE)
# ggsave(file.path(monocle_out_dir, str_glue("trajactory_Pseudotime_monocleHVGs.pdf")), height = 6, width = 8)
# pdf(file.path(monocle_out_dir, "barplot of cell props group by State - annotation_monocleHVGs.pdf"), width = 6, height = 5)
# cds %>%
#   Biobase::pData() %>%
#   group_by(State, subtype_reanno2 ) %>%
#   summarise(count = n()) %>%
#   group_by(State) %>%
#   mutate(freq = count / sum(count)) %>%
# 
#   ggplot(aes(x = State, y = freq)) +
#   geom_col(aes(fill = subtype_reanno2 ))
# dev.off()
# pdf(file.path(monocle_out_dir, "barplot of cell props group by State - clinical_type_monocleHVGs.pdf"), width = 6, height = 5)
# cds %>%
#   Biobase::pData() %>%
#   group_by(State, clinical_type) %>%
#   summarise(count = n()) %>%
#   group_by(State) %>%
#   mutate(freq = count / sum(count)) %>%
# 
#   ggplot(aes(x = State, y = freq)) +
#   geom_col(aes(fill = clinical_type))
# dev.off()
# pdf(file.path(monocle_out_dir, "barplot of cell props group by State - tumor_type_monocleHVGs.pdf"), width = 6, height = 5)
# cds %>%
#   Biobase::pData() %>%
#   group_by(State, tumor_type) %>%
#   summarise(count = n()) %>%
#   group_by(State) %>%
#   mutate(freq = count / sum(count)) %>%
#   
#   ggplot(aes(x = State, y = freq)) +
#   geom_col(aes(fill = tumor_type))
# dev.off()
# pdf(file.path(monocle_out_dir, "cell trajectory - State - split by annotation_monocleHVGs.pdf"), width = 12, height = 4)
# local({
#   Biobase::pData(cds)$subtype_reanno2 <- CAF[['subtype_reanno2']][rownames(Biobase:::pData(cds)), ]
# 
#   # cds2 <- cds[, !is.na(pData(cds)$subtype_reanno2)]
# 
#   plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
#     scale_color_manual(values = col_plan_2) +
#     facet_wrap(~ subtype_reanno2) +
#     theme(legend.position = "right")+
#     guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
# })
# dev.off()
# pdf(file.path(monocle_out_dir, "cell trajectory - State - split by clinical_type_monocleHVGs.pdf"), width = 12, height = 4)
# local({
#   Biobase::pData(cds)$clinical_type <- CAF[['clinical_type']][rownames(Biobase:::pData(cds)), ]
# 
#   # cds2 <- cds[, !is.na(pData(cds)$clinical_type)]
# 
#   plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
#     scale_color_manual(values = col_plan_2) +
#     facet_wrap(~ clinical_type) +
#     theme(legend.position = "right")+
#     guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
# })
# dev.off()




p_col_state <- cds %>% 
  Biobase::pData() %>%
  dplyr::select(clinical_type, State) %>%
  group_by(clinical_type, State) %>%
  mutate(clinical_type = as.character(clinical_type)) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = clinical_type, y = freq)) +
  geom_col(aes(fill = State)) +
  # scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
p_col_annotation <- cds %>% 
  Biobase::pData() %>%
  dplyr::select(clinical_type, State) %>%
  group_by(State, clinical_type) %>%
  mutate(clinical_type = as.character(clinical_type)) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = State, y = freq)) +
  geom_col(aes(fill = clinical_type)) +
  # scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  ggpubr::theme_pubr() +
  theme(
    legend.position = "right"
  )
p_col_state + p_col_annotation





cds <- orderCells(cds, root_state = 4)

# saveRDS(cds, file = file.path(monocle_out_dir, "cds_seuratHVGs.rds"))
saveRDS(cds, file = file.path(monocle_out_dir, "cds_clusterDEGs.rds"))
# saveRDS(cds, file = file.path(monocle_out_dir, "cds_monocleHVGs.rds"))




# plot ----------
# some ploting
# cds <- orderCells(cds, root_state = 4)
# plot_cell_trajectory(cds, color_by="Pseudotime", size=1,show_backbone=TRUE) 
# plot_cell_trajectory(cds, color_by="subtype_reanno2", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# plot_cell_trajectory(cds, color_by="clinical_type", size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_1)
# plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone=TRUE) + scale_color_manual(values = col_plan_2)
# plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "subtype_reanno2")
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "clinical_type", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)



# p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
# p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
# p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
# 
# p1=plot_cell_trajectory(cds, color_by = "CCL5")  + scale_color_gsea()

#monocle_out_dir2 <- file.path(monocle_out_dir, 'seuratHVGs')
monocle_out_dir2 <- file.path(monocle_out_dir, 'clusterDEGs')
#monocle_out_dir2 <- file.path(monocle_out_dir, 'monocleHVGs')
dir.create(monocle_out_dir2)

pdf(file.path(monocle_out_dir2, "cell trajector.pdf"), width = 8, height = 4)

  plot_cell_trajectory(cds,color_by="Pseudotime", cell_size =0.1, show_backbone=TRUE) + theme(legend.position = "right")
  
  plot_cell_trajectory(cds,color_by="clinical_type", cell_size =0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    theme(legend.position = "right")
  
  plot_cell_trajectory(cds,color_by="subtype_reanno2", cell_size =0.1, show_backbone=TRUE) + 
    scale_color_manual(values = col_plan_2) +
    theme(legend.position = "right") + 
    guides(color = guide_legend(ncol = 2,override.aes = list(size = 3)))
  
  plot_cell_trajectory(cds, color_by = "State",cell_size =1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    theme(legend.position = "right")

dev.off()


pdf(file.path(monocle_out_dir2, "cell trajectory - State - split by clinical_type.pdf"), width = 12, height = 4)
local({
  Biobase::pData(cds)$clinical_type <- all_data[['clinical_type']][rownames(Biobase:::pData(cds)), ]

  # cds2 <- cds[, !is.na(pData(cds)$clinical_type)]
  
  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ clinical_type) +
    theme(legend.position = "right")+ 
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()

pdf(file.path(monocle_out_dir2, "cell trajectory - State - split by annotation.pdf"), width = 12, height = 4)
local({
  Biobase::pData(cds)$subtype_reanno2 <- all_data[['subtype_reanno2']][rownames(Biobase:::pData(cds)), ]
  
  # cds2 <- cds[, !is.na(pData(cds)$subtype_reanno2)]
  
  plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1,show_backbone=TRUE) +
    scale_color_manual(values = col_plan_2) +
    facet_wrap(~ subtype_reanno2) +
    theme(legend.position = "right")+ 
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
})
dev.off()



pdf(file.path(monocle_out_dir2, "cell trajectory - tree.pdf"), width = 12, height = 8)
plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "subtype_reanno2") +
  scale_color_manual(values = col_plan_1) +
  facet_wrap(~ subtype_reanno2)
plot_complex_cell_trajectory(cds, x = 1, y = 2, color_by = "subtype_reanno2", cell_size = 0.5)+ scale_color_manual(values = col_plan_2)
dev.off()







# pdf(file.path(monocle_out_dir2, "Dimplot of sc_epi - clinical_type.pdf"), width = 8, height = 5)
# DimPlot(
#   CAF,
#   group.by = "annotation_full"
# )
# 
# DimPlot(
#   CAF,
#   group.by = "clinical_type"
# )
# dev.off()


# cell barplot 
p_col_cell_state <-
  cds %>% 
  Biobase::pData() %>%
  dplyr::select(annotation_full, State) %>%
  group_by(annotation_full, State) %>%
  mutate(annotation_full = as.character(annotation_full)) %>%
  # mutate(annotation_full = factor(annotation_full, levels = c(
  #   'S1 fibroblasts', 
  #   'Peri-cryptal S2 fibroblasts',
  #   'Immediate early-response (IER) fibroblasts',
  #   'S3 fibroblasts',
  #   'Myofibroblasts'
  # ))) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation_full, y = freq)) +
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

## TFF3 cells -> sub cluster ----
sc_epi_tff3 <- readRDS(file.path("epi_subcluster_tff3 goblet", "sc_epi.rds"))
# sc_epi_tff3 <- sc_epi

meta_annotation_full_df <-
  Biobase::pData(cds) %>%
  dplyr::select(annotation_full) %>%
  rownames_to_column("cell_id") %>%
  left_join(
    rownames_to_column(sc_epi_tff3[["annotation"]], "cell_id"),
    by = "cell_id"
  ) %>%
  mutate(annotation_full2 = ifelse(
    annotation_full == "TFF3 Goblet",
    as.character(annotation),
    as.character(annotation_full)
  )) %>%
  mutate(annotation_full2 = annotation_full2 %>% str_replace_na("TFF3 Goblet"))

Biobase::pData(cds)$annotation_full2 <- meta_annotation_full_df$annotation_full2

# export
pdf(file.path(monocle_out_dir, "cell trajectory - tff3 goblet subcluster - dim6.pdf"), width = 8, height = 4)

  plot_cell_trajectory(cds,color_by="annotation_full2", cell_size =0.1, show_backbone=TRUE) + 
    scale_color_manual(values = col_plan_2) +
    theme(legend.position = "right") + 
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))
  
  plot_complex_cell_trajectory(cds, x = 1, cell_size = 0.1, y = 2, color_by = "annotation_full2") +
    scale_color_manual(values = col_plan_1) +
    facet_wrap(~ annotation_full2) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  p_col_cell_state2 <-
    cds %>% 
    Biobase::pData() %>%
    dplyr::select(annotation_full2, State) %>%
    group_by(annotation_full2, State) %>%
    mutate(annotation_full2 = as.character(annotation_full2)) %>%
  
    summarise(count = n()) %>%
    mutate(freq = count/sum(count)) %>%
    
    ggplot(aes(x = annotation_full2, y = freq)) +
    geom_col(aes(fill = State)) +
    scale_fill_manual(values = col_plan_2) +
    scale_y_continuous(labels = scales::percent) +
    coord_flip() +
    labs(x = NULL, y = NULL) +
    ggpubr::theme_pubr() +
    theme(
      legend.position = "right"
    )
  
  print(p_col_cell_state2)
  
dev.off()


p_pseudotime_groupby_annotation <- 
  pData(cds) %>% 
  dplyr::filter(!is.na(annotation_full2)) %>%
  {
    levs <- summarise(group_by(., annotation_full2), vals = median(Pseudotime))
    levs_order <- arrange(levs, vals)[['annotation_full2']]
    mutate(., annotation_full2 = factor(annotation_full2, levels = levs_order))
  } %>%
  ggplot(aes(x = annotation_full2, y = Pseudotime)) +
  geom_boxplot() +
  theme_classic()

ggsave(
  p_pseudotime_groupby_annotation,
  filename = file.path(monocle_out_dir, "annotation order by pseudotime.pdf"),
  height = 4,
  width  = 8
)

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
gene_for_beam <- VariableFeatures(CAF)
BEAM_res <- BEAM(cds[gene_for_beam,], branch_point = 1, cores = 10) 

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
  # top_n(n = 100, dplyr::desc(qval)) %>%
  pull(gene_short_name) %>%
  as.character()

undebug(plot_genes_branched_heatmap)
p <- 
  cds[BEAM_genes,] %>%
  # cds[c('EVX2','CLEC2A','CNGA1', 'TCP10L2','RFX6'),] %>%
  plot_genes_branched_heatmap(
    branch_point = 1, 
    num_clusters = 4, # 3, 
    show_rownames = F, 
    return_heatmap = T
  )
ggsave(file.path(monocle_out_dir, "BEAM_heatmap_full.pdf"), p$ph_res, width = 6.5, height = 10)


plot_genes_branched_pseudotime(cds[c("CLCA1", "TFF3", "SLC26A3")], color_by = "State")
plot_genes_branched_pseudotime(cds[c("AREG", "EREG")], color_by = "State") + scale_y_log10()
plot_genes_branched_pseudotime(cds[c("MLPH", "PLA2G2A", "TNFRSF11A", "MUC4", "SLC28A2")], color_by = "State")

plot_genes_jitter(cds[c("AREG", "EREG")], color_by = "State")


# branch gene enrichment http://172.30.10.203:6601/graphics/plot_zoom_png?width=443&height=708
beam_genes_ls <- 
  p$ph_res$tree_row %>% cutree(k = 4) %>%
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


# de analysis of TFF3 Goblets of State --------
cds %>% dim
cds %>% exprs() %>% .[1:5, 1:5]

cds %>% pData() %>% head()

seu_cds <-
  CreateSeuratObject(
    counts = Biobase::exprs(cds),
    meta.data = Biobase::pData(cds)
  )

seu_cds <- seu_cds %>% NormalizeData()
Idents(seu_cds) <- "State"

seu_cds_markers_ls <- 
  2:4 %>%
  set_names(c("Tumor", "Normal", "TSA")) %>%
  map(function(state){
    sc <- seu_cds 
    sc$group2 <- str_c(sc$project, sc$group, sep = "_")
    
    sc <- sc[, sc$annotation_full == "TFF3 Goblet" & sc$State %in% as.character(2:4)]
    
    FindConservedMarkers(
      sc,
      ident.1 = state,
      ident.2 = NULL,
      grouping.var = "group2",
      min.cells.group = 10,
      logfc.threshold = 0
    )
  })

seu_cds_markers_ls <- 
  seu_cds_markers_ls %>%
  map(function(de){
    de %>%
      rownames_to_column("gene") %>%
      rowwise() %>%
      mutate(avg_log2FC = mean(c_across(ends_with("avg_log2FC")))) %>%
      mutate(pvalue = minimump_p_val)
  })

seu_cds_markers_sig_ls <-
  seu_cds_markers_ls %>%
  map(function(de){
    de %>% 
      ungroup() %>%
      dplyr::filter(avg_log2FC > 0.25) %>%
      slice_min(order_by = pvalue, n = 100)
  })

# export
seu_cds_markers_ls %>%
  openxlsx::write.xlsx(
    file = file.path(monocle_out_dir, "monocle de analysis of 3 State of TFF3 Goblet cells.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

seu_cds_markers_sig_ls %>%
  openxlsx::write.xlsx(
    file = file.path(monocle_out_dir, "significant gene of monocle de analysis of 3 State of TFF3 Goblet cells.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

expr_head_genes <- cds %>% exprs %>% rowMeans() %>% .[order(.)] %>% tail(100) %>% rev()

tff3_gene1 <- hub_df$gene_name[hub_df$module == "TFF3 Goblet1"]
tff3_gene3 <- hub_df$gene_name[hub_df$module == "TFF3 Goblet3"]

seu_cds_markers_sig_ls %>%
  map('gene') %>%
  map(~ intersect(.x, tff3_gene1))
# $Tumor
# [1] "PLA2G2A"
# 
# $Normal
# [1] "ZG16"    "IGHA1"   "MT2A"    "SELENOP" "LEFTY1"  "AGR2"   
# 
# $TSA
# [1] "OLFM4"   "PLA2G2A"

seu_cds_markers_sig_ls %>%
  map('gene') %>%
  map(~ intersect(.x, c(tff3_gene1, tff3_gene3)))
  

plot_genes_jitter(cds[c("OLFM4", "PLA2G2A", "REG4")], color_by = "State")


beam_genes_ls[as.character(3:4)] %>%
  map(~ intersect(.x, tff3_gene3)) 

plot_genes_branched_pseudotime(cds[c("PLA2G2A", "OLFM4", "REG4"), ], color_by = "State")

gene_beam_highexpr <-
  beam_genes_ls[3:4] %>%
  map(~ intersect(.x, names(expr_head_genes)))

DotPlot(
  all_data,
  features = gene_beam_highexpr,
  group.by = "annotation3"
)+
  RotatedAxis()

DotPlot(
  all_data,
  features = gene_beam_highexpr,
  group.by = "group"
) +
  RotatedAxis()

gene_beam_wgcna <-
  beam_genes_ls[3:4] %>%
  # map(~ intersect(.x, names(expr_head_genes))) %>%
  map(~ intersect(.x,c(tff3_gene1, tff3_gene3))) 

# gene_beam_wgcna <-
#   beam_genes_ls[3:4] %>%
#   map(~ intersect(.x, names(expr_head_genes))) %>%
#   map(~ intersect(.x,c(tff3_gene1, tff3_gene3))) 


# plot_genes_branched_pseudotime(
#   cds[c("IFITM3", "MYC", "YWHAZ", "ETS2", "HSP90AB1", "ZNF703", "SYNGR2", "SLC12A2", "FTL" ), ],
#   color_by = "State"
# )

all_data %>%
  DotPlot(
    features = gene_beam_wgcna$`3`,
    group.by = "annotation3"
  )
all_data %>%
  DotPlot(
    features = gene_beam_wgcna$`4`,
    group.by = "annotation3"
  )


DotPlot(
  all_data,
  features = c("APP", "CD74"),
  group.by = "annotation3"
) +
  RotatedAxis()

DotPlot(
  all_data,
  features = c("PECAM1", "CD38"),
  group.by = "annotation3"
) +
  RotatedAxis()

DotPlot(
  all_data,
  features = gene_beam_wgcna,
  group.by = "annotation3"
) +
  RotatedAxis()

DotPlot(
  all_data,
  features = gene_beam_wgcna,
  group.by = "group"
) +
  RotatedAxis()



gene_beam_wgcna

SpatialFeaturePlot(
  all_data,
  images = "TSA_233841.008um",
  features = "OLFM4"
) +
  scale_fill_gradientn(
    colours = Seurat:::SpatialColors(30),
    transform = scales::new_transform(name = "log2p1", transform = function(x)log2(x+1), inverse = function(x)2^x-1)
  )






SpatialFeaturePlot(
  all_data,
  images = "TSA_233841.008um",
  features = "OLFM4"
) +
  scale_fill_gradientn(
    colours = Seurat:::SpatialColors(30),
    transform = "log1p"
  )






coords_all_data_ls <-
  all_data@images %>% 
  map(GetTissueCoordinates)  
# range_x <- coords_all_data_ls %>% map('x') %>% map(range)


kde_ls <- 
  coords_all_data_ls %>% 
  map(function(x)(
    x %>%
      left_join(rownames_to_column(all_data[[]], "cell"), by = c('cell')) %>% 
      group_by(annotation3) %>%
      nest() %>%
      deframe()
  )) %>%
  map(function(df_ls){
    kde_df <- df_ls %>%
      map(function(x){
        kde <- MASS::kde2d(x$y, -x$x)
        data.frame(
          z = kde$z %>% c(),
          x = kde$x %>% rep(times = length(kde$y)),
          y = kde$y %>% rep(each  = length(kde$x))
        )
      }) %>%
        enframe("annotation3") %>%
        unnest(value)
    kde_df
  })
  

coords_TSA_df_ls <-
  coords_all_data_ls %>%
  map(function(x){
    x %>% left_join(rownames_to_column(all_data[[c('group', 'annotation3')]], "cell"), by = "cell")
  })

density_c1qc_ls <-
  coords_TSA_df_ls %>%
  map(function(df){
    df %>%
      # dplyr::slice_sample(prop = 0.5) %>%
      ggplot(aes(x = y, y = -x)) +
      geom_tile(aes(fill = group), height = 10, width = 10) +
      scale_fill_manual(values = c(
        Normal = "gray90", TSA = "gray75", Tumor = "gray50"
      )) +
      theme_void() +
      
      geom_tile(
        data = df %>% dplyr::filter(annotation3 == "C1QC Macrophage"),
        fill = "red",
        height = 50,
        width  = 50
      ) +
      
      geom_density2d(
        data = df %>% dplyr::filter(annotation3 == "C1QC Macrophage"),
        contour_var = "density", show.legend = FALSE
      )
  })
  

density_selenop_ls <-
  coords_TSA_df_ls %>%
  map(function(df){
    df %>%
      # dplyr::slice_sample(prop = 0.5) %>%
      ggplot(aes(x = y, y = -x)) +
      geom_tile(aes(fill = group), height = 10, width = 10) +
      scale_fill_manual(values = c(
        Normal = "gray90", TSA = "gray75", Tumor = "gray50"
      )) +
      theme_void() +
      
      geom_tile(
        data = df %>% dplyr::filter(annotation3 == "SELENOP Macrophage" ),
        fill = "red",
        height = 50,
        width  = 50
      ) +
      
      geom_density2d(
        data = df %>% dplyr::filter(annotation3 == "SELENOP Macrophage" ),
        contour_var = "density", show.legend = FALSE
      )
  })
  

  
  




















kde_ls$TSA_233841.008um %>% 
  dplyr::filter(annotation3 %in% c("C1QC Macrophage")) %>% # , "SELENOP Macrophage")) %>%
  ggplot(aes(x = x, y = y, z = z)) +
  geom_contour() +
  theme_void()



kde_ls$TSA_233841.008um %>% 
  dplyr::filter(annotation3 %in% c("C1QC Macrophage", "SELENOP Macrophage", "Enterocyte" ) ) %>% 
  ggplot(aes(x = x, y = y, z = z, color = annotation3)) +
  geom_contour() +
  theme_void() +
  facet_wrap(~ annotation3)











coords_all_data_ls %>% .[1] %>% 
  map(function(x){
    x %>%
      left_join(rownames_to_column(all_data[[]], "cell"), by = c('cell')) %>%
      
      dplyr::filter(annotation3 %in% c("C1QC Macrophage") ) %>% 
      
      ggplot(aes(x = y, y = -x)) +
      facet_wrap(~ project) +
      # stat_density_2d(
      #   geom="raster",
      #   aes(fill = after_stat(ndensity)),
      #   alpha=0.8,
      #   contour  = FALSE,
      #   n = 100,
      #   show.legend = T
      # )
      # geom_density2d(contour_var = "density", show.legend = T )
      geom_density2d(contour_var = "ndensity", show.legend = T )
      # stat_density_2d(geom = "point", aes(size = after_stat(density)), n = 20, contour = T)
  })


coords_all_data_ls$TSA_233841.008um %>%
  slice_sample(prop = 0.1) %>%
  ggplot(aes(x = y, y = -x)) +
  geom_point()




