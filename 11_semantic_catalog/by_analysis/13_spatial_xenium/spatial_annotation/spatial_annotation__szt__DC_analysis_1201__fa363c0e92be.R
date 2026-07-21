library(Seurat)
library(tidyverse)
library(BPCells)
#All_T subtype-----------------------
T_subtype1 <- readRDS("T_subtype/T_subtype.rds")
T_subtype2 <- readRDS("T_subtype/T_subtype2.rds")
T_subtype3 <- readRDS("T_subtype/third/T_subtype.rds")
dim(T_subtype3)

table(T_subtype1$subtype)
table(T_subtype2$subtype)
table(T_subtype3$subtype)

pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker - T_subtype1.pdf")), width = 20, height = 6)
DotPlot(T_subtype1, group.by = "subtype", features = major_marker_ls)+RotatedAxis()
DotPlot(T_subtype1, group.by = "subtype", features = major_marker_ls_given)+RotatedAxis()
DotPlot(T_subtype1, group.by = "subtype", features = major_marker_ls_CD8)+RotatedAxis()
DotPlot(T_subtype1, group.by = "RNA_snn_res.1", features = major_marker_ls)+RotatedAxis()
DotPlot(T_subtype1, group.by = "RNA_snn_res.1", features = major_marker_ls_given)+RotatedAxis()
DotPlot(T_subtype1, group.by = "RNA_snn_res.1", features = major_marker_ls_CD8)+RotatedAxis()
dev.off()

pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker - T_subtype2.pdf")), width = 20, height = 6)
DotPlot(T_subtype2, group.by = "subtype", features = major_marker_ls)+RotatedAxis()
DotPlot(T_subtype2, group.by = "subtype", features = major_marker_ls_given)+RotatedAxis()
DotPlot(T_subtype2, group.by = "subtype", features = major_marker_ls_CD8)+RotatedAxis()
DotPlot(T_subtype2, group.by = "RNA_snn_res.1", features = major_marker_ls)+RotatedAxis()
DotPlot(T_subtype2, group.by = "RNA_snn_res.1", features = major_marker_ls_given)+RotatedAxis()
DotPlot(T_subtype2, group.by = "RNA_snn_res.1", features = major_marker_ls_CD8)+RotatedAxis()
dev.off()

pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker - T_subtype3.pdf")), width = 20, height = 6)
DotPlot(T_subtype3, group.by = "subtype", features = major_marker_ls)+RotatedAxis()
DotPlot(T_subtype3, group.by = "subtype", features = major_marker_ls_given)+RotatedAxis()
DotPlot(T_subtype3, group.by = "subtype", features = major_marker_ls_CD8)+RotatedAxis()
DotPlot(T_subtype3, group.by = "RNA_snn_res.1", features = major_marker_ls)+RotatedAxis()
DotPlot(T_subtype3, group.by = "RNA_snn_res.1", features = major_marker_ls_given)+RotatedAxis()
DotPlot(T_subtype3, group.by = "RNA_snn_res.1", features = major_marker_ls_CD8)+RotatedAxis()
dev.off()


out_dir<-'T_subtype/Fourth'
dir.create(out_dir)


### dotplot -------
major_marker_ls <- c(
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  Endothelial = "CLDN5",
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  Epithelial  = "CEACAM6",
  Epithelial  = "CDH1",
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  "Smooth muscle cells" = "PLN",
  "Smooth muscle cells" = "MYLK",
  "Smooth muscle cells" = "MYH11",
  Immune = 'PTPRC',
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "CD3D",
  "T cells" = "CD4",
  "T cells" = "CD8A",
  "T cells" = "TRAC",
  "T cells" = "IL7R",
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  "Proliferating" = "MKI67",
  "Proliferating" = "TOP2A",
  "Myeloid" = "ITGAM",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC",
  "Macrophage" = 'CSF1R',
  "Macrophage" = "CXCL10",
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  "Mast cell" = "MS4A2",
  "Mast cell" = "KIT",
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  "Monocytes" = "FCN1",
  'Dendritic cell' = "LILRA4",  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "CD83",
  'Dendritic cell' = "ITGAX",
  'Dendritic cell' = "XCR1",
  'Dendritic cell' = "IDO1",
  'Neutrophi' = "S100A8",
  'Neutrophi' = "S100A9",
  'Neutrophi' = "CSF3R",
  'Neutrophi' = "FCGR3B"
)

major_marker_ls %>%
  str_to_title() %>%
  {names(.) <- names(major_marker_ls); .} -> major_marker_ls


# major_marker_ls <- c(
#   'CD8'='CD8A',
#   'CD8'='CD8B',
#   'CD8'='GZMB',
#   'CD8'='GZMA',
#   'CD8'='GZMK',
#   'CD8'='GZMH',
#   'CD8'='CCL5',
#   'Naive'='CCR6',
#   'Naive'='LEF1',
#   'Naive'='TCF7',
#   'Naive'='SELL',
#   'Naive'='CCR7',
#   'Naive'='YBX3',
#   'Naive'='TGFB1',
#   'Naive'='CCl5',
#   'Naive'='CXCR3',
#   'Memory'='GPR183',
#   'Memory'='IL7R',
#   'Memory'='FOS',
#   'Memory'='JUN',
#   'Memory'='KLRB1',
#   'TRM'='CD89',
#   'TRM'='CXCR8',
#   # 'TEM'='IFNG',
#   'TEM'='OCL4',
#   'TEM'='OCL5',
#   'TEM'='JAK3',
#   # 'TCM'='IDO1',
#   # 'TCM'='SKIV2L',
#   # 'TCM'='TRAF2',
#   # 'TCM'='SLC2A1',
#   # 'TCM'='CD40LG',
#   'TCM'='DUSP2',
#   'TCM'='LYAR',
#   'TCM'='ZNF683',
#   'TCM'='IL32',
#   'Effector'='FGFBP2',
#   'Effector'='TYROBP',
#   'Effector'='KLRC2',
#   'Effector'='CX3CR1',
#   'Effector'='FCGR3A',
#   'Effector'='GNLY',
#   'cytotoxic'='XCL1',
#   'cytotoxic'='IFNG',
#   'cytotoxic'='CCL4',
#   # 'cytotoxic'='CCL5',
#   'Exhausted'='CXCL13',
#   'Exhausted'='LAYN',
#   'Exhausted'='LAG3',
#   'Exhausted'='HAVCR2',
#   'Exhausted'='TNFRSF9',
#   'Exhausted'='CTLA4',
#   'Exhausted'='TIGIT',
#   'Exhausted'='TOX',
#   'Exhausted'='PDCD1',
#   'CD4'='FOXP3',
#   'CD4'='CD4',
#   'CD4'='TRDC',
#   'CD4'='IL2RA',
#   'CD4'='CD25',
#   'CD4'='STAT1',
#   'CD4'='IL2',
#   'Tfh'='CD200',
#   'Tfh'='CD83',
#   'Tfh'='CXCL3',
#   'Tfh'='ICOS',
#   'Tfh'='BATF',
#   'Tfh'='PDCD1',
#   "B cells" = "CD19",
#   "B cells" = "CD79A",
#   "B cells" = "MS4A1",
#   'proliferation'='MKI67',
#   'proliferation'='H4C3',
#   'proliferation'='TUBB',
#   'proliferation'='STMN1',
#   'proliferation'='TOP1',
#   'proliferation'='ANAPC11',
#   'proliferation'='TOP2A'
# )

major_marker_ls_given <- c(
  # "B cells" = "CD19",
  # "B cells" = "CD79A",
  # "B cells" = "MS4A1",
  'T' ='CD3E',
  'CD4'='CD4',
  'CD4'='IL7R',
  'CD8'='CD8A',
  'CD8'='CD8B',
  # 'CD4'='TRDC',
  # 'CD4'='IL2RA',
  # 'CD4'='CD25',
  # 'CD4'='STAT1',
  # 'CD4'='IL2',
  'CD4_Treg'='TNFRSF4',
  'CD4_Treg'='TNFRSF9',
  'CD4_Treg'='FOXP3',
  'CD4_Treg'='IL2RA',
  'CD4_Treg'='LAYN',
  'CD4_Treg'='CTLA4',
  'CD4_Treg'='TIGIT',
  'CD4_naive'='LEF1',
  'CD4_Tem'='GZMK',
  'CD4_Tem'='GZMA',
  'CD4_Tem'='NKG7',
  'CD4_Tem'='CCL4',
  'CD4_Tem'='CCL5',
  'CD4_Tem'='CAPG',
  'CD4_Th1'='ANXA1',
  'CD4_Th1'='CXCL13',
  'CD4_Th1'='BHLHE40',
  # 'CD4_Th1'='IFNG',
  'CD4_Tfh'='TOX2',
  'CD4_Tfh'='PDCD1',
  'CD4_Tfh'='CXCR5',
  'CD4_Th17'='IL17A',
  'CD4_Th17'='IL17F',
  'CD4_Th17'='IL23R',
  'CD4_Th17'='RORC',
  'CD4_Th17'='CCR6',
  # 'CD4_Tcm'='TCF7',
  # 'CD4_Tcm'='CCR7',
  'CD4_Tcm'='S1PR1',
  'CD4_Tcm'='TXK',
  'CD4_Tcm'='ICAM2',
  'CD4_Tcm'='CDC25B',
  'CD4_Tcm'='PTGER2',
  'CD8_Tex'='HAVCR2',
  'CD8_Tex'='LAG3',
  # 'CD8_Tex'='ENPD1',
  'CD8_Tex'='TOX',
  'CD8_Tex'='CD279',
  'CD8_Tem'='CD44',
  'CD8_Tem'='GZMB',
  'CD8_Tem'='PRF1',
  'CD8_Tem'='KLRG1',
  'CD8_Tem'='TNF',
  'CD8_Tem'='IFNG',
  # 'CD8_Tnaive'='LEF1',
  'CD8_Tnaive'='CCR7',
  'CD8_Tnaive'='TCF1',
  'CD8_Tnaive'='TCF7',
  'CD8_Tnaive'='CD45RA',
  'CD8_MAIT'='NCR3',
  'CD8_MAIT'='SLC4A10',
  'CD8_MAIT'='KLRB1',
  'CD8_Trm'='MYADM',
  'CD8_Trm'='CD69',
  'CD8_Trm'='ITGAE',
  'CD8_Trm'='CXCR6',
  'CD8_Trm'='ITGA1',
  'CD8_Trm'='ZNF683',
  'CD8_Tcm'='CD28',
  'CD8_Tcm'='GPR183'
)


major_marker_ls_CD8 <- c(  
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8_Tex'='HAVCR2',
  'CD8_Tex'='CTLA4',
  'CD8_Tex'='LAG3',
  'CD8_Tex'='TIGIT',
  'CD8_Tex'='ENPD1',
  'CD8_Tem'='GZMK',
  'CD8_Tem'='GZMA',
  'CD8_Tem'='IFNG',
  'CD8_Tem'='KLRG1',
  'CD8_Tem'='GZMB',
  'CD8_Tem'='PRF1',
  'CD8_Tem'='NKG7',
  'CD8_Tnaive'='SELL',
  'CD8_Tnaive'='LEF1',
  'CD8_Tnaive'='TCF7',
  'CD8_Tnaive'='CCR7',
  'CD8_MAIT'='NCR3',
  'CD8_MAIT'='SLC4A10',
  'CD8_MAIT'='KLRB1',
  'CD8_Trm'='MYADM',
  'CD8_Trm'='CD69',
  'CD8_Trm'='ITGAE',
  'CD8_Trm'='CXCR6',
  'CD8_Trm'='ITGA1',
  'CD8_Trm'='ZNF683',
  'CD8_Tcm'='IL7R',
  'CD8_Tcm'='CD28',
  'CD8_Tcm'='GPR183'
)

major_marker_ls %>%
  str_to_title() %>%
  {names(.) <- names(major_marker_ls); .} -> major_marker_ls

major_marker_ls_given %>%
  str_to_title() %>%
  {names(.) <- names(major_marker_ls_given); .} -> major_marker_ls_given

major_marker_ls_CD8 %>%
  str_to_title() %>%
  {names(.) <- names(major_marker_ls_CD8); .} -> major_marker_ls_CD8

# 
# major_marker_cd8 <- c(  
#   'CD8'='Cd8a',
#   'CD8_Effector'='Ccl5',
#   'CD8_Effector'='Prf1',
#   'CD8_TEM'='Ccr7',
#   'CD8_TEM'='Sell',
#   'CD8_TEM'='Il7r',
#   'CD8_Exhausted'='Lag3',
#   'CD8_Exhausted'='Pdcd1',
#   'CD8_Exhausted'='Ctla4',
#   'CD8_Exhausted'='Havcr2',
#   'CD8_Exhausted'='Tox'
# )
# 
# major_marker_cd4 <- c(  
#   'CD4'='Cd4',
#   'CD4_Regulatory'='Tigit',
#   'CD4_Regulatory'='Foxp3',
#   'CD4_naive'='Tcf7',
#   'CD4_naive'='Ccr7'
# )




##映射------
library(tibble)
library(dplyr)
meta.data.t <-
  T_subtype1[['subtype']] %>%
  rownames_to_column('cells') %>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation3 = ifelse((as.character(annotation)=='T cells'),as.character(subtype), as.character(annotation)))%>%
  dplyr::select(cells, annotation3) %>%
  column_to_rownames('cells')


all_data <- AddMetaData(all_data, meta.data.t['annotation3'])

t_cell_types <- c("Epithelial_like", "CD4_Th17", "CD4_Tcm", "Epithelial_Mast", 
                  "Proliferating_T", "Endothelial_Pericytes", "Plasma_like", 
                  "CD4_Treg", "B_like", "NK", "Endothelial_like", "CD8_Tem")

all_data$annotation3 <- case_when(
  all_data$annotation3 %in% t_cell_types ~ "T cells", 
  TRUE ~ as.character(all_data$annotation3)
)

all_data$annotation3 <- ifelse(
  all_data$annotation3 == "Epithelial_Dc_like", 
  "DC", 
  as.character(all_data$annotation3)
)



pdf(file.path("UMAP with DC annotation3.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation3", label = TRUE, reduction = 'umap', pt.size = 0.5,repel = TRUE, cols = ggsci::pal_igv()(28)) 
#DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


cols_plan <- c(
  "Duct_cells_1" = "gray50",
  "Endothelial" = "gray50",
  "Macrophage" = "gray50",
  "Duct_cells_2" = "gray50",
  "Fibroblasts" = "gray50",
  "Neutrophil" = "gray50",
  "T cells" = "gray50",
  "DC" = "#600473"
)

pdf(file.path("UMAP with DC annotation3 highlight.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation3", label = TRUE, reduction = 'umap', pt.size = 0.5,repel = TRUE, cols = cols_plan) 
#DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()








# cell props barplot  --------
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

all_data$sample <- 
  factor(all_data$sample, 
         levels = c("CON1", "CON2", "CON3", "IRE1_3", "IRE2_3", "IRE3_3", "IRE1_7", "IRE3_7", "IRE4_7"))

p_col_proj <- 
  all_data_meta %>%
  group_by(sample, annotation3) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation3)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p_col_group <-
  all_data_meta %>%
  group_by(group, annotation3) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation3)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p_box_group  <-
  all_data_meta %>%
  group_by(group, sample, annotation3) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation3, y = freq, fill = group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) +
  ggsci::scale_fill_igv() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = "cell_project_new.pdf",
  height = 6,
  width  = 10
)

ggsave(
  p_col_group,
  filename = "cell_group_new.pdf",
  height = 6,
  width  = 6
)

ggsave(
  p_box_group,
  filename = "cell_freq_group_box_new.pdf",
  height = 6,
  width  = 10
)








#Fibroblast subtype--------------------
Fibroblast_subtype <- readRDS("Fibroblast_subtype/Fibroblast_subtype.rds")
Fibroblast_sub2 <- readRDS("Fibroblast_subtype/Fibroblast_sub2.rds")
Fibroblast <- readRDS("Fibroblast_subtype/second/Fibroblast.rds")

dim(Fibroblast_subtype)
dim(Fibroblast_sub2)
dim(Fibroblast)

head(Fibroblast_subtype)
head(Fibroblast_sub2)
head(Fibroblast)

colnames(Fibroblast_subtype@meta.data)
colnames(Fibroblast_sub2@meta.data)
colnames(Fibroblast@meta.data)


table(Fibroblast_subtype$subtype, useNA = "ifany")
table(Fibroblast_sub2$subtype, useNA = "ifany")
table(Fibroblast$subtype, useNA = "ifany")





##映射------
library(tibble)
library(dplyr)
meta.data.fib <-
  Fibroblast_subtype[['subtype']] %>%
  rownames_to_column('cells') %>%
  right_join(
    all_data[['annotation3']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation4 = ifelse((as.character(annotation3)=='Fibroblasts'),as.character(subtype), as.character(annotation3)))%>%
  dplyr::select(cells, annotation4) %>%
  column_to_rownames('cells')


all_data <- AddMetaData(all_data, meta.data.fib['annotation4'])

fib_cell_types <- c("mCAF", "apCAF", "myCAF", "Immune_like", "Macrophage_like", "dCAF", "iCAF")

all_data$annotation4 <- case_when(
  all_data$annotation4 %in% fib_cell_types ~ "Fibroblasts", 
  TRUE ~ as.character(all_data$annotation4)
)

all_data$annotation4 <- ifelse(
  all_data$annotation4 == "DCs_like", 
  "DC", 
  as.character(all_data$annotation4)
)

pdf(file.path("UMAP with DC annotation4.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation4", label = TRUE, reduction = 'umap', pt.size = 0.5,repel = TRUE, cols = ggsci::pal_igv()(28)) 
#DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


cols_plan <- c(
  "Duct_cells_1" = "gray50",
  "Endothelial" = "gray50",
  "Macrophage" = "gray50",
  "Duct_cells_2" = "gray50",
  "Fibroblasts" = "gray50",
  "Neutrophil" = "gray50",
  "T cells" = "gray50",
  "DC" = "#600473"
)

pdf(file.path("UMAP with DC annotation4 highlight.pdf"), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation4", label = TRUE, reduction = 'umap', pt.size = 0.5,repel = TRUE, cols = cols_plan) 
#DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




# cell props barplot  --------
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

p_col_proj <- 
  all_data_meta %>%
  group_by(sample, annotation4) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = annotation4)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p_col_group <-
  all_data_meta %>%
  group_by(group, annotation4) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = annotation4)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

p_box_group  <-
  all_data_meta %>%
  group_by(group, sample, annotation4) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation4, y = freq, fill = group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) +
  ggsci::scale_fill_igv() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = "cell_project_annotation4.pdf",
  height = 6,
  width  = 10
)

ggsave(
  p_col_group,
  filename = "cell_group_annotation4.pdf",
  height = 6,
  width  = 6
)

ggsave(
  p_box_group,
  filename = "cell_freq_group_box_annotation4.pdf",
  height = 6,
  width  = 10
)










# Macrophage subcluster ---------
Macrophage_subtype <- readRDS("~/Macrophage_subtype/Macrophage_subtype.rds")
Macrophage_2 <- readRDS("~/Macrophage_subtype/Second/Macrophage_2.rds")

head(Macrophage_subtype)
head(Macrophage_2)

dim(Macrophage_subtype)
dim(Macrophage_2)

table(Macrophage_subtype$subtype)
table(Macrophage_2$subtype)


pdf(file.path("Dotplot with DC annotation4.pdf"), height = 6, width = 15)
DotPlot(all_data, group.by = "annotation4", dot.scale = 10, features = major_marker_ls) + RotatedAxis()
# DimPlot(all_data, group.by = "annotation4", label = TRUE, reduction = 'umap', pt.size = 0.5,repel = TRUE, cols = ggsci::pal_igv()(28)) 
dev.off()





all_data_filt <- all_data

cells_to_keep <- !(all_data_filt$sample == "CON2" & all_data_filt$group == "IRE_4")

all_data_filt <- subset(all_data, cells = which(cells_to_keep))

dim(all_data)
dim(all_data_filt)

table(all_data_filt$sample, all_data_filt$group)

saveRDS(all_data_filt, "all_data_DC_annotation4_1202.rds")





all_data$annotation5 <- ifelse(
  all_data$annotation4 == "DC", 
  "DC", 
  as.character(all_data$annotation2)
)

table(all_data$annotation5, all_data$annotation4, useNA = "always")

saveRDS(all_data, "all_data_DC_annotation5_1208.rds")






