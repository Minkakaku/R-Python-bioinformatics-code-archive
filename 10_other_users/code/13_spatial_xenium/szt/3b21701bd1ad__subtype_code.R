library(Seurat)
library(tidyverse)
library(BPCells)
#All_T subtype-----------------------
All_T<-subset(all_data,annotation=='T cells')
### PCA --------
DefaultAssay(All_T) <- "RNA"
All_T <- NormalizeData(All_T, scale.factor = median(All_T$nCount_RNA))
All_T <- FindVariableFeatures(All_T)
All_T <- ScaleData(All_T, vars.to.regress = 'nCount_RNA')

All_T <- RunPCA(All_T)

out_dir<-'T_subtype/third'
dir.create(out_dir)
### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(All_T, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  All_T[[]] %>%
  rownames_to_column("cell_id2") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(All_T)) %>%
  dplyr::pull("cell_id2")


pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(All_T[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(All_T[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(All_T[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(All_T[, suset_cells], dims = 16:20)
dev.off()
### dimension reduction --------
dims <- 12


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


All_T <- RunUMAP(All_T, dims = 1:dims,return.model = TRUE)



pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "project", label = TRUE) + NoLegend()
dev.off()


All_T <- All_T %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

All_T <- RunUMAP(
  All_T, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
All_T  <-
  FindNeighbors(
    All_T ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


All_T <-
  FindClusters(
    All_T,
    resolution = c(0.2, 0.4, 0.6, 0.8, 1)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
# DimPlot(All_T, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(All_T, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(All_T, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(All_T, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(All_T, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(All_T, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

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
  
  
  # Enterocyte  = "SLC26A3",
  # Enterocyte  = "CA1",
  # Enterocyte  = "GUCA2A",
  # Enterocyte  = "GUCA2B",
  # 
  # Goblet = "MUC2",
  # Goblet = "TFF3",
  # 
  # Enteroendocrine = "CHGA",
  
  # "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  # "Fibroblasts" = 'COL1A2',
  # "Fibroblasts" = 'COL3A1',
  # "Fibroblasts" = 'LUM',
  # "Fibroblasts" = 'FN1',
  # "Fibroblasts" = 'COL6A1',
  # "Fibroblasts" = 'FBLN1',
  # 
  # "Fibroblasts" = 'VIM',
  # "Fibroblasts" = 'THY',
  # "Fibroblasts" = 'S100A4',
  # "Fibroblasts" = 'PDPN',
  # "Fibroblasts" = 'FAP',
  
  
  # 
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  
  
  # "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "PLN",
  "Smooth muscle cells" = "MYLK",
  "Smooth muscle cells" = "MYH11",
  # "Smooth muscle cells" = "ACTA2",
  
  
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  # "Plasma" = "XBP1",
  
  # Plasmablast = "MZB1",
  # Plasmablast = "XBP1",
  
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
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  
  'Neutrophi' = "S100A8",
  'Neutrophi' = "S100A9",
  'Neutrophi' = "CSF3R",
  'Neutrophi' = "FCGR3B"
  
  # 'Chondroblast' = 'ACAN',
  # 'Chondroblast' = 'COL2A1',
  # 'Chondroblast' = 'SOX9',
  # 'Chondroblast' = 'FMOD',
  # 'Chondroblast' = 'COMP'
)



major_marker_ls %>%
  str_to_title() %>%
  {names(.) <- names(major_marker_ls); .} -> major_marker_ls





major_marker_ls <- c(
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8'='GZMB',
  'CD8'='GZMA',
  'CD8'='GZMK',
  'CD8'='GZMH',
  'CD8'='CCL5',
  'Naive'='CCR6',
  'Naive'='LEF1',
  'Naive'='TCF7',
  'Naive'='SELL',
  'Naive'='CCR7',
  'Naive'='YBX3',
  'Naive'='TGFB1',
  'Naive'='CCl5',
  'Naive'='CXCR3',
  'Memory'='GPR183',
  'Memory'='IL7R',
  'Memory'='FOS',
  'Memory'='JUN',
  'Memory'='KLRB1',
  'TRM'='CD89',
  'TRM'='CXCR8',
  # 'TEM'='IFNG',
  'TEM'='OCL4',
  'TEM'='OCL5',
  'TEM'='JAK3',
  # 'TCM'='IDO1',
  # 'TCM'='SKIV2L',
  # 'TCM'='TRAF2',
  # 'TCM'='SLC2A1',
  # 'TCM'='CD40LG',
  'TCM'='DUSP2',
  'TCM'='LYAR',
  'TCM'='ZNF683',
  'TCM'='IL32',
  'Effector'='FGFBP2',
  'Effector'='TYROBP',
  'Effector'='KLRC2',
  'Effector'='CX3CR1',
  'Effector'='FCGR3A',
  'Effector'='GNLY',
  'cytotoxic'='XCL1',
  'cytotoxic'='IFNG',
  'cytotoxic'='CCL4',
  # 'cytotoxic'='CCL5',
  'Exhausted'='CXCL13',
  'Exhausted'='LAYN',
  'Exhausted'='LAG3',
  'Exhausted'='HAVCR2',
  'Exhausted'='TNFRSF9',
  'Exhausted'='CTLA4',
  'Exhausted'='TIGIT',
  'Exhausted'='TOX',
  'Exhausted'='PDCD1',
  'CD4'='FOXP3',
  'CD4'='CD4',
  'CD4'='TRDC',
  'CD4'='IL2RA',
  'CD4'='CD25',
  'CD4'='STAT1',
  'CD4'='IL2',
  'Tfh'='CD200',
  'Tfh'='CD83',
  'Tfh'='CXCL3',
  'Tfh'='ICOS',
  'Tfh'='BATF',
  'Tfh'='PDCD1',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",


  'proliferation'='MKI67',
  'proliferation'='H4C3',
  'proliferation'='TUBB',
  'proliferation'='STMN1',
  'proliferation'='TOP1',
  'proliferation'='ANAPC11',
  # 'proliferation'='MCM7',
  # 'proliferation'='DNMT1',
  'proliferation'='TOP2A'
)

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
  'CD4_Th1'='IFNG',
  'CD4_Tfh'='TOX2',
  'CD4_Tfh'='PDCD1',
  'CD4_Tfh'='CXCR5',
  'CD4_Th17'='IL17A',
  'CD4_Th17'='IL17F',
  'CD4_Th17'='IL23R',
  'CD4_Th17'='RORC',
  'CD4_Th17'='CCR6',
  'CD4_Tcm'='TCF7',
  'CD4_Tcm'='CCR7',
  'CD4_Tcm'='S1PR1',
  'CD4_Tcm'='TXK',
  'CD4_Tcm'='ICAM2',
  'CD4_Tcm'='CDC25B',
  'CD4_Tcm'='PTGER2',
  'CD8_Tex'='HAVCR2',
  'CD8_Tex'='LAG3',
  'CD8_Tex'='CD279',
  'CD8_Tem'='GZMB',
  'CD8_Tem'='PRF1',
  'CD8_Tem'='KLRG1',
  'CD8_Tnaive'='SELL',
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

major_marker_ls_CD4 %>%
  str_to_title() %>%
  {names(.) <- names(major_marker_ls_CD4); .} -> major_marker_ls_CD4

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
# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 30, height = 6)
# DotPlot(All_T, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = major_marker_ls_CD4)+RotatedAxis()
# DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = major_marker_ls_CD8)+RotatedAxis()
# DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = major_marker_cd4)+RotatedAxis()
# DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = major_marker_cd8)+RotatedAxis()
dev.off()
# find markers -------
Idents(All_T) <- "RNA_snn_res.0.6"
all_markers <-
  FindAllMarkers(
    JoinLayers(All_T),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  # arrange(desc(pct.1)) %>%
  slice_head(n = 30) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

### annotation ----

cell_cluster <- c(
  '0'  = "CD4_Tcm",
  '1'  = "CD4_naive",
  '2'  = "CD4_Tcm",
  '3'  = "CD4_Tem",
  '4'  = "CD8_naive",
  '5'  = "CD8_Tex",
  '6'  = "CD4_Treg",
  '7'  = "CD4_Th17",
  '8'  = "CD4_Th1"
)


Idents(All_T) <- "RNA_snn_res.0.4"

All_T <- RenameIdents(All_T, cell_cluster)
All_T$subtype <- Idents(All_T)


pdf(file.path(out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 8, width = 10)
# DimPlot(All_T, group.by = "subtype", label = TRUE) + NoLegend()
DimPlot(All_T, group.by = "subtype", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf(file.path(out_dir,"dotplot celltype main_marker.pdf"), width = 30, height = 6)
DotPlot(All_T, group.by = "subtype", features = major_marker_ls_CD4, dot.scale = 10) + theme_90angle
# DotPlot(All_T, group.by = "subtype", features = major_marker_ls_CD8, dot.scale = 10) + theme_90angle
DotPlot(All_T, group.by = "subtype", features = major_marker, dot.scale = 10) + theme_90angle
dev.off()


### project data ------
All_T <- 
  All_T %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'harmonyumap
    dims = 1:dims, 
    
    refdata = list(subtype_full = "subtype")
  )


DefaultAssay(All_T) <- 'RNA'

pdf(file.path("/home/rstudio/dat_analysis_cycle3/subtype/All_T",str_glue("Dimplot dim{dims} annotation - projected.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "subtype_full", label = TRUE, reduction = 'full.harmonyumap.1') + NoLegend()
dev.off()

# saveRDS(All_T,"All_T_subtype.rds")
# all_data <- readRDS("CD8_subtype.rds")

#spatial plot------------
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")
All_T_meta <-
  All_T[[]] %>%
  rownames_to_column("cell_id")
cells_boundary_ls <- readRDS(file.path("R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls <- readRDS(file.path("R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", 
          "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", 
          "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879", 
          "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", 
          "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", "#ffc300", 
          "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){  #browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat2 <-
      cells_boundary %>%
      right_join(All_T_meta, by = "cell_id")
    
    p <- dat[which(dat$annotation != 'T cells'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray80")  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal()+
      theme_void() +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat2, aes(x = vertex_x, y = vertex_y, fill = subtype, group = cell_id)) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = plan)+
      theme_void()
    p
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("T_subtype/third/spatial_plot", str_glue("T subtype spatial plot {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(All_T_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })

all_cells <- unique(All_T_meta$subtype)
cell_cols <- structure(
  ggsci::pal_igv()(length(all_cells)),
  names = all_cells
)

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){  #browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray90")  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      theme_void()
    
    p_ls <-
      all_cells %>%
      set_names(., .) %>%
      map(function(cell){ #browser()
        dat_sub <- dat %>% dplyr::filter(annotation3 == .env[['cell']])
        p +
          geom_polygon(
            data = dat_sub,
            aes(group = cell_id), 
            fill = cell_cols[cell]
          )  +
          ggtitle(cell)
      })
    p_ls
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path("T_subtype/spatial_plot", str_glue("spatial plot {grp} celltype split.pdf")),
      height = 12,
      width  = 8
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })
# cell props barplot  --------

All_T_meta <-
  All_T[[]] %>%
  rownames_to_column("cell_id")

All_T_meta$project <- factor(All_T_meta$project,levels=c("GX_1","GYF_1","HXF_1","LS_1","MXH_1","GX_2","GYF_2","HXF_2","LS_2","MXH_2","GX_3","GYF_3","HXF_3","LS_3","MXH_3"))
p_col_proj <-
  All_T_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


p_box_group  <-
  all_data_meta %>%
  group_by(group, project, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation, y = freq, fill = group)) +
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
  filename = file.path("T_subtype", "cell_freq_group.pdf"),
  height = 6,
  width  = 6
)

ggsave(
  p_box_group,
  filename = file.path(dr_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 10
)

##映射------
library(tibble)
library(dplyr)
meta.data.t <-
  All_T[['subtype']] %>%
  rownames_to_column('cells')%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation)=='T cells'),as.character(subtype), as.character(annotation)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation2'])
pdf(file.path("mapping_annotation T umap.pdf"), height = 6, width = 12)
DimPlot(all_data, group.by = "annotation3", label = TRUE, reduction = 'harmonyumap',pt.size = 2,repel = TRUE) 
#DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


#macrophage subtype--------------------
Macrophage_out_dir <- "Macrophage_subtype/Second"
dir.create(Macrophage_out_dir)
Macrophage <- subset(all_data, annotation=='Macrophage')
### PCA --------
Macrophage <- NormalizeData(Macrophage, scale.factor = median(Macrophage$nCount_RNA))
Macrophage <- FindVariableFeatures(Macrophage)
Macrophage <- ScaleData(Macrophage, vars.to.regress = 'nCount_RNA')

Macrophage <- RunPCA(Macrophage)


### pca inspector ------
pdf(file.path(Macrophage_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  Macrophage[[]] %>%
  rownames_to_column("cell_id2") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(Macrophage)) %>%
  dplyr::pull("cell_id2")


pdf(file.path(Macrophage_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage[, suset_cells], dims = 16:20)
dev.off()
### dimension reduction --------
dims <- 20 # 18


library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)


Macrophage <- RunUMAP(Macrophage, dims = 1:dims,return.model = TRUE)



pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "project", label = TRUE) + NoLegend()
dev.off()


Macrophage <- Macrophage %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Macrophage <- RunUMAP(
  Macrophage, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Macrophage  <-
  FindNeighbors(
    Macrophage ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Macrophage <-
  FindClusters(
    Macrophage,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
# DimPlot(Macrophage, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(Macrophage, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### dotplot -------


major_marker <- c(
  M0 = "Csf1r",
  M0 = "Cd11b",
  M0 = "Cd68",
  M0 = "Cd14",
  M1 = 'Macro',
  M1 = 'Cxcl9',
  M1 = 'Cxcl10',
  M1 = 'Cxcl11',
  M1 = 'Nos2',
  M1 = 'Socs1',
  M2 = 'Cd206',
  M2 = 'Tgm2',
  M2 = 'Fizz1',
  M2 = 'Chil3',
  M2 = 'Arg1',
  M2 = 'Ccl22',
  M2 = 'Cd163'
)


major_marker <- c(
  M0 = "Csf1r",
  M0 = "Cd11b",
  M0 = "Cd68",
  M0 = "Cd14",
  M1 = 'Cd80',
  M1 = 'Cd86',
  M1 = 'Tnf',
  M1 = 'Nos2',
  M1 = 'Cxcl9',
  M1 = 'Cxcl10',
  
  M2 = 'Cd163',
  M2 = 'Cd206',
  M2 = 'Mrc1',
  M2 = 'Arg1',
  M2 = 'Chil3',
  M2 = 'Il10',
  M2 = 'Tgfb1'
)


major_marker %>%
  str_to_title() %>%
  {names(.) <- names(major_marker); .} -> major_marker


# dotplot
pdf(file.path(Macrophage_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 30, height = 6)
# DotPlot(B, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
# DotPlot(B, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.2", features = major_marker)+RotatedAxis()
# DotPlot(B, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
dev.off()



Idents(Macrophage)<-'RNA_snn_res.0.2'


cell_cluster <- c(
  '0'  = "M2",
  '1'  = "M0",
  '2'  = "M2",
  '3'  = "M1",
  '4'  = "M1",
  '5'  = "M2",
  '6'  = "M2"
)


Macrophage <- RenameIdents(Macrophage, cell_cluster)
Macrophage$subtype <- Idents(Macrophage)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf(file.path(Macrophage_out_dir, "dotplot celltype main_marker.pdf"), width = 30, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(Macrophage, group.by = "subtype", features = major_marker, dot.scale = 10) + theme_90angle
DotPlot(Macrophage, group.by = "subtype", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

#spatial plot------------
Macrophage_meta <-
  Macrophage[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){  #browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat2 <-
      cells_boundary %>%
      right_join(Macrophage_meta, by = "cell_id")
    
    p <- dat[which(dat$annotation2 != 'Macrophage'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray80")  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal()+
      theme_void() +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat2, aes(x = vertex_x, y = vertex_y, fill = subtype, group = cell_id)) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = plan)+
      theme_void()
    p
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("Macrophage_subtype/Second/spatial_plot", str_glue("Macrophage subtype spatial plot {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })
# 
# sp_p_ls <-
#   cells_boundary_ls %>%
#   imap(function(cells_boundary, project){ # browser()
#     image_wh_r <- image_wh_r_ls[[project]]
#     
#     dat <-
#       cells_boundary %>%
#       right_join(Macrophage_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = vertex_x, y = vertex_y)) +
#       geom_polygon(aes(fill = subtype, group = cell_id))  +
#       scale_y_reverse() +
#       coord_fixed(ratio = 1/image_wh_r) +
#       ggsci::scale_fill_igv() +
#       theme_void()
#     p
#   })

all_cells <- unique(Macrophage_meta$subtype)
cell_cols <- structure(
  ggsci::pal_igv()(length(all_cells)),
  names = all_cells
)

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){  #browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray60")  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      theme_void()
    
    p_ls <-
      all_cells %>%
      set_names(., .) %>%
      map(function(cell){ #browser()
        dat_sub <- dat %>% dplyr::filter(annotation3 == .env[['cell']])
        p +
          geom_polygon(
            data = dat_sub,
            aes(group = cell_id), 
            fill = cell_cols[cell]
          )  +
          ggtitle(cell)
      })
    p_ls
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path("Macrophage_subtype/spatial_plot", str_glue("spatial plot {grp} celltype split.pdf")),
      height = 12,
      width  = 8
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })
# cell props barplot  --------

Macrophage_meta <-
  Macrophage[[]] %>%
  rownames_to_column("cell_id")

Macrophage_meta$project <- factor(Macrophage_meta$project,levels=c("GX_1","GYF_1","HXF_1","LS_1","MXH_1","GX_2","GYF_2","HXF_2","LS_2","MXH_2","GX_3","GYF_3","HXF_3","LS_3","MXH_3"))
p_col_proj <-
  Macrophage_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1')


temp %>% 
  rownames_to_column('cell_id') %>% 
  left_join(Duct_cells_1[[]] %>%
            select(group) %>% 
            rownames_to_column('cell_id'),by = 'cell_id') %>%
  column_to_rownames('cell_id') %>%
  write.csv('vlnplot_data.csv')

p_box_group  <-
  all_data_meta %>%
  group_by(group, project, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation, y = freq, fill = group)) +
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
  filename = file.path("Macrophage_subtype", "cell_freq_group.pdf"),
  height = 6,
  width  = 6
)

### project data ------
B <- 
  B %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'harmonyumap
    dims = 1:dims, 
    
    refdata = list(subtype_full = "subtype")
  )


DefaultAssay(B) <- 'RNA'

pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} annotation - projected.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "subtype_full", label = TRUE, reduction = 'full.harmonyumap.1') + NoLegend()
dev.off()

# saveRDS(B,"B_subtype.rds")

##映射------
library(tibble)
library(dplyr)
meta.data.t <-
  Macrophage[['subtype']] %>%
  rownames_to_column('cells')%>%
  right_join(
    all_data[['annotation2']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation2)=='Macrophage'),as.character(subtype), as.character(annotation2)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation2'])
saveRDS(Macrophage,'Macrophage_subtype/Macrophage.rds')


#Fibroblast subtype--------------------
out_dir<-'Fibroblast_subtype/second'
dir.create(out_dir)
Fibroblast<-subset(all_data,annotation2=='Fibroblast')
### PCA --------
DefaultAssay(Fibroblast) <- "RNA"
Fibroblast <- FindVariableFeatures(Fibroblast)
Fibroblast <- ScaleData(Fibroblast, vars.to.regress = 'nCount_RNA')

Fibroblast <- RunPCA(Fibroblast)


### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Fibroblast, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  Fibroblast[[]] %>%
  rownames_to_column("cell_id2") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(Fibroblast)) %>%
  dplyr::pull("cell_id2")


pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Fibroblast[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Fibroblast[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Fibroblast[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Fibroblast[, suset_cells], dims = 16:20)
dev.off()
### dimension reduction --------
dims <- 15


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


Fibroblast <- RunUMAP(Fibroblast, dims = 1:dims,return.model = TRUE)



pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "project", label = TRUE) + NoLegend()
dev.off()


Fibroblast <- Fibroblast %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Fibroblast <- RunUMAP(
  Fibroblast, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Fibroblast  <-
  FindNeighbors(
    Fibroblast ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Fibroblast <-
  FindClusters(
    Fibroblast,
    resolution = c(0.2, 0.4, 0.6, 0.8, 0.5)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
# DimPlot(Fibroblast, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(Fibroblast, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(Fibroblast, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(Fibroblast, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(Fibroblast, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Fibroblast, group.by = "RNA_snn_res.0.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

### dotplot -------
major_marker_ls <- c(
  # Epithelial  = "EPCAM", 
  # # Epithelial  = "KRT5",
  # Epithelial  = "KRT8",
  # Epithelial  = "S100A8",
  # Epithelial  = "KRT18",
  # 
  # "Fibroblasts" = "Pdgfra",
  # "Fibroblasts" = "Col1a2", 
  # 
 'mCAF'='MMP11',
 'mCAF'='COL1A1',
 'mCAF'='POSTN',
 'mCAF'='COMP',
 'mCAF'='COL10A1',
 'mCAF'='COL11A1',
 'mCAF'='DCN',
 'mCAF'='CDH11',
 'mCAF'='LUM',
 'iCAF'='PLA2G2A',
 'iCAF'='CFD',
 'iCAF'='CD34',
 'iCAF'='CXCL14',
 'iCAF'='PI16',
 'iCAF'='PRG4',
 'iCAF'='APOD',
 'iCAF'='DPT',
 # 'iCAF'='CXCL12',
 'iCAF'='CXCL1',
 'iCAF'='IL6',
 'iCAF'='CCL2',
 'apCAF'='HLA-DRA',
 'apCAF'='HLA-DRB1',
 'apCAF'='CD74',
 'apCAF'='IGLC2',
 'apCAF'='IGKC',
 'apCAF'='CXCR4',
 'vCAF'='ACTA2',
 'vCAF'='MCAM',
 'vCAF'='RGS5',
 'vCAF'='MYH11',
 'vCAF'='ADIRF',
 'tCAF'='ENO1',
 'tCAF'='GAPDH',
 'tCAF'='NDRG1',
 # 'tCAF'='HSPH1',
 # 'tCAF'='VEGFA',
 'tCAF'='MME',
 'dCAF'='MKI67',
 'dCAF'='STMN1',
 'dCAF'='PTTG1',
 'dCAF'='H2AFZ',
 'dCAF'='CENPF',
 'infCAF'='IDO1',
 'infCAF'='CXCL9',
 'infCAF'='CXCL10',
 'infCAF'='CXCL11',
 'rCAF'='CCL19',
 'rCAF'='CCL21',
 'pericyte'='COL4A1',
 'pericyte'='COL4A2',
 'hsptcaf'='HSPH1',
 'hsptcaf'='VEGFA'
 
)
major_marker <- c(
  'myCAF'='MMP11',
  'myCAF'='MYL9',
  'myCAF'='POSTN',
  'myCAF'='TPM2',
  'myCAF'='TAGLN',
  'myCAF'='HOPX',
  'mCAF'='DCN',
  'mCAF'='CDH11',
  'mCAF'='LUM',
  'iCAF'='PLA2G2A',
  'iCAF'='CFD',
  'iCAF'='CD34',
  'iCAF'='CXCL14',
  'iCAF'='PI16',
  'iCAF'='PRG4',
  'iCAF'='APOD',
  'iCAF'='DPT',
  'iCAF'='CXCL12',
  'iCAF'='CXCL1',
  'iCAF'='IL6',
  'iCAF'='CCL2',
  'apCAF'='HLA-DRA',
  'apCAF'='HLA-DRB1',
  'apCAF'='CD74',
  'apCAF'='IGLC2',
  'apCAF'='IGKC',
  'apCAF'='CXCR4',
  'vCAF'='ACTA2',
  'vCAF'='MCAM',
  'vCAF'='RGS5',
  'vCAF'='MYH11',
  'vCAF'='ADIRF',
  'tCAF'='ENO1',
  'tCAF'='GAPDH',
  'tCAF'='NDRG1',
  # 'tCAF'='HSPH1',
  # 'tCAF'='VEGFA',
  'tCAF'='MME',
  'dCAF'='MKI67',
  'infCAF'='IDO1',
  'infCAF'='CXCL9',
  'infCAF'='CXCL10',
  'infCAF'='CXCL11',
  'rCAF'='CCL19',
  'rCAF'='CCL21'
  # 'hsptcaf'='HSPH1',
  # 'hsptcaf'='VEGFA'
)

major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  Endothelial = "CLDN5",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  Epithelial  = "CEACAM6",
  Epithelial  = "CDH1",

  
  # Enterocyte  = "SLC26A3",
  # Enterocyte  = "CA1",
  # Enterocyte  = "GUCA2A",
  # Enterocyte  = "GUCA2B",
  # 
  # Goblet = "MUC2",
  # Goblet = "TFF3",
  # 
  # Enteroendocrine = "CHGA",
  
  # "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  # "Fibroblasts" = 'COL1A2',
  # "Fibroblasts" = 'COL3A1',
  # "Fibroblasts" = 'LUM',
  # "Fibroblasts" = 'FN1',
  # "Fibroblasts" = 'COL6A1',
  # "Fibroblasts" = 'FBLN1',
  # 
  # "Fibroblasts" = 'VIM',
  # "Fibroblasts" = 'THY',
  # "Fibroblasts" = 'S100A4',
  # "Fibroblasts" = 'PDPN',
  # "Fibroblasts" = 'FAP',
  
  
  # 
  "Pericytes" = "PDGFRB",
  "Pericytes" = "RGS5",
  "Pericytes" = "MCAM",
  
  
  # "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "PLN",
  "Smooth muscle cells" = "MYLK",
  "Smooth muscle cells" = "MYH11",
  # "Smooth muscle cells" = "ACTA2",
  
  
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  # "Plasma" = "XBP1",
  
  # Plasmablast = "MZB1",
  # Plasmablast = "XBP1",
  
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "CD3D",
  "T cells" = "CD4",
  "T cells" = "CD8A",
  "T cells" = "TRAC",
  
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
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  
  'Neutrophi' = "S100A8",
  'Neutrophi' = "S100A9",
  'Neutrophi' = "CSF3R",
  'Neutrophi' = "FCGR3B"
  
  # 'Chondroblast' = 'ACAN',
  # 'Chondroblast' = 'COL2A1',
  # 'Chondroblast' = 'SOX9',
  # 'Chondroblast' = 'FMOD',
  # 'Chondroblast' = 'COMP'
)



major_marker %>%
  str_to_title() %>%
  {names(.) <- names(major_marker); .} -> major_marker


major_marker <- c(
  'myCAF'='HOPX',
  'myCAF'='POSTN',
  'myCAF'='TAGLN',
  'apCAF'='IGKC',
  'apCAF'='CXCR4',
  'iCAF'='PDGFRA',
  'iCAF'='CD34',
  'iCAF'='CXCL14',
  'iCAF'='CXCL12',
  'iCAF'='CXCL1',
  'iCAF'='CCL2'
 
)

major_marker %>%
  str_to_title() %>%
  {names(.) <- names(major_marker); .} -> major_marker

# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 20, height = 6)
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = major_marker_ls)+RotatedAxis()
DotPlot(Fibroblast, group.by = "RNA_snn_res.0.6", features = major_marker)+RotatedAxis()
dev.off()

# find markers -------
Idents(Fibroblast) <- "RNA_snn_res.0.4"
all_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )


all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  # arrange(desc(pct.1)) %>%
  slice_head(n = 30) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")



### annotation ----

cell_cluster <- c(
  '0'  = "apCAF",
  '1'  = "myCAF",
  '2'  = "iCAF",
  '3'  = "apCAF",
  '4'  = "iCAF",
  '5'  = "iCAF",
  '6'  = "myCAF",
  '7'  = "iCAF",
  '8'  = "iCAF",
  '9'  = "myCAF",
  '10'  = "apCAF"
)


Idents(Fibroblast) <- "RNA_snn_res.0.6"

Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype <- Idents(Fibroblast)


pdf(file.path(out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf("Fibroblast_subtype/second/dotplot celltype main_marker.pdf", width = 28, height = 6)
DotPlot(Fibroblast, group.by = "subtype", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(Fibroblast, group.by = "subtype", features = major_marker, dot.scale = 10) + theme_90angle
dev.off()


### project data ------
Fibroblast <- 
  Fibroblast %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",#pca.full
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap' harmonyumap
    dims = 1:dims, 
    
    refdata = list(subtype_full = "subtype")
  )


DefaultAssay(Fibroblast) <- 'RNA'

pdf(str_glue("Dimplot dim{dims} annotation - projected.pdf"), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype_full", label = TRUE, reduction = 'full.harmonyumap.1') + NoLegend()
dev.off()
DefaultAssay(Fibroblast) <- 'sketch'
pdf("dotplot celltype main_marker.pdf", width = 10, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(Fibroblast, group.by = "subtype_full", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()
# saveRDS(Fibroblast,"Fibroblast_subtype.rds")


Fibroblast <- 
  Fibroblast %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'harmonyumap
    dims = 1:dims, 
    
    refdata = list(subtype_full = "subtype")
  )


DefaultAssay(Fibroblast) <- 'RNA'

pdf(str_glue("Dimplot dim{dims} annotation - projected.pdf"), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype_full", label = TRUE, reduction = 'full.harmonyumap.1') + NoLegend()
dev.off()

# saveRDS(Fibroblast,"Fibroblast_subtype.rds")

#spatial plot------------
Fibroblast_meta <-
  Fibroblast[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){  #browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat2 <-
      cells_boundary %>%
      right_join(Fibroblast_meta, by = "cell_id")
    
    p <- dat[which(dat$annotation2 != 'Fibroblast'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray80")  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal()+
      theme_void() +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat2, aes(x = vertex_x, y = vertex_y, fill = subtype, group = cell_id)) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = plan)+
      theme_void()
    p
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("Fibroblast_subtype/second/spatial_plot", str_glue("Fibroblast subtype spatial plot {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })
# 
# sp_p_ls <-
#   cells_boundary_ls %>%
#   imap(function(cells_boundary, project){ # browser()
#     image_wh_r <- image_wh_r_ls[[project]]
#     
#     dat <-
#       cells_boundary %>%
#       right_join(Macrophage_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = vertex_x, y = vertex_y)) +
#       geom_polygon(aes(fill = subtype, group = cell_id))  +
#       scale_y_reverse() +
#       coord_fixed(ratio = 1/image_wh_r) +
#       ggsci::scale_fill_igv() +
#       theme_void()
#     p
#   })

all_cells <- unique(Fibroblast_meta$subtype)
cell_cols <- structure(
  ggsci::pal_igv()(length(all_cells)),
  names = all_cells
)

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){  #browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray60")  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      theme_void()
    
    p_ls <-
      all_cells %>%
      set_names(., .) %>%
      map(function(cell){ #browser()
        dat_sub <- dat %>% dplyr::filter(annotation3 == .env[['cell']])
        p +
          geom_polygon(
            data = dat_sub,
            aes(group = cell_id), 
            fill = cell_cols[cell]
          )  +
          ggtitle(cell)
      })
    p_ls
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path("Fibroblast_subtype/spatial_plot", str_glue("spatial plot {grp} celltype split.pdf")),
      height = 12,
      width  = 8
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })

##映射------
library(tibble)
library(data.table)
meta.data.t <-
  Fibroblast[['subtype']] %>%
  rownames_to_column('cells')%>%
  right_join(
    all_data[['annotation2']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation2)=='Fibroblasts'),as.character(subtype), as.character(annotation2)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation2'])
# meta.data <-
#   rbindlist(list(CD4T[['subtype_full']]%>%rownames_to_column('cells'),
#                  CD8T[['subtype_full']]%>%rownames_to_column('cells'),
#                  B[['subtype_full']]%>%rownames_to_column('cells'), 
#                  Fibroblast[['subtype_full']]%>%rownames_to_column('cells')))%>%as.data.frame()%>%
#   right_join(
#     all_data_cycle3[['annotation3']] %>% rownames_to_column('cells'),
#     by = "cells"
#   ) %>% 
#   mutate(
#     annotation4 = ifelse((as.character(annotation3)=='CD4-T'|as.character(annotation3)=='CD8-T'|as.character(annotation3)=='B'|as.character(annotation3)=='Epi0_B'|as.character(annotation3)=='CD83-dim B'|as.character(annotation3)=='CD83 B'|as.character(annotation3)=='Fibroblast'|as.character(annotation3)=='Epi0_Fibroblast'),as.character(subtype_full), as.character(annotation3)))%>%
#   dplyr::select(cells, annotation4) %>%
#   column_to_rownames('cells')
# 
# 
# all_data_cycle3 <- AddMetaData(all_data_cycle3, meta.data['annotation4'])

# cell props barplot  --------
Fibroblast_meta <-
  Fibroblast[[]] %>%
  rownames_to_column("cell_id")

Fibroblast_meta$project <- factor(Fibroblast_meta$project,levels=c("GX_1","GYF_1","HXF_1","LS_1","MXH_1","GX_2","GYF_2","HXF_2","LS_2","MXH_2","GX_3","GYF_3","HXF_3","LS_3","MXH_3"))
p_col_proj <-
  Fibroblast_meta %>%
  group_by(group, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj,
  filename = file.path("Fibroblast_subtype", "cell_freq_group.pdf"),
  height = 6,
  width  = 6
)


saveRDS(Fibroblast,'Fibroblast_subtype/Fibroblast.rds')



#Findmarkers -----------------
Idents(Fibroblast)<-'sketch_snn_res.0.6'
Fibroblast_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 1,
    min.pct = 0.1,
    only.pos = TRUE
  )

# export marker gene
All_T_markers %>%
  dplyr::select(gene, everything()) %>%
  write.csv(
    file = file.path(dr_outdir, "marker.csv"), 
    row.names = FALSE
  )

CD8T_de_genes <-
  CD8T_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  pull("gene")

Fibroblast_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(abs(avg_log2FC) > 0.75) %>%
  # arrange(desc(avg_log2FC)) %>%
  arrange(desc(pct.1),.by_group = TRUE) %>%
  slice_head(n =20) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

# spatial plot -----------------
sp_out_dir <- file.path(out_dir, 'spatial plot') # dir.create(sp_out_dir)

cells_boundary_ls <- readRDS("/home/rstudio/data_for_analysis/cells_boundary_ls.rds")
image_wh_r_ls     <- readRDS('/home/rstudio/data_for_analysis/image_wh_r_ls.rds')

all_data_meta <- all_data_cycle3[[]] %>% rownames_to_column('cell_id')

all_data_T <- all_data_T %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "CD4-Tmemory" | 
      as.character(annotation4) == "CD4-Tnaive" | 
      as.character(annotation4) == "CD8-Tem" | 
      as.character(annotation4) == "CD8-Tex" | 
      as.character(annotation4) == "CD8-Tnaive" | 
      as.character(annotation4) == "CD8-Trm" | 
      as.character(annotation4) == "Tfh" | 
      as.character(annotation4) == "Treg" | 
      as.character(annotation4) == "proliferating-T" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))
sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_T[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_T, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 "#5DB1DDFF",'grey90',"#802268FF","#6BD76BFF","#D595A7FF"))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

# all_data_T <- all_data_cycle3[[]] %>% rownames_to_column("cell_id")
all_data_T <- all_data_T %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "ASC" | 
    as.character(annotation4) == "B-memory" | 
      as.character(annotation4) == "B-naive" | 
      as.character(annotation4) == "B-GC" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))
sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_T[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_T, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF","#802268FF", '#D595A7', 'grey90'))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_T <- all_data_T %>%
  mutate(subtype3 = case_when(
    as.character(annotation4) == "apCAF" | 
      as.character(annotation4) == "iCAF" | 
      as.character(annotation4) == "mCAF"| 
      as.character(annotation4) == "vCAF"| 
      as.character(annotation4) == "rCAF"| 
      as.character(annotation4) == "Epithelial"  ~ as.character(annotation4),
    TRUE ~ "Others"
  ))
sp_p_ls_Fibro <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_T[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_T, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype3, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF",'grey90',"#466983FF","#BA6338FF"))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


# export image
sp_p_ls_Fibro %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

###score--------------------
M1<-c( 'NOS2',
       'IL12',
       'FCGR1A',
       'FCGR1B',
       'FCGR1C',
       'CD80',
       'CXCR10',
       'IL23',
       'CXCL9',
       'CXCL10',
       'CXCL11',
       'CD86',
       'IL1A',
       'IL1B',
        'IL6',
        'TNF',
       'MHCII',
       'CCL5',
       'IRF5',
       'IRF1',
       'CD40',
       'KYNU',
       'CCR7',
        'VCAN',
        'CD300E',
        'INHBA',
        'TNIP3')
M2<-c( 'ARG1',
       'IL10',
       'CD32',
       'CD163',
       'FCER2',
       'CD200R1',
       'PDCD1LG2',
       'CD274',
       'MARCO',
       'CSF1R',
       'MRC1',
       'IL1RN',
       'Il1R2',
       'IL4R',
       'CCL4',
       'CCL13',
       'CCL20',
       'CCL17',
       'CCL18',
       'CCL22',
       'CCL24',
       'LYVE1',
       'VEGFA',
       'VEGFB',
       'VEGFC',
       'VEGFD',
       'EGF',
       'CTSA',
       'CTSB',
       'CSTC',
       'CTSD',
       'TGFB1',
       'CD209',
       'CD163L1',
       'FOLR2')

pathway<-list(M1,M2)
names(pathway)<-c('M1','M2')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  Macrophage <- AddModuleScore(JoinLayers(Macrophage), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=50)
}

Macrophage <- Macrophage[[]] %>% rownames_to_column('cell_id')
p<-VlnPlot(Macrophage,features = c('M1','M2'),group.by='RNA_snn_res.0.4',pt.size = 0)

ggsave(
  p,
  filename = 'Macrophage_subtype/M1_M2 score.pdf',
  height = 6,
  width  = 8
)
all_data_T <- all_data_T %>%
  mutate(subtype4 = case_when(
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))
sp_p_ls_tumor <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_T[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_T, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype4, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#CE3D32FF",'grey90'))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


# export image
sp_p_ls_tumor %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })
all_data_meta <- all_data_cycle3[[]] %>% rownames_to_column('cell_id')
#HLA--------
# sp_out_dir <- "dat_analysis_cycle3/subtype/newplot"

all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "CD4-Tmemory" | 
      as.character(annotation4) == "CD4-Tnaive" | 
      as.character(annotation4) == "CD8-Tem" | 
      as.character(annotation4) == "CD8-Tex" | 
      as.character(annotation4) == "CD8-Tnaive" | 
      as.character(annotation4) == "CD8-Trm" | 
      as.character(annotation4) == "Tfh" | 
      as.character(annotation4) == "Treg" | 
      as.character(annotation4) == "proliferating-T" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_HLA <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 'grey90',"#802268FF","#6BD76BFF","#D595A7FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = HLA, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_HLA %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot T {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "B-memory" | 
      as.character(annotation4) == "B-naive" | 
      as.character(annotation4) == "B-GC" | 
      as.character(annotation4) == "ASC" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_HLA <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF",'grey90',"#802268FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = HLA, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_HLA %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot B_plasma {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
      as.character(annotation4) == "apCAF" | 
      as.character(annotation4) == "iCAF" | 
      as.character(annotation4) == "mCAF"| 
      as.character(annotation4) == "vCAF"| 
      as.character(annotation4) == "rCAF"| 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_HLA <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF",'grey90',"#802268FF","#466983FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = HLA, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_HLA %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot fibroblast {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

#tumor_cell_LN_metastasis----------
all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "CD4-Tmemory" | 
      as.character(annotation4) == "CD4-Tnaive" | 
      as.character(annotation4) == "CD8-Tem" | 
      as.character(annotation4) == "CD8-Tex" | 
      as.character(annotation4) == "CD8-Tnaive" | 
      as.character(annotation4) == "CD8-Trm" | 
      as.character(annotation4) == "Tfh" | 
      as.character(annotation4) == "Treg" | 
      as.character(annotation4) == "proliferating-T" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_tumor_cell_LN_metastasis <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 'grey90',"#802268FF","#6BD76BFF","#D595A7FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = tumor_cell_LN_metastasis, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-1,0,1,2))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_tumor_cell_LN_metastasis %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot T {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "B-memory" | 
      as.character(annotation4) == "B-naive" | 
      as.character(annotation4) == "B-GC" | 
      as.character(annotation4) == "ASC" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_tumor_cell_LN_metastasis <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF",'grey90',"#802268FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = tumor_cell_LN_metastasis, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_tumor_cell_LN_metastasis %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot B_plasma {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "apCAF" | 
      as.character(annotation4) == "iCAF" | 
      as.character(annotation4) == "mCAF"| 
      as.character(annotation4) == "vCAF"| 
      as.character(annotation4) == "rCAF"| 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_tumor_cell_LN_metastasis <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF",'grey90',"#802268FF","#466983FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = tumor_cell_LN_metastasis, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_tumor_cell_LN_metastasis %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot fibroblast {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })






all_data_meta2$Inflammation[which(all_data_meta2$subtype4=='Others')]<-NA
sp_p_ls_Inflammation <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta2[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta2, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = Inflammation, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      # scale_fill_manual(values=c("#CE3D32FF",'grey90'))+
      scale_fill_gradient(low = "#FFB800", high ="#CC00FF",na.value = 'grey90')+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })

sp_p_ls_Inflammation %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


#IFN-------------
all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "CD4-Tmemory" | 
      as.character(annotation4) == "CD4-Tnaive" | 
      as.character(annotation4) == "CD8-Tem" | 
      as.character(annotation4) == "CD8-Tex" | 
      as.character(annotation4) == "CD8-Tnaive" | 
      as.character(annotation4) == "CD8-Trm" | 
      as.character(annotation4) == "Tfh" | 
      as.character(annotation4) == "Treg" | 
      as.character(annotation4) == "proliferating-T" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_IFN <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 'grey90',"#802268FF","#6BD76BFF","#D595A7FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = IFN, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-1,0,1,2))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_IFN %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot T {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "B-memory" | 
      as.character(annotation4) == "B-naive" | 
      as.character(annotation4) == "B-GC" | 
      as.character(annotation4) == "ASC" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_IFN <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF",'grey90',"#802268FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = IFN, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_IFN %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot B_plasma {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "apCAF" | 
      as.character(annotation4) == "iCAF" | 
      as.character(annotation4) == "mCAF"| 
      as.character(annotation4) == "vCAF"| 
      as.character(annotation4) == "rCAF"| 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_IFN <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF",'grey90',"#802268FF","#466983FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = IFN, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_IFN %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot fibroblast {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

#Inflammation----------
all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "CD4-Tmemory" | 
      as.character(annotation4) == "CD4-Tnaive" | 
      as.character(annotation4) == "CD8-Tem" | 
      as.character(annotation4) == "CD8-Tex" | 
      as.character(annotation4) == "CD8-Tnaive" | 
      as.character(annotation4) == "CD8-Trm" | 
      as.character(annotation4) == "Tfh" | 
      as.character(annotation4) == "Treg" | 
      as.character(annotation4) == "proliferating-T" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_Inflammation <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 'grey90',"#802268FF","#6BD76BFF","#D595A7FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = Inflammation, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-1,0,1,2))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_Inflammation %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot T {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "B-memory" | 
      as.character(annotation4) == "B-naive" | 
      as.character(annotation4) == "B-GC" | 
      as.character(annotation4) == "ASC" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_Inflammation <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF",'grey90',"#802268FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = Inflammation, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_Inflammation %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot B_plasma {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "apCAF" | 
      as.character(annotation4) == "iCAF" | 
      as.character(annotation4) == "mCAF"| 
      as.character(annotation4) == "vCAF"| 
      as.character(annotation4) == "rCAF"| 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_Inflammation <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF",'grey90',"#802268FF","#466983FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = Inflammation, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_Inflammation %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot fibroblast {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })
#cancer_stem_cell_and_EMT-----------
all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "CD4-Tmemory" | 
      as.character(annotation4) == "CD4-Tnaive" | 
      as.character(annotation4) == "CD8-Tem" | 
      as.character(annotation4) == "CD8-Tex" | 
      as.character(annotation4) == "CD8-Tnaive" | 
      as.character(annotation4) == "CD8-Trm" | 
      as.character(annotation4) == "Tfh" | 
      as.character(annotation4) == "Treg" | 
      as.character(annotation4) == "proliferating-T" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_cancer_stem_cell_and_EMT <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 'grey90',"#802268FF","#6BD76BFF","#D595A7FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = cancer_stem_cell_and_EMT, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4,6,8))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_cancer_stem_cell_and_EMT %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot T {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "B-memory" | 
      as.character(annotation4) == "B-naive" | 
      as.character(annotation4) == "B-GC" | 
      as.character(annotation4) == "ASC" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_cancer_stem_cell_and_EMT <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF",'grey90',"#802268FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = cancer_stem_cell_and_EMT, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4,6,8))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_cancer_stem_cell_and_EMT %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot B_plasma {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "apCAF" | 
      as.character(annotation4) == "iCAF" | 
      as.character(annotation4) == "mCAF"| 
      as.character(annotation4) == "vCAF"| 
      as.character(annotation4) == "rCAF"| 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_cancer_stem_cell_and_EMT <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF",'grey90',"#802268FF","#466983FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = cancer_stem_cell_and_EMT, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4,6,8))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_cancer_stem_cell_and_EMT %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot fibroblast {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


pdf(file.path(str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 10)
DimPlot(all_data, group.by = "annotation2", label = TRUE, reduction = "umap",repel = TRUE,raster=FALSE,cols = ggsci::pal_igv()(28))
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


# CN barplot
data <- read.csv("CNs_15_all_cell_features.xls",sep='\t')

library(ggplot2)
library(tidyr)

library(purrr)

p<-
  data %>%
  dplyr::select(cell_id,sample,group,CNs_15) %>%
  group_by(sample,group,CNs_15) %>%
  summarise(count = n()) %>%
  mutate(proportion = count/sum(count)) %>%
  
  ggplot(aes(x = CNs_15, y = proportion, fill = group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0),alpha = 0.2) +
  
  ggsci::scale_fill_npg() +
  theme_classic() +
  theme(
    axis.text.x = element_text(size= 10,angle = 90,face="bold", hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size= 10,face="bold"),
    axis.title.x = element_text(size=12,face="bold"),
    axis.title.y = element_text(size=12,face="bold")
  )
p <- p + ggpubr::geom_pwc(label = "p.format", method = "wilcox.test")


if (length(unique(data$group)) == 3) {
  # 如果只有两组，添加显著性标记
  p <- p + ggpubr::geom_pwc(label = "p.signif", method = "wilcox.test")
}

print(p)



data %>%
  dplyr::select(cell_id,sample,group,CNs_10) %>%
  group_by(sample,group,CNs_10) %>%
  summarise(count = n()) %>%
  mutate(proportion = count/sum(count)) %>%
  ggpubr::compare_means(
    data = .,
    group.by = "CNs_10",
    formula = proportion ~ group
  ) %>% print(n=30)
