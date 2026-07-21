# subtype
# export: dat_analysis_cycle3/subtype/all_data_cycle3_1205.rds


library(Seurat)
library(tidyverse)
# library(BPCells)
#All_T subtype-----------------------
all_data_cycle3 <- readRDS('/home/rstudio/dat_analysis_cycle3/all_data_cycle3.rds')
All_T <- subset(
  all_data_cycle3,
  annotation2   == 'CD4T'| annotation2=='CD4T-Macrophage'|
    annotation2 == 'CD4 Endothelial'| annotation2=='Epi0_CD4T'|
    annotation2 == 'Epi0_CD83_T'| annotation2=='Epi0_T'|
    annotation2 == 'CD8T'
)

### PCA --------
DefaultAssay(All_T) <- "sketch"
All_T <- FindVariableFeatures(All_T)
All_T <- ScaleData(All_T, vars.to.regress = 'nCount_RNA')

All_T <- RunPCA(All_T)


### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(All_T, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  All_T[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(All_T)) %>%
  dplyr::pull("cell_id")


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
dims <- 20


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
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(All_T, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(All_T, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(All_T, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(All_T, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

### dotplot -------
major_marker_ls <- c(
  'CD8'='CD8A',
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
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  'Endothelial' = "PECAM1",
  'Endothelial' = "VWF",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC",
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "JCHAIN",
  'NK'='GNLY',
  'NK'='KLRC1',
  'NK'='FGFBP2',
  'NK'='NKG7',
  'NK'='KLRF1',
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  # "Fibroblasts" = "SFRP4", 
  # "Fibroblasts" = "PDGFRA", 
  # "Fibroblasts" = "CCDC80",
  # "Fibroblasts" = "COL1A1",
  
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "MYH11", 
  "Smooth muscle cells" = "ACTA2",
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",

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
# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(All_T, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(All_T, group.by = "sketch_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(All_T, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()

### annotation ----

cell_cluster <- c(
  '0'  = "CD4_T",
  '1'  = "CD8_T",
  '2'  = "Endothelial",
  '3'  = "Macrophage",
  '4'  = "B",
  '5'  = "CD4_T",
  '6'  = "Epithelial",
  '7'  = "NK",
  '8'  = "CD4_T",
  '9'  = "DCs",
  '10'  = "Mast",
  '11'  = "proliferating_T",
  '12'  = "Smooth muscle cells"
)


Idents(All_T) <- "sketch_snn_res.0.4"

All_T <- RenameIdents(All_T, cell_cluster)
All_T$subtype <- Idents(All_T)


pdf(str_glue("Dimplot dim{dims} annotation.pdf"), height = 6, width = 8)
DimPlot(All_T, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf("dotplot celltype main_marker.pdf", width = 10, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(All_T, group.by = "subtype", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
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

##映射------
library(tibble)
library(dplyr)
meta.data.t <-
  All_T[['subtype_full']] %>%
  rownames_to_column('cells')%>%
  right_join(
    all_data_cycle3[['annotation2']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation3 = ifelse((as.character(annotation2)=='CD4T'|as.character(annotation2)=='CD4 Endothelial'|as.character(annotation2)=='CD4T-Macrophage'|as.character(annotation2)=='Epi0_CD4T'|as.character(annotation2)=='Epi0_CD83_T'|as.character(annotation2)=='Epi0_T'|as.character(annotation2)=='CD8T'),as.character(subtype_full), as.character(annotation2)))%>%
  dplyr::select(cells, annotation3) %>%
  column_to_rownames('cells')


all_data_cycle3 <- AddMetaData(all_data_cycle3,meta.data.t['annotation3'])
pdf(file.path("mapping_annotation T umap.pdf"), height = 6, width = 12)
DimPlot(all_data_cycle3, group.by = "annotation3", label = TRUE, reduction = 'harmonyumap',pt.size = 2,repel = TRUE) 
#DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()




#CD8T subtype--------------------
CD8T<-subset(all_data_cycle3,annotation3=='CD8-T')
### PCA --------
DefaultAssay(CD8T) <- "sketch"
CD8T <- FindVariableFeatures(CD8T)
CD8T <- ScaleData(CD8T, vars.to.regress = 'nCount_RNA')

CD8T <- RunPCA(CD8T)


### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(CD8T, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  CD8T[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(CD8T)) %>%
  dplyr::pull("cell_id")


pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(CD8T[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(CD8T[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(CD8T[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(CD8T[, suset_cells], dims = 16:20)
dev.off()
### dimension reduction --------
dims <- 20


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


CD8T <- RunUMAP(CD8T, dims = 1:dims,return.model = TRUE)



pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(CD8T, group.by = "project", label = TRUE) + NoLegend()
dev.off()


CD8T <- CD8T %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

CD8T <- RunUMAP(
  CD8T, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(CD8T, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
CD8T  <-
  FindNeighbors(
    CD8T ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


CD8T <-
  FindClusters(
    CD8T,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(CD8T, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(CD8T, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(CD8T, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(CD8T, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(CD8T, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(CD8T, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(CD8T, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(CD8T, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(CD8T, group.by = "sub.cluster", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()



#subcluster
Idents(CD8T) <- "sketch_snn_res.0.6"
CD8T <- FindSubCluster(CD8T, cluster = 3, graph.name = "sketch_snn", algorithm = 1,resolution = 0.8)
Idents(CD8T) <- "sub.cluster"
### dotplot -------
major_marker_ls <- c(
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8'='TRAC',
  'Memory'='IL7R',
  'Memory'='FOS',
  'Memory'='JUN',
  'Memory'='KLRB1',
  'cytotoxic'='XCL1',
  # 'cytotoxic'='IFNG',
  # 'cytotoxic'='CCL4',
  'cytotoxic'='CCL5',
  'Exhausted'='CXCL13',
  'cytotoxic'='GZMB',
  # 'cytotoxic'='GZMK',
  # 'cytotoxic'='NKG7',
  # 'cytotoxic'='CST7',
  # 'cytotoxic'='PRF1',
  'Exhausted'='LAG3',
  'Exhausted'='TNFRSF9',
  'Exhausted'='CTLA4',
  'proliferation'='MKI67',
  'proliferation'='MCM7',
  'proliferation'='DNMT1',
  'proliferation'='TOP2A'
)
major_marker_ls <- c(
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8'='TRAC',
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
  'TEM'='GZMK',
  'TEM'='GZMA',
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
  'cytotoxic'='CCL5',
  'Exhausted'='CXCL13',
  'Exhausted'='LAYN',
  'Exhausted'='LAG3',
  'Exhausted'='HAVCR2',
  'Exhausted'='TNFRSF9',
  'Exhausted'='CTLA4',
  'Exhausted'='TIGIT',
  'Exhausted'='TOX',
  'Exhausted'='PDCD1',
  'proliferation'='MKI67',
  'proliferation'='H4C3',
  'proliferation'='TUBB',
  'proliferation'='STMN1',
  'proliferation'='TOP1',
  'proliferation'='ANAPC11',
  # 'proliferation'='MCM7',
  # 'proliferation'='DNMT1',
  'proliferation'='TOP2A',
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC",
  'Fat cell' = 'ADIPOQ',
  'Fat cell' = 'UCP1',
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1"
)


cd8t_markers <- c(
  'CD4',
  'CD8A',
  
  "CD28",
  "CCL4",
  "CCL5",
  "CCL3",
  "CCL20",
  "STAT1",
  "IL26",
  "IL32",
  "HIF1A",
  "TGFB1",
  
  CD8_Tex = "PDCD1",
  CD8_Tex = "HAVCR2",
  CD8_Tex = "CTLA4",
  CD8_Tex = "CXCL13",
  CD8_Tex = "LYAN",
  CD8_Tex = "LAG3",
  CD8_Tex = "TIGIT",
  CD8_Tex = "ENTPD1",
  
  CD8_Tem = "GZMK",
  CD8_Tem = "GZMA",
  CD8_Tem = "IFNG",
  CD8_Tem = "KLRG1",
  CD8_Tem = "GZMB",
  CD8_Tem = "RPF1",
  CD8_Tem = "NKG7",
  
  CD8_Trm = "MYADM",
  CD8_Trm = "CD69",
  CD8_Trm = "ITGAE",
  CD8_Trm = "CXCR6",
  CD8_Trm = "ITGA1",
  
  CD8_MAIT= "NCR3",
  CD8_MAIT= "SLC4A10",
  
  CD8_TNaive = "SELL",
  CD8_TNaive = "LEF1",
  CD8_TNaive = "TCF7",
  CD8_TNaive = "CCR7",
  
  CD8_Cycling = "MKI67",
  CD8_Cycling = "TOP2A"
)
FeaturePlot(CD8T,features = c(
  'H4C3'),reduction = 'harmonyumap')


FeaturePlot(CD8T,features = c(
  'IL23R','MUC1','TERT','ADIPOQ','CD274','INS','IPMK','PRF1','GZMA','GZMB','GZMK',
  'LAG3','CD200','TOX2','TOX','PRDM1',
  'NELL2','REG4','CD248','MGAT4A',
  'MALT1',"IDO1","SKIV2L",'TRAF2',
  "GFL1",'JAK3','PDCD1LG2','SELPLG'),reduction = 'harmonyumap')
# dotplot
pdf(file.path(outdir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(CD8T, group.by = "sketch_snn_res.0.2", features = cd8t_markers)
DotPlot(CD8T, group.by = "sketch_snn_res.0.4", features = cd8t_markers)
DotPlot(CD8T, group.by = "sketch_snn_res.0.6", features = cd8t_markers)+RotatedAxis()
DotPlot(CD8T, group.by = "sub.cluster", features = cd8t_markers)+RotatedAxis()
DotPlot(CD8T, group.by = "sketch_snn_res.0.8", features = cd8t_markers)+RotatedAxis()
dev.off()




### annotation ----

cell_cluster <- c(
  '0'  = "CD8_Tnaive",
  '1'  = "CD8_Tnaive",
  '2'  = "CD8_Tem",
  '3_0'  = "CD8_Tnaive",
  '3_1'  = "CD8_Tnaive",
  '3_2'  = "CD8_Trm",
  '3_3'  = "CD8_Tnaive",
  '3_4'  = "CD8_Tem",
  '3_5'  = "CD8_Tnaive",
  '4'  = "CD8_Tem",
  '5'  = "CD8_Tem",
  '6'  = "CD8_Tem",
  '7'  = "CD8_Tem",
  '8'  = "CD8_Tem",
  '9'  = "CD8_Tex",
  '10'  = "CD8_Tem",
  '11'  = "CD8_Tem",
  '12'  = "CD8_Tem"
  
)


Idents(CD8T) <- "sub.cluster"

CD8T <- RenameIdents(CD8T, cell_cluster)
CD8T$subtype <- Idents(CD8T)


pdf(str_glue("Dimplot dim{dims} annotation.pdf"), height = 6, width = 8)
DimPlot(CD8T, group.by = "subtype", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf("dotplot celltype main_marker.pdf", width = 10, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(CD8T, group.by = "subtype", features = cd8t_markers, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()


major_marker_ls <- c(
  'Naive'='CCR6',
  'Naive'='LEF1',
  'Naive'='TCF7',
  'Naive'='SELL',
  'Memory'='GPR183',
  'Memory'='IL7R',
  'Memory'='FOS',
  'Memory'='JUN',
  'Memory'='KLRB1',
  'Effector'='FGFBP2',
  'Effector'='TYROBP',
  'Effector'='KLRC2',
  'Effector'='CX3CR1',
  'Effector'='FCGR3A',
  'cytotoxic'='XCL1',
  'cytotoxic'='IFNG',
  'cytotoxic'='CCL4',
  'cytotoxic'='CCL5',
  'Exhausted'='CXCL13',
  'Exhausted'='LAYN',
  'Exhausted'='LAG3',
  'Exhausted'='HAVCR2',
  'Exhausted'='TNFRSF9',
  'Exhausted'='CTLA4',
  'Exhausted'='TIGIT',
  'Exhausted'='TOX',
  'Exhausted'='PDCD1'
)

### project data ------
CD8T <- 
  CD8T %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'harmonyumap
    dims = 1:dims, 
    
    refdata = list(subtype_full = "subtype")
  )

DefaultAssay(CD8T) <- 'RNA'

pdf(str_glue("Dimplot dim{dims} annotation - projected.pdf"), height = 6, width = 8)
DimPlot(CD8T, group.by = "subtype_full", label = TRUE, reduction = 'full.harmonyumap.1') + NoLegend()
dev.off()


pdf("dotplot celltype main_marker.pdf", width = 10, height = 6)
DotPlot(CD8T, group.by = "subtype_full", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()
# saveRDS(CD8T,"CD8_subtype.rds")
all_data <- readRDS("CD8_subtype.rds")







#CD4T subtype--------------------
CD4T <- subset(all_data_cycle3, annotation3 == 'CD4-T')
### PCA --------
DefaultAssay(CD4T) <- "sketch"
CD4T <- FindVariableFeatures(CD4T)
CD4T <- ScaleData(CD4T, vars.to.regress = 'nCount_RNA')

CD4T <- RunPCA(CD4T)


### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(CD4T, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  CD4T[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(CD4T)) %>%
  dplyr::pull("cell_id")


pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(CD4T[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(CD4T[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(CD4T[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(CD4T[, suset_cells], dims = 16:20)
dev.off()
### dimension reduction --------
dims <- 15


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


CD4T <- RunUMAP(CD4T, dims = 1:dims,return.model = TRUE)



pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(CD4T, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


CD4T <- CD4T %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

CD4T <- RunUMAP(
  CD4T, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(CD4T, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
CD4T  <-
  FindNeighbors(
    CD4T ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


CD4T <-
  FindClusters(
    CD4T,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(CD4T, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(CD4T, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(CD4T, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(CD4T, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(CD4T, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(CD4T, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(CD4T, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(CD4T, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


Idents(CD4T) <- "sketch_snn_res.0.4"
CD4T <- FindSubCluster(CD4T, cluster = 3, graph.name = "sketch_snn", algorithm = 1,resolution = 0.3)
Idents(CD4T) <- "sub.cluster"
DimPlot(CD4T, group.by = "sub.cluster", label = TRUE, reduction = 'harmonyumap') + NoLegend()

### dotplot -------
CD4T<-readRDS('./CD4_subtype.rds')
major_marker_ls <- c(
  'CD4'='CD4',
  # 'Endothelial' = "PECAM1",
  # 'Endothelial' = "VWF",
  # "Macrophage" = "CD68",
  # "Macrophage" = "CD163",
  # "Macrophage" = "C1QA",
  # "Macrophage" = "C1QC",
  # "B cells" = "CD19",
  # "B cells" = "CD79A",
  # "B cells" = "MS4A1",
  # "Plasma" = "JCHAIN",
  # 'Dendritic cell' = "LILRA4",  
  # 
  # 'Dendritic cell' = "CLEC9A",
  # 'Dendritic cell' = "IDO1",
  # 
  # 'Dendritic cell' = "CD1C",
  # 'Dendritic cell' = "FCER1A",
  # "Fibroblasts" = "SFRP4", 
  # "Fibroblasts" = "PDGFRA", 
  # "Fibroblasts" = "CCDC80",
  # "Fibroblasts" = "COL1A1",
  # 
  # "Smooth muscle cells" = "PDGFRB", 
  # "Smooth muscle cells" = "MYH11", 
  # "Smooth muscle cells" = "ACTA2",
  # "Mast cell" = "CPA3",
  # "Mast cell" = "CTSG",
  # "Mast cell" = "TPSAB1",
  CD4_Tem = "GZMK",
  CD4_Tem = "GZMA",
  CD4_Tem = "NKKG7",
  CD4_Tem = "CCL4",
  CD4_Tem = "CCL5",
  CD4_Tem = "CAPG",
  # CD4_Tcm = "MALT1",
  # CD4_Tcm = "PDK1",
  # CD4_Tcm = "SKIV2L",
  # CD4_Tcm = "IL23A",
  # CD4_Tcm = "CD28",
  # CD4_Tcm = "CD45RO",
  # CD4_Tcm = "CD62L",
  # CD4_Tcm = "CD44",
  'TCM'='CD27',
  'TCM'='CD69',
  'TCM'='CXCR4',
  'Naive'='SELL',
  'Naive'='LEF1',
  'Naive'='TCF7',
  'Naive'='S1PR1',
  # 'Naive'='CCR6',
  # 'Naive'='CCR7',
  'Naive'='YBX3',
  'Naive'='TGFB1',
  'Naive'='CCl5',
  'Naive'='CXCR3',
  'Treg'='FOXP3',
  'Treg'='CTLA4',
  'Treg'='CD25',
  'Treg'='IKZF2',
  'Treg'='IL2RA',
  'Treg'='TIGHT',
  'gdT'='TRDC',
  'gdT'='TRGC1',
  'gdT'='TRGC2',
  'Tfh'='CD200',
  'Tfh'='CD83',
  'Tfh'='CXCL3',
  'Tfh'='ICOS',
  'Tfh'='BATF',
  'Tfh'='PDCD1',
  'Th1'='IFNG',
  'Th1'='TNF',
  'Th2'='TBX21',
  'Th2'='GATA3',
  'Th17'='IL10RA',
  'Th17'='IL17',
  'Th17'='IL17A',
  'Th17'='IL17AF',
  'Th17'='IL17F',
  'Th17'='IL21',
  'Memory'='GPR183',
  'Memory'='IL7R',
  'Memory'='FOS',
  'Memory'='JUN',
  'Memory'='KLRB1',
  'TEM'='GFI1',
  'TEM'='SERPINB6',
  'TEM'='LRRC32',
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
cd4t_markers <- c(
  'CD4',
  'CD8A',
  
  CD4_Treg = "TNFRSF4",
  CD4_Treg = "TNFRSF9",
  CD4_Treg = "FOXP3",
  CD4_Treg = "IL2RA",
  CD4_Treg = "LAYN",
  CD4_Treg = "CTLA4",
  CD4_Treg = "TIGIT",
  CD4_Treg = "RTKN2",
  
  CD4_Tem = "GZMK",
  CD4_Tem = "GZMA",
  CD4_Tem = "NKKG7",
  CD4_Tem = "CCL4",
  CD4_Tem = "CCL5",
  CD4_Tem = "CAPG",
  
  CD4_Th17 = "IL17A",
  CD4_Th17 = "IL17F",
  
  "CD4_Th1-like" = "ANXA1",
  "CD4_Th1-like" = "CXCL13",
  "CD4_Th1-like" = "BHLHE40",
  "CD4_Th1-like" = "IFNG",
  
  CD4_Tfh = "TOX2",
  CD4_Tfh = "PDCD1",
  
  CD4_Tcm = "TCF7",
  CD4_Tcm = "CCR7",
  CD4_Tcm = "S1PR1",
  CD4_Tcm = "TXK",
  
  CD4_Cycling = "MKI67",
  CD4_Cycling = "TOP2A"
)
# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(CD4T, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(CD4T, group.by = "sketch_snn_res.0.4", features = major_marker_ls)
DotPlot(CD4T, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(CD4T, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()

### annotation ----

cell_cluster <- c(
  '0'  = "CD4_Tnaive",
  '1'  = "Treg",
  '2'  = "Tfh",
  '3'  = "CD4_Tmemory",
  '4'  = "CD4_Tnaive",
  '5'  = "Treg",
  '6'  = "CD4_Tnaive",
  '7'  = "CD4_Tnaive"
)


Idents(CD4T) <- "sketch_snn_res.0.4"

CD4T <- RenameIdents(CD4T, cell_cluster)
CD4T$subtype <- Idents(CD4T)


pdf(str_glue("Dimplot dim{dims} annotation.pdf"), height = 6, width = 8)
DimPlot(CD4T, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

DotPlot(CD4T, group.by = "sub.cluster", features = major_marker_ls, dot.scale = 10) + theme_90angle
pdf("dotplot celltype main_marker.pdf", width = 10, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(CD4T, group.by = "subtype", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()


### project data ------
CD4T <- 
  CD4T %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'harmonyumap
    dims = 1:dims, 
    
    refdata = list(subtype_full = "subtype")
  )


DefaultAssay(CD4T) <- 'RNA'

pdf(file.path("/home/rstudio/dat_analysis_cycle3/subtype/CD4T",str_glue("Dimplot dim{dims} annotation - projected.pdf")), height = 6, width = 8)
DimPlot(CD4T, group.by = "subtype_full", label = TRUE, reduction = 'full.harmonyumap.1') + NoLegend()
dev.off()

# saveRDS(CD4T,"CD4_subtype.rds")
# all_data <- readRDS("CD8_subtype.rds")

#B subtype--------------------
b_out_dir <- "dat_analysis_cycle3/subtype/B"
B <- subset(all_data_cycle3, annotation3=='B'|annotation3=='CD83-dim B'|annotation3=='CD83 B'|annotation3=='CD83 B'|annotation3=='Epi0_B'|annotation3=='Plasma')
### PCA --------
DefaultAssay(B) <- "sketch"
B <- FindVariableFeatures(B)
B <- ScaleData(B, vars.to.regress = 'nCount_RNA')

B <- RunPCA(B)


### pca inspector ------
pdf(file.path(b_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(B, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  B[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(B)) %>%
  dplyr::pull("cell_id")


pdf(file.path(b_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(B[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(b_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(B[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(b_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(B[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(b_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(B[, suset_cells], dims = 16:20)
dev.off()
### dimension reduction --------
dims <- 20 # 18


library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)


B <- RunUMAP(B, dims = 1:dims,return.model = TRUE)



pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "project", label = TRUE) + NoLegend()
dev.off()


B <- B %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

B <- RunUMAP(
  B, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
B  <-
  FindNeighbors(
    B ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


B <-
  FindClusters(
    B,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(B, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(B, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(B, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(B, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(B, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(B) <- "sketch_snn_res.0.4"
B <- FindSubCluster(B, cluster = 7, graph.name = "sketch_snn", algorithm = 1,resolution = 0.3)
Idents(B) <- "sub.cluster"
DimPlot(B, group.by = "sub.cluster", label = TRUE, reduction = 'harmonyumap') + NoLegend()

### dotplot -------
major_marker_ls <- c(
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "ASC" = "MZB1",
  "ASC" = "XBP1",
  # "ASC" = "CD38",
  "ASC" = "SDC1",
  "ASC" = "IGHG1",
  "ASC" = "IGHA1",
  'naive_B' = 'IGHM',
  'naive_B' = 'BANK1',
  'naive_B' = 'FCRL1',
  'naive_B' = 'HLA-DQB1',
  'naive_B' = 'HLA-DQA1',
  'naive_B' = 'HVCN1',
  'naive_B' = 'IGHD',
  'naive_B' = 'FCER2',
  'naive_B' = 'IL4R',
  'memory_B'= 'CD27',
  'memory_B'= 'AIM2',
  'memory_B'= 'TXNIP',
  'memory_B'= 'GPR183',
  'memory_B'= 'TNFRSF13B',
  'memory_B'= 'CD22',
  'memory_B'= 'FCRL2',
  'memory_B'= 'ADAM28',
  'memory_B'= 'IGLL5',
  'memory_B'= 'P2RX5',
  'memory_B'= 'BLK',
  # 'cyclingB'='MKI67',
  # 'cyclingB'='STMN1',
  # 'cyclingB'='HMGB2',
  # 'cyclingB'='BIRC5',
  # 'plasma'= 'IGHA1',
  # 'plasma'= 'IGHG1',
  # 'plasma'= 'IGHG3',
  # 'plasma'= 'IGHG4',
  'GC B'= 'S1PI2',
  'GC B'='BCL6',
  'GC B'='RGS13',
  'GC B'='TCL1A',
  'GC B'= 'LRMP',
  'GC B'= 'SUGCT',
  'GC B'= 'MME',
  'GC B'='AICDA'
)




# dotplot
pdf(file.path(b_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(B, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(B, group.by = "sketch_snn_res.0.4", features = major_marker_ls)
DotPlot(B, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(B, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()

pdf(file.path(b_out_dir, str_glue("dotplot dims{dims} major_marker_2.pdf")), width = 15, height = 6)
DotPlot(B, group.by = "sketch_snn_res.0.2", features = c(
  "CD3D",
  "CD3E",
  "TRAC",
  "CD4",
  "CD8A",
  "MS4A1",
  "CD79A",
  "JCHAIN",
  "EPCAM",
  "PECAM1",
  "COL1A1",
  "RGS5",
  "TAGLN",
  "LYZ",
  'CD74',
  'CD14',
  'FCN1',
  'CD68',
  'C1QA',
  'KIT'
))
dev.off()

# # 2: T, 3: macrophage
# B2 <- B %>% subset(cells = colnames(B)[! B$sketch_snn_res.0.2 %in% c(2, 3)] )
# 
# DimPlot(B2, reduction = "harmonyumap")
# 
# B2 <- ScaleData(B2, vars.to.regress = 'nCount_RNA')
# B2 <- RunPCA(B2)
# 
# ElbowPlot(B2, ndims = 30)
# 
# B2 <- RunUMAP(B2, dims = 1:dims, return.model = TRUE)
# 
# B2 <- B2 %>% 
#   harmony::RunHarmony(
#     dims.use = 1:dims,
#     group.by.vars = "project",
#     reduction.save = "harmony"
#     # theta = 3, # default 2, bigger value -- more diversity per cluster
#     # lambda = 0.5, # default, smaller value -- more integrated
#   )
# 
# B2 <- RunUMAP(
#   B2, 
#   dims = 1:dims, 
#   reduction = "harmony", 
#   reduction.name = "harmonyumap", 
#   reduction.key  = "harmonyUMAP_",
#   return.model = TRUE
# )
# 
# B2  <-
#   FindNeighbors(
#     B2 ,
#     reduction = reduction.use,
#     dims = 1:dims,
#     k.param = 20
#   )
# B2 <-
#   FindClusters(
#     B2,
#     resolution = c(0.2, 0.4, 0.6, 0.8)
#   )
# 
# 
# DimPlot(B2, reduction = "harmonyumap", group.by = "sketch_snn_res.0.2", label = TRUE)
# 
# DotPlot(B2, group.by = "sketch_snn_res.0.2", features = c(
#   "CD3D",
#   "CD3E",
#   "TRAC",
#   "MS4A1",
#   "CD79A",
#   "JCHAIN",
#   "EPCAM",
#   "PECAM1",
#   "COL1A1",
#   "RGS5",
#   "TAGLN",
#   "LYZ",
#   'CD74',
#   'CD14',
#   'FCN1',
#   'CD68',
#   'C1QA',
#   'KIT'
# ))
# 
# DotPlot(B2, group.by = "sketch_snn_res.0.2", features = major_marker_ls)

### annotation ----

# cell_cluster <- c(
#   '0'  = "B_naive",
#   '1'  = "B_memory",
#   '2'  = "B_naive",
#   '3'  = "B_naive",
#   '4'  = "GC_B",
#   '5'  = "B_naive",
#   '6'  = "B_naive",
#   '7'  = "B_naive",
#   '8'  = "B_naive",
#   '9'  = "B_naive"
# )
# Idents(B) <- "sketch_snn_res.0.4"

cell_cluster <- c(
  '0'  = "B_naive",
  '1'  = "ASC",
  '2'  = "CD4_T",
  '3'  = "Macrophage",
  '4'  = "B_GC",
  '5'  = "B_naive",
  '6'  = "B_naive",
  '7'  = "B_memory"
)

Idents(B) <- "sketch_snn_res.0.2"

B <- RenameIdents(B, cell_cluster)
B$subtype <- Idents(B)


pdf(file.path(b_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(B, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf(file.path(b_out_dir, "dotplot celltype main_marker.pdf"), width = 10, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(B, group.by = "subtype", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()


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





#Fibroblast subtype--------------------
Fibroblast<-subset(all_data_cycle3,annotation3=='Fibroblast'|annotation3=='Epi0_Fibroblast')
### PCA --------
DefaultAssay(Fibroblast) <- "sketch"
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
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(Fibroblast)) %>%
  dplyr::pull("cell_id")


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
dims <- 10


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
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(Fibroblast, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "sub.cluster", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()
Idents(Fibroblast) <- "sketch_snn_res.0.6"
Fibroblast <- FindSubCluster(Fibroblast, cluster = 9, graph.name = "sketch_snn", algorithm = 1,resolution = 0.8)
Idents(Fibroblast) <- "sub.cluster"
### dotplot -------
major_marker_ls <- c(
  Epithelial  = "EPCAM", 
  # Epithelial  = "KRT5",
  Epithelial  = "KRT8",
  Epithelial  = "S100A8",
  Epithelial  = "KRT18",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  
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
major_marker_ls <- c(
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
  'rCAF'='CCL21',
  'hsptcaf'='HSPH1',
  'hsptcaf'='VEGFA'
)

# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.4", features = major_marker_ls)
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.5", features = major_marker_ls)
DotPlot(Fibroblast, group.by = "sub.cluster", features = major_marker_ls)
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()

### annotation ----

cell_cluster <- c(
  '0'  = "rCAF",
  '1'  = "iCAF",
  '2'  = "mCAF",
  '3'  = "apCAF",
  '4'  = "vCAF",
  '5'  = "Epithelial",
  '6'  = "mCAF",
  '7'  = "mCAF",
  '8'  = "Epithelial",
  '9_0'  = "apCAF",
  '9_1'  = "iCAF",
  '9_2'  = "mCAF",
  '9_3'  = "mCAF",
  '9_4'  = "vCAF",
  '9_5'  = "iCAF",
  '10'  = "mast",
  '11'  = "Epithelial"
)


Idents(Fibroblast) <- "sub.cluster"

Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype <- Idents(Fibroblast)


pdf(str_glue("Dimplot dim{dims} annotation.pdf"), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf("dotplot celltype main_marker.pdf", width = 10, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(Fibroblast, group.by = "subtype", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
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


##映射----
library(tibble)
library(data.table)
meta.data <-
  bind_rows(list(CD4T[['subtype_full']] %>% rownames_to_column('cells'),
                 CD8T[['subtype_full']] %>% rownames_to_column('cells'),
                 B[['subtype_full']] %>% rownames_to_column('cells'), 
                 Fibroblast[['subtype_full']] %>% rownames_to_column('cells'))) %>%
  as.data.frame() %>%
  right_join(
    all_data_cycle3[['annotation3']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation4 = ifelse(
      (as.character(annotation3)=='CD4-T'|as.character(annotation3)=='CD8-T'|
         as.character(annotation3)=='B'|as.character(annotation3)=='Epi0_B'|
         as.character(annotation3)=='CD83-dim B'|as.character(annotation3)=='CD83 B'| as.character(annotation3)=='Plasma'| 
         as.character(annotation3)=='Fibroblast'|as.character(annotation3)=='Epi0_Fibroblast'),
      
      as.character(subtype_full),
      as.character(annotation3))
  ) %>%
  dplyr::select(cells, annotation4) %>%
  column_to_rownames('cells')


all_data_cycle3 <- AddMetaData(all_data_cycle3, meta.data['annotation4'])

# # replace B cell annotation 2024.11.28
# 2024.12.24 mapping code update.
# meta_b <- B[['subtype_full']] %>% rownames_to_column('cells')
# 
# meta_b_full <-
#   all_data_cycle3[['annotation4']] %>%
#   rownames_to_column('cells') %>%
#   left_join(meta_b, by = "cells") %>%
#   mutate(annotation4 = ifelse(
#     is.na(subtype_full),
#     annotation4,
#     subtype_full
#   )) %>%
#   dplyr::select(annotation4, cells) %>%
#   column_to_rownames("cells")
# 
# all_data_cycle3 <- AddMetaData(all_data_cycle3, meta_b_full)


pdf(file.path("mapping_annotation umap.pdf"), height = 6, width = 12)
DimPlot(all_data_cycle3, group.by = "annotation4", label = TRUE, reduction = 'harmonyumap',pt.size = 2,repel = TRUE) 
#DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# saveRDS(all_data_cycle3, "all_data_cycle3_newdata.rds")


# all_data <- readRDS('all_data_cycle3_newdata.rds')
# epithelial cells in substype should not annotate as epithelial cells
# use prefix: Fibro-like or Immune-like instead
all_data$annotation4[which(all_data$annotation2 == 'Epi0_Fibroblast' & all_data$annotation4 == 'Epithelial')] <- 'Fibro_like Epithelial'
all_data$annotation4[which(all_data$annotation2 == 'Epi0_CD4T'       & all_data$annotation3 == 'Epithelial')] <- 'Immune_like Epithelial'
all_data$annotation4[which(all_data$annotation2 == 'Epi0_CD83_T'     & all_data$annotatio3 == 'Epithelial')]  <- 'Immune_like Epithelial'
all_data$annotation4[which(all_data$annotation2 == 'Epi0_Epi0_T'     & all_data$annotation3 == 'Epithelial')] <- 'Immune_like Epithelial'


# 合并mast和Mast，CD4T改为CD4T-like, 2025.1.3
all_data$annotation4[all_data$annotation4 %in% c("mast", "Mast")] <- "Mast"
all_data$annotation4[all_data$annotation4 == 'CD4-T'] <- "CD4-T_like"



# saveRDS(all_data, file = file.path(out_dir, "subtype", "all_data_cycle3_1205.rds"))
saveRDS(all_data[[]], file = file.path(out_dir, "subtype", "all_data_cycle3_1205_meta.rds"))


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
HLA<-c('HLA-A','HLA-B','HLA-C','HLA-DRA','HLA-DRB1','HLA-DQA1','HLA-DQB1','TAP1','TAP2')
tumor_cell_LN_metastasis<-c('VEGFC','VEGFD','LYVE1','PROX1','CCL1','CCR8','FGF2','PROM1','HGF','MBTPS1')
IFN<-c('JAK1','JAK2','TYK2','IFNG','IFNB1','STAT1','STAT2','STAT3','IFNAR1','IFNGR1','ISG15','TMEM173','TLR6')
cancer_stem_cell_and_EMT<-c('CD44','CDH2','VIM')
Inflammation<-c('IL6','CXCL2','TNF','MAP3K8')


pathway<-list(HLA,tumor_cell_LN_metastasis,IFN,cancer_stem_cell_and_EMT,Inflammation)
names(pathway)<-c('HLA','tumor_cell_LN_metastasis','IFN','cancer_stem_cell_and_EMT','Inflammation')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  all_data_cycle3 <- AddModuleScore(JoinLayers(all_data_cycle3), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=50)
}

all_data_meta <- all_data_cycle3[[]] %>% rownames_to_column('cell_id')

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
#  sp_out_dir<-'/home/rstudio/dat_analysis_cycle3/subtype/newplot/HLA'

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
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF","#802268FF",'grey90'))+
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
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
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
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF","#802268FF",'grey90'))+
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
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF","#802268FF",'grey90'))+
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
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF","#802268FF",'grey90'))+
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
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF","#802268FF",'grey90'))+
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