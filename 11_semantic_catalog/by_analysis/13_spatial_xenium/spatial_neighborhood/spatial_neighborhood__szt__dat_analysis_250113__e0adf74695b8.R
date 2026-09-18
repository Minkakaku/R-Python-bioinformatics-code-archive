library(Seurat)
library(tidyverse)


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 10Gb



out_dir <- "dat_analysis_cycle3"


all_data <- readRDS(file.path(out_dir, "subtype", "all_data_cycle3_1205.rds"))
all_data_meta <- all_data[[]]

#1.  cell subcluster ------------------
#1.1 T细胞亚群中分出来的DCs合并到pDCs中 -------------------------
all_data$annotation5 <- all_data$annotation4 
all_data$annotation5[all_data$annotation5 == "DCs"] <- "pDCs"


#1.2 B细胞亚群中分出来CD4 T-like合并到CD4 T中分析，重新分亚群 -------
cd4_out_dir <- file.path(out_dir, "subtype", "cd4_cycle2_250113")
dir.create(cd4_out_dir)


sp_cd4t <- subset(all_data, annotation3 == 'CD4-T' | annotation5 == "CD4-T_like")

### PCA --------
DefaultAssay(sp_cd4t) <- "sketch"
sp_cd4t <- FindVariableFeatures(sp_cd4t)
sp_cd4t <- ScaleData(sp_cd4t, vars.to.regress = 'nCount_RNA')

sp_cd4t <- RunPCA(sp_cd4t)


### pca inspector ------
pdf(file.path(cd4_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(sp_cd4t, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  sp_cd4t[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sp_cd4t)) %>%
  dplyr::pull("cell_id")


pdf(file.path(cd4_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sp_cd4t[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(cd4_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sp_cd4t[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(cd4_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sp_cd4t[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(cd4_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sp_cd4t[, suset_cells], dims = 16:20)
dev.off()
### dimension reduction --------
dims <- 20

sp_cd4t <- RunUMAP(sp_cd4t, dims = 1:dims,return.model = TRUE)



pdf(file.path(cd4_out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(sp_cd4t, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


sp_cd4t <- sp_cd4t %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sp_cd4t <- RunUMAP(
  sp_cd4t, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(cd4_out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sp_cd4t, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
sp_cd4t  <-
  FindNeighbors(
    sp_cd4t ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


sp_cd4t <-
  FindClusters(
    sp_cd4t,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(cd4_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(sp_cd4t, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_cd4t, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_cd4t, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_cd4t, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sp_cd4t, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(sp_cd4t, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(sp_cd4t, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(sp_cd4t, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


# Idents(sp_cd4t) <- "sketch_snn_res.0.4"
# sp_cd4t <- FindSubCluster(sp_cd4t, cluster = 3, graph.name = "sketch_snn", algorithm = 1,resolution = 0.3)
# Idents(sp_cd4t) <- "sub.cluster"
# DimPlot(sp_cd4t, group.by = "sub.cluster", label = TRUE, reduction = 'harmonyumap') + NoLegend()

### dotplot -------
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
sp_cd4t_markers <- c(
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
pdf(file.path(cd4_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(sp_cd4t, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(sp_cd4t, group.by = "sketch_snn_res.0.4", features = major_marker_ls)
DotPlot(sp_cd4t, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(sp_cd4t, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()

### annotation ----

cell_cluster <- c(
  '0'  = "CD4_Tmemory",
  '1'  = "CD4T_Tcm",
  '2'  = "Treg",
  '3'  = "CD4_Tnaive",
  '4'  = "Tfh",
  '5'  = "proliferating-T",
  '6'  = "CD4_Tnaive",
  '7'  = "CD4_Tnaive",
  '8'  = "CD4_Tmemory"
)


Idents(sp_cd4t) <- "sketch_snn_res.0.4"

sp_cd4t <- RenameIdents(sp_cd4t, cell_cluster)
sp_cd4t$subtype <- Idents(sp_cd4t)


pdf(file.path(cd4_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(sp_cd4t, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

DotPlot(sp_cd4t, group.by = "sub.cluster", features = major_marker_ls, dot.scale = 10) + theme_90angle
pdf(file.path(cd4_out_dir, "dotplot celltype main_marker.pdf"), width = 14, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(sp_cd4t, group.by = "subtype", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()


### project data ------
sp_cd4t <- 
  sp_cd4t %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'harmonyumap
    dims = 1:dims, 
    
    refdata = list(subtype_full = "subtype")
  )


DefaultAssay(sp_cd4t) <- 'RNA'

pdf(file.path(cd4_out_dir, str_glue("Dimplot dim{dims} annotation - projected.pdf")), height = 6, width = 8)
DimPlot(sp_cd4t, group.by = "subtype_full", label = TRUE, reduction = 'full.harmonyumap.1') + NoLegend()
dev.off()

# saveRDS(sp_cd4t, file.path(cd4_out_dir, "CD4_subtype.rds"))







# 1.3 从B/T细胞亚群中分出的Macrophage与mono-like macrophage合并，对Macrophage分亚群 -------------
mye_out_dir <- file.path(out_dir, "subtype", "mye_cycle2_250113")
dir.create(mye_out_dir)

sp_mye <- subset(all_data, annotation5 %in% c("Macrophage", "Mono-like Macrophage"))

### PCA --------
DefaultAssay(sp_mye) <- "sketch"
sp_mye <- FindVariableFeatures(sp_mye)
sp_mye <- ScaleData(sp_mye, vars.to.regress = 'nCount_RNA')

sp_mye <- RunPCA(sp_mye)


### pca inspector ------
pdf(file.path(mye_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(sp_mye, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  sp_mye[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sp_mye)) %>%
  dplyr::pull("cell_id")


pdf(file.path(mye_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sp_mye[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(mye_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sp_mye[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(mye_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sp_mye[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(mye_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sp_mye[, suset_cells], dims = 16:20)
dev.off()
### dimension reduction --------
dims <- 20

sp_mye <- RunUMAP(sp_mye, dims = 1:dims,return.model = TRUE)



pdf(file.path(mye_out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(sp_mye, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


sp_mye <- sp_mye %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sp_mye <- RunUMAP(
  sp_mye, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(mye_out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sp_mye, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
sp_mye  <-
  FindNeighbors(
    sp_mye ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


sp_mye <-
  FindClusters(
    sp_mye,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(mye_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(sp_mye, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_mye, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_mye, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sp_mye, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(sp_mye, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(sp_mye, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(sp_mye, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(sp_mye, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


# Idents(sp_mye) <- "sketch_snn_res.0.4"
# sp_mye <- FindSubCluster(sp_mye, cluster = 3, graph.name = "sketch_snn", algorithm = 1,resolution = 0.3)
# Idents(sp_mye) <- "sub.cluster"
# DimPlot(sp_mye, group.by = "sub.cluster", label = TRUE, reduction = 'harmonyumap') + NoLegend()

### dotplot -------
sp_mye_markers <- c(
  'EPCAM',
  'CD3D',
  'CD4',
  'CD8A',
  'CD19',
  "MS4A1",
  "BANK1", # naive
  'NKG7',
  'NCAM1',
  'PECAM1',
  
  'IGHD',
  'CD27',
  'CD38',
  'CD24',
  
  'monocytes' = 'FCN1', 
  # 'monocytes' = 'S100A9', 
  'monocytes' = 'CD14', 
  'monocytes' = 'FCGR3A', 
  
  cDC1 = 'XCR1',
  cDC1 = 'CLEC9A',

  cDC2 = 'CD1C',
  cDC2 = 'CLEC10A',
  
  mast = "TPSAB1",
  mast = "CPA3",
  mast = "KIT",
  
  'macrophages' = 'CD68', 
  'macrophages' = 'CD163', 
  'macrophages' = 'C1QA', 
  'macrophages' = 'C1QB', 
  'macrophages' = 'MS4A7', 
  
  'macrophages' = 'NFKB2',  #'NLPR3', 
  'macrophages' = 'PLTP', 
  'macrophages' = 'IL1B', 
  'macrophages' = 'C1QC', 
  'macrophages' = 'SPP1', 
  
  'proliferation'='MKI67',
  'proliferation'='TOP2A'
)


# dotplot
pdf(file.path(mye_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(sp_mye, group.by = "sketch_snn_res.0.2", features = sp_mye_markers)
DotPlot(sp_mye, group.by = "sketch_snn_res.0.4", features = sp_mye_markers)
DotPlot(sp_mye, group.by = "sketch_snn_res.0.6", features = sp_mye_markers)
DotPlot(sp_mye, group.by = "sketch_snn_res.0.8", features = sp_mye_markers)
dev.off()

pdf(file.path(mye_out_dir, str_glue("dotplot dims{dims} B_cell_marker.pdf")), width = 25, height = 6)
DotPlot(
  sp_mye,
  group.by = "sketch_snn_res.0.2",
  features = c(
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

    'GC B'= 'S1PI2',
    'GC B'='BCL6',
    'GC B'='RGS13',
    'GC B'='TCL1A',
    'GC B'= 'LRMP',
    'GC B'= 'SUGCT',
    'GC B'= 'MME',
    'GC B'='AICDA'
  )
)
dev.off()

### annotation ----

mye_cell_cluster <- c(
  '0'  = "Macrophage",
  '1'  = "Macrophage",
  '2'  = "Macrophage",
  '3'  = "B-naive",
  '4'  = "cDC1",
  '5'  = "Monocyte",
  '6'  = "Mast"
)


Idents(sp_mye) <- "sketch_snn_res.0.2"

sp_mye <- RenameIdents(sp_mye, mye_cell_cluster)
sp_mye$subtype <- Idents(sp_mye)


pdf(file.path(mye_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(sp_mye, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# DotPlot(sp_mye, group.by = "sub.cluster", features = major_marker_ls, dot.scale = 10) + theme_90angle
pdf(file.path(mye_out_dir, "dotplot celltype main_marker.pdf"), width = 10, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(sp_mye, group.by = "subtype", features = sp_mye_markers, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()


### project data ------
sp_mye <- 
  sp_mye %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'harmonyumap
    dims = 1:dims, 
    
    refdata = list(subtype_full = "subtype")
  )


DefaultAssay(sp_mye) <- 'RNA'

pdf(file.path(mye_out_dir, str_glue("Dimplot dim{dims} annotation - projected.pdf")), height = 6, width = 8)
DimPlot(sp_mye, group.by = "subtype_full", label = TRUE, reduction = 'full.harmonyumap.1') + NoLegend()
dev.off()

# saveRDS(sp_mye, file.path(mye_out_dir, "macrophage_subtype.rds"))

## find markers -------
Idents(sp_mye) <- "sketch_snn_res.0.2"
all_markers <-
  FindAllMarkers(
    JoinLayers(sp_mye),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )


all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  pull("gene")

# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(sp_mye))[[1]]


all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

avg_expr[all_de_genes, ] %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_row = avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
    cluster_col = FALSE,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_col = NA
  )


# 1.4 cominbe all annotation ----------------------
meta_cd4t_mye_anno <-
  bind_rows(
    sp_cd4t[['subtype_full']],
    sp_mye[['subtype_full']]
  )


meta_cd4t_mye_anno_full <-
  all_data[['annotation5']] %>% 
  merge(meta_cd4t_mye_anno, by = "row.names") %>%
  mutate(annotation5 = ifelse(is.na(subtype_full), annotation5, subtype_full)) %>%
  column_to_rownames("Row.names") %>%
  dplyr::select(annotation5)


all_data <- all_data %>% AddMetaData(meta_cd4t_mye_anno_full)

saveRDS(all_data, file.path(out_dir, "subtype", "all_data_cycle3_250219.rds"))

# 1.5 plot ------
p_umap_all_cycle2 <-
  all_data %>%
  DimPlot(group.by = "annotation5", reduction = "harmonyumap", label = TRUE, repel = TRUE) +
  ggsci::scale_color_igv()

ggsave(
  p_umap_all_cycle2,
  filename = file.path(file.path(out_dir, "subtype", "umap plot of all cells - cd4t and macrophage.pdf")),
  height = 6,
  width  = 12
)

p_col_anno <-
  all_data[[c('group', 'annotation5')]] %>% 
  group_by(annotation5) %>%
  summarise(count = n()) %>%
  # group_by(group) %>%
  mutate(freq = count / sum(count)) %>%
  
  ggplot(aes(x = 1, fill = annotation5, y = freq)) +
  geom_col() +
  ggsci::scale_fill_igv()


all_data[[c('group', 'annotation5')]] %>% 
  group_by(annotation5) %>%
  summarise(count = n()) %>%
  openxlsx::write.xlsx(
    file.path(file.path(out_dir, "subtype", "cell number - cd4t and macrophage.xlsx")),
    asTable = TRUE
  )



