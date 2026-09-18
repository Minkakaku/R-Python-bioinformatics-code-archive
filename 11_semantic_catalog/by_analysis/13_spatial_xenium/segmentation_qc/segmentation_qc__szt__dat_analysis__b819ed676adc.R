library(Seurat)
library(tidyverse)
library(BPCells)


out_dir <- 'dat_analysis_cycle3'
dir.create(out_dir)


# read data --------------
all_data <- readRDS("all_data_merge.rds") 
anno_data <- read.csv(file.path(out_dir, 'annotation.csv')) 

all_data <-
  all_data %>%
  AddMetaData(
    all_data[['project']] %>%
      rownames_to_column('cell_id') %>%
      left_join(anno_data, by = 'project') %>%
      dplyr::select(cell_id, group, group2) %>%
      column_to_rownames('cell_id')
  ) 

all_data$project %>% 
  table %>% 
  enframe('project', 'count') %>%
  arrange(count) %>%
  as.data.frame()
# 2200,5900 - 370000

all_data <- all_data %>% 
  SketchData(
    ncells = 5000, # 5000 per layer
    method = 'LeverageScore', # 'LeverageScore',  
    sketched.assay = "sketch"
  )


stopifnot(DefaultAssay(all_data) == 'sketch')


# PCA --------
DefaultAssay(all_data) <- "sketch"
all_data <- FindVariableFeatures(all_data)
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')

all_data <- RunPCA(all_data)


## pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()

set.seed(1234)
suset_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 16:20)
dev.off()
# dimension reduction --------
dims <- 12


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


all_data <- RunUMAP(all_data, dims = 1:dims)



pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE) + NoLegend()
dev.off()


all_data <- all_data %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

all_data <- RunUMAP(
  all_data, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


# cluster -------
all_data <-
  FindNeighbors(
    all_data,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# dotplot -------
major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  # "Fibroblasts" = "ACTA2",
  
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "MYH11", 
  "Smooth muscle cells" = "ACTA2",
  
  "Breast glandular cells" = "SERPINA3", 
  "Breast glandular cells" = "DSC2",
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  # "Plasma" = "JCHAIN",
  
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "CD4",
  "T cells" = "CD8A",
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  "Myeloid" = "ITGAM",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A",
  
  
  'Fat cell' = 'ADIPOQ',
  'Fat cell' = 'UCP1'
)


# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()

# annotation ----

cell_cluster <- c(
  '0'  = "CD4T",
  '1'  = "CD83-dim B",
  '2'  = "Epithelial",
  '3'  = "Mono-like Macrophage",
  '4'  = "Fibroblast",
  '5'  = "Endothelial",
  '6'  = "CD8T",
  '7'  = "Plasma",
  '8'  = "B",
  '9'  = "CD4 Endothelial",
  '10' = "Smooth muscle cells",
  '11' = "CD4T-Macrophage",
  '12' = "Fat cells",
  '13' = "pDCs",
  '14' = "CD83 B"
)


Idents(all_data) <- "sketch_snn_res.0.8"

all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


# find markers -------
Idents(all_data) <- "sketch_snn_res.0.8"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
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
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(all_data))[[1]]


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

# project data ------
all_data <- 
  all_data %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "pca.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'
    dims = 1:dims, # <---- harmony只有13个dims
    
    refdata = list(annotation_full = "annotation")
  )


DefaultAssay(all_data) <- 'RNA'

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation - harmony - projected.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

# saveRDS(all_data, file.path(out_dir, "all_data_projected.rds"))
all_data <- readRDS(file.path(out_dir, "all_data_projected.rds"))

# all_data$annotation_full %>% table



# sub cluster --------
# all_data_diet <-
#   all_data %>%
#   DietSeurat(
#     layers = Layers(all_data, search = "counts")
#   )

sc_epi <- all_data[, all_data$annotation_full == 'Epithelial']
sc_epi <- DietSeurat(sc_epi)
sc_epi[['sketch']] <- NULL

# sc_epi <- 
#   sc_epi %>% 
#   FindVariableFeatures() %>%
#   NormalizeData()


# sc_epi[[c('project','group','group2', 'annotation_full')]] %>%
#   group_by(project,group,group2, annotation_full) %>%
#   summarise(count = n()) %>%
#   
#   ggplot(aes(x = project, y = count, fill = annotation_full)) +
#   geom_col() +
#   geom_label(aes(label = group2))

sc_epi <- sc_epi %>%
  FindVariableFeatures() %>%
  NormalizeData() %>%
  ScaleData(vars.to.regress = 'nCount_RNA')

## pca -----
sc_epi <- sc_epi %>% RunPCA()

pdf(file.path(out_dir, 'pca plot-epi.pdf'), height = 6, width = 8)
ElbowPlot(sc_epi, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_epi <- 
  sc_epi[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sc_epi) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(out_dir, 'epi dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi[, suset_cells_epi], dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'epi dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi[, suset_cells_epi], dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'epi dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi[, suset_cells_epi], dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'epi dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi[, suset_cells_epi], dims = 16:20)
dev.off()


## dimension reduction --------
epi_dims <- 7


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


sc_epi <- RunUMAP(sc_epi, dims = 1:epi_dims)


pdf(file.path(out_dir, str_glue("epi Dimplot dim{epi_dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_epi, group.by = "project", label = TRUE) + NoLegend()
dev.off()


sc_epi <- sc_epi %>% 
  harmony::RunHarmony(
    dims.use = 1:epi_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_epi <- RunUMAP(
  sc_epi, 
  dims = 1:epi_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(out_dir, str_glue("epi Dimplot dim{epi_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_epi, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'pca' # 'harmony' 'pca'

## cluster ------
sc_epi <-
  FindNeighbors(
    sc_epi,
    reduction = reduction.use,
    dims = 1:epi_dims,
    k.param = 20
  )

sc_epi <-
  FindClusters(
    sc_epi,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(out_dir, str_glue("epi Dimplot dim{epi_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(sc_epi, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_epi, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_epi, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(sc_epi, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


pdf(file.path(out_dir, str_glue("epi Dimplot dim{epi_dims} group2 - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(sc_epi, group.by = "group2")
dev.off()


## dotplot -------
major_marker <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  Tumor = 'CEACAM6',
  
  "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  # "Fibroblasts" = "ACTA2",
  
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "MYH11", 
  "Smooth muscle cells" = "ACTA2",
  
  "Breast glandular cells" = "SERPINA3", 
  "Breast glandular cells" = "DSC2",
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  # "Plasma" = "JCHAIN",
  
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "CD4",
  "T cells" = "CD8A",
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  "Myeloid" = "ITGAM",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A"
  
)


# dotplot
pdf(file.path(out_dir, str_glue("epi dotplot dims{epi_dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(sc_epi, group.by = "RNA_snn_res.0.2", features = major_marker)
DotPlot(sc_epi, group.by = "RNA_snn_res.0.4", features = major_marker)
DotPlot(sc_epi, group.by = "RNA_snn_res.0.6", features = major_marker)
DotPlot(sc_epi, group.by = "RNA_snn_res.0.8", features = major_marker)
dev.off()


### cluster 0 not a epithelial cluster -----

sc_epi0 <- sc_epi[, sc_epi$RNA_snn_res.0.2 == '0']


sc_epi0 <- 
  sc_epi0 %>%
  FindVariableFeatures() %>%
  NormalizeData() %>%
  ScaleData(vars.to.regress = 'nCount_RNA') %>%
  RunPCA()

pdf(file.path(out_dir, 'pca plot-epi0.pdf'), height = 6, width = 8)
ElbowPlot(sc_epi0, ndims = 30)
dev.off()

set.seed(1234)
suset_cells_epi0 <- 
  sc_epi0[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(sc_epi0) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(out_dir, 'epi0 dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi0[, suset_cells_epi0], dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'epi0 dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi0[, suset_cells_epi0], dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'epi0 dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi0[, suset_cells_epi0], dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'epi0 dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(sc_epi0[, suset_cells_epi0], dims = 16:20)
dev.off()

### dimension reduction --------
epi0_dims <- 7


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


sc_epi <- RunUMAP(sc_epi0, dims = 1:epi0_dims)


pdf(file.path(out_dir, str_glue("epi0 Dimplot dim{epi0_dims} project.pdf")), height = 6, width = 8)
DimPlot(sc_epi0, group.by = "project", label = TRUE) + NoLegend()
dev.off()


sc_epi0 <- sc_epi0 %>% 
  harmony::RunHarmony(
    dims.use = 1:epi0_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

sc_epi0 <- RunUMAP(
  sc_epi0, 
  dims = 1:epi0_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  # return.mcodel = TRUE
)


pdf(file.path(out_dir, str_glue("epi0 Dimplot dim{epi0_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_epi0, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use.epi0 <- 'harmony' # 'harmony' 'pca'

### cluster -----
sc_epi0 <-
  FindNeighbors(
    sc_epi0,
    reduction = reduction.use.epi0,
    dims = 1:epi0_dims,
    k.param = 20
  )

sc_epi0 <-
  FindClusters(
    sc_epi0,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(out_dir, str_glue("epi0 Dimplot dim{epi0_dims} cluster - {reduction.use.epi0}.pdf")), height = 6, width = 8)
DimPlot(sc_epi0, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi0, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi0, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi0, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_epi0, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi0, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi0, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(sc_epi0, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


pdf(file.path(out_dir, str_glue("epi0 Dimplot dim{epi0_dims} group2 - {reduction.use.epi0}.pdf")), height = 6, width = 8)
DimPlot(sc_epi0, group.by = "group2")
dev.off()


### dotplot ------

pdf(file.path(out_dir, str_glue("epi0 dotplot dims{epi0_dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(sc_epi0, group.by = "RNA_snn_res.0.2", features = major_marker)
DotPlot(sc_epi0, group.by = "RNA_snn_res.0.4", features = major_marker)
DotPlot(sc_epi0, group.by = "RNA_snn_res.0.6", features = major_marker)
DotPlot(sc_epi0, group.by = "RNA_snn_res.0.8", features = major_marker)
dev.off()


### annotation ----

cell_cluster_epi0 <- c(
  '0'  = "Epi0_Fibroblast",
  '1'  = "Epi0_B",
  '2'  = "Epi0_T",
  '3'  = "Epi0_CD83_T",
  '4'  = "Epi0_Mast cell",
  '5'  = "Epi0_CD4T"
)


Idents(sc_epi0) <- "RNA_snn_res.0.2"

sc_epi0 <- RenameIdents(sc_epi0, cell_cluster_epi0)
sc_epi0$annotation_epi0 <- Idents(sc_epi0)

pdf(file.path(out_dir, str_glue("epi0 Dimplot dim{epi0_dims} annotation - harmony.pdf")), height = 6, width = 8)
DimPlot(sc_epi0, group.by = "annotation_epi0", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(sc_epi0, group.by = "annotation_epi0", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


### combined epithelial cells to all_data ----
# 只将annotation_epi0更新至all_data
# sc_epi的其他群还是上皮细胞

sc_meta_epi0 <-
  sc_epi0[['annotation_epi0']] %>%
  rownames_to_column("cell_id")

sc_meta_all <- 
  all_data[['annotation_full']] %>% 
  rownames_to_column("cell_id") %>%
  left_join(sc_meta_epi0, by = "cell_id")

sc_meta_all$annotation2 <- sc_meta_all$annotation_epi0 %>% as.character()
sc_meta_all$annotation2[is.na(sc_meta_all$annotation2)] <-
  sc_meta_all$annotation_full[is.na(sc_meta_all$annotation2)] %>%
  as.character()

sc_meta_all <- sc_meta_all %>% column_to_rownames("cell_id")

# saveRDS(sc_meta_all, file = file.path(out_dir, "sc_meta_all.rds"))


### export epithelial count ---------
all_data2[[c("project", "group", "group2", "annotation2")]] %>%
  group_by(project, group, group2, annotation2) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = "annotation2", values_from = "count") %>%
  mutate(total_count = sum(c_across(B:pDCs), na.rm = T)) %>%
  dplyr::select(project, group, group2, Epithelial,Epi0_Fibroblast ,total_count) %>%
  openxlsx::write.xlsx(
    file.path(out_dir, "epithelial cells count.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

### average expression heatmap ------
avg_expr <- AverageExpression(all_data2)


# all_data annotation2 -------
all_data2 <- all_data %>% AddMetaData(sc_meta_all['annotation2'])


# saveRDS(all_data2, file = file.path(out_dir, 'all_data_cycle3.rds'))

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation - harmony - projected - annotation2.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation2", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

p_cell_freq_col <-
  all_data2[[]] %>%
  group_by(project, group, group2, annotation2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  arrange(group, group2) %>%
  as.data.frame() %>%
  mutate(project = factor(project, levels = unique(project))) %>%
  
  ggplot(aes()) +
  geom_col(aes(x = project, y = freq, fill = annotation2)) +
  ggsci::scale_fill_igv() +
  ggnewscale::new_scale_fill() +
  geom_col(aes(x = project, y = -0.2, fill = group), data = distinct(all_data2[[c('project', 'group', 'group2')]])) +
  ggnewscale::new_scale_fill() +
  geom_col(aes(x = project, y = -0.1, fill = group2), data = distinct(all_data2[[c('project', 'group', 'group2')]])) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  guides(
    fill = guide_legend(order = 3),
    fill_new = guide_legend(order = 2),
    fill_new_new = guide_legend(order = 1)
  )

ggsave(
  p_cell_freq_col,
  filename = file.path(out_dir, "cell proportion barplot - cycle3.pdf"),
  height = 6, 
  width  = 12 
)


# spatial plot -----------------
sp_out_dir <- file.path(out_dir, 'spatial plot') # dir.create(sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

all_data_meta <- all_data2[[]] %>% rownames_to_column('cell_id')

sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation2, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
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



# CEACAM6 expression -------
p_ceacam6_expr <-
  FeaturePlot(
    all_data2,
    features = "CEACAM6",
    split.by = "project",
    combine  = FALSE
  ) %>%
  patchwork::wrap_plots(ncol = 5, guides = "collect")


ggsave(
  filename = file.path(out_dir, "CEACAM6 expression.pdf"),
  height = 20,
  width  = 20
)






