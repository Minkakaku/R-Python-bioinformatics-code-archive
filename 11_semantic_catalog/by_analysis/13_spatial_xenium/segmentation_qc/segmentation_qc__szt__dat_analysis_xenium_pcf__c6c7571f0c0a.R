library(tidyverse)
library(Seurat)


dr_dir <- "Xenium_PCF_dimension_reduction"
dir.create(dr_dir)

# # library --------------
# options(java.parameters = '-Xmx10g') # options('java.parameters')
# library(rJava)
# library(RBioFormats)
# source("xenium.R")

# multi workers  -----
library(future)
plan("multisession", workers = 10, gc = TRUE)

options(future.globals.maxSize = 20 * 1024 ^ 3)

nbrOfFreeWorkers()

# data import ----------------

# xenium
all_data_projected <- readRDS("all_data_projected.rds")
all_data_xenium_ls <- SplitObject(all_data_projected, split.by = "project")

# pca
all_data <- readRDS("all_data.pcf.rds")
all_data_pcf_ls <- SplitObject(all_data, split.by = "project")

# match table
sample_match_table <-
   data.frame(
     xenium = c("ID1", "ID10", "ID100"),
     pcf    = c('sample1', 'sample2', 'sample3')
   )


# combined pcf data
all_data_xenium_ls2 <- 
  all_data_xenium_ls %>% 
  imap(function(xenium, name){
    pcf_name <- sample_match_table$pcf[sample_match_table$xenium == name]
    stopifnot(length(pcf_name) == 1)
    pcf <- all_data_pcf_ls[[pcf_name]]
    
    # cell_id prefix 
    xenium_prefix <- colnames(xenium)[1] %>% str_split_1("_") %>% .[1]
    pcf_prefix    <- colnames(pcf)[1]    %>% str_split_1("_") %>% .[1]
    
    pcf_name_temp <- colnames(pcf) %>% str_remove(str_glue("^{pcf_prefix}_"))
    pcf_name_temp <- str_c(xenium_prefix, pcf_name_temp, sep = "_")
    
    # rename pcf
    colnames(pcf) <- pcf_name_temp
    rownames(pcf) <- str_c("PR-", rownames(pcf))
    
    # intersect of cell_id
    inter_cells <- intersect(colnames(xenium), colnames(pcf))
    
    # construct counts
    xenium_count <- LayerData(xenium, layer = Layers(xenium, "counts"))
    pcf_count    <- LayerData(pcf,    layer = Layers(pcf, "counts"))
    
    new_xenium_count <- rbind(
      xenium_count[, inter_cells],
      pcf_count[, inter_cells]
    )
    
    # construct data
    xenium_data <- LayerData(xenium, layer = Layers(xenium, "data"))
    pcf_data    <- LayerData(pcf,    layer = Layers(pcf, "data"))
    
    new_xenium_data <- rbind(
      xenium_data[, inter_cells],
      pcf_data[, inter_cells]
    )
    
    # assay
    new_assay <- 
      CreateAssay5Object(
        counts = new_xenium_count,
        data   = new_xenium_data
      )
    
    xenium[["xenium_pcf"]] <- new_assay
    xenium
  })
  
# merge seurat
all_data_xenium_pcf <- merge(all_data_xenium_ls2[[1]], all_data_xenium_ls2[-1])
all_data_xenium_pcf[['sketch']] <- NULL
DefaultAssay(all_data_xenium_pcf) <- "xenium_pcf"

all_data_xenium_pcf <- all_data_xenium_pcf %>% FindVariableFeatures()

# sketch data --------------
all_data_xenium_pcf <-
  SketchData(
    all_data_xenium_pcf,
    ncells         = 10000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score.pcf",
    seed = 123L
  )


## dimension reduction - pca ------------

DefaultAssay(all_data_xenium_pcf) <- "sketch"
all_data_xenium_pcf <- FindVariableFeatures(all_data_xenium_pcf, nfeatures = 3000)
all_data_xenium_pcf <- ScaleData(all_data_xenium_pcf, vars.to.regress = 'nCount_RNA')

all_data_xenium_pcf <- RunPCA(all_data_xenium_pcf)

pdf(file.path(dr_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data_xenium_pcf, ndims = 30)
dev.off()


set.seed(1234)
suset_cells <- 
  all_data_xenium_pcf[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data_xenium_pcf) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(dr_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data_xenium_pcf[, suset_cells], dims = 1:5)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data_xenium_pcf[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data_xenium_pcf[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data_xenium_pcf[, suset_cells], dims = 16:20)
dev.off()


# dimension reduction --------
DefaultAssay(all_data_xenium_pcf) <- "sketch"


dims <- 10 # 20


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


all_data_xenium_pcf <- RunUMAP(all_data_xenium_pcf, dims = 1:dims, return.model = TRUE)


pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data_xenium_pcf, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


all_data_xenium_pcf <- all_data_xenium_pcf %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data_xenium_pcf <- RunUMAP(
  all_data_xenium_pcf,
  dims = 1:dims,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_", 
  return.model = TRUE
)

pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data_xenium_pcf, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'pca' # 'harmony' 'pca'

# cluster ------
all_data_xenium_pcf <-
  FindNeighbors(
    all_data_xenium_pcf,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

all_data_xenium_pcf <-
  FindClusters(
    all_data_xenium_pcf,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data_xenium_pcf, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data_xenium_pcf, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data_xenium_pcf, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data_xenium_pcf, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data_xenium_pcf, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data_xenium_pcf, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data_xenium_pcf, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data_xenium_pcf, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()
  

# # cell marker ----
major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  Hepatocytes = 'AHSG', 
  Hepatocytes = 'APOM',
  
  Cholangiocyte = 'SOX9',
  Cholangiocyte = 'CFTR',
  
  "Hepatic stellate cell" = 'CCBE1',
  "Hepatic stellate cell" = 'GRK5',
  "Hepatic stellate cell" = 'LBP',
  
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
  # "Smooth muscle cells" = "PLN", 
  # "Smooth muscle cells" = "MYLK", 
  # "Smooth muscle cells" = "MYH11", 
  # "Smooth muscle cells" = "ACTA2",
  
  
  # "Myefibroblast" = "PDGFRB", 
  "Myefibroblast" = "PLN", 
  "Myefibroblast" = "MYLK", 
  "Myefibroblast" = "MYH11", 
  "Myefibroblast" = "TAGLN", 
  # "Myefibroblast" = "ACTA2",
  
  # "Breast glandular cells" = "SERPINA3", 
  # "Breast glandular cells" = "DSC2",
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  "Plasma" = "MZB1",
  "Plasma" = "XBP1",
  
  # Plasmablast = "MZB1",
  # Plasmablast = "XBP1",
  
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "CD4",
  "T cells" = "CD8A",
  
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
  "Macrophage" = "CSF1R",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  "Mast cell" = "MS4A2",
  "Mast cell" = "KIT",
  
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
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
  
  # Gial = "MPZ"
)


fibro_markers <- c(
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "ACTA2",
  
  "Fibroblasts" = 'COL1A2',
  "Fibroblasts" = 'COL3A1',
  "Fibroblasts" = 'LUM',
  "Fibroblasts" = 'FN1',
  "Fibroblasts" = 'COL6A1',
  "Fibroblasts" = 'FBLN1',
  
  "Fibroblasts" = 'VIM',
  "Fibroblasts" = 'THY',
  "Fibroblasts" = 'S100A4',
  "Fibroblasts" = 'PDPN',
  "Fibroblasts" = 'FAP'
)

fibro_markers %in% rownames(all_data)

pcf_marker <- rownames(all_data_xenium_pcf[['xenium_pcf']])  %>% str_subset("^PR-")
DefaultAssay(all_data_xenium_pcf) <- "xenium_pcf"

# dotplot
pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} major_marker2.pdf")), width = 25, height = 6)

  DotPlot(
    all_data_xenium_pcf, 
    group.by = "sketch_snn_res.0.8", 
    features = major_marker_ls
  ) + 
  RotatedAxis()
  
  DotPlot(
    all_data_xenium_pcf, 
    group.by = "sketch_snn_res.0.8", 
    features = pcf_marker
  ) + 
  RotatedAxis() + 
  scale_color_gradientn(colours = c("white", "white", "blue"))
  
dev.off()

# annotation -----
# dim20
# xenium_pcf_cell_cluster <- c(
#   '0' = "Cholangiocyte",
#   '1' = "Cholangiocyte proliferating",
#   '2' = "Monocyte_Macrophage",
#   '3' = "T cells",
#   '4' = "Fibroblasts",
#   '5' = "Hepatocytes",
#   '6' = "Others",
#   '7' = "Endothelial",
#   '8' = "Plasma",
#   '9' = "Hepatocytes",
#   '10'= "Hepatic stellate cells",
#   '11'= "Others",
#   '12'= "Epithelial",
#   '13'= "Monocyte",
#   '14'= "Macrophage",
#   '15'= "Pericytes",
#   '16'= "Others",
#   '17'= "B cells",
#   '18'= "Epithelial",
#   '19'= "Natural killer cells"
# )

xenium_pcf_cell_cluster <- c(
  '0' = "Epithelial",
  '1' = "Epithelial proliferating",
  '2' = "T cells",
  '3' = "Others", # "proliferating", # 20250914 update
  '4' = "Myeloid cells",
  '5' = "Others",
  '6' = "Fibroblasts",
  '7' = "Myeloid cells",
  '8' = "Endothelial",
  '9' = "Fibroblasts",
  '10'= "B/Plasma",
  '11'= "Hepatocytes",
  '12'= "Cholangiocyte",
  '13'= "Hepatocytes",
  '14'= "Endothelial",
  '15'= "Hepatic stellate cell",
  '16'= "Epithelial proliferating",
  '17'= "Epithelial"
)

Idents(all_data_xenium_pcf) <- "sketch_snn_res.0.8"

all_data_xenium_pcf <- RenameIdents(all_data_xenium_pcf, xenium_pcf_cell_cluster)
all_data_xenium_pcf$annotation <- Idents(all_data_xenium_pcf)

pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation {reduction.use} 2.pdf")), height = 6, width = 8)
DimPlot(all_data_xenium_pcf, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data_xenium_pcf, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation {reduction.use} - split2.pdf")), height = 16, width = 18)
DimPlot(all_data_xenium_pcf, group.by = "annotation", label = TRUE, reduction = 'umap') + facet_wrap(~ annotation) + NoLegend() 
DimPlot(all_data_xenium_pcf, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + facet_wrap(~ annotation)  + NoLegend()
dev.off()

# find markers ---------
DefaultAssay(all_data_xenium_pcf) <- "sketch"
Idents(all_data_xenium_pcf) <- "sketch_snn_res.0.8"

all_markers <-
  FindAllMarkers(
    JoinLayers(all_data_xenium_pcf),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")



# project data ------
DefaultAssay(all_data_xenium_pcf) <- "sketch"
all_data_xenium_pcf <- 
  all_data_xenium_pcf %>% 
  ProjectData(
    assay = "xenium_pcf",
    full.reduction = "pca.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "pca",
    
    umap.model = "umap", 
    dims = 1:dims, 
    
    refdata = list(
      annotation_full = "annotation",
      RNA_snn_res.0.2 = 'sketch_snn_res.0.2',
      RNA_snn_res.0.4 = 'sketch_snn_res.0.4',
      RNA_snn_res.0.6 = 'sketch_snn_res.0.6',
      RNA_snn_res.0.8 = 'sketch_snn_res.0.8'
    )
  )

Idents(all_data_xenium_pcf) <- 'annotation_full'
DefaultAssay(all_data_xenium_pcf) <- "xenium_pcf"

pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")), height = 6, width = 8)
DimPlot(all_data_xenium_pcf, group.by = "annotation_full", label = TRUE, reduction = 'full.umap') + NoLegend()
dev.off()


pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2 projected.pdf")), width = 25, height = 6)
DotPlot(all_data_xenium_pcf, group.by = "annotation_full", features = major_marker_ls) + RotatedAxis()
dev.off()


saveRDS(all_data_xenium_pcf, "all_data_xenium_pcf_dim10.rds")


all_data_xenium_pcf2 <- 
  all_data_xenium_pcf[ ,
    !all_data_xenium_pcf$annotation_full %in% c("Monocyte−Macrophage", "Others") &
      all_data_xenium_pcf$annotation_full.score > 0.9
  ]


# sub cluster -------------

out_dir <- "subcluster_xenium_pcf"
dir.create(out_dir)


## B/Plasma -----
# sc_mye
# ndims_t
# mye_markers
# mye/
# all_markers_mye
# mye_cell_cluster
# annotation_mye
# mye annotation


sc_b <- all_data_xenium_pcf[, all_data_xenium_pcf$annotation_full == "B/Plasma"]

DefaultAssay(sc_b) <- 'xenium_pcf'

sc_b <- sc_normalize(sc_b, scale.factor = "median")
sc_b <- sc_pca(
  sc_b, 
  sc_out_dir = file.path(out_dir, "b/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_b <- 10

sc_b <- sc_dr(
  sc_b, 
  ndims = ndims_b,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "b/sc.2 dr"),
  dir.create.force = TRUE
)


sc_b <- sc_cluster(
  sc_b, 
  ndims = ndims_b,
  reduction =  "pca", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "b/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

b_markers <- c(
  'EPCAM',
  "CD3D",
  "CD3E",
  "CCDC80",
  "PECAM1",
  "CD68",
  "CD14",
  "LILRA4",
  "CD1C",
  "IDO1",
  "XCR1",
  "CLEC9A",
  "CLEC10A",
  
  "CD19",
  "MS4A1",
  "CD79A",
  
  'CD27',
  
  
  "TNFRSF17",
  "JCHAIN",
  "MZB1",
  "XBP1"
)

sc_b <- sc_dotplot(
  sc_b, 
  ndims = ndims_b,
  marker = b_markers,
  sc_out_dir = file.path(out_dir, "b/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sc_b <- sc_dotplot(
  sc_b, 
  ndims = ndims_b,
  marker = pcf_marker,
  sc_out_dir = file.path(out_dir, "b/sc.4 dotplot - pcf"), 
  dir.create.force = TRUE
)
sc_b <- sc_dotplot(
  sc_b, 
  ndims = ndims_b,
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "b/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

pdf(
  file.path(out_dir, "b/sc.4 dotplot - pcf", 'dotplot-pcf.pdf'),
  height = 6,
  width  = 12
)
  DotPlot(
    sc_b,
    features = pcf_marker,
    group.by = "xenium_pcf_snn_res.0.2"
  ) +
  RotatedAxis() + 
  scale_color_gradientn(
    colours = c("white", "white", "blue")
  )
dev.off()


DotPlot(
  sc_b,
  features = c("CCDC80", "EPCAM", "CD3E", "CD68", "CD14", "MKI67", "TOP2A", "MS4A1", "CD27", "MZB1", "XBP1"),
  group.by = "xenium_pcf_snn_res.0.2"
)



all_markers_b <- 
  sc_findmarker(
    sc_b,
    idents = "xenium_pcf_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "group"
  )

sc_findmarker(sc_b, all_markers = all_markers_b, sort_by = "avg")
sc_findmarker(sc_b, all_markers = all_markers_b, sort_by = "pct")



dir.create(file.path(out_dir, "b/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "b/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sc_b,
  markers = all_markers_b,
  idents = "xenium_pcf_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sc_mark_heatmap(
  sc_b,
  markers = all_markers_b,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

b_cell_cluster <- c(
  '0'  =  "Plasma",
  '1'  =  "B cells",
  '2'  =  "Plasma",
  '3'  =  "Plasma",
  '4'  =  "CAF-like plasma cells", # "Plasma",
  '5'  =  "Plasma proliferating", # "Plasma",
  '6'  =  "pDCs"
)



Idents(sc_b) <- "xenium_pcf_snn_res.0.2"
sc_b <- RenameIdents(sc_b, b_cell_cluster)
sc_b$annotation_b <- Idents(sc_b)

dir.create(file.path(out_dir, "b/sc.5 annotation"))
pdf(file.path(out_dir, 'b/sc.5 annotation', str_glue("Dimplot dim{ndims_b} b annotation .pdf")), height = 6, width = 8)
DimPlot(sc_b, group.by = "annotation_b", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_b, group.by = "annotation_b", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_b, file.path(out_dir, "sc_b.rds"))


## myeloid -----
# sc_mye
# ndims_t
# mye_markers
# mye/
# all_markers_mye
# mye_cell_cluster
# annotation_mye
# mye annotation

sc_mye <- all_data_xenium_pcf[, all_data_xenium_pcf$annotation_full == "Myeloid cells"]

DefaultAssay(sc_mye) <- 'xenium_pcf'

sc_mye <- sc_normalize(sc_mye, scale.factor = "median")
sc_mye <- sc_pca(
  sc_mye, 
  sc_out_dir = file.path(out_dir, "mye/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_mye <- 10

sc_mye <- sc_dr(
  sc_mye, 
  ndims = ndims_mye,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "mye/sc.2 dr"),
  dir.create.force = TRUE
)


sc_mye <- sc_cluster(
  sc_mye, 
  ndims = ndims_mye,
  reduction =  "pca", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "mye/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

mye_markers <- c(
  'EPCAM',
  "CD3D",
  "CD3E",
  "CCDC80",
  "PECAM1",
  
  "MS4A1",
  "MZB1",
  "XBP1",
  'CD27',
  
  "CD68",
  "CD14",
  "FCGR3A",
  "C1QA",
  "C1QC",
  'CSF1R', # macro
  
  'FCN1', # Mono
  
  'XCR1', # cDC1
  'IDO1', # cDC2
  'LAMP3', # mature DC
  'LILRA4', #pDC
  'CLEC9A', # cDC1
  
  'MS4A2', # mast cell
  'KIT',
  
  'CSF3R', # netutrophil
  'S100A8', # mono?
  'S100A9' # mono?
)

sc_mye <- sc_dotplot(
  sc_mye, 
  ndims = ndims_mye,
  marker = mye_markers,
  sc_out_dir = file.path(out_dir, "mye/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sc_mye <- sc_dotplot(
  sc_mye, 
  ndims = ndims_mye,
  marker = pcf_marker,
  sc_out_dir = file.path(out_dir, "mye/sc.4 dotplot-pcf"), 
  dir.create.force = TRUE
)
sc_mye <- sc_dotplot(
  sc_mye, 
  ndims = ndims_mye,
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "mye/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)


pdf(
  file.path(out_dir, "mye/sc.4 dotplot-pcf", 'dotplot-pcf.pdf'),
  height = 6,
  width  = 12
)
DotPlot(
  sc_mye,
  features = pcf_marker,
  group.by = "xenium_pcf_snn_res.0.2"
) +
  RotatedAxis() + 
  scale_color_gradientn(
    colours = c("white", "white", "blue")
  ) +
  ggtitle("xenium_pcf_snn_res.0.2")
dev.off()

pdf(
  file.path(out_dir, "mye/sc.4 dotplot", 'dotplot-some markers.pdf'),
  height = 4,
  width  = 8
)
  DotPlot(
    sc_mye,
    features = c("CCDC80", "EPCAM", "CD3E", "CD68", "CD14", "MKI67", "TOP2A"),
    group.by = "xenium_pcf_snn_res.0.2"
  ) +
  ggtitle("xenium_pcf_snn_res.0.2")
dev.off()

all_markers_mye <- 
  sc_findmarker(
    sc_mye,
    idents = "xenium_pcf_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "group"
  )

sc_findmarker(sc_mye, all_markers = all_markers_mye, sort_by = "pct")



dir.create(file.path(out_dir, "mye/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "mye/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sc_mye,
  markers = all_markers_mye,
  idents = "xenium_pcf_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sc_mark_heatmap(
  sc_mye,
  markers = all_markers_mye,
  idents = "xenium_pcf_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

mye_cell_cluster <- c(
  '0'  =  "Monocyte",
  '1'  =  "Macrophage",
  '2'  =  "Monocyte",
  '3'  =  "Epithelial-like myeloid cells", # "Monocyte" # 20250914 update
  '4'  =  "CAF-like myeloid cells", # "Monocyte"  # 20250914 update
  '5'  =  "Neutrophil",
  '6'  =  "Mast cells"
)



Idents(sc_mye) <- "xenium_pcf_snn_res.0.2"
sc_mye <- RenameIdents(sc_mye, mye_cell_cluster)
sc_mye$annotation_mye <- Idents(sc_mye)

dir.create(file.path(out_dir, "mye/sc.5 annotation"))
pdf(file.path(out_dir, 'mye/sc.5 annotation', str_glue("Dimplot dim{ndims_mye} mye annotation .pdf")), height = 6, width = 8)
DimPlot(sc_mye, group.by = "annotation_mye", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_mye, group.by = "annotation_mye", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_mye, file.path(out_dir, "sc_mye.rds"))


## t -----
# sc_t
# ndims_t
# t_markers
# t/
# all_markers_t
# t_cell_cluster
# annotation_t
# t annotation

sc_t <- all_data_xenium_pcf[, all_data_xenium_pcf$annotation_full == "T cells"]

DefaultAssay(sc_t) <- 'xenium_pcf'

sc_t <- sc_normalize(sc_t, scale.factor = "median")
sc_t <- sc_pca(
  sc_t, 
  sc_out_dir = file.path(out_dir, "t/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_t <- 20

sc_t <- sc_dr(
  sc_t, 
  ndims = ndims_t,
  return.model = TRUE,
  sc_out_dir = file.path(out_dir, "t/sc.2 dr"),
  dir.create.force = TRUE
)


sc_t <- sc_cluster(
  sc_t, 
  ndims = ndims_t,
  reduction =  "pca", # "harmony", "pca"
  sc_out_dir = file.path(out_dir, "t/sc.3 cluster"),
  dir.create.force = TRUE,
  
  harmony.reduction.name = "harmonyumap",
  resolution = c(0.2, 0.4, 0.6, 0.8)
)

t_markers <- c(
  'EPCAM',
  
  "CD3D",
  "CD3E",
  "TRAC",
  "IL7R",
  "CD4",
  "CD8A",
  
  "IL2RA",
  "FOXP3",
  
  'TRDC',
  
  'EOMES',
  'PRF1',
  'NKG7',
  "NCAM1",
  "KLRF1",
  "KLRD1",
  
  'TCF7', 'SELL', 'CCR7', #Tn
  'GZMK', # Tem
  'CXCR6', # Trm
  'KLRG1', # TEMRA
  'PDCD1', 'TOX', 'TIGIT', 'CTLA4', # Tex
  'TOX2', 'CXCR5', # Tfh
  
  "FAP",
  "PDPN",
  "CCDC80",
  'PDGFRA',
  "PECAM1",
  
  "MS4A1",
  "MZB1",
  "XBP1",
  'CD27',
  
  "CD68",
  "CD14",
  "FCGR3A",
  "C1QA",
  "C1QC",
  'CSF1R', # macro
  
  'FCN1', # Mono
  
  'XCR1', # cDC1
  'IDO1', # cDC2
  'LAMP3', # mature DC
  'LILRA4', #pDC
  'CLEC9A', # cDC1
  
  'MS4A2', # mast cell
  'KIT',
  
  'CSF3R', # netutrophil
  'S100A8', # mono?
  'S100A9', # mono?
  
  "MKI67",
  "TOP2A"
)

sc_t <- sc_dotplot(
  sc_t, 
  ndims = ndims_t,
  marker = t_markers,
  sc_out_dir = file.path(out_dir, "t/sc.4 dotplot"), 
  dir.create.force = TRUE
)
sc_t <- sc_dotplot(
  sc_t, 
  ndims = ndims_t,
  marker = pcf_marker,
  sc_out_dir = file.path(out_dir, "t/sc.4 dotplot-pcf"), 
  dir.create.force = TRUE
)
sc_t <- sc_dotplot(
  sc_t, 
  ndims = ndims_t,
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "t/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

pdf(
  file.path(out_dir, "t/sc.4 dotplot-pcf", 'dotplot-pcf.pdf'),
  height = 6,
  width  = 12
)
  DotPlot(
    sc_t,
    features = pcf_marker,
    group.by = "xenium_pcf_snn_res.0.8"
  ) +
  RotatedAxis() + 
  scale_color_gradientn(
    colours = c("white", "white", "blue")
  )
dev.off()

all_markers_t <- 
  sc_findmarker(
    sc_t,
    idents = "xenium_pcf_snn_res.0.6",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "group"
  )

sc_findmarker(sc_t, all_markers = all_markers_t)
sc_findmarker(sc_t, all_markers = all_markers_t, sort_by = "pct")



dir.create(file.path(out_dir, "t/sc.6 heatmap"), recursive = TRUE)
pdf(
  file.path(out_dir, "t/sc.6 heatmap", 'heatmap.pdf'),
  height = 4,
  width  = 8
)
sc_mark_heatmap(
  sc_t,
  markers = all_markers_t,
  idents = "xenium_pcf_snn_res.0.8",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sc_mark_heatmap(
  sc_t,
  markers = all_markers_t,
  idents = "xenium_pcf_snn_res.0.8",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

# t_cell_cluster <- c(
#   '0'  =  "CAF-like T cells",
#   '1'  =  "CD8 Temra1",
#   '2'  =  "CD4 Tn",
#   '3'  =  "myeloid-like T cells",
#   '4'  =  "epithelial-like cells",
#   '5'  =  "CD8 Tn",
#   '6'  =  "T proliferating",
#   '7'  =  "CD8 Temra2",
#   '8'  =  "CD8 Tex",
#   '9'  =  "gamma delta T cells", # "γδ T cells",
#   '10' =  "CD8 T",
#   '11' =  "CD4 Treg",
#   '12' =  "Endothelial-like T cells",
#   '13' =  "B-like cells",
#   '14' =  "Plasma-like cells",
#   '15' =  "CD4 Tfh"
# )

# t_cell_cluster <- c(
#   '0'  = "CAF-like T cells",
#   '1'  = "myeloid-like T cells",
#   '2'  = "CD8 Temra1",
#   '3'  = "CD8 T",
#   '4'  = "epithelial-like cells",
#   '5'  = "CD4 Tn",
#   '6'  = "gamma delta T cells", # "γδ T cells",
#   '7'  = "CD8 Temra2",
#   '8'  = "CD8 T",
#   '9'  = "T proliferating",
#   '10' = "CD8 Tex",
#   '11' = "Endothelial-like T cells",
#   '12' = "CD4 Treg",
#   '13' = "Plasma-like cells",
#   '14' = "B-like cells",
#   '15' = "CD4 Tfh"
# )

# t_cell_cluster <- c(
#   '0'  = "CAF-like T cells",
#   '1'  = "CD8 Temra",
#   '2'  = "CD4 Tn",
#   '3'  = "Macrophage",
#   '4'  = "Endothelial",
#   '5'  = "CD8 Tex",
#   '6'  = "B-like cells",
#   '7'  = "Epithelial-like cells",
#   '8'  = "B-like cells",
#   '9'  = "CD8 Temra",
#   '10' = "Macrophage",
#   '11' = "T proliferating",
#   '12' = "CD4 Treg",
#   '13' = "CD8 Tem",
#   '14' = "Gamma delta T cells",
#   '15' = "Others"
# )

# 20250914 update
t_cell_cluster <- c(
  '0'  = "CAF-like T cells",
  '1'  = "CD8 Temra",
  '2'  = "CD4 Tn",
  '3'  = "Monocyte",
  '4'  = "Endothelial",
  '5'  = "CD8 Tex",
  '6'  = "Others",
  '7'  = "Epithelial-like T cells",
  '8'  = "B-like cells",
  '9'  = "CD8 Temra",
  '10' = "Macrophage",
  '11' = "T proliferating",
  '12' = "CD4 Treg",
  '13' = "CD8 Tem",
  '14' = "Gamma delta T cells",
  '15' = "CD8T"
)
Idents(sc_t) <- "xenium_pcf_snn_res.0.8"
sc_t <- RenameIdents(sc_t, t_cell_cluster)
sc_t$annotation_t <- Idents(sc_t)

sc_t$annotation_t[is.na(sc_t$annotation_t)] <- "Others"

dir.create(file.path(out_dir, "t/sc.5 annotation"))
pdf(file.path(out_dir, 't/sc.5 annotation', str_glue("Dimplot dim{ndims_t} t annotation .pdf")), height = 6, width = 8)
DimPlot(sc_t, group.by = "annotation_t", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t, group.by = "annotation_t", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


# saveRDS(sc_t, file.path(out_dir, "sc_t.rds"))
saveRDS(sc_t, file.path(out_dir, "sc_t.rds")) # 20250914 save



# combined all annotation -------------
annotation_df_subcluster <-
  list(
    sc_b[['annotation_b']],
    sc_mye[['annotation_mye']],
    sc_t[['annotation_t']]
  ) %>%
  map(function(dat){
    dat %>%
      rownames_to_column("cell_id") %>%
      dplyr::rename(annotation = 2)
  }) %>%
  bind_rows()

annotation_df_combined_subcluster <-
  all_data_xenium_pcf[['annotation_full']] %>%
  rownames_to_column("cell_id") %>%
  left_join(annotation_df_subcluster, by = "cell_id") %>%
  mutate(annotation_full2 = ifelse(
    is.na(annotation),
    as.character(annotation_full),
    as.character(annotation)
  )) %>%
  dplyr::select(cell_id, annotation_full2) %>% 
  column_to_rownames("cell_id")

all_data_xenium_pcf <- 
  all_data_xenium_pcf %>% 
  AddMetaData(annotation_df_combined_subcluster)

# all_data_xenium_pcf$annotation_full2[all_data_xenium_pcf$annotation_full2 == "proliferating"] <- "Others"
# all_data_xenium_pcf$annotation_full[all_data_xenium_pcf$annotation_full == "proliferating"] <- "Others"


saveRDS(all_data_xenium_pcf, 'all_data_xenium_pcf-v1-final.rds')


  
# spatial plot ----------------
cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
sp_p_ls_xenium_pcf <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    all_data_meta <-
      all_data_xenium_pcf2[[]] %>%
      rownames_to_column("cell_id")
    
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_full, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })

# export image
sp_dir <- "spatial plot_xenium_pcf" 
dir.create(sp_dir)

sp_p_ls_xenium_pcf %>% 
  iwalk(function(p, grp){ 
    pdf(
      file.path(sp_dir, str_glue("spatial plot {grp} - v2.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })


# spatial plot - sub cluster ----------------
cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
sp_p_ls_xenium_pcf2 <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    all_data_meta <-
      all_data_xenium_pcf[[]] %>%
      rownames_to_column("cell_id")
    
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_full2, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })

# export image
sp_dir2 <- "spatial plot_xenium_pcf2" 
dir.create(sp_dir2)

sp_p_ls_xenium_pcf2 %>% 
  iwalk(function(p, grp){ 
    pdf(
      file.path(sp_dir2, str_glue("spatial plot {grp} - v2.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })










