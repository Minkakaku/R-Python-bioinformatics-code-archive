library(tidyverse)
library(Seurat)
library(patchwork)

source('sc_analysis_tools.R')

all_data <- readRDS("all_data_projected.rds")

all_data$annotation_full %>% table


out_dir <- "subcluster"
dir.create(out_dir)

# sub cluster -------------

# B/Plasma -----
# sc_mye
# ndims_t
# mye_markers
# mye/
# all_markers_mye
# mye_cell_cluster
# annotation_mye
# mye annotation


sc_b <- all_data[, all_data$annotation_full == "B/Plasma"]

DefaultAssay(sc_b) <- 'RNA'

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
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "b/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)

all_markers_b <- 
  sc_findmarker(
    sc_b,
    idents = "RNA_snn_res.0.2",
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE,
    # method = "conserved",
    grouping.var = "group"
  )

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
  idents = "RNA_snn_res.0.2",
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
  '4'  =  "Plasma",
  '5'  =  "Plasma",
  '6'  =  "Plasma",
  '7'  =  "pDCs"
)



Idents(sc_b) <- "RNA_snn_res.0.2"
sc_b <- RenameIdents(sc_b, b_cell_cluster)
sc_b$annotation_b <- Idents(sc_b)

dir.create(file.path(out_dir, "b/sc.5 annotation"))
pdf(file.path(out_dir, 'b/sc.5 annotation', str_glue("Dimplot dim{ndims_b} b annotation .pdf")), height = 6, width = 8)
DimPlot(sc_b, group.by = "annotation_b", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_b, group.by = "annotation_b", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_b, file.path(out_dir, "sc_b.rds"))


# myeloid -----
# sc_mye
# ndims_t
# mye_markers
# mye/
# all_markers_mye
# mye_cell_cluster
# annotation_mye
# mye annotation

sc_mye <- all_data[, all_data$annotation_full == "Myeloid cells"]

DefaultAssay(sc_mye) <- 'RNA'

sc_mye <- sc_normalize(sc_mye, scale.factor = "median")
sc_mye <- sc_pca(
  sc_mye, 
  sc_out_dir = file.path(out_dir, "mye/sc.1 pca"),
  dir.create.force = TRUE
)

ndims_mye <- 20

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
  marker = major_marker_for_subcluster,
  sc_out_dir = file.path(out_dir, "mye/sc.4 dotplot - major"), 
  dir.create.force = TRUE
)


all_markers_mye <- 
  sc_findmarker(
    sc_mye,
    idents = "RNA_snn_res.0.2",
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
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sc_mark_heatmap(
  sc_mye,
  markers = all_markers_mye,
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE,
  scale = "none",
  log_trans = TRUE
)

mye_cell_cluster <- c(
  '0'  =  "Macrophage",
  '1'  =  "Macrophage",
  '2'  =  "Epithelial-like myeloid cells",
  '3'  =  "Macrophage",
  '4'  =  "Plasma-like myeloid cells",
  '5'  =  "Neutrophil",
  '6'  =  "Macrophage", # "DC",
  '7'  =  "Mast cells"
)


Idents(sc_mye) <- "RNA_snn_res.0.2"
sc_mye <- RenameIdents(sc_mye, mye_cell_cluster)
sc_mye$annotation_mye <- Idents(sc_mye)

dir.create(file.path(out_dir, "mye/sc.5 annotation"))
pdf(file.path(out_dir, 'mye/sc.5 annotation', str_glue("Dimplot dim{ndims_mye} mye annotation .pdf")), height = 6, width = 8)
DimPlot(sc_mye, group.by = "annotation_mye", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_mye, group.by = "annotation_mye", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_mye, file.path(out_dir, "sc_mye.rds"))


# t -----
# sc_t
# ndims_t
# t_markers
# t/
# all_markers_t
# t_cell_cluster
# annotation_t
# t annotation

sc_t <- all_data[, all_data$annotation_full == "T cells"]

DefaultAssay(sc_t) <- 'RNA'

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
# sc_t <- sc_dotplot(
#   sc_t, 
#   ndims = ndims_t,
#   marker = major_marker_for_subcluster,
#   sc_out_dir = file.path(out_dir, "t/sc.4 dotplot - major"), 
#   dir.create.force = TRUE
# )

all_markers_t <- 
  sc_findmarker(
    sc_t,
    idents = "RNA_snn_res.0.8",
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
  idents = "RNA_snn_res.0.2",
  p_val_adj = 0.05,
  pct.1 = 0.15,
  n_head = 5,
  row_cluster = FALSE
)
dev.off()


sc_mark_heatmap(
  sc_t,
  markers = all_markers_t,
  idents = "RNA_snn_res.0.2",
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

t_cell_cluster <- c(
  '0'  =  "CAF-like T cells",
  '1'  =  "CD8 Temra",
  '2'  =  "CD4T",
  '3'  =  "myeloid-like T cells",
  '4'  =  "epithelial-like cells",
  '5'  =  "CD4T",
  '6'  =  "T proliferating",
  '7'  =  "CD8T",
  '8'  =  "CD8 Tex",
  '9'  =  "gamma delta T cells", # "γδ T cells",
  '10' =  "CD8T",
  '11' =  "CD4 Treg",
  '12' =  "Endothelial-like T cells",
  '13' =  "B-like cells",
  '14' =  "Plasma-like cells",
  '15' =  "CD4T"
)

Idents(sc_t) <- "RNA_snn_res.0.8"
sc_t <- RenameIdents(sc_t, t_cell_cluster)
sc_t$annotation_t <- Idents(sc_t)

dir.create(file.path(out_dir, "t/sc.5 annotation"))
pdf(file.path(out_dir, 't/sc.5 annotation', str_glue("Dimplot dim{ndims_t} t annotation .pdf")), height = 6, width = 8)
DimPlot(sc_t, group.by = "annotation_t", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(sc_t, group.by = "annotation_t", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(sc_t, file.path(out_dir, "sc_t.rds"))



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
  all_data[['annotation_full']] %>%
  rownames_to_column("cell_id") %>%
  left_join(annotation_df_subcluster, by = "cell_id") %>%
  mutate(annotation_full2 = ifelse(
    is.na(annotation),
    as.character(annotation_full),
    as.character(annotation)
  )) %>%
  dplyr::select(cell_id, annotation_full2) %>% 
  column_to_rownames("cell_id")

all_data <- 
  all_data %>% 
  AddMetaData(annotation_df_combined_subcluster)

# saveRDS(all_data, "all_data_projected_subcluster.rds")
saveRDS(all_data, "all_data_projected_subcluster_update250914.rds")

# spatial plot --------------


sp_p_ls2 <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]

    all_data_meta <-
      all_data[[]] %>%
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
sp_dir <- "spatial plot2" 
dir.create(sp_dir)

sp_p_ls2 %>% 
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



# export  annotation file












