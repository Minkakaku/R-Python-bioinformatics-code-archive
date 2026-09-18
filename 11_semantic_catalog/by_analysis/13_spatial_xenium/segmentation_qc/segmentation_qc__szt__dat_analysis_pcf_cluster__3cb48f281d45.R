library(tidyverse)
library(Seurat)

# multi workers  -----  
library(future)
plan("multisession", workers = 20, gc = TRUE)

options(future.globals.maxSize = 20 * 1024 ^ 3)

nbrOfFreeWorkers()


# set directory ---------
dr_outdir <- "pcf_reduction"
dir.create(dr_outdir)

# read data ---------
pcf_dir_ls <-
  list.files("pcf segment by xenium/", pattern = "alinged by xenium cell boundary", full.names = TRUE) %>% 
  str_subset("- 2")

names(pcf_dir_ls) <- 
  basename(pcf_dir_ls) %>%
  str_remove("\\s-.*$")

pcf_dat_ls <-
  pcf_dir_ls %>%
  map(function(path){ # browser()
    dat <- read_csv(path) %>% column_to_rownames("cell_id")
    dat <- dat[! colnames(dat) %in% "...1"]
    colnames(dat) <-
      colnames(dat) %>%
      str_split(":\\s") %>%
      map_chr(3)
    dat[is.na(dat)] <- 0
    dat
  })

pcf_sc_ls <-
  pcf_dat_ls %>%
  imap(function(dat, name){
    sc <-
      CreateSeuratObject(
        t(dat[, -1]),# 去DAPI，转置
        assay = "RNA", 
        project = name,
        min.cells = 1,
        min.features = 1
      )
    sc$project <- sc@project.name
    sc
  })


# rename cells 
pcf_sc_ls <-
  pcf_sc_ls %>% 
  imap(function(sc, name){
    RenameCells(sc, add.cell.id = name)
  })


# merge data -------
# all_data <- merge(visim_ls_pass[[1]], visim_ls_pass[-1])
all_data <- merge(pcf_sc_ls[[1]], pcf_sc_ls[-1])

# normalize for sketch data ------
# all_data <- NormalizeData(all_data, scale.factor = floor(mean(all_data$nCount_RNA)))
# all_data <- FindVariableFeatures(all_data)

counts_layer_name <- Layers(all_data, search = "counts")
# data_layer_name <-
#   counts_layer_name %>%
#   str_remove("^counts.") %>%
#   str_c("data", ., sep = ".")

for(idx in seq_along(counts_layer_name)){
  layer_counts <- counts_layer_name[idx]
  layer_data   <- counts_layer_name[idx] %>% str_remove("^counts.") %>% str_c("data", ., sep = ".")
  
  LayerData(all_data, layer = layer_data) <- as( log2( LayerData(all_data, layer = layer_counts) + 1), "dgCMatrix" )
}


VariableFeatures(all_data) <- rownames(all_data)


# # sketch data --------
# all_data <-
#   SketchData(
#     all_data,
#     ncells         = 50000L,
#     sketched.assay = "sketch",
#     method         = "LeverageScore",
#     var.name       = "leverage.score",
#     seed = 123L
#   )
# 
# sketch_cells_pass_qc <- all_data[['sketch']] %>% colnames() %>% .[. %in%  unlist(map(visim_ls_pass, colnames))]
# all_data[['sketch']]  <- all_data[['sketch']] %>% subset(cells = sketch_cells_pass_qc)

### dimension reduction
# DefaultAssay(all_data) <- "sketch"
# all_data <- FindVariableFeatures(all_data, nfeatures = 3000)

all_data <- ScaleData(
  all_data,
  vars.to.regress = c(
    'nCount_RNA'
  )
)




# all_data <- SCTransform(all_data, assay = "sketch")

all_data <- RunPCA(all_data, features = VariableFeatures(all_data))

pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()


set.seed(1234)
suset_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(dr_outdir, 'pca dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 1:5)
dev.off()


pdf(file.path(dr_outdir, 'pca dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 16:19)
dev.off()


# dimension reduction --------
dims <- 10

all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


all_data <- all_data %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data,
  dims = 1:dims,
  return.model = TRUE,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

# validate QC metric
FeaturePlot(
  all_data,
  features  = c('nCount_Spatial.008um', 'nFeature_Spatial.008um', 
                'percent.mt', 'log10GenesPerUMI'),
  reduction = "harmonyumap"
)


reduction.use <- 'harmony' # 'harmony' 'pca'
# cluster ------
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
    algorithm = 1, 
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()

dev.off()


# dotplot ---------
all_markers <- rownames(all_data)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} main_marker - default.pdf")), width = 20, height = 10)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", cols = c("white", "blue"), features = all_markers, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "RNA_snn_res.0.4", cols = c("white", "blue"), features = all_markers, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "RNA_snn_res.0.6", cols = c("white", "blue"), features = all_markers, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "RNA_snn_res.0.8", cols = c("white", "blue"), features = all_markers, dot.scale = 10) + theme_90angle
dev.off()

color_scale <- scale_color_gradientn(colours = c("white","white", "blue"))
pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} main_marker.pdf")), width = 20, height = 10)
DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = all_markers, dot.scale = 10) + theme_90angle + color_scale
DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = all_markers, dot.scale = 10) + theme_90angle + color_scale
DotPlot(all_data, group.by = "RNA_snn_res.0.6", features = all_markers, dot.scale = 10) + theme_90angle + color_scale
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = all_markers, dot.scale = 10) + theme_90angle + color_scale
dev.off()


pdf(file.path(dr_outdir, str_glue("feature plot dims{dims} {reduction.use} main_marker.pdf")), width = 8, height = 6)
FeaturePlot(
  all_data,
  reduction = "harmonyumap",
  features = all_markers,
  combine = FALSE
)
dev.off()


pdf(
  file.path(dr_outdir, str_glue("Vlnplot dims{dims} {reduction.use} main_marker.pdf")),
  width = 12,
  height = 20
)
VlnPlot(
  all_data, 
  features = rownames(all_data), 
  group.by = "RNA_snn_res.0.2", 
  ncol = 4,
  pt.size = 0
)
dev.off()

# annotation -----
# cell_cluster <- c(
#   '0' = "Fibroblast_macrophage cells",
#   '1' = "Others",
#   '2' = "Foxp3 cells",
#   '3' = "Others",
#   '4' = "CD8T",
#   '5' = "Monocyte",
#   '6' = "T_B cells",
#   '7' = "Endothelials",
#   '8' = "Others",
#   '9' = "proliferating cells",
#   '10'= "Fibroblasts",
#   '11'= "Macrophage",
#   '12'= "CD44 cells",
#   '13'= "HLA-I HLA-II cells",
#   '14'= "Multi"
# )
# 
# Idents(all_data) <- "RNA_snn_res.0.4"

cell_cluster <- c(
  '0' = "Tregs_C0",
  '1' = "Fibroblast_Macrophage",
  '2' = "Proliferating B cells",
  '3' = "Monocyte",
  '4' = "Endothelial",
  '5' = "Tregs_C5",
  '6' = "Fibroblast_C6",
  '7' = "Others",
  '8' = "CD4T_C8",
  '9' = "CD8T_C9",
  '10'= "Tregs_C10",
  '11'= "Fibroblast_Monocyte",
  '12'= "Proliferating cells",
  '13'= "HLA−A cells",
  '14'= "CD8T_C14",
  '15'= "Fibroblast_C15",
  '16'= "Macrophage",
  '17'= "Others",
  '18'= "B cells",
  '19'= "Others",
  '20'= "Others",
  '21'= "CD4T_C21",
  '22'= "Multi"
)

Idents(all_data) <- "RNA_snn_res.0.8"

all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation {reduction.use} 2.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation {reduction.use} - split2.pdf")), height = 16, width = 18)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + facet_wrap(~ annotation) + NoLegend() 
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + facet_wrap(~ annotation)  + NoLegend()
dev.off()

saveRDS(all_data, "all_data.pcf.rds")



# spatial plot ----------------
cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
sp_p_ls_pcf <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    all_data_meta <-
      all_data[[]] %>%
      rownames_to_column("cell_id")
    
    # sample1 -> 1
    # sample2 -> 2
    # sample3 -> 3
    all_data_meta$cell_id <- 
      all_data_meta$cell_id %>%
      str_remove("^sample")
    
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })

# export image
sp_dir <- "spatial plot_pcf" 
dir.create(sp_dir)

sp_p_ls_pcf %>% 
  iwalk(function(p, grp){ 
    pdf(
      file.path(sp_dir, str_glue("spatial plot {grp}.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })



