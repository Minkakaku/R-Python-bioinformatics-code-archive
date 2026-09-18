library(Seurat)
library(tidyverse)
library(spacexr)


rctd_out_dir <- 'RCTD'
dir.create(rctd_out_dir)

# all_data <- readRDS("all_data.full.annotation.rds")

#1. sketch data  
# NOT RUN. ROI data, only 50k cells
# visium_roi_ls %>%
#   map(function(visium){
#     # sketch the cortical subset of the Visium HD dataset
#     DefaultAssay(visium) <- "Spatial.008um"
#     visium <- FindVariableFeatures(visium)
#     visium <- SketchData(
#       object = cortex,
#       ncells = 50000,
#       method = "LeverageScore",
#       sketched.assay = "sketch"
#     )
#     
#     DefaultAssay(visium) <- "sketch"
#     visium <- ScaleData(visium)
#     visium <- RunPCA(visium, assay = "sketch", reduction.name = "pca.sketch", verbose = T)
#     visium <- FindNeighbors(visium, reduction = "pca.sketch", dims = 1:50)
#     visium <- RunUMAP(visium, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50, verbose = T)
#     
#     visium
#   })

#2. single data

# all_data_sc <- readRDS('BP_AD_Pso_20240930.rds')
# all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
# all_data_sc$annotation <- all_data_sc$celltype 
# 
# Idents(all_data_sc) <- "annotation"
# counts_sc  <- all_data_sc %>% JoinLayers() %>% LayerData(layer = "counts")
# cluster_sc <- as.factor(all_data_sc$annotation)
# nUMI_sc    <- all_data_sc$nCount_RNA


counts_sc    <- arrow::read_parquet(file.path(rctd_out_dir, "sc_rtcd_counts.parquet")) %>% column_to_rownames("index")
sc_rtcd_meta <- arrow::read_parquet(file.path(rctd_out_dir, "sc_rtcd_meta.parquet")) %>% column_to_rownames("index")

# 25 cells 重复
counts_sc %>% colnames() %>% str_remove("-.*$") %>% duplicated() %>% sum()

# raw data 
# D3 D5
sc_path <- str_glue("singlecell/FWSC20240315_cellranger/{c('3D', '5D')}/filtered_feature_bc_matrix")
names(sc_path) <- c('3D', '5D')

sc_counts_ls   <- sc_path %>% map(Read10X)
sc_inter_cells <- sc_counts_ls %>% map(colnames) %>% purrr::reduce(intersect) 
sc_counts_ls   <- sc_counts_ls %>% map(~ .x[, ! colnames(.x) %in% sc_inter_cells])

sc_counts <- sc_counts_ls %>% map(as.matrix) %>% do.call(cbind, .)
stopifnot(all(!duplicated(colnames(sc_counts))))

# subset cells
subset_cells <- intersect(colnames(sc_counts), rownames(sc_rtcd_meta))

counts_sc_used    <- sc_counts[, subset_cells]
nUMI_sc           <- counts_sc_used %>% colSums()
cluster_sc        <- as.factor(sc_rtcd_meta[subset_cells, ]$cell_type_sub)
names(cluster_sc) <- subset_cells

cells_lt_25 <- (table(cluster_sc) < 25) %>% .[.] %>% names()

sc_cell_filter_idx <-  ! cluster_sc %in% c( "Unknown", cells_lt_25 )

# create the RCTD reference object
reference <- Reference(
  counts_sc_used[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
rctd_query_ls <-
  visim_ls %>%
  map(function(visium){# browser()
    counts_hd <- visium %>% LayerData(layer = "counts")
    coords    <- visium %>% GetTissueCoordinates() %>% .[, 1:2]
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })

#4. run RCTD
visim_ls[[1]] %>% .$nCount_Spatial.008um  %>% density  %>% plot
abline(v = 10)
table( visim_ls[[1]]$nCount_Spatial.008um > 10 )



visim_ls[[2]] %>% .$nCount_Spatial.008um  %>% density  %>% plot
abline(v = 10)
table( visim_ls[[2]]$nCount_Spatial.008um > 10 )



rctd_result_ls<- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(
      query, reference, 
      max_cores = 10,
      UMI_min   = 10
      # UMI_max
    )
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

# saveRDS(rctd_result_ls, 'rctd_result_ls.rds')

#5. addMetaData
visium_rctd_ls <-
  visim_ls %>%
  imap(function(visium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    visium <- AddMetaData(visium, metadata = rctd@results$results_df)
    
    visium$first_type <- as.character(visium$first_type)
    visium$first_type[is.na(visium$first_type)] <- "Other" # "Unknown"
    visium
  })

#6. spatial plot 
# SpatialFeaturePlot(
#   visium_roi_ls[[1]], 
#   image.scale = 'hires', 
#   pt.size.factor = 5,
#   shape = 22,
#   alpha = c(0.5, 0.8),
#   features = "nFeature_Spatial.008um"
# )

image_scale_ls <-
  visium_rctd_ls %>%
  # visium_roi_ls %>%
  map(GetTissueCoordinates) %>%
  map('x') %>%
  map(range) %>%
  map(diff) %>%
  map(log) %>%
  {map(., function(x, max_val)x/max_val, max_val = max(unlist(.)))}

pdf(
  file.path(rctd_out_dir, 'spatial plot of RCTD.pdf'), 
  height = 15, 
  width = 15
)
visium_rctd_ls %>%
  imap(function(visium, name){
    if(name == "HC") name_fix <- "NML" else name_fix <- name
    
    # visium$first_type[visium$first_type == "Unknown"] <- NA
    # visium$first_type[visium$first_type == "Other"] <- NA
    
    SpatialDimPlot(
      visium,
      group.by       = 'first_type',
      image.scale    = 'hires',
      pt.size.factor = 1 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()


# merge
visium_rctd <- merge(visium_rctd_ls[[1]], visium_rctd_ls[-1])
saveRDS(visium_rctd, file.path(rctd_out_dir, "visium_rctd.rds"))

# export csv
visium_rctd %>%
  GetTissueCoordinates() %>%
  write.csv(file = file.path(rctd_out_dir, "TissueCoordinates.csv"), row.names = FALSE)


visium_rctd[[]] %>%
  write.csv(file = file.path(rctd_out_dir, "rctd predicted table.csv"), row.names = FALSE)

