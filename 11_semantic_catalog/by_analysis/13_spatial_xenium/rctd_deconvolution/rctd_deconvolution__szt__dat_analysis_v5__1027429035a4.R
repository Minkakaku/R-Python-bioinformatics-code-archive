library(tidyverse)
library(Seurat)


# out_dir <- "output" 
# dir.create(out_dir)

# all_data <- readRDS("all_data_projected_v2.rds")
all_data <- readRDS("all_data_projected_v2_hi_res.rds")

# 1 read data -------
sample_df <- 
  # read.csv("data/Annoation.csv")
  data.frame(
    sample = c('2310163', '2305735', '2304480', '2304434'),
    group  = c('HC', 'BP1', 'BP2', 'BP3')
  )

file_list <- 
  list.files(
    # "/home/ztt/spatial/visium_HD/Human_Skin_FFPE/spaceranger_output",
    "/mnt/NexData_file/home/ztt/backup_203/spatial/visium_HD/Human_Skin_FFPE/spaceranger_output/",
    pattern = paste(sample_df$sample, collapse = "|"), 
    full.names = TRUE
  ) %>%
  str_subset("\\..*$", negate = TRUE)

stopifnot(nrow(file_list) == nrow(sample_df))

# TODO: fix bugs
# Samples order in file_list may not identical to sample_df
# names(file_list) <- sample_df$group # basename(file_list)
names(file_list) <- sample_df$group[match(basename(file_list), sample_df$sample)]




visim_ls <-
  file_list %>% 
  imap(function(path, name){
    path <- file.path(path, "outs")
    
    Read10X_Image <-
      function (image.dir, image.name = "tissue_lowres_image.png", 
                assay = "Spatial", slice = "slice1", filter.matrix = TRUE) {
        image <- png::readPNG(source = file.path(image.dir, image.name))
        scale.factors <- Read10X_ScaleFactors(filename = file.path(image.dir, 
                                                                   "scalefactors_json.json"))
        coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(image.dir, 
                                                                         "*tissue_positions*parquet")), filter.matrix)
        fov <- CreateFOV(coordinates[, c("imagerow", "imagecol")], 
                         type = "centroids", radius = scale.factors[["spot"]], 
                         assay = assay, key = Key(slice, quiet = TRUE))
        visium.fov <- new(Class = "VisiumV2", boundaries = fov@boundaries, 
                          molecules = fov@molecules, assay = fov@assay, key = fov@key, 
                          image = image, scale.factors = scale.factors)
        return(visium.fov)
      }
    
    Load10X_Spatial(
      path, 
      bin.size = 8, 
      slice = name, 
      image = Read10X_Image(file.path(path, "binned_outputs/square_008um/spatial"), assay = 'Spatial.008um')
    )
  })

# add sample_id or sample_name
visim_ls <-
  visim_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })


#2  RCTD -----------------------
library(spacexr)
rctd_out_dir <- 'RCTD - public singlecell' 
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
# all_data_sc <- readRDS('BP_AD_Pso-v2.rds')
# all_data_sc <- readRDS('BP_AD_Pso_20240930.rds')
# all_data_sc <- readRDS('BP_AD_Pso_20241028.rds')
load('/mnt/NexData_file/home/szt/lr/20250512_FWSP20241358_skin/Filter_Basicinfo.rData')
all_data_sc <- pbmc_filt

all_data_sc <- all_data_sc[, !all_data_sc$celltype %in% c("Malignant cells", "other")]

all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
all_data_sc$annotation <- all_data_sc$celltype 

Idents(all_data_sc) <- "annotation"
counts_sc  <- all_data_sc %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation)
nUMI_sc    <-  all_data_sc$nCount_RNA

sc_cell_filter_idx <-  cluster_sc != "Unknown"

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
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
rctd_result_ls<- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(
      query, reference, max_cores = 10,
      UMI_min = 50
      # UMI_max
    )
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

saveRDS(rctd_result_ls, file.path(rctd_out_dir, 'rctd_result_ls.rds'))

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


all_data <-
  visium_rctd_ls %>%
  imap(function(visium, name){
    visium$project <- name
    visium <- visium %>% RenameCells(add.cell.id = name)
    visium
  }) %>%
  {merge(.[[1]], .[-1])}

saveRDS(all_data, 'all_data_rctd - v5.rds')


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
    visium$first_type[visium$first_type == "Other"] <- NA
    
    SpatialDimPlot(
      visium,
      group.by       = 'first_type',
      image.scale    = 'lowres',
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()
