library(Seurat)
library(tidyverse)
library(spacexr)


rctd_out_dir <- 'RCTD_20260627'
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

visim_ls <- readRDS("data_analysis_cycle2/visim_ls.rds")

#2. single data

# 要使用“RNA_snn_res.0.05”中的0，1，2，3四个细胞群
all_data_sc <- readRDS('data/hast_data.rds')

all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
all_data_sc$annotation <- all_data_sc$RNA_snn_res.0.05

Idents(all_data_sc) <- "annotation"
counts_sc  <- all_data_sc %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation)
nUMI_sc    <- all_data_sc$nCount_RNA


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
umi_min <- 30
visim_ls[[1]] %>% .$nCount_Spatial.008um  %>% density  %>% plot
abline(v = umi_min)
table( visim_ls[[1]]$nCount_Spatial.008um > umi_min )



visim_ls[[2]] %>% .$nCount_Spatial.008um  %>% density  %>% plot
abline(v = umi_min)
table( visim_ls[[2]]$nCount_Spatial.008um > umi_min )



rctd_result_ls<- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(
      query, reference, 
      max_cores = 10,
      UMI_min   = umi_min ### 
      # UMI_max
    )
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    saveRDS(RCTD, file.path(rctd_out_dir, str_glue("rctd_{name}.rds")))
    
    RCTD
  })

# saveRDS(rctd_result_ls, file.path(rctd_out_dir, 'rctd_result_ls.rds'))


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
all_celltype_rctd <- visium_rctd_ls$AZFb$first_type %>% unique() %>% str_sort()
color_celltype <- structure(
  names = all_celltype_rctd,
  ggsci::pal_igv()(length(all_celltype_rctd))
)

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
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values = color_celltype) +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()

# split sample 
visium_rctd_ls %>%
  imap(function(visium, name){
    if(name == "HC") name_fix <- "NML" else name_fix <- name
    
    pdf(
      file.path(rctd_out_dir, str_glue('spatial plot of RCTD - {name}.pdf')), 
      height = 15, 
      width = 15
    )
      p <-
        SpatialDimPlot(
          visium,
          group.by       = 'first_type',
          image.scale    = 'hires',
          pt.size.factor = 3 * image_scale_ls[[name]],
          shape          = 22,
          # alpha          = 0.5
        ) +
        ggtitle(name_fix) +
        # ggsci::scale_fill_igv() +
        scale_fill_manual(values = color_celltype) +
        guides(fill = guide_legend(override.aes = list(size = 5)))
      print(p)
    dev.off()
  })


# split sample and celltype
visium_rctd_ls %>%
  imap(function(visium, name){
    if(name == "HC") name_fix <- "NML" else name_fix <- name
    
    pdf(
      file.path(rctd_out_dir, str_glue('spatial plot of RCTD - {name} - split celltype.pdf')), 
      height = 15, 
      width = 15
    )
      p <- 
        visium$first_type %>%
        unique() %>%
        map(function(x){
          color_celltype[names(color_celltype) != x] <-  NA
          color_celltype[names(color_celltype) == x] <- 'blue'
          
          SpatialDimPlot(
            visium,
            group.by       = 'first_type',
            image.scale    = 'hires',
            pt.size.factor = 3 * image_scale_ls[[name]],
            shape          = 22,
            # alpha          = 0.5
          ) +
          ggtitle(name_fix) +
          # ggsci::scale_fill_igv() +
          scale_fill_manual(values = color_celltype, na.value = NA) +
          guides(fill = guide_legend(override.aes = list(size = 5)))
        })
      print(p)
    dev.off()
  })



# # merge
# visium_rctd <- merge(visium_rctd_ls[[1]], visium_rctd_ls[-1])
# saveRDS(visium_rctd, file.path(rctd_out_dir, "visium_rctd.rds"))
# 
# # export csv
# visium_rctd %>%
#   GetTissueCoordinates() %>%
#   write.csv(file = file.path(rctd_out_dir, "TissueCoordinates.csv"), row.names = FALSE)
# 
# 
# visium_rctd[[]] %>%
#   write.csv(file = file.path(rctd_out_dir, "rctd predicted table.csv"), row.names = FALSE)

