library(Banksy)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)


# HD示例数据：/mnt/NexData_file/home/ztt/test_data/Visium_HD/Colorectal_Cancer/output/data

out_dir <- "banksy"
dir.create(out_dir)


meta_data_anno <- read.csv('adata_celltype_obs.csv')

all_data <- readRDS("adata_filter.rds")
all_data[['RNA']] <- as(all_data[['RNA']], 'Assay5')

all_data$x <- Embeddings(all_data[['spatial']])[, 'SPATIAL_1'] 
all_data$y <- Embeddings(all_data[['spatial']])[, 'SPATIAL_2'] 

# # coords 平移
# rng_coord <-
#   all_data[[]] %>%
#   group_by(sample) %>%
#   mutate(
#     xmin = min(x),
#     xmax = max(x),
#     rng = xmax - xmin
#   )
# 
# all_data$xx <- all_data$x - rng_coord$xmin + as.numeric(factor(all_data$sample)) * 1.2 * max(rng_coord$rng)
# 
# data.frame(
#   all_data$xx,
#   all_data$y
# ) %>%
#   slice_sample(n = 10000) %>%
#   plot

# banksy ----------
all_data <- 
  all_data %>% 
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 1000) 


all_data <- 
  RunBanksy(all_data,
            lambda = 0.8, 
            verbose = TRUE,
            assay = "RNA", # "Spatial.008um", 
            slot = "data",
            features = "variable",
            k_geom = 50,
            
            dimx = "x", # use coord 
            dimy = "y",
            
            group = "sample"
  )

DefaultAssay(all_data) <- "BANKSY"

# data_own <- SeuratWrappers:::get_data(all_data, assay = "RNA", slot = 'data', features = 'variable', verbose = F)
# locs <- SeuratWrappers:::get_locs(all_data, dimx = "x", dimy = "y", dimz = NULL,
#                                   ndim = 2, data_own = data_own, group = 'sample', verbose = F)
# 
# 
# locs %>%
#   slice_sample(n = 10000) %>%
#   plot





# Run PCA and UMAP -----------
all_data <- RunPCA(
  all_data, 
  assay = 'BANKSY', 
  features = rownames(all_data[['BANKSY']]), 
  npcs = 30
)

# umap
ElbowPlot(all_data, ndims = 30)
all_data <- RunUMAP(all_data, dims = 1:20)


# harmony
all_data <- all_data %>%
  harmony::RunHarmony(
    dims.use = 1:20,
    group.by.vars = "sample",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data,
  dims = 1:20,
  return.model = TRUE,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)




# Clustering
all_data <- FindNeighbors(all_data, dims = 1:20, reduction = "harmony")
all_data <- FindClusters(all_data, resolution = c(0.2, 0.4, 0.6, 0.8))



DimPlot(
  all_data,
  group.by = "BANKSY_snn_res.0.2"
)


cluster_lt_1000 <- all_data$BANKSY_snn_res.0.2 %>% table %>% {names(.)[. < 1000]}

all_data2 <- all_data[, !all_data$BANKSY_snn_res.0.2 %in% cluster_lt_1000]



DimPlot(
  all_data2,
  group.by = "BANKSY_snn_res.0.2"
)

DimPlot(
  all_data2,
  group.by = "sample"
)


saveRDS(all_data, 'all_data_banksy.rds')

# plot  ---------------
all_data_banksy <- readRDS('all_data_banksy.rds')
# export pdf
pdf(
  file.path(out_dir, "umap plot of banksy.pdf"),
  height = 6,
  width = 8
)
DimPlot(
  all_data_banksy,
  group.by = "BANKSY_snn_res.0.2"
)
DimPlot(
  all_data_banksy,
  group.by = "sample"
)
dev.off()


# spatial plot
annotation_col <- c("BANKSY_snn_res.0.2")

cells_position_ls <-
  all_data_banksy[[c("x", "y", "sample")]] %>%
  rownames_to_column('cell_id') %>%
  dplyr::rename(
    SPATIAL_1 = x,
    SPATIAL_2 = y
  ) %>%
  group_by(sample) %>%
  nest() %>%
  deframe()

plot_dat <-
  FetchData(
    all_data_banksy, 
    vars = c(annotation_col)
  ) %>%
  rownames_to_column('cell_id')


p_sp_ct_ls <-
  cells_position_ls %>% # .[1] %>%
  imap(function(coords, sam){
    dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
    
    height <- 15
    width <- 15
    p <-
      dat2 %>%
      ggplot(aes(x = SPATIAL_1, SPATIAL_2)) +
      geom_tile(
        aes(fill = .data[[annotation_col]]),
        height = height, 
        width = width,
        show.legend = c(fill = TRUE, alpha = FALSE)
      ) +
      theme_void() +
      coord_equal() +
      scale_y_reverse() +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values = allcolour) +
      theme(panel.background = element_rect(fill = "black", color = NA)) +
      ggtitle(sam)
    p
    
  })


pdf(
  file.path(out_dir, "spatial plot of celltype in each sample.pdf"),
  height = 8,
  width = 8
)
p_sp_ct_ls %>% walk(print)
dev.off()


