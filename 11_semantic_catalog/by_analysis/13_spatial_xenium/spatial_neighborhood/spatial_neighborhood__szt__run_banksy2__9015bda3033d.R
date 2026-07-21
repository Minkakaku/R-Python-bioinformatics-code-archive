library(Banksy)
library(Seurat)
library(SeuratWrappers)
library(tidyverse)


# HD示例数据：/mnt/NexData_file/home/ztt/test_data/Visium_HD/Colorectal_Cancer/output/data


meta_data_anno <- read.csv('adata_celltype_obs.csv')

all_data <- readRDS("adata_filter.rds")
all_data[['RNA']] <- as(all_data[['RNA']], 'Assay5')
all_data[['RNA']] <- split(all_data[['RNA']], f = all_data$sample)

all_data$x <- Embeddings(all_data[['spatial']])[, 'SPATIAL_1'] 
all_data$y <- Embeddings(all_data[['spatial']])[, 'SPATIAL_2'] 


all_data <-
  all_data %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000)


## sketch data --------
# undebug(SketchData)
all_data <-
  SketchData(
    all_data,
    ncells         = 50000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    seed = 123L
  )

# sketch_cells_pass_qc <- all_data[['sketch']] %>% colnames() %>% .[. %in%  unlist(map(visim_ls_pass, colnames))]
# all_data[['sketch']]  <- all_data[['sketch']] %>% subset(cells = sketch_cells_pass_qc)

### dimension reduction
DefaultAssay(all_data) <- "sketch"
all_data <- FindVariableFeatures(all_data, nfeatures = 3000)
all_data <- ScaleData(
  all_data,
  vars.to.regress = c(
    'nCount_RNA'
  )
)

all_data_sketch <- CreateSeuratObject(all_data[['sketch']], meta.data = all_data[[]])
# all_data <- JoinLayers(all_data)

all_data_sketch <- 
  RunBanksy(all_data_sketch,
            lambda = 0.8, 
            verbose = TRUE,
            assay = "RNA", # "Spatial.008um", 
            slot = "data",
            features = "variable",
            k_geom = 50,
            
            dimx = "x",
            dimy = "y",
            
            group = "sample"
  )

# Run PCA and UMAP
all_data[['BANKSY']] <- all_data_sketch[['BANKSY']]
all_data <- 
  all_data %>%
  AddMetaData(all_data_sketch[[c('staggered_sdimx', 'staggered_sdimy')]])

DefaultAssay(all_data) <- "BANKSY"

all_data <- RunPCA(
  all_data_sketch, 
  assay = 'BANKSY', 
  features = rownames(all_data[['BANKSY']]), 
  npcs = 30
)

# umap
ElbowPlot(all_data, ndims = 30)
all_data <- RunUMAP(all_data, dims = 1:10)


# harmony
all_data <- all_data %>%
  harmony::RunHarmony(
    dims.use = 1:10,
    group.by.vars = "sample",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data,
  dims = 1:10,
  return.model = TRUE,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)




# Clustering
all_data <- FindNeighbors(all_data, dims = 1:10, reduction = "harmony")
all_data <- FindClusters(all_data, resolution = c(0.2, 0.4, 0.6, 0.8))



DimPlot(
  all_data,
  group.by = "BANKSY_snn_res.0.2"
)

DimPlot(
  all_data,
  group.by = "sample"
)



# all_data <- 
#   all_data %>% 
#   ProjectData(
#     assay = "RNA",
#     full.reduction = "pca.full",
#     
#     sketched.assay = "BANKSY",
#     sketched.reduction = "harmony",
#     
#     umap.model = "harmonyumap", 
#     dims = 1:10, 
#     
#     refdata = list(
#       # annotation_full = "annotation",
#       # annotation2_full = "annotation2",
#       BANKSY_snn_res.0.2_full = "BANKSY_snn_res.0.2",
#       BANKSY_snn_res.0.4_full = "BANKSY_snn_res.0.4",
#       BANKSY_snn_res.0.6_full = "BANKSY_snn_res.0.6",
#       BANKSY_snn_res.0.8_full = "BANKSY_snn_res.0.8"
#     )
#   )

# DefaultAssay(all_data) <- "RNA"

# 
# 
# # annotation
# # all_data$annotation2_full2 <- all_data$annotation2_full
# all_data$BANKSY_snn_res.0.2_full[all_data$BANKSY_snn_res.0.2_full < 0.9] <- NA
# 


# 
# cluster_lt_1000 <- all_data$BANKSY_snn_res.0.2 %>% table %>% {names(.)[. < 1000]}
# 
# all_data2 <- all_data[, !all_data$BANKSY_snn_res.0.2 %in% cluster_lt_1000]
# 
# 
# 
# DimPlot(
#   all_data2,
#   group.by = "BANKSY_snn_res.0.2"
# )
# 
# DimPlot(
#   all_data2,
#   group.by = "sample"
# )


saveRDS(all_data, 'all_data_banksy_use_sketch.rds')







