#xenium数据-----
all_data <- readRDS("R_OUTPUT/all_data.rds")

library(tidyverse)
library(Seurat)
getwd()
library(patchwork)
options(java.parameters = '-Xmx10g')  # options('java.parameters')
library(rJava)
library(RBioFormats)
library(ggpubr)
library(SeuratDisk)
library(ggrepel)
library(loupeR)

# result_dir <- 'result' 
# dir.create(result_dir)


# RCTD ------------------
library(spacexr)
rctd_out_dir <- 'RCTD' 
dir.create(rctd_out_dir)


##2. singlecell data处理----
# all_data_rctd <- readRDS("all_data_sc_dim15_anno.rds")
# all_data_rctd <- all_data
sc_data_sub <- readRDS("/mnt/NexData_file/home/sh/Projects/AE2508274589-01/01.Mission/mapping_raw/Main_mapping.rds")
sc_data_sub@meta.data %>% head()
subtype <- sc_data_sub[["sub_cell_type"]]

sc_data <- readRDS("/mnt/NexData_file/home/sh/Projects/AE2508274589-01/01.Mission/Main_mapping.rds")
sc_data <- AddMetaData(sc_data,subtype["sub_cell_type"])

anno_counts <- table(sc_data@meta.data$sub_cell_type)
print(anno_counts)

# 从 SCT 提取原始计数矩阵
counts_matrix <- GetAssayData(sc_data, "RNA","counts")  

sc_data$celltype <- sc_data$sub_cell_type 
Idents(sc_data) <- "celltype"
cluster_sc <- as.factor(sc_data$celltype)
nUMI_sc    <- sc_data$total_counts

# create the RCTD reference object
reference <- Reference(
  counts_matrix, 
  droplevels(cluster_sc), 
  nUMI_sc
)


## 3. create the RCTD query object-----
all_data_sp_rctd_ls <- all_data %>% 
  JoinLayers() %>% 
  SplitObject('sample')

cells_position_ls <- readRDS("~/R_OUTPUT/data_for_analysis/cells_position_ls.rds")

rctd_query_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    counts_hd <- xenium %>% LayerData(layer = "counts")
    coords    <- 
      cells_position_ls[[name]] %>% 
      dplyr::select(
        cell_id,
        x = x_centroid, 
        y = y_centroid
      ) %>%
      column_to_rownames('cell_id')
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })
saveRDS(rctd_query_ls,"RCTD/rctd_query_ls.rds")


##4. run RCTD-----
rctd_result_ls <- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 30, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

saveRDS(rctd_result_ls, 'RCTD/rctd_result_ls.rds')

#rctd_result_ls <- readRDS("RCTD/rctd_result_ls.rds")

##5. addMetaData-----
all_data_sp_rctd_anno_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)
    
    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })

saveRDS(all_data_sp_rctd_anno_ls,'RCTD/rctd.rds')

all_data_sp_rctd_anno_ls$low_1_1@meta.data$first_type %>% unique()

#all_data_sp_rctd_anno_ls <- readRDS("rctd.rds")


#6. dimplot 
# all_data_sp_anno <-
#   merge(
#     all_data_sp_rctd_anno_ls[[1]],
#     all_data_sp_rctd_anno_ls[-1]
#   )

meta_rctd_result <-
  rctd_result_ls %>%
  map(~ .@results$results_df) %>%
  bind_rows()

all_data_rctd <-
  all_data %>% 
  AddMetaData(meta_rctd_result)

saveRDS(all_data_rctd,"RCTD/all_data_rctd.rds")

all_data_rctd <- JoinLayers(all_data_rctd)

all_data_rctd %>% DimPlot(reduction = "umap", group.by = "first_type", label = TRUE) + NoLegend()

all_data_rctd$first_type %>% table()

ttt <- all_data_rctd$first_type %>% is.na()



