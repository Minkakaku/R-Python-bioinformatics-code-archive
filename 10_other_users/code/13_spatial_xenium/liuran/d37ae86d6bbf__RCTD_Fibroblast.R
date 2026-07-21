library(tidyverse)
library(Seurat)
library(future)

plan("multisession", workers = 16, gc = TRUE)
options(future.globals.maxSize = 50 * 1024 ^ 3)

all_data <- readRDS("/mnt/NexData_file/home/liuran/projects/3_Xenium/XE2512247768-08_human_skin/data_for_analysis/all_data_rctd_260716.rds")
all_data_sp_rctd <- all_data[, all_data@meta.data$RCTD == "Fibroblast"]

cells_position_ls <- readRDS("/mnt/NexData_file/home/liuran/projects/3_Xenium/XE2512247768-08_human_skin/data_for_analysis/cells_position_ls.rds")
cells_boundary_ls <- readRDS("/mnt/NexData_file/home/liuran/projects/3_Xenium/XE2512247768-08_human_skin/data_for_analysis/cells_boundary_ls.rds")

cells_boundary_ls_raw <- readRDS("/mnt/NexData_file/home/liuran/projects/3_Xenium/XE2512247768-08_human_skin/data_for_analysis/cells_boundary_ls_raw.rds")
cells_position_ls_raw <- readRDS("/mnt/NexData_file/home/liuran/projects/3_Xenium/XE2512247768-08_human_skin/data_for_analysis/cells_position_ls_raw.rds")

### RCTD ------------------
library(spacexr)
rctd_out_dir <- "/mnt/NexData_file/home/liuran/projects/3_Xenium/XE2512247768-08_human_skin/260717-RCTD_Fibroblast"
dir.create(rctd_out_dir)


### 2.single data ------
all_data_sc <- readRDS("/mnt/NexData_file/home/liuran/projects/3_Xenium/XE2512247768-08_human_skin/scRNA annotation/Ps_Fibroblast_annotation/all_data_sc.rds")
colnames(all_data_sc@meta.data)

all_data_sc$annotation <- all_data_sc$celltype

Idents(all_data_sc) <- "annotation"
counts_sc <- all_data_sc %>%
  JoinLayers() %>%
  LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation)
nUMI_sc <- all_data_sc$nCount_RNA

sc_cell_filter_idx <- !cluster_sc == "other"

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx],
  droplevels(cluster_sc[sc_cell_filter_idx]),
  nUMI_sc[sc_cell_filter_idx]
)

### 3.create the RCTD query object ------
all_data_sp_rctd@meta.data$sample %>% unique()
all_data_sp_rctd@meta.data$group %>% unique()

all_data_sp_rctd_ls <- all_data_sp_rctd %>%
  JoinLayers() %>%
  SplitObject("sample")

rctd_query_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name) { # browser()
    counts <- xenium %>% LayerData(layer = "counts")
    coords <-
      cells_position_ls[[name]] %>%
      dplyr::select(
        cell_id,
        x = x_centroid,
        y = y_centroid
      ) %>%
      column_to_rownames("cell_id")
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts, colSums(counts))
    query
  })

### 4.run RCTD ------
rctd_result_ls <-
  rctd_query_ls %>%
  imap(function(query, name) {
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 10, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

saveRDS(rctd_result_ls, file.path(rctd_out_dir, "rctd_result_ls.rds"))


### 5.addMetaData ------
all_data_sp_rctd_anno_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name) { # browser()
    rctd <- rctd_result_ls[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)
    
    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })


### 6.dimplot ------
meta_rctd_result <-
  rctd_result_ls %>%
  map(~ .@results$results_df) %>%
  bind_rows()

all_data_sp_rctd_anno <-
  all_data_sp_rctd %>%
  AddMetaData(meta_rctd_result)

all_data_sp_rctd_anno$first_type <- as.character(all_data_sp_rctd_anno$first_type)
all_data_sp_rctd_anno$first_type[is.na(all_data_sp_rctd_anno$first_type)] <- "Unknown"

saveRDS(all_data_sp_rctd_anno, file.path(rctd_out_dir, "all_data_sp_rctd_anno.rds"))
