

all_data <- readRDS("~/R_OUTPUT/all_data_4slides_sample_filtered_reUMAP.rds")
# qs2::qs_save(all_data, 'spatialEcoTyper/all_data_4slides_sample_filtered_reUMAP.qs2')

all_data2 <- all_data[, all_data$annotation != 'Skeletal muscle cell']
all_data2$annotation %>% table()
# qs2::qs_save(all_data, 'spatialEcoTyper/all_data_4slides_sample_filtered_reUMAP_2.qs2')

cells_position_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_position_ls_4slides_filtered.rds')

coord_mat <- 
  cells_position_ls %>% 
  bind_rows() %>% 
  filter(cell_id %in% Cells(all_data2)) %>% 
  select(cell_id, x_centroid, y_centroid) %>% 
  rename(spatial_1 = x_centroid, spatial_2 = y_centroid) %>% 
  column_to_rownames('cell_id') %>% 
  as.matrix()

all_data2[["spatial"]] <- CreateDimReducObject(
  embeddings = coord_mat,
  assay = "RNA",        # 与你的空间 assay 名称对应
  key = "spatial_"
)

all_data2@reductions$spatial

# /mnt/NexData_file/home/hhq/envs/r4.2
# qs::qsave(all_data2, 'spatialEcoTyper/all_data_4slides_sample_filtered_reUMAP_2.qs')
