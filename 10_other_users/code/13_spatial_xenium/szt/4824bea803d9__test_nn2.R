library(Seurat)
library(tidyverse)


file_list <- 
  list.files(
    "/home/ztt/spatial/visium_HD/FWSP20240131",
    pattern = paste(sample_df$sample, collapse = "|"), 
    full.names = TRUE
  ) %>%
  str_subset("_[no|old|remove]", negate = TRUE)

names(file_list) <- sample_df$group # basename(file_list)


visim_ls_16 <-
  file_list %>% 
  imap(function(path, name){
    path <- file.path(path, "outs")
    Load10X_Spatial(path, bin.size = 16, slice = name)
  })


# meta.data
meta_181 <- readRDS("NN_mac_T/181_all_RCTD_results_df_241015.rds")

visium_181 <- visim_ls_16$Tumor
visium_181 <- visium_181 %>% AddMetaData(meta_181)

# coords
coords_df <- visium_181@images$Tumor.016um %>% GetTissueCoordinates()
coords_df_add_celltype <- 
  coords_df %>% left_join(
    rownames_to_column(visium_181[['first_type']], "cell"),
    by = "cell"
  )

# microns/pixels
# 16: 16um bins
scale_factor <- 16 / visium_181@images$Tumor.016um@scale.factors$spot 



# NN -----------

# 只是取最近邻细胞，k为10即可
coords_query <- coords_df_add_celltype %>% dplyr::filter(first_type == "B_cell") 

cn_nn <- RANN::nn2(
  coords_df_add_celltype[c('x', 'y')], 
  query = coords_query[c('x', 'y')],
  k = 10
)


# cn_nn$nn.idx是最近邻细胞的index，k为10，因此有10列，分别是距离由近到远的cell index
# 由于query cells也在all cells里面，因此query cells的最近邻一定是自身
# 所以要找真正的最近邻细胞，请使用第二列

query_idx <- which(coords_df_add_celltype$cell %in% coords_query$cell)
all(cn_nn$nn.idx[, 1] == query_idx) # query cells index和最近邻细胞index相同


# NN cells是第二列
nn_idx <- cn_nn$nn.idx[, 2]
nn_cells <- coords_df_add_celltype[nn_idx, ] # B cells周边其他细胞的cell_id等信息


nn_cells_dist <- cn_nn$nn.dists[, 2]
nn_cells_dist_um <- nn_cells_dist * scale_factor


cell_type <- as.character(coords_df_add_celltype$first_type[nn_idx])
nn_cells_dist_by_cells <- tapply(nn_cells_dist, cell_type, mean) # B cells周边其他细胞的平均距离 - 像素值
nn_cells_dist_um_by_cells <- tapply(nn_cells_dist_um, cell_type, mean) # B cells周边其他细胞的平均距离 - um


# barplot距离示意图
# 其他细胞的distance的距离比较接近，0 - 1 scale之后绘图
nn_cells_dist_um_by_cells %>%
  enframe("Cell", "distance") %>%
  mutate(distance_scale = (distance - min(distance))/diff(range(distance)) ) %>%
  ggplot(aes(x = Cell, y = distance_scale, fill = Cell)) +
  geom_col()









