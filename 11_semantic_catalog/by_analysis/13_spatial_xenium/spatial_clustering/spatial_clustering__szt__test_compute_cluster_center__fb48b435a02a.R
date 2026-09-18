library(Seurat)
library(tidyverse)
library(patchwork)


list.files("compute_cluster_center/")

cc_df <- readRDS("compute_cluster_center/radius_barcode_coord_df.rds")

center_cc_df <- 
  cc_df %>%
  dplyr::filter(radius %in% c("radius_0um", "radius_25um") )

p1 <-
  cc_df %>%
  ggplot(aes(x = x, y = y, color = radius)) +
  geom_point()


dbscan::kNNdistplot(cc_df[1:2], minPts = 3)


cluster_dbscan <-
  center_cc_df[c('x', 'y')] %>%
  dbscan::dbscan(
    eps = 100, # <------------------ 测试不同的值，找到合适的分群
    minPts = 3, 
    # weights = NULL, 
    # borderPoints = TRUE
  )

cluster_df <- 
  data.frame(
    cell = center_cc_df$cell,
    cluster = as.character(cluster_dbscan$cluster)
  )


p2 <- cc_df %>%
  left_join(cluster_df, by = "cell") %>%
  
  ggplot(aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.01, show.legend = F)


# 检查是否是合适的分群
p1 | p2



# 计算边界
coords_center_df <-
  center_cc_df %>%
  left_join(cluster_df, by = "cell") %>%
  group_by(cluster) %>%
  summarise(
    x = mean(x),
    y = mean(y)
  )

coords_center_df %>%
  ggplot(aes(x = x, y = y)) +
  geom_point() 

# 检查中心点位置  
p2 + 
  # ggplot() +
  geom_point(
    data = coords_center_df,
    aes(x = x, y = y),
    shape = 16,
    size  = 2,
    color = 'black',
    show.legend = FALSE
  )


p1 + 
  # ggplot() +
  geom_point(
    data = coords_center_df,
    aes(x = x, y = y),
    shape = 16,
    size  = 2,
    color = 'black',
    show.legend = FALSE
  )
