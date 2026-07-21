library(tidyverse)
library(Seurat)

out_dir_1110 <- "dat_analysis_251110"
dir.create(out_dir_1110)


all_data <- readRDS("~/temp_data/spatial_all_data.rds")

ratio_for_dim <- 
  all_data@images$DZZ_Saline.008um %>% 
  GetTissueCoordinates() %>% 
  .[1:2] %>% 
  map(range) %>%
  map(diff) %>% 
  {.$y / .$x}

# SpatialDimPlot(
#   all_data,
#   images = "DZZ_Saline.008um",
#   group.by = "annotation_full"
# ) +
#   theme_void()

#1. UMAP plot ----------

p_dim <- 
  DimPlot(
    all_data,
    group.by = "annotation_full",
    reduction = "full.harmonyumap",
    label = TRUE
  ) +
  coord_fixed(ratio = ratio_for_dim) +
  NoLegend() +
  ggsci::scale_color_igv()


# export
ggsave(
  p_dim,
  filename = file.path(out_dir_1110, "umap plot of all cells.pdf"),
  height = 6,
  width = 8
)

#2. feature plot ------------

p_feature <- 
  FeaturePlot(
    all_data,
    features = "Calb1",
    reduction = "full.harmonyumap",
    label = TRUE
  ) +
  coord_fixed(ratio = ratio_for_dim) 

# export
ggsave(
  p_feature,
  filename = file.path(out_dir_1110, "feature plot of Calb1.pdf"),
  height = 6,
  width = 8
)





