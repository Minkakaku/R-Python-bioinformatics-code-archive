visium <- Load10X_Spatial("/home/ztt/spatial/visium_HD/FWSP20240131/838375-1-181/outs" )

library(tidyverse)
library(Seurat)

visium <- Load10X_Spatial("path_to_visium_data" )


coord_8um  <- visium@images$slice1.008um %>% GetTissueCoordinates()
coord_16um <- visium@images$slice1.016um %>% GetTissueCoordinates()

nn <-
  RANN::nn2(
    data  = coord_16um[1:2],
    query = coord_8um[1:2]
  )

coords_match_table <-
  data.frame(
    id_16um = coord_16um$cell[nn$nn.idx[,1]],
    id_8um  = coord_8um$cell,
    distance = nn$nn.dists[, 1]
  )


coords_match_table %>% group_by(id_16um) %>% summarise(n()) %>% .[[2]] %>% table
# 1      2      3      4 
# 137    243    132 173186
