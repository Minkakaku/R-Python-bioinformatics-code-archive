library(tidyverse)
library(Seurat)
library(RANN)


visium <- visim_ls$Tumor


# radius 50um
radius <- 50

# coords
coords_df <- visium@images$Tumor.008um %>% GetTissueCoordinates()

# microns/pixels
# 8: 8um bins
scale_factor <- 8 / visium@images$Tumor.008um@scale.factors$spot 


# max nn point to find
max_k <- ceiling( (50 / (8/2) ) ^ 2 ) + 10

# compute nn cells
res <-
  RANN::nn2(
    coords_df[c("x", "y")], # all cells
    query = coords_df[c("x", "y")] %>% dplyr::slice_sample(n = 100000), # query cells, 仅选择了100,000个细胞作为测试
    k = max_k,
    searchtype = 'radius',
    radius = 50/scale_factor # 50um -> pixels
  )

k_enough <- sum(apply(res$nn.idx == 0, 1, any))/nrow(res$nn.idx )
if(k_enough < 1) stop(str_glue("k is {max_k} and should be larger, then re-run RANN::nn2"))


# change  index to spot name
# res$nn.idx is cell index of NN
# when nn.idx is 0, it's not a cell

cells <- c("NA", coords_df$cell)
nn.idx <- res$nn.idx + 1

nn.name <- cells[nn.idx]
nn.name <- matrix(
  nn.name, 
  nrow = nrow(nn.idx)
)



# nn.name就是每个query细胞50um内的最近邻细胞ID
# row: query cells
# col: all nn cells






