library(tidyverse)
library(Seurat)
library(sp)

# 91 16um 
visim_91_16 <-
  {
    path <- file.path("/home/ztt/spatial/visium_HD/FWSP20240131/838375-2-91_16" , "outs")
    Load10X_Spatial(path, bin.size = 16, slice = 'Handover_16')
  }

# microns_per_pixel
# 16: 数据是16um bins
mpp <- 16 / visim_91_16@images$Handover_16.016um@scale.factors$spot

# tumor boundary
# tumor_boundary仅用于重新圈选right area，method1中未使用
tumor_boundary <- read.csv("border_center.csv")
coords <- GetTissueCoordinates(visim_91_16)
tumor_boundary_coords <- coords[tumor_boundary$Barcode, ]


# method1 ---------------------
# LOUPE中根据tumor boundary重新圈选包含tumor boundary的右侧区域
# 读入这个右侧后，从总区域减去，则剩余部分为左侧区域
# 对tumor boundary进行100um的最近邻计算
# 则所有最近邻点组成了100um的浸润区域


# left and right area
right_area <- read.csv("838375-2-91_16_group.csv")
coords_right_area <- coords[coords$cell %in% right_area$Barcode, ]
coords_left_area  <- coords[! coords$cell %in% right_area$Barcode, ]

dim(coords_right_area)
dim(coords_left_area)

# plot and inspection
plot(coords_right_area$x, coords_right_area$y)
coords_right_area[chull(coords_right_area$x, coords_right_area$y), c('x', 'y')] %>% plot


# find tumor boundary
# 根据left area和right area找到绘制的tumor boundary
# 寻找left和right中的最近邻点，以20um（大概一个细胞，根据效果调整）为限度，
# 如果没有nn点则不是边界点。
# 找到left和right中的边界点，取均值即是最终的边界点

# dat: left_area
# query: right_area
nn_left_right <- 
  RANN::nn2(
    coords_left_area[1:2],
    coords_right_area[1:2],
    k = 10,
    searchtype = 'radius',
    radius = 20 / mpp
  )

non_boundary_idx <- apply(nn_left_right$nn.idx == 0, 1, all)
if(sum(non_boundary_idx) == nrow(coords_right_area)){
  print("radius is small, no boundary point.")
}

left_boundary_idx   <- nn_left_right$nn.idx[! non_boundary_idx, 1]
left_boundary       <- coords_left_area[left_boundary_idx, ]
right_boundary_idx  <- ! non_boundary_idx
right_boundary      <- coords_right_area[right_boundary_idx, ]

# comuted boudary coords
coords_df <-
  data.frame(
    x = (left_boundary$x + right_boundary$x) / 2,
    y = (left_boundary$y + right_boundary$y) / 2
  )

# validate tumor boundary
# 若边界未能找到，则调整radius的值
ggplot() +
  geom_point(
    # left area
    data = coords_left_area %>% dplyr::slice_sample(prop = 0.1), # 仅绘制10%的点，加快绘制速度
    aes(x = y, y = -x),
    color = "gray"
  ) +
  geom_point(
    # right area
    data = coords_right_area %>% dplyr::slice_sample(prop = 0.1), # 仅绘制10%的点，加快绘制速度
    aes(x = y, y = -x),
    color = "gray50"
  ) +
  
  geom_point(
    # comute left_boundary
    data = left_boundary,
    aes(x = y, y = -x),
    color = "red",
    size = 0.1
  ) + 
  geom_point(
    # compute right_boundary
    data = right_boundary,
    aes(x = y, y = -x),
    color = "blue",
    size = 0.1
  ) +
  geom_point(
    # the final tumor boundary
    data = coords_df,
    aes(x = y, y = -x),
    color = "green",
    size = 0.1
  ) +
  theme_void()

  



# right_area

# reduceFindNN find all (k) nn points in a certain radius,
# k would be increased when k is not enough for a certain radius
reduceFindNN <- function(dat, k = 10, query = NULL, searchtype = "radius", radius = 100, mpp = 1, max_k = 2000){
  if(is.null(query)){
    query <- dat
  }
  
  stopifnot(all(c('x', 'y') %in% colnames(dat)))
  stopifnot(all(c('x', 'y') %in% colnames(query)))
  
  dat   <- dat[c('x', 'y')]
  query <- query[c('x', 'y')]
  
  nn_dat <- RANN::nn2(
    dat, 
    query = query,
    k = k,
    radius = radius / mpp # radius, um -> pixel
  )
  
  k_enough <- sum(apply(nn_dat$nn.idx == 0, 1, any))/nrow(nn_dat$nn.idx)
  cat(str_glue("k is {k} and is enough for {k_enough * 100}% cells.\n\n"))
  
  while(k < max_k && k_enough < 0.99){
    
    if(k < 700){
      k <- k + 100
    }else{
      k <- k + 200
    }
    
    nn_dat <- RANN::nn2(
      dat, 
      query = query,
      searchtype = searchtype,
      k = k,
      radius = radius / mpp # radius, um -> pixel
    )
    
    k_enough <- sum(apply(nn_dat$nn.idx == 0, 1, any))/nrow(nn_dat$nn.idx)
    cat(str_glue("k is {k} and is enough for {k_enough * 100}% cells.\n\n"))
    
  }
 return(nn_dat) 
}

# 100um --
dat_nn <-
  reduceFindNN(
    dat   = coords_right_area, # %>% rbind(tumor_boundary_coords) %>% distinct(),
    query = coords_df, # tumor_boundary_coords,
    radius = 200,
    mpp = mpp,
    k = 100,
    max_k = 1000
  )

# all nn points is 100um infiltration points
points_idx_100um <- dat_nn$nn.idx %>% c() %>% .[. != 0] %>% unique

# 100um infiltration area
coords_100um     <- coords_right_area[points_idx_100um, ] 

# 200um  --
dat_nn <-
  reduceFindNN(
    dat   = coords_right_area, # %>% rbind(tumor_boundary_coords) %>% distinct(),
    query = coords_df, # tumor_boundary_coords,
    radius = 200,
    mpp = mpp,
    k = 100,
    max_k = 1000
  )

points_idx_200um <- dat_nn$nn.idx %>% c() %>% .[. != 0] %>% unique
coords_200um     <- coords_right_area[points_idx_200um, ] 

# validate 100um, 200um area
coords %>%
  dplyr::slice_sample(prop = 0.1) %>%
  ggplot(aes(x = y, y = -x)) +
  geom_point() +
  geom_point(
    data = coords_200um,
    color = "yellow",
    alpha = 0.5,
    size = 0.5
  ) +
  geom_point(
    data = coords_100um,
    color = "blue",
    alpha = 0.5,
    size = 0.5
  ) +
  geom_point(
    data = coords_df, # tumor_boundary_coords,
    color = "red",
    size = 0.5
  ) + 
  theme_void()


# method2 ---------------------
# coords平移
# 效果一般，平移方向不太好寻找


# 100 micron -> pixels
pixel100 <- 100 / mpp


# translation direction
tumor_bodundary_range <-
  tumor_boundary_coords %>% 
  .[c('x', 'y')] %>%
  map(range) %>%
  map_dbl(diff)

# y / x
theta <- atan(tumor_bodundary_range['y']/tumor_bodundary_range['x'])

# coords translation

coords_100um_v2 <-
  data.frame(
    x = tumor_boundary_coords$x + pixel100 * sin(theta),
    y = tumor_boundary_coords$y + pixel100 * cos(theta)
  )


# validate 100um boundary

coords %>%
  dplyr::slice_sample(prop = 0.1) %>%
  ggplot(aes(x = y, y = -x)) +
  geom_point() +
  geom_point(
    data = tumor_boundary_coords,
    color = "red"
  ) + 
  geom_point(
    data = coords_100um_v2,
    color = "blue"
  )
  theme_void()



# plot(coords$y, -coords$x)
# points(tumor_boundary_coords$y, -tumor_boundary_coords$x, col = "red")
# points(coords_100um$y, -coords_100um$x, col = "blue")
# 







