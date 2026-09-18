# tumor boundary: polygon to a line
nn_tumor_boundary <-
  RANN::nn2(
    tumor_boundary_coords[1:2],
    tumor_boundary_coords[1:2],
    k = 10
  )
compute_jaccard_index <- function(x, y){
  length(intersect(x, y))/length(union(x, y))
}

jaccard_idx_vector <- c()
for(i in seq.int(nrow(nn_tumor_boundary$nn.idx))){
  inter_len <- apply(nn_tumor_boundary$nn.idx, 1, function(x) length(intersect(x, nn_tumor_boundary$nn.idx[i, ])))
  union_len <- apply(nn_tumor_boundary$nn.idx, 1, function(x) length(union(x, nn_tumor_boundary$nn.idx[i, ])))
  
  jaccard_idx_vector <- c(jaccard_idx_vector, inter_len / union_len)
}

jaccard_idx_mat <- 
  matrix(
    jaccard_idx_vector,
    nrow = nrow(nn_tumor_boundary$nn.idx)
  )


tumor_boundary_coords_tidy <- 
  jaccard_idx_mat %>%
  apply(1, function(x) which(x > 0.4), simplify = FALSE) %>%
  sapply(function(x) c(
    x = mean(tumor_boundary_coords[x, 'x']),
    y = mean(tumor_boundary_coords[x, 'y'])
  )) %>%
  t

xspline(
  tumor_boundary_coords_tidy[, 'x'],
  tumor_boundary_coords_tidy[, 'y'],
  shape = 1,
  open = FALSE,
  draw = FALSE
)

plot(tumor_boundary_coords_tidy, cex = 0.2)



# find_last_nn_points
# the last nn point would be the boundary point.
find_last_nn_points <- function(x){
  stopifnot(is.matrix(x))
  
  res <- apply(x, 1, function(xx){
    idx <- xx == 0
    
    if(any(idx)){
      idx2 <- which(idx)
      idx2 <- idx2[1] - 1
      if(idx2 == 0) return(0) # have no nn points
      return(xx[idx2])
    }else{
      return(xx[length(xx)])
    }
  })
  
  res
}
# debug(reduceFindNN)


# left and right area
he_bbox <- bbox(as.matrix(coords[c("x", "y")]))

tumor_boundary_endpoint <-
  tumor_boundary_coords[c('x', 'y')] %>%
  map_dfc(range) %>%
  t %>%
  {colnames(.) <- c("min", "max"); .} %>%
  as.matrix()

tumor_boundary_centerpoint <-
  tumor_boundary_endpoint %>%
  apply(1, mean)

# thread line
# HE: coords_x是纵坐标，coords_y是横坐标
slope_tumor_boundary <- diff(tumor_boundary_endpoint['x', ]) / diff(tumor_boundary_endpoint['y', ])
intercept_tumor_boundary <- tumor_boundary_centerpoint[1] - slope_tumor_boundary * tumor_boundary_centerpoint[2]

# thread line: 
# coords_x = slope_tumor_boundary * coords_y + intercept_tumor_boundary
# coords_x = 2.477732 * coords_y - 515.5326

query_point_position <-
  function(x){
    stopifnot(length(x) == 2)
    
    # coords_x_on_line <- 2.477732 * x[2] - 515.5326
    coords_x_on_line <- slope_tumor_boundary * x[2] + intercept_tumor_boundary
    coords_x_query   <- x[1]
    
    if(coords_x_query > coords_x_on_line) {
      return("position1")
    }else{
      return("position2")
    }
  }


bbox2rect <- function(x){
  stopifnot(all(rownames(x) == c("x", "y")) )
  stopifnot(all(colnames(x) == c("min", "max")) )
  
  x <- as.matrix(x)
  
  data.frame(
    x = c(x['x', "min"], x['x', "max"], x['x', "min"], x['x', "max"]),
    y = c(x['y', "min"], x['y', "min"], x['y', "max"], x['y', "max"])
  )
}

# left and right boundary point
he_bbox_rect <- he_bbox %>% bbox2rect() 
he_bbox_position <- apply(he_bbox_rect, 1, query_point_position)





# HE的纵坐标是翻转的
if(sign(slope_tumor_boundary) > 0){
  print("tumor_boundary趋势1： 左上 -> 右下")
}else{
  print("tumor_boundary趋势2： 左下 -> 右上")
}



# HE绘图坐标系
#      |
# -----+--------------------------------------> x
#      |      ymin            ymax
#      |      +------------------+
#      |  xmin|                  |
#      |      |                  |
#      |      |                  |
#      |  xmax|                  |
#      |      +------------------+
#      |
#    y |
#      V


# left and right full area

coords_area_1 <-
  rbind(
    tumor_boundary_coords[c("x", "y")],
    he_bbox_rect[he_bbox_position == "position1", ]
  )

coords_area_2 <-
  rbind(
    tumor_boundary_coords[c("x", "y")],
    he_bbox_rect[he_bbox_position == "position2", ]
  )


# validate two areas
coords %>%
  ggplot(aes(x = y, y = -x)) +
  geom_point() +
  # geom_polygon(data = coords_area_1, fill = "red") +
  # geom_polygon(data = coords_area_2, fill = "blue") +
  geom_point(data = tumor_boundary_coords, color = "white")

ggplot() +
  geom_point(data = tumor_boundary_coords, aes(x = y, y = -x), color = "white")



