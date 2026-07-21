den_df <- sc_n$nCount_Spatial.008um %>% density() %>% {data.frame(x = .$x, y = .$y)} 
den_df_delta <- 
  data.frame(
    delda_x = den_df$x[-1] - den_df$x[-nrow(den_df)],
    delda_y = den_df$y[-1] - den_df$y[-nrow(den_df)]
  )

d_df <- den_df_delta$delda_y/den_df_delta$delda_x

plot(den_df$x[-1], d_df)



# 测试Dimplot，颠倒及镜像

p <- visium_roi_ls_v2[[1]] %>% SpatialPlot2(alpha = 0.5, features = "CD4", image.scale = "hires", flip_y = T, flip_xy = F)
p <- visium_roi_ls[[1]] %>% SpatialPlot2(alpha = 0.5, features = "CD4", image.scale = "hires", flip_y = T, flip_xy = F, plot_he = FALSE)

p <- visium_roi_ls[[1]] %>% SpatialPlot(alpha = 0, image.scale = "hires")

p <- visium_roi_ls_v2[[1]] %>%
  SpatialPlot_manual(
    image = NULL,
    image.alpha = 1,
    image.scale = "hires",
    pt.alpha = NULL,
    crop = TRUE,
    
    slot = 'data',
    group.by = NULL,
    features = 'CD4' %>% c("CD8A", "CD3D"),
    
    flip_x = FALSE,
    flip_y = FALSE,
    flip_xy = FALSE,
    plot_he = TRUE,
    
    alpha = c(0, 1),
    
    pt.size.factor = 10,
    default_shape = 20,
    stroke = NA,
    col.by = NULL,
    alpha.by = NULL,
    na.value = 'grey50' 
  ) +
  # scale_shape_manual(values = c(
  #   'CD4'  = 21,
  #   "CD8A" = 22, 
  #   "CD3D" = 23
  # )) +
  guides(shape = guide_legend(override.aes = list(size = 5)))

# p2 <- p + scale_y_reverse()
p  # + coord_flip()


visium_roi_ls_v2[[1]] %>%
  SpatialPlot_manual(
    image.scale = "hires",
    features = 'CD4' %>% c("CD8A", "CD3D"),
    
    flip_x = FALSE,
    flip_y = FALSE,
    flip_xy = FALSE,
    plot_he = TRUE,
    alpha = c(0, 1),
    pt.size.factor = 10,
  ) +
  scale_fill_gradient(low = "#FFFFFF11", high = "red", trans = 'log1p') +
  guides(shape = guide_legend(override.aes = list(size = 5)))



g <- p %>% ggplotGrob()
gg <- g %>% gtable::gtable_filter(pattern = "panel") %>% .$grobs %>% .[[1]]

grid.newpage()
gg$children[[1]] %>% grid.draw()

gg$children[[3]] %>% grid.draw()


# p + 
#   # coord_flip()
#   scale_y_reverse()
#   # scale_y_continuous(transform = "reverse")



ll <- coords$x

0.04545455 0.95454545

-ll %>% {. - min(.)} %>% {. + min(ll)}







