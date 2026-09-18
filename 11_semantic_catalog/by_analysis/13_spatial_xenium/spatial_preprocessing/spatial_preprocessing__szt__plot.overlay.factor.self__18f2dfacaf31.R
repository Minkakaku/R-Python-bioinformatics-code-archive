
get_lsgi_lin_df <- function(
    info, 
    scf = 1,
    sel.factors = NULL,
    r_squared_thresh = 0.6,
    minimum.fctr = 20
){
  lin.res.df <- get.ind.rsqrs(info)
  lin.res.df <- na.omit(lin.res.df)
  lin.res.df <- lin.res.df[lin.res.df$rsquared > r_squared_thresh,
  ]
  if (!is.null(sel.factors)) {
    lin.res.df <- lin.res.df[lin.res.df$fctr %in% sel.factors,
    ]
  }
  lin.res.df <- lin.res.df %>%
    group_by(fctr) %>%
    filter(n() >= minimum.fctr) %>%
    ungroup()
  
  lin.res.df$Xend <- lin.res.df$X + lin.res.df$vx.u
  lin.res.df$Yend <- lin.res.df$Y + lin.res.df$vy.u
  
  lin.res.df$X <- lin.res.df$X * scf
  lin.res.df$Xend <- lin.res.df$Xend * scf
  lin.res.df$Y <- lin.res.df$Y * scf
  lin.res.df$Yend <- lin.res.df$Yend * scf
  
  lin.res.df
}


plot.overlay.factor.self <- function(
    visium,
    info, 
    img = 'slice1',
    scale = c("lowres", "hires"),
    scf = 1,
    sel.factors = NULL,
    r_squared_thresh = 0.6,
    minimum.fctr = 20
){
  require(Seurat)
  require(EBImage)
  require(LSGI)
  
  scale <- match.arg(scale)
  
  lin_res_df <- get_lsgi_lin_df(
    info,
    scf              = scf,
    sel.factors      = sel.factors,
    r_squared_thresh = r_squared_thresh,
    minimum.fctr     = minimum.fctr
  )
  
  p <- 
    lin_res_df %>%
    as.data.frame %>%
    ggplot() +
    geom_segment(aes(x = X, y = Y, xend = Xend, yend = Yend, color = fctr, fill = NULL), 
                 linewidth = 0.4,
                 arrow = arrow(length = unit(0.1,"cm"))) +
    scale_color_brewer(palette = "Paired") +
    # scale_fill_gradient(low='lightgrey', high='navy')
    # +
    ggpubr::theme_transparent() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      axis.text.x = element_text(face = "bold", color = "black", size = 12, angle = 0, hjust = 1), 
      axis.text.y = element_text(face = "bold",color = "black", size = 12, angle = 0),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(colour = "black")
    )
  
  # image cropped
  stopifnot(img %in% Images(visium))
  img <- visium@images[[img]]
  coords_lowres <- img %>% GetTissueCoordinates(scale = scale)
  coords_crop <-
    coords_lowres[1:2] %>%
    map(range) %>%
    map(round) %>%
    map(~ .x[1]:.x[2])
  
  img_cropped <- img@image[coords_crop$x, coords_crop$y, ]
  
  ## image rotation
  
  eb_img <- img_cropped %>% as.Image()
  colorMode(eb_img) <- "Color"
  
  # display(eb_img) 
  
  eb_img <-
    eb_img %>% 
    rotate(90)  
  # display(eb_img) 
  
  
  # img_limits <- lin_res_df[c('X', 'Y')] %>% map(range)
  
  p_combined <-
    p +
    patchwork::inset_element(
      grid::rasterGrob(imageData(eb_img), height = unit(x = 1, units = "npc"), width = unit(x = 1, units = "npc")), 
      left   = 0, # img_limits$X[1],
      bottom = 0, # img_limits$Y[1],
      right  = 1, # img_limits$X[2],
      top    = 1, # img_limits$Y[1],
      align_to = "panel",
      on_top   = F,
      clip     = TRUE,
      ignore_tag = FALSE
    ) 
  p_combined
}



