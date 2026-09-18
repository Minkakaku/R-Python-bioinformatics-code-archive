
# 1. fix bug in ggplot2 ----------------------------------------------------------------
# Error in UseMethod("rescale") : "rescale"没有适用于"AsIs"目标对象的方法

# print.ggplot: original code
# just make sure function ggplot_build which is manually defined in .GlobalEnv can be called

print.ggplot <- 
function (x, newpage = is.null(vp), vp = NULL, ...) {
  # cat("print called...\n\n")
  require(grid)
  
  set_last_plot(x)
  if (newpage) 
    grid.newpage()
  grDevices::recordGraphics(requireNamespace("ggplot2", quietly = TRUE),
                            list(),
                            getNamespace("ggplot2"))
  data <- ggplot_build(x)
  gtable <- ggplot_gtable(data)
  if (is.null(vp)) {
    grid.draw(gtable)
  }
  else {
    if (is.character(vp)) 
      seekViewport(vp)
    else pushViewport(vp)
    grid.draw(gtable)
    upViewport()
  }
  if (isTRUE(getOption("BrailleR.VI")) && rlang::is_installed("BrailleR")) {
    print(asNamespace("BrailleR")$VI(x))
  }
  invisible(x)
}

# ggplot_build.ggplot: fixing code inside
ggplot_build.ggplot <- 
function (plot){
  # cat("ggplot_build.ggplot called...\n\n")
  
  plot <- ggplot2:::plot_clone(plot)
  if (length(plot$layers) == 0) {
    plot <- plot + geom_blank()
  }
  layers <- plot$layers
  layer_data <- lapply(layers, function(y) y$layer_data(plot$data))
  scales <- plot$scales
  by_layer <- function(f) {
    out <- vector("list", length(data))
    for (i in seq_along(data)) {
      out[[i]] <- f(l = layers[[i]], d = data[[i]])
    }
    out
  }
  data <- layer_data
  data <- by_layer(function(l, d) l$setup_layer(d, plot))
  layout <- ggplot2:::create_layout(plot$facet, plot$coordinates)
  data <- layout$setup(data, plot$data, plot$plot_env)
  data <- by_layer(function(l, d) l$compute_aesthetics(d, plot))
  # AsIs出现在下一条命令，原因未知
  data <- lapply(data, ggplot2:::scales_transform_df, scales = scales)
  scale_x <- function() scales$get_scales("x")
  scale_y <- function() scales$get_scales("y")
  layout$train_position(data, scale_x(), scale_y())
  data <- layout$map_position(data)
  data <- by_layer(function(l, d) l$compute_statistic(d, layout))
  data <- by_layer(function(l, d) l$map_statistic(d, plot))
  ggplot2:::scales_add_missing(plot, c("x", "y"), plot$plot_env)
  data <- by_layer(function(l, d) l$compute_geom_1(d))
  data <- by_layer(function(l, d) l$compute_position(d, layout))
  layout$reset_scales()
  layout$train_position(data, scale_x(), scale_y())
  layout$setup_panel_params()
  data <- layout$map_position(data)
  npscales <- scales$non_position_scales()
  if (npscales$n() > 0) {
    lapply(data, ggplot2:::scales_train_df, scales = npscales)
    
    # fix bug
    for(i in seq_along(data)){
      is_AsIs <- vapply(data[[i]], function(x) "AsIs" %in% class(x), FUN.VALUE = logical(1))
      if(sum(is_AsIs) > 0){
        col_with_bug <- which(is_AsIs)
        for(j in col_with_bug){
          data[[i]][[j]] <- unclass(data[[i]][[j]])
        }
      }
    }
    
    # 报错在此处，由于data中的AsIs对象的存在
    data <- lapply(data, ggplot2:::scales_map_df, scales = npscales)
  }
  data <- by_layer(function(l, d) l$compute_geom_2(d))
  data <- by_layer(function(l, d) l$finish_statistics(d))
  data <- layout$finish_data(data)
  plot$labels$alt <- get_alt_text(plot)
  structure(list(data = data, layout = layout, plot = plot), 
            class = "ggplot_built")
}

#2. fix bug cause by cowplot ---------------------------------------------------
# According to call stack, make sure ggplot_build.ggplot defined manually can be called
# ggplotGrob <- as_grob.gg <- as_gtable.gg
# class gg can raise the priority of corresponding method

ggplotGrob <- 
function (x) {
  ggplot_gtable(ggplot_build(x))
}

as_grob.gg <-
function (plot, device = NULL) {
  if (is.null(device)) {
    device <- cowplot:::null_dev_env$current
  }
  cur_dev <- grDevices::dev.cur()
  device(width = 6, height = 6)
  null_dev <- grDevices::dev.cur()
  on.exit({
    grDevices::dev.off(null_dev)
    if (cur_dev > 1) grDevices::dev.set(cur_dev)
  })
  ggplotGrob(plot)
}

as_gtable.gg <-
function (plot){
  grob <- cowplot::as_grob(plot)
  cowplot::as_gtable(grob)
}
