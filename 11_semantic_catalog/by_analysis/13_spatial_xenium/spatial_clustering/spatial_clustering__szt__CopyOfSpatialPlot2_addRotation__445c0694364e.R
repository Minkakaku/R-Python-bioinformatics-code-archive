# For plotting the tissue image
#' @importFrom ggplot2 ggproto Geom aes ggproto_parent alpha draw_key_point
#' @importFrom grid unit gpar editGrob pointsGrob viewport gTree addGrob grobName
#'
#'

library(ggplot2)
library(grid)
library(EBImage)

GeomSpatial2 <- ggproto(
  "GeomSpatial2",
  Geom,
  required_aes = c("x", "y"),
  extra_params = c("na.rm", "image", "image.alpha", "image.scale", "crop", "filp_x", "flip_y","flip_xy", "plot_he"),
  default_aes = aes(
    default_shape = 21,
    shape = NA,
    colour = "black",
    point.size.factor = 1.0,
    fill = NA,
    alpha = NA,
    stroke = NA
  ),
  setup_data = function(self, data, params) { # browser()
    
    data <- ggproto_parent(Geom, self)$setup_data(data, params)
    # We need to flip the image as the Y coordinates are reversed
    if(params$flip_xy || params$rotate_90 || params$rotate_r90){
      data$x = max(data$x) - data$x + min(data$x)
    }else{
      data$y = max(data$y) - data$y + min(data$y)
    }
    
    data
  },
  draw_key = draw_key_point,
  draw_panel = function(
    data, panel_scales, coord, image, image.alpha, image.scale, 
    crop, flip_x, flip_y, flip_xy, plot_he,
    rotate_r90, rotate_90
  ) { 
    
    # browser()
    
    if(sum(c(flip_x, flip_y, flip_xy, rotate_r90, rotate_90)) > 1){
      warning(str_glue(
        "Only one of image transformation can be set: flip_x, flip_y, flip_xy, rotate_r90, rotate_90"
      ))
    }
    
    
    img.dim <- dim(image)
    
    
    # This should be in native units, where
    # Locations and sizes are relative to the x- and yscales for the current viewport.
    if (!crop && plot_he) {
      y.transform <- c(0, img.dim[[1]]) - panel_scales$y.range
      data$y <- data$y + sum(y.transform)
      panel_scales$x$continuous_range <- c(0, img.dim[[2]])
      panel_scales$y$continuous_range <- c(0, img.dim[[1]])
      panel_scales$y.range <- c(0, img.dim[[1]])
      panel_scales$x.range <- c(0, img.dim[[2]])
    }
    
    if(flip_xy || rotate_90 || rotate_r90){
      z <- coord$transform(
        data.frame(x = c(0, img.dim[[1]]), y = c(0, img.dim[[2]])),
        panel_scales
      )
    }else{
      z <- coord$transform(
        data.frame(x = c(0, img.dim[[2]]), y = c(0, img.dim[[1]])),
        panel_scales
      )
    }
    
    # Flip X/Y axis for image
    if(flip_xy){
      z$x <- -rev(z$x) + 1
    }else if(flip_x || flip_y){
      z$y <- -rev(z$y) + 1
    }else if(rotate_90){
      z$x <- -rev(z$x) + 1
      z$y <- -rev(z$y) + 1
    }else if(rotate_r90){
      
    }else{
      z$y <- -rev(z$y) + 1
    }
    
    eb_img <- EBImage::Image(image@image)
    img <- Seurat:::GetImage.VisiumV2(image)
    
    cat('flip_x: ', flip_x, "\n")
    cat('flip_y: ', flip_y, "\n")
    cat('flip_xy: ', flip_xy, "\n")
    cat('rotate_90: ', rotate_90, "\n")
    cat('rotate_r90: ', rotate_r90, "\n")
    
    if(flip_y && plot_he) {
      z$y <- -rev(z$y) + 1
      eb_img <- eb_img %>% EBImage::flop() # y: 函数是水平翻转
      img$raster <- as.raster(imageData(eb_img))
    }
    if(flip_x && plot_he) {
      z$x <- -rev(z$x) + 1
      eb_img <- eb_img %>% EBImage::flip() # x: 函数是垂直翻转
      img$raster <- as.raster(imageData(eb_img))
    }
    
    if(flip_xy && plot_he){
      new_image <- func_flip_xy(eb_img)
      img$raster <- as.raster(new_image)
    }
    
    if(rotate_90 && plot_he){
      
      new_image <- func_rotate_90(eb_img)
      img$raster <- as.raster(new_image)
    }
    
    if(rotate_r90 && plot_he){
      
      new_image <- func_rotate_r90(eb_img)
      img$raster <- as.raster(new_image)
    }
    
    wdth <- z$x[2] - z$x[1]
    hgth <- z$y[2] - z$y[1]
    vp <- viewport(
      x = unit(x = z$x[1], units = "npc"),
      y = unit(x = z$y[1], units = "npc"),
      width = unit(x = wdth, units = "npc"),
      height = unit(x = hgth, units = "npc"),
      just = c("left", "bottom")
    )
    
    spot.size <- Seurat:::Radius.VisiumV2(object = image, scale = image.scale)
    
    # translate data to shape
    if (is.character(data$shape)) {
      data$shape <- ggplot2::translate_shape_string(data$shape)
    }
    
    coords <- coord$transform(data, panel_scales)
    
    img <- editGrob(grob = img, vp = vp)
    
    # flip coords of points
    if(flip_x || rotate_r90){
      coords$x <- max(coords$x) - coords$x + min(coords$x)
    }
    if(flip_y || rotate_90){
      coords$y <- max(coords$y) - coords$y + min(coords$y)
    }
    # if(rotate_90){
    #   coords$y <- max(coords$y) - coords$y + min(coords$y)
    # }
    # if(rotate_r90){
    #   coords$x <- max(coords$x) - coords$x + min(coords$x)
    # }
    
    if(n_distinct(data$shape) > 1){
      pch <- coords$shape # turn to shape 21 - 25
    }else{
      pch <- data$default_shape
    }
    
    pts <- pointsGrob(
      x = coords$x,
      y = coords$y,
      pch = pch,
      size = unit(spot.size, "npc") * data$point.size.factor * 10,
      gp = gpar(
        # col = alpha(colour = coords$colour, alpha = coords$alpha),
        col = alpha(colour = coords$fill, alpha = coords$alpha),
        fill = alpha(colour = coords$fill, alpha = coords$alpha),
        lwd = coords$stroke)
    )
    vp <- viewport()
    gt <- gTree(vp = vp)
    if (image.alpha > 0 && plot_he) {
      if (image.alpha != 1) {
        img$raster = as.raster(
          x = matrix(
            data = alpha(colour = img$raster, alpha = image.alpha),
            nrow = nrow(x = img$raster),
            ncol = ncol(x = img$raster),
            byrow = TRUE)
        )
      }
      gt <- addGrob(gTree = gt, child = img)
    }
    gt <- addGrob(gTree = gt, child = pts)
    # Replacement for ggname
    gt$name <- grobName(grob = gt, prefix = 'geom_spatial2')
    return(gt)
  }
)


# flip xy 
func_flip_xy <- function(eb_img){
  
  
  old_image <- imageData(eb_img)
  dim_img <- dim(old_image)
  
  stopifnot(all(dim_img > 0))
  
  new_image <- array(0, dim = dim_img[c(2, 1, 3)] )
  
  i_seq <- seq.int(dim_img[1])
  i_seq_rev <- max(i_seq) - i_seq + min(i_seq)
  for(i in i_seq){
    
    j_seq <- seq.int(dim_img[2])
    j_seq_rev <- max(j_seq) - j_seq + min(j_seq)
    
    for(j in j_seq){
      
      new_image[j, i, ] <- old_image[i_seq_rev[i], j_seq_rev[j], ]
    }
  }
  
  new_image
}

# rotate -90 degree
func_rotate_r90 <- function(eb_img){
  
  
  old_image <- imageData(eb_img)
  dim_img <- dim(old_image)
  
  stopifnot(all(dim_img > 0))
  
  new_image <- array(0, dim = dim_img[c(2, 1, 3)] )
  
  i_seq <- seq.int(dim_img[1])
  i_seq_rev <- max(i_seq) - i_seq + min(i_seq)
  for(i in i_seq){
    
    j_seq <- seq.int(dim_img[2])
    j_seq_rev <- max(j_seq) - j_seq + min(j_seq)
    
    for(j in j_seq){
      
      new_image[j, i, ] <- old_image[i, j_seq_rev[j], ]
    }
  }
  
  new_image
}

# rotate 90 degree
func_rotate_90 <- function(eb_img){
  
  
  old_image <- imageData(eb_img)
  dim_img <- dim(old_image)
  
  stopifnot(all(dim_img > 0))
  
  new_image <- array(0, dim = dim_img[c(2, 1, 3)] )
  
  i_seq <- seq.int(dim_img[1])
  i_seq_rev <- max(i_seq) - i_seq + min(i_seq)
  for(i in i_seq){
    
    j_seq <- seq.int(dim_img[2])
    j_seq_rev <- max(j_seq) - j_seq + min(j_seq)
    
    for(j in j_seq){
      
      new_image[j, i, ] <- old_image[i_seq_rev[i], j, ]
    }
  }
  
  new_image
}





# influenced by: https://stackoverflow.com/questions/49475201/adding-tables-to-ggplot2-with-facet-wrap-in-r
# https://ggplot2.tidyverse.org/articles/extending-ggplot2.html
#' @importFrom ggplot2 layer
#'
#'
library(ggplot2)
geom_spatial2 <-  function(
    mapping = NULL,
    data = NULL,
    image = image,
    image.alpha = image.alpha,
    image.scale = image.scale,
    crop = crop,
    flip_x = FALSE,
    flip_y = FALSE,
    flip_xy = FALSE,
    rotate_r90 = FALSE,
    rotate_90 = FALSE,
    plot_he = TRUE,
    stat = "identity",
    position = "identity",
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE,
    ...
) {
  layer(
    geom = GeomSpatial2,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      image = image,
      image.alpha = image.alpha,
      image.scale = image.scale,
      crop = crop,
      plot_he = plot_he,
      flip_x = flip_x,
      flip_y = flip_y,
      flip_xy = flip_xy,
      rotate_r90 = rotate_r90,
      rotate_90 = rotate_90,
      ...
    )
  )
}



#
library(ggplot2)
library(tibble)
SingleSpatialPlot2 <- function(
    data,
    image,
    cols = NULL,
    image.alpha = 1,
    image.scale = "lowres",
    pt.alpha = NULL,
    crop = TRUE,
    
    flip_x = FALSE,
    flip_y = FALSE,
    flip_xy = FALSE,
    plot_he = TRUE,
    
    pt.size.factor = NULL,
    shape = 21,
    stroke = NA,
    col.by = NULL,
    alpha.by = NULL,
    cells.highlight = NULL,
    cols.highlight = c('#DE2D26', 'grey50'),
    geom = c('spatial', 'interactive', 'poly'),
    na.value = 'grey50'
) {
  
  
  browser()
  geom <- match.arg(arg = geom)
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find '", col.by, "' in data, not coloring", call. = FALSE, immediate. = TRUE)
    col.by <- NULL
  }
  col.by <- col.by %iff% paste0("`", col.by, "`")
  alpha.by <- alpha.by %iff% paste0("`", alpha.by, "`")
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(
      cells.highlight = cells.highlight,
      cells.all = rownames(x = data),
      sizes.highlight = pt.size.factor,
      cols.highlight = cols.highlight[1],
      col.base = cols.highlight[2]
    )
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- 'highlight'
    levels(x = data$ident) <- c(order, setdiff(x = levels(x = data$ident), y = order))
    data <- data[order(data$ident), ]
  }
  
  if(flip_xy){
    plot <- ggplot(data = data, aes_string(
      x = colnames(x = data)[1],
      y = colnames(x = data)[2],
      fill = col.by,
      alpha = alpha.by
    ))
  }else{
    plot <- ggplot(data = data, aes_string(
      x = colnames(x = data)[2],
      y = colnames(x = data)[1],
      fill = col.by,
      alpha = alpha.by
    ))
  }
  
  
  plot <- switch(
    EXPR = geom,
    'spatial' = {
      if (is.null(x = pt.alpha)) {
        plot <- plot + geom_spatial2(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          image.scale = image.scale,
          crop = crop,
          flip_x = flip_x,
          flip_y = flip_y,
          flip_xy = flip_xy,
          plot_he = plot_he,
          shape = shape,
          stroke = stroke,
        )
      } else {
        plot <- plot + geom_spatial2(
          point.size.factor = pt.size.factor,
          data = data,
          image = image,
          image.alpha = image.alpha,
          image.scale = image.scale,
          crop = crop,
          flip_x = flip_x,
          flip_y = flip_y,
          flip_xy = flip_xy,
          plot_he = plot_he,
          shape = shape,
          stroke = stroke,
          alpha = pt.alpha
        )
      }
      plot + coord_fixed() + theme(aspect.ratio = 1)
    },
    'interactive' = {
      plot + geom_spatial_interactive(
        data = tibble(
          grob = list(
            GetImage(
              object = image,
              mode = 'grob'
            )
          )
        ),
        mapping = aes_string(grob = 'grob'),
        x = 0.5,
        y = 0.5
      ) +
        geom_point(mapping = aes_string(color = col.by)) +
        xlim(0, ncol(x = image)) +
        ylim(nrow(x = image), 0) +
        coord_cartesian(expand = FALSE)
    },
    'poly' = {
      data$cell <- rownames(x = data)
      data[, c('x', 'y')] <- NULL
      data <- merge(
        x = data,
        y = GetTissueCoordinates(
          object = image,
          qhulls = TRUE,
          scale = image.scale
        ),
        by = "cell"
      )
      plot + geom_polygon(
        data = data,
        mapping = aes_string(fill = col.by, group = 'cell')
      ) + coord_fixed() + theme_cowplot()
      
    },
    stop("Unknown geom, choose from 'spatial' or 'interactive'", call. = FALSE)
  )
  if (!is.null(x = cells.highlight)) {
    plot <- plot + scale_fill_manual(values = cols.highlight)
  }
  if (!is.null(x = cols) && is.null(x = cells.highlight)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || cols %in% rownames(x = RColorBrewer::brewer.pal.info))) {
      scale <- scale_fill_brewer(palette = cols, na.value = na.value)
    } else if (length(x = cols) == 1 && (cols %in% c('alphabet', 'alphabet2', 'glasbey', 'polychrome', 'stepped'))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), palette = cols)
      scale <- scale_fill_manual(values = colors, na.value = na.value)
    } else {
      cols <- cols[names(x = cols) %in% data[[gsub(pattern = '`', replacement = "", x = col.by)]]]
      scale <- scale_fill_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + NoAxes() + theme(panel.background = element_blank())
  return(plot)
}


#' Visualize spatially  
#' 
library(ggplot2)
library(patchwork)
library(scales)
SpatialPlot2 <- function(
    object,
    group.by = NULL,
    features = NULL,
    images = NULL,
    cols = NULL,
    image.alpha = 1,
    image.scale = "lowres",
    crop = TRUE,
    slot = 'data',
    keep.scale = "feature",
    min.cutoff = NA,
    max.cutoff = NA,
    cells.highlight = NULL,
    cols.highlight = c('#DE2D26', 'grey50'),
    facet.highlight = FALSE,
    label = FALSE,
    label.size = 5,
    label.color = 'white',
    label.box = TRUE,
    repel = FALSE,
    ncol = NULL,
    combine = TRUE,
    pt.size.factor = 1.6,
    alpha = c(1, 1),
    shape = 21,
    stroke = NA,
    flip_x = FALSE,
    flip_y = FALSE,
    flip_xy = FALSE,
    plot_he = TRUE,
    interactive = FALSE,
    do.identify = FALSE,
    identify.ident = NULL,
    do.hover = FALSE,
    information = NULL
) {
  if (isTRUE(x = do.hover) || isTRUE(x = do.identify)) {
    warning(
      "'do.hover' and 'do.identify' are deprecated as we are removing plotly-based interactive graphics, use 'interactive' instead for Shiny-based interactivity",
      call. = FALSE,
      immediate. = TRUE
    )
    interactive <- TRUE
  }
  if (!is.null(x = group.by) & !is.null(x = features)) {
    stop("Please specific either group.by or features, not both.")
  }
  images <- images %||% Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) == 0) {
    images <- Images(object = object)
  }
  if (length(x = images) < 1) {
    stop("Could not find any spatial image information")
  }
  
  # Check keep.scale param for valid entries
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature", "all"))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  
  cells <- unique(Seurat:::CellsByImage(object, images = images, unlist = TRUE))
  if (is.null(x = features)) {
    if (interactive) {
      return(ISpatialDimPlot(
        object = object,
        image = images[1],
        image.scale = image.scale,
        group.by = group.by,
        alpha = alpha
      ))
    }
    group.by <- group.by %||% 'ident'
    object[['ident']] <- Idents(object = object)
    data <- object[[group.by]]
    data <- data[cells,,drop=F]
    for (group in group.by) {
      if (!is.factor(x = data[, group])) {
        data[, group] <- factor(x = data[, group])
      }
    }
  } else {
    if (interactive) {
      return(ISpatialFeaturePlot(
        object = object,
        feature = features[1],
        image = images[1],
        image.scale = image.scale,
        slot = slot,
        alpha = alpha
      ))
    }
    data <- FetchData(
      object = object,
      vars = features,
      cells = cells,
      layer = slot,
      clean = FALSE
    )
    features <- colnames(x = data)
    # Determine cutoffs
    min.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        return(ifelse(
          test = is.na(x = cutoff),
          yes = min(data[, feature]),
          no = cutoff
        ))
      },
      cutoff = min.cutoff,
      feature = features
    )
    max.cutoff <- mapply(
      FUN = function(cutoff, feature) {
        return(ifelse(
          test = is.na(x = cutoff),
          yes = max(data[, feature]),
          no = cutoff
        ))
      },
      cutoff = max.cutoff,
      feature = features
    )
    check.lengths <- unique(x = vapply(
      X = list(features, min.cutoff, max.cutoff),
      FUN = length,
      FUN.VALUE = numeric(length = 1)
    ))
    if (length(x = check.lengths) != 1) {
      stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    # Apply cutoffs
    data <- sapply(
      X = 1:ncol(x = data),
      FUN = function(index) {
        data.feature <- as.vector(x = data[, index])
        min.use <- SetQuantile(cutoff = min.cutoff[index], data.feature)
        max.use <- SetQuantile(cutoff = max.cutoff[index], data.feature)
        data.feature[data.feature < min.use] <- min.use
        data.feature[data.feature > max.use] <- max.use
        return(data.feature)
      }
    )
    colnames(x = data) <- features
    rownames(x = data) <- cells
  }
  features <- colnames(x = data)
  colnames(x = data) <- features
  rownames(x = data) <- cells
  facet.highlight <- facet.highlight && (!is.null(x = cells.highlight) && is.list(x = cells.highlight))
  if (do.hover) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning(
        "'do.hover' requires only one image, using image ",
        images,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (length(x = features) > 1) {
      features <- features[1]
      type <- ifelse(test = is.null(x = group.by), yes = 'feature', no = 'grouping')
      warning(
        "'do.hover' requires only one ",
        type,
        ", using ",
        features,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (facet.highlight) {
      warning(
        "'do.hover' requires no faceting highlighted cells",
        call. = FALSE,
        immediate. = TRUE
      )
      facet.highlight <- FALSE
    }
  }
  if (facet.highlight) {
    if (length(x = images) > 1) {
      images <- images[1]
      warning(
        "Faceting the highlight only works with a single image, using image ",
        images,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    ncols <- length(x = cells.highlight)
  } else {
    ncols <- length(x = images)
  }
  plots <- vector(
    mode = "list",
    length = length(x = features) * ncols
  )
  
  # Get max across all features
  if (!(is.null(x = keep.scale)) && keep.scale == "all") {
    max.feature.value <- max(apply(data, 2, function(x) max(x, na.rm = TRUE)))
  }
  
  for (i in 1:ncols) {
    plot.idx <- i
    image.idx <- ifelse(test = facet.highlight, yes = 1, no = i)
    image.use <- object[[images[[image.idx]]]]
    
    coordinates <- GetTissueCoordinates(
      object = image.use,
      scale = image.scale
    )
    
    highlight.use <- if (facet.highlight) {
      cells.highlight[i]
    } else {
      cells.highlight
    }
    for (j in 1:length(x = features)) {
      cols.unset <- is.factor(x = data[, features[j]]) && is.null(x = cols)
      if (cols.unset) {
        cols <- hue_pal()(n = length(x = levels(x = data[, features[j]])))
        names(x = cols) <- levels(x = data[, features[j]])
      }
      
      # Get feature max for individual feature
      if (!(is.null(x = keep.scale)) && keep.scale == "feature" && !inherits(x = data[, features[j]], what = "factor") ) {
        max.feature.value <- max(data[, features[j]])
      }
      
      plot <- SingleSpatialPlot2(
        data = cbind(
          coordinates,
          data[rownames(x = coordinates), features[j], drop = FALSE]
        ),
        image = image.use,
        image.scale = image.scale,
        image.alpha = image.alpha,
        col.by = features[j],
        cols = cols,
        alpha.by = if (is.null(x = group.by)) {
          features[j]
        } else {
          NULL
        },
        pt.alpha = if (!is.null(x = group.by)) {
          alpha[j]
        } else {
          NULL
        },
        geom = if (inherits(x = image.use, what = "STARmap")) {
          'poly'
        } else {
          'spatial'
        },
        cells.highlight = highlight.use,
        cols.highlight = cols.highlight,
        pt.size.factor = pt.size.factor,
        shape = shape,
        stroke = stroke,
        crop = crop,
        
        flip_x = flip_x,
        flip_y = flip_y,
        flip_xy = flip_xy,
        plot_he = plot_he
      )
      
      # browser()
      if (is.null(x = group.by)) {
        plot <- plot +
          scale_fill_gradientn(
            name = features[j],
            colours = Seurat:::SpatialColors(n = 100)
          ) +
          theme(legend.position = 'top') +
          scale_alpha(range = alpha) +
          guides(alpha = "none")
      } else if (label) {
        plot <- LabelClusters(
          plot = plot,
          id = ifelse(
            test = is.null(x = cells.highlight),
            yes = features[j],
            no = 'highlight'
          ),
          geom = if (inherits(x = image.use, what = "STARmap")) {
            'GeomPolygon'
          } else {
            'GeomSpatial2'
          },
          repel = repel,
          size = label.size,
          color = label.color,
          box = label.box,
          position = "nearest"
        )
      }
      if (j == 1 && length(x = images) > 1 && !facet.highlight) {
        plot <- plot +
          ggtitle(label = images[[image.idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
      }
      if (facet.highlight) {
        plot <- plot +
          ggtitle(label = names(x = cells.highlight)[i]) +
          theme(plot.title = element_text(hjust = 0.5)) +
          NoLegend()
      }
      
      # Plot multiple images depending on keep.scale
      if (!(is.null(x = keep.scale)) && !inherits(x = data[, features[j]], "factor")) {
        plot <- suppressMessages(plot & scale_fill_gradientn(colors = Seurat:::SpatialColors(n = 100), limits = c(NA, max.feature.value)))
      }
      
      plots[[plot.idx]] <- plot
      plot.idx <- plot.idx + ncols
      if (cols.unset) {
        cols <- NULL
      }
    }
  }
  
  if (combine) {
    if (!is.null(x = ncol)) {
      return(wrap_plots(plots = plots, ncol = ncol))
    }
    if (length(x = images) > 1) {
      return(wrap_plots(plots = plots, ncol = length(x = images)))
    }
    return(wrap_plots(plots = plots))
  }
  return(plots)
}



## spatial plot
SpatialPlot_manual <- function(
    object,
    image = NULL,
    image.alpha = 1,
    image.scale = "lowres",
    pt.alpha = NULL,
    crop = TRUE,
    
    slot = 'data',
    group.by = NULL,
    features = NULL,
    
    flip_x = FALSE,
    flip_y = FALSE,
    flip_xy = FALSE,
    rotate_r90 = FALSE,
    rotate_90  = FALSE,
    plot_he    = TRUE,
    
    pt.size.factor = NULL,
    default_shape = 21,
    alpha = c(1, 1),
    stroke = NA,
    col.by = NULL,
    alpha.by = NULL,
    na.value = 'grey50'    
){
  # browser()
  
  if(sum(c(flip_x, flip_y, flip_xy, rotate_r90, rotate_90)) > 1){
    warning(str_glue(
      "Only one of image transformation can be set: flip_x, flip_y, flip_xy, rotate_r90, rotate_90"
    ))
  } 
  
  if(!is.null(features)){
    data <- FetchData(
      object = object,
      vars  = features,
      layer = slot,
      clean = FALSE
    )
    
    if(is.null(col.by)) col.by <- features
    if(is.null(alpha.by)) alpha.by <- features
    
  }else if(!is.null(group.by)){
    data <- object[[group.by]]
    
    if(is.null(col.by)) col.by <- group.by
    # if(is.null(alpha.by)) alpha.by <- NULL
  }else {
    data <- NULL
  }
  
  
  if(is.null(image)) {
    image <-  Images(object)[1]
  }
  
  if(is.character(image)){
    image.use <- object[[image]]
  }else{
    image.use <- image
  }
  
  
  coordinates <- GetTissueCoordinates(
    object = image.use,
    scale  = image.scale
  )
  
  if(!is.null(data)){
    data <- cbind(
      coordinates,
      data[rownames(coordinates), , drop = FALSE]
    )
  }else{
    data <- coordinates
  }
  
  
  # multi features
  if(length(features) > 1){
    data <- data %>%
      pivot_longer(
        any_of(features),
        names_to  = "gene",
        values_to = "expression"
      )
    shape.by <- "gene"
    col.by   <- 'expression'
    alpha.by <- 'expression'
  }else{
    shape.by <- NULL
  }
  
  
  if(flip_xy || rotate_90 || rotate_r90){
    plot <- ggplot(data = data, aes_string(
      x = colnames(x = data)[1],
      y = colnames(x = data)[2],
      fill = col.by,
      alpha = alpha.by,
      shape = shape.by
    ))
  }else{
    plot <- ggplot(data = data, aes_string(
      x = colnames(x = data)[2],
      y = colnames(x = data)[1],
      fill = col.by,
      alpha = alpha.by,
      shape = shape.by
    ))
  }
  
  if (is.null(pt.alpha)) {
    plot <- plot + geom_spatial2(
      point.size.factor = pt.size.factor,
      data = data,
      image = image.use,
      image.alpha = image.alpha,
      image.scale = image.scale,
      crop = crop,
      flip_x = flip_x,
      flip_y = flip_y,
      flip_xy = flip_xy,
      rotate_r90 = rotate_r90,
      rotate_90 = rotate_90,
      plot_he = plot_he,
      default_shape = default_shape,
      stroke = stroke,
    )
  } else {
    plot <- plot + geom_spatial2(
      point.size.factor = pt.size.factor,
      data = data,
      image = image.use,
      image.alpha = image.alpha,
      image.scale = image.scale,
      crop = crop,
      flip_x = flip_x,
      flip_y = flip_y,
      flip_xy = flip_xy,
      rotate_r90 = rotate_r90,
      rotate_90 = rotate_90,
      plot_he = plot_he,
      default_shape = default_shape,
      stroke = stroke,
      alpha = pt.alpha
    )
  }
  
  # browser()
  
  # return
  plot + 
    NoAxes() +
    coord_fixed() +
    scale_fill_gradientn(
      # name = features[j],
      colours = Seurat:::SpatialColors(n = 100)
    ) +
    theme(legend.position = 'top') +
    scale_alpha(range = alpha) +
    
    guides(alpha = "none") +
    
    theme(panel.background = element_blank()) +
    theme(aspect.ratio = 1)
}


