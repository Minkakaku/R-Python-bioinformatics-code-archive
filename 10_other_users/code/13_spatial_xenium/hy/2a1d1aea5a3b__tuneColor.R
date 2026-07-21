# use for produce and change color

# 自动生成颜色
makeColors <- function(x,
                       palette = "Dark2",
                       continuous_palette = "Oranges",
                       continuous = FALSE,
                       detect_continuous= TRUE,
                       rainbow_alpha = 0.6,
                       tuneHue = FALSE,
                       random = FALSE,
                       seed = 1234){
  
  xx <- unique(x) %>% sort()
  
  set.seed(seed + sum(nchar(as.character(xx))))
  
  if(random){
    continuous_palette <- sample(
      c('Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges',
        'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 
        'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd'),
      size = 1, replace = FALSE
    )
    
    palette <- sample(
      c('Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'),
      size = 1, replace = FALSE
    )
  }
  
  if(detect_continuous && !is.factor(x)){
    continuous <- all(str_detect(xx, "^[0-9]+$"))
  }
  
  if(is.numeric(x)){
    continuous <- TRUE
  }else{
    if(is.factor(x)){
      xx <- levels(x)
    }
  }
  
  if(continuous){
    color <- suppressWarnings(RColorBrewer::brewer.pal(length(xx), name = continuous_palette))
    if(length(color) < length(xx)){
      color <- colorRampPalette(colors = color)(length(xx))
    }else{
      color <- color[seq_len(length(xx))]
    }
    names(color) <- xx
  }else{
    color <- suppressWarnings(RColorBrewer::brewer.pal(length(xx), name = palette))
    if(length(color) < length(xx)){
      color <- rainbow(length(xx)) %>% setColorAlpha(alpha = rainbow_alpha) 
      if(tuneHue) color <- setColorHue(color)
    }else{
      color <- color[seq_len(length(xx))]
    }
    names(color) <- xx
  }
  
  return(color)
}

# 调整颜色alpha
setColorAlpha <- function(x, alpha = 0.6){
  # color: "#123456"
  x %>% col2rgb() %>% 
    apply(2, function(x) x/255) %>%  
    apply(2, function(x) rgb(x[1], x[2], x[3], alpha = 0.6))
}

# 微调颜色的色调
setColorHue <- function(x, offset = 0.1,  random = TRUE, seed = 1234, alpha = 1){
  # color: "#123456"
  set.seed(seed)
  xx <- x %>% col2rgb() %>% apply(2, function(x) rgb2hsv(x[1], x[2], x[3]))
  
  if(random){
    xx[1, ] <- xx[1, ] + rnorm(length(x), mean = 0.5, sd = 0.25)
  }else{
    xx[1, ] <- xx[1, ] + offset
  }
  
  # 0-1
  if(min(xx[1, ] < 0)){
    xx[1, ] <- xx[1, ] - min(xx[1, ])
  }
  
  if(max(xx[1, ]) > 1){
    xx[1, ] <- xx[1, ]/max(xx[1, ])
  }
  
  # return
  xxx <- xx %>% apply(2, function(x) hsv(x[1], x[2], x[3], alpha = alpha))
  
  xxx
}

# 展示颜色
checkColors <- function(x){
  # x: "#123456"
  xx <- rep(1, length(x))
  barplot(xx, col = x)
}