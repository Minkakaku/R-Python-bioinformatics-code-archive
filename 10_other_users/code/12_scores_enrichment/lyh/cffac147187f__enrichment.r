# 代码来源：https://www.r2omics.cn/
# 加载必要库
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

# 核心绘图函数
circosEnrichmentPlot <- function(df,
                                 topN = 6,
                                 classCol = c("#f7cb16", "#65c3fc", "#bfe046", "#bdb5e3", "#fbccca", "#54beaa"),
                                 ifLog = FALSE,
                                 type = "base") {

  df.top  <- df %>%
    group_by(Class) %>%
    slice_min(PValue, n = topN, with_ties = FALSE) %>%
    arrange(desc(Ratio), .by_group = TRUE) %>%
    ungroup()

  main.col  <- classCol[as.numeric(as.factor(df.top$Class))]

  # 创建四个数据层
  df1 <- data.frame(TermID  = df.top$TermID,  start = 0, end = max(df.top$bg_pro_num,  na.rm  = TRUE))

  # 颜色映射
  p_max <- round(max(df.top$Enrichment))  + 1
  color_assign <- colorRamp2(0:p_max, colorRampPalette(c('#fee5d9', '#fb6a4a'))(p_max + 1))

  df2 <- data.frame(
    TermID = df.top$TermID,
    start = 0,
    end = df.top$bg_term_num,
    bg_term_num = df.top$bg_term_num,
    bg_term_num_col = color_assign(df.top$Enrichment)
  )

  # 根据类型处理前景数据
  if (type == "base") {
    df3 <- data.frame(
      TermID = df.top$TermID,
      start = 0,
      end = df.top$fg_term_num,
      fg_term_num = df.top$fg_term_num,
      color = "#ba55d3"
    )
  } else if (type == "sig") {
    tempLong <- if (ifLog) log10(max(df.top$bg_pro_num,  na.rm  = TRUE) + 1) else max(df.top$bg_pro_num,  na.rm  = TRUE)

    df3 <- bind_rows(
      data.frame(
        TermID = df.top$TermID,
        start = 0,
        end = df.top$Up  / (df.top$Up  + df.top$Down)  * tempLong,
        count = df.top$Up,
        color = "#69115dA1"
      ),
      data.frame(
        TermID = df.top$TermID,
        start = df.top$Up  / (df.top$Up  + df.top$Down)  * tempLong,
        end = tempLong,
        count = df.top$Down,
        color = "#838bc5"
      )
    ) %>% mutate(count = ifelse(count == 0, "", count))
  } else {
    stop("type must be 'base' or 'sig'")
  }

  # 比例数据
  df4 <- data.frame(
    TermID = df.top$TermID,
    start = 0,
    end = max(df.top$bg_pro_num,  na.rm  = TRUE),
    ratio = df.top$Ratio  / max(df.top$Ratio,  na.rm  = TRUE) * 10,
    col = main.col
  )

  # 坐标轴对数转换
  if (ifLog) {
    df1 <- mutate(df1, end = log10(end + 1))
    df2 <- mutate(df2, end = log10(end + 1))
    if (type == "base") df3 <- mutate(df3, end = log10(end + 1))
    df4 <- mutate(df4, end = log10(end + 1))
  }

  # 绘图设置
  par(omi = c(0.1, 0.1, 0.1, 1.5))
  circos.par(track.margin  = c(0.01, 0.01))

  # 初始化轨道
  circos.genomicInitialize(df1,  plotType = "none")

  # 第一轨道：分类标签
  circos.trackPlotRegion(
    ylim = c(0, 1),
    panel.fun  = function(x, y) {
      sector.index  <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(mean(xlim),  mean(ylim), sector.index,
                  cex = 0.8, facing = "bending.inside",  niceFacing = TRUE)
    },
    track.height  = 0.08,
    bg.border  = NA,
    bg.col  = main.col
  )

  # 坐标轴
  if (!ifLog) {
    for (si in get.all.sector.index())  {
      circos.axis(
        h = "top",
        labels.cex  = 0.6,
        sector.index  = si,
        track.index  = 1,
        major.at  = pretty(c(0, max(df1$end, na.rm  = TRUE)), n = 5),
        labels.facing  = "clockwise"
      )
    }
  } else {
    for (si in get.all.sector.index())  {
      circos.axis(
        h = "top",
        labels.cex  = 0.6,
        sector.index  = si,
        track.index  = 1,
        major.at  = 0:5,
        labels = c(0, 10, 100, 1000, 10000),
        labels.facing  = "clockwise"
      )
    }
  }

  # 第二轨道：背景基因数
  circos.genomicTrack(
    df2,
    ylim = c(0, 1),
    track.height  = 0.1,
    bg.border  = "white",
    panel.fun  = function(region, value, ...) {
      circos.genomicRect(region,  value, ytop = 0, ybottom = 1,
                         col = value[, 2], border = NA, ...)
      circos.genomicText(region,  value, y = 0.5, labels = value[, 1],
                         adj = c(0.5, 0.5), cex = 0.8, ...)
    }
  )

  # 第三轨道：前景基因数
  circos.genomicTrack(
    df3,
    ylim = c(0, 1),
    track.height  = 0.1,
    bg.border  = "white",
    panel.fun  = function(region, value, ...) {
      circos.genomicRect(region,  value, ytop = 0, ybottom = 1,
                         col = value[, 2], border = NA, ...)
      circos.genomicText(region,  value, y = 0.5, labels = value[, 1],
                         cex = 0.9, adj = c(0.5, 0.5), ...)
    }
  )

  # 第四轨道：富集比例
  circos.genomicTrack(
    df4,
    ylim = c(0, 10),
    track.height  = 0.35,
    bg.border  = "white",
    bg.col  = "grey90",
    panel.fun  = function(region, value, ...) {
      cell.xlim  <- get.cell.meta.data("cell.xlim")
      cell.ylim  <- get.cell.meta.data("cell.ylim")
      for (j in 1:9) {
        y <- cell.ylim[1]  + (cell.ylim[2]  - cell.ylim[1])  / 10 * j
        circos.lines(cell.xlim,  c(y, y), col = "#FFFFFF", lwd = 0.3)
      }
      circos.genomicRect(region,  value, ytop = 0, ybottom = value[, 1],
                         col = value[, 2], border = NA, ...)
    }
  )

  circos.clear()

  # 添加图例
  if (type == "base") {
    middle.legend  <- Legend(
      labels = c('Number of Genes', 'Number of Select', 'Rich Factor(0-1)'),
      type = "points",
      pch = c(15, 15, 17),
      legend_gp = gpar(col = c('pink', '#BA55D3', main.col[1])),
      title = "",
      nrow = 3,
      size = unit(3, "mm")
    )
  } else {
    middle.legend  <- Legend(
      labels = c('Number of Genes', 'Number of Up', 'Number of Down', 'Rich Factor(0-1)'),
      type = "points",
      pch = c(15, 15, 15, 17),
      legend_gp = gpar(col = c('pink', '#69115dA1', '#838bc5', main.col[1])),
      title = "",
      nrow = 4,
      size = unit(3, "mm")
    )
  }

  circle_size <- unit(1, "snpc")
  draw(middle.legend,  x = circle_size * 0.42)

  # 主图例
  main.legend  <- Legend(
    labels = unique(df.top$Class),
    type = "points",
    pch = 15,
    legend_gp = gpar(col = classCol),
    title = "Class",
    nrow = 3,
    size = unit(3, "mm"),
    grid_height = unit(5, "mm"),
    grid_width = unit(5, "mm")
  )

  # P值图例
  logp.legend  <- Legend(
    col_fun = colorRamp2(
      round(seq(0, p_max, length.out  = 6), 0),
      colorRampPalette(c('#fee5d9', '#fb6a4a'))(6)
    ),
    legend_height = unit(3, 'cm'),
    labels_gp = gpar(fontsize = 8),
    title_position = 'topleft',
    title = '-Log10(Pvalue)'
  )

  lgd <- packLegend(main.legend,  logp.legend)
  draw(lgd, x = circle_size * 0.85, y = circle_size * 0.55, just = "left")
}