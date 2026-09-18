.libPaths('/mnt/NexData_file/home/wyt/software/minicoda/envs/R_envs/lib/R/library')
library(circlize)
colors_36 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", 
  "#FFFF33", "#F781BF", "#999999", "#E5D8BD", "#B3CDE3", "#CCEBC5", 
  "#FED9A6", "#FBB4AE", "#8DD3C7", "#BEBADA", "#80B1D3", "#B3DE69", 
  "#FCCDE5", "#BC80BD", "#FFED6F", "#8DA0CB", "#E78AC3", "#E5C494", 
  "#CCCCCC", "#FB9A99", "#800080", "#CAB2D6", "#6A3D9A","#E7300A",
  "#66C2A5","#FC8D62","#A6D854", "#FFD92F", "#E7298A","#008B8B", "#FF00FF" 
)
colors_46 <- c(colors_36,
               "#2F4F4F",
               "#556B2F",
               "#4682B4",
               "#D2691E",
               "#5F9EA0",
               "#DA70D6",
               "#6495ED",
               "#DC143C",
               "#00CED1")
colors <-  colors_46
circos_celltype_plot <- function(data,source.use,bottom_value,right_value1,right_value2,title_value){
  
  if(dev.cur() == 1) {
    dev.new()  # 如果没有图形设备，创建一个
  }
  
  set.seed(123)
  merged_data_sub <- subset(data,source %in% source.use)
  cell_types <- unique(c(merged_data_sub$source, merged_data_sub$target))
  cols <- colors
  names(cols) <- cell_types
  link_colors <- ifelse(merged_data_sub$prob_diff > 0, "#E41A1C", "#377EB8")
  names(link_colors) <- as.character(merged_data_sub$source)
  # 2. 设置图形参数
  par(mar = c(2, 2, 3, 13))  # 右边留出图例空间
  
  # 3. 创建和弦图，隐藏扇区文字标签
  
  chordDiagram(merged_data_sub[, c("source", "target", "prob_diff")],
               grid.col = cols,
               annotationTrack = "grid")
  
  actual_order <- get.all.sector.index()
  names(colors) <- actual_order
  n_cells <- length(actual_order)
  n_first_col <- ceiling(n_cells / 2)
  col1 <- actual_order[1:n_first_col]
  col2 <- actual_order[(n_first_col + 1):n_cells]
  circos.clear()
  circos.par(gap.after = 5)
  
  chordDiagram(merged_data_sub[, c("source", "target", "prob_diff")],
               order = actual_order,
               directional = 1,
               direction.type = "arrows",
               link.arr.type = "big.arrow",
               grid.col = colors[actual_order],
               col = link_colors, 
               transparency = 0.5,
               reduce = -1,  
               annotationTrack = c("grid"),
               annotationTrackHeight = c(0.03, 0.05),
               preAllocateTracks = list(
                 list(track.height = 0.1),
                 list(track.height = 0.1)
               )
  )
  
  circos.track(track.index = 2, panel.fun = function(x, y) {
    sector.index = CELL_META$sector.index
    if (grepl(":", sector.index)) {
      gene = strsplit(sector.index, ":")[[1]][2]
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], gene,
                  facing = "clockwise", niceFacing = TRUE,
                  adj = c(0, 0.5), cex = 0.6, col = "darkblue")
    }
  }, bg.border = NA)
  
  title(main = paste(source.use,"Cell-Cell Communication Chord Diagram"),
        line = title_value,      # 适中距离
        cex.main = 0.9,
        col.main = "black",
        adj = 0.55)
  # 第一列图例
  legend("right",
         legend = col1,
         fill = colors[col1],
         border = "black",
         title = "Cell Types",
         cex = 0.7,
         xpd = TRUE,
         inset = c(right_value1, 0),  # 向右移动
         bty = "n",
         x.intersp = 0.3,
         y.intersp = 0.6)
  
  # 第二列图例（稍微向右移动）
  legend("right",
         legend = col2,
         fill = colors[col2],
         border = "black",
         title = "",
         cex = 0.7,
         xpd = TRUE,
         inset = c(right_value2, 0),  # 更向右移动
         bty = "n",
         x.intersp = 0.3,
         y.intersp = 0.6)
  
  legend("bottom",
         legend = c("increase", "decrease"),
         fill = c("#E41A1C", "#377EB8"),
         border = "black",
         title = "",
         cex = 0.7,
         xpd = TRUE,
         inset = c(0,bottom_value),  # 负值向下移动
         bty = "n",
         horiz = TRUE,
         x.intersp = 0.5,
         y.intersp = 1)
  circos.clear()
}

circos_gene_plot <- function(data,source.use = NULL,bottom_value,right_value,num,title_value,signaling =  NULL ){
  
  if(dev.cur() == 1) {
    dev.new()  # 如果没有图形设备，创建一个
  }
  
  set.seed(123)
  
  if (is.null(signaling) && is.null(source.use)) {
    stop("错误：必须指定 'signaling' 或 'source.use' 参数之一")
  }
  
  # 设置图形参数
  par(mar = c(2, 2, 3, 10))
  if (is.null(signaling)){
    num_data <- subset(merged_data, source == source.use ) %>% 
      group_by(target) %>%
      count(target) %>%
      summarise(count=n()) %>%
      arrange(desc(count)) %>%
      slice_head(n = num) %>%
      pull(target)
    
    merged_data_sub <- subset(merged_data, source == source.use & target %in% num_data)
  } else{
    merged_data_sub <- subset(merged_data, pathway_name == signaling )
  }
  
  link_colors <- ifelse(merged_data_sub$prob_diff > 0, "#E41A1C", "#377EB8")
  names(link_colors) <-as.character(merged_data_sub$source)
  
  cell_types <- unique(c(merged_data_sub$source, merged_data_sub$target))
  df <- merged_data_sub
  
  # 给数据添加ID
  df$id <- 1:nrow(df)
  
  # 处理配体标签（CellChat风格：添加空格）
  ligand.uni <- unique(df$ligand)
  for (i in 1:length(ligand.uni)) {
    df.i <- df[df$ligand == ligand.uni[i], ]
    source.uni <- unique(df.i$source)
    for (j in 1:length(source.uni)) {
      df.i.j <- df.i[df.i$source == source.uni[j], ]
      df.i.j$ligand <- paste0(df.i.j$ligand, paste(rep(" ", j - 1), collapse = ""))
      df$ligand[df$id %in% df.i.j$id] <- df.i.j$ligand
    }
  }
  
  # 处理受体标签
  receptor.uni <- unique(df$receptor)
  for (i in 1:length(receptor.uni)) {
    df.i <- df[df$receptor == receptor.uni[i], ]
    target.uni <- unique(df.i$target)
    for (j in 1:length(target.uni)) {
      df.i.j <- df.i[df.i$target == target.uni[j], ]
      df.i.j$receptor <- paste0(df.i.j$receptor, paste(rep(" ", j - 1), collapse = ""))
      df$receptor[df$id %in% df.i.j$id] <- df.i.j$receptor
    }
  }
  
  # 排序
  df.ordered.source <- df[with(df, order(source, -prob_diff)), ]
  df.ordered.target <- df[with(df, order(target, -prob_diff)), ]
  
  # 确定扇区顺序
  order.source <- unique(df.ordered.source[, c("ligand", "source")])
  order.target <- unique(df.ordered.target[, c("receptor", "target")])
  order.sector <- c(order.source$ligand, order.target$receptor)
  
  # 设置细胞类型颜色
  if (!exists("scPalette")) {
    scPalette <- function(n) {
      library(RColorBrewer)
      if (n <= 8) {
        return(brewer.pal(max(3, n), "Set2"))
      } else if (n <= 12) {
        return(brewer.pal(n, "Set3"))
      } else {
        return(rainbow(n))
      }
    }
  }
  
  color.use <- scPalette(length(cell_types))
  names(color.use) <- cell_types
  
  # ========== 关键修正：扇形颜色设置 ==========
  # 1. 配体扇形颜色 = 来源细胞颜色
  grid.col.ligand <- color.use[as.character(order.source$source)]
  names(grid.col.ligand) <- order.source$ligand  # 用配体名称命名
  
  # 2. 受体扇形颜色 = 目标细胞颜色
  grid.col.receptor <- color.use[as.character(order.target$target)]
  names(grid.col.receptor) <- order.target$receptor  # 用受体名称命名
  
  # 3. 合并所有扇形颜色
  grid.col <- c(grid.col.ligand, grid.col.receptor)
  
  # 准备绘图数据
  df.plot <- df.ordered.source[, c("ligand", "receptor", "prob_diff")]
  names(df.plot) <- c("ligand", "receptor", "value")
  
  # 清除之前的布局
  circos.clear()
  
  # 设置circos参数
  circos.par(
    start.degree = 90,
    gap.degree = 1,
    track.margin = c(0.01, 0.05),
    points.overflow.warning = FALSE
  )
  
  # 绘制和弦图
  chordDiagram(df.plot, order = order.sector, col = link_colors, 
               grid.col = grid.col, transparency = 0.5, link.border = NA, 
               directional = 1, direction.type = c("diffHeight", "arrows"), 
               link.arr.type = "big.arrow", annotationTrack = "grid", 
               annotationTrackHeight = c(0.03), 
               preAllocateTracks = list(track.height = 0.1),  # 固定值
               small.gap = 1, big.gap = 10, link.visible = TRUE, 
               scale = FALSE, link.target.prop = TRUE, 
               reduce = 0)
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.name = get.cell.meta.data("sector.index")
    # 去除空格显示原始基因名
    gene_name <- sub("\\s+$", "", sector.name)
    circos.text(CELL_META$xcenter, 
                CELL_META$ylim[1], 
                gene_name,
                facing = "clockwise",
                niceFacing = TRUE,
                adj = c(0, 0.5),
                cex = 0.4)
  }, bg.border = NA)
  
  # 添加图例
  legend("right",
         legend = names(color.use),
         fill = color.use,
         border = NA,
         bty = "n",
         xpd = TRUE,
         inset = c(right_value, 0),
         cex = 0.8)
  
  legend("bottom",
         legend = c("increase", "decrease"),
         fill = c("#E41A1C", "#377EB8"),
         border = "black",
         title = "",
         cex = 0.7,
         xpd = TRUE,
         inset = c(bottom_value,0),  # 负值向下移动
         bty = "n",
         horiz = TRUE,
         x.intersp = 0.5,
         y.intersp = 1)
  
  if (is.null(signaling)){
    title(main = paste(source.use,"to major Cells Communication"),
         line = title_value,      # 适中距离
         cex = 0.15,
         col.main = "black")
  } else{
    title(main = paste(signaling,"to Cells - Cells Communication"),
          line = title_value,      # 适中距离
          cex = 0.15,
          col.main = "black")
  }
}

barplot_celltype <- function(cellchat.pos,cellchat.other,celltype_label){
  pos.net.source <- subsetCommunication(cellchat.pos) %>% 
    filter(source == celltype_label) %>%
    count(target) %>% 
    rename(celltype = target, ncount = n) %>%
    mutate(group = "pos") %>%
    mutate(class = "source")
  
  other.net.source <- subsetCommunication(cellchat.other)  %>% 
    filter(source == celltype_label) %>%
    count(target) %>% 
    rename(celltype = target, ncount = n) %>%
    mutate(group = "other") %>%
    mutate(class = "source")
  
  pos.net.target <- subsetCommunication(cellchat.pos) %>% 
    filter(target == celltype_label) %>%
    count(source) %>% 
    rename(celltype = source, ncount = n) %>%
    mutate(group = "pos") %>%
    mutate(class = "target")
  
  other.net.target <- subsetCommunication(cellchat.other)  %>% 
    filter(target == celltype_label) %>%
    count(source) %>% 
    rename(celltype = source, ncount = n) %>%
    mutate(group = "other") %>%
    mutate(class = "target")
  
  net <- rbind(pos.net.source,other.net.source,pos.net.target,other.net.target)
  
  
  p <- ggplot(net, aes(x = celltype, y = ncount, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = ncount), 
              position = position_dodge(width = 0.8), 
              vjust = -0.3, size = 3, fontface = "bold") +  # 加粗字体
    facet_wrap(~ class, scales = "free_y") +
    scale_fill_manual(values = c("other" = "#1f77b4", "pos" = "#ff7f0e")) +
    # 添加以下美化部分
    labs(
      title = paste0(celltype_label),
      x = "celltype",
      y = "count",
      fill = "group"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "top"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  
  
  return(p)
}

#circos_gene_plot(source.use = i,title_value = -1,bottom_value = 0,right_value = -0.3,num = 10,data = merged_data)
#circos_celltype_plot(data=merged_data,source.use = "B",right_value2=-0.5,right_value1 =-0.15,bottom_value=0.2,title_value = -7.5)
#circos_gene_plot(signaling = "ADGRE",title_value = -1,bottom_value = 0,right_value = -0.3,num = 10,data = merged_data)
#p <- barplot_celltype(cellchat.pos = cellchat.pos,cellchat.other = cellchat.other,celltype_label = "Tumor")


cellchat_filtered <- function(cellchat.1,cellchat.2){
  celltypes_to_keep <- intersect(levels(cellchat.1@idents), levels(cellchat.2@idents))
  keep_idx <- which(cellchat.1@idents %in% celltypes_to_keep)
  remove_cell_idx <- which(!levels(cellchat.1@idents) %in% celltypes_to_keep)
  
  cellchat.filtered <- cellchat.1
  cellchat.filtered@data.signaling <- cellchat.filtered@data.signaling[,keep_idx]
  cellchat.filtered@data <- cellchat.filtered@data[,keep_idx]
  cellchat.filtered@idents <- cellchat.filtered@idents[keep_idx]
  cellchat.filtered@idents <- factor(cellchat.filtered@idents,levels = celltypes_to_keep)
  cellchat.filtered@meta <- cellchat.pos@meta[keep_idx, ]
  
  cellchat.filtered@images$coordinates <- cellchat.filtered@images$coordinates[keep_idx,]                            
  cellchat.filtered@images$distance <- cellchat.filtered@images$distance[celltypes_to_keep,celltypes_to_keep]                                             
  cellchat.filtered@net$prob <- cellchat.filtered@net$prob[celltypes_to_keep, celltypes_to_keep,]
  cellchat.filtered@net$pval <- cellchat.filtered@net$pval[celltypes_to_keep,celltypes_to_keep,]
  cellchat.filtered@net$count <- cellchat.filtered@net$count[celltypes_to_keep,celltypes_to_keep] 
  cellchat.filtered@net$weight <-  cellchat.filtered@net$weight[celltypes_to_keep,celltypes_to_keep] 
  cellchat.filtered@var.features$features.info <- subset(cellchat.filtered@var.features$features.info,clusters %in% celltypes_to_keep)
  if (length(remove_cell_idx) != 0){
    cellchat.filtered@netP$prob <- cellchat.filtered@netP$prob[-remove_cell_idx, -remove_cell_idx,]
    pathway_list <- names(cellchat.filtered@netP$centr)
    cellchat.filtered@netP$centr[pathway_list] <- 
      map(cellchat.filtered@netP$centr[pathway_list], function(centr_data) {
        vector_attrs <- names(centr_data)[sapply(centr_data, is.vector)]
        for(attr in vector_attrs) {
          centr_data[[attr]] <- centr_data[[attr]][-remove_cell_idx]
        }
        return(centr_data)
      })
  }
  return(cellchat.filtered)
}
