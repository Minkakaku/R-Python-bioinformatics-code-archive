rm(list = ls())

library(cowplot)
library(Seurat)
library(RColorBrewer)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)   # 新增，用于热图
library(igraph)           # 新增，用于网络图

options(stringsAsFactors = FALSE)

outdir <- "/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/core/make_band/cellchat/in/test"
inrds  <- "/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/core/make_band/cellchat/in/test/adata_in.rds"
species <- 'human'

setwd(outdir)

# 创建输出子目录用于存放差异分析结果
out_dir_cellchat <- file.path(outdir, "CellChat_diff_analysis")
dir.create(out_dir_cellchat, showWarnings = FALSE)

# ------------------------- 1. 自定义可视化函数（来自新脚本） -------------------------

# 1.1 差异相互作用网络图
netVisual_diffInteraction <- function(object, comparison = c(1,2), 
                                      measure = c("count", "weight", "count.merged", "weight.merged"), 
                                      color.use = NULL, 
                                      color.edge = c('#b2182b','#2166ac'), 
                                      title.name = NULL, 
                                      sources.use = NULL, 
                                      targets.use = NULL, 
                                      remove.isolate = FALSE, 
                                      top = 1,
                                      weight.scale = FALSE, 
                                      vertex.weight = 20, 
                                      vertex.weight.max = NULL, 
                                      vertex.size.max = 15, 
                                      vertex.label.cex=1,
                                      vertex.label.color= "black",
                                      edge.weight.max = NULL, 
                                      edge.width.max=8, 
                                      alpha.edge = 0.6, 
                                      label.edge = FALSE,
                                      edge.label.color='black',
                                      edge.label.cex=0.8,
                                      edge.curved=0.2,
                                      shape='circle',
                                      layout=in_circle(), 
                                      margin=0.2,
                                      arrow.width=1,
                                      arrow.size = 0.2
){
  options(warn = -1)
  measure <- match.arg(measure)
  
  obj1_raw <- object@net[[comparison[1]]][[measure]]
  obj2_raw <- object@net[[comparison[2]]][[measure]]
  
  # all labels, missing labels set to 0
  add_data_square <- function(dat, all_labels){
    diff_label <- all_labels %>% setdiff(rownames(dat))
    if(length(diff_label) ==0) return(dat)
    mat_row <- matrix(0, nrow = length(diff_label), ncol = ncol(dat),
                      dimnames = list(diff_label, colnames(dat)))
    dat2 <- rbind(dat, mat_row)
    mat_col <- matrix(0, nrow = nrow(dat2), ncol = length(diff_label),
                      dimnames = list(rownames(dat2), diff_label))
    dat3 <- cbind(dat2, mat_col)
    dat3 <- dat3[all_labels, all_labels]
    dat3
  }
  
  all_label <- union(rownames(obj1_raw), rownames(obj2_raw))
  obj1 <- obj1_raw %>% add_data_square(all_label)
  obj2 <- obj2_raw %>% add_data_square(all_label)
  
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) title.name = "Differential number of interactions"
  } else {
    if (is.null(title.name)) title.name = "Differential interaction strength"
  }
  net <- net.diff
  
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) sources.use <- rownames(net.diff)[sources.use]
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) targets.use <- rownames(net.diff)[targets.use]
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  
  net[abs(net) < stats::quantile(abs(net), probs = 1-top, na.rm= T)] <- 0
  
  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
  
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$weight, digits = 1)
  }
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1],color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }
  
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle = rep(0, nrow(edge.start))
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica")
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}

# 1.2 差异相互作用热图
netVisual_heatmap <- function(object, comparison = c(1,2), measure = c("count", "weight"), signaling = NULL, slot.name = c("netP", "net"), color.use = NULL, color.heatmap = c("#2166ac","#b2182b"),
                              title.name = NULL, width = NULL, height = NULL, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                              sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, row.show = NULL, col.show = NULL){
  
  require(ComplexHeatmap)
  
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1_raw <- object@net[[comparison[1]]][[measure]]
    obj2_raw <- object@net[[comparison[2]]][[measure]]
    
    add_data_square <- function(dat, all_labels){
      diff_label <- all_labels %>% setdiff(rownames(dat))
      if(length(diff_label) ==0) return(dat)
      mat_row <- matrix(0, nrow = length(diff_label), ncol = ncol(dat),
                        dimnames = list(diff_label, colnames(dat)))
      dat2 <- rbind(dat, mat_row)
      mat_col <- matrix(0, nrow = nrow(dat2), ncol = length(diff_label),
                        dimnames = list(rownames(dat2), diff_label))
      dat3 <- cbind(dat2, mat_col)
      dat3 <- dat3[all_labels, all_labels]
      dat3
    }
    
    all_label <- union(rownames(obj1_raw), rownames(obj2_raw))
    obj1 <- obj1_raw %>% add_data_square(all_label)
    obj2 <- obj2_raw %>% add_data_square(all_label)
    
    net.diff <- obj2 - obj1
    
    if (measure == "count") {
      if (is.null(title.name)) title.name = "Differential number of interactions"
    } else if (measure == "weight") {
      if (is.null(title.name)) title.name = "Differential interaction strength"
    }
    legend.name = "Relative values"
  } else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[,,signaling]
      if (is.null(title.name)) title.name = paste0(signaling, " signaling network")
      legend.name <- "Communication Prob."
    } else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) title.name = "Number of interactions"
      } else if (measure == "weight") {
        if (is.null(title.name)) title.name = "Interaction strength"
      }
      legend.name <- title.name
    }
  }
  
  net <- net.diff
  
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) sources.use <- rownames(net.diff)[sources.use]
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) targets.use <- rownames(net.diff)[targets.use]
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }
  
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[ ,col.show]
    color.use <- color.use[col.show]
  }
  
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp2(c(min(mat), 0, max(mat)), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = 2), 0, round(max(mat, na.rm = T), digits = 2))
  } else {
    if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp2(c(min(mat), max(mat)), color.heatmap)
    } else {
      color.heatmap.use = colorRamp2(seq(min(mat), max(mat), length.out = 100), 
                                     colorRampPalette(brewer.pal(9, color.heatmap[1]))(100))
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = 2), round(max(mat, na.rm = T), digits = 2))
  }
  
  df <- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE, gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE, gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  } else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
                bottom_annotation = col_annotation, left_annotation = row_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows, cluster_columns = cluster.rows,
                row_names_side = "left", row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), column_names_gp = gpar(fontsize = font.size),
                column_title = title.name, column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90,
                row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = font.size.title), row_title_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"), title_position = "leftcenter-rot",
                                            border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, "mm")))
  return(ht1)
}

# 1.3 修正版的netAnalysis_signalingRole_heatmap（解决原函数bug）
netAnalysis_signalingRole_heatmap_fixed <- function (object, signaling = NULL, 
                                                     pattern = c("outgoing", "incoming", "all"), 
                                                     slot.name = "netP", color.use = NULL, 
                                                     color.heatmap = "BuGn", title = NULL, width = 10, height = 8, 
                                                     font.size = 8, font.size.title = 10, cluster.rows = FALSE, 
                                                     cluster.cols = FALSE) 
{
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
  dimnames(outgoing) <- list(levels(object@idents), names(centr))
  dimnames(incoming) <- dimnames(outgoing)
  for (i in 1:length(centr)) {
    outgoing[, i] <- centr[[i]]$outdeg
    incoming[, i] <- centr[[i]]$indeg
  }
  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  }
  else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  }
  else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  }
  else {
    title <- paste0(paste0(legend.name, " signaling patterns"), " - ", title)
  }
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
  mat[mat == 0] <- NA
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap))))(99)
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                         length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                      show_annotation_name = FALSE)
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  }
  else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                      round(max(mat, na.rm = T), digits = 1))
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = "Relative strength", bottom_annotation = col_annotation, 
                top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), width = unit(width, 
                                                                           "cm"), height = unit(height, "cm"), column_title = title, 
                column_title_gp = gpar(fontsize = font.size.title), 
                column_names_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                                                   fontface = "plain"), title_position = "leftcenter-rot", 
                                                                   border = NA, at = legend.break, legend_height = unit(20, 
                                                                                                                        "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, 
                                                                                                                                                                                 "mm")))
  return(ht1)
}

# ------------------------- 2. 原脚本流程：读取数据、拆分、构建CellChat对象 -------------------------

# 读入 Seurat 对象（已包含 counts 矩阵）
Seuset <- readRDS(inrds)

# 按分组拆分
cellchat.list <- SplitObject(Seuset, split.by = "group_final")
meta.list <- lapply(cellchat.list, function(x) x@meta.data)

# 分组名称（必须和 split.by 一致）
group <- as.character(unique(Seuset$group_final))

# 创建 CellChat 对象：先标准化，再传入（CellChat 自动取 data 矩阵）
cellchat.objet.list <- lapply(group, FUN = function(x) {
  sub_seu <- NormalizeData(cellchat.list[[x]], 
                           normalization.method = "LogNormalize",
                           scale.factor = 10000, 
                           verbose = FALSE)
  createCellChat(sub_seu, meta = meta.list[[x]], group.by = "celltype")
})
names(cellchat.objet.list) <- group

# 运行 CellChat 分析的自定义函数
runCellchat <- function(cellchat, species = "human", use.ident = "celltype", 
                        filter.cells = 10, name_ct) {
  print(name_ct)
  cellchat <- setIdent(cellchat, ident.use = use.ident)
  
  if (species == "human") {
    CellChatDB <- CellChatDB.human  
  } else {
    CellChatDB <- CellChatDB.mouse  
  }
  
  write.table(CellChatDB$interaction, 'cellchatDB.xls', sep = "\t")
  showDatabaseCategory(CellChatDB)
  
  cellchat@DB <- subsetDB(CellChatDB, 
                          search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")) 
  cellchat <- subsetData(cellchat)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat@idents <- droplevels(cellchat@idents, 
                                exclude = setdiff(levels(cellchat@idents), 
                                                  unique(cellchat@idents)))
  
  cellchat <- computeCommunProb(cellchat, type = "triMean", trim = 0)
  cellchat <- filterCommunication(cellchat, min.cells = filter.cells)
  cellchat <- computeCommunProbPathway(cellchat)
  
  df.net <- subsetCommunication(cellchat)
  relation <- paste0(df.net[,1], "->", df.net[,2])
  df.net <- cbind(relation, df.net)
  colnames(df.net)[1] <- "relation"
  write.table(df.net, paste0(name_ct, ".net_lr.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.table(df.netp, paste0(name_ct, ".netPathways_lr.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)
  return(cellchat)
}

# 对每个分组运行 CellChat
cellchat.objet.list.new <- lapply(group, FUN = function(x) {
  runCellchat(cellchat.objet.list[[x]], species = species, name_ct = x)
})
names(cellchat.objet.list.new) <- group

# 检查各分组的细胞类型数是否一致
idents.list <- lapply(cellchat.objet.list.new, function(x) length(levels(x@idents)))
if (length(unique(unlist(idents.list))) > 1) {
  maxgroup <- names(idents.list)[idents.list == max(unlist(idents.list))]
  group.max <- levels(cellchat.objet.list.new[[maxgroup]]@idents)
  cellchat.objet.list.new <- lapply(cellchat.objet.list.new, function(x) {
    liftCellChat(x, group.max)
  })
}

# 合并 CellChat 对象（用于比较）
cellchat <- mergeCellChat(cellchat.objet.list.new, 
                          add.names = names(cellchat.objet.list.new))
celltype <- unique(cellchat@meta$celltype_final) %>% as.character()

save(cellchat, file = file.path(out_dir_cellchat, "cellchat_object.RData"))
save(cellchat.objet.list.new, file = file.path(out_dir_cellchat, "cellchat_object_list.RData"))

# ------------------------- 3. 原脚本中的基本对比图 -------------------------
# 绘图辅助函数（保留原样）
plot_save <- function(filename, plt, dpi = 300, w = 7, h = 7) {
  ggsave(paste0(filename, '.pdf'), plt, height = h, width = w)
}

fig_save <- function(filename, plt, dpi = 300, w = 7, h = 7) {
  pdf(paste0(filename, '.pdf'), height = h, width = w)
  draw(plt)
  dev.off()
}

recor_save <- function(filename, plt, dpi = 300, w = 7, h = 7) {
  pdf(paste0(filename, '.pdf'), height = h, width = w)
  replayPlot(plt)
  dev.off()
}

# 比较两组互作数量与强度
gg1 <- compareInteractions(cellchat, show.legend = FALSE, measure = "count")
plot_save(file.path(out_dir_cellchat, "Number of inferred interactions.barplot"), gg1)

gg2 <- compareInteractions(cellchat, show.legend = FALSE, measure = "weight")
plot_save(file.path(out_dir_cellchat, "Interaction strength.barplot"), gg2)

# ------------------------- 4. 新增：差异网络、热图、气泡图、排名图、中心性热图等 -------------------------

# 获取所有分组对
all_combns_ls <- combn(seq_along(group), 2, simplify = FALSE)

# 4.1 差异网络图（netVisual_diffInteraction）
for (idx in all_combns_ls) {
  grp_pair <- group[idx]
  pdf(file.path(out_dir_cellchat, paste0("diff_network_", grp_pair[1], "_vs_", grp_pair[2], ".pdf")),
      height = 10, width = 10)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, comparison = idx,
                            title.name = "Differential number of interactions",
                            margin = c(0.2, 0.2, 0.5, 0.2), arrow.size = 0.6)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, comparison = idx, measure = "weight",
                            title.name = "Differential interaction strength",
                            margin = c(0.2, 0.2, 0.5, 0.2), arrow.size = 0.6)
  dev.off()
}

# 4.2 差异热图（netVisual_heatmap）
for (idx in all_combns_ls) {
  grp_pair <- group[idx]
  pdf(file.path(out_dir_cellchat, paste0("diff_heatmap_", grp_pair[1], "_vs_", grp_pair[2], ".pdf")),
      height = 10, width = 10)
  draw(netVisual_heatmap(cellchat, comparison = idx))
  draw(netVisual_heatmap(cellchat, comparison = idx, measure = "weight"))
  dev.off()
}

# 4.3 差异气泡图（netVisual_bubble）需要指定部分细胞类型，这里使用全部细胞类型或自定义
# 获取所有细胞类型
all_celltypes <- levels(cellchat@idents$joint)  # 合并后的细胞类型
# 示例：选取前几个作为发送者和接收者，可根据实际修改
sources <- all_celltypes[1:min(5, length(all_celltypes))]
targets <- all_celltypes
if (length(sources) > 0 && length(targets) > 0) {
  for (idx in all_combns_ls) {
    grp_pair <- group[idx]
    pdf(file.path(out_dir_cellchat, paste0("diff_bubble_", grp_pair[1], "_vs_", grp_pair[2], ".pdf")),
        height = 10, width = 12)
    print(netVisual_bubble(cellchat, sources.use = sources, targets.use = targets, 
                           comparison = idx, angle.x = 45))
    dev.off()
  }
}

# 4.4 差异排名图（rankNet）
for (idx in all_combns_ls) {
  grp_pair <- group[idx]
  pdf(file.path(out_dir_cellchat, paste0("rankNet_", grp_pair[1], "_vs_", grp_pair[2], ".pdf")),
      height = 8, width = 12)
  p1 <- rankNet(cellchat, mode = "comparison", measure = "count", comparison = idx,
                stacked = TRUE, do.stat = TRUE) + ggtitle("Relative information flow - count")
  p2 <- rankNet(cellchat, mode = "comparison", measure = "weight", comparison = idx,
                stacked = TRUE, do.stat = TRUE) + ggtitle("Relative information flow - weight")
  p3 <- rankNet(cellchat, mode = "comparison", measure = "count", comparison = idx,
                stacked = FALSE, do.stat = TRUE) + ggtitle("Information flow - count")
  p4 <- rankNet(cellchat, mode = "comparison", measure = "weight", comparison = idx,
                stacked = FALSE, do.stat = TRUE) + ggtitle("Information flow - weight")
  print(p1 + p2 + p3 + p4 + plot_layout(ncol = 2))
  dev.off()
}

# 4.5 各分组内信号对细胞群贡献的热图（outgoing/incoming/all）
# 注意：需要确保cellchat.objet.list.new中每个对象都已计算中心性（runCellchat中已包含）
for (i in seq_along(cellchat.objet.list.new)) {
  name <- names(cellchat.objet.list.new)[i]
  obj <- cellchat.objet.list.new[[i]]
  pathways <- obj@netP$pathways
  if (length(pathways) == 0) next
  
  # outgoing
  p_out <- netAnalysis_signalingRole_heatmap_fixed(obj, pattern = "outgoing", signaling = pathways,
                                                   title = name, width = 15, height = 15)
  fig_save(file.path(out_dir_cellchat, paste0(name, "_signaling_role_outgoing")), p_out, w = 15, h = 15)
  # incoming
  p_in <- netAnalysis_signalingRole_heatmap_fixed(obj, pattern = "incoming", signaling = pathways,
                                                  title = name, width = 15, height = 15)
  fig_save(file.path(out_dir_cellchat, paste0(name, "_signaling_role_incoming")), p_in, w = 15, h = 15)
  # all
  p_all <- netAnalysis_signalingRole_heatmap_fixed(obj, pattern = "all", signaling = pathways,
                                                   title = name, width = 15, height = 15)
  fig_save(file.path(out_dir_cellchat, paste0(name, "_signaling_role_all")), p_all, w = 15, h = 15)
}

# 4.6 差异交互强度散点图（netAnalysis_signalingRole_scatter）
num.link <- sapply(cellchat.objet.list.new, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
scatter_plots <- lapply(seq_along(cellchat.objet.list.new), function(i) {
  netAnalysis_signalingRole_scatter(cellchat.objet.list.new[[i]], 
                                     title = names(cellchat.objet.list.new)[i],
                                     weight.MinMax = range(num.link))
})
# 拼合散点图
combined_scatter <- wrap_plots(scatter_plots, ncol = 2)
plot_save(file.path(out_dir_cellchat, "diff_interaction_strength_scatter"), combined_scatter, w = 12, h = 8)

message("所有差异分析可视化已完成，结果保存在：", out_dir_cellchat)
