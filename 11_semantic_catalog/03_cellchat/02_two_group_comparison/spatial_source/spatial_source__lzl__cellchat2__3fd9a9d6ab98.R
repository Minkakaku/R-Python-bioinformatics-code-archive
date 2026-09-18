library(tidyverse)
library(Seurat)
library(CellChat)

# define function 
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
  
  # shared_label = intersect(rownames(obj1_raw), rownames(obj2_raw))
  # obj1 = obj1_raw[shared_label, shared_label]
  # obj2 = obj2_raw[shared_label, shared_label]
  
  
  # all labels, missing labels set to 0
  add_data_square <- function(dat, all_labels){
    diff_label <- all_labels %>% setdiff(rownames(dat))
    
    if(length(diff_label) ==0) return(dat)
    
    # row
    mat_row <- matrix(
      0, 
      nrow = length(diff_label),
      ncol = ncol(dat),
      dimnames = list(
        diff_label,
        colnames(dat)
      )
    )
    dat2 <- rbind(dat, mat_row)
    
    # column
    mat_col <- matrix(
      0, 
      nrow = nrow(dat2),
      ncol = length(diff_label),
      dimnames = list(
        rownames(dat2),
        diff_label
      )
    )
    dat3 <- cbind(dat2, mat_col)
    
    dat3 <- dat3[all_labels, all_labels]
    dat3
  }
  
  all_label <- union(rownames(obj1_raw), rownames(obj2_raw))
  obj1 <- obj1_raw %>% add_data_square(all_label)
  obj2 <- obj2_raw %>% add_data_square(all_label)
  
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  } else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
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
  
  # ERROR: Length of new attribute value must be 1 or 45
  # igraph::V(g)$size<-vertex.weight
  
  # Fix below
  igraph::V(g)$size<-vertex.weight
  
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  #igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1],color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, alpha.edge)
  
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }
  
  
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    # add this line
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
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}


netVisual_heatmap <- function(object, comparison = c(1,2), measure = c("count", "weight"), signaling = NULL, slot.name = c("netP", "net"), color.use = NULL, color.heatmap = c("#2166ac","#b2182b"),
                              title.name = NULL, width = NULL, height = NULL, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                              sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, row.show = NULL, col.show = NULL){
  
  require(ComplexHeatmap)
  
  # obj1 <- object.list[[comparison[1]]]
  # obj2 <- object.list[[comparison[2]]]
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    # obj1 <- object@net[[comparison[1]]][[measure]]
    # obj2 <- object@net[[comparison[2]]][[measure]]
    
    obj1_raw <- object@net[[comparison[1]]][[measure]]
    obj2_raw <- object@net[[comparison[2]]][[measure]]
    
    # shared_label = intersect(rownames(obj1_raw), rownames(obj2_raw))
    # obj1 = obj1_raw[shared_label, shared_label]
    # obj2 = obj2_raw[shared_label, shared_label]
    
    # all labels, missing labels set to 0
    add_data_square <- function(dat, all_labels){
      diff_label <- all_labels %>% setdiff(rownames(dat))
      
      if(length(diff_label) ==0) return(dat)
      
      # row
      mat_row <- matrix(
        0, 
        nrow = length(diff_label),
        ncol = ncol(dat),
        dimnames = list(
          diff_label,
          colnames(dat)
        )
      )
      dat2 <- rbind(dat, mat_row)
      
      # column
      mat_col <- matrix(
        0, 
        nrow = nrow(dat2),
        ncol = length(diff_label),
        dimnames = list(
          rownames(dat2),
          diff_label
        )
      )
      dat3 <- cbind(dat2, mat_col)
      
      dat3 <- dat3[all_labels, all_labels]
      dat3
    }
    
    all_label <- union(rownames(obj1_raw), rownames(obj2_raw))
    obj1 <- obj1_raw %>% add_data_square(all_label)
    obj2 <- obj2_raw %>% add_data_square(all_label)
    
    
    net.diff <- obj2 - obj1
    
    if (measure == "count") {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    } else if (measure == "weight") {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  } else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[,,signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    } else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      } else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  
  net <- net.diff
  
  
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
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
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), 0, round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
    # color.heatmap.use = colorRamp3(c(seq(min(mat), -(max(mat)-min(max(mat)))/9, length.out = 4), 0, seq((max(mat)-min(max(mat)))/9, max(mat), length.out = 4)), RColorBrewer::brewer.pal(n = 9, name = color.heatmap))
  } else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), color.heatmap)
    } else if (length(color.heatmap) == 1) {
      color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
  }
  # col_fun(as.vector(mat))
  
  df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))
  
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
  
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  } else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
                bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
                cluster_rows = cluster.rows,cluster_columns = cluster.rows,
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                # width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, #at = colorbar.break,
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
  )
  #  draw(ht1)
  return(ht1)
}

# define default-------------
out_dir_cellchat <- "cellchat_dir2" 
dir.create(out_dir_cellchat)

# read data  -------
all_data <- readRDS('/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/data/all_data.rds')
table(all_data$annotation)
Idents(all_data) <- "annotation"
all_data <- JoinLayers(all_data)

Macrophage_meta <- Macrophage_meta[,-42]
Tcell_meta


#all_data_cn_meta <- readRDS("/mnt/NexData_file/home/lzl/20260515_XE2510316140-03/CN_analysis/all_data_meta_p_cn.rds")

all_data_immnue_meta <- rbind(Macrophage_meta,Tcell_meta)

length(unique(all_data_immnue_meta$cell_id_sample))

table(all_data_immnue_meta$subtype)
immnue_idx <- all_data_immnue_meta$cell_id_sample
# only used 50,000 cells
set.seed(12345)
subset_cell_idx <- 
  all_data_meta_all %>%
  rownames_to_column("cell_id") %>%
  # group_by(CN2, group, annotation2) %>%
  slice_sample(prop = signif(50000/nrow(all_data_meta_all)))

all_data_meta <- all_data_meta_all[subset_cell_idx$cell_id, ]
all_data_norm <- all_data_norm_all[, subset_cell_idx$cell_id]

###
all_data_meta_all     <- all_data[[]] %>% .[immnue_idx, ]

all_data_meta_all <- merge(
  all_data_meta_all,
  all_data_immnue_meta %>% dplyr::select(cell_id_sample, subtype),  # 只传需要的列
  by = "cell_id_sample",
  all.x = TRUE
)
#all_data_meta_all <- na.omit(all_data_meta_all)

all_data_norm_all     <- LayerData(all_data, 'data') %>% .[, immnue_idx]
rownames(all_data_meta_all) <- all_data_meta_all$cell_id_sample
all_data_meta_all <- all_data_meta_all[colnames(all_data_norm_all), ]

stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all))) 

all_data_meta <- all_data_meta_all
all_data_norm <- all_data_norm_all

all_data$group_3 %>% table()

all_grp_ls <-
  all_data_meta %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, group_3) %>%
  group_by(group_3) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.name = FALSE)

names(all_grp_ls)
# using CellChat db --------
all_cellchat_ls <-
  all_grp_ls %>%
  map(function(cells){ # browser()
    
    meta <- all_data_meta[cells, ] 
    dat  <- all_data_norm[, cells]
    
    # table
    cells <- all_data_meta_all$subtype %>% unique()
    if(!all(cells %in% meta$subtype)){
      missing_cells <- cells %>% setdiff(unique(meta$subtype))
      missing_meta <-
        all_data_meta_all %>%
        rownames_to_column('cell_id') %>%
        dplyr::filter(subtype %in% missing_cells) %>%
        group_by(subtype) %>%
        slice_head(n = 1) %>%
        column_to_rownames("cell_id")
      meta <- rbind(meta, missing_meta)
      print(missing_meta)
      dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta)])
    }
    
    meta$samples <- factor(meta$group)
    
    cellchat <- createCellChat(object = dat, meta = meta, group.by = "subtype")
    
    # cellchat db
    CellChatDB.use <- 
      CellChatDB.human %>% 
      subsetDB(
        search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), 
        key    = "annotation"
      )
    
    cellchat@DB <- CellChatDB.use
    
    
    # preprocess
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # optinal - project PPI data
    # cellchat = projectData(cellchat, PPI.human)
    cellchat <- smoothData(cellchat, adj = PPI.human)
    
    # Compute the communication probability
    cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, population.size = TRUE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    # Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    
    # aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    # cellchat <- netAnalysis_computeCentrality(cellchat)
    
    cellchat
  })


saveRDS(all_cellchat_ls, file.path(out_dir_cellchat, "all_cellchat_ls.rds"))

# export ligand receptor table
all_cellchat_ls  %>%
  iwalk(function(cellchat, name){ 
    df.net <- tryCatch(
      cellchat %>% subsetCommunication() %>% as.data.frame(),
      error = function(e) NULL
    )
    if(!is.null(df.net)){
      df.net$relation <-  paste0(df.net[[1]], "->", df.net[[2]])
      write.table(
        df.net,
        file      = file.path(out_dir_cellchat, str_glue("net_lr - {name}.txt")),
        sep       = "\t",
        quote     = FALSE,
        row.names = FALSE,
        col.names = TRUE
      )
    }
    
    # network
    df.netp <- tryCatch(
      subsetCommunication(cellchat, slot.name="netP"),
      error = function(e) NULL
    )
    if(!is.null(df.netp)){
      write.table(
        df.netp,
        file      = file.path(out_dir_cellchat, str_glue("netPathways_lr - {name}.txt")),
        sep       = "\t",
        quote     = FALSE,
        row.names = FALSE,
        col.names = TRUE
      )
    }
  })



# merge cellchat
idents.list <- map(all_cellchat_ls, ~ levels(.x@idents))
idents.list.len <- idents.list %>% map_dbl(length)

if(any(idents.list.len[1] != idents.list.len)){
  
  idents.max      <- idents.list[[which.max(idents.list.len)]]
  all_cellchat_ls <- all_cellchat_ls %>% map(function(x) liftCellChat(x, idents.max))
}



all_cellchat_new_ls <-
  all_cellchat_ls %>%
  imap(function(cellchat, name){
    cellchat <- tryCatch(netAnalysis_computeCentrality(cellchat), error = function(e) NULL)
    if(is.null(cellchat)){
      cat(str_glue("group {name} cannot run netAnalysis_computeCentrality and removed!\n\n"))
    }
    cellchat
  }) %>% 
  discard(is.null)  

groups   <- names(all_cellchat_new_ls)
cellchat <- mergeCellChat(all_cellchat_new_ls, add.names = groups, cell.prefix = TRUE)


saveRDS(all_cellchat_new_ls, file.path(out_dir_cellchat, "all_cellchat_new_ls.rds")) 
saveRDS(cellchat, file.path(out_dir_cellchat, "cellchat.rds")) 


all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)


## comprare interation ------
pdf(
  file.path(out_dir_cellchat, str_glue("comprare interation plot.pdf")),
  height = 4,
  width  = 6
)
compareInteractions(cellchat, show.legend = FALSE, measure = "count")
compareInteractions(cellchat, show.legend = FALSE, measure = "weight")
dev.off()


## network plot ------------
# #统计细胞和细胞之间通信的数量（有多少个配体-受体对）和强度（概率）
# cellchat <- aggregateNet(cellchat)
all_cellchat_new_ls$C@netP$pathways
all_cellchat_new_ls %>%
  iwalk(function(cellchat, name){
    pdf(
      file.path(out_dir_cellchat, str_glue("network plot - {name}.pdf")),
      height = 8,
      width  = 16
    )
    
    #计算每种细胞各有多少个
    groupSize <- as.numeric(table(cellchat@idents))
    
    par(mfrow = c(1, 2), xpd=TRUE)
    netVisual_circle(
      cellchat@net$count,
      vertex.weight = groupSize,
      weight.scale = T, 
      label.edge= F, 
      title.name = "Number of interactions"
    )
    netVisual_circle(
      cellchat@net$weight, 
      vertex.weight = groupSize,
      weight.scale = T, 
      label.edge= F,
      title.name = "Interaction weights/strength"
    )
    
    dev.off()
  })


unique(cellchat@meta$datasets)
## diff network plot --------
all_combns_ls %>%
  walk(function(x){
    grp <- groups[x]
    pdf(
      file.path(out_dir_cellchat, str_glue("diff network plot - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
      height = 10,
      width  = 10
    )
    # par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(
      cellchat,
      weight.scale  = T, 
      comparison    = x,
      title.name    = "Differential number of interactions",
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    netVisual_diffInteraction(
      cellchat,
      weight.scale  = T,
      comparison    = x,
      measure       = "weight", 
      title.name    = "Differential interaction strength",
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    dev.off()
  })


all_combns_ls[[1]]
## diff heatmap plot --------
all_combns_ls %>%
  walk(function(x){
    grp <- groups[x]
    
    pdf(
      file.path(out_dir_cellchat, str_glue("diff heatmap - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
      height = 10,
      width  = 10
    )
    print(netVisual_heatmap(cellchat, comparison = x))                    # count
    print(netVisual_heatmap(cellchat, comparison = x, measure = "weight")) # weight
    dev.off()
  })
## diff rank plot --------
all_combns_ls %>%
  walk(function(x){
    grp <- groups[x]
    pdf(
      file.path(out_dir_cellchat, str_glue("diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
      height = 6,
      width  = 6
    )
    # par(mfrow = c(1,2), xpd=TRUE)
    p1 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "count",
      sources.use = NULL,
      targets.use = NULL, 
      comparison  = x,
      stacked     = T,
      do.stat     = TRUE
    ) +
      ggtitle("Relative information flow - count")
    print(p1)
    
    p2 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "weight",
      sources.use = NULL,
      targets.use = NULL, 
      comparison  = x,
      stacked       = T,
      do.stat       = TRUE
    ) +
      ggtitle("Relative information flow - weight")
    print(p2)
    
    p3 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "count",
      sources.use = NULL,
      targets.use = NULL, 
      comparison  = x,
      stacked     = F,
      do.stat     = TRUE
    ) +
      ggtitle("Information flow - count")
    print(p3)
    
    p4 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "weight",
      sources.use = NULL,
      targets.use = NULL, 
      comparison  = x,
      stacked       = F,
      do.stat       = TRUE
    ) +
      ggtitle("Information flow - weight")
    print(p4)
    
    dev.off()
  })


## diff incoming and outgoing cells --------
library(future)
plan("sequential")
# outgoing
pdf(
  file.path(out_dir_cellchat, str_glue("contribution of signals to all groups - outgoing.pdf")),
  height = 10, 
  width  = 6
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){
    pathways <- cellchat@netP$pathways
    p <- cellchat %>%
      netAnalysis_computeCentrality() %>%
      netAnalysis_signalingRole_heatmap(
        pattern   = "outgoing", 
        signaling = pathways,
        title     = name,
        width     = 10,
        height    = length(pathways)/3
      )
    print(p)
  })
dev.off()

# incoming 
pdf(
  file.path(out_dir_cellchat, str_glue("contribution of signals to all groups - incoming.pdf")),
  height = 10, 
  width  = 6
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){
    pathways <- cellchat@netP$pathways
    p <- cellchat %>%
      netAnalysis_computeCentrality() %>%
      netAnalysis_signalingRole_heatmap(
        pattern   = "incoming", 
        signaling = pathways,
        title     = name,
        width     = 10, 
        height    = length(pathways)/3
      )
    print(p)
  })
dev.off()

# all 
pdf(
  file.path(out_dir_cellchat, str_glue("contribution of signals to all groups - all.pdf")),
  height = 10, 
  width  = 6
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){
    pathways <- cellchat@netP$pathways
    p <- cellchat %>%
      netAnalysis_computeCentrality() %>%
      netAnalysis_signalingRole_heatmap(
        pattern   = "all", 
        signaling = pathways,
        title     = name,
        width     = 10, 
        height    = length(pathways)/3
      )
    print(p)
  })
dev.off()

## diff interaction strength -----

pdf(
  file.path(out_dir_cellchat, str_glue("diff interaction strength - all groups.pdf")),
  height = 5,
  width  = 10
)
num.link <- sapply(all_cellchat_new_ls, function(x) {rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
all_cellchat_new_ls %>%
  imap(function(cellchat, name){
    cellchat %>% 
      netAnalysis_computeCentrality() %>%
      netAnalysis_signalingRole_scatter(
        title = name,
        weight.MinMax = range(num.link)
      )
  }) %>%
  patchwork::wrap_plots()
dev.off()


# NOT RUN.
# using specific cells or moleculars

# pathway network - signaling Role -----
signaling_pathway_ls <- list(
  Tumor = c("APP", "CEACAM", "CLDN"),
  TSA = c("CEACAM", "COLLAGEN")
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){
    
    pdf(
      file.path(out_dir_cellchat, str_glue("pathway network - signaling role - {name}.pdf")),
      height = 6, 
      width  = 6
    )
    p <- cellchat %>%
      netAnalysis_computeCentrality()  %>%
      netAnalysis_signalingRole_network(
        signaling = signaling_pathway_ls[[name]],
        # width = 8, 
        # height = 2.5,
        font.size = 10
      )
    print(p)
    dev.off()
  })



# diff molecule ---------
pdf(
  file.path(out_dir_cellchat, str_glue("diff molecule A vs B.pdf")),
  height = 4,
  width  = 6
)

# some cells to plot
all_cells <- cellchat@idents$joint %>% levels()
c("CN1", "CN2", "CN5", "CN7","CN8") %>%
  #all_cells %>%
  map(function(cell){
    netAnalysis_signalingChanges_scatter(
      cellchat, 
      idents.use = cell
    )
  })

dev.off()


# pathway plot----------------
CellChatDB.human$interaction %>%
  dplyr::filter(
    ligand %in% c("CD80", "CD86", "NECTIN", 'MHC-II') |
      receptor  %in% c("CD80", "CD86", "NECTIN", 'MHC-II')
  )


ecm_remodeling <- c(
  "LAMININ",      # 基底膜结构
  "COLLAGEN",     # 基质沉积/纤维化
  "JAM",          # 细胞间粘附/紧密连接
  "TENASCIN",     # 基质重塑/侵袭
  "PERIOSTIN",    # 纤维化修复
  "MPZ",          # 粘附/髓鞘样
  "GAP",          # 细胞骨架/粘附
  "ICAM",         # 细胞间粘附（双重功能）
  "CDH5",         # 血管内皮粘附
  "PECAM1",       # 内皮细胞粘附
  "SEMA4",        # 轴突导向/细胞迁移
  "PTPRM",        # 酪氨酸磷酸酶/粘附
  "THBS",         # 血栓反应蛋白/基质（双重功能）
  "ADGRE",        # 粘附GPCR
  "HSPG",         # 硫酸乙酰肝素/基质
  "GDF",          # 生长分化/基质
  "TGFB"          # 纤维化主驱动
)

inflammation <- c(
  "CSF",          # 巨噬细胞招募/活化
  "CXCL",         # 趋化因子/免疫招募
  "PROS",         # 前列腺素/炎症
  "GAS",          # 生长停滞/免疫调节
  "CD276",        # 免疫检查点
  "CD99",         # 凋亡/免疫调节
  "GRN",          # 颗粒蛋白/炎症
  "ACTIVIN",      # 激活素/免疫
  "THY1",         # 间充质/炎症
  "KIT",          # 干细胞因子/炎症
  "ANGPT",        # 血管生成（炎症相关）
  "VEGF"          # 血管生成（缺氧/炎症）
)
vertex.receiver <- c("Extravillous Trophoblast")
cell2id_a <- 
  all_cellchat_ls$`T`@net$pro %>% 
  rownames() %>%
  factor(., levels = .) %>%
  {structure(
    as.numeric(.),
    names = as.character(.)
  )}


netVisual_aggregate(
  all_cellchat_ls$`T`,
  signaling = ecm_remodeling[1],  
  #vertex.receiver = all_cells[3],
  layout = "hierarchy"
)

all_cellchat_new_ls$`T`@netP$pathways

netVisual_aggregate(
  object = all_cellchat_new_ls$`C`,
  signaling = "CSF",       # 你要看的通路名
  layout = "circle"       # 三选一：circle / hierarchy / chord
)

pdf(
  file.path(out_dir_cellchat, str_glue("heatmap CXCL.pdf")),
  height = 6,
  width  = 6
)
netVisual_heatmap(
  object = all_cellchat_new_ls$`T`,
  signaling = 'CXCL',  
  color.heatmap = "Reds"
)
dev.off()
## 
netAnalysis_contribution(all_cellchat_ls$A, signaling = pathways.show[2])


netVisual_individual(
  all_cellchat_ls$B, 
  signaling = pathways.show[4]
)


## dotplot
# 1
pdf(
  file.path(out_dir_cellchat, "CD4_macrophage sourced signaling molecular.pdf"),
  height = 6,
  width = 10
)
all_patheways <- 
  cellchat@netP %>%
  map('pathways') %>%
  unname() %>% 
  do.call(intersect, .)
all_patheways %>%
  iwalk(
    ~ tryCatch(
      netVisual_bubble(
        cellchat,
        signaling  = .x,
        sources.use = 'CD4T-Macrophage',
        comparison  = c(1,2),
        title.name  = .x,
        remove.isolate = T
      ) %>%
        print(),
      error = function(e) {}
    )
  )

dev.off()

# 2
pdf(
  file.path(out_dir_cellchat, "CD4_macrophage targeted signaling molecular.pdf"),
  height = 6,
  width = 10
)
all_patheways <- 
  cellchat@netP %>%
  map('pathways') %>%
  unname() %>% 
  do.call(intersect, .)
all_patheways %>%
  iwalk(
    ~ tryCatch(
      netVisual_bubble(
        cellchat,
        signaling  = .x,
        targets.use = 'CD4T-Macrophage',
        comparison  = c(1,2),
        title.name  = .x,
        remove.isolate = T
      ) %>%
        print(),
      error = function(e) {}
    )
  )

dev.off()


# 3
pdf(
  file.path(out_dir_cellchat, "Epithelial sourced signaling molecular.pdf"),
  height = 6,
  width = 10
)
all_patheways <- 
  cellchat@netP %>%
  map('pathways') %>%
  unname() %>% 
  do.call(intersect, .)
all_patheways %>%
  iwalk(
    ~ tryCatch(
      netVisual_bubble(
        cellchat,
        signaling  = .x,
        sources.use = 'Epithelial',
        comparison  = c(1,2),
        title.name  = .x,
        remove.isolate = T
      ) %>%
        print(),
      error = function(e) {}
    )
  )

dev.off()

# 4
pdf(
  file.path(out_dir_cellchat, "Epithelial targeted signaling molecular.pdf"),
  height = 6,
  width = 10
)
all_patheways <- 
  cellchat@netP %>%
  map('pathways') %>%
  unname() %>% 
  do.call(intersect, .)
all_patheways %>%
  iwalk(
    ~ tryCatch(
      netVisual_bubble(
        cellchat,
        signaling  = .x,
        targets.use = 'Epithelial',
        comparison  = c(1,2),
        title.name  = .x,
        remove.isolate = T
      ) %>%
        print(),
      error = function(e) {}
    )
  )

dev.off()

plotGeneExpression(all_cellchat_ls$A, signaling = "NECTIN")
plotGeneExpression(all_cellchat_ls$B, signaling = "NECTIN")

plotGeneExpression(all_cellchat_ls$A, signaling = "CXCL")
plotGeneExpression(all_cellchat_ls$B, signaling = "CXCL")

plotGeneExpression(all_cellchat_ls$A, signaling = "PVR")
plotGeneExpression(all_cellchat_ls$B, signaling = "PVR")


plotGeneExpression(all_cellchat_ls$A, signaling = "COLLAGEN")
plotGeneExpression(all_cellchat_ls$B, signaling = "COLLAGEN")

plotGeneExpression(all_cellchat_ls$A, signaling = "MHC-II")
plotGeneExpression(all_cellchat_ls$B, signaling = "MHC-II")

# # chord plot
# all_cellchat_ls %>%
#   iwalk(function(cellchat, name){
#     
#     pdf(
#       file.path(out_dir_cellchat, str_glue("chord plot - {name}.pdf")),
#       height = 10,
#       width  = 10
#     )
#     tryCatch(netVisual_chord_gene(cellchat, sources.use = "Epithelial cells", targets.use = NULL), error = function(x) x)
#     tryCatch(netVisual_chord_gene(cellchat, sources.use = NULL, targets.use = "Epithelial cells"), error = function(x) x)
#     dev.off()
#     
#   })


# gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
# gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
# gg1 + gg2

# using ChellChat and CellPhone db --------
db_full <- readRDS("CellChatDB.human_user.merged.rds")
all_cellchat_ls2 <-
  all_grp_ls %>%
  map(function(cells){# browser()
    
    meta <- all_data_meta[cells, ] 
    dat  <- all_data_norm[, cells]
    
    meta$samples <- factor(meta$project)
    
    cellchat <- createCellChat(object = dat, meta = meta, group.by = "final_annotation")
    
    # cellchat db
    CellChatDB.use <- 
      db_full %>% 
      subsetDB(
        search = db_full$interaction$annotation %>% unique %>% setdiff("Non-protein Signaling"), 
        key    = "annotation"
      )
    
    cellchat@DB <- CellChatDB.use
    
    
    # preprocess
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # Compute the communication probability
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    # Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    
    # aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    
    cellchat
  })


# saveRDS(all_cellchat_ls2, "all_cellchat_ls_all_silde2.rds")

# dir.create("cellchat_result")
all_cellchat_ls2 %>%
  iwalk(function(cellchat, name){
    
    pdf(
      file.path("cellchat_result", str_glue("chord plot - tumor - all slide - {name} - 2.pdf")),
      height = 10,
      width  = 10
    )
    tryCatch(netVisual_chord_gene(cellchat, sources.use = "Epithelial cells", targets.use = NULL), error = function(x) x)
    tryCatch(netVisual_chord_gene(cellchat, sources.use = NULL, targets.use = "Epithelial cells"), error = function(x) x)
    dev.off()
    
  })

###########
groupSize <- as.numeric(table(cellchat@idents))
groupSize <- 10
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = paste0("Interactions with ", "CN7"),
  # 高亮目标细胞
  color.use = ifelse(levels(cellchat@idents) == "CN7", "red", "grey80"),
  vertex.label.color = ifelse(levels(cellchat@idents) == "CN7", "red", "black")
)

all_cellchat_new_ls %>%
  iwalk(function(cellchat, name){
    
    # 提取细胞类型
    cell_types <- levels(cellchat@idents)
    
    # 计算 groupSize（修复版）
    if (is.list(cellchat@idents)) {
      idents_vec <- cellchat@idents$joint
    } else {
      idents_vec <- cellchat@idents
    }
    groupSize <- as.numeric(table(idents_vec))
    
    # 创建命名颜色向量（关键修复！）
    cols <- setNames(rep("grey80", length(cell_types)), cell_types)
    cols["CN8"] <- "#E63946"  # 红色高亮
    
    # 标签颜色同理
    label_cols <- setNames(rep("black", length(cell_types)), cell_types)
    label_cols["CN8"] <- "#E63946"
    
    pdf(
      file.path(out_dir_cellchat, str_glue("network plot CN8 - {name}.pdf")),
      height = 8,
      width  = 16
    )
    
    par(mfrow = c(1, 2), xpd = TRUE)
    
    # Count 网络
    netVisual_circle(
      cellchat@net$count,
      vertex.weight = groupSize,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = paste0("Interactions with CN8 (", name, ")"),
      color.use = cols,
      vertex.label.color = label_cols
    )
    
    # Weight 网络
    netVisual_circle(
      cellchat@net$weight,
      vertex.weight = groupSize,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = paste0("Interactions with CN8 (", name, ")"),
      color.use = cols,
      vertex.label.color = label_cols
    )
    
    dev.off()
    
    message(name, " CN8 网络完成")
  })


all_cellchat_new_ls %>%
  iwalk(function(cellchat, name){
    
    # 提取细胞类型
    cell_types <- levels(cellchat@idents)
    
    # 计算 groupSize（修复版）
    if (is.list(cellchat@idents)) {
      idents_vec <- cellchat@idents$joint
    } else {
      idents_vec <- cellchat@idents
    }
    groupSize <- as.numeric(table(idents_vec))
    
    # 创建命名颜色向量（关键修复！）
    cols <- setNames(rep("grey80", length(cell_types)), cell_types)
    cols["Hofbauer Cell"] <- "#E63946"  # 红色高亮
    
    # 标签颜色同理
    label_cols <- setNames(rep("black", length(cell_types)), cell_types)
    label_cols["Hofbauer Cell"] <- "#E63946"
    
    pdf(
      file.path(out_dir_cellchat, str_glue("network plot Hofbauer Cell - {name}.pdf")),
      height = 8,
      width  = 16
    )
    
    par(mfrow = c(1, 2), xpd = TRUE)
    
    # Count 网络
    netVisual_circle(
      cellchat@net$count,
      vertex.weight = groupSize,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = paste0("Interactions with Hofbauer Cell (", name, ")"),
      color.use = cols,
      vertex.label.color = label_cols
    )
    
    # Weight 网络
    netVisual_circle(
      cellchat@net$weight,
      vertex.weight = groupSize,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = paste0("Interactions with Hofbauer Cell (", name, ")"),
      color.use = cols,
      vertex.label.color = label_cols
    )
    
    dev.off()
    
    message(name, " EVT 网络完成")
  })
####################
pdf(
  file.path(out_dir_cellchat, str_glue("STAB1 Macro bubble.pdf")),
  height = 8,
  width  = 8
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = "STAB1+ Macro",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "STAB1+ Macro",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = "STAB1+ Macro",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "STAB1+ Macro",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = "STAB1+ Macro",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "STAB1+ Macro",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = "STAB1+ Macro",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "STAB1+ Macro",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
dev.off()
########
####################
pdf(
  file.path(out_dir_cellchat, str_glue("MMP9 Macro bubble.pdf")),
  height = 8,
  width  = 8
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = "MMP9+ Macro",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "MMP9+ Macro",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = "MMP9+ Macro",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "MMP9+ Macro",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = "MMP9+ Macro",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "MMP9+ Macro",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = "MMP9+ Macro",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "MMP9+ Macro",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
dev.off()

#########
pdf(
  file.path(out_dir_cellchat, str_glue("CD8+TEMRA bubble.pdf")),
  height = 8,
  width  = 8
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = "CD8+TEMRA",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD8+TEMRA",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = "CD8+TEMRA",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD8+TEMRA",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = "CD8+TEMRA",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD8+TEMRA",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = "CD8+TEMRA",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD8+TEMRA",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
dev.off()
#########
pdf(
  file.path(out_dir_cellchat, str_glue("CD8+TEM bubble.pdf")),
  height = 8,
  width  = 8
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = "CD8+TEM",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD8+TEM",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = "CD8+TEM",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD8+TEM",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = "CD8+TEM",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD8+TEM",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = "CD8+TEM",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD8+TEM",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
dev.off()
#########
pdf(
  file.path(out_dir_cellchat, str_glue("CD4+naive bubble.pdf")),
  height = 8,
  width  = 8
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = "CD4+naive",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD4+naive",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = "CD4+naive",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD4+naive",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = "CD4+naive",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD4+naive",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = "CD4+naive",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "CD4+naive",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
dev.off()
#########
pdf(
  file.path(out_dir_cellchat, str_glue("Th1 bubble.pdf")),
  height = 8,
  width  = 8
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = "Th1",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "Th1",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = "Th1",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "Th1",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = "Th1",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "Th1",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = "Th1",  # EVT 作为发送方
  targets.use = NULL,                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = NULL,  # EVT 作为发送方
  targets.use = "Th1",                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
dev.off()
#######
table(all_data_immnue_meta$subtype)
pdf(
  file.path(out_dir_cellchat, str_glue("immnue bubble 1.pdf")),
  height = 8,
  width  = 16
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = c("STAB1+ Macro","MMP9+ Macro"),  # EVT 作为发送方
  targets.use = c("CD8+TEM","CD8+TEMRA","Gamma-delta_T","CD4+naive","Treg","Th1"),                         # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[1]],
  sources.use = c("CD8+TEM","CD8+TEMRA","Gamma-delta_T","CD4+naive","Treg","Th1"),    # EVT 作为发送方
  targets.use =  c("STAB1+ Macro","MMP9+ Macro"),                       # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[1],
  font.size.title = 12
)


netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = c("STAB1+ Macro","MMP9+ Macro"),  # EVT 作为发送方
  targets.use = c("CD8+TEM","CD8+TEMRA","Gamma-delta_T","CD4+naive","Treg","Th1"),    
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[2]],
  sources.use = c("CD8+TEM","CD8+TEMRA","Gamma-delta_T","CD4+naive","Treg","Th1"),    # EVT 作为发送方
  targets.use =  c("STAB1+ Macro","MMP9+ Macro"),                       # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[2],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = c("STAB1+ Macro","MMP9+ Macro"),  # EVT 作为发送方
  targets.use = c("CD8+TEM","CD8+TEMRA","Gamma-delta_T","CD4+naive","Treg","Th1"), 
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[3]],
  sources.use = c("CD8+TEM","CD8+TEMRA","Gamma-delta_T","CD4+naive","Treg","Th1"),    # EVT 作为发送方
  targets.use =  c("STAB1+ Macro","MMP9+ Macro"),                       # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[3],
  font.size.title = 12
)

netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = c("STAB1+ Macro","MMP9+ Macro"),  # EVT 作为发送方
  targets.use = c("CD8+TEM","CD8+TEMRA","Gamma-delta_T","CD4+naive","Treg","Th1"), 
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
netVisual_bubble(
  all_cellchat_new_ls[[4]],
  sources.use = c("CD8+TEM","CD8+TEMRA","Gamma-delta_T","CD4+naive","Treg","Th1"),    # EVT 作为发送方
  targets.use =  c("STAB1+ Macro","MMP9+ Macro"),                       # 所有接收方
  signaling = NULL,                           # 所有信号通路
  comparison = NULL,                          # 单组分析
  angle.x = 45,                               # X轴标签角度
  font.size = 8,
  title.name = names(all_cellchat_new_ls)[4],
  font.size.title = 12
)
dev.off()
########
netVisual_bubble(
  all_cellchat_new_ls$`T`,
  sources.use = c("Extravillous Trophoblast"),
  targets.use =  c("Vascular Endothelial Cell(VEC)"),
  #signaling = c("CXCL"),  # 包含 MMP12 的通路
  comparison = NULL,
  angle.x = 45
)
netVisual_bubble(
  all_cellchat_new_ls$`C`,
  sources.use = c("Extravillous Trophoblast"),
  targets.use =  c("Vascular Endothelial Cell(VEC)"),
  #signaling = c("CXCL"),  # 包含 MMP12 的通路
  comparison = NULL,
  angle.x = 45
)
unique(all_data$annotation)
file.remove("/mnt/NexData_file/home/lzl/20260415_XE2602040886-02/cellchat_dir/EVT as target.pdf")

#####
genes_of_interest <- c("CXCL16", "CXCR6", "ANGPTL4", "CDH5")
expr_data <- FetchData(all_data, vars = genes_of_interest, layer = "data")

# 合并 metadata
plot_df <- all_data@meta.data %>%
  dplyr::select(annotation, group) %>%
  mutate(cell_id = rownames(all_data@meta.data)) %>%
  left_join(
    expr_data %>% 
      as.data.frame() %>%
      mutate(cell_id = rownames(expr_data)),
    by = "cell_id"
  )


# ============================================
# Step 2: 计算各细胞类型的表达统计
# ============================================

summary_stats <- plot_df %>%
  pivot_longer(cols = all_of(genes_of_interest), 
               names_to = "gene", 
               values_to = "expression") %>%
  group_by(annotation, gene) %>%
  summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    median_expr = median(expression, na.rm = TRUE),
    pct_expressed = sum(expression > 0) / n() * 100,  # 表达比例
    .groups = "drop"
  )

print(summary_stats)
genes_present <- colnames(plot_df)[colnames(plot_df) %in% genes_of_interest]
plot_df_long <- plot_df %>%
  dplyr::select(-cell_id) %>%
  pivot_longer(
    cols = all_of(genes_present),
    names_to = "gene",
    values_to = "expression"
  )

p1 <- ggplot(plot_df_long, aes(x = annotation, y = expression, fill = annotation)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(title = "Gene Expression by Cell Type")
p1

p2 <- ggplot(plot_df_long, aes(x = annotation, y = expression, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.3) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression: Control vs Treatment")
p2

all_data_C <- subset(all_data, subset = group == "C")
all_data_T <- subset(all_data, subset = group == "T")

# C组 DotPlot
out_dir <- "cellchat_dir"
pdf(
  file.path(out_dir, "Dotplot special L-R gene.pdf"),
  height = 6,
  width = 10
)
DotPlot(all_data_C, 
        features = genes_of_interest, 
        group.by = "annotation") +
  coord_flip() +
  ggtitle("C") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# T组 DotPlot
DotPlot(all_data_T, 
        features = genes_of_interest, 
        group.by = "annotation") +
  coord_flip() +
  ggtitle("T") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
dev.off()
p_C
p_T
##
target_genes <- c("CXCL16", "CXCR6")
rownames(expr_matrix)
annotation_col <- 'annotation'
HRAS_position_ls <- map(cells_position_ls, function(sample_df) {
  
  # 3.1 从 Seurat 对象中提取当前样本细胞的表达量
  # 注意：这里假设 sample_df$cell_id 中的 ID 与 all_data 的 colnames 完全匹配
  # FetchData 会自动处理 cell_id 的匹配
  
  # 提取表达矩阵子集
  gene_expr_matrix <- FetchData(all_data, 
                                vars = c(annotation_col,target_genes), 
                                cells = sample_df$cell_id)
  
  # 3.2 将提取的矩阵转换为 data.frame 并加上 cell_id 列以便合并
  # rownames(gene_expr_matrix) 就是 cell_id
  expr_df <- gene_expr_matrix %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "cell_id")
  
  # 3.3 将表达量数据 (expr_df) 左连接合并到原始坐标数据 (sample_df) 中
  # 使用 left_join 保证原始坐标数据的顺序和行数不变
  final_df <- left_join(sample_df, expr_df, by = "cell_id")
  
  return(final_df)
})

p_hras_contour_ls <- 
  HRAS_position_ls %>% 
  imap(function(coords, sam){
    
    # --- 1. 准备当前样本的绘图数据 ---
    # coords 就是当前样本的 HRAS_position_ls[[i]] (即之前的 dat)
    all_cells <- coords %>% remove_missing()
    hras_positive_cells <- all_cells %>% filter(CXCR6 > 0)
    
    # 如果该样本没有阳性细胞，直接返回空，防止后续报错
    if (nrow(hras_positive_cells) == 0) return(NULL)
    
    # --- 2. 计算当前样本的二维核密度估计 (KDE) ---
    kde_result <- MASS::kde2d(
      x = hras_positive_cells$SPATIAL_1, 
      y = -hras_positive_cells$SPATIAL_2, 
      n = 100            
    )
    
    # 转换为长格式数据框
    kde_df <- data.frame(
      x = rep(kde_result$x, times = length(kde_result$y)),
      y = rep(kde_result$y, each = length(kde_result$x)),
      z = as.vector(kde_result$z)
    )
    
    # --- 3. 绘制当前样本的叠加图 ---
    p <- ggplot() +
      # 【第一层：背景】细胞类型注释底图
      geom_tile(
        data = all_cells,
        aes(x = SPATIAL_1, y = -SPATIAL_2, fill = annotation), # Y轴同步翻转
        height = 8,
        width = 8
      ) +
      scale_fill_manual(values = cell_colors) +
      
      # 【第二层：前景】HRAS 阳性细胞的密度等高线
      geom_contour(data = kde_df, 
                   aes(x = x, y = y, z = z, color = after_stat(level)), 
                   size = 0.8) +
      
      # 设置配色、主题与标签
      scale_color_viridis_c(option = "C") +
      theme_void() +
      coord_equal() +
      labs(fill = "Celltype", color = "Density") +
      theme(
        panel.background = element_rect(fill = "black", color = NA), # 保持黑色背景
        legend.key.size = unit(0.5, 'cm'),       
        legend.text = element_text(size = 10),   
        legend.title = element_text(size = 12)   
      ) +
      ggtitle(sam)+ # 加上样本名称作为标题
      guides(
        fill = guide_legend(order = 1),       # 离散图例用 guide_legend
        color = guide_colorbar(order = 2)     # 连续图例用 guide_colorbar
      )
    
    return(p)
  })
pdf(
  file.path(out_dir, "contour_plot_CXCR6 in each sample.pdf"),
  height = 6,
  width = 10
)
p_hras_contour_ls %>% walk(print)
dev.off()

###
some_genes <- c("CXCL16","CXCR6")
plot_dat <-
  FetchData(
    all_data, 
    vars = c(annotation_col,some_genes)
  ) %>%
  rownames_to_column('cell_id')
gene_colors <- c(
  "CXCL16" = "#38B000",    # 墨绿
  "CXCR6" = "#C77DFF" 
)
p_sp_gene_ls <-
  cells_position_ls %>%
  imap(function(coords, sam){
    some_genes %>%
      setNames(., .) %>%
      map(function(x){
        dat2 <- coords %>% left_join(plot_dat, by = "cell_id")
        stopifnot(x %in% colnames(dat2))
        
        # quantile cutoff
        q_cutoff <- 0.9
        n <- 0
        while( near(quantile(dat2[[x]], q_cutoff), 0) ){
          q_cutoff <- q_cutoff + 9/10^(n + 2)
          print(q_cutoff)
          n <- n + 1
          if(n > 3) break
        }
        
        if(n <= 3){
          dat2[[x]] <- pmin(dat2[[x]], quantile(dat2[[x]], q_cutoff))
        }
        
        val_rng <- dat2[[x]] %>% range()
        height <- 8
        width <- 8
        
        # 获取当前基因的颜色
        gene_color <- gene_colors[x]
        
        p <- dat2 %>%
          ggplot(aes(x = SPATIAL_1, y = SPATIAL_2)) +
          geom_tile(fill = "gray60", height = height, width = width, alpha = 0.3) + 
          geom_tile(
            aes(fill = .data[[x]], alpha = .data[[x]]),
            height = height, 
            width = width
          ) +
          theme_void() +
          coord_equal() +
          scale_y_reverse() +
          # 使用基因指定颜色：白→基因色
          scale_fill_gradient(low = "white", high = gene_color, limits = val_rng) +
          scale_alpha_continuous(range = c(0.25, 1), guide = "none") +
          theme(panel.background = element_rect(fill = "black", color = NA)) +
          ggtitle(str_glue("{sam} - {x}"))
        
        p
      })
  })

# 导出
p_sp_gene_ls %>%
  iwalk(function(p_ls, sam){
    pdf(
      file.path(out_dir, str_glue("spatial plot of cellchat special gene in sample {sam}.pdf")),
      height = 8,
      width = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })
