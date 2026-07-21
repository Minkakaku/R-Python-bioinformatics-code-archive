library(tidyverse)
library(Seurat)
library(CellChat)


# set future -------------
library(future)
plan(multisession, workers = 1)
options(future.globals.maxSize = 20 * 1024 ^ 3) # 20Gb




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
out_dir_cellchat <- "cellchat_dir_sample" 
dir.create(out_dir_cellchat)

# read data  -------
all_data_sc2 <- readRDS("all_data_sc2.rds")

all_data_meta_all     <- all_data_sc2[[]]  %>% dplyr::select(-cell_id)
all_data_norm_all     <- LayerData(JoinLayers(all_data_sc2), 'data')  
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all)))


# only used 50,000 cells
set.seed(12345)
subset_cell_idx <- 
  all_data_meta_all %>%
  rownames_to_column("cell_id") %>%
  # group_by(CN2, group, annotation2) %>%
  slice_sample(prop = signif(20000/nrow(all_data_meta_all)))

all_data_meta <- all_data_meta_all[subset_cell_idx$cell_id, ]
all_data_norm <- all_data_norm_all[, subset_cell_idx$cell_id]

all_grp_ls <-
  all_data_meta %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, Sample) %>%
  group_by(Sample) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.name = FALSE)




# using CellChat db --------
all_cellchat_ls <-
  all_grp_ls %>%
  map(function(cells){ # browser()
    
    meta <- all_data_meta[cells, ] 
    dat  <- all_data_norm[, cells]
    
    # table
    cells <- all_data_meta_all$annotation2 %>% unique()
    if(!all(cells %in% meta$annotation2)){
      missing_cells <- cells %>% setdiff(unique(meta$annotation2))
      missing_meta <-
        all_data_meta_all %>%
        rownames_to_column('cell_id') %>%
        dplyr::filter(annotation2 %in% missing_cells) %>%
        group_by(annotation2) %>%
        slice_head(n = 1) %>%
        column_to_rownames("cell_id")
      meta <- rbind(meta, missing_meta)
      print(missing_meta)
      dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta), drop = FALSE])
    }
    
    meta$samples <- factor(meta$project)
    
    cellchat <- createCellChat(object = dat, meta = meta, group.by = "annotation2")
    
    # cellchat db
    CellChatDB.use <- 
      # CellChatDB.human %>% 
      CellChatDB.mouse %>%
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


# all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)
# groups
# 1 P7_Hyperoxia 
# 2 P14_Normoxia 
# 3 P3_Hyperoxia 
# 4 P14_Hyperoxia 
# 5 P3_Normoxia 
# 6 P7_Normoxia

# Hyer 3 1 4
# Norm 5 6 2
all_combns_ls <- list(
  c(3, 1), c(3, 4), c(1 ,4),
  c(5, 6), c(5, 2), c(6, 2),
  c(5, 3), c(6, 1), c(2, 4)
)

## comprare interation ------
pdf(
  file.path(out_dir_cellchat, str_glue("comprare interation plot.pdf")),
  height = 4,
  width  = 6
)
compareInteractions(cellchat, show.legend = FALSE, measure = "count")
compareInteractions(cellchat, show.legend = FALSE, measure = "weight")
dev.off()



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



## diff heatmap plot --------
all_combns_ls %>%
  walk(function(x){
    grp <- groups[x]
    
    pdf(
      file.path(out_dir_cellchat, str_glue("diff heatmap - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
      height = 10,
      width  = 10
    )
    print(netVisual_heatmap(cellchat))
    print(netVisual_heatmap(cellchat, measure = "weight"))
    dev.off()
  })


## diff rank plot --------
all_combns_ls %>%
  walk(function(x){
    grp <- groups[x]
    pdf(
      file.path(out_dir_cellchat, str_glue("diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
      height = 14,
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
# outgoing
pdf(
  file.path(out_dir_cellchat, str_glue("contribution of signals to cell groups - outgoing.pdf")),
  height = 20, 
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
  file.path(out_dir_cellchat, str_glue("contribution of signals to cell groups - incoming.pdf")),
  height = 20, 
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
  file.path(out_dir_cellchat, str_glue("contribution of signals to cell groups - all.pdf")),
  height = 20, 
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
  height = 5 * 2,
  width  = 5 * 3
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


## 

netAnalysis_signalingChanges_scatter(
  cellchat,
  comparison = c(3, 1), # c(4, 2),
  idents.use = c("Capillary 1 cells") #, "Capillary 2 cells")
)



# NOT RUN.
# using specific cells or moleculars

# diff signaling in cells ------
cells_to_used <- c(
  "Capillary 1 cells",
  "Capillary 2 cells",
  'ILC2',
  'Dlk1 AT2', 'Areg AT2',
  'Acoxl AT2', 'Phgdh AT2',
  "Peribronchial fibroblasts",
  "Alveolar fibroblasts"
)

all_combns_ls %>%
  walk(function(x){# browser()
    grp <- groups[x]
    
    pdf(
      file.path(out_dir_cellchat, str_glue("diff signaling in selected cells - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
      height = 6,
      width  = 8
    )
    for(cell in cells_to_used){
      tryCatch({
        p <- netAnalysis_signalingChanges_scatter(
          cellchat,
          comparison = x,
          idents.use = cell
        )
        p <- p + patchwork::plot_annotation(title = cell)
        print(p)
      },
      error = function(e)print(e$message)
    )
    }
    dev.off()
  })




# diff signaling pathway ----- 

# pathway_ls <-
#   list(
#     "P14_Normoxia" = c('RELN', 'NRXN', 'VTN', 'NTS'),
#     "P14_Hyperoxia"  = c(
#       'TWEAK', 'BST2', 'HH', 'LIFR', 'CD34', 'NGL', 'TRAIL', 'LAIR1',
#       'CD80', 'IL6', 'PD-L1', 'BAFF', 'HGF', 'RANKL', 'CX3C', 'NKG2D',
#       'PDL2', 'CD137', 'NEGR'
#     )
#   )

pathway_ls <-
  list(
    # "P7_Normoxia"  = c('KIT','LIFR','CD6','VWF','EGF','PLAU', 'SLITRK'),
    # "P7_Hyperoxia" = c('LT', 'PD-L1', 'CTSG', 'IL17', 'VEGI', 'NEGR')
    
    
    P3_Hyperoxia = c("EGF", "HGF"),
    P3_Normoxia = c("EGF", "HGF"),
    P7_Hyperoxia = c("EGF", "HGF"),
    P7_Normoxia = c("EGF", "HGF"),
    P14_Hyperoxia = c("EGF", "HGF"),
    P14_Normoxia = c("EGF", "HGF")
  )

# path_name <- pathway_ls %>% names() %>% str_c(collapse = ' vs ')
# pathway_ls %>%
#    iwalk(function(pathway, name){
#      cellchat <- all_cellchat_new_ls[[name]]
#      pdf(
#        file.path(out_dir_cellchat, str_glue("{path_name} selected signaling pathway network - {name}.pdf")),
#        height = 5,
#        width  = 10
#      )
#      
#      idents <- cellchat@idents %>% levels()
#      print(
#        cellchat %>% 
#          netAnalysis_signalingRole_network(
#            signaling = pathway,
#            width  = length(idents) / 2,
#            height = 2
#          )
#      )
#      
#      dev.off()
#    })

all_cellchat_new_ls %>%
  iwalk(function(cellchat, name){
    pdf(
      # file.path(out_dir_cellchat, str_glue("{path_name} selected signaling pathway network - {name}.pdf")),
      file.path(out_dir_cellchat, str_glue("selected signaling pathway network - {name} - EGF HGF.pdf")),
      height = 5,
      width  = 10
    )

    idents <- cellchat@idents %>% levels()
    tryCatch(
      print(
        cellchat %>%
          netAnalysis_signalingRole_network(
            signaling = pathway_ls[[name]] %>% .[. %in% cellchat@netP$pathways ],
            width  = length(idents) / 2,
            height = 2
          )
      ),
      error = function(x){
        print(str_glue("{name} - error of {x$message}"))
      }
    )
    

    dev.off()
  })

netAnalysis_computeCentrality
undebug(netAnalysis_signalingRole_network)

all_cellchat_new_ls$P3_Hyperoxia %>% 
  netAnalysis_signalingRole_network(
    signaling = c(
      'ADGRA',
      'SEMA7',
      'LAIR1',
      'CD6',
      'PLAU'
    ),
    width  = 15,
    height = 2
  )


all_cellchat_new_ls$P3_Hyperoxia %>% 
  netAnalysis_signalingRole_network(
    signaling = c(
      'ADGRA',
      'SEMA7',
      'LAIR1',
      'CD6',
      'PLAU',
      'FASLG',
      'TULP',
      'CD96',
      'PCDH',
      'BMP',
      'PVR'
    ),
    width  = 15,
    height = 2
  )



# dotplot -------
all_data %>%
  DotPlot(
    # features = c('Areg', 'Hbegf', 'Egfr'),
    features = c('Tnfsf14',	'Ltbr'),
    group.by = "annotation2",
    split.by = "Oxygen"
  )





# diff molecule ---------
pdf(
  file.path(out_dir_cellchat, str_glue("diff molecule A vs B.pdf")),
  height = 4,
  width  = 6
)

# some cells to plot
# all_cells <- cellchat@idents$joint %>% levels()
c("CD4T-Macrophage", "CD4T", "CD8T", "Epithelial") %>%
  # all_cells %>%
  map(function(cell){
    netAnalysis_signalingChanges_scatter(
      cellchat, 
      idents.use = cell
    )
  })

dev.off()


netAnalysis_signalingChanges_scatter(
  cellchat, 
  idents.use = "Mier1 T cells"
)




# pathway plot----------------
CellChatDB.human$interaction %>%
  dplyr::filter(
    ligand %in% c("CD80", "CD86", "NECTIN", 'MHC-II') |
      receptor  %in% c("CD80", "CD86", "NECTIN", 'MHC-II')
  )

pathways.show <- c("CD80", "CD86", "NECTIN", 'MHC-II')
vertex.receiver <- c("CD4T-Macrophage", "CD4T", "CD8T", "Epithelial")

cell2id_a <- 
  all_cellchat_ls$A@net$pro %>% 
  rownames() %>%
  factor(., levels = .) %>%
  {structure(
    as.numeric(.),
    names = as.character(.)
  )}


netVisual_aggregate(
  all_cellchat_ls$A,
  signaling = pathways.show[3],  
  vertex.receiver = cell2id_a[vertex.receiver],
  layout = "hierarchy"
)


netVisual_aggregate(
  all_cellchat_ls$A,
  signaling = pathways.show[3], 
  layout = "chord"
)


netVisual_heatmap(
  all_cellchat_ls$A,
  signaling = pathways.show[3],  
  color.heatmap = "Reds"
)


netVisual_aggregate(
  all_cellchat_ls$P7_Hyperoxia,
  signaling = 'COLLAGEN',  
  vertex.receiver = 1:10, # cell2id_a[vertex.receiver],
  layout = "hierarchy"
)

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

## plot ohters -------------------

## dotplot of pathway -------------
p_dot_light <-
  DotPlot(
    all_data_sc2,
    features = c('Tnfsf14', 'Ltbr'),
    group.by = "annotation2"
  ) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  )

p_dot_light2 <-
  DotPlot(
    all_data_sc2[, all_data_sc2$annotation2 == "ILC2"],
    features = c('Tnfsf14', 'Ltbr'),
    group.by = "Sample" # "Oxygen" # Sample
  ) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  ) +
  ggtitle('ILC2')

DotPlot(
  all_data_sc2,
  features = c('Tnfsf14', 'Ltbr'),
  group.by = "Sample" # "Oxygen" # Sample
) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  ) 

# Vwf	ITGAV_ITGB3
p_dot_vwf <-
  DotPlot(
    all_data_sc2,
    features = c('Vwf', 'Itgav', 'Itgb3'),
    group.by = "annotation2"
  ) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  )

p_dot_vwf2 <-
  DotPlot(
    all_data_sc2[, all_data_sc2$annotation2 == "ILC2"],
    features = c('Vwf', 'Itgav', 'Itgb3'),
    group.by = "Sample" # "Oxygen" # Sample
  ) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  ) +
  ggtitle('ILC2')


DotPlot(
  all_data_sc2,
  features = c('Vwf', 'Itgav', 'Itgb3'),
  group.by = "Sample" # "Oxygen" # Sample
) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  ) 

pdf(
  file.path(out_dir_cellchat, str_glue("dotplot - LIGHT pathway.pdf")),
  height = 8,
  width  = 6
)
p_dot_light %>% print()
p_dot_light2 %>% print()
dev.off()

pdf(
  file.path(out_dir_cellchat, str_glue("dotplot - VWF pathway.pdf")),
  height = 8,
  width  = 6
)
p_dot_vwf %>% print()
p_dot_vwf2 %>% print()
dev.off()

## network plot of some cells -----
c(
  "ILC2"
) %>% 
  map(function(cell){
    pdf(
      file.path(out_dir_cellchat, str_glue("circle plot - {cell}.pdf")),
      height = 5,
      width  = 8
    )
    on.exit(dev.off())
    
    all_cellchat_new_ls %>% 
      imap(function(cellchat, name){ # browser()
        count_dat  <- cellchat@net$count
        weight_dat <- cellchat@net$weight
        
        count_dat2  <- matrix(0, nrow = nrow(count_dat),  ncol = ncol(count_dat), dimnames = dimnames(count_dat))
        weight_dat2 <- matrix(0, nrow = nrow(weight_dat), ncol = ncol(weight_dat), dimnames = dimnames(weight_dat))
        
        count_dat2_source  <- count_dat2_target  <- count_dat2
        weight_dat2_source <- weight_dat2_target <- weight_dat2
        
        count_dat2_source[cell, ] <- count_dat[cell, ]
        count_dat2_target[, cell] <- count_dat[, cell]
        
        weight_dat2_source[cell, ] <- weight_dat[cell, ]
        weight_dat2_target[, cell] <- weight_dat[, cell]
        
        # par
        opar <- par(no.readonly = TRUE)
        par(mar = c(2,2,5,2), xpd = TRUE)
        on.exit(par(opar))
        
        # source count
        print(netVisual_circle(
          count_dat2_source,
          vertex.weight = cellchat@idents %>% table %>% as.numeric(), 
          weight.scale  = T, 
          label.edge    = F, 
          title.name    = str_glue("{name} - {cell} - Number of interactions")
        ))
        
        # source weigth
        print(netVisual_circle(
          count_dat2_target,
          vertex.weight = cellchat@idents %>% table %>% as.numeric(), 
          weight.scale  = T, 
          label.edge    = F, 
          title.name    = str_glue("{name} - {cell} - Number of interactions")
        ))
        
        # target count 
        print(netVisual_circle(
          weight_dat2_source,
          vertex.weight = cellchat@idents %>% table %>% as.numeric(), 
          weight.scale  = T, 
          label.edge    = F, 
          title.name    = str_glue("{name} - {cell} - Interaction weights/strength")
        ))
        
        # target weight
        print(netVisual_circle(
          weight_dat2_target,
          vertex.weight = cellchat@idents %>% table %>% as.numeric(), 
          weight.scale  = T, 
          label.edge    = F, 
          title.name    = str_glue("{name} - {cell} - Interaction weights/strength")
        ))
        
      })
  })


## bubble plot ---------------
all_lr_dat <-
  list.files(out_dir_cellchat, pattern = "net_lr", full.names = TRUE) %>%
  setNames(., str_remove(str_remove(basename(.), "^.*\\s"), "\\..*$") ) %>%
  map(read.table, sep = "\t", header = TRUE) %>%
  enframe(name = "cluster") %>%
  unnest("value")




# p_ct_cellchat_light <-
#   all_lr_dat %>%
#   # mutate(cluster = factor(cluster, levels = c("Ctrl", "CD", "CD-F"))) %>%
#   dplyr::filter(source %in% c("ILC2")) %>%
#   dplyr::filter(pathway_name  %in% c("LIGHT")) %>%
#   
#   
#   ggplot(aes(x = relation, y = interaction_name)) +
#   geom_point(aes(size = pval, color = prob)) +
#   facet_wrap(~ cluster) +
#   scale_color_gradientn(colors = Seurat:::SpatialColors(30), breaks = function(x) x, labels = c("min", "max")) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
#     plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
#   )


p_ct_cellchat_LIGHT_VWF <-
  all_lr_dat %>%
  # mutate(cluster = factor(cluster, levels = c("Ctrl", "CD", "CD-F"))) %>%
  dplyr::filter(target %in% c("ILC2") | source %in% c("ILC2")) %>%
  dplyr::filter(pathway_name  %in% c("LIGHT", "VWF")) %>%
  
  
  ggplot(aes(x = relation, y = interaction_name)) +
  geom_point(aes(size = pval, color = prob)) +
  facet_wrap(~ cluster, scales = "free_x") +
  scale_color_gradientn(colors = Seurat:::SpatialColors(30), breaks = function(x) x, labels = c("min", "max")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )

# pdf(
#   file.path(out_dir_cellchat, str_glue("bubble plot - LIGHT pathway.pdf")),
#   height = 4,
#   width  = 8
# )
# p_ct_cellchat_light %>% print()
# dev.off()

pdf(
  file.path(out_dir_cellchat, str_glue("bubble plot - LIGHT - VWF pathway.pdf")),
  height = 4,
  width  = 8
)
p_ct_cellchat_LIGHT_VWF %>% print()
dev.off()



## de analysis -------------
de_tab_AEC <-
  all_data_sc2 %>%
  {Idents(.) <- "annotation2"; .} %>% 
  JoinLayers() %>% 
  FindMarkers(
    ident.1 = 'Hyperoxia',
    ident.2 = 'Normoxia',
    group.by = "Oxygen",
    subset.ident = 'Arterial endothelial cells'
  )

de_tab_VEC <-
  all_data_sc2 %>%
  {Idents(.) <- "annotation2"; .} %>% 
  JoinLayers() %>% 
  FindMarkers(
    ident.1 = 'Hyperoxia',
    ident.2 = 'Normoxia',
    group.by = "Oxygen",
    subset.ident = 'Venous endothelial cells'
  )


p_volcano_AEC <- 
  de_tab_AEC %>%
  rownames_to_column('gene') %>%
  mutate(`-log10(padj)` = -log10(p_val_adj)) %>%  
  mutate(Group = case_when(
    avg_log2FC >  1 & p_val_adj < 0.05 ~ "Pos",
    avg_log2FC < -1 & p_val_adj < 0.05 ~ "Neg",
    TRUE ~ "non-sig"
  )) %>%
  ggplot(aes(x = avg_log2FC, y = `-log10(padj)`)) +
  geom_point(aes(color = Group), show.legend = F) +
  ggrepel::geom_text_repel(
    data = ~ slice_min(dplyr::filter(., Group == "Pos"), p_val_adj, n = 10),
    aes(label = gene)
  )+
  ggrepel::geom_text_repel(
    data = ~ slice_min(dplyr::filter(., Group == "Neg"), p_val_adj, n = 10),
    aes(label = gene)
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_color_manual(values = c(
    Pos = "red",
    Neg = "blue",
    "non-sig" = "gray"
  )) +
  theme_bw() +
  ggtitle("Hyperoxia vs Normoxia - AEC")


p_volcano_VEC <- 
  de_tab_VEC %>%
  rownames_to_column('gene') %>%
  mutate(`-log10(padj)` = -log10(p_val_adj)) %>%  
  mutate(Group = case_when(
    avg_log2FC >  1 & p_val_adj < 0.05 ~ "Pos",
    avg_log2FC < -1 & p_val_adj < 0.05 ~ "Neg",
    TRUE ~ "non-sig"
  )) %>%
  ggplot(aes(x = avg_log2FC, y = `-log10(padj)`)) +
  geom_point(aes(color = Group), show.legend = F) +
  ggrepel::geom_text_repel(
    data = ~ slice_min(dplyr::filter(., Group == "Pos"), p_val_adj, n = 10),
    aes(label = gene)
  )+
  ggrepel::geom_text_repel(
    data = ~ slice_min(dplyr::filter(., Group == "Neg"), p_val_adj, n = 10),
    aes(label = gene)
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_color_manual(values = c(
    Pos = "red",
    Neg = "blue",
    "non-sig" = "gray"
  )) +
  theme_bw() +
  ggtitle("Hyperoxia vs Normoxia - VEC")




pdf(
  file.path(out_dir_cellchat, str_glue("volcano plot - VEC and AEC.pdf")),
  height = 5,
  width  = 6
)
p_volcano_AEC %>% print()
p_volcano_VEC %>% print()
dev.off()

# dotplot 

DotPlot(
  all_data_sc2[, all_data_sc2$annotation2 == "Arterial endothelial cells"],
  features = c(
    'Ctla2a', 'Phlda3', 'Sfrp2', 'Olfm2', 'Plvap',
    'Xist', 'Mgll', 'Sparcl1', 'Mctp1', 'Grrp1'
  ),
  group.by = "Sample" # "Oxygen" # Sample
) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  ) +
  ggtitle('Arterial endothelial cells')


DotPlot(
  all_data_sc2[, all_data_sc2$annotation2 == "Venous endothelial cells"],
  features = c(
    'Rps27l', 'Ccnd1', 'Bax', 'Cdkn1a', 'Tinagl1',
    'Gpihbp1', 'Rnf144a', 'Aqp1', 'Slfn5', 'Cd36'
  ),
  group.by = "Sample" # "Oxygen" # Sample
) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  ) +
  ggtitle('Venous endothelial cells')



# dotplot of IL33 -----------
DotPlot(
  all_data_sc2[, all_data_sc2$annotation2 == "Arterial endothelial cells"],
  features = c(
    'Ctla2a', 'Phlda3', 'Sfrp2', 'Olfm2', 'Plvap',
    'Xist', 'Mgll', 'Sparcl1', 'Mctp1', 'Grrp1'
  ),
  group.by = "Sample" # "Oxygen" # Sample
) +
  scale_color_gradientn(
    # colours = colorRampPalette(c())(50)
    colours = Seurat:::SpatialColors(50)
  ) +
  ggtitle('Arterial endothelial cells')


