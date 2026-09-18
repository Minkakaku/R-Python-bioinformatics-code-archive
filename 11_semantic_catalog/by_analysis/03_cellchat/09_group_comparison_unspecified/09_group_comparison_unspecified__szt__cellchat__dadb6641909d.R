library(tidyverse)
library(Seurat)
library(CellChat)
library(NMF)      # 关键依赖，提供矩阵计算功能
library(ggalluvial) # 用于流式图绘制，影响内部数据处理


# define default-------------
out_dir_cellchat <- "cellchat_dir" 
dir.create(out_dir_cellchat)

# read data  -------
all_data <- readRDS("all_data.rds") # YOUR XENIUM DATA

all_data@meta.data$group.y %>% unique()
# "N0_CA","N1_CA","Ni_CA","Nmi_CA"
# "N0_LN","N1_LN","Ni_LN","Nmi_LN"

all_data_meta_all     <- all_data[[]]
all_data_norm_all     <- LayerData(JoinLayers(all_data), 'data')
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all)))

# only used 50,000 cells
#set.seed(12345)
subset_cell_idx <- 
  all_data_meta_all %>%
  filter(group.y %in% c("N1_CA","N1_LN")) %>% 
  rownames_to_column("cell_id") #%>%
  #slice_sample(prop = signif(500000/nrow(all_data_meta_all)))

all_data_meta <- all_data_meta_all[subset_cell_idx$cell_id, ]
all_data_norm <- all_data_norm_all[, subset_cell_idx$cell_id]

all_grp_ls <-
  all_data_meta %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, group.y) %>%
  group_by(group.y) %>%
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
    cells <- all_data_meta_all$annotation %>% unique()
    if(!all(cells %in% meta$annotation)){
      missing_cells <- cells %>% setdiff(unique(meta$annotation))
      missing_meta <-
        all_data_meta_all %>%
        rownames_to_column('cell_id') %>%
        dplyr::filter(annotation %in% missing_cells) %>%
        group_by(annotation2) %>%
        slice_head(n = 1) %>%
        column_to_rownames("cell_id")
      meta <- rbind(meta, missing_meta)
      print(missing_meta)
      dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta)])
    }
    
    meta$samples <- factor(meta$group.y)
    
    cellchat <- createCellChat(object = dat, meta = meta, group.by = "annotation")
    
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


saveRDS(all_cellchat_ls, file.path(out_dir_cellchat, "all_cellchat_ls2.rds"))

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


# all_cellchat_new_ls <-
#   all_cellchat_ls %>%
#   imap(function(cellchat, name){
#     cellchat <- tryCatch(netAnalysis_computeCentrality(cellchat), error = function(e) NULL)
#     if(is.null(cellchat)){
#       cat(str_glue("group {name} cannot run netAnalysis_computeCentrality and removed!\n\n"))
#     }
#     cellchat
#   }) %>% 
#   discard(is.null)  

all_cellchat_new_ls <- all_cellchat_ls %>%
  imap(function(cellchat, name) {
    range <- NULL  # 显式定义局部变量
    tryCatch(
      { cellchat <- netAnalysis_computeCentrality(cellchat) },
      error = function(e) {
        message("Error in group ", name, ": ", e$message)
        NULL
      }
    )
  })


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
      # comparison  = x,
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
      # comparison  = x,
      stacked       = F,
      do.stat       = TRUE
    ) +
      ggtitle("Information flow - weight")
    print(p3)
    
    dev.off()
  })


## diff incoming and outgoing cells --------
# outgoing
pdf(
  file.path(out_dir_cellchat, str_glue("contribution of signals to cell groups - outgoing.pdf")),
  height = 6, 
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
  height = 6, 
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
  height = 6, 
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





#1、所有细胞大类循环各时期分析----
library(tidyverse)
library(Seurat)
library(CellChat)
library(NMF)      # 关键依赖，提供矩阵计算功能
library(ggalluvial) # 用于流式图绘制，影响内部数据处理


# define default-------------
# out_dir_cellchat <- "cellchat_dir" 
# dir.create(out_dir_cellchat)

# read data  -------
all_data <- readRDS("all_data.rds") # YOUR XENIUM DATA

all_data@meta.data$group.y %>% unique()
# "N0_CA","N1_CA","Ni_CA","Nmi_CA"
# "N0_LN","N1_LN","Ni_LN","Nmi_LN"

group_pairs <- list(
  c("N0_CA", "N0_LN"),
  c("N1_CA", "N1_LN"),
  c("Ni_CA", "Ni_LN"),
  c("Nmi_CA", "Nmi_LN")
)

all_data_meta_all     <- all_data[[]]
all_data_norm_all     <- LayerData(JoinLayers(all_data), 'data')
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all)))


# 循环处理每个组对
for (pair in group_pairs){
  group1 <- pair[1]
  group2 <- pair[2]
  
  # 创建动态输出目录（格式：N0_CA_vs_N0_LN）
  pair_name <- paste0("Cellchat_", group1, "_vs_", group2)
  out_dir_pair <- file.path(out_dir_cellchat, pair_name)
  dir.create(out_dir_pair, showWarnings = FALSE, recursive = TRUE)
  
  cat("正在分析组对:", group1, "vs", group2, "\n")
  
  # only used 50,000 cells
  #set.seed(12345)
  subset_cell_idx <- 
    all_data_meta_all %>%
    filter(group.y %in% c(group1, group2)) %>% 
    rownames_to_column("cell_id") #%>%
  #slice_sample(prop = signif(500000/nrow(all_data_meta_all)))
  
  all_data_meta <- all_data_meta_all[subset_cell_idx$cell_id, ]
  all_data_norm <- all_data_norm_all[, subset_cell_idx$cell_id]
  
  all_grp_ls <-
    all_data_meta %>%
    rownames_to_column("cell_id") %>%
    dplyr::select(cell_id, group.y) %>%
    group_by(group.y) %>%
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
      cells <- all_data_meta_all$annotation %>% unique()
      if(!all(cells %in% meta$annotation)){
        missing_cells <- cells %>% setdiff(unique(meta$annotation))
        missing_meta <-
          all_data_meta_all %>%
          rownames_to_column('cell_id') %>%
          dplyr::filter(annotation %in% missing_cells) %>%
          group_by(annotation2) %>%
          slice_head(n = 1) %>%
          column_to_rownames("cell_id")
        meta <- rbind(meta, missing_meta)
        print(missing_meta)
        dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta)])
      }
      
      meta$samples <- factor(meta$group.y)
      
      cellchat <- createCellChat(object = dat, meta = meta, group.by = "annotation")
      
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
  
  
  saveRDS(all_cellchat_ls, file.path(out_dir_pair, "all_cellchat_ls2.rds"))
  
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
          file      = file.path(out_dir_pair, str_glue("net_lr - {name}.txt")),
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
          file      = file.path(out_dir_pair, str_glue("netPathways_lr - {name}.txt")),
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
  
  
  all_cellchat_new_ls <- all_cellchat_ls %>%
    imap(function(cellchat, name) {
      range <- NULL  # 显式定义局部变量
      tryCatch(
        { cellchat <- netAnalysis_computeCentrality(cellchat) },
        error = function(e) {
          message("Error in group ", name, ": ", e$message)
          NULL
        }
      )
    }) 
  
  groups   <- names(all_cellchat_new_ls)
  cellchat <- mergeCellChat(all_cellchat_new_ls, add.names = groups, cell.prefix = TRUE)
  
  
  saveRDS(all_cellchat_new_ls, file.path(out_dir_pair, "all_cellchat_new_ls.rds")) 
  saveRDS(cellchat, file.path(out_dir_pair, "cellchat.rds")) 
  
  
  all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)
  
  
  ## comprare interation ------
  pdf(
    file.path(out_dir_pair, str_glue("comprare interation plot.pdf")),
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
        file.path(out_dir_pair, str_glue("diff network plot - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
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
        file.path(out_dir_pair, str_glue("diff heatmap - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
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
        file.path(out_dir_pair, str_glue("diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
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
        # comparison  = x,
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
        # comparison  = x,
        stacked       = F,
        do.stat       = TRUE
      ) +
        ggtitle("Information flow - weight")
      print(p3)
      
      dev.off()
    })
  
  
  ## diff incoming and outgoing cells --------
  # outgoing
  pdf(
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - outgoing.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - incoming.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - all.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("diff interaction strength - all groups.pdf")),
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
}



L2 <- c("B cells","pDCs","CD4_Tfh","CD8_Tem","Naive_T","T_like_Macrophage","Tfh","B_like_Macrophage","cycling macrophage", "Treg","CD8_naive","proliferating T","M2","Plasma", "M1","MergeDC","CD4_naive", "DC_like_Macrophage", "Epithelial_like_Macrophage","CD8_MAIT","Fibroblast_like_Macrophage","Mast_like_Macrophage","Th2")
#2.1、上皮和免疫大类循环各时期分析----
library(tidyverse)
library(Seurat)
library(CellChat)
library(NMF)      # 关键依赖，提供矩阵计算功能
library(ggalluvial) # 用于流式图绘制，影响内部数据处理


# define default-------------
# out_dir_cellchat <- "cellchat_dir" 
# dir.create(out_dir_cellchat)

# read data  -------
all_data <- readRDS("all_data.rds") # YOUR XENIUM DATA

all_data@meta.data$group.y %>% unique()
# "N0_CA","N1_CA","Ni_CA","Nmi_CA"
# "N0_LN","N1_LN","Ni_LN","Nmi_LN"

group_pairs <- list(
  c("N0_CA", "N0_LN"),
  c("N1_CA", "N1_LN"),
  c("Ni_CA", "Ni_LN"),
  c("Nmi_CA", "Nmi_LN")
)

all_data_meta_all     <- all_data[[]]
all_data_norm_all     <- LayerData(JoinLayers(all_data), 'data')
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all)))

all_data_meta_all$cell_id2 <-  all_data_meta_all$cell_id
all_data_meta_all$cell_id <- NULL

# 循环处理每个组对
for (pair in group_pairs){
  group1 <- pair[1]
  group2 <- pair[2]
  
  # 创建动态输出目录（格式：N0_CA_vs_N0_LN）
  pair_name <- paste0("Cellchat_", group1, "_vs_", group2)
  out_dir_pair <- file.path(out_dir_cellchat, pair_name)
  dir.create(out_dir_pair, showWarnings = FALSE, recursive = TRUE)
  
  cat("正在分析组对:", group1, "vs", group2, "\n")
  
  # only used 50,000 cells
  #set.seed(12345)
  subset_cell_idx <- 
    all_data_meta_all %>%
    filter(annotation %in% c("B cells","pDCs","T cells","Macrophage","Plasma","MergeDC","Epithelial")) %>% 
    filter(group.y %in% c(group1, group2)) %>% 
    rownames_to_column("cell_id") #%>%
  #slice_sample(prop = signif(500000/nrow(all_data_meta_all)))
  
  all_data_meta <- all_data_meta_all[subset_cell_idx$cell_id, ]
  
  
  all_data_meta$annotation <- all_data_meta$annotation %>% as.character()
  #levels(all_data_meta$annotation) %>% factor(all_data_meta$annotation %>% unique())
  all_data_meta$annotation <- all_data_meta$annotation %>% as.factor()
  all_data_meta$annotation2 <- all_data_meta$annotation2 %>% as.character()
  all_data_meta$annotation2 <- all_data_meta$annotation2 %>% as.factor()
  
  
  all_data_norm <- all_data_norm_all[, subset_cell_idx$cell_id]
  
  all_grp_ls <-
    all_data_meta %>%
    rownames_to_column("cell_id") %>%
    dplyr::select(cell_id, group.y) %>%
    group_by(group.y) %>%
    nest() %>%
    deframe() %>%
    map(unlist, use.name = FALSE)
  
  
  # using CellChat db --------
  all_cellchat_ls <-
    all_grp_ls %>%
    map(function(cells){ #browser()
      
      meta <- all_data_meta[cells, ] 
      dat  <- all_data_norm[, cells]
      
      # table
      cells <- all_data_meta$annotation %>% unique()
      if(!all(cells %in% meta$annotation)){
        missing_cells <- cells %>% setdiff(unique(meta$annotation))
        missing_meta <-
          all_data_meta %>%
          rownames_to_column('cell_id') %>%
          dplyr::filter(annotation %in% missing_cells) %>%
          group_by(annotation) %>%
          slice_head(n = 1) %>%
          column_to_rownames("cell_id")
        meta <- rbind(meta, missing_meta)
        print(missing_meta)
        dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta)])
      }
      
      meta$samples <- factor(meta$group.y)
      ###强制更新因子水平,
      meta$annotation2 <- droplevels(meta$annotation2)
      
      cellchat <- createCellChat(object = dat, meta = meta, group.by = "annotation")
      
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
  
  
  saveRDS(all_cellchat_ls, file.path(out_dir_pair, "all_cellchat_ls2.rds"))
  
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
          file      = file.path(out_dir_pair, str_glue("net_lr - {name}.txt")),
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
          file      = file.path(out_dir_pair, str_glue("netPathways_lr - {name}.txt")),
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
  
  
  all_cellchat_new_ls <- all_cellchat_ls %>%
    imap(function(cellchat, name) {
      range <- NULL  # 显式定义局部变量
      tryCatch(
        { cellchat <- netAnalysis_computeCentrality(cellchat) },
        error = function(e) {
          message("Error in group ", name, ": ", e$message)
          NULL
        }
      )
    }) 
  
  groups   <- names(all_cellchat_new_ls)
  cellchat <- mergeCellChat(all_cellchat_new_ls, add.names = groups, cell.prefix = TRUE)
  
  
  saveRDS(all_cellchat_new_ls, file.path(out_dir_pair, "all_cellchat_new_ls.rds")) 
  saveRDS(cellchat, file.path(out_dir_pair, "cellchat.rds")) 
  
  
  all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)
  
  
  ## comprare interation ------
  pdf(
    file.path(out_dir_pair, str_glue("comprare interation plot.pdf")),
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
        file.path(out_dir_pair, str_glue("diff network plot - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
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
        file.path(out_dir_pair, str_glue("diff heatmap - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
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
        file.path(out_dir_pair, str_glue("diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
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
        # comparison  = x,
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
        # comparison  = x,
        stacked       = F,
        do.stat       = TRUE
      ) +
        ggtitle("Information flow - weight")
      print(p3)
      
      dev.off()
    })
  
  
  ## diff incoming and outgoing cells --------
  # outgoing
  pdf(
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - outgoing.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - incoming.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - all.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("diff interaction strength - all groups.pdf")),
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
}




#2.2、上皮和免疫亚类循环各时期分析----
library(tidyverse)
library(Seurat)
library(CellChat)
library(NMF)      # 关键依赖，提供矩阵计算功能
library(ggalluvial) # 用于流式图绘制，影响内部数据处理


# define default-------------
# out_dir_cellchat <- "cellchat_dir" 
# dir.create(out_dir_cellchat)

# read data  -------
all_data <- readRDS("all_data.rds") # YOUR XENIUM DATA

all_data@meta.data$group.y %>% unique()
# "N0_CA","N1_CA","Ni_CA","Nmi_CA"
# "N0_LN","N1_LN","Ni_LN","Nmi_LN"

group_pairs <- list(
  #c("N0_CA", "N0_LN"),
  #c("N1_CA", "N1_LN"),
  c("Ni_CA", "Ni_LN"),
  c("Nmi_CA", "Nmi_LN")
)

all_data_meta_all     <- all_data[[]]
all_data_norm_all     <- LayerData(JoinLayers(all_data), 'data')
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all)))

all_data_meta_all$cell_id2 <-  all_data_meta_all$cell_id
all_data_meta_all$cell_id <- NULL

# 循环处理每个组对
for (pair in group_pairs){
  group1 <- pair[1]
  group2 <- pair[2]
  
  # 创建动态输出目录（格式：N0_CA_vs_N0_LN）
  pair_name <- paste0("Cellchat_", group1, "_vs_", group2)
  out_dir_pair <- file.path(out_dir_cellchat, pair_name)
  dir.create(out_dir_pair, showWarnings = FALSE, recursive = TRUE)
  
  cat("正在分析组对:", group1, "vs", group2, "\n")
  
  # only used 50,000 cells
  #set.seed(12345)
  subset_cell_idx <- 
    all_data_meta_all %>%
    filter(annotation %in% c("B cells","pDCs","T cells","Macrophage","Plasma","MergeDC","Epithelial")) %>% 
    filter(group.y %in% c(group1, group2)) %>% 
    rownames_to_column("cell_id") #%>%
  #slice_sample(prop = signif(500000/nrow(all_data_meta_all)))
  
  all_data_meta <- all_data_meta_all[subset_cell_idx$cell_id, ]
  
  
  all_data_meta$annotation <- all_data_meta$annotation %>% as.character()
  #levels(all_data_meta$annotation) %>% factor(all_data_meta$annotation %>% unique())
  all_data_meta$annotation <- all_data_meta$annotation %>% as.factor()
  all_data_meta$annotation2 <- all_data_meta$annotation2 %>% as.character()
  all_data_meta$annotation2 <- all_data_meta$annotation2 %>% as.factor()
  
  
  all_data_norm <- all_data_norm_all[, subset_cell_idx$cell_id]
  
  all_grp_ls <-
    all_data_meta %>%
    rownames_to_column("cell_id") %>%
    dplyr::select(cell_id, group.y) %>%
    group_by(group.y) %>%
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
      cells <- all_data_meta$annotation %>% unique()
      if(!all(cells %in% meta$annotation)){
        missing_cells <- cells %>% setdiff(unique(meta$annotation))
        missing_meta <-
          all_data_meta %>%
          rownames_to_column('cell_id') %>%
          dplyr::filter(annotation %in% missing_cells) %>%
          group_by(annotation) %>%
          slice_head(n = 1) %>%
          column_to_rownames("cell_id")
        meta <- rbind(meta, missing_meta)
        print(missing_meta)
        dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta)])
      }
      
      meta$samples <- factor(meta$group.y)
      ###强制更新因子水平,
      meta$annotation2 <- droplevels(meta$annotation2)
      
      cellchat <- createCellChat(object = dat, meta = meta, group.by = "annotation2")
      
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
  
  
  saveRDS(all_cellchat_ls, file.path(out_dir_pair, "all_cellchat_ls2.rds"))
  
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
          file      = file.path(out_dir_pair, str_glue("net_lr - {name}.txt")),
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
          file      = file.path(out_dir_pair, str_glue("netPathways_lr - {name}.txt")),
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
  
  
  all_cellchat_new_ls <- all_cellchat_ls %>%
    imap(function(cellchat, name) {
      range <- NULL  # 显式定义局部变量
      tryCatch(
        { cellchat <- netAnalysis_computeCentrality(cellchat) },
        error = function(e) {
          message("Error in group ", name, ": ", e$message)
          NULL
        }
      )
    }) 
  
  groups   <- names(all_cellchat_new_ls)
  cellchat <- mergeCellChat(all_cellchat_new_ls, add.names = groups, cell.prefix = TRUE)
  
  
  saveRDS(all_cellchat_new_ls, file.path(out_dir_pair, "all_cellchat_new_ls.rds")) 
  saveRDS(cellchat, file.path(out_dir_pair, "cellchat.rds")) 
  
  
  all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)
  
  
  ## comprare interation ------
  pdf(
    file.path(out_dir_pair, str_glue("comprare interation plot.pdf")),
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
        file.path(out_dir_pair, str_glue("diff network plot - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
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
        file.path(out_dir_pair, str_glue("diff heatmap - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
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
        file.path(out_dir_pair, str_glue("diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
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
        # comparison  = x,
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
        # comparison  = x,
        stacked       = F,
        do.stat       = TRUE
      ) +
        ggtitle("Information flow - weight")
      print(p3)
      
      dev.off()
    })
  
  
  ## diff incoming and outgoing cells --------
  # outgoing
  pdf(
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - outgoing.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - incoming.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - all.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("diff interaction strength - all groups.pdf")),
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
}





#3、免疫免疫大类循环各时期分析----
library(tidyverse)
library(Seurat)
library(CellChat)
library(NMF)      # 关键依赖，提供矩阵计算功能
library(ggalluvial) # 用于流式图绘制，影响内部数据处理


# define default-------------
# out_dir_cellchat <- "cellchat_dir" 
# dir.create(out_dir_cellchat)

# read data  -------
all_data <- readRDS("all_data.rds") # YOUR XENIUM DATA

all_data@meta.data$group.y %>% unique()
# "N0_CA","N1_CA","Ni_CA","Nmi_CA"
# "N0_LN","N1_LN","Ni_LN","Nmi_LN"

group_pairs <- list(
  c("N0_CA", "N0_LN"),
  c("N1_CA", "N1_LN"),
  c("Ni_CA", "Ni_LN"),
  c("Nmi_CA", "Nmi_LN")
)

all_data_meta_all     <- all_data[[]]
all_data_norm_all     <- LayerData(JoinLayers(all_data), 'data')
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all)))

all_data_meta_all$cell_id2 <-  all_data_meta_all$cell_id
all_data_meta_all$cell_id <- NULL

# 循环处理每个组对
for (pair in group_pairs){
  group1 <- pair[1]
  group2 <- pair[2]
  
  # 创建动态输出目录（格式：N0_CA_vs_N0_LN）
  pair_name <- paste0("Cellchat_", group1, "_vs_", group2)
  out_dir_pair <- file.path(out_dir_cellchat, pair_name)
  dir.create(out_dir_pair, showWarnings = FALSE, recursive = TRUE)
  
  cat("正在分析组对:", group1, "vs", group2, "\n")
  
  # only used 50,000 cells
  #set.seed(12345)
  subset_cell_idx <- 
    all_data_meta_all %>%
    filter(annotation %in% c("B cells","pDCs","T cells","Macrophage","Plasma","MergeDC")) %>% 
    filter(group.y %in% c(group1, group2)) %>% 
    rownames_to_column("cell_id") #%>%
  #slice_sample(prop = signif(500000/nrow(all_data_meta_all)))
  
  all_data_meta <- all_data_meta_all[subset_cell_idx$cell_id, ]
  
  
  all_data_meta$annotation <- all_data_meta$annotation %>% as.character()
  #levels(all_data_meta$annotation) %>% factor(all_data_meta$annotation %>% unique())
  all_data_meta$annotation <- all_data_meta$annotation %>% as.factor()
  all_data_meta$annotation2 <- all_data_meta$annotation2 %>% as.character()
  all_data_meta$annotation2 <- all_data_meta$annotation2 %>% as.factor()
  
  
  all_data_norm <- all_data_norm_all[, subset_cell_idx$cell_id]
  
  all_grp_ls <-
    all_data_meta %>%
    rownames_to_column("cell_id") %>%
    dplyr::select(cell_id, group.y) %>%
    group_by(group.y) %>%
    nest() %>%
    deframe() %>%
    map(unlist, use.name = FALSE)
  
  
  # using CellChat db --------
  all_cellchat_ls <-
    all_grp_ls %>%
    map(function(cells){ #browser()
      
      meta <- all_data_meta[cells, ] 
      dat  <- all_data_norm[, cells]
      
      # table
      cells <- all_data_meta$annotation %>% unique()
      if(!all(cells %in% meta$annotation)){
        missing_cells <- cells %>% setdiff(unique(meta$annotation))
        missing_meta <-
          all_data_meta %>%
          rownames_to_column('cell_id') %>%
          dplyr::filter(annotation %in% missing_cells) %>%
          group_by(annotation) %>%
          slice_head(n = 1) %>%
          column_to_rownames("cell_id")
        meta <- rbind(meta, missing_meta)
        print(missing_meta)
        dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta)])
      }
      
      meta$samples <- factor(meta$group.y)
      ###强制更新因子水平,
      meta$annotation2 <- droplevels(meta$annotation2)
      
      cellchat <- createCellChat(object = dat, meta = meta, group.by = "annotation")
      
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
  
  
  saveRDS(all_cellchat_ls, file.path(out_dir_pair, "all_cellchat_ls2.rds"))
  
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
          file      = file.path(out_dir_pair, str_glue("net_lr - {name}.txt")),
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
          file      = file.path(out_dir_pair, str_glue("netPathways_lr - {name}.txt")),
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
  
  
  all_cellchat_new_ls <- all_cellchat_ls %>%
    imap(function(cellchat, name) {
      range <- NULL  # 显式定义局部变量
      tryCatch(
        { cellchat <- netAnalysis_computeCentrality(cellchat) },
        error = function(e) {
          message("Error in group ", name, ": ", e$message)
          NULL
        }
      )
    }) 
  
  groups   <- names(all_cellchat_new_ls)
  cellchat <- mergeCellChat(all_cellchat_new_ls, add.names = groups, cell.prefix = TRUE)
  
  
  saveRDS(all_cellchat_new_ls, file.path(out_dir_pair, "all_cellchat_new_ls.rds")) 
  saveRDS(cellchat, file.path(out_dir_pair, "cellchat.rds")) 
  
  
  all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)
  
  
  ## comprare interation ------
  pdf(
    file.path(out_dir_pair, str_glue("comprare interation plot.pdf")),
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
        file.path(out_dir_pair, str_glue("diff network plot - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
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
        file.path(out_dir_pair, str_glue("diff heatmap - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
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
        file.path(out_dir_pair, str_glue("diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
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
        # comparison  = x,
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
        # comparison  = x,
        stacked       = F,
        do.stat       = TRUE
      ) +
        ggtitle("Information flow - weight")
      print(p3)
      
      dev.off()
    })
  
  
  ## diff incoming and outgoing cells --------
  # outgoing
  pdf(
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - outgoing.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - incoming.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("contribution of signals to cell groups - all.pdf")),
    height = 6, 
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
    file.path(out_dir_pair, str_glue("diff interaction strength - all groups.pdf")),
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
}

