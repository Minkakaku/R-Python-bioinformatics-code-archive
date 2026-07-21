library(tidyverse)
library(Seurat)
library(CellChat)


# define default-------------
out_dir_cellchat <- "cellchat_dir" 
dir.create(out_dir_cellchat)

# read data  -------
all_data <- readRDS("data/adata_CN15.rds")
all_data$project <- all_data$sample
all_data$annotation <- all_data$celltype

all_data <- all_data %>% NormalizeData()
all_data <- all_data %>% RenameCells(add.cell.id = 'cell')
all_data$cell_id <- colnames(all_data)


# CN15
cn15_idx <- all_data$cell_id[all_data$CN == 'CN15']
all_data_meta_all     <- all_data[[]] %>% .[cn15_idx, ]
all_data_norm_all     <- LayerData(all_data, 'data') %>% .[, cn15_idx]
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all)))
 

# only used 50,000 cells
set.seed(12345)
subset_cell_idx <- 
  all_data_meta_all %>%
  slice_sample(n = 50000)

all_data_meta <- all_data_meta_all[subset_cell_idx$cell_id, ]
all_data_norm <- all_data_norm_all[, subset_cell_idx$cell_id]

all_grp_ls <-
  all_data_meta %>%
  dplyr::select(cell_id, group) %>%
  group_by(group) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.name = FALSE)

# NS have no cells(only  6 cells)
all_grp_ls <- all_grp_ls['KS']


# using CellChat db --------
all_cellchat_ls <-
  all_grp_ls %>%
  map(function(cells){ # browser()
    
    meta <- all_data_meta[cells, ] 
    dat  <- all_data_norm[, cells]

    meta$annotation <- meta$annotation %>% droplevels()
    meta$samples <- factor(meta$project)
    
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

# only have 1 groups, not merged cellchat obj
# cellchat <- mergeCellChat(all_cellchat_new_ls, add.names = groups, cell.prefix = TRUE)
cellchat <- all_cellchat_new_ls$KS

saveRDS(all_cellchat_new_ls, file.path(out_dir_cellchat, "all_cellchat_new_ls.rds")) 
saveRDS(cellchat, file.path(out_dir_cellchat, "cellchat.rds")) 


all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)


### ONLY ONLE GROUP, NOT RUN ###
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


### ONLY ONLE GROUP, NOT RUN ###
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


### ONLY ONLE GROUP, NOT RUN ###
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

### ONLY ONLE GROUP, NOT RUN ###
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


