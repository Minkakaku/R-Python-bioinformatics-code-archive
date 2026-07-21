library(tidyverse)
library(Seurat)
library(CellChat)
source('cellchat_fix.R')


# define default-------------
out_dir_cellchat <- "cellchat_dir" 
dir.create(out_dir_cellchat)
out_dir_cellchat <- "cellchat_dir/Treg_20um"
dir.create(out_dir_cellchat)

all_data <- readRDS('all_data_remove_Tumor3_annotation_full_re.rds')
Idents(all_data) <- "annotation_full_re"

res_ls <- readRDS('Treg_proximity_result/20um/res_ls_20um.rds')
data_meta_pos_ls <- readRDS('Treg_proximity_result/data_meta_pos_ls.rds')
neighbor_cellid_ls <- 
  res_ls %>% 
  imap(function(res, sample){
    index2cellid <- as.character(data_meta_pos_ls[[sample]][['cell_id']])[t(res$nn.idx)]
  })
all_neighbor_cellid <- unlist(neighbor_cellid_ls)


Treg_beighbor_20um <- subset(JoinLayers(all_data), cells = all_neighbor_cellid)

all_data <- Treg_beighbor_20um

all_data_meta_all <- all_data[[]]
all_data_norm_all <- LayerData(all_data, 'data')
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all))) 

all_data_meta <- all_data_meta_all
all_data_norm <- all_data_norm_all

# all_data_meta_all %>% group_by(sample) %>% summarise(count = n()) %>% mutate(prop = count/sum(count))
# all_data_meta %>% group_by(sample) %>% summarise(count = n()) %>% mutate(prop = count/sum(count))

all_data$group %>% table()
all_grp_ls <- 
  all_data_meta %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, group) %>%
  group_by(group) %>%
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
    cells <- all_data_meta_all$annotation_full_re %>% unique()
    if(!all(cells %in% meta$annotation_full_re)){
      missing_cells <- cells %>% setdiff(unique(meta$annotation_full_re))
      missing_meta <-
        all_data_meta_all %>%
        rownames_to_column('cell_id') %>%
        dplyr::filter(annotation_full_re %in% missing_cells) %>%
        group_by(annotation_full_re) %>%
        slice_head(n = 1) %>%
        column_to_rownames("cell_id")
      meta <- rbind(meta, missing_meta)
      print(missing_meta)
      dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta)])
    }
    
    meta$samples <- factor(meta$group)
    
    cellchat <- createCellChat(object = dat, meta = meta, group.by = "annotation_full_re")
    
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
  imap(function(cellchat, name){ # browser()
    cellchat <- tryCatch(netAnalysis_computeCentrality(cellchat), error = function(e) NULL)
    if(is.null(cellchat)){
      cat(str_glue("group {name} cannot run netAnalysis_computeCentrality and removed!\n\n"))
    }
    cellchat
  }) %>% 
  purrr::discard(is.null)  

groups   <- names(all_cellchat_new_ls)
cellchat <- mergeCellChat(all_cellchat_new_ls, add.names = groups, cell.prefix = TRUE)


saveRDS(all_cellchat_new_ls, file.path(out_dir_cellchat, "all_cellchat_new_ls.rds")) 
saveRDS(cellchat, file.path(out_dir_cellchat, "cellchat.rds")) 


all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)
all_combns_ls <- c(all_combns_ls[1], all_combns_ls[3])




# plot ---------------
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
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T, 
      comparison    = x,
      title.name    = "Differential number of interactions",
      vertex.label.cex = 0.8, 
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T,
      comparison    = x,
      measure       = "weight", 
      title.name    = "Differential interaction strength",
      vertex.label.cex = 0.8, 
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    dev.off()
  })


# CD4 Treg
all_combns_ls %>%
  walk(function(x){
    grp <- groups[x]
    pdf(
      file.path(out_dir_cellchat, str_glue("CD4 Treg - diff network plot - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
      height = 10,
      width  = 10
    )
    # par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T, 
      comparison    = x,
      title.name    = "Differential number of interactions",
      targets.use = 'CD4 Treg', 
      vertex.label.cex = 0.8, 
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T,
      comparison    = x,
      measure       = "weight", 
      title.name    = "Differential interaction strength",
      targets.use = 'CD4 Treg', 
      vertex.label.cex = 0.8, 
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T, 
      comparison    = x,
      title.name    = "Differential number of interactions",
      sources.use = 'CD4 Treg', 
      vertex.label.cex = 0.8, 
      margin        = c(0.2, 0.2, 0.5, 0.2),
      arrow.size    = 0.6
    )
    netVisual_diffInteraction_fix(
      cellchat,
      weight.scale  = T,
      comparison    = x,
      measure       = "weight", 
      title.name    = "Differential interaction strength",
      sources.use = 'CD4 Treg', 
      vertex.label.cex = 0.8, 
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


# Treg-source
all_combns_ls %>% 
  walk(function(x){
    grp <- groups[x]
    pdf(
      file.path(out_dir_cellchat, str_glue("Treg-source diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
      height = 6,
      width  = 6
    )
    # par(mfrow = c(1,2), xpd=TRUE)
    p1 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "count",
      sources.use = 'CD4 Treg',
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
      sources.use = 'CD4 Treg',
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
      sources.use = 'CD4 Treg',
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
      sources.use = 'CD4 Treg',
      targets.use = NULL, 
      comparison  = x,
      stacked       = F,
      do.stat       = TRUE
    ) +
      ggtitle("Information flow - weight")
    print(p4)
    
    dev.off()
  })

# Treg-target
all_combns_ls %>% 
  walk(function(x){
    grp <- groups[x]
    pdf(
      file.path(out_dir_cellchat, str_glue("Treg-target diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
      height = 6,
      width  = 6
    )
    # par(mfrow = c(1,2), xpd=TRUE)
    p1 <- rankNet(
      cellchat,
      mode        = "comparison",
      measure     = "count",
      sources.use = NULL,
      targets.use = 'CD4 Treg', 
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
      targets.use = 'CD4 Treg', 
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
      targets.use = 'CD4 Treg', 
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
      targets.use = 'CD4 Treg', 
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
  height = 10, 
  width  = 6
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){ # browser()
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
  height = 10, 
  width  = 6
)

all_cellchat_new_ls %>% 
  iwalk(function(cellchat, name){
    pathways <- cellchat@netP$pathways
    
    tryCatch({
      ## YOUR CODE
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
    },
    error = function(x) {
      print(x$message)
      NULL
    }
    )
    
  })
dev.off()

# all 
pdf(
  file.path(out_dir_cellchat, str_glue("contribution of signals to cell groups - all.pdf")),
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
  width  = 15
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
