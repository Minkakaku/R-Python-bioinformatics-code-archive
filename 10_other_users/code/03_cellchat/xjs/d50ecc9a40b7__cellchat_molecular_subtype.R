library(tidyverse)
library(Seurat)
library(CellChat)
source('cellchat_fix.R')

# define default-------------
out_dir_cellchat <- "cellchat_result_molecular_subtype" 
dir.create(out_dir_cellchat)

all_data <- readRDS('/R_OUTPUT/all_data_4slides_sample_filtered_celltype_filtered_reUMAP.rds')
Idents(all_data) <- "annotation_full"

all_data_meta_all <- all_data[[]]
all_data_norm_all <- LayerData(JoinLayers(all_data), 'data')
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all))) 

# all_data_meta <- all_data_meta_all
# all_data_norm <- all_data_norm_all

set.seed(12345)
subset_cell_idx <- 
  all_data_meta_all %>%
  rownames_to_column("cell_id") %>%
  slice_sample(prop = signif(50000/nrow(all_data_meta_all)))

all_data_meta <- all_data_meta_all[subset_cell_idx$cell_id, ]
all_data_norm <- all_data_norm_all[, subset_cell_idx$cell_id]


# all_data_meta_all %>% group_by(sample) %>% summarise(count = n()) %>% mutate(prop = count/sum(count))
# all_data_meta %>% group_by(sample) %>% summarise(count = n()) %>% mutate(prop = count/sum(count))

all_data$molecular_subtype %>% table()
all_grp_ls <- 
  all_data_meta %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id, molecular_subtype) %>%
  group_by(molecular_subtype) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.name = FALSE)

names(all_grp_ls)
# all_grp_ls <- all_grp_ls[c(3, 2, 1)] # <-------------------


# using CellChat db --------
all_cellchat_ls <-
  all_grp_ls %>%
  map(function(cells){ # browser()
    
    meta <- all_data_meta[cells, ] 
    dat  <- all_data_norm[, cells]
    
    # table
    cells <- all_data_meta_all$annotation_full %>% unique()
    if(!all(cells %in% meta$annotation_full)){
      missing_cells <- cells %>% setdiff(unique(meta$annotation_full))
      missing_meta <-
        all_data_meta_all %>%
        rownames_to_column('cell_id') %>%
        dplyr::filter(annotation_full %in% missing_cells) %>%
        group_by(annotation_full) %>%
        slice_head(n = 1) %>%
        column_to_rownames("cell_id")
      meta <- rbind(meta, missing_meta)
      print(missing_meta)
      dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta), drop=FALSE])
    }
    
    meta$samples <- factor(meta$molecular_subtype)
    
    cellchat <- createCellChat(object = dat, meta = meta, group.by = "annotation_full")
    
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
# all_combns_ls <- c(all_combns_ls[1], all_combns_ls[3])




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
      file.path(out_dir_cellchat, str_glue("diff network plot - {grp[2]} vs {grp[1]} - red high, blue low.pdf")),
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


### CD4 Treg ----
# all_combns_ls %>%
#   walk(function(x){
#     grp <- groups[x]
#     pdf(
#       file.path(out_dir_cellchat, str_glue("CD4 Treg - diff network plot - {grp[2]} vs {grp[1]} - red high, blue low.pdf")),
#       height = 10,
#       width  = 10
#     )
#     # par(mfrow = c(1,2), xpd=TRUE)
#     netVisual_diffInteraction_fix(
#       cellchat,
#       weight.scale  = T, 
#       comparison    = x,
#       title.name    = "Differential number of interactions",
#       targets.use = 'CD4 Treg', 
#       vertex.label.cex = 0.8, 
#       margin        = c(0.2, 0.2, 0.5, 0.2),
#       arrow.size    = 0.6
#     )
#     netVisual_diffInteraction_fix(
#       cellchat,
#       weight.scale  = T,
#       comparison    = x,
#       measure       = "weight", 
#       title.name    = "Differential interaction strength",
#       targets.use = 'CD4 Treg', 
#       vertex.label.cex = 0.8, 
#       margin        = c(0.2, 0.2, 0.5, 0.2),
#       arrow.size    = 0.6
#     )
#     netVisual_diffInteraction_fix(
#       cellchat,
#       weight.scale  = T, 
#       comparison    = x,
#       title.name    = "Differential number of interactions",
#       sources.use = 'CD4 Treg', 
#       vertex.label.cex = 0.8, 
#       margin        = c(0.2, 0.2, 0.5, 0.2),
#       arrow.size    = 0.6
#     )
#     netVisual_diffInteraction_fix(
#       cellchat,
#       weight.scale  = T,
#       comparison    = x,
#       measure       = "weight", 
#       title.name    = "Differential interaction strength",
#       sources.use = 'CD4 Treg', 
#       vertex.label.cex = 0.8, 
#       margin        = c(0.2, 0.2, 0.5, 0.2),
#       arrow.size    = 0.6
#     )
#     dev.off()
#   })


## diff heatmap plot --------
all_combns_ls %>% 
  walk(function(x){
    grp <- groups[x]
    
    pdf(
      file.path(out_dir_cellchat, str_glue("diff heatmap - {grp[2]} vs {grp[1]} - red high, blue low.pdf")),
      height = 10,
      width  = 10
    )
    print(netVisual_heatmap(cellchat, comparison = x))
    print(netVisual_heatmap(cellchat, comparison = x, measure = "weight"))
    dev.off()
  })


## diff rank plot --------
all_combns_ls %>% 
  walk(function(x){
    x <- c(x[2], x[1])
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


### Treg-source ---------
# all_combns_ls %>% 
#   walk(function(x){
#     x <- c(x[2], x[1])
#     grp <- groups[x]
#     pdf(
#       file.path(out_dir_cellchat, str_glue("Treg-source diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
#       height = 6,
#       width  = 6
#     )
#     # par(mfrow = c(1,2), xpd=TRUE)
#     p1 <- rankNet(
#       cellchat,
#       mode        = "comparison",
#       measure     = "count",
#       sources.use = 'CD4 Treg',
#       targets.use = NULL, 
#       comparison  = x,
#       stacked     = T,
#       do.stat     = TRUE
#     ) +
#       ggtitle("Relative information flow - count")
#     print(p1)
#     
#     p2 <- rankNet(
#       cellchat,
#       mode        = "comparison",
#       measure     = "weight",
#       sources.use = 'CD4 Treg',
#       targets.use = NULL, 
#       comparison  = x,
#       stacked       = T,
#       do.stat       = TRUE
#     ) +
#       ggtitle("Relative information flow - weight")
#     print(p2)
#     
#     p3 <- rankNet(
#       cellchat,
#       mode        = "comparison",
#       measure     = "count",
#       sources.use = 'CD4 Treg',
#       targets.use = NULL, 
#       comparison  = x,
#       stacked     = F,
#       do.stat     = TRUE
#     ) +
#       ggtitle("Information flow - count")
#     print(p3)
#     
#     p4 <- rankNet(
#       cellchat,
#       mode        = "comparison",
#       measure     = "weight",
#       sources.use = 'CD4 Treg',
#       targets.use = NULL, 
#       comparison  = x,
#       stacked       = F,
#       do.stat       = TRUE
#     ) +
#       ggtitle("Information flow - weight")
#     print(p4)
#     
#     dev.off()
#   })

### Treg-target --------------
# all_combns_ls %>% 
#   walk(function(x){
#     x <- c(x[2], x[1])
#     grp <- groups[x]
#     pdf(
#       file.path(out_dir_cellchat, str_glue("Treg-target diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
#       height = 6,
#       width  = 6
#     )
#     # par(mfrow = c(1,2), xpd=TRUE)
#     p1 <- rankNet(
#       cellchat,
#       mode        = "comparison",
#       measure     = "count",
#       sources.use = NULL,
#       targets.use = 'CD4 Treg', 
#       comparison  = x,
#       stacked     = T,
#       do.stat     = TRUE
#     ) +
#       ggtitle("Relative information flow - count")
#     print(p1)
#     
#     p2 <- rankNet(
#       cellchat,
#       mode        = "comparison",
#       measure     = "weight",
#       sources.use = NULL,
#       targets.use = 'CD4 Treg', 
#       comparison  = x,
#       stacked       = T,
#       do.stat       = TRUE
#     ) +
#       ggtitle("Relative information flow - weight")
#     print(p2)
#     
#     p3 <- rankNet(
#       cellchat,
#       mode        = "comparison",
#       measure     = "count",
#       sources.use = NULL,
#       targets.use = 'CD4 Treg', 
#       comparison  = x,
#       stacked     = F,
#       do.stat     = TRUE
#     ) +
#       ggtitle("Information flow - count")
#     print(p3)
#     
#     p4 <- rankNet(
#       cellchat,
#       mode        = "comparison",
#       measure     = "weight",
#       sources.use = NULL,
#       targets.use = 'CD4 Treg', 
#       comparison  = x,
#       stacked       = F,
#       do.stat       = TRUE
#     ) +
#       ggtitle("Information flow - weight")
#     print(p4)
#     
#     dev.off()
#   })




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




## chordDiagram 4 `Treg 2 Fibroblast` -----------------
# celltype <- all_data$annotation_full %>% unique()
# 
# df_ct_T <-
#   read.table('cellchat_dir2/whole_slide/net_lr - T.txt', header = T, sep = '\t') %>%
#   dplyr::filter(source == 'CD4 Treg') %>%
#   dplyr::filter(target %in% c('Fibroblast')) %>%
#   select(source, target, prob, interaction_name_2)
# df_ct_T$interaction_name_2 <- str_c(df_ct_T$interaction_name_2, ' (Tumor)')
# 
# df_ct_CP <-
#   read.table('cellchat_dir2/whole_slide/net_lr - CP.txt', header = T, sep = '\t') %>%
#   dplyr::filter(source == 'CD4 Treg') %>%
#   dplyr::filter(target %in% c('Fibroblast')) %>%
#   select(source, target, prob, interaction_name_2)
# df_ct_CP$interaction_name_2 <- str_c(df_ct_CP$interaction_name_2, ' (CP)')
# 
# # set.seed(1234)
# df_ct_NA <-
#   data.frame(
#     source = celltype,
#     target = 'Fibroblast',
#     prob = 0,
#     interaction_name_2 = 'NA'
#   )  %>%
#   dplyr::filter(! source %in% c("CD4 Treg", "Fibroblast")) %>%
#   dplyr::filter(
#     str_detect(source, "^[AB]") |
#       str_detect(source, "^[DEF]")
#   ) 
# 
# 
# df_ct_NA$prob <- min(df_ct_T$prob, df_ct_CP$prob)/2
# # df_ct <- rbind(df_ct_T, df_ct_CP, df_ct_NA)
# df_ct <- rbind(df_ct_NA[1:2, ], df_ct_T, df_ct_CP, df_ct_NA[3:nrow(df_ct_NA), ])
# 
# df_ct$prob <- log1p( df_ct$prob/min(df_ct$prob) )
# # df_ct$source <- factor(df_ct$source)
# # df_ct <- arrange(df_ct, source)
# 
# # color alpha
# alpha_val <-  0.5
# 
# # pathway and pathway color
# col_map <- structure(
#   # ggsci::pal_igv()(length(unique(df_ct$interaction_name_2))),
#   allcolour[1:length(unique(df_ct$interaction_name_2))],
#   names = unique(df_ct$interaction_name_2)
# )
# col_map['NA'] <- alpha('white', 0)
# 
# # linke color
# link.colors <-
#   df_ct$interaction_name_2 %>%
#   fct_recode(!!! structure(names(col_map), names = col_map)) %>%
#   alpha(alpha = alpha_val)
# 
# # plot
# pdf(file.path(out_dir_cellchat, 'chordDiagram - Treg 2 Fibroblast - group.pdf'), height = 8, width = 12)
# # par(mar = c(5, 4, 4, 2))
# par(mar = c(2, 2, 2, 10))
# circlize::circos.par()
# 
# circlize::chordDiagram(
#   df_ct,
#   
#   grid.col = ggsci::pal_igv()(7),
#   # linke color
#   col = link.colors,
#   
#   # transparency
#   transparency = 1-alpha_val,
#   
#   # link is a big arrow
#   directional = 1,
#   direction.type = c("diffHeight", "arrows"),
#   link.arr.length = 0.1,
#   # link.target.prop = 0.05,
#   # target.prop.height = circlize::mm_h(5),
#   # diffHeight = circlize::mm_h(5),
#   # link.arr.type = "triangle",
#   link.arr.type = "big.arrow",
#   big.gap = 15,
#   small.gap = 3
# )
# 
# # add legend
# legend(
#   # 'topright',
#   x = 1.2,
#   y = 0.5,
#   legend = names(col_map[2:3]),
#   fill = alpha(col_map[2:3], alpha_val),
#   bty = "o",
#   xpd = TRUE
# )
# 
# circlize::circos.clear()
# dev.off()




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
