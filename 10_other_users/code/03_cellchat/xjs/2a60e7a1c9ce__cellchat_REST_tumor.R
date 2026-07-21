library(tidyverse)
library(Seurat)
library(CellChat)


# define default-------------
out_dir_cellchat <- "cellchat_dir" 
out_dir_cellchat <- file.path(out_dir_cellchat, 'REST tumor')
dir.create(out_dir_cellchat)

# read data  -------
all_data <- readRDS("all_data_annotation_filtered.rds")
all_data$annotation_full %>% table()
all_data <- subset(all_data, annotation_full %in% c('CD8 Tem', 'M1 macrophage', 'REST_tumor cell'))
# all_data_cn_meta <- readRDS("CN_result/all_data_meta_cn.rds")
# all_data_cn_meta$annotation_full %>% table()

all_data_meta_all <- 
  all_data[[]] # %>% 
  # rownames_to_column('cell_id') %>% 
  # left_join(all_data_cn_meta[c('cell_id', 'CN')], by = 'cell_id') %>% 
  # column_to_rownames('cell_id')
all_data_meta_all %>% head()
all_data_norm_all <- LayerData(JoinLayers(all_data), 'data')
stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all)))

# # CN
# CN_idx <- all_data_cn_meta$cell_id[all_data_cn_meta$CN == c('CN1', 'CN2', 'CN7')]
# all_data_meta_all     <- all_data_meta_all %>% .[CN_idx, ]
# all_data_norm_all     <- all_data_norm_all %>% .[, CN_idx]
# stopifnot(all(rownames(all_data_meta_all) == colnames(all_data_norm_all))) 

# # only used 50,000 cells
# set.seed(12345)
# subset_cell_idx <- 
#   all_data_meta_all %>%
#   rownames_to_column("cell_id") %>%
#   slice_sample(prop = signif(50000/nrow(all_data_meta_all)))
# 
# all_data_meta <- all_data_meta_all[subset_cell_idx$cell_id, ]
# all_data_norm <- all_data_norm_all[, subset_cell_idx$cell_id]

all_data_meta <- all_data_meta_all
all_data_norm <- all_data_norm_all

all_grp_ls <-
  all_data_meta %>%
  rownames_to_column("cell_id") %>%
  dplyr::select(cell_id) # %>%
  # group_by(CN) %>%
  # nest() %>%
  # deframe() %>%
  map(unlist, use.name = FALSE)

# all_data_meta$CN %>% table()
names(all_grp_ls) <- 'whole_slide'


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
      dat <- cbind(dat, all_data_norm_all[, rownames(missing_meta)])
    }
    
    meta$samples <- factor(meta$sample)
    
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

groups <- names(all_cellchat_new_ls)
# cellchat <- mergeCellChat(all_cellchat_new_ls, add.names = groups, cell.prefix = TRUE)
cellchat <- all_cellchat_new_ls$whole_slide

saveRDS(all_cellchat_new_ls, file.path(out_dir_cellchat, "all_cellchat_new_ls.rds")) 
saveRDS(cellchat, file.path(out_dir_cellchat, "cellchat.rds")) 


# all_combns_ls <- groups %>% seq_along() %>% combn(2) %>% apply(2, c, simplify = FALSE)


# ## comprare interation ------
# pdf(
#   file.path(out_dir_cellchat, str_glue("comprare interation plot.pdf")),
#   height = 4,
#   width  = 6
# )
# compareInteractions(cellchat, show.legend = FALSE, measure = "count")
# compareInteractions(cellchat, show.legend = FALSE, measure = "weight")
# dev.off()

# ## diff network plot --------
# all_combns_ls %>%
#   walk(function(x){
#     grp <- groups[x]
#     pdf(
#       file.path(out_dir_cellchat, str_glue("diff network plot - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
#       height = 10,
#       width  = 10
#     )
#     # par(mfrow = c(1,2), xpd=TRUE)
#     netVisual_diffInteraction(
#       cellchat,
#       weight.scale  = T, 
#       comparison    = x,
#       title.name    = "Differential number of interactions",
#       margin        = c(0.2, 0.2, 0.5, 0.2),
#       arrow.size    = 0.6
#     )
#     netVisual_diffInteraction(
#       cellchat,
#       weight.scale  = T,
#       comparison    = x,
#       measure       = "weight", 
#       title.name    = "Differential interaction strength",
#       margin        = c(0.2, 0.2, 0.5, 0.2),
#       arrow.size    = 0.6
#     )
#     dev.off()
#   })
# 
# 
# 
# ## diff heatmap plot --------
# all_combns_ls %>%
#   walk(function(x){
#     grp <- groups[x]
#     
#     pdf(
#       file.path(out_dir_cellchat, str_glue("diff heatmap - {grp[1]} vs {grp[2]} - red high, blue low.pdf")),
#       height = 10,
#       width  = 10
#     )
#     print(netVisual_heatmap(cellchat))
#     print(netVisual_heatmap(cellchat, measure = "weight"))
#     dev.off()
#   })
# 
# 
# ## diff rank plot --------
# all_combns_ls %>%
#   walk(function(x){
#     grp <- groups[x]
#     pdf(
#       file.path(out_dir_cellchat, str_glue("diff rank plot - {grp[1]} vs {grp[2]}.pdf")),
#       height = 6,
#       width  = 6
#     )
#     # par(mfrow = c(1,2), xpd=TRUE)
#     p1 <- rankNet(
#       cellchat,
#       mode        = "comparison",
#       measure     = "count",
#       sources.use = NULL,
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
#       sources.use = NULL,
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
#       sources.use = NULL,
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
#       sources.use = NULL,
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


## diff incoming and outgoing cells --------
# outgoing
pdf(
  file.path(out_dir_cellchat, str_glue("contribution of signals to cell groups.pdf")),
  height = 8, 
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
  height = 8, 
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
  width  = 5
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



# ## diff signal gene -------------------------------------------------------------
# signal_ls <- c("CXCL", "SEMA4", "CDH", "ADGRE", "CD23", "IL16", "COLLAGEN")
# 
# pdf(file.path(out_dir_cellchat, 'T2 signal net.pdf'), width = 10, height = 5)
# signal_ls %>% 
#   map(function(signal){
#     CellChat::netAnalysis_signalingRole_network(
#       all_cellchat_new_ls$T2,
#       signaling = signal,
#       height = 5,
#       width = 15
#     )
#   })
# dev.off()
# 
# pdf(file.path(out_dir_cellchat, 'dotplot_CXCL_by_clinical_type.pdf'), width = 4, height = 5)
# all_data %>%
#   DotPlot(
#     features  = c('CXCR4',	'CXCL12'),
#     group.by = "clinical_type"
#   )
# dev.off()
# 
# pdf(file.path(out_dir_cellchat, 'dotplot_CXCL_by_annotation_full.pdf'), width = 6, height = 8)
# all_data %>%
#   DotPlot(
#     features  = c('CXCR4',	'CXCL12'),
#     group.by = "annotation_full"
#   )
# dev.off()
# 
# 
# p_dot_ls_splitby_clinical <-
#   all_data$clinical_type %>% 
#   unique() %>%
#   setNames(., .) %>%
#   map(function(x){
#     sc <- all_data[, all_data$clinical_type == x]
#     DotPlot(sc, features = c("CXCR4", 'CXCL12'), group.by = "annotation_full") + ggtitle(x)
#   })
# 
# p_dot_ls_splitby_clinical %>%
#   iwalk(function(p, name){ # browser()
#     pdf(file.path(out_dir_cellchat, str_glue('dotplot_CXCL_by_annotation_full_{name}.pdf')), width = 6, height = 8)
#     print(p)
#     dev.off()
#   })
# 
# 
# 
# 
# clinical_type_ls <- c('T1', 'T2', 'T3', 'T4')
# names(clinical_type_ls) <- clinical_type_ls
# 
# p_SEMA4D_by_clinical <- 
#   clinical_type_ls %>% 
#   map(function(clinical_typeV){
#     dat <- subset(all_data, clinical_type == clinical_typeV)
#     
#     DotPlot(dat, features = c('SEMA4D', 'PLXNB2', 'MET', 'ERBB2'), group.by = 'annotation_full')
#   })
# 
# pdf(file.path(out_dir_cellchat, str_glue('dotplot_SEMA4D_by_annotatation6_split.pdf')), width = 8, height = 8)
# p_SEMA4D_by_clinical %>% 
#   iwalk(function(p, name){
#     print(p+ggtitle(name))
#   })
# dev.off()
# 
# pdf(file.path(out_dir_cellchat, str_glue('dotplot_SEMA4D_by_annotatation6_split_all.pdf')), width = 35, height = 8)
# DotPlot(all_data, features = c('SEMA4D', 'PLXNB2', 'MET', 'ERBB2'), group.by = 'annotation_full', split.by = 'clinical_type', cols = c('red', 'green', 'yellow', 'blue')) + 
#   coord_flip() + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
# dev.off()
# 
# 
# 
# 
# all_data_meta <- readRDS(file.path(cn_out_dir, 'all_data_meta_cn.rds'))
# p_col_proj <-
#   all_data_meta %>%
#   group_by(sample, CN2) %>%
#   summarise(count = n()) %>%
#   mutate(freq = count/sum(count)) %>%
#   
#   ggplot(aes(x = sample, y = freq, fill = CN2)) + 
#   geom_col() +
#   ggsci::scale_fill_igv() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   )
# 
# 
# # export
# ggsave(
#   p_col_proj,
#   filename = file.path(cn_out_dir, "CN_freq_sample.pdf"),
#   height = 6,
#   width  = 8
# )
# 
# p_col_proj <-
#   all_data_meta %>%
#   group_by(clinical_type, CN2) %>%
#   summarise(count = n()) %>%
#   mutate(freq = count/sum(count)) %>%
#   
#   ggplot(aes(x = clinical_type, y = freq, fill = CN2)) + 
#   geom_col() +
#   ggsci::scale_fill_igv() +
#   theme(
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
#   )
# 
# 
# # export
# ggsave(
#   p_col_proj,
#   filename = file.path(cn_out_dir, "CN_freq_clinical_type.pdf"),
#   height = 6,
#   width  = 8
# )


# 绘制整体弦图
pdf(file.path(out_dir_cellchat, 'chord plot.pdf'), width = 9, height = 8)
netVisual_chord_gene(cellchat, 
                     sources.use = NULL,  # 所有发送细胞
                     targets.use = NULL,  # 所有接收细胞
                     lab.cex = 0.8,       # 标签大小
                     transparency = 0.4,       # 弦的透明度
                     legend.pos.y = 30)   # 图例位置
dev.off()
