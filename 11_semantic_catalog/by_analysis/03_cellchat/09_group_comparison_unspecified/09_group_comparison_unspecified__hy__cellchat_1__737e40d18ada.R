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
options(stringsAsFactors = FALSE)

outdir <- "/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/core/make_band/cellchat/out"
inrds  <- "/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/core/make_band/cellchat/out/adata_out.rds"
species <- 'human'

setwd(outdir)

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

save(cellchat, file = "cellchat_object.RData")
save(cellchat.objet.list.new, file = "cellchat_object_list.RData")

# 绘图辅助函数
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
gg1
plot_save("Number of inferred interactions.barplot", gg1)

gg2 <- compareInteractions(cellchat, show.legend = FALSE, measure = "weight")
plot_save("Interaction strength.barplot", gg2)
