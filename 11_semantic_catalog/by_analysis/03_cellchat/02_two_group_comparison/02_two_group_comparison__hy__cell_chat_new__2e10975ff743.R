rm(list = ls())
library(cowplot)
library(Seurat)
# library(dior)
library(RColorBrewer)

# install.packages('NMF')
# devtools::install_github("jokergoo/circlize")
# devtools::install_github("jokergoo/ComplexHeatmap")
# devtools::install_github("sqjin/CellChat")
library(CellChat)
library(ggplot2)
# BiocManager::install("ggalluvial")
library(ggalluvial)
library(svglite)
library(patchwork)
#library(SingleR)
#library(getopt)
library(CellChat)
options(stringsAsFactors = FALSE)

rm(list = ls())
setwd("/home/hy/Project/FWSP20251095/Redo/subcluster/sub/cellchat/CN6")
Seuset1 = readRDS("/home/hy/Project/FWSP20251095/Redo/subcluster/sub/cellchat/CN6.rds")
selected_types <- c("CD4+ CD2+ T","SOCS3+ M1","SOCS3+ Fibro","SOCS3+ CTMCs","THY1+ FB","SIGLEC1+ Macrophages","MARCO+ Marcophages","FAMs","CTMs_activated","M2")
Seuset<- subset(Seuset1, subset = celltype_final %in% selected_types)
cellchat_sp = SplitObject(Seuset,split.by = "group")

meta1 = cellchat_sp[["disease"]]@meta.data
meta2 = cellchat_sp[["normal"]]@meta.data

cellchat_A <- createCellChat(object = cellchat_sp[["disease"]], meta = meta1, group.by = "celltype_final")
cellchat_B <- createCellChat(object = cellchat_sp[["normal"]], meta = meta2, group.by = "celltype_final")

runCellchat = function(cellchat,name_ct){
  # showDatabaseCategory(CellChatDB)###作者提供了可视化的代码，可以看到该数据库中“Secreted Signaling”占比过半
  cellchat <- setIdent(cellchat, ident.use = "celltype_final") # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
  CellChatDB = CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB,search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")) # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
  cellchat@DB <- CellChatDB.use
  
  unique(cellchat@idents)
  cellchat@data.signaling 
  cellchat <- subsetData(cellchat)
  # future::plan("multiprocess", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
  cellchat <- identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
  
  # 检查当前idents的因子水平
  print(levels(cellchat@idents))

# 删除未使用的因子水平
  cellchat@idents <- droplevels(cellchat@idents)

# 验证清理后的水平
  print(levels(cellchat@idents))
  # project gene expression data onto PPI network (optional)
  #cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat,type="triMean",trim=0,raw.use = TRUE)###这应该是核心的一行代码
  cellchat <- filterCommunication(cellchat, min.cells = 5)##过滤掉小于10个细胞的胞间通讯网络
  cellchat <- computeCommunProbPathway(cellchat)
  
  df.net <- subsetCommunication(cellchat)
  relation <- paste0(df.net[,1],"->",df.net[,2])
  df.net = cbind(relation,df.net)
  colnames(df.net)[1]="relation"
  write.table(df.net,paste0(name_ct,".","net_lr.txt"),sep="\t",quote=F,row.names=F,col.names=T)
  #保存信号通路细胞通讯网络
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat,slot.name="netP")
  write.table(df.netp,paste0(name_ct,".","netPathways_lr.txt"),sep="\t",quote=F,row.names=F,col.names=T)
  
  cellchat <- aggregateNet(cellchat)
  cellchat = netAnalysis_computeCentrality(cellchat)
  return(cellchat)
}
cellchat_A = runCellchat(cellchat_A,"disease")
cellchat_B = runCellchat(cellchat_B,"normal")
object.list = list(disease=cellchat_A,normal=cellchat_B)

#common_celltypes <- Reduce(intersect, list(
#  levels(cellchat_A@idents), 
#  levels(cellchat_B@idents)  # 添加更多对象...
#))
#common_pathways <- Reduce(intersect, list(
#  names(cellchat_A@netP), 
#  names(cellchat_B@netP)  # 若合并多个对象，需列出所有对象
#))
#
#for (obj in list(cellchat_A, cellchat_B)) {
#  obj@idents <- factor(obj@idents, levels = common_celltypes)
#  obj@netP <- obj@netP[common_pathways]
#  obj@netP[common_pathways]  # 调整网络维度
#}

# 1. 统一细胞类型和通路
common_celltypes <- intersect(levels(cellchat_A@idents), levels(cellchat_B@idents))
common_pathways <- intersect(names(cellchat_A@netP), names(cellchat_B@netP))

# 2. 更新对象（强制修改）
cellchat_A@idents <- factor(cellchat_A@idents, levels = common_celltypes)
cellchat_B@idents <- factor(cellchat_B@idents, levels = common_celltypes)
cellchat_A@netP <- cellchat_A@netP[common_pathways]
cellchat_B@netP <- cellchat_B@netP[common_pathways]

# 3. 合并对象（命名必须匹配！）
merged_obj <- mergeCellChat(
  object.list = list(disease = cellchat_A, normal = cellchat_B),
  add.names = c("disease", "normal"),
  cell.prefix = FALSE  # 避免细胞名添加前缀
)


#merged_obj <- mergeCellChat(
#  object.list = list(cellchat_A, cellchat_B),
#  add.names = c("disease", "normal"),  # 命名各组
#  cell.prefix = TRUE                 # 禁用细胞名前缀
#)

cellchat = merged_obj

#cellchat = mergeCellChat(object.list,add.names = names(object.list))

cellchat
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),measure = "count")
gg1
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg2
gg = gg1 + gg2

ggsave("compareInteractions.png",width = 10,height = 6)
dev.off()

pdf("diffInteraction.pdf",width = 10,height = 6)
par(mfrow = c(1,2),xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T,measure = "count")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()


gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg = gg1 + gg2
pdf("heatmap.pdf",width = 10,height = 6)
gg
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()


library(ComplexHeatmap)

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", 
                                        signaling = c("TGFb","CXCL","CCL","THBS","VEGF"),
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing",
                                        signaling = c("TGFb","CXCL","CCL","THBS","VEGF"),
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", 
                                        signaling = c("TGFb","CXCL","CCL","THBS","VEGF"),
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming",
                                        signaling = c("TGFb","CXCL","CCL","THBS","VEGF"),
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

netVisual_bubble(cellchat, sources.use = 3, targets.use = c(6,11), 
                 min.dataset = 2,
                 remove.isolate = T,
                 comparison = c(1, 2), angle.x = 45)

pos.dataset = "disease"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat,
                                       group.dataset = "datasets",
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name,
                                       only.pos = FALSE,
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 0.05)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat,
                              net = net,
                              datasets = "disease",
                              sources.use = 3,
                              ligand.logFC = 0.2,
                              receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat,
                                net = net, 
                                datasets = "normal",
                                sources.use = 3,
                                ligand.logFC = -0.1,
                                receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,
                        sources.use = 3, targets.use = c(6,11),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
                        sources.use = 3, targets.use = c(6,11),
                        comparison = c(1, 2), 
                        angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2



group.cellType <- c(rep("Fibroblast", 3), rep("Monocytic", 3))
group.cellType <- factor(group.cellType, levels = c("Fibroblast", "Monocytic"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 8, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","weight", "weight.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 8, title.name = paste0("weight of interactions - ", names(object.list)[i]))
}
dev.off()


# gg1 = netVisual_heatmap(cellchat_A,color.heatmap = c("red2","white"))
# gg2 = netVisual_heatmap(cellchat_B,color.heatmap = c("red2","white"))
# gg1+gg2
# weight.max = getMaxWeight(object.list, attribute = c("group","count"))
# par(mfrow=c(1,2),xpd=TRUE)
# for (i in 1:2) {
#   netVisual_circle(object.list[[i]]@net$count,weight.scale = T,label.edge = F,edge.weight.max = weight.max[2],
#                    edge.width.max = 12,title.name = paste0("Number of interactions - ",names(object.list)[i]))
# }
# dev.off()
# 
# mat = cellchat@net$disease$count
# par(mfrow=c(3,4),xpd=TRUE)
# for (i in 1:nrow(mat)) {
#   mat2 = matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
#   mat2[i,] = mat[i,]
#   netVisual_circle(mat2,vertex.weight = 11,weight.scale = T,edge.weight.max = max(mat),title.name = rownames(mat)[i])
# }
# dev.off()

