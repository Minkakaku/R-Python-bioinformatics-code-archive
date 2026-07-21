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
library('getopt')
library(CellChat)
options(stringsAsFactors = FALSE)

rm(list = ls())
gc()
setwd("/mnt/NexData_file/home/zwl/Project/FWSC20241172/require5/cellchat/B/")
# ###########################################exchange h5ad to rds############################################
library(sceasy)
library(Seurat)
library(reticulate)

Seuset <- readRDS(file = '../../monocle/B_raw.rds') #B细胞数量太少没结果；最好每组每类>50
Seuset <- FindVariableFeatures(Seuset,nfeatures = 2000)
Seuset$cell_type = paste0('C',Seuset$leiden)

Seurat::DimPlot(Seuset,reduction = 'umap',group.by = c('group','cell_type'))
table(Seuset$group,Seuset$cell_type)
#########################################save image function##############################
plot_save = function(filename,plt,dpi=300,w=7,h=7){
  ggsave(paste0(filename,'.pdf'),plt,height = h,width = w)
}
library(ComplexHeatmap)
fig_save = function(filename,plt,dpi=300,w=7,h=7){
  pdf(paste0(filename,'.pdf'),height = h,width = w)
  draw(plt)
  dev.off()
}

cellchat_sp = SplitObject(Seuset,split.by = "group")

meta1 = cellchat_sp[["EXP"]]@meta.data
meta2 = cellchat_sp[["CON"]]@meta.data

cellchat_A <- createCellChat(object = cellchat_sp[["EXP"]], meta = meta1, group.by = "cell_type")
cellchat_B <- createCellChat(object = cellchat_sp[["CON"]], meta = meta2, group.by = "cell_type")

runcellchat = function (cellchat,name_ct){
  # showDatabaseCategory(CellChatDB)###作者提供了可视化的代码，可以看到该数据库中“Secreted Signaling”占比过半
  cellchat <- setIdent(cellchat, ident.use = "cell_type") # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
  CellChatDB = CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB,search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")) # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
  cellchat@DB <- CellChatDB.use
  
  unique(cellchat@idents)
  cellchat@data.signaling 
  cellchat <- subsetData(cellchat)
  # future::plan("multiprocess", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
  cellchat <- identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用 
  # project gene expression data onto PPI network (optional)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat,type="triMean",trim=0,raw.use = TRUE)###这应该是核心的一行代码
  cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络
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
cellchat_A = runcellchat(cellchat_A,"EXP")
cellchat_B = runcellchat(cellchat_B,"CON")

# #CON少一个细胞类型
 group.new = levels(cellchat_A@idents)
 cellchat_B <- liftCellChat(cellchat_B, group.new)

object.list = list(EXP=cellchat_A,CON=cellchat_B)
cellchat = mergeCellChat(object.list,add.names = names(object.list))

save(cellchat,file = "cellchat_object.RData")
save(object.list,file = "cellchat_object_list.RData")

cellchat
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),measure = "count")
gg1
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg2
gg = gg1 + gg2

ggsave("compareInteractions.png",width = 10,height = 6,gg)

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
pdf(file = paste0('cell-cell communication_count_net.pdf'))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()


library(ComplexHeatmap)

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", 
                                        signaling = pathway.union,
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i], width = 6, 
                                        height = length(pathway.union)/4)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing",
                                        signaling = pathway.union,
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i+1], width = 6, 
                                        height = length(pathway.union)/4)
fig_save("EXP VS CON_Outcoming compare.heatmap",plt = draw(ht1 + ht2, ht_gap = unit(0.5, "cm")),
         dpi = 300,w=10,h=16)


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", 
                                        signaling = pathway.union,
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i], width = 6, 
                                        height = length(pathway.union)/4)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming",
                                        signaling = object.list[[i+1]]@netP[["pathways"]],
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i+1], width = 6, 
                                        height = length(pathway.union)/4)
fig_save("EXP VS CON_Incoming compare.heatmap_EXP",plt = ht1,dpi = 300,w=10,h=16)
fig_save("EXP VS CON_Incoming compare.heatmap_CON",plt = ht2,dpi = 300,w=10,h=16)


pos.dataset = "EXP"
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
                              datasets = "EXP",
                              #sources.use = 3,
                              ligand.logFC = 0.2,
                              receptor.logFC = 0.2)
net.down <- subsetCommunication(cellchat,
                                net = net, 
                                datasets = "CON",
                                #sources.use = 3,
                                ligand.logFC = -0.1,
                                receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,
                        #sources.use = 3, targets.use = c(6,11),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
                        #sources.use = 3, targets.use = c(6,11),
                        comparison = c(1, 2), 
                        angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2
ggsave(filename = './regulated_up_signaling_in_EXP.png',w=12,h=8,gg1,dpi = 200)
ggsave(filename = './regulated_down_signaling_in_EXP.png',w=12,h=8,gg2)


gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
plot_save(paste0("EXP VS CON_Relative information flow.weight stacked"),gg1,w=6,h=10)
plot_save(paste0("EXP VS CON_Relative information flow.weight"),gg2,w=6,h=10)

gg1 <- rankNet(cellchat, mode = "comparison", measure = "count", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "count", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
plot_save(paste0("EXP VS CON_Relative information flow.count stacked"),gg1,w=6,h=10)
plot_save(paste0("EXP VS CON_Relative information flow.count"),gg2,w=6,h=10)

