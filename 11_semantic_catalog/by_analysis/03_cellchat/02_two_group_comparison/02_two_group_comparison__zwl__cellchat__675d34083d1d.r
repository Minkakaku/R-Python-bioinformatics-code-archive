rm(list = ls())
gc()
library(cowplot)
library(Seurat)
library(RColorBrewer)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(getopt)
library(CellChat)
options(stringsAsFactors = FALSE)
outdir = "/mnt/NexData_file/home/zwl/Project/FWSC20251792/require1/cellchat/"
setwd(outdir)

#read rds########################################
library(sceasy)
library(Seurat)
library(reticulate)

sceasy::convertFormat('../raw.h5ad', from="anndata", to="seurat",outFile='raw.rds')
Seuset <- readRDS(file = "/mnt/NexData_file/home/zwl/Project/FWSC20251792/require1/cellchat/raw.rds")
Seuset <- NormalizeData(Seuset,normalization.method = "LogNormalize",scale.factor = 1e4) 
Seuset <- FindVariableFeatures(Seuset)
Seuset <- ScaleData(Seuset)
Seurat::DimPlot(Seuset,reduction = 'umap',group.by = c('cell_type','group'))

pos.dataset1 = "EDA1"
ref.dataset = "Control1"
species = 'mouse'



# begin##############################################

compare.list = factor(c(pos.dataset1,ref.dataset),
                      levels = c(ref.dataset,pos.dataset1))

cellchat.list = SplitObject(Seuset,split.by = "group")
meta.list = lapply(cellchat.list,FUN = function(x){
  return(x@meta.data)
})
group = as.character(unique(Seuset$group))
cellchat.objet.list = lapply(group,FUN = function(x){
  createCellChat(cellchat.list[[x]],meta = meta.list[[x]],group.by = "cell_type")
})
names(cellchat.objet.list) = group

runCellchat = function(cellchat,species = "human",use.ident = "cell_type",filter.cells = 10,name_ct){
  print(name_ct)
  cellchat = setIdent(cellchat, ident.use = use.ident) # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
  if (species=="human") {
    CellChatDB = CellChatDB.human  
  }else{
    CellChatDB = CellChatDB.mouse  
  }
  write.table(CellChatDB$interaction,'cellchatDB.xls',sep = "\t")
  showDatabaseCategory(CellChatDB)
  # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
  cellchat@DB = subsetDB(CellChatDB,search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")) 
  cellchat = subsetData(cellchat)
  # future::plan("multiprocess", workers = 4)
  cellchat = identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
  cellchat = identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
  # project gene expression data onto PPI network (optional)
  if(species=="human"){
    cellchat = projectData(cellchat, PPI.human)
  }else{
    cellchat = projectData(cellchat, PPI.mouse)
  }
  # group.new = levels(cellchat@idents)
  # cell.chat = liftCellChat(cellchat,group.new)
  cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))
  cellchat = computeCommunProb(cellchat,type="triMean",trim = 0)###这应该是核心的一行代码
  cellchat = filterCommunication(cellchat, min.cells = filter.cells)##过滤掉小于10个细胞的胞间通讯网络
  cellchat = computeCommunProbPathway(cellchat)
  df.net = subsetCommunication(cellchat)
  relation = paste0(df.net[,1],"->",df.net[,2])
  df.net = cbind(relation,df.net)
  colnames(df.net)[1]="relation"
  write.table(df.net,paste0(name_ct,".","net_lr.txt"),sep="\t",quote=F,row.names=F,col.names=T)
  #保存信号通路细胞通讯网络
  cellchat = computeCommunProbPathway(cellchat)
  df.netp = subsetCommunication(cellchat,slot.name="netP")
  write.table(df.netp,paste0(name_ct,".","netPathways_lr.txt"),sep="\t",quote=F,row.names=F,col.names=T)
  cellchat = aggregateNet(cellchat)
  cellchat = netAnalysis_computeCentrality(cellchat)
  return(cellchat)
}

cellchat.objet.list.new = lapply(group,FUN = function(x){
  runCellchat(cellchat.objet.list[[x]],species = species,name_ct=x)
})

names(cellchat.objet.list.new) = group
idents.list = lapply(cellchat.objet.list.new,FUN = function(x){
  length(levels(x@idents))
})


# 这里判断是否包含的细胞类型数目一致不一致的话选择里面数量最多的，其他的看齐
if (length(unique(unlist(idents.list)))>1) {
  maxgroup = names(idents.list)[idents.list == max(unique(unlist(idents.list)))]
  group.max = levels(cellchat.objet.list.new[[maxgroup]]@idents)
  cellchat.objet.list.new = lapply(cellchat.objet.list.new, function(x){
    liftCellChat(x, group.max)
  })
  
}

cellchat = mergeCellChat(cellchat.objet.list.new,add.names = names(cellchat.objet.list.new))
celltype = unique(cellchat@meta$cell_type) %>% as.character()
save(cellchat,file = "cellchat_object.RData")
save(cellchat.objet.list.new,file = "cellchat_object_list.RData")

object.list = cellchat.objet.list.new
########################四个类器官的细胞-细胞通讯网络的总数和相互作用强度柱状图#######################################
plot_save = function(filename,plt,dpi=300,w=7,h=7){
  ggsave(paste0(filename,'.pdf'),plt,height = h,width = w)
  
}
fig_save = function(filename,plt,dpi=300,w=7,h=7){
  pdf(paste0(filename,'.pdf'),height = h,width = w)
  draw(plt)
  dev.off()
}
recor_save = function(filename,plt,dpi=300,w=7,h=7){
  pdf(paste0(filename,'.pdf'),height = h,width = w)
  replayPlot(plt)
  dev.off()
}



################细胞间通信网络图和热图
compare.list = c(2,1)
cellchat
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),measure = "count")
gg1
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg2
gg = gg1 + gg2

ggsave("compareInteractions.png",width = 10,height = 6,gg)
dev.off()

pdf("diffInteraction.pdf",width = 10,height = 6)
par(mfrow = c(1,2),xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T,comparison = compare.list,measure = "count")
netVisual_diffInteraction(cellchat, weight.scale = T,comparison = compare.list, measure = "weight")
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
fig_save("EDA VS Control_Outcoming compare.heatmap",plt = draw(ht1 + ht2, ht_gap = unit(0.5, "cm")),
         dpi = 300,w=10,h=16)


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", 
                                        signaling = pathway.union,
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i], width = 6, 
                                        height = length(pathway.union)/4)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming",
                                        signaling = pathway.union,
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i+1], width = 6, 
                                        height = length(pathway.union)/4)
fig_save("EDA VS Control_Incoming compare.heatmap",plt = draw(ht1 + ht2, ht_gap = unit(0.5, "cm")),
         dpi = 300,w=10,h=16)



pos.dataset = "EDA1"
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
                              datasets = pos.dataset,
                              #sources.use = 3,
                              ligand.logFC = 0.2,
                              receptor.logFC = 0.2)
net.down <- subsetCommunication(cellchat,
                                net = net, 
                                datasets = pos.dataset,
                                #sources.use = 3,
                                ligand.logFC = -0.1,
                                receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,
                        sources.use = 'Neutrophil',
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
                        sources.use = 'Neutrophil',
                        comparison = c(1, 2), 
                        angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2
ggsave(filename = './regulated_up_signaling_in_EDA1.png',w=12,h=8,gg1)
ggsave(filename = './regulated_down_signaling_in_EDA1.png',w=12,h=8,gg2)


gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
plot_save(paste0("EDA VS Control_Relative information flow.weight stacked"),gg1,w=6,h=10)
plot_save(paste0("EDA VS Control_Relative information flow.weight"),gg2,w=6,h=10)

gg1 <- rankNet(cellchat, mode = "comparison", measure = "count", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "count", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
plot_save(paste0("EDA VS Control_Relative information flow.count stacked"),gg1,w=6,h=10)
plot_save(paste0("EDA VS Control_Relative information flow.count"),gg2,w=6,h=10)


pathways.show <- c('IFN-II') #TGFb,'PD-L1'抑制
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, 
                               color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
