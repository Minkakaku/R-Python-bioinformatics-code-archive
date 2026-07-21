library(cowplot)
library(Seurat)
# library(dior)
library(RColorBrewer)
library(cowplot)
library(RColorBrewer)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(getopt)
options(stringsAsFactors = FALSE)

rm(list = ls())
gc()
setwd('/home/zwl/Project/FW2024-225-G/cellchat/')


# ###########################################exchange h5ad to rds############################################
library(sceasy)
library(Seurat)
library(reticulate)

Seuset <- readRDS(file = './HSPC_redsc.rds')
Seuset <- FindVariableFeatures(Seuset,nfeatures = 2000)
metadata <-read.csv(file = './3.ALL_subtype_metadata.csv',row.names = 1)

Seuset$cell_type = metadata[colnames(Seuset),'sub_celltype2']
Seuset$group = Seuset$disease

Seurat::DimPlot(Seuset,reduction = 'umap',group.by = 'cell_type')



###########################################begin################################################
cellchat_sp = SplitObject(Seuset,split.by = "group")

meta1 = cellchat_sp[["iMCD"]]@meta.data
meta2 = cellchat_sp[["NC"]]@meta.data

cellchat_A <- createCellChat(object = cellchat_sp[["iMCD"]], meta = meta1, group.by = "cell_type")
cellchat_B <- createCellChat(object = cellchat_sp[["NC"]], meta = meta2, group.by = "cell_type")

runCellchat = function(cellchat,name_ct){
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
  cellchat <- identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
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
cellchat_A = runCellchat(cellchat_A,"iMCD")
cellchat_B = runCellchat(cellchat_B,"NC")
object.list = list(iMCD=cellchat_A,NC=cellchat_B)
cellchat = mergeCellChat(object.list,add.names = names(object.list))
save(cellchat,file = "cellchat_object.RData")
save(object.list,file = "cellchat_object_list.RData")


cellchat
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),measure = "count")
gg1
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg2
gg = gg1 + gg2

ggsave("compareInteractions.png",width = 5,height = 4,gg)
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
fig_save = function(filename,plt,dpi=300,w=7,h=7){
  pdf(paste0(filename,'.pdf'),height = h,width = w)
  draw(plt)
  dev.off()
}

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", 
                                        signaling = c("TGFb","CXCL","CCL","THBS","IL6","IFN-I",'IFN-II'),
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing",
                                        signaling = c("TGFb","CXCL","CCL","THBS","IL6","IFN-I",'IFN-II'),
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i+1], width = 5, height = 6)
fig_save("iMCD VS NC_Outcoming compare.heatmap",plt = draw(ht1 + ht2, ht_gap = unit(0.5, "cm")),
         dpi = 300,w=10,h=6)

pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", 
                                        signaling = c("TGFb","CXCL","CCL","THBS","IL6","IFN-I",'IFN-II','IFN-III'),
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming",
                                        signaling = c("TGFb","CXCL","CCL","THBS","IL6","IFN-I",'IFN-II','IFN-III'),
                                        # cluster.rows = TRUE,
                                        # cluster.cols = TRUE,
                                        title = names(object.list)[i+1], width = 5, height = 6)
fig_save("iMCD VS NC_Incoming compare.heatmap",plt = draw(ht1 + ht2, ht_gap = unit(0.5, "cm")),
         dpi = 300,w=10,h=6)


# pos.dataset = "iMCD"
# features.name = pos.dataset
# cellchat <- identifyOverExpressedGenes(cellchat,
#                                        group.dataset = "datasets",
#                                        pos.dataset = pos.dataset,
#                                        features.name = features.name,
#                                        only.pos = FALSE,
#                                        thresh.pc = 0.1,
#                                        thresh.fc = 0.1,
#                                        thresh.p = 0.05)
# net <- netMappingDEG(cellchat, features.name = features.name)
# net.up <- subsetCommunication(cellchat,
#                               net = net,
#                               datasets = "iMCD",
#                               sources.use = 'MEP',
#                               targets.use = c("Megakaryocyte",'Erythrocyte'),
#                               ligand.logFC = 0.2,
#                               receptor.logFC = NULL)
# net.down <- subsetCommunication(cellchat,
#                                 net = net,
#                                 datasets = "iMCD",
#                                 sources.use = 'MEP',
#                                 targets.use = c("Megakaryocyte",'Erythrocyte'),
#                                 ligand.logFC = -0.2,
#                                 receptor.logFC = NULL)
# gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
# gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
# pairLR.use.up = net.up[, "interaction_name", drop = F]
# gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,
#                         sources.use = 3, targets.use = c(6,11),
#                         comparison = c(1, 2),
#                         angle.x = 90,
#                         remove.isolate = T,
#                         title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
# #> Comparing communications on a merged object
# pairLR.use.down = net.down[, "interaction_name", drop = F]
# gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
#                         sources.use = 3, targets.use = c(6,11),
#                         comparison = c(1, 2),
#                         angle.x = 90,
#                         remove.isolate = T,
#                         title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
# #> Comparing communications on a merged object
# gg1 + gg2
# ggsave(filename = './regulated_signaling_in_iMCD.png',w=10,h=8,gg1 + gg2)


gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
ggsave(filename = './relative_information_flow.png',gg1,w=6,h=6)
ggsave(filename = './information_flow.png',gg2,w=6,h=6)
########################################################special pathway#######################################################
p = netVisual_bubble(cellchat, sources.use = 'MEP',
                 targets.use = c("Megakaryocyte",'Erythrocyte'),
                 comparison = c(1, 2), angle.x = 45,
                 title.name = paste0('singnal communication'),
                 thresh = 0.05)
ggsave(filename = 'Singnal_communication_P001_ALL.png',w=6,h=6,p)



################没弄好
infla_gene = readxl::read_xlsx(path = '../05.geneset/inflammatory genes.xlsx') %>% data.frame()
All_gene = NULL
for (i in 1:4) {
  gene = infla_gene[,i][infla_gene[,i] %in% rownames(cellchat@data.signaling)]
  All_gene = c(All_gene,gene)
}

p = netVisual_bubble(cellchat, sources.use = 'MEP',
                     targets.use = c("Megakaryocyte",'Erythrocyte'),
                     comparison = c(1, 2), angle.x = 45,
                     signaling = All_gene,
                     title.name = paste0('singnal communication'),
                     thresh = 0.05)





