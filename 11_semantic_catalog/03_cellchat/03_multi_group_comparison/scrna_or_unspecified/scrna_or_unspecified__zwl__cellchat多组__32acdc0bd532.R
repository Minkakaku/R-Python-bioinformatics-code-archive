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
outdir = "/mnt/NexData_file/home/zwl/Project/FWAE20250109/require3/cellchat"
species = 'human'

pos.dataset2 = 'Elderly'
pos.dataset1 = "Middle"
ref.dataset = "Young"



setwd(outdir)


compare.list = factor(c(pos.dataset2,pos.dataset1,ref.dataset),
                      levels = c(ref.dataset,pos.dataset1,pos.dataset2))

#read rds########################################
Seuset <- readRDS(file = '/mnt/NexData_file/home/zwl/Project/FWAE20250109/rds_h5ad/Main/P20121501.diff_PRO.rds')
Seuset <- FindVariableFeatures(Seuset,nfeatures = 2000)
Seuset$type = Seuset$group
Seuset$group = do::Replace(Seuset$sample,pattern = c("B1_1:Young","B1_2:Young","B2_1:Young","B2_2:Young",
                                                     "B3_1:Middle","B3_2:Middle","B4_1:Middle","B4_2:Middle",
                                                     "B5_1:Elderly","B5_2:Elderly","B6_1:Elderly",
                                                     "P1_1:Young","P1_2:Young","P2_1:Young","P2_2:Young",
                                                     "P3_1:Middle","P3_2:Middle","P4_1:Middle","P4_2:Middle",
                                                     "P5_1:Elderly","P5_2:Elderly","P6_1:Elderly"))
Seurat::DimPlot(Seuset,reduction = 'umap',group.by = c('cluster','group','type'))

Seuset$cell_type = Seuset$cluster

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
compare.list = c(2,1)
gg1 = compareInteractions(cellchat, show.legend = F,measure = "count")
gg1
plot_save(filename = "Number of inferred interactions.barplot",plt = gg1)
gg2 = compareInteractions(cellchat, show.legend = F, measure = "weight")
plot_save("Interaction strength.barplot",gg2)


################control-HCT/control-HT/control-OHP/OHP-HCT/HT-HCT 5种两类器官的细胞间通信网络图和热图
library(ComplexHeatmap)
for (i in 1:length(group)) {
  for (j in (1:length(group))[(1:length(group)) !=i]) {
    compare.list = c(i,j)
    
    # cellchat@netP$
    gg3 = netVisual_diffInteraction(cellchat,  comparison =compare.list,
                                    weight.scale = T,measure = "count")
    dev.off()
    recor_save(paste0(group[i],"VS",group[j],"_Differential number of interactions.net"),gg3)
    gg4 = netVisual_diffInteraction(cellchat, comparison = compare.list,
                                    weight.scale = T, measure = "weight")
    recor_save(paste0(group[i],"VS",group[j],"_Differential interaction strength.net"),gg4)
    gg5 = netVisual_heatmap(cellchat, comparison = compare.list)
    fig_save(paste0(group[i],"VS",group[j],"_Differential number of interactions.heatmap"),gg5)
    gg6 = netVisual_heatmap(cellchat, measure = "weight", comparison = compare.list)
    fig_save(paste0(group[i],"VS",group[j],"_Differential interaction strength.heatmap"),gg6)
    
    gg7 = rankNet(cellchat,
                  mode = "comparison",
                  measure = "weight",
                  sources.use = NULL,
                  targets.use = NULL, comparison =compare.list ,
                  stacked = T,
                  do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Relative information flow.weight"),gg7)
    gg8 = rankNet(cellchat,
                  mode = "comparison",
                  measure = "weight",
                  sources.use = NULL,
                  targets.use = NULL, comparison = compare.list,
                  stacked = F,
                  do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Information flow.weight"),gg8)
    
    gg9 = rankNet(cellchat,
                  mode = "comparison",
                  measure = "count",
                  sources.use = NULL,
                  targets.use = NULL,
                  comparison =compare.list,
                  stacked = T,
                  do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Relative information flow.count"),gg9)
    gg10 = rankNet(cellchat,
                   mode = "comparison",
                   measure = "count",
                   sources.use = NULL,
                   targets.use = NULL, comparison = compare.list,
                   stacked = F,
                   do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Information flow.count"),gg10)
    
    
    
    dif_pathway = as.character(gg8$data[c(1:3,(nrow(gg8$data)-3):nrow(gg8$data)),1])
    
    # 比较群体间的传入传出区别
    library(ComplexHeatmap)
    pathway.list = lapply(cellchat@netP, function(x){
      unique(x$pathways)
    }) %>% unlist() %>% unique()
    # group.max = 
    ht1 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "outgoing", 
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[i], 
                                            width = 10, height = length(pathway.list)/3)
    ht2 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "outgoing",
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[j],
                                            width = 10, height = length(pathway.list)/3)
    fig_save(paste0(group[i],"VS",group[j],'_Outgoing compare.heatmap'),plt = draw(ht1 + ht2, ht_gap = unit(0.5, "cm")),
             dpi = 300,w=10,h=length(pathway.list)/3)
    
    ht3 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "incoming", 
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[i], 
                                            width = 10, height = length(pathway.list)/3)
    ht4 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "incoming",
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[j], 
                                            width = 10, height = length(pathway.list)/3)
    fig_save(paste0(group[i],"VS",group[j],'_Incoming compare.heatmap'),plt = draw(ht3 + ht4, ht_gap = unit(0.5, "cm")),
             dpi = 300,w=10,h=length(pathway.list)/3)
    
    ht5 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "all", 
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[i], 
                                            width = 10, height = length(pathway.list)/3)
    ht6 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "all",
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[j], 
                                            width = 10, height = length(pathway.list)/3)
    fig_save(paste0(group[i],"VS",group[j],'_All compare.heatmap'),plt = draw(ht5 + ht6, ht_gap = unit(0.5, "cm")),
             dpi = 300,w=10,h=length(pathway.list)/3)
    
    
  }
}





##### control-HCT/control-HT/control-OHP/OHP-HCT/HT-HCT 5种两类器官的配体-受体对分析
# 特异上调的下调的配受体关系对
for (i in 1:3) {
  compare.list_ALL = list(c('Young','Middle'),c('Young','Elderly'),c('Middle','Elderly'))
  compare.list = compare.list_ALL[[i]]
  subcellchat = mergeCellChat(cellchat.objet.list.new[compare.list],
                              add.names = names(cellchat.objet.list.new[compare.list]))
  features.name = compare.list[2]
  cellchat_dif = identifyOverExpressedGenes(subcellchat,
                                            group.dataset = "datasets",
                                            pos.dataset = compare.list[2], 
                                            features.name = features.name,
                                            only.pos = FALSE,
                                            thresh.pc = 0.1, 
                                            thresh.fc = 0.1, 
                                            thresh.p = 0.05)
  net = netMappingDEG(cellchat_dif, features.name = features.name)
  net.up = subsetCommunication(cellchat_dif,
                               net = net,
                               datasets = compare.list[2],
                               # sources.use = 3,
                               ligand.logFC = 0.1,
                               receptor.logFC = 0.1)
  
  gene.up = extractGeneSubsetFromPair(net.up, cellchat_dif)
  pairLR.use.up = net.up[, "interaction_name", drop = F]
  gg1 <- netVisual_bubble(cellchat_dif, pairLR.use = pairLR.use.up,
                          # sources.use = 3, targets.use = c(6,11),
                          comparison = c(1, 2),  
                          angle.x = 90, 
                          remove.isolate = F,
                          title.name = paste0("Up-regulated signaling in ", compare.list[2]))
  ggsave(filename = paste0("Up-regulated sig signaling in ", compare.list[1],' vs ',compare.list[2],'.png'),gg1,dpi = 200,w=15,h=8)
  #> Comparing communications on a merged object
  
  tryCatch(
    {
      net.down = subsetCommunication(cellchat_dif,
                                     net = net, 
                                     datasets = ref.dataset,
                                     # sources.use = 3,
                                     ligand.logFC = -0.1,
                                     receptor.logFC = -0.1)
      gene.down = extractGeneSubsetFromPair(net.down, cellchat_dif)
      pairLR.use.down = net.down[, "interaction_name", drop = F]
      gg2 <- netVisual_bubble(cellchat_dif, pairLR.use = pairLR.use.down,
                              comparison = c(1, 2), 
                              angle.x = 90, 
                              remove.isolate = FALSE,
                              title.name = paste0("Down-regulated signaling in ", compare.list[2]))
      ggsave(filename = paste0("Down-regulated sig signaling in ", compare.list[1],' vs ',compare.list[2],'.png'),gg2,dpi = 200,w=15,h=8)
      return(gg2)},
    error = function(e) {
      message(paste("netVisual_bubble Error：", conditionMessage(e)))
      return(NULL)
    }
  )
}




######################4个类器官单独的细胞间通信网络图和热图###############################
for (i in c('Young','Middle','Elderly')){
  subcellchat = cellchat.objet.list.new[[i]]
  subcellchat <- aggregateNet(subcellchat)
  groupSize <- as.numeric(table(subcellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  pdf(file = paste0('Sample_',i,' cell-cell communication_net.pdf'))
  netVisual_circle(subcellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  dev.off()
  p2 = netVisual_heatmap(subcellchat, color.heatmap = "Reds",title.name = "Number of interactions")
  fig_save(filename = paste0('Sample_',i,' cell-cell communication_heatmap.pdf'),p2)
}

######################4个类器官单独的配体-受体对分析###############################

for (i in c('Young','Middle','Elderly')){
  subcellchat = cellchat.objet.list.new[[i]]
  pathways.show = subcellchat@netP$pathways
  
  pdf(file = paste0('Sample_',i,' contribution of each ligand-receptor.pdf'),h=15,w=5)
  netAnalysis_contribution(subcellchat, signaling = pathways.show)
  dev.off()
}





