rm(list = ls())

library(cowplot)

library(Seurat)

library(RColorBrewer)

library(CellChat)

library(ggplot2)

library(ggalluvial)

library(svglite)

library(patchwork)

library(CellChat)

options(stringsAsFactors = FALSE)

outdir = "/mnt/NexData_file/home/hy/Project/Project/liuxiaoxiao/xenium/RCTD/TMA/combine/subset_mpc/cell_chat/out_a_and_in_d"

inrds = "/mnt/NexData_file/home/hy/Project/Project/liuxiaoxiao/xenium/RCTD/TMA/combine/subset_mpc/cell_chat/out_a_and_in_d/adata_out_a_and_in_d.rds"

species = 'human'

# pos.dataset = "iMCD"

# ref.dataset = "NC"







setwd(outdir)




# 检查函数是否存在
#exists("projectData")
# 查找函数属于哪个包
#find("projectData")
#compare.list = factor(c(pos.dataset,ref.dataset),levels = c(ref.dataset,pos.dataset))

Seuset = readRDS(inrds)

cellchat.list = SplitObject(Seuset,split.by = "site_band")

meta.list = lapply(cellchat.list,FUN = function(x){

  return(x@meta.data)

})

group = as.character(unique(Seuset$site_band))

cellchat.objet.list = lapply(group,FUN = function(x){

  createCellChat(cellchat.list[[x]],meta = meta.list[[x]],group.by = "annotation")

})

names(cellchat.objet.list) = group



runCellchat = function(cellchat,species = "human",use.ident = "annotation",filter.cells = 10,name_ct){

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

 # if(species=="human"){

  #  cellchat = projectData(cellchat, PPI.human)

  #}else{

   # cellchat = projectData(cellchat, PPI.mouse)

  #}

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

celltype = unique(cellchat@meta$celltype_final) %>% as.character()

save(cellchat,file = "cellchat_object.RData")

save(cellchat.objet.list.new,file = "cellchat_object_list.RData")





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
