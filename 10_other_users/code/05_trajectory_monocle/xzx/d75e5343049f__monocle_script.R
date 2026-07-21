
setwd("/mnt/NexData_file/home/xzx/project/AE2604082792-01/shouhou2")

rm(list = ls())

dir.create("./monocle")
setwd('/mnt/NexData_file/home/xzx/project/AE2604082792-01/shouhou2/monocle')

library(monocle)
library(Seurat)
library(ggplot2)
#library(argparse)
library(dplyr)
library(tidyverse)
library(ggsci)
library(patchwork)
library(stringr)
library(getopt)


#读取参数和文件

rdsfile ="/mnt/NexData_file/home/xzx/project/AE2604082792-01/shouhou2/Epithelium/celltype_fin/adata_epi_anno_count.rds" 
geneFile = NULL                        #/t分割符的文本
lambda = 8                             #数值
RootState = NULL                       #起点位置state数值
Reverse = NULL                           #逻辑值
celllist = NULL
cluster = NULL
RandomCells = NULL
cellsPercent = 0.1
minCellsPerCluster = 1000
outpath = "./"
dir.create(paste0(outpath,"Pseudotime"),recursive = T)
setwd(paste0(outpath,"Pseudotime"))

# 函数名: Psepre
# 参数: seuset - Seurat对象 celllists -抽取的细胞列表 cluster -选择细胞类型 RandomCells -随机抽取
# 返回值: seurat - Seurat对象
# 功能: 对传入的seuset进行处理并返回
Psepre= function(seuset,celllist,cluster,RandomCells,cellsPercent,minCellsPerCluster){
  # if(!is.null(celllist)){
  #   newcluster = read.table(celllist,sep="\t",header=F)
  #   realcluster = newcluster[is.finite(match(newcluster[,1],colnames(seuset))),]
  #   seuset = subset(seuset,cells=as.character(realcluster[,1]))
  #   seuset = SetIdent(object = seuset, cells= as.character(realcluster[,1]),
  #                     as.character(realcluster[,2]))
  # }else if(!is.null(cluster)){
  #   clusters = unlist(strsplit(cluster, "[,]"))
  #   celllist = names(seuset@active.ident[is.finite(match(seuset@active.ident,clusters))])
  #   seuset = subset(seuset,cells=celllist)
  # }
  if(!is.null(RandomCells)){
    newcelllist = randomSample(seuset@active.ident,cellsPercent = 0.1,minCellsPerCluster = 1000)
    seuset = subset(seuset,cells=as.character(newcelllist))
    write.table(data.frame(Cell = colnames(seuset),Cluster=seuset@active.ident),
                file=paste0("RandomCells.txt"),sep="\t",row.names=F,quot=F)
  }
  seuset <- AddMetaData(object = seuset, metadata = seuset@active.ident, col.name = "Cluster")
  return(seuset)
}

# 函数名: setNum
# 参数: number - 传入细胞数 cellsPercent -抽取的比例 minCellsPerCluster -最小细胞数
# 返回值: finalNum - 需要抽取的细胞数量
# 功能: 对传入的细胞数以及参数进行处理并返回最终抽取细胞数量
setNum = function(number,cellsPercent,minCellsPerCluster){
  selectNum = round(number*cellsPercent)
  if(selectNum>=minCellsPerCluster){
    finalNum = selectNum
  }else{
    finalNum = min(number,minCellsPerCluster)
  }
  return(finalNum)
}
# 函数名: setNum
# 参数: cluster - 细胞类型定义 cellsPercent -抽取的比例 minCellsPerCluster -最小细胞数
# 返回值: finalNum - 需要抽取的细胞数量
# 功能: 对传入的细胞数以及参数进行处理并返回最终抽取细胞数量
randomSample = function(cluster,size=NULL,cellsPercent,minCellsPerCluster){
  samples = c()
  id = unique(cluster)
  for(i in id){
    clusterInfo = cluster[cluster==i]
    if(is.null(size)){
      finalNum = setNum(length(clusterInfo),cellsPercent,minCellsPerCluster)
    }else{
      finalNum = size
    }
    insample = sample(names(clusterInfo),finalNum,replace = FALSE)
    samples = c(samples,insample)
  }
  return(samples)
}
# 函数名: ImportSeurat4CDS
# 参数: otherCDS - Seurat对象
# 返回值: monocle_cds - 转换后的Monocle对象
# 功能: 将Seurat对象转换为Monocle对象，以便在Monocle中进行进一步的分析和可视化。
ImportSeurat4CDS = function(otherCDS){
  # 获取原始数据，并判断是否为稀疏矩阵，不是则进行转换
  # 官方文档推荐加快计算速度
  # otherCDS = seuset
  data = GetAssayData(otherCDS, slot = "counts", assay = "RNA")
  if (class(data) == "data.frame") {
    data = as(as.matrix(data),"sparseMatrix")
  }
  # 获取meta.data,用于monocle定义数据
  pd = tryCatch(
    {
      pd = new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    },error = function(e){
      pData = data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd = new("AnnotatedDataFrame", data = pData)
      message("This Seurat object doesn't provide any meta data")
      pd
    }
  )
  # raw meta不一致以meta为主，用于monocle的featureData
  if (length(setdiff(colnames(data),rownames(pd))) >0) {
    data = data[,rownames(pd)]
  }
  fd = data.frame(gene_short_name = row.names(data),row.names = row.names(data))
  fd = new("AnnotatedDataFrame",data = fd)
  # 这里使用VGAM包的函数定义传入数据的分布类型，其中negbinomial.size()以及negbinomial()适用于UMI数据前者速度更快
  # 通常情况下将会运行第一种，FPKM值适用tobit
  # if (all(data == floor(data))) {
  expressionFamily <- negbinomial.size()
  # }else if(any(data < 0)){
  #   expressionFamily <- uninormal()
  # }else{
  # expressionFamily <- tobit()
  # }
  # 创建monocle使用object，expressionFamily表达相应的VGAM函数
  monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd,
                                lowerDetectionLimit = 0, expressionFamily = expressionFamily)
  # 如果传入的Seurat包含monocle的数据进行清除
  if ("Monocle" %in% names(otherCDS@misc)) {
    otherCDS@misc$Monocle@auxClusteringData$seurat = NULL
    otherCDS@misc$Monocle@auxClusteringData$scran = NULL
    monocle_cds = otherCDS@misc$Monocle
    mist_list = otherCDS
  }else{
    mist_list = otherCDS
  }
  # 其他信息存储
  monocle_cds@auxClusteringData$seurat <- mist_list
  return(monocle_cds)
}
# 函数名: PsePlot
# 参数: CDS - monocle_cds对象
# 返回值: 画图
# 功能: 画各种图
PsePlot = function(seuCDS){
  gg =  plot_cell_trajectory(seuCDS, color_by = "State")+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.State.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.State.pdf"),limitsize = FALSE)
  # gg =  plot_cell_trajectory(seuCDS, color_by = "Cluster")+
  #   theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
  #         axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  # ggsave(paste0("Pseudotime.Cluster.png"),limitsize = FALSE)
  # ggsave(paste0("Pseudotime.Cluster.pdf"),limitsize = FALSE)
  gg =  plot_cell_trajectory(seuCDS, color_by = "Pseudotime")+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.Pseudotime.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.Pseudotime.pdf"),limitsize = FALSE)
  # gg =  plot_cell_trajectory(seuCDS, color_by = "Cluster") + 
  #   facet_wrap(~Cluster, nrow = floor(sqrt(length(unique(seuCDS$Cluster)))))+
  #   theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
  #         axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  # ggsave(paste0("Pseudotime.ClusterEach.png"),limitsize = FALSE)
  # ggsave(paste0("Pseudotime.ClusterEach.pdf"),limitsize = FALSE)
  gg =  plot_cell_trajectory(seuCDS, color_by = "State") +
    facet_wrap(~State, nrow = floor(sqrt(length(unique(seuCDS$State)))))+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.StateEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.StateEach.pdf"),limitsize = FALSE)
  # 增加的特定列的绘图
  gg =  plot_cell_trajectory(seuCDS, color_by = "cell_type")+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.cell_type.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.cell_type.pdf"),limitsize = FALSE)
  
  gg =  plot_cell_trajectory(seuCDS, color_by = "inferCNV")+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.inferCNV.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.inferCNV.pdf"),limitsize = FALSE)
  
  gg =  plot_cell_trajectory(seuCDS, color_by = "leiden")+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.leiden.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.leiden.pdf"),limitsize = FALSE)
  
  gg =  plot_cell_trajectory(seuCDS, color_by = "cell_type") +
    facet_wrap(~cell_type, nrow = floor(sqrt(length(unique(seuCDS$cell_type)))))+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.CellTypeEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.CellTypeEach.pdf"),limitsize = FALSE)
  
  gg =  plot_cell_trajectory(seuCDS, color_by = 'inferCNV') +
    facet_wrap(~inferCNV, nrow = floor(sqrt(length(unique(seuCDS$inferCNV)))))+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.inferCNVEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.inferCNVEach.pdf"),limitsize = FALSE)
  
  gg =  plot_cell_trajectory(seuCDS, color_by = 'leiden') +
    facet_wrap(~leiden, nrow = floor(sqrt(length(unique(seuCDS$leiden)))))+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime.leidenEach.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime.leidenEach.pdf"),limitsize = FALSE)
  
  
  
  gg = plot_complex_cell_trajectory(seuCDS, color_by = 'inferCNV', 
                                    show_branch_points = T, cell_size = 0.5, 
                                    cell_link_size = 0.3)+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime_inferCNV.tree.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime_inferCNV.tree.pdf"),limitsize = FALSE)
  gg = plot_complex_cell_trajectory(seuCDS, color_by = 'cell_type', 
                                    show_branch_points = T, cell_size = 0.5, 
                                    cell_link_size = 0.3)+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime_Celltype.tree.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime_Celltype.tree.pdf"),limitsize = FALSE)
  gg = plot_complex_cell_trajectory(seuCDS, color_by = 'leiden', 
                                    show_branch_points = T, cell_size = 0.5, 
                                    cell_link_size = 0.3)+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime_leiden.tree.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime_leiden.tree.pdf"),limitsize = FALSE)
  gg = plot_complex_cell_trajectory(seuCDS, color_by = 'State', 
                                    show_branch_points = T, cell_size = 0.5, 
                                    cell_link_size = 0.3)+
    theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
  ggsave(paste0("Pseudotime_state.tree.png"),limitsize = FALSE)
  ggsave(paste0("Pseudotime_state.tree.pdf"),limitsize = FALSE)
}
# 函数名: PsePlotSample
# 参数: CDS - monocle_cds对象
# 返回值: 画样本图
# 功能: 画样本图
PsePlotSample = function(seuCDS){
  if(("sample"%in% colnames(seuCDS@phenoData@data))|("orig.ident"%in% colnames(seuCDS@phenoData@data))){
    if("orig.ident"%in% colnames(seuCDS@phenoData@data)){
      seuCDS$Sample = as.factor(seuCDS$orig.ident)
    }else {
      seuCDS$Sample = as.factor(seuCDS$sample)     
    }
    gg =  plot_cell_trajectory(seuCDS, color_by = "Sample")+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
    ggsave(paste0("Pseudotime.Sample.png"),limitsize = FALSE)
    ggsave(paste0("Pseudotime.Sample.pdf"),limitsize = FALSE)
  }
  if("group" %in% colnames(seuCDS@phenoData@data)){
    gg =  plot_cell_trajectory(seuCDS, color_by = "group")+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
    ggsave(paste0("Pseudotime.group.png"),limitsize = FALSE)
    ggsave(paste0("Pseudotime.group.pdf"),limitsize = FALSE)
  }
  if("Sample"%in% colnames(seuCDS@phenoData@data)){
    gg =  plot_cell_trajectory(seuCDS, color_by = "Sample") + facet_wrap(~Sample, nrow = floor(sqrt(length(unique(seuCDS$Sample)))))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
    ggsave(paste0("Pseudotime.SampleEach.png"),limitsize = FALSE)
    ggsave(paste0("Pseudotime.SampleEach.pdf"),limitsize = FALSE)
  }
  if("group" %in% colnames(seuCDS@phenoData@data)){
    gg =  plot_cell_trajectory(seuCDS, color_by = "group") + facet_wrap(~group, nrow = floor(sqrt(length(unique(seuCDS$group)))))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
    ggsave(paste0("Pseudotime.groupEach.png"),limitsize = FALSE)
    ggsave(paste0("Pseudotime.groupEach.pdf"),limitsize = FALSE)
  }
  if("Sample"%in% colnames(seuCDS@phenoData@data)){
    gg = plot_complex_cell_trajectory(seuCDS, color_by = 'Sample', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Sample, nrow = floor(sqrt(length(unique(seuCDS$Sample))))) + scale_size(range = c(0.2, 0.2))+theme(legend.position="right",axis.text.x = element_text(angle = 30, hjust = 1))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
    ggsave(paste0("SampleEach.tree.png"),limitsize = FALSE)
    ggsave(paste0("SampleEach.tree.pdf"),limitsize = FALSE)
  }
  if("group" %in% colnames(seuCDS@phenoData@data)){
    gg = plot_complex_cell_trajectory(seuCDS, color_by = 'group', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~group, nrow = floor(sqrt(length(unique(seuCDS$Sample))))) + scale_size(range = c(0.2, 0.2))+theme(legend.position="right",axis.text.x = element_text(angle = 30, hjust = 1))+theme(legend.position="right",axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
    ggsave(paste0("groupEach.tree.png"),limitsize = FALSE)
    ggsave(paste0("groupEach.tree.pdf"),limitsize = FALSE)
  }
}

PseDenPlot = function(forplot){
  if("Sample" %in%colnames(forplot)){
    if(length(unique(forplot$Sample))>1){
      nrow = floor(sqrt(length(unique(forplot$Sample))))
      ncol = ceiling(length(unique(forplot$Sample))/nrow)
      gg = ggplot(forplot,aes(x=Pseudotime,fill=Sample))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~Sample,nrow = nrow)
      ggsave("PseudotimeDensity.SampleEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE,bg = 'white')
      ggsave("PseudotimeDensity.SampleEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
      gg = ggplot(forplot,aes(x=Pseudotime,fill=Sample))+geom_density(alpha=0.5)+theme_cowplot()
      ggsave("PseudotimeDensity.Sample.png",gg,width= 8,height=5,limitsize = FALSE,bg = 'white')
      ggsave("PseudotimeDensity.Sample.pdf",gg,width= 8,height=5,limitsize = FALSE)
    }else{
      # plotSampledata=ggplot_build(gg)$data[[1]]
      # write.table(plotSampledata,file="../DCPlotForVis/orig.ident.DCPlot.txt",sep="\t",row.names=F,quot=F)
      nrow = floor(sqrt(length(unique(forplot$leiden))))
      ncol = ceiling(length(unique(forplot$leiden))/nrow)
      gg = ggplot(forplot,aes(x=Pseudotime,fill=leiden))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~Cluster,nrow = nrow)
      ggsave("PseudotimeDensity.leidenEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE,bg = 'white')
      ggsave("PseudotimeDensity.leidenEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
      gg = ggplot(forplot,aes(x=Pseudotime,fill=leiden))+geom_density(alpha=0.5)+theme_cowplot()
      ggsave("PseudotimeDensity.leiden.png",gg,width= 8,height=5,limitsize = FALSE,bg = 'white')
      ggsave("PseudotimeDensity.leiden.pdf",gg,width= 8,height=5,limitsize = FALSE)
      # plotdata=ggplot_build(gg)$data[[1]]
      # write.table(plotdata,file="../DCPlotForVis/Cluster.DCPlot.txt",sep="\t",row.names=F,quot=F)   
    }
  }
  if("group" %in%colnames(forplot)){
    nrow = floor(sqrt(length(unique(forplot$group))))
    ncol = ceiling(length(unique(forplot$group))/nrow)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=group))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~group,nrow = nrow)
    ggsave("PseudotimeDensity.groupEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE,bg = 'white')
    ggsave("PseudotimeDensity.groupEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=group))+geom_density(alpha=0.5)+theme_cowplot()
    ggsave("PseudotimeDensity.group.png",gg,width= 8,height=5,limitsize = FALSE,bg = 'white')
    ggsave("PseudotimeDensity.group.pdf",gg,width= 8,height=5,limitsize = FALSE)
    
    gg = ggplot(forplot,aes(x=Pseudotime))+geom_density(alpha=0.5,fill="blue")+theme_cowplot()
    ggsave("PseudotimeDensity.png",gg,width= 8,height=5,limitsize = FALSE,bg = 'white')
    ggsave("PseudotimeDensity.pdf",gg,width= 8,height=5,limitsize = FALSE)
  }
  if("cell_type" %in%colnames(forplot)){
    nrow = floor(sqrt(length(unique(forplot$cell_type))))
    ncol = ceiling(length(unique(forplot$cell_type))/nrow)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=cell_type))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~cell_type,nrow = nrow)
    ggsave("PseudotimeDensity.cell_typeEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE,bg = 'white')
    ggsave("PseudotimeDensity.cell_typeEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=cell_type))+geom_density(alpha=0.5)+theme_cowplot()
    ggsave("PseudotimeDensity.cell_type.png",gg,width= 8,height=5,limitsize = FALSE,bg = 'white')
    ggsave("PseudotimeDensity.cell_type.pdf",gg,width= 8,height=5,limitsize = FALSE)
  }
  if("inferCNV" %in%colnames(forplot)){
    nrow = floor(sqrt(length(unique(forplot$inferCNV))))
    ncol = ceiling(length(unique(forplot$inferCNV))/nrow)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=inferCNV))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~inferCNV,nrow = nrow)
    ggsave("PseudotimeDensity.inferCNVEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE,bg = 'white')
    ggsave("PseudotimeDensity.inferCNVEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=inferCNV))+geom_density(alpha=0.5)+theme_cowplot()
    ggsave("PseudotimeDensity.inferCNV.png",gg,width= 8,height=5,limitsize = FALSE,bg = 'white')
    ggsave("PseudotimeDensity.inferCNV.pdf",gg,width= 8,height=5,limitsize = FALSE)
  }
  
  if("leiden" %in%colnames(forplot)){
    nrow = floor(sqrt(length(unique(forplot$leiden))))
    ncol = ceiling(length(unique(forplot$leiden))/nrow)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=leiden))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~leiden,nrow = nrow)
    ggsave("PseudotimeDensity.leidenEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE,bg = 'white')
    ggsave("PseudotimeDensity.leidenEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=leiden))+geom_density(alpha=0.5)+theme_cowplot()
    ggsave("PseudotimeDensity.leiden.png",gg,width= 8,height=5,limitsize = FALSE,bg = 'white')
    ggsave("PseudotimeDensity.leiden.pdf",gg,width= 8,height=5,limitsize = FALSE)
  }
  
  
  if("State" %in%colnames(forplot)){
    nrow = floor(sqrt(length(unique(forplot$State))))
    ncol = ceiling(length(unique(forplot$State))/nrow)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=State))+geom_density(alpha=0.5)+theme_cowplot()+facet_wrap(~State,nrow = nrow)
    ggsave("PseudotimeDensity.StateEach.png",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE,bg = 'white')
    ggsave("PseudotimeDensity.StateEach.pdf",gg,width= 4*ncol,height=2.5*nrow,limitsize = FALSE)
    gg = ggplot(forplot,aes(x=Pseudotime,fill=State))+geom_density(alpha=0.5)+theme_cowplot()
    ggsave("PseudotimeDensity.State.png",gg,width= 8,height=5,limitsize = FALSE,bg = 'white')
    ggsave("PseudotimeDensity.State.pdf",gg,width= 8,height=5,limitsize = FALSE)
  }
  
}

################################################################################
seuset = rdsfile

# mainfun = function(seuset,geneFile,lambda,RootState,Reverse,celllist,cluster,RandomCells,
#                    cellsPercent,minCellsPerCluster){

seuset = readRDS(seuset)
seuset@active.ident = seuset$inferCNV
seuset = Psepre(seuset,celllist,cluster,RandomCells,cellsPercent,minCellsPerCluster)
print("Psepre done")
# 若没有传入基因选择列表则这里获取seuset对象里面的vargene进行使用
if(!is.null(geneFile)){
  genelist = unique(read.table(geneFile,sep="\t",header=T)[,1])
}else{
  seuset=FindVariableFeatures(seuset,selection.method ='vst',nfeatures =2000)
  genelist = VariableFeatures(seuset)
  # print(length(genelist))
  # genelist = genelist[1:100]
}
# 获取lambda参数
lambda = lambda*ncol(seuset)
# 创建monocle对象,标记用于聚类的基因，标化数据，并进行降维,排序，计算拟时间
# estimateSizeFactors 用于裱花细胞间的差异 lamda越大branch数越小
# reduceDimension max_components降维空间维数
seuCDS = ImportSeurat4CDS(seuset) %>% setOrderingFilter(genelist) %>%
  estimateSizeFactors() %>% 
  reduceDimension(reduction_method = "DDRTree",
                  lambda=lambda,
                  # auto_param_selection = F,    #细胞数量较低的时候需要有这行代码 Error in kmeans(t(Z), K, centers = centers) : initial centers are not distinct
                  # ncenter = 2
  ) %>% orderCells()
print("orderCells done")

# 修改起点
# if(!is.null(RootState) | !is.null(Reverse)){
#   seuCDS = orderCells(seuCDS,root_state=RootState,reverse=Reverse)
# }
library(reshape2)
table = table(seuCDS$inferCNV,seuCDS$State) %>% as.data.frame() %>% 
  dcast(Var1~Var2)
colnames(table)[1] = "Cluster"
table %>% write.table("Pseudotime.Statistics.txt",sep = "\t",quote = F)
seuCDS@phenoData@data %>% data.frame(CellName=row.names(seuCDS@phenoData@data)) %>%
  write.table("Pseudotime.Summary.txt",sep = "\t",quote = F,row.names = F)
PsePlot(seuCDS)
if (length(unique(seuCDS$sample))>1|length(unique(seuCDS$orig.ident))>1) {
  PsePlotSample(seuCDS)
}
if("sample"%in% colnames(seuCDS@phenoData@data)){
  forplot = data.frame(Sample=seuCDS$sample,cell_type = seuCDS$cell_type, inferCNV = seuCDS$inferCNV,group=seuCDS$group,Pseudotime=seuCDS$Pseudotime,State = seuCDS$State)
  PseDenPlot(forplot)
}
#   if("orig.ident"%in% colnames(seuCDS@phenoData@data)){
#   forplot = data.frame(Sample=seuCDS$orig.ident,Cluster=seuCDS$Cluster,Pseudotime=seuCDS$Pseudotime)
#   PseDenPlot(forplot)
#   }
#   if("group"%in% colnames(seuCDS@phenoData@data)){
#   forplot = data.frame(group=seuCDS$group,Cluster=seuCDS$Cluster,Pseudotime=seuCDS$Pseudotime)
#   PseDenPlot(forplot)
#   }
# 储存monocle rds
saveRDS(seuCDS, file = paste0("Pseudotime.monocle.rds"))
# }

# mainfun(rdsFile,geneFile = geneFile,lambda = lambda,RootState = RootState,
#         Reverse = Reverse,celllist = celllist,cluster = cluster,
#         RandomCells = RandomCells,cellsPercent = cellsPercent,
#         minCellsPerCluster = minCellsPerCluster)

################################################################################
pseudotime_df = seuCDS@phenoData@data[,c("sample","group","inferCNV","leiden","cell_type","Size_Factor","Pseudotime","State")]
df = data.frame(Cells=rownames(pseudotime_df),pseudotime_df)
write.table(df,"monocle.metadata.xls",sep = "\t",quote = F,row.names = F)


cds <- seuCDS
outpath = '/mnt/NexData_file/home/xzx/project/AE2604082792-01/shouhou2/monocle/'
for(branch_point in 1:((length(levels(cds@phenoData@data[,"State"]))-1)/2)){
  outbeam = paste0(outpath,"/BEAM_analysis.Node",branch_point)
  print(paste0("Analysis Beam Node ",branch_point))
  if(!dir.exists(outbeam)){
    dir.create(outbeam)
  }
  
  BEAM_res = BEAM(cds,branch_point = branch_point,cores = 4, progenitor_method = "duplicate")
  BEAM_res = BEAM_res[,c("gene_short_name","pval","qval")]
  BEAM_res = BEAM_res[order(BEAM_res$qval),]
  write.csv(BEAM_res,paste0(outbeam,"/BEAM_res.csv"),row.names = F)
  selectgenes = row.names(subset(BEAM_res,qval<0.05))
  p = plot_genes_branched_heatmap(cds[selectgenes[1:min(100,length(selectgenes))],],
                                  branch_point = branch_point,
                                  num_clusters = 4,
                                  cores=1,
                                  use_gene_short_name = T,
                                  show_rownames = T,return_heatmap = T)
  ggsave(paste0(outbeam,"/BEAM_heatmap.png"),p$ph_res,width=7,height = 8)
  ggsave(paste0(outbeam,"/BEAM_heatmap.pdf"),p$ph_res,width=7,height = 8)
  p=plot_genes_branched_pseudotime(cds[selectgenes[1:min(20,length(selectgenes))],],branch_point = branch_point,color_by = "State",ncol = 4)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.png"),p,width=15,height = 12)
  ggsave(paste0(outbeam,"/BEAM_pseudotime_top20.pdf"),p,width=15,height = 12)
}

fData(cds)$num_cells_expressed <- rowSums(exprs(cds) > 0)
ordergene = cds@featureData@data$gene_short_name[cds@featureData@data$use_for_ordering]
diff_test_res <- differentialGeneTest(cds[ordergene,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- diff_test_res %>% dplyr::filter(num_cells_expressed >10)%>% #防止小类干扰
  dplyr::arrange(qval) %>% head(100) %>% row.names() 

write.csv(file = './Pseudotime_ordergene.csv',diff_test_res)
p1 = plot_pseudotime_heatmap(cds[sig_gene_names,],
                             num_clusters = 5,
                             cores = 1,
                             show_rownames = T,return_heatmap = TRUE)
p1
ggsave(filename = "./monocle.top100gene_cortime.png",p1,h=8,w=6,dpi = 300)
ggsave(filename = "./monocle.top100gene_cortime.pdf",p1,h=8,w=6,dpi = 300)


################################################################################
### 绘制不同state的比例图
# 1. 提取 Monocle 对象的 State 信息
# pData(cds) 提取的是 metadata 数据框
monocle_states <- pData(cds)$State

# 2. 为提取出的 State 命名（确保能和 Seurat 匹配）
names(monocle_states) <- colnames(cds)

# 3. 将 State 添加到 Seurat 对象的 metadata 中
# Seurat 会自动根据 names(monocle_states) 匹配细胞 ID
seuset[["Monocle_State"]] <- monocle_states[colnames(seuset)]

# 4. 验证是否成功
DimPlot(seuset, group.by = "Monocle_State", label = TRUE)

# 1. 从Seurat对象提取元数据
meta_data <- seuset@meta.data  # sc 是你的Seurat对象

# 2. 统计：每个细胞类型 → 各state的数量 & 百分比
plot_df <- meta_data %>%
  group_by(cell_type, Monocle_State) %>%  # 按细胞类型+state分组
  summarise(count = n(), .groups = "drop") %>%  # 计数
  group_by(Monocle_State) %>%
  mutate(percent = count / sum(count) * 100)  # 计算占比(%)

# 3. 绘制百分比堆叠柱状图
ggplot(plot_df, 
       aes(x = Monocle_State,           # X轴：细胞类型
           y = percent,             # Y轴：百分比
           fill = factor(cell_type))) +  # 填充色：monocle状态
  geom_col(position = "fill", color = "black", size = 0.2) +  # 堆叠柱
  scale_y_continuous(labels = scales::percent) +  # Y轴显示%
  # 配色（可自定义）
  scale_fill_brewer(palette = "Paired") +
  # 标签与主题
  labs(x = "State", 
       y = "Proportion of Monocle State") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold")
  )
ggsave(plot = last_plot(),"../celltype_monoclestate_proportion.png",device = 'png',width = 9,height = 9,dpi = 300,bg = 'white')
ggsave(plot = last_plot(),"../celltype_monoclestate_proportion.pdf",device = 'pdf',width = 9,height = 9,dpi = 300,bg = 'white')


# 2. 统计：每个细胞类型 → 各state的数量 & 百分比
plot_df <- meta_data %>%
  group_by(inferCNV, Monocle_State) %>%  # 按细胞类型+state分组
  summarise(count = n(), .groups = "drop") %>%  # 计数
  group_by(Monocle_State) %>%
  mutate(percent = count / sum(count) * 100)  # 计算占比(%)

# 3. 绘制百分比堆叠柱状图
ggplot(plot_df, 
       aes(x = Monocle_State,           # X轴：细胞类型
           y = percent,             # Y轴：百分比
           fill = factor(inferCNV))) +  # 填充色：monocle状态
  geom_col(position = "fill", color = "black", size = 0.2) +  # 堆叠柱
  scale_y_continuous(labels = scales::percent) +  # Y轴显示%
  # 配色（可自定义）
  scale_fill_brewer(palette = "Paired") +
  # 标签与主题
  labs(x = "State", 
       y = "Proportion of Monocle State") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold")
  )
ggsave(plot = last_plot(),"../inferCNV_monoclestate_proportion.png",device = 'png',width = 9,height = 9,dpi = 300,bg = 'white')
ggsave(plot = last_plot(),"../inferCNV_monoclestate_proportion.pdf",device = 'pdf',width = 9,height = 9,dpi = 300,bg = 'white')


# 2. 统计：每个细胞类型 → 各state的数量 & 百分比
plot_df <- meta_data %>%
  group_by(leiden, Monocle_State) %>%  # 按细胞类型+state分组
  summarise(count = n(), .groups = "drop") %>%  # 计数
  group_by(Monocle_State) %>%
  mutate(percent = count / sum(count) * 100)  # 计算占比(%)

# 3. 绘制百分比堆叠柱状图
ggplot(plot_df, 
       aes(x = Monocle_State,           # X轴：细胞类型
           y = percent,             # Y轴：百分比
           fill = factor(leiden))) +  # 填充色：monocle状态
  geom_col(position = "fill", color = "black", size = 0.2) +  # 堆叠柱
  scale_y_continuous(labels = scales::percent) +  # Y轴显示%
  # 配色（可自定义）
  scale_fill_brewer(palette = "Paired") +
  # 标签与主题
  labs(x = "State", 
       y = "Proportion of Monocle State") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold")
  )
ggsave(plot = last_plot(),"../leiden_monoclestate_proportion.png",device = 'png',width = 9,height = 9,dpi = 300,bg = 'white')
ggsave(plot = last_plot(),"../leiden_monoclestate_proportion.pdf",device = 'pdf',width = 9,height = 9,dpi = 300,bg = 'white')


# 2. 统计：每个细胞类型 → 各state的数量 & 百分比
plot_df <- meta_data %>%
  group_by(group, Monocle_State) %>%  # 按细胞类型+state分组
  summarise(count = n(), .groups = "drop") %>%  # 计数
  group_by(Monocle_State) %>%
  mutate(percent = count / sum(count) * 100)  # 计算占比(%)

# 3. 绘制百分比堆叠柱状图
ggplot(plot_df, 
       aes(x = Monocle_State,           # X轴：细胞类型
           y = percent,             # Y轴：百分比
           fill = factor(group))) +  # 填充色：monocle状态
  geom_col(position = "fill", color = "black", size = 0.2) +  # 堆叠柱
  scale_y_continuous(labels = scales::percent) +  # Y轴显示%
  # 配色（可自定义）
  scale_fill_brewer(palette = "Paired") +
  # 标签与主题
  labs(x = "State", 
       y = "Proportion of Monocle State") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold")
  )
ggsave(plot = last_plot(),"../group_monoclestate_proportion.png",device = 'png',width = 9,height = 9,dpi = 300,bg = 'white')
ggsave(plot = last_plot(),"../group_monoclestate_proportion.pdf",device = 'pdf',width = 9,height = 9,dpi = 300,bg = 'white')



# # 2. 统计：每个细胞类型 → 各state的数量 & 百分比
# plot_df <- meta_data %>%
#   group_by(sample, Monocle_State) %>%  # 按细胞类型+state分组
#   summarise(count = n(), .groups = "drop") %>%  # 计数
#   group_by(Monocle_State) %>%
#   mutate(percent = count / sum(count) * 100)  # 计算占比(%)
# 
# # 3. 绘制百分比堆叠柱状图
# ggplot(plot_df, 
#        aes(x = Monocle_State,           # X轴：细胞类型
#            y = percent,             # Y轴：百分比
#            fill = factor(sample))) +  # 填充色：monocle状态
#   geom_col(position = "fill", color = "black", size = 0.2) +  # 堆叠柱
#   scale_y_continuous(labels = scales::percent) +  # Y轴显示%
#   # 配色（可自定义）
#   scale_fill_brewer(palette = "Paired") +
#   # 标签与主题
#   labs(x = "State", 
#        y = "Proportion of Monocle State") +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
#     axis.title = element_text(face = "bold", size = 12),
#     legend.title = element_text(face = "bold")
#   )
# ggsave(plot = last_plot(),"../sample_monoclestate_proportion.png",device = 'png',width = 9,height = 9,dpi = 300,bg = 'white')
# ggsave(plot = last_plot(),"../sample_monoclestate_proportion.pdf",device = 'pdf',width = 9,height = 9,dpi = 300,bg = 'white')
# 

save.image("../monocleDT.Rdata")

################################################################################
# 将 Pseudotime 信息添加到 Seurat 对象的 meta.data 中
seuset@meta.data$Pseudotime <- cds@phenoData@data$Pseudotime


# 提取 UMAP 和 Pseudotime 信息
umap_data <- FetchData(seuset, vars = c("UMAP_1", "UMAP_2", "Pseudotime"))

# 绘制 Pseudotime 映射到 UMAP 图
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Pseudotime)) +
  geom_point(size = 0.2) +
  scale_color_gradientn(colors = c('blue', 'cyan', 'green', 'yellow', 'orange', 'red')) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))
ggsave(plot = last_plot(),"/mnt/NexData_file/home/xzx/project/AE2604082792-01/shouhou2/monocle/UMAP_pseudotime.pdf",device = 'pdf',width = 9,height = 9,dpi = 300,bg = 'white')
ggsave(plot = last_plot(),"/mnt/NexData_file/home/xzx/project/AE2604082792-01/shouhou2/monocle/UMAP_pseudotime.png",device = 'png',width = 9,height = 9,dpi = 300,bg = 'white')


if("sample"%in% colnames(seuCDS@phenoData@data)){
  forplot = data.frame(Sample=seuCDS$sample,cell_type = seuCDS$cell_type, inferCNV = seuCDS$inferCNV,group=seuCDS$group,Pseudotime=seuCDS$Pseudotime,leiden = seuCDS$leiden)
  PseDenPlot(forplot)
}



seuset@meta.data$Monocle_State <- cds@phenoData@data$Pseudotime


# 提取 UMAP 和 Pseudotime 信息
umap_data <- FetchData(seuset, vars = c("UMAP_1", "UMAP_2", "Monocle_State"))

# 绘制 Pseudotime 映射到 UMAP 图
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = Monocle_State)) +
  geom_point(size = 0.2) +
  #scale_color_gradientn(colors = c('blue', 'cyan', 'green', 'yellow', 'orange', 'red')) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"))
ggsave(plot = last_plot(),"/mnt/NexData_file/home/xzx/project/AE2604082792-01/shouhou2/monocle/UMAP_pseudotime.pdf",device = 'pdf',width = 9,height = 9,dpi = 300,bg = 'white')
ggsave(plot = last_plot(),"/mnt/NexData_file/home/xzx/project/AE2604082792-01/shouhou2/monocle/UMAP_pseudotime.png",device = 'png',width = 9,height = 9,dpi = 300,bg = 'white')


if("sample"%in% colnames(seuCDS@phenoData@data)){
  forplot = data.frame(Sample=seuCDS$sample,cell_type = seuCDS$cell_type, inferCNV = seuCDS$inferCNV,group=seuCDS$group,Pseudotime=seuCDS$Pseudotime,leiden = seuCDS$leiden)
  PseDenPlot(forplot)
}


p = plot_genes_branched_heatmap(cds['SOD2',],
                                branch_point = branch_point,
                                num_clusters = 4,
                                cores=1,
                                use_gene_short_name = T,
                                show_rownames = T,return_heatmap = T)
p=plot_genes_branched_pseudotime(cds[selectgenes[1:min(20,length(selectgenes))],],
                                 branch_point = branch_point,color_by = "State",ncol = 4)
p