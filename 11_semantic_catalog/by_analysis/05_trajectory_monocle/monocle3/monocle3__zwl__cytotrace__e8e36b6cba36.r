rm(list = ls())
#.libPaths("/home/wzb/R/x86_64-pc-linux-gnu-library/4.2")
library(Seurat)
library(getopt)
library(stringr)
library(Matrix)
# BiocManager::install("sva")
# devtools::install_local("/home/wzb/pkg/CytoTRACE_0.3.3.tar.gz")

library(CytoTRACE)
library(ggplot2)
library(cowplot)
command = matrix(c('rds','r',1,"character",
                   'outpath','o',1,"character",
                   'plotgeneNum','n',1,"numeric",
                   'colorTab','C',2,"character",
                   'celllist','c',2,"character")
                 ,byrow=TRUE, ncol=4)
args=getopt(command)
options(stringsAsFactors = FALSE)

rdsFile = args$rds;
outpath = args$outpath;
plotgeneNum = args$plotgeneNum
if(is.null(args$celllist)){celllist=NULL
}else {
  celllist = args$celllist
}
if(is.null(args$colorTab)){colorTab=NULL
}else {
  colorTab = args$colorTab
}

rdsFile = 'F:/FW2024-285/售后2/pseudotime/Stem_epc.rds'
outpath ='G:/FW2024-348/FW2024-348_10/售后4/monocle3/result/CytoTRACE'
plotgeneNum = 20



changeColor = FALSE
#seuset = readRDS(rdsFile)
seuset=adata
Idents(seuset)='cell_type'
seuset= subset(seuset, idents = c("Hepatic stellate cells", "Fibroblasts","Hepatocytes","Cholangiocyte"))
seuset= subset(seuset, idents = c("Hepatic stellate cells", "Fibroblasts","Hepatocytes","Cholangiocyte"))



Idents(seuset)='group'
seuset= subset(seuset, idents = c("KO-ccl4-010","KO-ccl4-488","WT-ccl4-4"))
write.table(seuset@meta.data,file = '../meta.txt', sep = "\t", quote = FALSE, row.names = TRUE)
new_metadata_df <- read.table("../meta1.txt", header = TRUE, sep = "\t", row.names = 1)
seuset@meta.data <- new_metadata_df


#Idents(seuset) = seuset@meta.data$cell_type
new.cluster.ids=c('KO','KO','WT-ccl4-4') 
names(new.cluster.ids) <- levels(seuset)
seuset <- RenameIdents(seuset, new.cluster.ids)
setwd(outpath)
# 绘制umap图
#seuset=scdata
# seuset@assays$RNA@key 
metadataplot = function (object, features, dims = c(1, 2),pt.size = 1,cols=c("lightgrey", "blue"),reduction = "UMAP"){
  dims = paste0(Key(object = object[[reduction]]), dims)
  data = FetchData(object = object, vars = c(dims,features))
  gg = ggplot(data=data,
              mapping = aes_string(x = dims[1],y=dims[2],color=features))+
    geom_point(size=pt.size)+
    scale_color_gradientn(colours = cols)+
    labs(color = NULL,title = features)+ 
    theme_cowplot()+
    theme(plot.title = element_text(hjust = 0.5,size = 30),
          axis.title.x=element_text(size=25),
          axis.title.y=element_text(size=25),
          axis.text.x=element_text(size=15),
          axis.text.y=element_text(size=15))
  return(gg)}

mainfun = function(seuset,plotgeneNum,celllist,colorTab){
  # 抽取细胞重新定义类型绘制图片
  if(!is.null(celllist)){
    celllist = read.table(celllist,sep="\t",header=F,comment.char = "")
    realcluster = celllist[is.finite(match(celllist[,1],colnames(seuset))),]
    seuset = subset(seuset,cells=as.character(realcluster[,1]))
    seuset = SetIdent(object = seuset, cells= as.character(realcluster[,1]), as.character(realcluster[,2]))
  }
  
  # 对细胞类型定义定义颜色
  if(!is.null(colorTab)){
    colorlist = read.table(colorTab,sep="\t",header=F,comment.char = "")
    realColor = colorlist[is.finite(match(colorlist[,1],levels(seuset))),]
    if(nrow(realColor)!= length(levels(seuset@active.ident))){
      changeColor = FALSE
      print("ColorList is Error, uses default colors.")}else{
        type = as.character(realColor[,1])
        color = as.character(realColor[,2])
        names(color) = type
        changeColor = TRUE}}
  
  if("SCT" %in% names(seuset@assays)){
    DefaultAssay(seuset) = "SCT"}else{
      DefaultAssay(seuset) = "RNA"}
  
  # 抽取原始数据
  data = GetAssayData(seuset, slot = "counts")
  data = as.matrix(data[,colnames(seuset)])
  results = CytoTRACE(data,ncores=1)
  saveRDS(results,"CytoTRACE.result.rds")
  #results=readRDS('./CytoTRACE/CytoTRACE.result.rds')
  
  CytoValue = results$CytoTRACE
  Cytodata = data.frame(Cell = names(CytoValue),CytoTRACE=CytoValue)
  seuset = AddMetaData(seuset,Cytodata)
  colnames(Cytodata) = c("Cell","CytoTRACE_value")
  write.table(Cytodata,file="CytoTRACE.value.txt",sep="\t",row.names=F,quot=F)
  # tmp = read.csv("CytoTRACE.value.txt",sep = "\t")
  # row.names(tmp) = tmp[,1]
  # tmp = tmp[,-1,drop=F]
  # seuset$CytoTRACE = tmp
  CytoGene = results$cytoGenes
  CytoGene = sort(CytoGene,decreasing = T)
  Genedata = data.frame(Gene = names(CytoGene),GenesCorrelated=CytoGene)
  write.table(Genedata,file="CytoTRACE.gene.txt",sep="\t",row.names=F,quot=F)
  if(!is.null(seuset@reductions$umap)){
    gg = metadataplot(seuset, features="CytoTRACE", cols = rainbow(8)[6:1], reduction = "umap")
    gg = gg + labs(color = "Predicted\norder")
    ggsave("CytoTRACE.umapPlot.png",gg,width = 9,height = 8)
    ggsave("CytoTRACE.umapPlot.pdf",gg,width = 9,height = 8)
  }
  
  if(!is.null(seuset@reductions$tsne)){
    gg = metadataplot(seuset, features="CytoTRACE", cols = rainbow(8)[6:1], reduction = "tsne")
    gg = gg + labs(color = "Predicted\norder")
    ggsave("CytoTRACE.tsnePlot.png",gg,width = 9,height = 8)
    ggsave("CytoTRACE.tsnePlot.pdf",gg,width = 9,height = 8)
  }
  
  forplot = data.frame(seuset@meta.data$CytoTRACE,seuset@active.ident)
  colnames(forplot) = c("CytoTRACE","Cluster")
  median = tapply(forplot$CytoTRACE,forplot$Cluster,median)
  median = sort(median,decreasing = T)
  order = names(median)
  forplot$Cluster = factor(forplot$Cluster,levels = order)
  write.table(forplot,"CytoTRACE.Rmd.txt",sep ="\t",quote=F)
  gg = ggplot(data = forplot, mapping = aes(x= Cluster,y = CytoTRACE, col=Cluster,fill=Cluster))+geom_boxplot(alpha=0.5,outlier.size = -1)+theme_cowplot()+theme(legend.position = "None")+geom_jitter(height = 0,size=1,width = 0.1)+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if(changeColor){
    plotcolor = color[order]
    gg = gg+scale_color_manual(values = plotcolor)+scale_fill_manual(values = plotcolor)
  }
  ggsave("CytoTRACE.boxPlot.png",gg,width = 6+0.4*length(levels(seuset)),height = 8)
  ggsave("CytoTRACE.boxPlot.pdf",gg,width = 6+0.4*length(levels(seuset)),height = 8)
  
  forplot = rbind(head(Genedata,plotgeneNum),tail(Genedata,plotgeneNum))
  forplot$TF = forplot$GenesCorrelated<=0
  limits = ceiling(max(abs(forplot$GenesCorrelated))*100)/100
  forplot$Gene = factor(forplot$Gene,levels = rev(as.character(forplot$Gene)))
  gg = ggplot(data=forplot,mapping = aes(y=GenesCorrelated,x=Gene,fill=TF))+geom_col()+ylab("Correlation with CytoTRACE")+ylim(-limits,limits)+theme_cowplot()+coord_flip()+theme(legend.position = "none")
  ggsave("CytoTRACE.geneBarPlot.png",gg,width = 7,height = 0.8*plotgeneNum)
  ggsave("CytoTRACE.geneBarlot.pdf",gg,width = 7,height = 0.8*plotgeneNum)
  
}

mainfun(seuset,plotgeneNum,celllist,colorTab)
