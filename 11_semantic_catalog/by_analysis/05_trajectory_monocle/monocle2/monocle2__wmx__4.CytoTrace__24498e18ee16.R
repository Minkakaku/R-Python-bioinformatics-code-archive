rm(list = ls())
library(Seurat)
library("CytoTRACE",lib.loc="/mnt/NexData_file/share/scRNA-tools/envs/pip_CytoTRACE/lib/R/library/")
library("HiClimR",lib.loc="/mnt/NexData_file/share/scRNA-tools/envs/pip_CytoTRACE/lib/R/library/")
library("ccaPP",lib.loc="/mnt/NexData_file/share/scRNA-tools/envs/pip_CytoTRACE/lib/R/library/")
library("nnls",lib.loc="/mnt/NexData_file/share/scRNA-tools/envs/pip_CytoTRACE/lib/R/library/")
library(stringr)
library(Matrix)
library(ggplot2)
library(cowplot)

rdsFile <-"/mnt/NexData_file/home/wmx/Project/SC2603132229-01/shouhou7/data/Neu_fixed.rds"
sce <- readRDS(rdsFile)
setwd('/mnt/NexData_file/home/wmx/Project/SC2603132229-01/shouhou7/Neu_Pseudotime/')
output <- '//mnt/NexData_file/home/wmx/Project/SC2603132229-01/shouhou7/Neu_Pseudotime/Cytotrace/'
dir.create(output)


# cytotrace analysis ------
data = GetAssayData(sce, slot = "counts")
data = as.matrix(data[,colnames(sce)])
data = na.omit(data)
results = CytoTRACE(data,ncores=1,enableFast = TRUE)
# save results
saveRDS(results,paste0(output,"/CytoTRACE.result.rds"))
CytoValue = results$CytoTRACE
Cytodata = data.frame(Cell = names(CytoValue),CytoTRACE=CytoValue)
sce = AddMetaData(sce,Cytodata)
colnames(Cytodata) = c("Cell","CytoTRACE_value")
write.table(Cytodata,file=paste0(output,"/CytoTRACE_value.txt"),sep="\t",row.names=F,quot=F)

CytoGene = results$cytoGenes
CytoGene = sort(CytoGene,decreasing = T)
Genedata = data.frame(Gene = names(CytoGene),GenesCorrelated=CytoGene)
write.table(Genedata,file=paste0(output,"/CytoTRACE_gene.txt"),sep="\t",row.names=F,quot=F)


# plot -----
# fuction reduction plot
metadataplot = function (object, features, dims = c(1, 2),pt.size = 1,cols=c("lightgrey", "blue"),reduction = "umap"){
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


if(!is.null(sce@reductions$umap)){
  gg = metadataplot(sce, features="CytoTRACE", cols = rainbow(8)[6:1], reduction = "umap")
  gg = gg + labs(color = "Predicted\norder")
  ggsave(paste0(output,"CytoTRACE.umapPlot.png"),gg,width = 9,height = 8,bg = "white")
  ggsave(paste0(output,"CytoTRACE.umapPlot.pdf"),gg,width = 9,height = 8)
}

# boxplot celltype/cluster
forplot = data.frame(sce@meta.data$CytoTRACE,sce@meta.data$cell_type)
colnames(forplot) = c("CytoTRACE","Cluster")
median = tapply(forplot$CytoTRACE,forplot$Cluster,median)
median = sort(median,decreasing = T)
order = names(median)
forplot$Cluster = factor(forplot$Cluster,levels = order)
write.table(forplot,paste0(output,"CytoTRACE.Rmd.txt"),sep ="\t",quote=F)
gg = ggplot(data = forplot, mapping = aes(x= Cluster,y = CytoTRACE, col=Cluster,fill=Cluster)) +
  geom_boxplot(alpha=0.5,outlier.size = -1) +
  theme_cowplot()+theme(legend.position = "None") +
  geom_jitter(height = 0,size=1,width = 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(output,"CytoTRACE.boxPlot.png"),gg,width = 6+0.4*length(levels(sce)),height = 8,bg = "white")
ggsave(paste0(output,"CytoTRACE.boxPlot.pdf"),gg,width = 6+0.4*length(levels(sce)),height = 8)

# gene barplot
forplot = rbind(head(Genedata,10),tail(Genedata,10))
forplot$TF = forplot$GenesCorrelated<=0
limits = ceiling(max(abs(forplot$GenesCorrelated))*100)/100
forplot$Gene = factor(forplot$Gene,levels = rev(as.character(forplot$Gene)))
gg = ggplot(data=forplot,mapping = aes(y=GenesCorrelated,x=Gene,fill=TF)) +
  geom_col()+ylab("Correlation with CytoTRACE")+ylim(-limits,limits) +
  theme_cowplot()+coord_flip()+theme(legend.position = "none")
ggsave(paste0(output,"CytoTRACE.geneBarPlot.png"),gg,width = 7,height = 0.8*10,bg = "white")
ggsave(paste0(output,"CytoTRACE.geneBarlot.pdf"),gg,width = 7,height = 0.8*10)

forplot = data.frame(sce@meta.data$CytoTRACE,sce@meta.data$group)
colnames(forplot) = c("CytoTRACE","Cluster")
median = tapply(forplot$CytoTRACE,forplot$Cluster,median)
median = sort(median,decreasing = T)
order = names(median)
forplot$Cluster = factor(forplot$Cluster,levels = order)
write.table(forplot,paste0(output,"CytoTRACE.Rmd.txt"),sep ="\t",quote=F)
gg = ggplot(data = forplot, mapping = aes(x= Cluster,y = CytoTRACE, col=Cluster,fill=Cluster)) +
  geom_boxplot(alpha=0.5,outlier.size = -1) +
  theme_cowplot()+theme(legend.position = "None") +
  geom_jitter(height = 0,size=1,width = 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# a = wilcox.test(CytoTRACE~Cluster, data = forplot ) 
# print(paste("wilcox test result ", sep = " ", "pvalue =" ,a$p.value, ";"))
# # > print(paste("wilcox test result ", sep = " ", "pvalue =" ,a$p.value, ";"))
# # [1] "wilcox test result  pvalue = 3.81636348081416e-10 ;"
# 
# ggsave(paste0(output,"CytoTRACE.boxPlot.group.png"),gg,width = 6+0.4*length(levels(seuset)),height = 8)
# ggsave(paste0(output,"CytoTRACE.boxPlot.group.pdf"),gg,width = 6+0.4*length(levels(seuset)),height = 8)
