library(Seurat)
library(corrplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(ggpmisc)

adata <- readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/06.Mission/Myeloid cells/Myeloid.rds")
MIF<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
#转换成list
gene<-as.list(MIF)
names(gene)<-c('MIF')
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[18]<-'MIF'
outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/06.Mission/AddModuleScore/"
levels(adata$cell_type)

for (pathway in c('MIF')){
	print(pathway)
p <- VlnPlot(adata,features=pathway,pt.size = 0)+scale_fill_manual(values = c('DUSP1 Classical monocyte'="#1F77B4",'ISG15 Classical monocyte'="#ff7f0e",'Non-classical monocyte'="#d62728",'VCAN Classical monocyte'="#aa40fc",'cDCs'="#8c564b"))
ggsave(paste0(outdir,"Myeloid-",pathway,"-VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,"Myeloid-",pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c('NC'='#1f77b4','Negative'='#ff7f0e','Positive'='#2ca02c','Rheumatic Diseases'='#d62728'))
ggsave(paste0(outdir,"Myeloid_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Myeloid_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Myeloid_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c('NC'='#1f77b4','Negative'='#ff7f0e','Positive'='#2ca02c','Rheumatic Diseases'='#d62728'))
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Myeloid_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Myeloid_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 20,height = 10)
ggsave(paste0(outdir,"Myeloid_group_sig_ana-",pathway,"-VlnPlot.png"),width = 20,height = 10)

mydata<-FetchData(adata,vars=c("UMAP_1","UMAP_2",'group',pathway))
a4<-ggplot(mydata[which(mydata$group=="NC"),],aes(x=UMAP_1,y=UMAP_2,colour=MIF)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="NC")
a1<-ggplot(mydata[which(mydata$group=="Negative"),],aes(x=UMAP_1,y=UMAP_2,colour=MIF)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Negative")
a2<-ggplot(mydata[which(mydata$group=="Positive"),],aes(x=UMAP_1,y=UMAP_2,colour=MIF)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Positive")
a3<-ggplot(mydata[which(mydata$group=="Rheumatic Diseases"),],aes(x=UMAP_1,y=UMAP_2,colour=MIF)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Rheumatic Diseases")
a<-a4+a1+a2+a3
ggsave(paste0(outdir,"Myeloid-",pathway,"-UMAP.pdf"),width = 13,height = 12)
ggsave(paste0(outdir,"Myeloid-",pathway,"-UMAP.png"),width = 13,height = 12)

}