library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/all_data_sc_doublet.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)

geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/05.Mission/Enrich/enrich_select.csv",sep = ",",head=TRUE)
#转换成list
gene<-as.list(geneset)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
    name='select_geneset'
)

colnames(adata@meta.data)
colnames(adata@meta.data)[23]<-'select_geneset_Score'
p <- VlnPlot(adata,features='select_geneset_Score',pt.size = 0)+scale_fill_manual(values = c(Monocyte='#5050FFFF',Neutrophil='#CE3D32FF',CD4T='#749B58FF',CD8T="#F0E685FF","T proliferating"="#466983FF",NK="#BA6338FF",B="#5DB1DDFF"))
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/AddModuleScore/Main.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/AddModuleScore/Main.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(adata,features='select_geneset_Score',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/AddModuleScore/Main_group.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/AddModuleScore/Main_group.VlnPlot.png",width = 15,height = 10)

p <- get_sig_ana(p,select_geneset_Score~split,group.by='ident',prefix="Main_group",map_sig_level=T,only_sig = F,method = "t.test")+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
get_sig_ana(p$data,select_geneset_Score~split,group.by='ident',prefix="Main_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/AddModuleScore/",map_sig_level=T,only_sig = F,method = "t.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/AddModuleScore/Main_group_sig_ana.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/AddModuleScore/Main_group_sig_ana.VlnPlot.png",width = 15,height = 10)
# data<-FetchData(adata,vars =c("annotation","select_geneset_Score"))
# p<-ggplot(data,aes(cell_type,select_geneset_Score))
# p+geom_boxplot()+theme_bw()+RotatedAxis()

# mydata<-FetchData(adata,vars=c("umap_1","umap_2","select_geneset_Score"))
# a<-ggplot(mydata,aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
mydata<-FetchData(adata,vars=c("umap_1","umap_2","select_geneset_Score","group"))
a1<-ggplot(mydata[which(mydata$group=="LAHS"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="LAHS")
a2<-ggplot(mydata[which(mydata$group=="NLAHS"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="NLAHS")
a3<-ggplot(mydata[which(mydata$group=="NC"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="NC")
a<-a1+a2+a3
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/AddModuleScore/Main.UMAP.pdf",width = 18,height = 6)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/AddModuleScore/Main.UMAP.png",width = 18,height = 6)


groupfile<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/02.Mission/script/groupfile.txt")
colnames(groupfile)<-c("project","sample","group")
outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/AddModuleScore/"
cluster <- c("mono")
adata<-readRDS(paste0("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/02.Mission/data/sc_mono.rds"))
adata <- JoinLayers(adata)
adata$project<- ifelse(adata$project == "YLY", "YLY0701", adata$project)
adata@meta.data <- adata@meta.data  %>% rownames_to_column(var = "cell_id")%>% left_join(groupfile,by="project") %>% column_to_rownames(var = "cell_id")
Samples<-c("HLH01-LAHS-B","HLH01-LAHS-B-Post","HLH03-LAHS-T","HLH03-LAHS-T-Post","HLH07-NLAHS","HLH07-NLAHS-Post","HC01","HC02")
barcodes<-rownames(adata@meta.data[which(adata$sample %in% Samples),])
adata<-subset(adata,cells = barcodes)
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
    name='select_geneset'
)

colnames(adata@meta.data)
colnames(adata@meta.data)[16]<-'select_geneset_Score'
p <- VlnPlot(adata,features='select_geneset_Score',pt.size = 0)+scale_fill_manual(values = c("CD14 Monocyte"='#F8766D',"FCGR3A Monocyte"='#00BFC4',"Dendritic cell"='#749B58FF',"TPX2 Monocyte"="#C77CFF"))
ggsave(paste0(outdir,cluster,".VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,cluster,".VlnPlot.png"),width = 10,height = 10)

p <- get_sig_ana(p,select_geneset_Score~ident,group.by='ident',prefix=cluster,map_sig_level=T,only_sig = F,reset_y=T,method = "t.test")+scale_fill_manual(values = c("CD14 Monocyte"='#F8766D',"FCGR3A Monocyte"='#00BFC4',"Dendritic cell"='#749B58FF',"TPX2 Monocyte"="#C77CFF"))
get_sig_ana(p$data,select_geneset_Score~ident,group.by='ident',prefix=cluster,outdir=outdir,map_sig_level=T,only_sig = F,reset_y=T,method = "t.test")
ggsave(paste0(outdir,cluster,"_sig_ana.VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,cluster,"_sig_ana.VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(adata,features='select_geneset_Score',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
ggsave(paste0(outdir,cluster,"_group.VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,cluster,"_group.VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,select_geneset_Score~split,group.by='ident',prefix=paste0(cluster,"_group"),map_sig_level=T,only_sig = F,method = "t.test")+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
get_sig_ana(p$data,select_geneset_Score~split,group.by='ident',prefix=paste0(cluster,"_group"),outdir=outdir,map_sig_level=T,only_sig = F,method = "t.test")
ggsave(paste0(outdir,cluster,"_group_sig_ana.VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,cluster,"_group_sig_ana.VlnPlot.png"),width = 15,height = 10)
# data<-FetchData(adata,vars =c("annotation","select_geneset_Score"))
# p<-ggplot(data,aes(cell_type,select_geneset_Score))
# p+geom_boxplot()+theme_bw()+RotatedAxis()

mydata<-FetchData(adata,vars=c("umap_1","umap_2","select_geneset_Score","group"))
a1<-ggplot(mydata[which(mydata$group=="LAHS"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="LAHS")
a2<-ggplot(mydata[which(mydata$group=="NLAHS"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="NLAHS")
a3<-ggplot(mydata[which(mydata$group=="NC"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="NC")
a<-a1+a2+a3
ggsave(paste0(outdir,cluster,".UMAP.pdf"),width = 18,height = 6)
ggsave(paste0(outdir,cluster,".UMAP.png"),width = 18,height = 6)
