library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(readxl)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/mapping_raw/result_mapping_clean.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
# adata <- subset(adata,subset=(cell_type %in% c("B cells", "Plasma cells")))

geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC20250829/04.Mission/script/geneset.csv",sep = "\t",head=TRUE)
#转换成list
gene<-as.list(geneset)
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[18:25]<-names(gene)

outdir <- "/mnt/NexData_file/home/sh/Projects/FWSC20250829/04.Mission/AddModuleScore/"
dir.create(outdir)

for (pathway in names(gene)){
	print(pathway)
p <- VlnPlot(adata,features=pathway,pt.size = 0) + geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)
ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.pdf"),width = 27,height = 10)
ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.png"),width = 27,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0) + geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.pdf"),width = 25,height = 10)
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.png"),width = 25,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test") + geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 25,height = 10)
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.png"),width = 25,height = 10)

# mydata<-FetchData(adata,vars=c("UMAP_1","UMAP_2","geneset","group"))
# a1<-ggplot(mydata[which(mydata$group=="Negative"),],aes(x=UMAP_1,y=UMAP_2,colour=geneset)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Negative")
# a2<-ggplot(mydata[which(mydata$group=="Positive"),],aes(x=UMAP_1,y=UMAP_2,colour=geneset)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Positive")
# a3<-ggplot(mydata[which(mydata$group=="Rheumatic Diseases"),],aes(x=UMAP_1,y=UMAP_2,colour=geneset)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Rheumatic Diseases")
# a<-a1+a2+a3
# ggsave(paste0(outdir,"Main-",pathway,"-UMAP.pdf"),width = 18,height = 6)
# ggsave(paste0(outdir,"Main-",pathway,"-UMAP.png"),width = 18,height = 6)

}