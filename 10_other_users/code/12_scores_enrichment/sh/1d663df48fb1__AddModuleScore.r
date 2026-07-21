library(Seurat)
library(corrplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
# library(ggpmisc)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC20251753/01.Mission/mapping_raw/Main_mapping.rds")
levels(adata)

geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC20251753/02.Mission/script/KEGG.csv",sep = "\t",head=TRUE,strip.white = TRUE)
#转换成list
gene<-as.list(unique(geneset))
gene<-lapply(gene, function(x) x[x != ""])

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[15:32]<-names(gene)
outdir <- "/mnt/NexData_file/home/sh/Projects/FWSC20251753/02.Mission/AddModuleScore/"
dir.create(outdir)

for (pathway in names(gene)){
	print(pathway)
p <- VlnPlot(adata,features=pathway,pt.size = 0) + geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)
ggsave(paste0(outdir,pathway,"-VlnPlot.pdf"),width = 14,height = 10)
ggsave(paste0(outdir,pathway,"-VlnPlot.png"),width = 14,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0) + geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave(paste0(outdir,pathway,"-group-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,pathway,"-group-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=pathway,map_sig_level=T,only_sig = F,method = "wilcox.test")
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=pathway,outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)

mydata<-FetchData(adata,vars=c("UMAP_1","UMAP_2","group",pathway))
colnames(mydata)[4]<-'pathway'
a1<-ggplot(mydata[which(mydata$group=="L"),],aes(x=UMAP_1,y=UMAP_2,colour=pathway)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c('#d3d3d3','#195e6c','#eb2d4a','#f1681e','#ffeee5')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'right',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="L")
a2<-ggplot(mydata[which(mydata$group=="NL"),],aes(x=UMAP_1,y=UMAP_2,colour=pathway)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c('#d3d3d3','#195e6c','#eb2d4a','#f1681e','#ffeee5')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'right',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="NL")
a<-a1+a2
ggsave(paste0(outdir,pathway,"-UMAP.pdf"),width = 13,height = 6)
ggsave(paste0(outdir,pathway,"-UMAP.png"),width = 13,height = 6)

}