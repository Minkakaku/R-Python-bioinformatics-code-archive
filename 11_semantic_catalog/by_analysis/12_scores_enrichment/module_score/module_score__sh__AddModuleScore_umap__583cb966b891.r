library(Seurat)
library(corrplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(ggpmisc)

MIF<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
IL6<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/il6_gene.csv",sep = ",",head=FALSE)
gene<-list(c(MIF$V1),c(IL6$V1))

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/script/Main_mapping_clean.rds")
outdir<-"/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/fig3_C_score_umap/"
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[17]<-'MIF'
colnames(adata@meta.data)[18]<-'IL6'

for (pathway in c('MIF','IL6')){
	print(pathway)

mydata<-FetchData(adata,vars=c("UMAP_1","UMAP_2",pathway,"group"))
a1<-ggplot(mydata[which(mydata$group=="N"),],aes(x=UMAP_1,y=UMAP_2,colour=MIF)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c("royalblue","royalblue","white","firebrick3","firebrick3")) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="N")
a2<-ggplot(mydata[which(mydata$group=="P"),],aes(x=UMAP_1,y=UMAP_2,colour=MIF)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c("royalblue","royalblue","white","firebrick3","firebrick3")) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="P")
a3<-ggplot(mydata[which(mydata$group=="AID"),],aes(x=UMAP_1,y=UMAP_2,colour=MIF)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c("royalblue","royalblue","white","firebrick3","firebrick3")) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="AID")
a<-a1+a2+a3
ggsave(paste0(outdir,"Main-",pathway,"-UMAP.pdf"),width = 18,height = 6)
ggsave(paste0(outdir,"Main-",pathway,"-UMAP.png"),width = 18,height = 6)
ggsave(paste0(outdir,"Main-",pathway,"-UMAP.tiff"),width = 18,height = 6)
}