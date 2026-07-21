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
Idents(adata)<-adata$raw_cell_type

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

outdir <- "/mnt/NexData_file/home/sh/Projects/FWSC20250829/05.Mission/AddModuleScore/"
dir.create(outdir)

for (pathway in names(gene)){
	print(pathway)
p <- VlnPlot(adata,features=pathway,pt.size = 0) + geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)+scale_fill_manual(values = c('B cells'='#1f77b4','Mast cells'='#ff7f0e','Monocytes'='#279e68','NK cells'='#d62728','Plasmacytoid dendritic cells'='#8c564b','Platelets'='#e377c2','T cells'='#b5bd61'))
ggsave(paste0(outdir,pathway,"-VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0) + geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)+scale_fill_manual(values = c('case'='#1f77b4','control'='#ff7f0e'))
ggsave(paste0(outdir,"group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test") + geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)

mydata<-FetchData(adata,vars=c("UMAP_1","UMAP_2",pathway,"group"))
colnames(mydata)[3]<-'geneset'
a1<-ggplot(mydata[which(mydata$group=="case"),],aes(x=UMAP_1,y=UMAP_2,colour=geneset)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c('#d3d3d3','#195e6c','#eb2d4a','#f1681e','#ffeee5')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="case")
a2<-ggplot(mydata[which(mydata$group=="control"),],aes(x=UMAP_1,y=UMAP_2,colour=geneset)) + 
		geom_point(size=1) + 
		scale_color_gradientn(colours = c('#d3d3d3','#195e6c','#eb2d4a','#f1681e','#ffeee5')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'right',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="control")
a<-a1+a2
ggsave(paste0(outdir,pathway,"-UMAP.pdf"),width = 13,height = 6)
ggsave(paste0(outdir,pathway,"-UMAP.png"),width = 13,height = 6)

mydata<-FetchData(adata,vars=c("raw_cell_type",pathway,"group"))
mean_scores <- mydata %>%
  group_by(raw_cell_type, group) %>%
  summarise(across(everything(), mean))

score_select<-as.data.frame(mean_scores[,3])
score_select <- matrix(score_select[,1], nrow = 7, ncol = 2, byrow = TRUE)
row.names(score_select)<-c('B cells','Mast cells','Monocytes','NK cells','Plasmacytoid dendritic cells','Platelets','T cells')
colnames(score_select)<-c("case","control")
pdf(paste0(outdir,pathway,"-heatmap.pdf"),width = 6,height = 8)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()
png(paste0(outdir,pathway,"-heatmap.png"),width = 360,height = 540)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()

}