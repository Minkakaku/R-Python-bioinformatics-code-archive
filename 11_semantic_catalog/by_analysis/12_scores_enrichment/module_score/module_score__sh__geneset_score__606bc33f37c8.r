library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(corrplot)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/all_data_anno.rds")
# adata$group<-factor(adata$group)
# adata$sample<-factor(adata$sample)

# geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/NLAHS_LAHS-NC/Venn_result_select.csv",sep = ",",head=TRUE)
geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/NLAHS_LAHS-NC/Venn_result_select.csv",sep = ",",head=FALSE)
#转换成list
gene<-as.list(geneset)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
    name='select_geneset'
)

colnames(adata@meta.data)
colnames(adata@meta.data)[21]<-'select_geneset_Score'

# Idents(adata)<-adata$annotation
p <- VlnPlot(adata,features='select_geneset_Score',pt.size = 0)+scale_fill_manual(values = c(Monocyte='#5050FFFF',Neutrophil='#CE3D32FF',CD4T='#749B58FF',CD8T="#F0E685FF","T proliferating"="#466983FF",NK="#BA6338FF",B="#5DB1DDFF"))
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Main.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Main.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(adata,features='select_geneset_Score',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(Before='#1f77b4',Post='#ff7f0e'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Main_group.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Main_group.VlnPlot.png",width = 15,height = 10)

p <- get_sig_ana(p,select_geneset_Score~split,group.by='ident',prefix="Main_group",map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(Before='#1f77b4',Post='#ff7f0e'))
get_sig_ana(p$data,select_geneset_Score~split,group.by='ident',prefix="Main_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Main_group_sig_ana.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Main_group_sig_ana.VlnPlot.png",width = 15,height = 10)

bdata<-subset(adata,idents='Monocyte')
Idents(bdata)<-bdata$group
p <- VlnPlot(bdata,features='select_geneset_Score',pt.size = 0)+scale_fill_manual(values = c(Before='#1f77b4',Post='#ff7f0e'))+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Monocyte_group.VlnPlot.pdf",width = 7,height = 7)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Monocyte_group.VlnPlot.png",width = 7,height = 7)

p <- get_sig_ana(p,select_geneset_Score~ident,group.by='ident',prefix="Monocyte_group",map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)+scale_fill_manual(values = c(Before='#1f77b4',Post='#ff7f0e'))
get_sig_ana(p$data,select_geneset_Score~ident,group.by='ident',prefix="Monocyte_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Monocyte_group_sig_ana.VlnPlot.pdf",width = 7,height = 9)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Monocyte_group_sig_ana.VlnPlot.png",width = 7,height = 9)
# data<-FetchData(adata,vars =c("annotation","select_geneset_Score"))
# p<-ggplot(data,aes(cell_type,select_geneset_Score))
# p+geom_boxplot()+theme_bw()+RotatedAxis()

# mydata<-FetchData(adata,vars=c("umap_1","umap_2","select_geneset_Score"))
# a<-ggplot(mydata,aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(values=seq(0,1,0.2),colours = c("#375093","#fbe2ab",'#d9950b','#f1681e')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))
mydata<-FetchData(adata,vars=c("umap_1","umap_2","select_geneset_Score","group"))
a1<-ggplot(mydata[which(mydata$group=="Before"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c("#375093","#fbe2ab",'#d9950b','#f1681e')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Before")
a2<-ggplot(mydata[which(mydata$group=="Post"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c("#375093","#fbe2ab",'#d9950b','#f1681e')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Post")
a<-a1+a2
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Main.UMAP.pdf",width = 12,height = 6)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore/Main.UMAP.png",width = 12,height = 6)


library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(pheatmap)
adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main/all_data_anno.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)

geneset1<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/pathway.csv",sep = "\t",head=TRUE)
geneset2<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/cytokine_genes.csv",sep = "\t",head=FALSE)

#转换成list
gene<-c(as.list(geneset1),as.list(geneset2))

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100 #默认值是100
)

colnames(adata@meta.data)
colnames(adata@meta.data)[21]<-'HALLMARK_INFLAMMATORY_RESPONSE'
colnames(adata@meta.data)[22]<-'HALLMARK_GLYCOLYSIS'
colnames(adata@meta.data)[23]<-'CYTOKINE_RECEPTOR_INTERACTION'
colnames(adata@meta.data)[24]<-'CYTOKINE'
mydata<-FetchData(adata,vars=c("cell_type","group","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_GLYCOLYSIS","CYTOKINE_RECEPTOR_INTERACTION","CYTOKINE"))
mean_scores <- mydata %>%
  group_by(cell_type, group) %>%
  summarise(across(everything(), mean))

score_select<-mean_scores$CYTOKINE_RECEPTOR_INTERACTION
score_select <- matrix(score_select, nrow = 7, ncol = 2, byrow = TRUE)
row.names(score_select)<-c('Monocyte','Neutrophil','CD4T','CD8T',"T proliferating",'NK','B')
colnames(score_select)<-c('Before','Post')
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore_heatmap/CYTOKINE_RECEPTOR_INTERACTION.pdf",width = 6,height = 8)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "row",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore_heatmap/CYTOKINE_RECEPTOR_INTERACTION.png",width = 360,height = 540)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "row",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()

write.table(mean_scores,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/AddModuleScore_heatmap/CYTOKINE_RECEPTOR_INTERACTION.csv", sep = ",", append = TRUE, row.names = FALSE)

