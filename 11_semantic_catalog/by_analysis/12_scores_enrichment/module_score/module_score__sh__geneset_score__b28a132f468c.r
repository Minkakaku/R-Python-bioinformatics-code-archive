library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(corrplot)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Main2/all_data_anno.rds")

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
colnames(adata@meta.data)[22]<-'select_geneset_Score'

# Idents(adata)<-adata$annotation
p <- VlnPlot(adata,features='select_geneset_Score',pt.size = 0)+scale_fill_manual(values = c(Monocyte='#5050FFFF',Neutrophil='#CE3D32FF',CD4T='#749B58FF',CD8T="#F0E685FF","T proliferating"="#466983FF",NK="#BA6338FF",B="#5DB1DDFF"))
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Main.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Main.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(adata,features='select_geneset_Score',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(Before='#1f77b4',Post='#ff7f0e'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Main_group.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Main_group.VlnPlot.png",width = 15,height = 10)

p <- get_sig_ana(p,select_geneset_Score~split,group.by='ident',prefix="Main_group",map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(Before='#1f77b4',Post='#ff7f0e'))
get_sig_ana(p$data,select_geneset_Score~split,group.by='ident',prefix="Main_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Main_group_sig_ana.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Main_group_sig_ana.VlnPlot.png",width = 15,height = 10)

bdata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/Monocyte2/Monocyte_anno.rds")
bdata<-AddModuleScore(
    object=bdata,
    features=gene,
    ctr1=100, #默认值是100
    name='select_geneset'
)

colnames(bdata@meta.data)
colnames(bdata@meta.data)[22]<-'select_geneset_Score'

# Idents(adata)<-adata$annotation
p <- VlnPlot(bdata,features='select_geneset_Score',pt.size = 0)+scale_fill_manual(values = c("Non-classical monocyte"='#f8766d',"PTTG1 classical monocyte"='#b79f00',"ZFP36L1 classical monocyte"='#00ba38',"S100A8 classical monocyte"="#00bfc4",'Classical monocyte'='#619cff','Dendritic cell'='#f564e3'))+ theme(plot.margin = unit(c(1, 1, 1, 3), "lines"))  # 上、右、下、左
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Monocyte.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Monocyte.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(bdata,features='select_geneset_Score',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(Before='#1f77b4',Post='#ff7f0e'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Monocyte_group.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Monocyte_group.VlnPlot.png",width = 15,height = 10)

p <- get_sig_ana(p,select_geneset_Score~split,group.by='ident',prefix="Monocyte_group",map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(Before='#1f77b4',Post='#ff7f0e'))
get_sig_ana(p$data,select_geneset_Score~split,group.by='ident',prefix="Monocyte_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Monocyte_group_sig_ana.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Monocyte_group_sig_ana.VlnPlot.png",width = 15,height = 10)


mydata<-FetchData(adata,vars=c("UMAP_1","UMAP_2","select_geneset_Score","group"))
a1<-ggplot(mydata[which(mydata$group=="Before"),],aes(x=UMAP_1,y=UMAP_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c("#375093","#fbe2ab",'#d9950b','#f1681e')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Before")
a2<-ggplot(mydata[which(mydata$group=="Post"),],aes(x=UMAP_1,y=UMAP_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c("#375093","#fbe2ab",'#d9950b','#f1681e')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'right',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Post")
a<-a1+a2
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Main.UMAP.pdf",width = 13,height = 6)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Main.UMAP.png",width = 13,height = 6)

mydata<-FetchData(bdata,vars=c("umap_1","umap_2","select_geneset_Score","group"))
a1<-ggplot(mydata[which(mydata$group=="Before"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c("#375093","#fbe2ab",'#d9950b','#f1681e')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Before")
a2<-ggplot(mydata[which(mydata$group=="Post"),],aes(x=umap_1,y=umap_2,colour=select_geneset_Score)) + 
		geom_point(size=1) + 
		scale_color_gradientn(values=seq(0,1,0.2),colours = c("#375093","#fbe2ab",'#d9950b','#f1681e')) + 
		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'right',
    panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Post")
a<-a1+a2
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Monocyte.UMAP.pdf",width = 13,height = 6)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore/Monocyte.UMAP.png",width = 13,height = 6)

library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(pheatmap)

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
colnames(adata@meta.data)[23]<-'HALLMARK_INFLAMMATORY_RESPONSE'
colnames(adata@meta.data)[24]<-'HALLMARK_GLYCOLYSIS'
colnames(adata@meta.data)[25]<-'CYTOKINE_RECEPTOR_INTERACTION'
colnames(adata@meta.data)[26]<-'CYTOKINE'
mydata<-FetchData(adata,vars=c("cell_type","group","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_GLYCOLYSIS","CYTOKINE_RECEPTOR_INTERACTION","CYTOKINE"))
mean_scores <- mydata %>%
  group_by(cell_type, group) %>%
  summarise(across(everything(), mean))

score_select<-mean_scores$CYTOKINE
score_select <- matrix(score_select, nrow = 7, ncol = 2, byrow = TRUE)
row.names(score_select)<-c('Monocyte','Neutrophil','CD4T','CD8T',"T proliferating",'NK','B')
colnames(score_select)<-c('Before','Post')
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore_heatmap/CYTOKINE.pdf",width = 6,height = 8)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore_heatmap/CYTOKINE.png",width = 360,height = 540)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()

write.table(mean_scores,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore_heatmap/Main-heatmap.csv", sep = ",", append = TRUE, row.names = FALSE)

score_select<-as.data.frame(mean_scores[1:2,3:6])
rownames(score_select)<-c('Before','Post')
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore_heatmap/Monocyte.pdf",width = 7,height = 7)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "90",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore_heatmap/Monocyte.png",width = 540,height = 540)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "90",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()


bdata<-AddModuleScore(
    object=bdata,
    features=gene,
    ctr1=100 #默认值是100
)

colnames(bdata@meta.data)
colnames(bdata@meta.data)[23]<-'HALLMARK_INFLAMMATORY_RESPONSE'
colnames(bdata@meta.data)[24]<-'HALLMARK_GLYCOLYSIS'
colnames(bdata@meta.data)[25]<-'CYTOKINE_RECEPTOR_INTERACTION'
colnames(bdata@meta.data)[26]<-'CYTOKINE'
mydata<-FetchData(bdata,vars=c("cell_type","group","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_GLYCOLYSIS","CYTOKINE_RECEPTOR_INTERACTION","CYTOKINE"))
mean_scores <- mydata %>%
  group_by(cell_type, group) %>%
  summarise(across(everything(), mean))

score_select<-mean_scores$CYTOKINE_RECEPTOR_INTERACTION
score_select <- matrix(score_select, nrow = 6, ncol = 2, byrow = TRUE)
row.names(score_select)<-c("Non-classical monocyte","PTTG1 classical monocyte","ZFP36L1 classical monocyte","S100A8 classical monocyte",'Classical monocyte','Dendritic cell')
colnames(score_select)<-c('Before','Post')
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore_heatmap/Monocyte-CYTOKINE_RECEPTOR_INTERACTION.pdf",width = 6,height = 8)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore_heatmap/Monocyte-CYTOKINE_RECEPTOR_INTERACTION.png",width = 360,height = 540)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()

write.table(mean_scores,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/AddModuleScore_heatmap/Monocyte-heatmap.csv", sep = ",", append = TRUE, row.names = FALSE)

