library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(corrplot)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/all_data_sc_doublet_clean.rds")

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
colnames(adata@meta.data)[25]<-'select_geneset_Score'

Idents(adata)<-adata$annotation
p <- VlnPlot(adata,features='select_geneset_Score',pt.size = 0)+scale_fill_manual(values = c(Monocyte='#5050FFFF',Neutrophil='#CE3D32FF',CD4T='#749B58FF',CD8T="#F0E685FF","T proliferating"="#466983FF",NK="#BA6338FF",B="#5DB1DDFF"))
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/2_f_g_AddModuleScore/Main.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/2_f_g_AddModuleScore/Main.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(adata,features='select_geneset_Score',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/2_f_g_AddModuleScore/Main_group.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/2_f_g_AddModuleScore/Main_group.VlnPlot.png",width = 15,height = 10)

p <- get_sig_ana(p,select_geneset_Score~split,group.by='ident',prefix="Main_group",map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
get_sig_ana(p$data,select_geneset_Score~split,group.by='ident',prefix="Main_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/2_f_g_AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/2_f_g_AddModuleScore/Main_group_sig_ana.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/2_f_g_AddModuleScore/Main_group_sig_ana.VlnPlot.png",width = 15,height = 10)


adata_sub<-subset(adata,subset=(annotation %in% c("Monocyte")))
Idents(adata_sub)<-'group'
p <- VlnPlot(adata_sub,features='select_geneset_Score',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)+theme(
    legend.position = "none",                    # 去掉图例
		axis.text.x  = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank(),              # 可选：去掉横坐标标题
    # 加宽提琴间距：通过调整宽度实现（数值越小间距越大）
    strip.spacing.x = unit(0.8, "cm")            # 分组条间隔（可选）
  )
p <- get_sig_ana(p,select_geneset_Score~ident,group.by='ident',prefix="Monocyte_group",map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/2_f_g_AddModuleScore/Monocyte_group_sig_ana.VlnPlot.pdf",width = 7,height = 11)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/2_f_g_AddModuleScore/Monocyte_group_sig_ana.VlnPlot.png",width = 7,height = 11)

adata_sub<-subset(adata,subset=(annotation %in% c("Neutrophil")))
Idents(adata_sub)<-'group'
p <- VlnPlot(adata_sub,features='select_geneset_Score',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)+theme(
    legend.position = "none",                    # 去掉图例
		axis.text.x  = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank(),              # 可选：去掉横坐标标题
    # 加宽提琴间距：通过调整宽度实现（数值越小间距越大）
    strip.spacing.x = unit(0.8, "cm")            # 分组条间隔（可选）
  )
p <- get_sig_ana(p,select_geneset_Score~ident,group.by='ident',prefix="Neutrophil_group",map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/2_f_g_AddModuleScore/Neutrophil_group_sig_ana.VlnPlot.pdf",width = 7,height = 11)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/2_f_g_AddModuleScore/Neutrophil_group_sig_ana.VlnPlot.png",width = 7,height = 11)