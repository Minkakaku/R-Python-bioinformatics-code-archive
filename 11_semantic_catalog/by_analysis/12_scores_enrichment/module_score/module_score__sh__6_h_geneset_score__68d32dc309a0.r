library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(pheatmap)
# library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyr)
library(stringr)
# library(GSVA)
library(patchwork)
# library(limma)
# library(clusterProfiler)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(ggpubr)


genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/script/new_geneset.csv",sep = "\t",head=TRUE)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/sc_mono_anno.rds")
# adata<-subset(adata,cell_type %in% c("Dendritic cell","Classical monocyte","Non-classical monocyte"))
gene<-c(genesets$GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION)
gene <- gene[gene != ""]
gene<-list(gene)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
    name='select_geneset'
)

colnames(adata@meta.data)
colnames(adata@meta.data)[14]<-'Antigen_processing_and_presentation'
# colnames(adata@meta.data)[15]<-'GOBP_CELL_ACTIVATION'
write.table(adata@meta.data,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_h_geneset_score/monocyte-all_cells_features.xls",sep="\t")
# Idents(adata)<-adata$annotation
# p <- VlnPlot(adata,features='GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION',pt.size = 0)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_h_geneset_score/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all.VlnPlot.pdf",width = 10,height = 10)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_h_geneset_score/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(adata,features='Antigen_processing_and_presentation',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)+theme(
    legend.position = "none",                    # 去掉图例
		axis.text.x  = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank()
  )
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_h_geneset_score/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all_group.VlnPlot.pdf",width = 15,height = 10)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_h_geneset_score/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all_group.VlnPlot.png",width = 15,height = 10)

p <- get_sig_ana(p,Antigen_processing_and_presentation~split,group.by='ident',prefix="Antigen_processing_and_presentation",map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
# get_sig_ana(p$data,as.formula(Antigen_processing_and_presentation~split),group.by='ident',prefix="Antigen_processing_and_presentation",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_h_geneset_score/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_h_geneset_score/Antigen_processing_and_presentation_all_group_sig_ana.VlnPlot.pdf",width = 21,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/20.Mission/6_h_geneset_score/Antigen_processing_and_presentation_all_group_sig_ana.VlnPlot.png",width = 21,height = 10)