library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(stringr)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/Fibroblasts.rds")
# adata <- JoinLayers(adata)
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)

geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/script/Fib.csv",sep = "\t",head=TRUE)
#转换成list
gene<-as.list(geneset)
gene <- lapply(gene, function(x) na.omit(Filter(function(y) y != "", x)))

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
    name='Fib_geneset'
)

colnames(adata@meta.data)
colnames(adata@meta.data)[15]<-'GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY'
colnames(adata@meta.data)[16]<-'GOBP_CANONICAL_WNT_SIGNALING_PATHWAY'
colnames(adata@meta.data)[17]<-'GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY'
outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/AddModuleScore/"

median_scores <- adata@meta.data %>%
  group_by(cell_type) %>%
  summarise(median_score = median(GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY, na.rm = TRUE)) %>%
  arrange(desc(median_score))
levels(adata)<-median_scores$cell_type
p <- VlnPlot(adata,features='GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY',pt.size = 0)+geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)+ labs(title = 'GOBP_TRANSFORMING_GROWTH_FACTOR \n _BETA_RECEPTOR_SIGNALING_PATHWAY')
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/AddModuleScore/GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/AddModuleScore/GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY.VlnPlot.png",width = 10,height = 10)

p <- get_sig_ana(p,GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY~ident,group.by='ident',prefix='GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY',map_sig_level=T,only_sig = T,reset_y=T)+ labs(title = 'GOBP_TRANSFORMING_GROWTH_FACTOR \n _BETA_RECEPTOR_SIGNALING_PATHWAY')
get_sig_ana(p$data,GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY~ident,group.by='ident',prefix='GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY',outdir="/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/AddModuleScore/",map_sig_level=T,only_sig = T,reset_y=T)
ggsave(paste0(outdir,"GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY_sig_ana.VlnPlot.pdf"),width = 10,height = 20)
ggsave(paste0(outdir,"GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY_sig_ana.VlnPlot.png"),width = 10,height = 20)


median_scores <- adata@meta.data %>%
  group_by(cell_type) %>%
  summarise(median_score = median(GOBP_CANONICAL_WNT_SIGNALING_PATHWAY, na.rm = TRUE)) %>%
  arrange(desc(median_score))
levels(adata)<-median_scores$cell_type
p <- VlnPlot(adata,features='GOBP_CANONICAL_WNT_SIGNALING_PATHWAY',pt.size = 0)+geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/AddModuleScore/GOBP_CANONICAL_WNT_SIGNALING_PATHWAY.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/AddModuleScore/GOBP_CANONICAL_WNT_SIGNALING_PATHWAY.VlnPlot.png",width = 10,height = 10)

p <- get_sig_ana(p,GOBP_CANONICAL_WNT_SIGNALING_PATHWAY~ident,group.by='ident',prefix='GOBP_CANONICAL_WNT_SIGNALING_PATHWAY',map_sig_level=T,only_sig = T,reset_y=T)
get_sig_ana(p$data,GOBP_CANONICAL_WNT_SIGNALING_PATHWAY~ident,group.by='ident',prefix='GOBP_CANONICAL_WNT_SIGNALING_PATHWAY',outdir="/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/AddModuleScore/",map_sig_level=T,only_sig = T,reset_y=T)
ggsave(paste0(outdir,"GOBP_CANONICAL_WNT_SIGNALING_PATHWAY_sig_ana.VlnPlot.pdf"),width = 10,height = 20)
ggsave(paste0(outdir,"GOBP_CANONICAL_WNT_SIGNALING_PATHWAY_sig_ana.VlnPlot.png"),width = 10,height = 20)


median_scores <- adata@meta.data %>%
  group_by(cell_type) %>%
  summarise(median_score = median(GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY, na.rm = TRUE)) %>%
  arrange(desc(median_score))
levels(adata)<-median_scores$cell_type
p <- VlnPlot(adata,features='GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY',pt.size = 0)+geom_boxplot(width=0.1,col="black",fill="white",outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/AddModuleScore/GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/AddModuleScore/GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY.VlnPlot.png",width = 10,height = 10)

p <- get_sig_ana(p,GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY~ident,group.by='ident',prefix='GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY',map_sig_level=T,only_sig = T,reset_y=T)
get_sig_ana(p$data,GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY~ident,group.by='ident',prefix='GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY',outdir="/mnt/NexData_file/home/sh/Projects/FWSC20250326/03.Mission/AddModuleScore/",map_sig_level=T,only_sig = T,reset_y=T)
ggsave(paste0(outdir,"GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY_sig_ana.VlnPlot.pdf"),width = 10,height = 20)
ggsave(paste0(outdir,"GOBP_REGULATION_OF_NOTCH_SIGNALING_PATHWAY_sig_ana.VlnPlot.png"),width = 10,height = 20)