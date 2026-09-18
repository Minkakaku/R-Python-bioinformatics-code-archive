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
library(clusterProfiler)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")
library(ggpubr)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/TandNK/sc_TandNK_anno.rds")
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
# Idents(adata) <- "annotation"
expr <- AverageExpression(adata, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] #选取非零基因
expr <- as.matrix(expr)
head(expr)
outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/GSVA/"
dir.create(outdir)

genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/script/new_geneset.csv",sep = "\t",head=TRUE)

gene<-c(genesets$T.cell.exhaustion)
gene <- gene[gene != ""]
gene<-list(gene)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
    name='select_geneset'
)

colnames(adata@meta.data)
colnames(adata@meta.data)[23]<-'GOBP_MITOCHONDRIAL_GENE_EXPRESSION'
colnames(adata@meta.data)[23]<-'T.cell.exhaustion'
# Idents(adata)<-adata$annotation
p <- VlnPlot(adata,features='T.cell.exhaustion',pt.size = 0)+scale_fill_manual(values = c('CD4_Tcm'="#A6CEE3","CD4_Tem"="#3F8EAA","CD4_TNaive"="#79C360","CD4_Treg"="#B89B74","CD8_Tem"="#E52829","CD8_Trm"="#FDB762","CD8_Tex"="#ED8F47","CD8_MAIT"="#9471B4","T proliferating"="#DDD399","NK"="#B15928"))
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/T.cell.exhaustion.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/T.cell.exhaustion.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(adata,features='T.cell.exhaustion',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/T.cell.exhaustion_group.VlnPlot.pdf",width = 18,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/T.cell.exhaustion_group.VlnPlot.png",width = 18,height = 10)

p <- get_sig_ana(p,T.cell.exhaustion~split,group.by='ident',prefix="T.cell.exhaustion_group",map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
get_sig_ana(p$data,T.cell.exhaustion~split,group.by='ident',prefix="T.cell.exhaustion_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/T.cell.exhaustion_group_sig_ana.VlnPlot.pdf",width = 18,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/T.cell.exhaustion_group_sig_ana.VlnPlot.png",width = 18,height = 10)

df<-p$data
df$cell_type = factor(df$ident,levels = c('CD4_Tcm',"CD4_Tem","CD4_TNaive","CD4_Treg","CD8_Tem","CD8_Trm","CD8_Tex","CD8_MAIT","T proliferating","NK"))
gg = ggplot(data=df,mapping=aes(x=cell_type,y=df$`T.cell.exhaustion`))+
  geom_jitter(aes(colour  = cell_type),
              width = 0.1,  size=2.5, stroke=0.5, shape=16,alpha=10) +
  geom_violin(alpha=0.65, linewidth=1, width=0.5,adjust = 0.9) +
  geom_boxplot(width=0.3, linewidth=1,outlier.color = NA, fill=NA) +
  scale_colour_manual(values=c("#A6CEE3","#3F8EAA","#79C360","#B89B74", "#E52829", "#FDB762", "#ED8F47", "#9471B4","#DDD399","#B15928"))+
  stat_compare_means(method = "wilcox.test",ref.group=".all.",  label = "p.format") +
  # stat_compare_means(method = "wilcox.test", comparisons = combn(celltype, 2, simplify = FALSE), label = "p.format") +
  theme_bw()+theme(panel.grid=element_blank(),axis.text.x=element_text(color='black', angle=45, hjust=1))+
  geom_hline(yintercept  = mean(df$`T.cell.exhaustion`),
             linetype = "dashed", color = "red4", size = 1)+
  labs(x="",y = "T.cell.exhaustion")
# +
#    facet_wrap(~ celltype, nrow = 2)  # 按照 cell_type 分面
gg
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/T.cell.exhaustion-celltype.pdf",gg,width = 12,height = 5)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/T.cell.exhaustion-celltype.png",gg,width = 12,height = 5)



adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/sc_mono_anno.rds")
# adata<-subset(adata,cell_type %in% c("Dendritic cell","Classical monocyte","Non-classical monocyte"))

gene<-c(genesets$GOBP_CELL_ACTIVATION)
gene <- gene[gene != ""]
gene<-list(gene)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
    name='select_geneset'
)

colnames(adata@meta.data)
colnames(adata@meta.data)[14]<-'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION'
colnames(adata@meta.data)[15]<-'GOBP_CELL_ACTIVATION'
write.table(adata@meta.data,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/monocyte-all_cells_features.xls",sep="\t")
# Idents(adata)<-adata$annotation
p <- VlnPlot(adata,features='GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION',pt.size = 0)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(adata,features='GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all_group.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all_group.VlnPlot.png",width = 15,height = 10)

p <- get_sig_ana(p,GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION~split,group.by='ident',prefix="GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_group",map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
get_sig_ana(p$data,GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION~split,group.by='ident',prefix="GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all_group_sig_ana.VlnPlot.pdf",width = 15,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_all_group_sig_ana.VlnPlot.png",width = 15,height = 10)

# p <- VlnPlot(adata,features='GOBP_CELL_ACTIVATION',pt.size = 0)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_CELL_ACTIVATION.VlnPlot.pdf",width = 7,height = 10)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_CELL_ACTIVATION.VlnPlot.png",width = 7,height = 10)

# p <- VlnPlot(adata,features='GOBP_CELL_ACTIVATION',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_CELL_ACTIVATION_group.VlnPlot.pdf",width = 10,height = 10)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_CELL_ACTIVATION_group.VlnPlot.png",width = 10,height = 10)

# p <- get_sig_ana(p,GOBP_CELL_ACTIVATION~split,group.by='ident',prefix="GOBP_CELL_ACTIVATION_group",map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
# get_sig_ana(p$data,GOBP_CELL_ACTIVATION~split,group.by='ident',prefix="GOBP_CELL_ACTIVATION_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_CELL_ACTIVATION_group_sig_ana.VlnPlot.pdf",width = 10,height = 10)
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0922_AddModuleScore/GOBP_CELL_ACTIVATION_group_sig_ana.VlnPlot.png",width = 10,height = 10)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/TandNK/sc_TandNK_anno.rds")
genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/script/new_geneset.csv",sep = "\t",head=TRUE)

gene<-c(genesets$T.cell.exhaustion)
gene <- gene[gene != ""]
gene<-list(gene)

Idents(adata)<-adata$annotation
adata_sub<-subset(adata,idents=c("NK"),invert=TRUE)
adata_sub<-AddModuleScore(
    object=adata_sub,
    features=gene,
    ctr1=100, #默认值是100
    name='select_geneset'
)

colnames(adata_sub@meta.data)
colnames(adata_sub@meta.data)[23]<-'T.cell.exhaustion'

p <- VlnPlot(adata_sub,features='T.cell.exhaustion',pt.size = 0)+scale_fill_manual(values = c(CD4T='#749B58FF',CD8T="#F0E685FF","T proliferating"="#466983FF"))
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/T.cell.exhaustion.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/T.cell.exhaustion.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(adata_sub,features='T.cell.exhaustion',split.by = 'group',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/T.cell.exhaustion_group.VlnPlot.pdf",width = 18,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/T.cell.exhaustion_group.VlnPlot.png",width = 18,height = 10)

p <- get_sig_ana(p,T.cell.exhaustion~split,group.by='ident',prefix="T.cell.exhaustion_group",map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
get_sig_ana(p$data,T.cell.exhaustion~split,group.by='ident',prefix="T.cell.exhaustion_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/T.cell.exhaustion_group_sig_ana.VlnPlot.pdf",width = 18,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/T.cell.exhaustion_group_sig_ana.VlnPlot.png",width = 18,height = 10)

df<-p$data
df$cell_type = factor(df$ident,levels = c('CD4T',"CD8T","T proliferating"))
gg = ggplot(data=df,mapping=aes(x=cell_type,y=df$`T.cell.exhaustion`))+
  geom_jitter(aes(colour  = cell_type),
              width = 0.1,  size=2.5, stroke=0.5, shape=16,alpha=10) +
  geom_violin(alpha=0.65, linewidth=1, width=0.5,adjust = 0.9) +
  geom_boxplot(width=0.3, linewidth=1,outlier.color = NA, fill=NA) +
  scale_colour_manual(values=c("#749B58FF","#F0E685FF","#466983FF"))+
  stat_compare_means(method = "wilcox.test",ref.group=".all.",  label = "p.format") +
  # stat_compare_means(method = "wilcox.test", comparisons = combn(celltype, 2, simplify = FALSE), label = "p.format") +
  theme_bw()+theme(panel.grid=element_blank(),axis.text.x=element_text(color='black', angle=45, hjust=1))+
  geom_hline(yintercept  = mean(df$`T.cell.exhaustion`),
             linetype = "dashed", color = "red4", size = 1)+
  labs(x="",y = "T.cell.exhaustion")
# +
#    facet_wrap(~ celltype, nrow = 2)  # 按照 cell_type 分面
gg
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/T.cell.exhaustion-celltype.pdf",gg,width = 7,height = 5)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/T.cell.exhaustion-celltype.png",gg,width = 7,height = 5)


genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/script/new_geneset.csv",sep = "\t",head=TRUE)

gene<-c(genesets$NK.cell.exhaustion)
gene <- gene[gene != ""]
gene<-list(gene)

adata_sub<-subset(adata,idents=c("NK"))
Idents(adata_sub)<-adata_sub$group
adata_sub<-AddModuleScore(
    object=adata_sub,
    features=gene,
    ctr1=100, #默认值是100
    name='select_geneset'
)

colnames(adata_sub@meta.data)
colnames(adata_sub@meta.data)[23]<-'NK.cell.exhaustion'

p <- VlnPlot(adata_sub,features='NK.cell.exhaustion',pt.size = 0)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/NK.cell.exhaustion_group.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/NK.cell.exhaustion_group.png",width = 10,height = 10)

p <- get_sig_ana(p,NK.cell.exhaustion~ident,group.by='ident',prefix="NK.cell.exhaustion_group",map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))+scale_y_continuous(limits=c(0,2.2))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
# get_sig_ana(p$data,NK.cell.exhaustion~ident,group.by='ident',prefix="NK.cell.exhaustion_group",outdir="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/NK.cell.exhaustion_group_sig_ana.VlnPlot.pdf",width = 10,height = 11)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/0928_AddModuleScore/NK.cell.exhaustion_group_sig_ana.VlnPlot.png",width = 10,height = 11)