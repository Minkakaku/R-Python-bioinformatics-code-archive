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

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251249/03.Mission/mapping_raw/Main_mapping.rds")
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
# Idents(adata) <- "annotation"
adata<-subset(adata,raw_cell_type %in% c("T cells","B cells"))
genesets<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251249/03.Mission/script/new_geneset.csv",sep = "\t",head=TRUE)

gene1<-c(genesets$GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION)
gene1 <- gene1[gene1 != ""]
gene2<-c(genesets$T.cell.exhaustion)
gene2 <- gene2[gene2 != ""]
gene3<-c(genesets$cytokine)
gene3 <- gene3[gene3 != ""]
gene<-list(gene1,gene2,gene3)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
    name='select_geneset'
)

colnames(adata@meta.data)
colnames(adata@meta.data)[16]<-'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION'
colnames(adata@meta.data)[17]<-'T.cell.exhaustion'
colnames(adata@meta.data)[18]<-'Cytokine'
Idents(adata)<-adata$annotation
p <- VlnPlot(adata,features='Cytokine',pt.size = 0)+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA) #+scale_fill_manual(values = c('CD4_Tcm'="#A6CEE3","CD4_Tem"="#3F8EAA","CD4_TNaive"="#79C360","CD4_Treg"="#B89B74","CD8_Tem"="#E52829","CD8_Trm"="#FDB762","CD8_Tex"="#ED8F47","CD8_MAIT"="#9471B4","T proliferating"="#DDD399","NK"="#B15928"))
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251249/03.Mission/AddModuleScore/Cytokine.VlnPlot.pdf",width = 10,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251249/03.Mission/AddModuleScore/Cytokine.VlnPlot.png",width = 10,height = 10)

p <- VlnPlot(adata,features='Cytokine',split.by = 'group',pt.size = 0)+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA) #+scale_fill_manual(values = c(LAHS='#1f77b4',NLAHS='#ff7f0e',NC='#279e68'))
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251249/03.Mission/AddModuleScore/Cytokine_group.VlnPlot.pdf",width = 18,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251249/03.Mission/AddModuleScore/Cytokine_group.VlnPlot.png",width = 18,height = 10)

p <- get_sig_ana(p,Cytokine~split,group.by='ident',prefix="Cytokine_group",map_sig_level=T,only_sig = F,method = "wilcox.test")
get_sig_ana(p$data,Cytokine~split,group.by='ident',prefix="Cytokine_group",outdir="/mnt/NexData_file/home/sh/Projects/FWAE20251249/03.Mission/AddModuleScore/",map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251249/03.Mission/AddModuleScore/Cytokine_group_sig_ana.VlnPlot.pdf",width = 18,height = 10)
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251249/03.Mission/AddModuleScore/Cytokine_group_sig_ana.VlnPlot.png",width = 18,height = 10)

df<-p$data
df$cell_type = factor(df$ident,levels = c("B cells","Cd4T","Cd8T","GDT","MAIT","NK Cells","Plasma cells","Treg","pDC"))
gg = ggplot(data=df,mapping=aes(x=cell_type,y=df$`Cytokine`))+
  geom_jitter(aes(colour  = cell_type),
              width = 0.1,  size=2.5, stroke=0.5, shape=16,alpha=10) +
  geom_violin(alpha=0.65, linewidth=1, width=0.5,adjust = 0.9) +
  geom_boxplot(width=0.3, linewidth=1,outlier.color = NA, fill=NA) +
  # scale_colour_manual(values=c("#A6CEE3","#3F8EAA","#79C360","#B89B74", "#E52829", "#FDB762", "#ED8F47", "#9471B4","#DDD399","#B15928"))+
  stat_compare_means(method = "wilcox.test",ref.group=".all.",  label = "p.format") +
  # stat_compare_means(method = "wilcox.test", comparisons = combn(celltype, 2, simplify = FALSE), label = "p.format") +
  theme_bw()+theme(panel.grid=element_blank(),axis.text.x=element_text(color='black', angle=45, hjust=1))+
  geom_hline(yintercept  = mean(df$`Cytokine`),
             linetype = "dashed", color = "red4", size = 1)+
  labs(x="",y = "Cytokine")
# +
#    facet_wrap(~ celltype, nrow = 2)  # 按照 cell_type 分面
gg
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251249/03.Mission/AddModuleScore/Cytokine-celltype.pdf",gg,width = 12,height = 5)
ggsave("/mnt/NexData_file/home/sh/Projects/FWAE20251249/03.Mission/AddModuleScore/Cytokine-celltype.png",gg,width = 12,height = 5)