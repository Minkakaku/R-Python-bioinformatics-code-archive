library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(readxl)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/02.Mission/mapping_raw/Main_mapping.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
# adata <- subset(adata,subset=(cell_type %in% c("B cells", "Plasma cells")))
levels(adata$group)<-c('N','P','AID')

MIF<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
IL6<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/il6_gene.csv",sep = ",",head=TRUE)
#转换成list
gene<-c(as.list(MIF),as.list(IL6))
names(gene)<-c('MIF','IL6')
outdir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/fig3_B_vlnplot/"

bdata<-subset(adata,idents = "Plasma cells")
Idents(bdata)<-bdata$sub_cell_type
bdata<-AddModuleScore(
    object=bdata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(bdata@meta.data)
colnames(bdata@meta.data)[16]<-'MIF'
colnames(bdata@meta.data)[17]<-'IL6'

for (pathway in c('MIF','IL6')){
	print(pathway)
p <- VlnPlot(bdata,features=pathway,pt.size = 0)+scale_fill_manual(values = c("IgA Plasma"='#1f77b4',"IgG Plasma C1"='#ff7f0e',"IgG Plasma C2"='#279e68',"IgM Plasma C1"='#d62728',"IgM Plasma C2"='#aa40fc'))
ggsave(paste0(outdir,"Plasmacells-",pathway,"-VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,"Plasmacells-",pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(bdata,features=pathway,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c('N'='#1f77b4','P'='#ff7f0e','AID'='#279e68'))
ggsave(paste0(outdir,"Plasmacells_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Plasmacells_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Plasmacells_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(N='#1f77b4',P='#ff7f0e','AID'='#279e68'))
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Plasmacells_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Plasmacells_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Plasmacells_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)


p <- VlnPlot(bdata,features=pathway,group.by = 'group',pt.size = 0)+scale_fill_manual(values = c('N'='#1f77b4','P'='#ff7f0e','AID'='#279e68'))+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9), linewidth=1,outlier.color = NA, fill=NA)
ggsave(paste0(outdir,"Plasma_group-",pathway,"-VlnPlot.pdf"),width = 8,height = 10)
ggsave(paste0(outdir,"Plasma_group-",pathway,"-VlnPlot.png"),width = 8,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~ident')),group.by='ident',prefix=paste0("Plasma_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)+scale_fill_manual(values = c(N='#1f77b4',P='#ff7f0e','AID'='#279e68'))+scale_y_continuous(limits=c(-0.2,0.4)) #+scale_y_continuous(limits=c(-0.25,1.1))
# get_sig_ana(p$data,as.formula(paste0(pathway,'~ident')),group.by='ident',prefix=paste0("Plasma_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Plasma_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 10,height = 15)
ggsave(paste0(outdir,"Plasma_group_sig_ana-",pathway,"-VlnPlot.png"),width = 10,height = 15)
}