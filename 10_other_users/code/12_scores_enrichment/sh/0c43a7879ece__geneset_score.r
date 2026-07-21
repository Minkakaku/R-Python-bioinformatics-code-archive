library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(readxl)
source("/mnt/NexData_file/share/scRNA-tools/diffgene_enrich/get_sig_ana.R")

#########################
adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/script/Main_mapping_clean.rds")
levels(adata$cell_type)<-c("Myeloid cells","Myeloid cells","Myeloid cells","Plasma cells","Plasma cells","Plasma cells","Plasma cells","Plasma cells","Megakaryocytes","B cells","B cells","B cells","NK cells","NKT cells","Myeloid cells","pDCs","T cells","Myeloid cells")
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
Idents(adata)<-adata$cell_type

MIF<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
gene<-c(as.list(MIF))

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[17]<-'MIF'
outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/19.Mission/result6/"

for (pathway in c('MIF')){
	print(pathway)
# p <- VlnPlot(adata,features=pathway,pt.size = 0)
# ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
# ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0)
# ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
# ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"1.Main_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"1.Main_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)
}
# metadata<-adata@meta.data
# file_path<-paste0(outdir,"Main_pathway.csv")
# write.table(metadata, file = file_path, sep = ",", append = TRUE, row.names = FALSE)



adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_A_umap/Main_mapping.rds")
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
Idents(adata)<-adata$raw_cell_type
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[16]<-'MIF'
for (pathway in c('MIF')){
	print(pathway)
p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"2.Main_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"2.Main_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)
}