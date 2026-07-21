library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(readxl)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/mapping_raw/Main_mapping.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
# adata <- subset(adata,subset=(cell_type %in% c("B cells", "Plasma cells")))

cytokine<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/cytokine_genes.csv",sep = ",",head=FALSE)
inflammatory<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/inflammatory_genes.csv",sep = ",",head=FALSE)
mTOR_pathway<-c('MTOR','RPTOR','RICTOR','RHEB','TSC1','TSC2')
IL6<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/il6_gene.csv",sep = ",",head=TRUE)
#转换成list
gene<-c(as.list(cytokine),as.list(inflammatory),list(mTOR_pathway),IL6)
names(gene)<-c('cytokine','inflammatory','mTOR_pathway','KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION')
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[16]<-'cytokine'
colnames(adata@meta.data)[17]<-'inflammatory'
colnames(adata@meta.data)[18]<-'mTOR_pathway'
colnames(adata@meta.data)[19]<-'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION'
outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/AddModuleScore/"

Idents(adata)<-adata$raw_cell_type
for (pathway in c('cytokine','inflammatory','mTOR_pathway','KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION')){
	print(pathway)
p <- VlnPlot(adata,features=pathway,pt.size = 0)+scale_fill_manual(values = c("B cells"="#ea7070","Basophils"="#fdc4b6","Erythroid cells"="#e59572","Megakaryocyte"="#2694ab","Myeloid cells"="#96ceb4","Neutrophils"="#ffeead","Plasma cells"="#d9534f","T/NK"="#ffad60"))
ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c('HD'='#E76BF3','P'='#00BF7D','P-after'='#A3A500'))
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c('HD'='#E76BF3','P'='#00BF7D','P-after'='#A3A500'))
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)
}
metadata<-adata@meta.data
file_path<-paste0(outdir,"Main_pathway.csv")
write.table(metadata, file = file_path, sep = ",", append = TRUE, row.names = FALSE)


adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/mapping_raw/Main_mapping.rds")
Idents(adata)<-adata$raw_cell_type
bdata<-subset(adata,idents = "Myeloid cells")
Idents(bdata)<-bdata$sub_cell_type
bdata<-AddModuleScore(
    object=bdata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(bdata@meta.data)
colnames(bdata@meta.data)[16]<-'cytokine'
colnames(bdata@meta.data)[17]<-'inflammatory'
colnames(bdata@meta.data)[18]<-'mTOR_pathway'
colnames(bdata@meta.data)[19]<-'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION'

for (pathway in c('cytokine','inflammatory','mTOR_pathway','KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION')){
	print(pathway)
p <- VlnPlot(bdata,features=pathway,pt.size = 0)+scale_fill_manual(values = c("DUSP1 Classical monocyte"="#1f77b4","ISG15 Classical monocyte"="#ff7f0e","LRMDA Classical monocyte"="#279e68","Non-classical monocyte"="#d62728","VCAN Classical monocyte"="#aa40fc","cDCs"="#8c564b","pDCS"="#e377c2"))
ggsave(paste0(outdir,"Myeloid-",pathway,"-VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,"Myeloid-",pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(bdata,features=pathway,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c('HD'='#E76BF3','P'='#00BF7D','P-after'='#A3A500'))
ggsave(paste0(outdir,"Myeloid_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Myeloid_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Myeloid_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c('HD'='#E76BF3','P'='#00BF7D','P-after'='#A3A500'))
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Myeloid_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Myeloid_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Myeloid_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)
}
metadata<-bdata@meta.data
file_path<-paste0(outdir,"Myeloid_pathway.csv")
write.table(metadata, file = file_path, sep = ",", append = TRUE, row.names = FALSE)

# adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/mapping_raw/Main_mapping.rds")
# Idents(adata)<-adata$raw_cell_type
bdata<-subset(adata,sub_cell_type %in% c("IgA Plasma","IgG Plasma C1","IgG Plasma C2","IgM Plasma"))
Idents(bdata)<-bdata$sub_cell_type
bdata<-AddModuleScore(
    object=bdata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(bdata@meta.data)
colnames(bdata@meta.data)[16]<-'cytokine'
colnames(bdata@meta.data)[17]<-'inflammatory'
colnames(bdata@meta.data)[18]<-'mTOR_pathway'
colnames(bdata@meta.data)[19]<-'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION'

for (pathway in c('cytokine','inflammatory','mTOR_pathway','KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION')){
	print(pathway)
p <- VlnPlot(bdata,features=pathway,pt.size = 0)+scale_fill_manual(values = c("IgA Plasma"="#ff7f0e","IgG Plasma C1"="#279e68","IgG Plasma C2"="#d62728","IgM Plasma"="#aa40fc"))
ggsave(paste0(outdir,"Plasmacells-",pathway,"-VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,"Plasmacells-",pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(bdata,features=pathway,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c('HD'='#E76BF3','P'='#00BF7D','P-after'='#A3A500'))
ggsave(paste0(outdir,"Plasmacells_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Plasmacells_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Plasmacells_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c('HD'='#E76BF3','P'='#00BF7D','P-after'='#A3A500'))
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Plasmacells_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Plasmacells_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Plasmacells_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)
}
metadata<-bdata@meta.data
file_path<-paste0(outdir,"Plasma_pathway.csv")
write.table(metadata, file = file_path, sep = ",", append = TRUE, row.names = FALSE)