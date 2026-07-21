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

cytokine<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/cytokine_genes.csv",sep = ",",head=FALSE)
inflammatory<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/inflammatory_genes.csv",sep = ",",head=FALSE)
#转换成list
gene<-c(as.list(cytokine),as.list(inflammatory))
names(gene)<-c('cytokine','inflammatory')
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[16]<-'cytokine'
colnames(adata@meta.data)[17]<-'inflammatory'
outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/AddModuleScore/"

for (pathway in c('cytokine','inflammatory')){
	print(pathway)
p <- VlnPlot(adata,features=pathway,pt.size = 0)+scale_fill_manual(values = c('B cells'='#1f77b4','Erythroid cells'='#ff7f0e','HSPCs'='#279e68','Megakaryocytes'='#d62728','Myeloid cells'='#aa40fc','NK cells'='#8c564b','Neutrophils'='#e377c2','Plasma cells'='#b5bd61','T cells'='#17becf'))
ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c('Negative'='#1f77b4','Positive'='#ff7f0e','Rheumatic Diseases'='#279e68'))
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(Negative='#1f77b4',Positive='#ff7f0e','Rheumatic Diseases'='#279e68'))
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)

# mydata<-FetchData(adata,vars=c("UMAP_1","UMAP_2","cytokine","group"))
# a1<-ggplot(mydata[which(mydata$group=="Negative"),],aes(x=UMAP_1,y=UMAP_2,colour=cytokine)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Negative")
# a2<-ggplot(mydata[which(mydata$group=="Positive"),],aes(x=UMAP_1,y=UMAP_2,colour=cytokine)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),legend.position = 'none',
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Positive")
# a3<-ggplot(mydata[which(mydata$group=="Rheumatic Diseases"),],aes(x=UMAP_1,y=UMAP_2,colour=cytokine)) + 
# 		geom_point(size=1) + 
# 		scale_color_gradientn(values=seq(0,1,0.2),colours = c('#ffeee5','#d9950b','#f1681e','#ce3d32','#eb2d4a')) + 
# 		theme_bw() + theme(panel.border=element_blank(),panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),axis.line=element_line(colour="black"))+labs(title="Rheumatic Diseases")
# a<-a1+a2+a3
# ggsave(paste0(outdir,"Main-",pathway,"-UMAP.pdf"),width = 18,height = 6)
# ggsave(paste0(outdir,"Main-",pathway,"-UMAP.png"),width = 18,height = 6)

}


bdata<-subset(adata,idents = "B cells")
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

for (pathway in c('cytokine','inflammatory')){
	print(pathway)
p <- VlnPlot(bdata,features=pathway,pt.size = 0)+scale_fill_manual(values = c("Memory B C1"='#1f77b4',"Memory B C2"='#ff7f0e',"Naive B"='#279e68'))
ggsave(paste0(outdir,"Bcells-",pathway,"-VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,"Bcells-",pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(bdata,features=pathway,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c('Negative'='#1f77b4','Positive'='#ff7f0e','Rheumatic Diseases'='#279e68'))
ggsave(paste0(outdir,"Bcells_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Bcells_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Bcells_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(Negative='#1f77b4',Positive='#ff7f0e','Rheumatic Diseases'='#279e68'))
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Bcells_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Bcells_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Bcells_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)

}


bdata<-subset(adata,idents = "Plasma cells")
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

for (pathway in c('cytokine','inflammatory')){
	print(pathway)
p <- VlnPlot(bdata,features=pathway,pt.size = 0)+scale_fill_manual(values = c("IgA Plasma"='#1f77b4',"IgG Plasma C1"='#ff7f0e',"IgG Plasma C2"='#279e68',"IgM Plasma C1"='#d62728',"IgM Plasma C2"='#aa40fc'))
ggsave(paste0(outdir,"Plasmacells-",pathway,"-VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,"Plasmacells-",pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(bdata,features=pathway,split.by = 'group',pt.size = 0)+scale_fill_manual(values = c('Negative'='#1f77b4','Positive'='#ff7f0e','Rheumatic Diseases'='#279e68'))
ggsave(paste0(outdir,"Plasmacells_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Plasmacells_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Plasmacells_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+scale_fill_manual(values = c(Negative='#1f77b4',Positive='#ff7f0e','Rheumatic Diseases'='#279e68'))
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Plasmacells_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Plasmacells_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Plasmacells_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)

}