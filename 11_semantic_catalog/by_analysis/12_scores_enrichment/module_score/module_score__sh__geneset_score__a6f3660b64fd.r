library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(readxl)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_A_umap/Main_mapping.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
# adata <- subset(adata,subset=(cell_type %in% c("B cells", "Plasma cells")))
Idents(adata)<-adata$raw_cell_type
cytokine<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/cytokine_genes.csv",sep = ",",head=FALSE)
inflammatory<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/inflammatory_genes.csv",sep = ",",head=FALSE)
mTOR_pathway<-c('MTOR','RPTOR','RICTOR','RHEB','TSC1','TSC2')
IL6<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/il6_gene.csv",sep = ",",head=FALSE)
IL6_Signaling_Score<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/script/IL6_Signaling_Score.csv",sep = ",",head=TRUE)
#转换成list
gene<-c(as.list(cytokine),as.list(inflammatory),list(mTOR_pathway),IL6,IL6_Signaling_Score)
names(gene)<-c('cytokine','inflammatory','mTOR_pathway','IL6','IL6_Signaling_Score')
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
colnames(adata@meta.data)[19]<-'IL6'
colnames(adata@meta.data)[20]<-'IL6_Signaling_Score'
outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_E_scoreplot/Main_new/"

for (pathway in c('cytokine','inflammatory','mTOR_pathway','IL6','IL6_Signaling_Score')){
	print(pathway)
p <- VlnPlot(adata,features=pathway,pt.size = 0)
ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0)
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)
}
metadata<-adata@meta.data
file_path<-paste0(outdir,"Main_pathway.csv")
write.table(metadata, file = file_path, sep = ",", append = TRUE, row.names = FALSE)

adata$cell_type_group<-paste0(adata$raw_cell_type,"_",adata$group)
p<-DotPlot(adata, features = IL6$V1,group.by = 'cell_type_group',dot.scale = 12,scale = FALSE) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#375093','#e8edf1','#831a21')) #颜色
ggsave(paste0(outdir,'BIOCARTA_IL6_PATHWAY_DotPlot.pdf'),width = 12,height = 12,bg = 'white')
ggsave(paste0(outdir,'BIOCARTA_IL6_PATHWAY_DotPlot.png'),width = 12,height = 12,bg = 'white')
p<-DotPlot(adata, features = IL6_Signaling_Score$IL6_Signaling_Score,group.by = 'cell_type_group',dot.scale = 12,scale = FALSE) + 
		coord_flip() + #翻转
		theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
		guides(size = guide_legend("Percent Expression") )+ #legend
		scale_color_gradientn(colours = c('#375093','#e8edf1','#831a21')) #颜色
ggsave(paste0(outdir,'IL6_Signaling_Score_DotPlot.pdf'),width = 12,height = 12,bg = 'white')
ggsave(paste0(outdir,'IL6_Signaling_Score_DotPlot.png'),width = 12,height = 12,bg = 'white')


MIF<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/mif_receptor.csv",sep = ",",head=FALSE)
IL6<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/il6_gene.csv",sep = ",",head=FALSE)
gene<-list(c(MIF$V1),c(IL6$V1))

bdata<-subset(adata,sub_cell_type %in% c("IgA Plasma","IgG Plasma C1","IgG Plasma C2","IgM Plasma"))
Idents(bdata)<-bdata$sub_cell_type
bdata<-AddModuleScore(
    object=bdata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(bdata@meta.data)
colnames(bdata@meta.data)[20]<-'MIF'
colnames(bdata@meta.data)[21]<-'IL6'
outdir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_E_scoreplot/Plasma/"
for (pathway in c('MIF','IL6')){
	print(pathway)
p <- VlnPlot(bdata,features=pathway,pt.size = 0)+scale_fill_manual(values = c("IgA Plasma"="#8595e1","IgG Plasma C1"="#b5bbe3","IgG Plasma C2"="#e6afb9","IgM Plasma"="#e07b91"))
ggsave(paste0(outdir,"Plasmacells-",pathway,"-VlnPlot.pdf"),width = 10,height = 10)
ggsave(paste0(outdir,"Plasmacells-",pathway,"-VlnPlot.png"),width = 10,height = 10)

p <- VlnPlot(bdata,features=pathway,split.by = 'group',pt.size = 0)
ggsave(paste0(outdir,"Plasmacells_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Plasmacells_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Plasmacells_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Plasmacells_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Plasmacells_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Plasmacells_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)
}

for (pathway in c('MIF')){
	print(pathway)
p <- VlnPlot(bdata,features=pathway,group.by = 'group',pt.size = 0)+geom_boxplot(width=0.05,col="black",fill="white",outlier.shape = NA)
ggsave(paste0(outdir,"Plasma_group-",pathway,"-VlnPlot.pdf"),width = 7,height = 10)
ggsave(paste0(outdir,"Plasma_group-",pathway,"-VlnPlot.png"),width = 7,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~ident')),group.by='ident',prefix=paste0("Plasma_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test",reset_y=T)+ylim(-0.3,1.5)
# get_sig_ana(p$data,as.formula(paste0(pathway,'~ident')),group.by='ident',prefix=paste0("Plasma_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Plasma_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 7,height = 15)
ggsave(paste0(outdir,"Plasma_group_sig_ana-",pathway,"-VlnPlot.png"),width = 7,height = 15)
}
metadata<-bdata@meta.data
file_path<-paste0(outdir,"Plasma_pathway.csv")
write.table(metadata, file = file_path, sep = ",", append = TRUE, row.names = FALSE)



#########################
adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/script/Main_mapping_clean.rds")
levels(adata$cell_type)<-c("Myeloid cells","Myeloid cells","Myeloid cells","Plasma cells","Plasma cells","Plasma cells","Plasma cells","Plasma cells","Megakaryocytes","B cells","B cells","B cells","NK cells","NKT cells","Myeloid cells","pDCs","T cells","Myeloid cells")
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)

Idents(adata)<-adata$cell_type
cytokine<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/cytokine_genes.csv",sep = ",",head=FALSE)
inflammatory<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/inflammatory_genes.csv",sep = ",",head=FALSE)
mTOR_pathway<-c('MTOR','RPTOR','RICTOR','RHEB','TSC1','TSC2')
IL6<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/04.Mission/script/il6_gene.csv",sep = ",",head=FALSE)
#转换成list
gene<-c(as.list(cytokine),as.list(inflammatory),list(mTOR_pathway),IL6)
names(gene)<-c('cytokine','inflammatory','mTOR_pathway','IL6')
adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[17]<-'cytokine'
colnames(adata@meta.data)[18]<-'inflammatory'
colnames(adata@meta.data)[19]<-'mTOR_pathway'
colnames(adata@meta.data)[20]<-'IL6'
outdir <- "/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/addplot_scoreplot/"

for (pathway in c('cytokine','inflammatory','mTOR_pathway','IL6')){
	print(pathway)
p <- VlnPlot(adata,features=pathway,pt.size = 0)
ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0)
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Main_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Main_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)
}
metadata<-adata@meta.data
file_path<-paste0(outdir,"Main_pathway.csv")
write.table(metadata, file = file_path, sep = ",", append = TRUE, row.names = FALSE)