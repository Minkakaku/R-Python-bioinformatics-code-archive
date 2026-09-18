library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(readxl)
source("/home/sh/sc_pipline/geneset_score/get_sig_ana.R")

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/SC2509265318-02/02.Mission/Chondrocyte/Chondrocyte.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)
# adata <- subset(adata,subset=(cell_type %in% c("B cells", "Plasma cells")))
# Idents(adata)<-adata$raw_cell_type
genset<-read.table("/mnt/NexData_file/home/sh/Projects/SC2509265318-02/02.Mission/script/genset.csv",sep = "\t",head=FALSE)
genset$V2 <- strsplit(genset$V2, ",")
gene <- setNames(genset$V2, genset$V1)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100, #默认值是100
     name='select_geneset'
)
colnames(adata@meta.data)
colnames(adata@meta.data)[15]<-'Aging_GeneSet'
outdir <- "/mnt/NexData_file/home/sh/Projects/SC2509265318-02/02.Mission/sc_score/"

for (pathway in c('Aging_GeneSet')){
	print(pathway)
p <- VlnPlot(adata,features=pathway,pt.size = 0)+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave(paste0(outdir,"Chondrocyte-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Chondrocyte-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- VlnPlot(adata,features=pathway,split.by = 'group',pt.size = 0)+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
ggsave(paste0(outdir,"Chondrocyte_group-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Chondrocyte_group-",pathway,"-VlnPlot.png"),width = 15,height = 10)

p <- get_sig_ana(p,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Chondrocyte_group-",pathway),map_sig_level=T,only_sig = F,method = "wilcox.test")+geom_boxplot(width=0.1,col="black",position = position_dodge(width = 0.9),outlier.shape = NA)
get_sig_ana(p$data,as.formula(paste0(pathway,'~split')),group.by='ident',prefix=paste0("Chondrocyte_group-",pathway),outdir=outdir,map_sig_level=T,only_sig = F,method = "wilcox.test")
ggsave(paste0(outdir,"Chondrocyte_group_sig_ana-",pathway,"-VlnPlot.pdf"),width = 15,height = 10)
ggsave(paste0(outdir,"Chondrocyte_group_sig_ana-",pathway,"-VlnPlot.png"),width = 15,height = 10)
}
metadata<-adata@meta.data
file_path<-paste0(outdir,"Chondrocyte_pathway.csv")
write.table(metadata, file = file_path, sep = ",", append = TRUE, row.names = FALSE)



# adata$cell_type_group<-paste0(adata$raw_cell_type,"_",adata$group)
# p<-DotPlot(adata, features = IL6$V1,group.by = 'cell_type_group',dot.scale = 12,scale = FALSE) + 
# 		coord_flip() + #翻转
# 		theme(panel.grid = element_blank(), 
#         axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
# 		guides(size = guide_legend("Percent Expression") )+ #legend
# 		scale_color_gradientn(colours = c('#375093','#e8edf1','#831a21')) #颜色
# ggsave(paste0(outdir,'BIOCARTA_IL6_PATHWAY_DotPlot.pdf'),width = 12,height = 12,bg = 'white')
# ggsave(paste0(outdir,'BIOCARTA_IL6_PATHWAY_DotPlot.png'),width = 12,height = 12,bg = 'white')
# p<-DotPlot(adata, features = IL6_Signaling_Score$IL6_Signaling_Score,group.by = 'cell_type_group',dot.scale = 12,scale = FALSE) + 
# 		coord_flip() + #翻转
# 		theme(panel.grid = element_blank(), 
#         axis.text.x=element_text(angle = 45, hjust = 1,vjust=1))+ #轴标签labs(x=NULL,y=NULL) + 
# 		guides(size = guide_legend("Percent Expression") )+ #legend
# 		scale_color_gradientn(colours = c('#375093','#e8edf1','#831a21')) #颜色
# ggsave(paste0(outdir,'IL6_Signaling_Score_DotPlot.pdf'),width = 12,height = 12,bg = 'white')
# ggsave(paste0(outdir,'IL6_Signaling_Score_DotPlot.png'),width = 12,height = 12,bg = 'white')