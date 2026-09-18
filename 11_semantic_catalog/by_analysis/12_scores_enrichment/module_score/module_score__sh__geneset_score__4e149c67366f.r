library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(pheatmap)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/sc_mono_anno.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)

geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/pathway.csv",sep = "\t",head=TRUE)
geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWAE20251216/03.Mission/script/cytokine_genes.csv",sep = "\t",head=TRUE)

#转换成list
gene<-as.list(geneset)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100 #默认值是100
)

colnames(adata@meta.data)
colnames(adata@meta.data)[14]<-'HALLMARK_INFLAMMATORY_RESPONSE'
colnames(adata@meta.data)[15]<-'HALLMARK_GLYCOLYSIS'
colnames(adata@meta.data)[16]<-'CYTOKINE_RECEPTOR_INTERACTION'
colnames(adata@meta.data)[14]<-'CYTOKINE'
mydata<-FetchData(adata,vars=c("cell_type","group","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_GLYCOLYSIS","CYTOKINE_RECEPTOR_INTERACTION"))
mean_scores <- mydata %>%
  group_by(cell_type, group) %>%
  summarise(across(everything(), mean))

score_select<-mean_scores$CYTOKINE
score_select <- matrix(score_select, nrow = 6, ncol = 3, byrow = TRUE)
row.names(score_select)<-c("Non-classical monocyte","PTTG1 classical monocyte","ZFP36L1 classical monocyte","S100A8 classical monocyte","Classical monocyte","Dendritic cell")
colnames(score_select)<-c("LAHS","NLAHS","NC")
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/AddModuleScore/CYTOKINE.pdf",width = 7,height = 8)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/AddModuleScore/CYTOKINE.png",width = 480,height = 540)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()

write.table(mean_scores,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/AddModuleScore/CYTOKINE_AddModuleScore.csv", sep = ",", append = TRUE, row.names = FALSE)


# data <- cor(mydata[,3:5], method = "pearson")
# cor_scores <- mydata %>%
#   group_by(cell_type, group) %>%
#   summarise(cor_HALLMARK_INFLAMMATORY_RESPONSE_HALLMARK_GLYCOLYSIS = cor(HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_GLYCOLYSIS),
#             cor_HALLMARK_INFLAMMATORY_RESPONSE_CYTOKINE_RECEPTOR_INTERACTION = cor(HALLMARK_INFLAMMATORY_RESPONSE, CYTOKINE_RECEPTOR_INTERACTION),
#             cor_HALLMARK_GLYCOLYSIS_CYTOKINE_RECEPTOR_INTERACTION = cor(HALLMARK_GLYCOLYSIS, CYTOKINE_RECEPTOR_INTERACTION))
# write.table(cor_scores,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/AddModuleScore/cor_scores.csv", sep = ",", append = TRUE, row.names = FALSE)

# write.table(mydata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/12.Mission/AddModuleScore/all_scores.csv", sep = ",", append = TRUE, row.names = FALSE)
