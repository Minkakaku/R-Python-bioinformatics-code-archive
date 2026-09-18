library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
library(ggplot2)
library(pheatmap)

adata<-readRDS("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/Monocyte/sc_mono.rds")
# adata <- JoinLayers(adata)
# adata <- subset(adata,subset=(sample %in% c("LXZ", "ZGP0408","ZCX0503","YLY0701", "ZLY0717","FJH0806-W","ZM0814","HXL0726","Y0726")))
adata$group<-factor(adata$group)
adata$sample<-factor(adata$sample)

geneset<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/script/pathway.csv",sep = "\t",head=TRUE)
#转换成list
gene<-as.list(geneset)

adata<-AddModuleScore(
    object=adata,
    features=gene,
    ctr1=100 #默认值是100
)

colnames(adata@meta.data)
colnames(adata@meta.data)[16]<-'HALLMARK_INFLAMMATORY_RESPONSE'
colnames(adata@meta.data)[17]<-'HALLMARK_GLYCOLYSIS'
colnames(adata@meta.data)[18]<-'CYTOKINE_RECEPTOR_INTERACTION'
mydata<-FetchData(adata,vars=c("annotation","group","HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_GLYCOLYSIS","CYTOKINE_RECEPTOR_INTERACTION"))
mean_scores <- mydata %>%
  group_by(annotation, group) %>%
  summarise(across(everything(), mean))
new_row <- data.frame(annotation="TPX2 Monocyte",group="NC",HALLMARK_INFLAMMATORY_RESPONSE= NA,HALLMARK_GLYCOLYSIS= NA,CYTOKINE_RECEPTOR_INTERACTION= NA)
mean_scores <- rbind(mean_scores, new_row)

score_select<-mean_scores$CYTOKINE_RECEPTOR_INTERACTION
score_select <- matrix(score_select, nrow = 4, ncol = 3, byrow = TRUE)
row.names(score_select)<-c("CD14 Monocyte","FCGR3A Monocyte","Dendritic cell","TPX2 Monocyte")
colnames(score_select)<-c("LAHS","NLAHS","NC")
pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/AddModuleScore/CYTOKINE_RECEPTOR_INTERACTION.pdf",width = 7,height = 8)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/AddModuleScore/CYTOKINE_RECEPTOR_INTERACTION.png",width = 480,height = 540)
pheatmap::pheatmap(score_select, show_colnames = T,cluster_cols=F,cluster_rows=F,
                   scale = "none",angle_col = "45",
                   color = c(colorRampPalette(colors = c("#375093","#fbe2ab","#9D000A"))(50)))
dev.off()

write.table(mean_scores,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/AddModuleScore/AddModuleScore.csv", sep = ",", append = TRUE, row.names = FALSE)


data <- cor(mydata[,3:5], method = "pearson")
cor_scores <- mydata %>%
  group_by(annotation, group) %>%
  summarise(cor_HALLMARK_INFLAMMATORY_RESPONSE_HALLMARK_GLYCOLYSIS = cor(HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_GLYCOLYSIS),
            cor_HALLMARK_INFLAMMATORY_RESPONSE_CYTOKINE_RECEPTOR_INTERACTION = cor(HALLMARK_INFLAMMATORY_RESPONSE, CYTOKINE_RECEPTOR_INTERACTION),
            cor_HALLMARK_GLYCOLYSIS_CYTOKINE_RECEPTOR_INTERACTION = cor(HALLMARK_GLYCOLYSIS, CYTOKINE_RECEPTOR_INTERACTION))
write.table(cor_scores,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/AddModuleScore/cor_scores.csv", sep = ",", append = TRUE, row.names = FALSE)

write.table(mydata,"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/10.Mission/AddModuleScore/all_scores.csv", sep = ",", append = TRUE, row.names = FALSE)


# x_col <- 'HIF1A'
# y_cols <- setdiff(gene$x, x_col)
# for (y_col in y_cols) {
#   # 绘制散点图
#   p <- ggplot(genedata, aes_string(x = x_col, y = y_col)) +
#     geom_point(size = 2, alpha = 0.7, shape= 1, color = "#82b4ad") +  # 散点图
#     geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +  # 添加线性拟合线
#     labs(title = paste(x_col, "vs", y_col, "Cor"),
#          x = x_col,
#          y = y_col) +
#     theme_bw() +  # 使用简洁主题
#     theme(plot.title = element_text(hjust = 0.5),  # 居中标题
#           panel.grid.major = element_blank(),  # 添加网格线
#           panel.grid.minor = element_blank())  # 去掉次要网格线
  
#   # 打印图形
#   ggsave(paste("/data/NAS02/sh/Projects/FWSC2024-192/07.Mission/cor/","scatter_", y_col, ".png", sep = ""), plot = p, width = 6, height = 4)
#   ggsave(paste("/data/NAS02/sh/Projects/FWSC2024-192/07.Mission/cor/","scatter_", y_col, ".pdf", sep = ""), plot = p, width = 6, height = 4)

# }