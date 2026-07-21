library(dplyr)
library(ggplot2)
library(RColorBrewer)

GO_result<-read.table("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/07.Mission/Fibroblasts/clusterADC/GO_Biological_Process_2023.Human.enrichr.reports.txt",sep="\t",header=TRUE,fill = TRUE,quote = "",fileEncoding = "UTF-8")
# GO_result<-GO_result[1:10,]
top10<- GO_result
#2.论文画法：
top10$Description <- gsub("\\s*\\([^)]*\\)", "", top10$Term)
top10 <- top10[top10$Description %in% c('Molecular Function Activator Activity','Vasculature Development','Leukocyte Activation','Hemopoiesis','Central Nervous System Development','Angiogenesis','Calcium Ion Binding','Positive Regulation Of Programmed Cell Death','Morphogenesis Of An Epithelium','Regulation Of Angiogenesis','Maintenance Of Location','Humoral Immune Response','Developmental Maturation','Positive Regulation Of Angiogenesis','Sprouting Angiogenesis','Positive Regulation Of MAPK Cascade'),]
top10$Description <- factor(top10$Description, levels = rev(top10$Description)) 
top10$Count <- sub("/.*", "", top10$Overlap)
top10$Count <- as.character(top10$Count)
# mytheme2<- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
top10$text_x <- rep(0.03,7) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
# 2. 作图
p2 <- ggplot(data = top10,aes(x = Odds.Ratio, y = Description)) +
  geom_bar(aes(fill = P.value),stat = "identity", width = 0.8, alpha = 0.7) +        # 加边框
  scale_fill_distiller(palette = "YlOrRd", direction = -1) + # 色带倒过来
  labs(x = "Odds.Ratio", y = "Pathway", title = "GOBP enrichment barplot") +
  geom_text(aes(x = text_x, label = Description),hjust = 0,fontface = "bold",size = 7) +                   # 文本加粗
  theme_classic() +
	theme(
    axis.title.x = element_text(size = 25), # X轴标题字体大小
		axis.title.y = element_text(size = 25), # Y轴标题字体大小
		axis.text.x = element_text(size = 15), # X轴标签字体大小
		axis.text.y =  element_blank(),          # 坐标轴标题
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 图边框
    legend.text = element_text(size = 25),        # 图例文字大小
    legend.title = element_text(size = 25),       # 图例标题大小
    legend.key.size = unit(1.2, "lines"),         # 色块大小
    legend.position = c(0.95, 0.75),              # 右上角 (x,y) 坐标，0-1 比例
    legend.justification = c("right", "top")      # 对齐方式
  )

ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/10.Mission/Figure3/GOBP_barplot_select.pdf",width = 12,height = 8,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/10.Mission/Figure3/GOBP_barplot_select.png",width = 12,height = 8,bg='white')


GO_result<-read.table("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/07.Mission/Fibroblasts/clusterADC/KEGG_2021_Human.Human.enrichr.reports.txt",sep="\t",header=TRUE,fill = TRUE,quote = "",fileEncoding = "UTF-8")
# GO_result<-GO_result[1:10,]
top10<- GO_result
#2.论文画法：
top10$Description <- gsub("\\s*\\([^)]*\\)", "", top10$Term)
top10 <- top10[top10$Description %in% c('Glycolysis / Gluconeogenesis','AMPK signaling pathway','HIF-1 signaling pathway','Metabolism of xenobiotics by cytochrome P450','Ribosome','MAPK signaling pathway','Cytokine-cytokine receptor interaction'),]
top10$Description <- factor(top10$Description, levels = rev(top10$Description)) 
top10$Count <- sub("/.*", "", top10$Overlap)
top10$Count <- as.character(top10$Count)
# mytheme2<- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
top10$text_x <- rep(0.03,7) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
# 2. 作图
p2 <- ggplot(data = top10,aes(x = Odds.Ratio, y = Description)) +
  geom_bar(aes(fill = P.value),stat = "identity", width = 0.8, alpha = 0.7) +        # 加边框
  scale_fill_distiller(palette = "YlOrRd", direction = -1) + # 色带倒过来
  labs(x = "Odds.Ratio", y = "Pathway", title = "KEGG enrichment barplot") +
  geom_text(aes(x = text_x, label = Description),hjust = 0,fontface = "bold",size = 7) +                   # 文本加粗
  theme_classic() +
	theme(
    axis.title.x = element_text(size = 25), # X轴标题字体大小
		axis.title.y = element_text(size = 25), # Y轴标题字体大小
		axis.text.x = element_text(size = 15), # X轴标签字体大小
		axis.text.y =  element_blank(),          # 坐标轴标题
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 图边框
    legend.text = element_text(size = 25),        # 图例文字大小
    legend.title = element_text(size = 25),       # 图例标题大小
    legend.key.size = unit(1.2, "lines"),         # 色块大小
    legend.position = c(0.95, 0.75),              # 右上角 (x,y) 坐标，0-1 比例
    legend.justification = c("right", "top")      # 对齐方式
  )

ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/10.Mission/Figure3/KEGG_barplot_select.pdf",width = 12,height = 8,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/10.Mission/Figure3/KEGG_barplot_select.png",width = 12,height = 8,bg='white')