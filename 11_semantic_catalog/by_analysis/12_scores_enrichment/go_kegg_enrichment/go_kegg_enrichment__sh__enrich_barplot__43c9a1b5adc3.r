library(dplyr)
library(ggplot2)
library(RColorBrewer)

#载入差异基因富集结果表：
#这里使用自己的富集结果表即可
# load( "GO_diff.Rdata")
# View(GO_result) #clusterProfiler富集结果表
# GO_result<- mutate(GO_result, RichFactor= Count / as.numeric(sub( "/\\d+", "", BgRatio)))

GO_result<-read.table("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/14.Mission/Enrich/genelist_enrich.go.csv",sep=",",header=TRUE)
GO_result<-GO_result[1:10,]
##计算Fold Enrichment（富集倍数）：

# GO_result$FoldEnrichment <- apply(GO_result,1,function(x){
# GeneRatio<- eval(parse(text = x[ "GeneRatio"]))
# BgRatio<- eval(parse(text = x[ "BgRatio"]))
# foldEnrichment<- round(GeneRatio/BgRatio,2)
# foldEnrichment
# })

# head(GO_result$RichFactor)
# head(GO_result$FoldEnrichment)
# colnames(GO_result)

#提取结果前Top10绘图(或自定义所需pathway绘图)：

# top10<- GO_result[1:10,]
top10<- GO_result
#指定绘图顺序（转换为因子）：

# top10$pathway <- factor(top10$Description,levels = rev(top10$Description))

######
#1.常规画法：
#自定义主题：
mytheme<- theme(axis.title = element_text(size = 13),
axis.text = element_text(size = 11),
plot.title = element_text(size = 14,
hjust= 0.5,
face= "bold"),
legend.title = element_text(size = 13),
legend.text = element_text(size = 11))

#Top10富集数目条形图：

p<- ggplot(data = top10,
aes(x = Count, y = pathway, fill = -log10(pvalue)))+
scale_fill_distiller(palette = "RdPu",direction = 1) +
geom_bar(stat = "identity", width = 0.8) +
labs(x = "Number of Gene", y = "pathway", title = "GO enrichment barplot") +
theme_bw()+ mytheme
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/enrich_barplot/barplot.pdf",width = 10,height = 10,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/enrich_barplot/barplot.png",width = 10,height = 10,bg='white')

#Top10富集因子条形图（挑选所需数据列绘图即可）：
# p1<- ggplot(data = top10,
# aes(x = -log10(pvalue), y = pathway, fill = RichFactor)) +
# scale_fill_distiller(palette = "RdPu",direction = 1) +
# geom_bar(stat = "identity", width = 0.8) +
# labs(x = "Number of Gene", y = "pathway", title = "GO enrichment barplot") +
# theme_bw()+ mytheme
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/enrich_barplot/barplot_top10.pdf",width = 10,height = 10,bg='white')
# ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/enrich_barplot/barplot_top10.png",width = 10,height = 10,bg='white')


#2.论文画法：
top10$Description <- factor(top10$Description, levels = rev(top10$Description)) 
# mytheme2<- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
top10$text_x <- rep(0.03,10) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
# 2. 作图
p2 <- ggplot(data = top10,aes(x = Count, y = Description)) +
  geom_bar(aes(fill = p.adjust),stat = "identity", width = 0.8, alpha = 0.7) +        # 加边框
  scale_fill_distiller(palette = "YlOrRd", direction = -1) + # 色带倒过来
  labs(x = "Count", y = "Pathway", title = "GO enrichment barplot") +
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

ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/enrich_barplot/barplot_top10.pdf",width = 16,height = 10,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/enrich_barplot/barplot_top10.png",width = 16,height = 10,bg='white')