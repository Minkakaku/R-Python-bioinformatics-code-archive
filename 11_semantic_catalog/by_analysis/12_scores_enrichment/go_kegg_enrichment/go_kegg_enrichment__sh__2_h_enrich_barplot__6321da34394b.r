library(dplyr)
library(ggplot2)
library(RColorBrewer)

#载入差异基因富集结果表：
#这里使用自己的富集结果表即可
# load( "GO_diff.Rdata")
# View(GO_result) #clusterProfiler富集结果表
# GO_result<- mutate(GO_result, RichFactor= Count / as.numeric(sub( "/\\d+", "", BgRatio)))

GO_result<-read.table("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/07.Mission/Epithelial/clusterADC/GO_Biological_Process_2023.Human.enrichr.reports.txt",sep="\t",header=TRUE,fill = TRUE,quote = "",fileEncoding = "UTF-8")
# GO_result<-GO_result[1:10,]
top10<- GO_result
#2.论文画法：
top10$Description <- gsub("\\s*\\([^)]*\\)", "", top10$Term)
top10 <- top10[top10$Description %in% c('Actin Filament Organization','Establishment Of Protein Localization To Extracellular Region','Organic Acid Transport','Extracellular Structure Organization','Extracellular Matrix Organization','External Encapsulating Structure Organization','Amino Acid Metabolic Process','Peptide Secretion','Peptide Hormone Secretion','Cellular Ketone Metabolic Process','Alpha−Amino Acid Metabolic Process','Carbohydrate Transport','Biogenic Amine Metabolic Process','Dicarboxylic Acid Metabolic Process','Extracellular Matrix Disassembly'),]
top10$Description <- factor(top10$Description, levels = rev(top10$Description)) 
top10$Count <- sub("/.*", "", top10$Overlap)
top10$Count <- as.character(top10$Count)
# mytheme2<- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
top10$text_x <- rep(0.03,8) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
# 2. 作图
p2 <- ggplot(data = top10,aes(x = Odds.Ratio, y = Description)) +
  geom_bar(aes(fill = Adjusted.P.value),stat = "identity", width = 0.8, alpha = 0.7) +        # 加边框
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

ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/GOBP_barplot_select.pdf",width = 15,height = 12,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/GOBP_barplot_select.png",width = 15,height = 12,bg='white')


top10<- GO_result[1:10,]
#2.论文画法：
top10$Description <- gsub("\\s*\\([^)]*\\)", "", top10$Term)
top10$Description <- factor(top10$Description, levels = rev(top10$Description)) 
top10$Count <- sub("/.*", "", top10$Overlap)
top10$Count <- as.character(top10$Count)
# mytheme2<- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
top10$text_x <- rep(0.03,10) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
# 2. 作图
p2 <- ggplot(data = top10,aes(x = Odds.Ratio, y = Description)) +
  geom_bar(aes(fill = Adjusted.P.value),stat = "identity", width = 0.8, alpha = 0.7) +        # 加边框
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

ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/GOBP_barplot_top10.pdf",width = 15,height = 12,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/GOBP_barplot_top10.png",width = 15,height = 12,bg='white')


GO_result<-read.table("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/07.Mission/Epithelial/clusterADC/GO_Cellular_Component_2023.Human.enrichr.reports.txt",sep="\t",header=TRUE,fill = TRUE,quote = "",fileEncoding = "UTF-8")
# GO_result<-GO_result[1:10,]
top10<- GO_result
#2.论文画法：
top10$Description <- gsub("\\s*\\([^)]*\\)", "", top10$Term)
top10 <- top10[top10$Description %in% c('Actin Cytoskeleton','Collagen−Containing Extracellular Matrix','Lytic Vacuole Membrane','Lysosomal Membrane','Cell Projection Membrane','Motile Cilium','Basal Part Of Cell','Basolateral Plasma Membrane','Ciliary Plasm','Axoneme','Basement Membrane','Actin Filament Bundle','Peroxisomal Matrix','Microbody Lumen','M Band'),]
top10$Description <- factor(top10$Description, levels = rev(top10$Description)) 
top10$Count <- sub("/.*", "", top10$Overlap)
top10$Count <- as.character(top10$Count)
# mytheme2<- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
top10$text_x <- rep(0.03,7) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
# 2. 作图
p2 <- ggplot(data = top10,aes(x = Odds.Ratio, y = Description)) +
  geom_bar(aes(fill = Adjusted.P.value),stat = "identity", width = 0.8, alpha = 0.7) +        # 加边框
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

ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/GOCC_barplot_select.pdf",width = 15,height = 12,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/GOCC_barplot_select.png",width = 15,height = 12,bg='white')


top10<- GO_result[1:10,]
#2.论文画法：
top10$Description <- gsub("\\s*\\([^)]*\\)", "", top10$Term)
top10$Description <- factor(top10$Description, levels = rev(top10$Description)) 
top10$Count <- sub("/.*", "", top10$Overlap)
top10$Count <- as.character(top10$Count)
# mytheme2<- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
top10$text_x <- rep(0.03,10) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
# 2. 作图
p2 <- ggplot(data = top10,aes(x = Odds.Ratio, y = Description)) +
  geom_bar(aes(fill = Adjusted.P.value),stat = "identity", width = 0.8, alpha = 0.7) +        # 加边框
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

ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/GOCC_barplot_top10.pdf",width = 15,height = 12,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure2/GOCC_barplot_top10.png",width = 15,height = 12,bg='white')



GO_result<-read.table("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/07.Mission/Fibroblasts/clusterADC/GO_Biological_Process_2023.Human.enrichr.reports.txt",sep="\t",header=TRUE,fill = TRUE,quote = "",fileEncoding = "UTF-8")
# GO_result<-GO_result[1:10,]
top10<- GO_result
#2.论文画法：
top10$Description <- gsub("\\s*\\([^)]*\\)", "", top10$Term)
top10 <- top10[top10$Description %in% c('Molecular Function Activator Activity','Vasculature Development','Leukocyte Activation','Hemopoiesis','Central Nervous System Development','Angiogenesis','Calcium Ion Binding','Positive Regulation Of Programmed Cell Death','Morphogenesis Of An Epithelium','Regulation Of Angiogenesis','Maintenance Of Location','Humoral Immune Response','Developmental Maturation','Positive Regulation Of Angiogenesis','Sprouting Angiogenesis'),]
top10$Description <- factor(top10$Description, levels = rev(top10$Description)) 
top10$Count <- sub("/.*", "", top10$Overlap)
top10$Count <- as.character(top10$Count)
# mytheme2<- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
top10$text_x <- rep(0.03,6) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
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

ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure3/GOBP_barplot_select.pdf",width = 12,height = 8,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure3/GOBP_barplot_select.png",width = 12,height = 8,bg='white')


top10<- GO_result[1:10,]
#2.论文画法：
top10$Description <- gsub("\\s*\\([^)]*\\)", "", top10$Term)
top10$Description <- factor(top10$Description, levels = rev(top10$Description)) 
top10$Count <- sub("/.*", "", top10$Overlap)
top10$Count <- as.character(top10$Count)
# mytheme2<- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
top10$text_x <- rep(0.03,10) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
# 2. 作图
p2 <- ggplot(data = top10,aes(x = Odds.Ratio, y = Description)) +
  geom_bar(aes(fill = Adjusted.P.value),stat = "identity", width = 0.8, alpha = 0.7) +        # 加边框
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

ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure3/GOBP_barplot_top10.pdf",width = 15,height = 12,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure3/GOBP_barplot_top10.png",width = 15,height = 12,bg='white')


GO_result<-read.table("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure4/enrich/AAEs_Fibroblasts_ADC-SCC_enrich.go.csv",sep=",",header=TRUE)
# GO_result<-GO_result[1:10,]
top10<- GO_result
#2.论文画法：
# top10$Description <- gsub("\\s*\\([^)]*\\)", "", top10$Term)
top10 <- top10[top10$Description %in% c('vasculature development','signaling receptor activator activity','peptidyl−tyrosine phosphorylation','peptidyl−tyrosine modification','angiogenesis','positive regulation of MAPK cascade','gland development','response to growth factor','regulation of cell development','head development','immune response','Golgi apparatus','regulation of vasculature development','tube formation','sprouting angiogenesis','presynaptic membrane','granulocyte migration','cell recognition','regulation of viral process','adaptive thermogenesis'),]
top10$Description <- factor(top10$Description, levels = rev(top10$Description)) 
# top10$Count <- sub("/.*", "", top10$Overlap)
# top10$Count <- as.character(top10$Count)
# mytheme2<- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
top10$text_x <- rep(0.03,5) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
# 2. 作图
p2 <- ggplot(data = top10,aes(x = Count, y = Description)) +
  geom_bar(aes(fill = p.adjust),stat = "identity", width = 0.8, alpha = 0.7) +        # 加边框
  scale_fill_distiller(palette = "YlOrRd", direction = -1) + # 色带倒过来
  labs(x = "Count", y = "Pathway", title = "GOBP enrichment barplot") +
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

ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure4/enrich/GOBP_barplot_select.pdf",width = 15,height = 10,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure4/enrich/GOBP_barplot_select.png",width = 15,height = 10,bg='white')


top10<- GO_result[1:10,]
#2.论文画法：
# top10$Description <- gsub("\\s*\\([^)]*\\)", "", top10$Term)
top10$Description <- factor(top10$Description, levels = rev(top10$Description)) 
# top10$Count <- sub("/.*", "", top10$Overlap)
# top10$Count <- as.character(top10$Count)
# mytheme2<- mytheme + theme(axis.text.y = element_blank()) #先在自定义主题中隐去y轴文本标签显示
top10$text_x <- rep(0.03,10) #新增一列重复数组，使绘图时文本标签能从固定位置开始；
# 2. 作图
p2 <- ggplot(data = top10,aes(x = Count, y = Description)) +
  geom_bar(aes(fill = p.adjust),stat = "identity", width = 0.8, alpha = 0.7) +        # 加边框
  scale_fill_distiller(palette = "YlOrRd", direction = -1) + # 色带倒过来
  labs(x = "Count", y = "Pathway", title = "GOBP enrichment barplot") +
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
ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure4/enrich/GOBP_barplot_top10.pdf",width = 15,height = 12,bg='white')
ggsave("/mnt/NexData_file/home/sh/Projects/AE2507033447-07/08.Mission/Figure4/enrich/GOBP_barplot_top10.png",width = 15,height = 12,bg='white')