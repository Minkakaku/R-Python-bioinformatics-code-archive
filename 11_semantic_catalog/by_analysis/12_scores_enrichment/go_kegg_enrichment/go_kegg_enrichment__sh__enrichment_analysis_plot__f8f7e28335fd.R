library(ggplot2)
library(stringr)
#GO可视化之BP
BP <- read.delim("GO可视化文件之BP.txt")
#柱状图
#tiff("Up_barplot.tiff", height = 15, width = 30, units = "cm", res = 600)
pdf('GO_barplot之BP.pdf',width=8,height=5)
p <- ggplot(data = BP)+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=-log10(pvalue)),stat='identity')+
  scale_fill_gradient(expression(-log["10"](pvalue)),low="blue",high="red")+
  ylab("")+
  xlab("Gene Count")+ # 设置坐标轴名称
  #facet_grid(ONTOLOGY~., scales ="free")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 12,color = "black", face = "bold", vjust = 0.5))+
  scale_y_discrete(labels=function(x)str_wrap(x,width = 55))+
  theme(axis.title.x = element_text(size = 11,color = "black", face = "bold", vjust = 0.5))+
  theme(legend.title = element_text(face = 'bold'))
p
dev.off()

#GO可视化之CC
CC <- read.delim("GO可视化文件之CC.txt")
#柱状图
#tiff("Up_barplot.tiff", height = 15, width = 30, units = "cm", res = 600)
pdf('GO_barplot之CC.pdf',width=8,height=5)
p <- ggplot(data = CC)+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=-log10(pvalue)),stat='identity')+
  scale_fill_gradient(expression(-log["10"](pvalue)),low="blue",high="red")+
  ylab("")+
  xlab("Gene Count")+ # 设置坐标轴名称
  #facet_grid(ONTOLOGY~., scales ="free")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 12,color = "black", face = "bold", vjust = 0.5))+
  scale_y_discrete(labels=function(x)str_wrap(x,width = 55))+
  theme(axis.title.x = element_text(size = 11,color = "black", face = "bold", vjust = 0.5))+
  theme(legend.title = element_text(face = 'bold'))
p
dev.off()

#GO可视化之MF
MF <- read.delim("GO可视化文件之MF.txt")
#柱状图
#tiff("Up_barplot.tiff", height = 15, width = 30, units = "cm", res = 600)
pdf('GO_barplot之MF.pdf',width=8,height=5)
p <- ggplot(data = MF)+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=-log10(pvalue)),stat='identity')+
  scale_fill_gradient(expression(-log["10"](pvalue)),low="blue",high="red")+
  ylab("")+
  xlab("Gene Count")+ # 设置坐标轴名称
  #facet_grid(ONTOLOGY~., scales ="free")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 12,color = "black", face = "bold", vjust = 0.5))+
  scale_y_discrete(labels=function(x)str_wrap(x,width = 55))+
  theme(axis.title.x = element_text(size = 11,color = "black", face = "bold", vjust = 0.5))+
  theme(legend.title = element_text(face = 'bold'))
p
dev.off()

#KEGG气泡图
kk <- read.delim("KK可视化文件.txt")
kk$Description<-factor(kk$Description,levels = rev(kk$Description))
#tiff("KK_dotplot.tiff", height = 15, width = 30, units = "cm", res = 600)
pdf('KEGG_dotplot.pdf',width=11,height=5)
p <- ggplot(kk,aes(x=reorder(Description,GeneRatio),y=GeneRatio))+
  geom_point(aes(color=-log10(pvalue),size=Count),alpha=0.8)+
  scale_color_gradient(expression(-log["10"](pvalue)),low = "blue",high = "red")+
  scale_size_continuous(range=c(4,9))+
  labs(x='',y='GeneRatio')+
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(size = 12,color = "black", face = "bold", vjust = 0.5))+
  #  scale_y_discrete(labels=function(x)str_wrap(x,width = 100))+
  theme(axis.title.x = element_text(size = 11,color = "black", face = "bold", vjust = 0.5))+
  theme(legend.title = element_text(face = 'bold'))
p
dev.off()
