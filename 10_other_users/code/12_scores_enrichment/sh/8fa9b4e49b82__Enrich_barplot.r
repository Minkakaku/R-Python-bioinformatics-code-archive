library(ggplot2)
library(stringr)

A <- read.csv("/mnt/NexData_file/home/sh/Projects/FWSC20250829/06.Mission/case-control/Non-classical Monocytes_case-control_enrich.go.csv",sep=",",header = T)
A <- A[1:10,]
p <- ggplot(A,aes(Count,Description)) + 
	geom_bar(aes(y=reorder(Description,Count),x=Count,fill=p.adjust),stat='identity') + 
	scale_fill_gradient(low="#CE3D32FF",high="#9467BD") + 
	theme_bw() + 
	theme(panel.grid = element_blank(),
		panel.border = element_rect(color = "black", size = 2, fill = NA),
		axis.ticks.y = element_blank(),
		plot.title = element_text(hjust = 0.5, size = 10),
		axis.title.x = element_text(size = 20), # X轴标题字体大小
		axis.title.y = element_text(size = 20), # Y轴标题字体大小
		axis.text.x = element_text(size = 15,angle=90,hjust = 1,vjust=0.5), # X轴标签字体大小
		axis.text.y = element_text(size = 15), # Y轴标签字体大小
		plot.margin=unit(x=c(top.mar=0.2,right.mar=0.2,bottom.mar=0.2,left.mar=0.2), units="inches"),
		legend.key.size = unit(0.5, 'inches'),
		legend.title = element_text(size=20),
		legend.text = element_text(size=15)) + 
		scale_y_discrete(labels=function(x) str_wrap(x, width=35))

pdf("/mnt/NexData_file/home/sh/Projects/FWSC20250829/06.Mission/case-control/Enrich_barplot.pdf", height = 9, width = 12)
print(p)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC20250829/06.Mission/case-control/Enrich_barplot.png", height = 630, width = 840)
print(p)
dev.off()

A <- read.csv("/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/Myeloid/cluster3/GO_Biological_Process_2023.Human.enrichr.reports.txt",sep="\t",header = T)
A$Term <- sub("\\(.*", "", A$Term)
A <- A[order(-A$Odds.Ratio), ][1:10, c(2,5,8)]
colnames(A)<-c('Description','p.adjust','Odds.Ratio')
p <- ggplot(A,aes(Odds.Ratio,Description)) + 
	geom_bar(aes(y=reorder(Description,Odds.Ratio),x=Odds.Ratio,fill=p.adjust),stat='identity') + 
	scale_fill_gradient(low="#CE3D32FF",high="#9467BD") + 
	theme_bw() + 
	theme(panel.grid = element_blank(),
		panel.border = element_rect(color = "black", size = 2, fill = NA),
		axis.ticks.y = element_blank(),
		plot.title = element_text(hjust = 0.5, size = 10),
		axis.title.x = element_text(size = 20), # X轴标题字体大小
		axis.title.y = element_text(size = 20), # Y轴标题字体大小
		axis.text.x = element_text(size = 15,angle=90,hjust = 1,vjust=0.5), # X轴标签字体大小
		axis.text.y = element_text(size = 15), # Y轴标签字体大小
		plot.margin=unit(x=c(top.mar=0.2,right.mar=0.2,bottom.mar=0.2,left.mar=0.2), units="inches"),
		legend.key.size = unit(0.5, 'inches'),
		legend.title = element_text(size=20),
		legend.text = element_text(size=15)) + 
		scale_y_discrete(labels=function(x) str_wrap(x, width=35))

pdf("/mnt/NexData_file/home/sh/Projects/FWSC20250829/06.Mission/case-control/cluster-Enrich_barplot.pdf", height = 9, width = 12)
print(p)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC20250829/06.Mission/case-control/cluster-Enrich_barplot.png", height = 630, width = 840)
print(p)
dev.off()