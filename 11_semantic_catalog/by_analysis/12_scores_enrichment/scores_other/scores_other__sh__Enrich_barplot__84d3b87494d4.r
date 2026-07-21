library(ggplot2)
library(stringr)

A <- read.csv("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/script/genelist_enrich.go_select.csv",sep="\t",header = T)
p <- ggplot(A,aes(Count,Description)) + 
	geom_bar(aes(y=reorder(Description,Count),x=Count,fill=p.adjust),stat='identity') + 
	scale_fill_gradient(low="#CE3D32FF",high="#d9950b") + 
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

pdf("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/Enrich_plot/Enrich_barplot.pdf", height = 9, width = 12)
print(p)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/Enrich_plot/Enrich_barplot.png", height = 630, width = 840)
print(p)
dev.off()