library(ggplot2)
library(stringr)
library(clusterProfiler)

outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC20252616/04.Mission/Enrich_plot/"
A <- read.csv("/mnt/NexData_file/home/sh/Projects/FWSC20252616/04.Mission/script/Myeloid_HHD-CTRL_enrich.KEGG_select.csv",sep="\t",header = T)
cluster_diff <- c("Myeloid")
diffgroup <- c("HHD/CTRL")
diffgroup <- strsplit(diffgroup, "/")

n=1
m=1
cluster_diff_sub<-cluster_diff[n]
print(cluster_diff_sub)
diff1<-diffgroup[[m]][1]
diff2<-diffgroup[[m]][2]
print(paste0(diff1,"vs",diff2))
data = read.table(paste0("/mnt/NexData_file/home/sh/Projects/FWSC20252616/04.Mission/HHD-CTRL-Myeloid/",cluster_diff_sub,"_",diff1,"-",diff2, "_diffgenes_significant.csv"),sep = ",",head=TRUE)
colnames(data) <- c("gene_name", "PValue", "log2FC", "pct_nz_group" , "pct_nz_reference", "FDR")
#自行选择筛选自己感兴趣的基因集
# gene_set <- rownames(pro_result[which(pro_result$FDR < 0.01 & abs(pro_result$log2FC) >= 2),])
gene_set <- data$gene_name
gene_symbol <- bitr(geneID = gene_set,  #感兴趣的基因集
                    fromType="SYMBOL",   #输入ID的类型
                    toType=c("ENTREZID"),   #输出ID的类型，可为多个
                    OrgDb="org.Hs.eg.db")  #物种注释数据库
gene_KEGG <- gene_symbol[,2]

KEGG<- enrichKEGG(gene = gene_KEGG,   #基因列表(同GO) 
                  organism = "hsa",  #物种
                  keyType = "kegg",  #指定的基因ID类型，默认为kegg
                  minGSSize = 3, 
                  maxGSSize = 500,
                  pvalueCutoff = 0.05,  
                  pAdjustMethod = "BH",
                  qvalueCutoff = 1)
write.table(KEGG, paste0(outdir,cluster_diff_sub,"_",diff1,"-",diff2, "_enrich.KEGG.csv"), sep = ',', row.names = FALSE, quote = TRUE)

KEGG@result$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", KEGG@result$Description)

showCategory<-A$Description
p <- dotplot(KEGG,  #GO富集分析结果
        x = "GeneRatio",  #横坐标,默认GeneRation,也可以为Count
        color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
        showCategory = showCategory,  #展示前20个点，默认为10个
        size = NULL,  #点的大小
        title = "KEGG_dotplot",  #设置图片的标题
				label_format = 100)
p + guides(colour = "Count")
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC20252616/04.Mission/Enrich_plot/Myeloid-Enrich_dotplot.pdf",width = 12,height = 5)
ggsave("/mnt/NexData_file/home/sh/Projects/FWSC20252616/04.Mission/Enrich_plot/Myeloid-Enrich_dotplot.png",width = 12,height = 5)

p<-barplot(KEGG,  #GO富集分析结果
        x = "Count",  #横坐标,默认Count,也可以为GeneRation
        color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
        showCategory = showCategory,  #展示前20个，默认为10个
        size = NULL,  
        title = "KEGG_barplot",  #设置图片的标题
				label_format = 80,orderCategory= FALSE)

ggsave(paste0(outdir,"Myeloid-Enrich_barplot.pdf"),width = 12,height = 5)
ggsave(paste0(outdir,"Myeloid-Enrich_barplot.png"),width = 12,height = 5)


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

pdf("/mnt/NexData_file/home/sh/Projects/FWSC20252616/04.Mission/Enrich_plot/M1-Enrich_barplot_order.pdf", height = 9, width = 12)
print(p)
dev.off()
png("/mnt/NexData_file/home/sh/Projects/FWSC20252616/04.Mission/Enrich_plot/M1-Enrich_barplot_order.png", height = 630, width = 840)
print(p)
dev.off()