# conda activate diffgene
#加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
# gene_symbol <- bitr(geneID = gene_set,  #感兴趣的基因集
#                     fromType="ENSEMBL",   #输入ID的类型
#                     toType=c("SYMBOL", "ENTREZID"),   #输出ID的类型，可为多个
#                     OrgDb="org.Hs.eg.db")  #物种注释数据库
# gene <- gene_symbol[,3]
# setwd("/data/NAS02/sh/Projects/FWSC20241287/02.Mission/Enrich/")
outdir <- paste0("G:/FWSC20242121-05/售后1/DEG_FunctionEnrichment")
dir.create(outdir)
cluster_diff <- list('Endothelial cells',
                     'SMC',
                     'Maternal_M1',
                     'Maternal_M2',
                     'Trophoblasts')
diffgroup <- c("HC/PTB")
diffgroup <- strsplit(diffgroup, "/")

for(n in 1:length(cluster_diff)){
	cluster_diff_sub<-cluster_diff[[n]]
	print(cluster_diff_sub)
	for(m in 1:length(diffgroup)){
		diff1<-diffgroup[[m]][1]
		diff2<-diffgroup[[m]][2]
		print(paste0(diff1,"vs",diff2))
data = read.table(paste0("G:/FWSC20242121-05/售后1/DEG_FunctionEnrichment/",cluster_diff_sub,"/",diff1,"-",diff2, "_diffgenes_significant.csv"),sep = ",",head=TRUE)
colnames(data) <- c("gene_name", "PValue", "log2FC", "pct_nz_group" , "pct_nz_reference", "FDR")
#自行选择筛选自己感兴趣的基因集
# gene_set <- rownames(pro_result[which(pro_result$FDR < 0.01 & abs(pro_result$log2FC) >= 2),])
gene_set <- data$gene_name

GO_ALL <- enrichGO(gene = gene_set,  #基因列表(转换的ID)
               keyType = "SYMBOL",  #指定的基因ID类型，默认为ENTREZID
               OrgDb=org.Hs.eg.db,  #物种对应的org包
               ont = "ALL",   #CC细胞组件，MF分子功能，BP生物学过程
               pvalueCutoff = 0.05,  #p值阈值
               pAdjustMethod = "BH",  #多重假设检验校正方式
               minGSSize = 3,   #注释的最小基因集，默认为10
               maxGSSize = 500,  #注释的最大基因集，默认为500
               qvalueCutoff = 1,  #q值阈值
               readable = TRUE)  #基因ID转换为基因名
write.table(GO_ALL, paste0(outdir,"/",cluster_diff_sub,"/",diff1,"-",diff2, "_enrich.go.csv"), sep = ',', row.names = FALSE, quote = TRUE)

if(nrow(GO_ALL)>0){
dotplot(GO_ALL,  #GO富集分析结果
        x = "GeneRatio",  #横坐标,默认GeneRation,也可以为Count
        color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
        showCategory = 10,  #展示前20个点，默认为10个
        size = NULL,  #点的大小
        title = "GO_dotplot",  #设置图片的标题
        split='ONTOLOGY',
				label_format = 80)+
				facet_grid(ONTOLOGY~.,scale="free")

ggsave(paste0(outdir,"/",cluster_diff_sub,"/",diff1,"-",diff2, ".GO_dotplot.pdf"),width = 12,height = 9)
ggsave(paste0(outdir,"/",cluster_diff_sub,"/",diff1,"-",diff2, ".GO_dotplot.png"),width = 12,height = 9)

barplot(GO_ALL,  #GO富集分析结果
        x = "Count",  #横坐标,默认Count,也可以为GeneRation
        color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
        showCategory = 10,  #展示前20个，默认为10个
        size = NULL,  
        title = "GO_barplot",  #设置图片的标题
        split='ONTOLOGY',
				label_format = 80)+
				facet_grid(ONTOLOGY~.,scale="free")

ggsave(paste0(outdir,"/",cluster_diff_sub,"/",diff1,"-",diff2, ".GO_barplot.pdf"),width = 12,height = 9)
ggsave(paste0(outdir,"/",cluster_diff_sub,"/",diff1,"-",diff2, ".GO_barplot.png"),width = 12,height = 9)
}

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
write.table(KEGG, paste0(outdir,"/",cluster_diff_sub,"/",diff1,"-",diff2, "_enrich.KEGG.csv"), sep = ',', row.names = FALSE, quote = TRUE)

KEGG@result$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", KEGG@result$Description)
if(nrow(KEGG)>0){
p <- dotplot(KEGG,  #GO富集分析结果
        x = "GeneRatio",  #横坐标,默认GeneRation,也可以为Count
        color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
        showCategory = 10,  #展示前20个点，默认为10个
        size = NULL,  #点的大小
        title = "KEGG_dotplot",  #设置图片的标题
				label_format = 100)
p + guides(colour = "Count")
ggsave(paste0(outdir,"/",cluster_diff_sub,"/",diff1,"-",diff2, ".KEGG_dotplot.pdf"),width = 12,height = 5)
ggsave(paste0(outdir,"/",cluster_diff_sub,"/",diff1,"-",diff2, ".KEGG_dotplot.png"),width = 12,height = 5)

barplot(KEGG,  #GO富集分析结果
        x = "Count",  #横坐标,默认Count,也可以为GeneRation
        color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
        showCategory = 10,  #展示前20个，默认为10个
        size = NULL,  
        title = "KEGG_barplot",  #设置图片的标题
				label_format = 80)

ggsave(paste0(outdir,"/",cluster_diff_sub,"/",diff1,"-",diff2, ".KEGG_barplot.pdf"),width = 12,height = 5)
ggsave(paste0(outdir,"/",cluster_diff_sub,"/",diff1,"-",diff2, ".KEGG_barplot.png"),width = 12,height = 5)
		}
	}
}