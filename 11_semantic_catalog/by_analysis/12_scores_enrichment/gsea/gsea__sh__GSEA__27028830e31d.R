library(enrichplot)
library(clusterProfiler)
library(DOSE)
# library(ReactomePA)
library(data.table)
library(ggplot2)
library(org.Mm.eg.db)

# 读取GMT文件（使用更稳健的方法）
gmt_file <- "/mnt/NexData_file/home/sh/Projects/SC2512227671-05/02.Mission/script/lipid_pathway.gmt"
gene_sets <- read.gmt(gmt_file)  # 使用clusterProfiler的内置函数读取

outdir <- "/mnt/NexData_file/home/sh/Projects/SC2512227671-05/02.Mission/GSEA/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cluster_diff <- c("PT1_Slc5a12","PT2_Cmss1","PT3_Gpc6","PT4_Kap","PT5_Anks1b","PT6_Gramd1b","PT7_Atp11a")
diffgroup <- c("K/CF","KFD/CF","KFG/CF","KFZ/CF","KFD/K","KFG/K","KFZ/K")
diffgroup <- strsplit(diffgroup, "/")

for(n in 1:length(cluster_diff)){
	cluster_diff_sub<-cluster_diff[[n]]
	print(cluster_diff_sub)
	for(m in 1:length(diffgroup)){
		diff1<-diffgroup[[m]][1]
		diff2<-diffgroup[[m]][2]
		print(paste0(diff1,"vs",diff2))

degs = read.table(paste0("/mnt/NexData_file/home/sh/Projects/SC2512227671-05/01.Mission/diffgene_cell_type/",cluster_diff_sub,"_",diff1,"-",diff2, "_diffgenes.csv"),sep = ",",head=TRUE)

  comp_name <- paste0(cluster_diff_sub,"_",diff1,"-",diff2)
  cat("\nProcessing:", comp_name, "\n")
  
  # 读取差异表达基因
  # degs <- read.table(file, sep = ',',header=TRUE)
  
  # 准备排序的基因列表
  gene_list <- degs$avg_log2FC
  names(gene_list) <- degs$gene
  gene_list <- na.omit(gene_list)  # 移除NA值
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # 运行GSEA（关键修正：使用TERM2GENE参数）
	set.seed(123)
  gsea_res <- GSEA(
    geneList = gene_list,
    TERM2GENE = gene_sets,  # 正确参数传递方式
    minGSSize = 3,
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,  # 开启详细输出便于调试
    eps = 0         # 避免p值=0的情况
  )
  
  # 保存完整GSEA结果
  result_file <- file.path(outdir, paste0(comp_name, "_gsea_results.csv"))
  write.csv(gsea_res@result, result_file, row.names = FALSE)
  
  # 检查目标基因集是否存在
  # target_set <- "Oxidative phosphorylation"
  # if (target_set %in% gsea_res$ID) {
    # 绘制GSEA图
	geneSet<-gsea_res@result$ID[which(gsea_res@result$pvalue<=0.05)]
	if(length(geneSet)>0){
    p <- gseaplot2(gsea_res, 
                   geneSetID = geneSet, 
                   title = 'lipid geneset',
                   pvalue_table = TRUE,
                   base_size = 11)
    
    # 保存图片
    plot_file <- file.path(outdir, paste0(comp_name, "_gsea_geneset.pdf"))
    ggsave(plot_file, p, width = 10, height = 6)
    plot_file <- file.path(outdir, paste0(comp_name, "_gsea_geneset.png"))
    ggsave(plot_file, p, width = 10, height = 6)
	}
  # } else {
    # warning(paste("Target gene set not found:", target_set))
  # }
	}
}














