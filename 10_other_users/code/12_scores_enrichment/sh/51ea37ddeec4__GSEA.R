library(enrichplot)
library(clusterProfiler)
library(DOSE)
# library(ReactomePA)
library(data.table)
library(ggplot2)
library(org.Mm.eg.db)

# 读取GMT文件（使用更稳健的方法）
gmt_file <- "/mnt/NexData_file/home/sh/Projects/SC2507163761-01/05.Mission/script/genset.gmt"
gene_sets <- read.gmt(gmt_file)  # 使用clusterProfiler的内置函数读取

outdir <- "/mnt/NexData_file/home/sh/Projects/SC2507163761-01/05.Mission/GSEA/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cluster_diff <- c("M1","M2a","M2b","M2c","M2d")
for(n in 1:length(cluster_diff)){
	comp_name<-cluster_diff[[n]]
	print(comp_name)
	degs <- read.table(paste0("/mnt/NexData_file/home/sh/Projects/SC2507163761-01/05.Mission/diffgene_cell_type/",comp_name,"_UDA-CMC_diffgenes.csv"),sep = ",",head=TRUE)
	# degs <- read.table("/mnt/NexData_file/home/sh/Projects/SC2507163761-01/04.Mission/diffgene-cell_type/Proximal tubule cells_UDA-CMC_diffgenes.csv",sep = ",",head=TRUE)
# for (file in diff_files) {
  # comp_name <- gsub("_diffgenes.csv", "", basename(file))
  # cat("\nProcessing:", comp_name, "\n")
  
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
    pvalueCutoff = 0.05,
    verbose = TRUE,  # 开启详细输出便于调试
    eps = 0,         # 避免p值=0的情况
    seed = 123       # 设置随机种子保证可重复性
  )
  
  # 保存完整GSEA结果
  result_file <- file.path(outdir, paste0(comp_name, "_gsea_results.csv"))
  write.csv(gsea_res@result, result_file, row.names = FALSE)
  tmp <- as.data.frame(gsea_res)
  # 检查目标基因集是否存在
  # target_set <- "Oxidative phosphorylation"
  # if (target_set %in% gsea_res$ID) {
    # 绘制GSEA图
		for(i in gsea_res@result$ID){
			p <- gseaplot2(gsea_res, 
                   geneSetID = i, 
                   title = i,
                   pvalue_table = TRUE,
                   base_size = 11)
    
    # 保存图片
    plot_file <- file.path(outdir, paste0(comp_name, "_",i,".pdf"))
		ggsave(plot_file, p, width = 10, height = 6)
		plot_file <- file.path(outdir, paste0(comp_name, "_",i,".png"))
    ggsave(plot_file, p, width = 10, height = 6)
		}

  # } else {
    # warning(paste("Target gene set not found:", target_set))
  # }
# }
}