library(enrichplot)
library(clusterProfiler)
library(DOSE)
# library(ReactomePA)
library(data.table)
library(ggplot2)
library(org.Mm.eg.db)

# 读取GMT文件（使用更稳健的方法）
gmt_file <- "/mnt/NexData_file/home/sh/Projects/AE2507033447-08/19.Mission/script/KEGG_metabolism_nc.gmt"
gene_sets <- read.gmt(gmt_file)  # 使用clusterProfiler的内置函数读取

# outdir <- "/mnt/NexData_file/home/sh/Projects/AE2507033447-08/19.Mission/CD4_Tfh_pcr_A_GSEA/"
# outdir <- "/mnt/NexData_file/home/sh/Projects/AE2507033447-08/19.Mission/CD8_Tm_C2_pcr_A_GSEA/"
outdir <- "/mnt/NexData_file/home/sh/Projects/AE2507033447-08/19.Mission/KRT16_Tumor_GSEA/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# diff_files <- c("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD4T/groupDEG/CD4_Tfh/pcr_B_vs_pcr_A/clusterpcr_A/differentially_expressed_genes.xls")
# diff_files <- c("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD8T/groupDEG/CD8_Tm_C2/pcr_B_vs_pcr_A/clusterpcr_A/differentially_expressed_genes.xls")
diff_files <- c("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/Epithelial/clusterKRT16+ Tumor/differentially_expressed_genes.xls")
for (file in diff_files) {
  # comp_name <- gsub("_diffgenes.csv", "", basename(file))
	# comp_name <- 'CD4_Tfh_pcr_A'
	# comp_name <- 'CD8_Tm_C2_pcr_A'
	comp_name <- 'KRT16_Tumor'
  cat("\nProcessing:", comp_name, "\n")
  
  # 读取差异表达基因
  degs <- read.table(file, sep = '\t',header=TRUE)
  
  # 准备排序的基因列表
  # gene_list <- degs$avg_log2FC
  # names(gene_list) <- degs$gene
  # gene_list <- na.omit(gene_list)  # 移除NA值
  # gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- degs$logfoldchanges
  names(gene_list) <- degs$names
  gene_list <- na.omit(gene_list)  # 移除NA值
  gene_list <- sort(gene_list, decreasing = TRUE)
  # 运行GSEA（关键修正：使用TERM2GENE参数）
	set.seed(123)
  gsea_res <- GSEA(
    geneList = gene_list,
    TERM2GENE = gene_sets,  # 正确参数传递方式
    minGSSize = 1,
    maxGSSize = 1000,
    pvalueCutoff = 1,
    verbose = TRUE,  # 开启详细输出便于调试
    eps = 0,         # 避免p值=0的情况
    seed = 123       # 设置随机种子保证可重复性
  )
  
  # 保存完整GSEA结果
  result_file <- file.path(outdir, paste0(comp_name, "_gsea_results.csv"))
  write.csv(gsea_res@result, result_file, row.names = FALSE)
  
	all_sets <- gsea_res$ID
	# all_sets <- c("Angiogenesis")
# 循环处理每个基因集
for (i in seq_along(all_sets)) {
  target_set <- all_sets[i]
  
  # 绘制GSEA图
  p <- gseaplot2(gsea_res, 
                 geneSetID = target_set,
                 title = paste0(target_set, " (NES=", round(gsea_res$NES[gsea_res$ID == target_set], 2), ")"),
                 pvalue_table = TRUE,
                 base_size = 11)
  
  # 保存图片
  plot_file <- file.path(outdir, paste0(comp_name, "_gsea_", target_set, ".pdf"))
  ggsave(plot_file, p, width = 10, height = 6)
  plot_file <- file.path(outdir, paste0(comp_name, "_gsea_", target_set, ".png"))
  ggsave(plot_file, p, width = 10, height = 6)
  message(paste("Saved:", plot_file))
}
}