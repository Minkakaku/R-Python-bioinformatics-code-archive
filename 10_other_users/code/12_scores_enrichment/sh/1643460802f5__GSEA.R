library(enrichplot)
library(clusterProfiler)
library(DOSE)
# library(ReactomePA)
library(data.table)
library(ggplot2)
library(org.Mm.eg.db)

# 读取GMT文件（使用更稳健的方法）
gmt_file <- "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/16.Mission/script/geneset.gmt"
gene_sets <- read.gmt(gmt_file)  # 使用clusterProfiler的内置函数读取

outdir <- "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/GSEA/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

diff_files <- c("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/Main-diffgene/Monocyte_Post-Before_diffgenes.csv","/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/Monocyte-diffgene/Classical monocyte_Post-Before_diffgenes.csv","/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/Monocyte-diffgene/Dendritic cell_Post-Before_diffgenes.csv","/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/Monocyte-diffgene/Non-classical monocyte_Post-Before_diffgenes.csv","/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/Monocyte-diffgene/PTTG1 classical monocyte_Post-Before_diffgenes.csv","/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/Monocyte-diffgene/S100A8 classical monocyte_Post-Before_diffgenes.csv","/mnt/NexData_file/home/sh/Projects/FWSC2024-192/17.Mission/Monocyte-diffgene/ZFP36L1 classical monocyte_Post-Before_diffgenes.csv")
for (file in diff_files) {
  comp_name <- gsub("_diffgenes.csv", "", basename(file))
  cat("\nProcessing:", comp_name, "\n")
  
  # 读取差异表达基因
  degs <- read.table(file, sep = ',',header=TRUE)
  
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
    minGSSize = 5,
    maxGSSize = 1000,
    pvalueCutoff = 0.9,
    verbose = TRUE,  # 开启详细输出便于调试
    eps = 0,         # 避免p值=0的情况
    seed = 123       # 设置随机种子保证可重复性
  )
  
  # 保存完整GSEA结果
  result_file <- file.path(outdir, paste0(comp_name, "_gsea_results.csv"))
  write.csv(gsea_res@result, result_file, row.names = FALSE)
  
  # 检查目标基因集是否存在
  # target_set <- "Oxidative phosphorylation"
  # if (target_set %in% gsea_res$ID) {
    # 绘制GSEA图
    p <- gseaplot2(gsea_res, 
                   geneSetID = 'geneset', 
                   title = 'geneset',
                   pvalue_table = TRUE,
                   base_size = 11)
    
    # 保存图片
    plot_file <- file.path(outdir, paste0(comp_name, "_gsea_geneset.pdf"))
    ggsave(plot_file, p, width = 10, height = 6)
  # } else {
    # warning(paste("Target gene set not found:", target_set))
  # }
}















