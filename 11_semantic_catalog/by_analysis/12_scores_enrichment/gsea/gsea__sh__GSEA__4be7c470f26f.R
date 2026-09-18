library(enrichplot)
library(clusterProfiler)
library(DOSE)
# library(ReactomePA)
library(data.table)
library(ggplot2)
library(org.Mm.eg.db)

# 读取GMT文件（使用更稳健的方法）
gmt_file <- "/mnt/NexData_file/home/sh/Projects/FWAE20242338/02.Mission/script/h.all.v2024.1.Hs.symbols.gmt"
gene_sets <- read.gmt(gmt_file)  # 使用clusterProfiler的内置函数读取

outdir <- "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/22.Mission/GSEA/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

diff_files <- c("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/21.Mission/Monocyte_LAHS-NLAHS_diffgenes.csv")
# "/mnt/NexData_file/home/sh/Projects/FWSC2024-192/21.Mission/Monocyte_LAHS-NLAHS_diffgenes_significant.csv"
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
    minGSSize = 3,
    maxGSSize = 500,
    pvalueCutoff = 1,
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
                   geneSetID = 'HALLMARK_GLYCOLYSIS', 
                   title = 'HALLMARK_GLYCOLYSIS',
                   pvalue_table = TRUE,
                   base_size = 11)
    
    # 保存图片
    plot_file <- file.path(outdir, paste0(comp_name, "_gsea_HALLMARK_GLYCOLYSIS.pdf"))
    ggsave(plot_file, p, width = 10, height = 6)
		plot_file <- file.path(outdir, paste0(comp_name, "_gsea_HALLMARK_GLYCOLYSIS.png"))
    ggsave(plot_file, p, width = 10, height = 6)
  # } else {
    # warning(paste("Target gene set not found:", target_set))
  # }
	p <- dotplot(gsea_res,  #GO富集分析结果
        x = "NES",  #横坐标,默认GeneRation,也可以为Count
        color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
        showCategory = 10,  #展示前20个点，默认为10个
        size = NULL,  #点的大小
        title = "GSEA_dotplot",  #设置图片的标题
				label_format = 100)
p + guides(colour = "NES")
ggsave(paste0(outdir,comp_name, "_GSEA_dotplot.pdf"),width = 10,height = 5)
ggsave(paste0(outdir,comp_name, "_GSEA_dotplot.png"),width = 10,height = 5)
}















