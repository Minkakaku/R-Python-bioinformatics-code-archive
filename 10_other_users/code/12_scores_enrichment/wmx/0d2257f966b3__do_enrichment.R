.libPaths("/mnt/NexData_file/share/scRNA-tools/envs/pip_diffgene_enrich/lib/R/library")
library(optparse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)

# ==== 参数设置 ====
option_list <- list(
  make_option(c("-f", "--file"), type = "character", help = "输入差异基因 xls 文件"),
  make_option(c("-s", "--species"), type = "character", default = "human", help = "物种: human 或 mouse"),
  make_option(c("-o", "--outdir"), type = "character", help = "输出结果目录")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ==== 路径设置 ====
input_file <- opt$file
species <- opt$species
outdir <- opt$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==== 读取数据 ====
data <- read.table(input_file, sep = "\t", header = TRUE, check.names = FALSE)
gene_set <- unique(data[, 1])  # 默认使用第一列为基因 symbol

# ==== 物种数据库选择 ====
if (species == "human") {
  OrgDb <- org.Hs.eg.db
  org <- "hsa"
} else if (species == "mouse") {
  OrgDb <- org.Mm.eg.db
  org <- "mmu"
} else {
  stop("错误：不支持的物种，仅支持 'human' 或 'mouse'")
}

# ==== GO 富集分析 ====
GO_ALL <- enrichGO(
  gene = gene_set,
  keyType = "SYMBOL",
  OrgDb = OrgDb,
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 3,
  maxGSSize = 500,
  qvalueCutoff = 1,
  readable = TRUE
)
write.csv(GO_ALL@result, file = file.path(outdir, "GO_enrich_result.csv"), row.names = FALSE)

if (nrow(GO_ALL) > 0) {
  p1 <- dotplot(GO_ALL, x = "GeneRatio", color = "p.adjust", showCategory = 10, split = "ONTOLOGY", label_format = 80) +
    facet_grid(ONTOLOGY ~ ., scale = "free")
  ggsave(file.path(outdir, "GO_dotplot.pdf"), p1, width = 12, height = 9)
  ggsave(file.path(outdir, "GO_dotplot.png"), p1, width = 12, height = 9)

  p2 <- barplot(GO_ALL, x = "Count", color = "p.adjust", showCategory = 10, split = "ONTOLOGY", label_format = 80) +
    facet_grid(ONTOLOGY ~ ., scale = "free")
  ggsave(file.path(outdir, "GO_barplot.pdf"), p2, width = 12, height = 9)
  ggsave(file.path(outdir, "GO_barplot.png"), p2, width = 12, height = 9)
}

# ==== KEGG 富集分析 ====
gene_symbol <- bitr(geneID = gene_set, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
gene_KEGG <- unique(gene_symbol$ENTREZID)

KEGG <- enrichKEGG(
  gene = gene_KEGG,
  organism = org,
  keyType = "kegg",
  minGSSize = 3,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  use_internal_data=TRUE
)
KEGG@result$Description <- sub(" - .*", "", KEGG@result$Description)
write.csv(KEGG@result, file = file.path(outdir, "KEGG_enrich_result.csv"), row.names = FALSE)

if (nrow(KEGG) > 0) {
  p3 <- dotplot(KEGG, x = "GeneRatio", color = "p.adjust", showCategory = 10, label_format = 100)
  ggsave(file.path(outdir, "KEGG_dotplot.pdf"), p3, width = 12, height = 5)
  ggsave(file.path(outdir, "KEGG_dotplot.png"), p3, width = 12, height = 5)

  p4 <- barplot(KEGG, x = "Count", color = "p.adjust", showCategory = 10, label_format = 80)
  ggsave(file.path(outdir, "KEGG_barplot.pdf"), p4, width = 12, height = 5)
  ggsave(file.path(outdir, "KEGG_barplot.png"), p4, width = 12, height = 5)
}
