rm(list = ls())
# 对data绘制火山图
# 筛选特征基因（pvals_adj<0.05,logfoldchanges > 1）进行hallmark数据库富集分析
library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)

Degs <- "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/Degs/Epith_Cluster1_8_marker.xls"
outdir <- "/mnt/NexData_file/home/wmx/Project/AE2507033447-08/shouhou13/Hallmark_enrich"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

data <- read.delim2(Degs, header = TRUE, check.names = FALSE)

colnames(data)
head(data)

data <- data %>%
  mutate(
    logfoldchanges = as.numeric(logfoldchanges),
    pvals_adj = as.numeric(pvals_adj),
    pvals = as.numeric(pvals),
    sig = case_when(
      pvals_adj < 0.05 & logfoldchanges > 1  ~ "Up",
      pvals_adj < 0.05 & logfoldchanges < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    neg_log10_fdr = -log10(pvals_adj)
  )

p_volcano <- ggplot(data, aes(x = logfoldchanges, y = neg_log10_fdr)) +
  geom_point(aes(color = sig), size = 1.2, alpha = 0.75) +
  scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey75")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(
    title = "Volcano plot",
    x = "logfoldchanges",
    y = "-log10(FDR)",
    color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(outdir, "volcano_plot.pdf"),
  plot = p_volcano,
  width = 7,
  height = 6
)

ggsave(
  filename = file.path(outdir, "volcano_plot.png"),
  plot = p_volcano,
  width = 7,
  height = 6,
  dpi = 300
)

deg_up <- data %>%
  filter(pvals_adj < 0.05, logfoldchanges > 1) %>%
  pull(var = names) %>%
  unique() %>%
  na.omit()

write.table(
  deg_up,
  file = file.path(outdir, "up_genes_for_hallmark.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

cat("Up genes:", length(deg_up), "\n")

m_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

hallmark_enrich <- enricher(
  gene = deg_up,
  TERM2GENE = m_df,
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  pAdjustMethod = "BH"
)

hallmark_res <- as.data.frame(hallmark_enrich)

write.csv(
  hallmark_res,
  file = file.path(outdir, "hallmark_enrichment_up_genes.csv"),
  row.names = FALSE
)

if (nrow(hallmark_res) > 0) {
  
  p_dot <- dotplot(
    hallmark_enrich,
    showCategory = min(20, nrow(hallmark_res))
  ) +
    ggtitle("Hallmark enrichment of up-regulated genes") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    filename = file.path(outdir, "hallmark_enrichment_dotplot.pdf"),
    plot = p_dot,
    width = 8,
    height = 12
  )
  
  ggsave(
    filename = file.path(outdir, "hallmark_enrichment_dotplot.png"),
    plot = p_dot,
    width = 8,
    height = 12,
    dpi = 300
  )
}

cat("Done. Results saved to:", outdir, "\n")
