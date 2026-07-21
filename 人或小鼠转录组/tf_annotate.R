log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)

# Input
## DEG results
DEG <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, check.names=FALSE)
## Cutoff setup
tfdb <- read.table(snakemake@params[["db"]], sep="\t", header=TRUE)

# Output
output_tf <- snakemake@output[[1]]

need_DEG <- DEG[, c("gene","log2FoldChange","pvalue","padj","regulate")]

# Merge
tf_results_df <- merge(need_DEG, tfdb, by.x="gene", by.y="Ensembl") %>% 
  select(gene, Symbol, log2FoldChange, padj, regulate, Family) %>%
  arrange(padj)
colnames(tf_results_df) <- c("Gene_ID", "Symbol", "log2FoldChange", "padj", "Regulated", "Family")
write.table(tf_results_df, file=output_tf, row.names=FALSE, sep='\t')