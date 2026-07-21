log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(reshape2)
library(ggplot2)
library(dplyr)

select_params <- read.table(snakemake@input[['Recommend']], header=TRUE, check.names=FALSE)

p_type <- snakemake@params[["p_type"]]
FC_cutoff <- snakemake@params[["FC_cutoff"]]
if (p_type == "auto") {
  p_type <- select_params[1,"p_type"]
} else {
  p_type <- snakemake@params[["p_type"]]
}

if (FC_cutoff == "auto") {
  FC_cutoff <- as.numeric(select_params[1,"fc"])
} else {
  FC_cutoff <- as.numeric(snakemake@params[["FC_cutoff"]])
}

log2FC_cutoff <- log2(FC_cutoff)

p_cutoff <- snakemake@params[["p_cutoff"]]

tool <- snakemake@params[["tool"]]
tsv_output <- snakemake@output[[1]]
png_output <- snakemake@output[[2]]
pdf_output <- snakemake@output[[3]]

diff_stats_df <- data.frame(contrast = character(), 
                            all = numeric(),
                            up = numeric(), 
                            down = numeric(), 
                            threshold = character()) 

if (p_type == "padj") {
  for (f in snakemake@input[['EGs']]) {
    f_df <- read.table(file=f, header=TRUE, sep="\t", check.names=FALSE)
    contrast <- gsub("_diffexp.xls", "", basename(f))
    f_df$significance  <- as.factor(ifelse(is.na(f_df$padj), 'NOT', 
                                           ifelse(f_df$padj < p_cutoff & f_df$log2FoldChange > log2FC_cutoff, 'UP', 
                                                  ifelse(f_df$padj < p_cutoff & f_df$log2FoldChange < -log2FC_cutoff, 'DOWN', 'NOT'))))
    gene_up <- nrow(f_df[f_df$significance=="UP",])
    gene_down <- nrow(f_df[f_df$significance=="DOWN",])
    gene_all <- gene_up + gene_down
    threshold <- paste0("padj < ", p_cutoff, " & |FoldChange| > ", FC_cutoff)
    diff_stats_df[nrow(diff_stats_df) + 1, ] <- c(contrast, gene_all, gene_up, gene_down, threshold)
  }
} else if (p_type == "pvalue") {
  for (f in snakemake@input[['EGs']]) {
    f_df <- read.table(file=f, header=TRUE, sep="\t", check.names=FALSE)
    contrast <- gsub("_diffexp.xls", "", basename(f))
    f_df$significance  <- as.factor(ifelse(is.na(f_df$pvalue), 'NOT', 
                                           ifelse(f_df$pvalue < p_cutoff & f_df$log2FoldChange > log2FC_cutoff, 'UP', 
                                                  ifelse(f_df$pvalue < p_cutoff & f_df$log2FoldChange < -log2FC_cutoff, 'DOWN', 'NOT'))))
    gene_up <- nrow(f_df[f_df$significance=="UP",])
    gene_down <- nrow(f_df[f_df$significance=="DOWN",])
    gene_all <- gene_up + gene_down
    threshold <- paste0("pvalue < ", p_cutoff, " & |FoldChange| > ", FC_cutoff)
    diff_stats_df[nrow(diff_stats_df) + 1, ] <- c(contrast, gene_all, gene_up, gene_down, threshold)
  }
}


write.table(diff_stats_df, file=tsv_output, row.names=FALSE, sep='\t')

plot_df <- diff_stats_df %>% 
  select(-threshold)

plot_df <- melt(plot_df, id.vars=c("contrast"))

plot_df$variable <- factor(plot_df$variable, levels = c("all", "down", "up"))
plot_df$value <- as.numeric(plot_df$value)

color <- c("all" = "#E41A1C",
           "down" = "#377EB8",
           "up" = "#4DAF4A")

p <- ggplot(data=plot_df, aes(x=contrast, y=value, fill=variable)) +
  geom_bar(stat = "identity", 
           position = 'dodge',
           width = 0.6) + 
  geom_text(aes(label=value), 
            position = position_dodge(width=0.6), 
            size = 4, 
            vjust=-0.5) + 
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold",
                                  hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  labs(x = "Contrast Group", 
       y = "DEG Counts",
       title = threshold) +
  scale_fill_manual(values = color) + 
  guides(fill = guide_legend(title="Type"))

ggsave(filename = png_output, 
       plot = p,
       width = 10,
       height = 8)
ggsave(filename = pdf_output, 
       plot = p,
       width = 10,
       height = 8)
