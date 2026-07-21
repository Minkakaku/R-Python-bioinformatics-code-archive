log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(ggplot2)
# library(pheatmap)
# library(circlize)
# library(ComplexHeatmap)
library(pheatmap)

# Input
## DEG results
DEG <- read.table(snakemake@input[[1]], sep="\t", header=TRUE, row.names=1, check.names=FALSE)
## Cutoff setup
select_params <- read.table(snakemake@input[[3]], header=TRUE, check.names=FALSE)

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

## Contrast
contrast <- gsub("_diffexp.xls", "", basename(snakemake@input[[1]]))
## Group
coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)
coldata <- coldata["group"]
coldata <- coldata[coldata$group %in% snakemake@params[["contrast"]], , drop=F]
coldata <- coldata[order(row.names(coldata)), , drop=F]
## Count table
count <- read.table(snakemake@input[[2]], sep="\t", header=TRUE, row.names=1, check.names=FALSE)[row.names(coldata)]

# Output
output_volcano_png <- snakemake@output[["volcano_png"]]
output_volcano_pdf <- snakemake@output[["volcano_pdf"]]
output_heatmap_png <- snakemake@output[["heatmap_png"]]
output_heatmap_pdf <- snakemake@output[["heatmap_pdf"]]

# Filter DESeq2 results
if (p_type == "padj") {
  need_DEG <- DEG[,c("log2FoldChange","pvalue","padj")]
  need_DEG <- need_DEG[!is.na(need_DEG$padj),]
  need_DEG$significance  <- as.factor(ifelse(is.na(need_DEG$padj), 'Normal', 
                                             ifelse(need_DEG$padj < p_cutoff & need_DEG$log2FoldChange > log2FC_cutoff, 'Up',
                                                    ifelse(need_DEG$padj < p_cutoff & need_DEG$log2FoldChange < -log2FC_cutoff, 'Down', 'Normal'))))
} else if (p_type == "pvalue") {
  need_DEG <- DEG[,c("log2FoldChange","pvalue","padj")]
  need_DEG <- need_DEG[!is.na(need_DEG$pvalue),]
  need_DEG$significance  <- as.factor(ifelse(is.na(need_DEG$pvalue), 'Normal', 
                                             ifelse(need_DEG$pvalue < p_cutoff & need_DEG$log2FoldChange > log2FC_cutoff, 'Up',
                                                    ifelse(need_DEG$pvalue < p_cutoff & need_DEG$log2FoldChange < -log2FC_cutoff, 'Down', 'Normal'))))
}

# Volcano plot

volcano_df <- need_DEG
up_n = nrow(volcano_df[volcano_df$significance == "Up",])
normal_n = nrow(volcano_df[volcano_df$significance == "Normal",])
down_n = nrow(volcano_df[volcano_df$significance == "Down",])

if (up_n == 0 & down_n == 0) {
  volcano_df$significance <- factor(volcano_df$significance, levels = c("Normal"))
  levels(volcano_df$significance) <- c(paste0("Normal (", normal_n, ")"))
  color_pal <- c("grey")
} else if (up_n != 0 & down_n == 0) {
  volcano_df$significance <- factor(volcano_df$significance, levels = c("Up","Normal"))
  levels(volcano_df$significance) <- c(paste0("Up (", up_n, ")"),
                                       paste0("Normal (", normal_n, ")"))
  color_pal <- c('red','grey')
} else if (up_n == 0 & down_n != 0) {
  volcano_df$significance <- factor(volcano_df$significance, levels = c("Normal", "Down"))
  levels(volcano_df$significance) <- c(paste0("Normal (", normal_n, ")"),
                                       paste0("Down (", down_n, ")"))
  color_pal <- c('grey','blue')
} else {
  volcano_df$significance <- factor(volcano_df$significance, levels = c("Up", "Normal", "Down"))
  levels(volcano_df$significance) <- c(paste0("Up (", up_n, ")"),
                                       paste0("Normal (", normal_n, ")"),
                                       paste0("Down (", down_n, ")"))
  color_pal <- c('red','grey','blue')
}

if (p_type == "padj") {
  g <- ggplot(data = volcano_df, 
              aes(x=log2FoldChange, y=-log10(padj), color=significance)) +
    geom_point(alpha=0.4, size=1) +
    theme_bw() +
    xlab("log2(FoldChange)") + 
    ylab("-log10(padj)") +
    ggtitle(contrast) + 
    scale_colour_manual(paste0("padj < ", p_cutoff, "\n", "|FoldChange| > ", FC_cutoff),
                        values=color_pal) +
    geom_vline(xintercept=c(-log2FC_cutoff, log2FC_cutoff), lty=4, col="grey", lwd=0.8) +
    geom_hline(yintercept=-log10(p_cutoff), lty=4, col="grey", lwd=0.8) +
    theme(plot.title = element_text(face = "bold", hjust=0.5), 
          plot.margin =unit(c(2,2,2,2), 'lines'),
          legend.position ="right", 
          axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
} else if (p_type == "pvalue") {
  g <- ggplot(data = volcano_df, 
              aes(x=log2FoldChange, y=-log10(pvalue), color=significance)) +
    geom_point(alpha=0.4, size=1) +
    theme_bw() +
    xlab("log2(FoldChange)") + 
    ylab("-log10(pvalue)") +
    ggtitle(contrast) +
    scale_colour_manual(paste0("pvalue < ", p_cutoff, "\n", "|FoldChange| > ", FC_cutoff),
                        values=color_pal) +
    geom_vline(xintercept=c(-log2FC_cutoff, log2FC_cutoff), lty=4, col="grey", lwd=0.8) +
    geom_hline(yintercept=-log10(p_cutoff), lty=4, col="grey", lwd=0.8) +
    theme(plot.title = element_text(face = "bold", hjust=0.5), 
          plot.margin =unit(c(2,2,2,2), 'lines'),
          legend.position ="right", 
          axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
}

ggsave(g, filename=output_volcano_png, width=10, height=8)
ggsave(g, filename=output_volcano_pdf, width=10, height=8)

# Heatmap plot
gene_up <- rownames(need_DEG[need_DEG$significance=="Up",])
gene_down <- rownames(need_DEG[need_DEG$significance=="Down",])
cg <- c(gene_up, gene_down)

if (length(cg) >= 3) {
  if (snakemake@params[["tool"]] == "DESeq2") {
    dat <- log2(count + 1)
    n <- dat[cg,]
    exp <- apply(n, 1, scale)
    rownames(exp) <- colnames(n)
    exp <- t(exp)
    pheatmap::pheatmap(exp,
                       show_colnames = T,
                       show_rownames = F,
                       #scale = "row",
                       fontsize = 7 ,
                       angle_col = 45,
                       cluster_cols = T,
                       annotation_col = coldata,
                       filename = output_heatmap_png)
    pheatmap::pheatmap(exp,
                       show_colnames = T,
                       show_rownames = F,
                       #scale = "row",
                       fontsize = 7 ,
                       angle_col = 45,
                       cluster_cols = T,
                       annotation_col = coldata,
                       filename = output_heatmap_pdf)
  } else if (snakemake@params[["tool"]] == "edgeR") {
    dat <- count
    n <- dat[cg,]
    exp <- apply(n, 1, scale)
    rownames(exp) <- colnames(n)
    exp <- t(exp)
    pheatmap::pheatmap(exp,
                       show_colnames = T,
                       show_rownames = F,
                       #scale = "row",
                       fontsize = 7 ,
                       angle_col = 45,
                       cluster_cols = T,
                       annotation_col = coldata,
                       filename = output_heatmap_png)
    pheatmap::pheatmap(exp,
                       show_colnames = T,
                       show_rownames = F,
                       #scale = "row",
                       fontsize = 7 ,
                       angle_col = 45,
                       cluster_cols = T,
                       annotation_col = coldata,
                       filename = output_heatmap_pdf)
  }
} else {
  p <- ggplot() +
    theme_void() +
    geom_text(aes(0,0,label='N/A'), size=15) +
    xlab(NULL) #optional, but safer in case another theme is applied later
  ggsave(p, filename=output_heatmap_png, width=10, height=8)
  ggsave(p, filename=output_heatmap_pdf, width=10, height=8)
}


