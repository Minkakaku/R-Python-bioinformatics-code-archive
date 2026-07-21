# R script to generate heatmap based on beta distance matrix among samples
# by Yewei on 2024.4.16
# Output: {opts$outdir}/{opts$prefix}.heatmap.{png,pdf}

library(optparse)
library(dplyr)
library(pheatmap)
library(glue)
library(reshape2)

# Parse arguments
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="beta diversity distance matrix"),
    make_option(c("-m", "--metadata"), type="character", default=NULL,
                help="Metadata file for all samples"),
    make_option(c("-e", "--contrast"), type="character", default=NULL,
                help="Contrast with groups to be included"),
    make_option(c("-d", "--dist"), type="character", default="braycurtis",
                help="Beta distance (braycurtis, jaccard, weighted_unifrac, unweighted_unifrac) [default %default]"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="Group variable name in metadata [default %default]"),
    make_option(c("-s", "--sample"), type="character", default="sample",
                help="Sample variable name in metadata [default %default]"),
    make_option(c("-o", "--outdir"), type="character", default=".",
                help="Output directory [default %default]"),
    make_option(c("-f", "--font"), type="integer", default=6,
                help="Font size for sample names (6 does a good job for sample size less than 30) [default %default]"),
    make_option(c("-c", "--cluster"), type="logical", action = "store_true", default = FALSE,
                help="Cluster columns and rows in heatmap [default %default]"),
    make_option(c("-p", "--prefix"), type="character", default="output",
                help="Output file name prefix [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

# Read beta diversity distance matrix
dist_matrix_df <- read.table(opts$input, header=T, check.names=F, row.names=1, sep="\t", quote="")

# Read contrast groups
group_list <- unlist(strsplit(opts$contrast, "-vs-"))

# Beta distance type
if (opts$dist == "braycurtis") {
  dist_type <- "Bray-Curtis Distance"
} else if (opts$dist == "jaccard") {
  dist_type <- "Jaccard Distance"
} else if (opts$dist == "weighted_unifrac") {
  dist_type <- "Weighted UniFrac Distance"
} else if (opts$dist == "unweighted_unifrac") {
  dist_type <- "Unweighted UniFrac Distance"
}
          
# Plot
meta_df <- read.table(opts$metadata, header=T, row.names=opts$sample, sep="\t") %>%
  dplyr::select(opts$group) %>%
  dplyr::rename(group = opts$group) %>%
  dplyr::filter(group %in% group_list) %>%
  dplyr::mutate(sample = rownames(.)) %>%
  dplyr::arrange(group, sample) %>%
  dplyr::select(-sample)
# Get sample order
meta_df_tmp <- meta_df %>%
  dplyr::arrange(factor(group, levels = group_list))
sample_order <- rownames(meta_df_tmp)

plot_df <- dist_matrix_df[sample_order, sample_order]

meta_df$group <- factor(meta_df$group, levels = group_list)

pheatmap::pheatmap(plot_df,
                   main = glue("Heatmap - {dist_type}"),
                   show_colnames = T,
                   show_rownames = T,
                   fontsize_row = opts$font,
                   fontsize_col = opts$font,
                   angle_col = 45,
                   cluster_rows = opts$cluster,
                   cluster_cols = opts$cluster,
                   annotation_col = meta_df,
                   filename = file.path(opts$outdir, glue("{opts$prefix}.heatmap.pdf")))
pheatmap::pheatmap(plot_df,
                   show_colnames = T,
                   show_rownames = T,
                   main = glue("Heatmap - {dist_type}"),
                   fontsize_row = opts$font,
                   fontsize_col = opts$font,
                   angle_col = 45,
                   cluster_rows = opts$cluster,
                   cluster_cols = opts$cluster,
                   annotation_col = meta_df,
                   filename = file.path(opts$outdir, glue("{opts$prefix}.heatmap.png")))
