# R script to generate heatmap based on taxonomic/functional abundance matrix
# by Yewei on 2024.1.4
# Output: {opts$outdir}/top{opts$topn}.{opts$level}.heatmap.{png,pdf}

library(optparse)
library(dplyr)
library(pheatmap)
library(glue)
library(reshape2)
library(tidyverse)

# Parse arguments
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Matrix file of taxonomic abundace for all samples"),
    make_option(c("-m", "--metadata"), type="character", default=NULL,
                help="Metadata file for all samples"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="Group variable name in metadata [default %default]"),
    make_option(c("-s", "--sample"), type="character", default="sample",
                help="Sample variable name in metadata [default %default]"),
    make_option(c("-l", "--level"), type="character", default="Species",
                help="Taxonomic level (Kingdom, Phylum, Class, Order, Family, Genus, Species) or Functional database type (KEGG_KO, KEGG_level1, KEGG_level2, KEGG_level3, CAZy, COG_category) [default %default]"),
    make_option(c("-o", "--outdir"), type="character", default=".",
                help="Output directory [default %default]"),
    make_option(c("-t", "--topn"), type="integer", default=0,
                help="Top N items with highest abundance to present ('0' to display all) [default %default]"),
    make_option(c("-f", "--font"), type="integer", default=6,
                help="Font size for sample names (6 does a good job for sample size less than 30) [default %default]"),
    make_option(c("-c", "--cluster_cols"), type="logical", action = "store_true", default = FALSE,
                help="Cluster columns in heatmap [default %default]"),
    make_option(c("-e", "--contrast"), type="character", default="all",
                help="Display samples in contrast groups ('all' or 'A-vs-B') [default %default]"),
    make_option(c("-u", "--unassigned"), type="logical", action="store_true", default=FALSE,
                help="If there is an 'Unassigned' row to present [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

# Read meta data
meta_df <- read.table(opts$metadata, header=T, row.names=opts$sample, sep="\t") %>%
  dplyr::select(opts$group) %>%
  dplyr::rename(group = opts$group) %>%
  dplyr::mutate(sample = rownames(.)) %>%
  dplyr::arrange(group, sample) %>%
  dplyr::select(-sample)
if (opts$contrast != "all") {
  group_list <- unlist(strsplit(opts$contrast, "-vs-"))
  meta_df <- meta_df %>%
    dplyr::filter(group %in% group_list)
  # Get sample order
  meta_df_tmp <- meta_df %>%
    dplyr::arrange(factor(group, levels = group_list))
  sample_order <- rownames(meta_df_tmp)
} else {
  sample_order <- rownames(meta_df)
}

# Read abundance matrix (sum up for duplicate taxonomy names)
abundance_df <- read.table(opts$input, header=T, check.names=F, sep="\t", quote="") %>%
  dplyr::select(opts$level, rownames(meta_df)) %>%
  dplyr::group_by(across(all_of(opts$level))) %>%
  dplyr::summarise_all(sum) %>%
  column_to_rownames(var=opts$level)

# Remove rows with all counts equal zero otherwise pheatmap will throw an error (hclust)
abundance_df <- abundance_df[rowSums(abundance_df) != 0,]

# Sort and select topN rows
abundance_sort_df <- abundance_df[order(rowMeans(abundance_df), decreasing = T),]

if (! opts$unassigned) {
  abundance_sort_df <- abundance_sort_df[rownames(abundance_sort_df)!="Unassigned",]
}

if (opts$topn == 0) {
  abundance_top_df <- abundance_sort_df
} else if (nrow(abundance_sort_df) >= opts$topn) {
  abundance_top_df <- abundance_sort_df[1:opts$topn,]
} else {
  abundance_top_df <- abundance_sort_df
}

# Plot
plot_df <- abundance_top_df[, sample_order]
pheatmap::pheatmap(plot_df,
                   show_colnames = T,
                   show_rownames = T,
                   scale = "row",
                   fontsize = 6,
                   fontsize_row = 6,
                   fontsize_col = opts$font,
                   angle_col = 45,
                   cluster_rows = T,
                   cluster_cols = opts$cluster_cols,
                   annotation_col = meta_df,
                   filename = file.path(opts$outdir, glue("{opts$contrast}.{opts$level}.heatmap.pdf")))
pheatmap::pheatmap(plot_df,
                   show_colnames = T,
                   show_rownames = T,
                   scale = "row",
                   fontsize = 6,
                   fontsize_row = 6,
                   fontsize_col = opts$font,
                   angle_col = 45,
                   cluster_rows = T,
                   cluster_cols = opts$cluster_cols,
                   annotation_col = meta_df,
                   filename = file.path(opts$outdir, glue("{opts$contrast}.{opts$level}.heatmap.png")))
