# R script to generate venn plot of ASVs count distribution
# by Yewei on 2024.4.24
# Output: {opts$outdir}/{opts$prefix}.ASVs_count.barplot.{png,pdf,xls}

library(optparse)
library(dplyr)
library(glue)
library(reshape2)
library(venn)

# Parse arguments
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Merged ASVs count table of all samples"),
    make_option(c("-m", "--metadata"), type="character", default=NULL,
                help="Metadata file for all samples"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="Group variable name in metadata [default %default]"),
    make_option(c("-s", "--sample"), type="character", default="sample",
                help="Sample variable name in metadata [default %default]"),
    make_option(c("-c", "--contrast"), type="character", default="all",
                help="Display samples in contrast groups ('all' or 'A-vs-B') [default %default]"),
    make_option(c("-o", "--outdir"), type="character", default=".",
                help="Output directory [default %default]"),
    make_option(c("-p", "--prefix"), type="character", default="output",
                help="Output file name prefix [default %default]"),
    make_option(c("-e", "--ellipse"), type="logical", action="store_true", default=FALSE,
                help="Whether to force the shape to an ellipse, where possible [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

# Convert contrast to group list
if (opts$contrast != "all") {
  group_list <- unlist(strsplit(opts$contrast, "-vs-"))
}

# Read metadata
if (opts$contrast == "all") {
  meta_df <- read.table(opts$meta, header=T, check.names=F, row.names=NULL, sep="\t") %>%
    dplyr::rename("sample" = opts$sample, "group" = opts$group) %>%
    dplyr::select(sample, group)
  group_list <- unique(meta_df$group)
} else {
  meta_df <- read.table(opts$meta, header=T, check.names=F, row.names=NULL, sep="\t") %>%
    dplyr::rename("sample" = opts$sample, "group" = opts$group) %>%
    dplyr::select(sample, group) %>%
    dplyr::filter(group %in% group_list)
}

# Read Merged ASVs count table
asv_count_df <- read.table(opts$input, check.names=F, header=T, sep="\t")
asv_count_df <- asv_count_df %>%
  dplyr::select(colnames(asv_count_df)[1], meta_df$sample)

# Generate barplot
plot_df <- asv_count_df
for (group in group_list) {
  samples <- meta_df[meta_df$group==group, "sample"]
  sums <- rowSums(asv_count_df[samples])
  plot_df[,group] <- sums
}

plot_df <- plot_df %>%
  dplyr::select(colnames(plot_df)[1], group_list)

venn_list <- list()
for (group in group_list) {
  venn_list <- append(venn_list, list(plot_df[plot_df[,group]>0, 1]))
}
names(venn_list) <- group_list

# Output
dir.create(file.path(opts$outdir))
write.table(plot_df,
            file=file.path(opts$outdir, glue("{opts$prefix}.ASVs_count.venn.xls")),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

pdf(file = file.path(opts$outdir, glue("{opts$prefix}.ASVs_count.venn.pdf")),
    width = 10,
    height = 10)
venn(venn_list,
     zcolor='style',
     ilabels = "counts",
     opacity = 0.3,
     box = F,
     ilcs = 1.5,
     sncs = 1.5,
     ellipse = opts$ellipse)
dev.off()

png(filename = file.path(opts$outdir, glue("{opts$prefix}.ASVs_count.venn.png")),
    width = 10,
    height = 10,
    units = "in",
    res=100)
venn(venn_list,
     zcolor='style',
     ilabels = "counts",
     opacity = 0.3,
     box = F,
     ilcs = 1.5,
     sncs = 1.5,
     ellipse = opts$ellipse)
dev.off()