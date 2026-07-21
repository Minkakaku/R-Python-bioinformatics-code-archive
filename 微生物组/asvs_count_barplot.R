# R script to generate barplot of ASVs count distribution
# by Yewei on 2024.4.24
# Output: {opts$outdir}/{opts$prefix}.ASVs_count.barplot.{png,pdf,xls}

library(optparse)
library(dplyr)
library(ggplot2)
library(glue)
library(reshape2)
library(RColorBrewer)
library(scales)

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
    make_option(c("-f", "--font"), type="integer", default=12,
                help="Font size for sample names (12 does a good job for sample size less than 30) [default %default]"),
    make_option(c("-w", "--plot_width"), type="numeric", default=10,
                help="Plot width [default %default]"),
    make_option(c("-e", "--plot_height"), type="numeric", default=8,
                help="Plot height [default %default]"),
    make_option(c("-b", "--bar_width"), type="numeric", default=0.9,
                help="Bar width [default %default]")
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
  sample_order <- meta_df$sample
} else {
  meta_df <- read.table(opts$meta, header=T, check.names=F, row.names=NULL, sep="\t") %>%
    dplyr::rename("sample" = opts$sample, "group" = opts$group) %>%
    dplyr::select(sample, group) %>%
    dplyr::filter(group %in% group_list)
  # Get sample order
  sample_order <- meta_df %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(sample, .by_group=TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(factor(group, levels = group_list)) %>%
    dplyr::pull(sample)
}

# Read Merged ASVs count table
asv_count_df <- read.table(opts$input, check.names=F, header=T, sep="\t")
asv_count_df <- asv_count_df %>%
  dplyr::select(colnames(asv_count_df)[1], meta_df$sample)

# Generate barplot
plot_df <- asv_count_df[rowSums(asv_count_df[2:ncol(asv_count_df)])>0,]
plot_df <- melt(plot_df, id.vars=colnames(asv_count_df)[1], value.name="count")
plot_df <- plot_df %>%
  dplyr::filter(count != 0)
count_total <- length(unique(plot_df[[1]]))
plot_df <- plot_df %>%
  dplyr::rename(Sample = variable) %>%
  dplyr::group_by(Sample) %>%
  summarise(ASVs_count = n())
plot_df <- plot_df %>%
  add_row(Sample = "Total", ASVs_count = count_total)
plot_df$Sample <- factor(plot_df$Sample, levels=c(sample_order, "Total"))

max <- 500 * ceiling(max(plot_df$ASVs_count)/500)

p <- ggplot(plot_df, aes(y=ASVs_count, x=Sample)) +
  geom_bar(stat="identity", fill="#377EB8", width=opts$bar_width) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = opts$font),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = "bold")) +
  labs(title = glue("ASVs Count Distribution: {opts$contrast}"),
       x = "Sample ID",
       y = "ASVs Count") +
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     limits = c(0, max),
                     breaks = seq(0,max,500),
                     labels = label_comma())

# Output
dir.create(file.path(opts$outdir))
write.table(plot_df,
            file=file.path(opts$outdir, glue("{opts$prefix}.ASVs_count.barplot.xls")),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
ggsave(filename = file.path(opts$outdir, glue("{opts$prefix}.ASVs_count.barplot.png")),
       device = "png",
       plot = p,
       width = opts$plot_width,
       height = opts$plot_height)
ggsave(filename = file.path(opts$outdir, glue("{opts$prefix}.ASVs_count.barplot.pdf")),
       device = "pdf",
       plot = p,
       width = opts$plot_width,
       height = opts$plot_height)
