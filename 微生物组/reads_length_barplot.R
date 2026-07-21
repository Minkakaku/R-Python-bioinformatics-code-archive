# R script to generate barplot of reads length distribution
# by Yewei on 2024.4.23
# Output: {opts$outdir}/{sample}/{sample}.reads_length.{png,pdf,xls}

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
    make_option(c("-l", "--length"), type="character", default=NULL,
                help="ASVs name and length table with tab separator"),
    make_option(c("-o", "--outdir"), type="character", default=".",
                help="Output directory [default %default]"),
    make_option(c("-m", "--max"), type="integer", default=490,
                help="The maxium length value to show in barplot (will be rounded to a number divisible by ten) [default %default]"),
    make_option(c("-n", "--min"), type="integer", default=350,
                help="Reads with length smaller than this value will be merged into one category (will be rounded to a number divisible by ten) [default %default]"),
    make_option(c("-w", "--width"), type="integer", default=10,
                help="Length increment for each bar [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

# Read Merged ASVs count table
asv_count_df <- read.table(opts$input, check.names=F, header=T, sep="\t")
sample_list <- colnames(asv_count_df)[-1]

# Read ASVs name and length table
asv_length_df <- read.table(opts$length, header=F, sep="\t", check.names=F, quote="") %>%
  dplyr::rename(ASV_ID = V1, length = V2)

# Merge two tables
asv_count_length_df <- merge(asv_count_df, asv_length_df, by="ASV_ID")

# Generate barplot for each sample
max_val <- 10 * ceiling(opts$max / 10)
min_val <- 10 * floor(opts$min / 10)

for (sample in sample_list) {
  plot_df <- asv_count_length_df %>%
    dplyr::select(ASV_ID, sample, length) %>%
    dplyr::rename(count = sample) %>%
    dplyr::filter(count != 0) %>%
    dplyr::mutate(length_cat = cut(length, c(0, seq(min_val,max_val,opts$width)), right=FALSE)) %>%
    dplyr::group_by(length_cat) %>%
    dplyr::summarise(sum = sum(count))

  levels <- levels(plot_df$length_cat)
  for (l in levels) {
    if (! l %in% plot_df$length_cat) {
      plot_df <- rbind(plot_df, c(l, 0))
    }
  }
  plot_df$sum <- as.numeric(plot_df$sum)
  plot_df <- plot_df %>%
    dplyr::arrange(length_cat) %>%
    dplyr::rename(Read_Length = length_cat, Read_Count = sum)

  p <- ggplot(plot_df, aes(y=Read_Count, x=Read_Length)) +
    geom_bar(stat="identity", fill="#377EB8") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(face = "bold")) +
    labs(title = glue("Read Length Distribution: {sample}"),
         x = "Read Length (bp)",
         y = "Read Counts") +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       trans = 'log10',
                       limits = c(10^0, 10^5),
                       breaks = c(0,10^(1:5)),
                       labels = label_comma())

  # Output
  dir.create(file.path(opts$outdir, sample))
  write.table(plot_df,
              file=file.path(opts$outdir, sample, glue("{sample}.reads_length.xls")),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  ggsave(filename = file.path(opts$outdir, sample, glue("{sample}.reads_length.png")),
         device = "png",
         plot = p,
         width = 10,
         height = 8)
  ggsave(filename = file.path(opts$outdir, sample, glue("{sample}.reads_length.pdf")),
         device = "pdf",
         plot = p,
         width = 10,
         height = 8)
}
