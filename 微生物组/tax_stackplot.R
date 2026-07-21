# R script to generate stack barplots based on relative taxonomic abundance matrix for AmpliSeq workflow
# by Yewei on 2024.4.9
# Output: {opts$outdir}/{opts$contrast}.{opts$level}.barplot.{png,pdf,xls}

library(optparse)
library(dplyr)
library(ggplot2)
library(glue)
library(reshape2)
library(RColorBrewer)

# Parse arguments
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Matrix file of relative taxonomic abundace for all samples"),
    make_option(c("-m", "--metadata"), type="character", default=NULL,
                help="Metadata file for all samples"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="Group variable name in metadata [default %default]"),
    make_option(c("-s", "--sample"), type="character", default="sample",
                help="Sample variable name in metadata [default %default]"),
    make_option(c("-l", "--level"), type="character", default="Species",
                help="Taxonomic level (Kingdom, Phylum, Class, Order, Family, Genus, Species) [default %default]"),
    make_option(c("-o", "--outdir"), type="character", default=".",
                help="Output directory [default %default]"),
    make_option(c("-t", "--topn"), type="integer", default=10,
                help="Top N items with highest abundance to present [default %default]"),
    make_option(c("-e", "--others"), type="logical", action="store_true", default=FALSE,
                help="If there is an 'Others' row other than Top N items to present [default %default]"),
    make_option(c("-f", "--font"), type="integer", default=12,
                help="Font size for sample names (12 does a good job for sample size less than 30) [default %default]"),
    make_option(c("-c", "--contrast"), type="character", default="all",
                help="Display samples in contrast groups ('all' or 'A-vs-B') [default %default]"),
    make_option(c("-u", "--unassigned"), type="logical", action="store_true", default=FALSE,
                help="If there is an 'Unassigned' row to present [default %default]"),
    make_option(c("-b", "--bar_width_sample"), type="numeric", default=0.9,
                help="Bar width in sample-level stack barplot [default %default]"),
    make_option(c("-d", "--bar_width_group"), type="numeric", default=0.5,
                help="Bar width in group-level stack barplot [default %default]"),
    make_option(c("-j", "--plot_width_sample"), type="numeric", default=10,
                help="Plot width in sample-level stack barplot [default %default]"),
    make_option(c("-k", "--plot_height_sample"), type="numeric", default=8,
                help="Plot height in sample-level stack barplot [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

# Read meta data
meta_df <- read.table(opts$metadata, check.names=F, header=T, sep="\t") %>%
  dplyr::select(opts$sample, opts$group) %>%
  dplyr::rename(sample = opts$sample, group = opts$group)
if (opts$contrast != "all") {
  group_list <- unlist(strsplit(opts$contrast, "-vs-"))
  meta_df <- meta_df %>%
    dplyr::filter(group %in% group_list)
  # Get sample order
  sample_order <- meta_df %>%
    dplyr::group_by(group) %>%
    dplyr::arrange(sample, .by_group=TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(factor(group, levels = group_list)) %>%
    dplyr::pull(sample)
} else {
  sample_order <- meta_df$sample
  group_list <- unique(meta_df$group)
}

# Read abundance matrix (sum up for duplicate taxonomy names)
abundance_df <- read.table(opts$input, header=T, sep="\t", check.names=F, quote="")
abundance_df <- abundance_df %>%
  dplyr::select(opts$level, meta_df$sample) %>%
  dplyr::group_by(across(all_of(opts$level))) %>%
  dplyr::summarise_all(sum)

# Sort based on rowMeans and select the topN rows
abundance_sort_df <- abundance_df[order(rowMeans(abundance_df[, 2:length(colnames(abundance_df))]), decreasing = T),]

if (opts$unassigned) {
  if (nrow(abundance_sort_df) >= opts$topn) {
    abundance_top_df <- abundance_sort_df[1:opts$topn,]
  } else {
    abundance_top_df <- abundance_sort_df
  }

  if (opts$others) {
    others_line <- c("Others", as.vector(1 - colSums(abundance_sort_df[abundance_sort_df[,1]=="Unassigned",][, 2:ncol(abundance_sort_df)]) - colSums(abundance_top_df[, 2:ncol(abundance_top_df)])))
    abundance_top_df <- rbind(abundance_top_df, others_line)
  }

  abundance_top_df <- rbind(abundance_top_df, abundance_sort_df[abundance_sort_df[,1]=="Unassigned",])
} else {
  if (nrow(abundance_sort_df) - 1 >= opts$topn) {
    abundance_top_df <- abundance_sort_df[abundance_sort_df[,1]!="Unassigned",][1:opts$topn,]
  } else {
    abundance_top_df <- abundance_sort_df[abundance_sort_df[,1]!="Unassigned",]
  }

  if (opts$others) {
    others_line <- c("Others", as.vector(1 - colSums(abundance_top_df[, 2:ncol(abundance_top_df)])))
    abundance_top_df <- rbind(abundance_top_df, others_line)
  }
}

order_vector <- rev(abundance_top_df[[1]])

################
# Sample-level #
################

# Data

write.table(abundance_top_df, 
            file=file.path(opts$outdir, glue("{opts$contrast}.{opts$level}.barplot.xls")),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
          
# Plot
abundance_top_long_df <- melt(abundance_top_df, id.vars = opts$level) %>% 
  rename(sample = variable)
abundance_top_long_df$value <- as.numeric(abundance_top_long_df$value)
plot_df <- merge(x=abundance_top_long_df, 
                 y=meta_df, 
                 by.x="sample", 
                 by.y="sample",
                 all.x=TRUE)
plot_df[[opts$level]] <- factor(plot_df[[opts$level]], levels=order_vector)

plot_df$sample <- factor(plot_df$sample, levels=sample_order)

# Color palette
color_count <- nrow(abundance_top_df)
color_pal <- brewer.pal(n = 12, "Paired")
if (color_count <= 12) {
  get_palette = colorRampPalette(color_pal[1:color_count])
} else {
  get_palette = colorRampPalette(color_pal)
}

p <- ggplot(plot_df, aes(fill=!!sym(opts$level), y=value, x=sample)) +
  geom_bar(position="stack", stat="identity", width=opts$bar_width_sample) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = opts$font, hjust = 1),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = "bold")) +
  labs(x = "Sample", 
       y = "Relative Abundance") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values=rev(get_palette(color_count)))

ggsave(filename = file.path(opts$outdir, glue("{opts$contrast}.{opts$level}.barplot.png")),
       device = "png",
       plot = p,
       width = opts$plot_width_sample,
       height = opts$plot_height_sample)
ggsave(filename = file.path(opts$outdir, glue("{opts$contrast}.{opts$level}.barplot.pdf")),
       device = "pdf",
       plot = p,
       width = opts$plot_width_sample,
       height = opts$plot_height_sample)

###############
# Group-level #
###############

# Data

abundance_top_group_df <- abundance_top_df
for (i in 2:ncol(abundance_top_group_df)) {
  abundance_top_group_df[i] <- as.numeric(abundance_top_group_df[[i]])
}

for (g in group_list) {
  samples <- meta_df[meta_df$group==g, "sample"]
  abundance_top_group_df[g] <- rowMeans(abundance_top_group_df[, samples])
}

abundance_top_group_df <- abundance_top_group_df[, c(opts$level, group_list)]

write.table(abundance_top_group_df,
            file=file.path(opts$outdir, glue("{opts$contrast}.{opts$level}.group.barplot.xls")),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Plot
plot_df <- melt(abundance_top_group_df, id.vars = opts$level) %>%
  rename(group = variable)
plot_df$value <- as.numeric(plot_df$value)

plot_df[[opts$level]] <- factor(plot_df[[opts$level]], levels=order_vector)

plot_df$group <- factor(plot_df$group, levels=group_list)

# Color palette
color_count <- nrow(abundance_top_group_df)
color_pal <- brewer.pal(n = 12, "Paired")
if (color_count <= 12) {
  get_palette = colorRampPalette(color_pal[1:color_count])
} else {
  get_palette = colorRampPalette(color_pal)
}

p <- ggplot(plot_df, aes(fill=!!sym(opts$level), y=value, x=group)) +
  geom_bar(position="stack", stat="identity", width=opts$bar_width_group) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = "bold")) +
  labs(x = "Group",
       y = "Relative Abundance") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values=rev(get_palette(color_count)))

ggsave(filename = file.path(opts$outdir, glue("{opts$contrast}.{opts$level}.group.barplot.png")),
       device = "png",
       plot = p,
       width = 10,
       height = 8)
ggsave(filename = file.path(opts$outdir, glue("{opts$contrast}.{opts$level}.group.barplot.pdf")),
       device = "pdf",
       plot = p,
       width = 10,
       height = 8)
