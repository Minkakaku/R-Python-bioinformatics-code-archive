# R script to perform Anosim analysis based on beta distance matrix among samples
# by Yewei on 2024.4.16
# Output: {opts$outdir}/{opts$prefix}.Anosim.{png,pdf,xls}

library(optparse)
library(dplyr)
library(vegan)
library(glue)
library(RColorBrewer)

# Parse arguments
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="beta diversity distance matrix"),
    make_option(c("-m", "--metadata"), type="character", default=NULL,
                help="Metadata file for all samples"),
    make_option(c("-c", "--contrast"), type="character", default=NULL,
                help="Contrast in groups [default %default]"),
    make_option(c("-d", "--dist"), type="character", default="braycurtis",
                help="Beta distance type (braycurtis, jaccard, weighted_unifrac, unweighted_unifrac) [default %default]"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="Group variable name in metadata [default %default]"),
    make_option(c("-s", "--sample"), type="character", default="sample",
                help="Sample variable name in metadata [default %default]"),
    make_option(c("-o", "--outdir"), type="character", default=".",
                help="Output directory [default %default]"),
    make_option(c("-p", "--prefix"), type="character", default="output",
                help="Output file name prefix [default %default]"),
    make_option(c("-e", "--perms"), type="integer", default=9999,
                help="Number of permutations performed in Anosim [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

# Read beta distance matrix
dist_matrix_df <- read.table(opts$input, header=T, row.names=1, check.names=F, sep="\t", quote="")
dist_matrix_df <- as.data.frame(t(dist_matrix_df))

# Read metadata
meta_df <- read.table(opts$meta, header=T, check.names=F, row.names=NULL, sep="\t") %>% 
  dplyr::rename("sample" = opts$sample, "group" = opts$group) %>%
  dplyr::select(sample, group)
rownames(meta_df) <- meta_df$sample
meta_df <- meta_df %>% 
  dplyr::select(-sample)

# Select samples in two contrast groups
groups <- unlist(strsplit(opts$contrast, split="-vs-"))
meta_df <- meta_df[meta_df$group %in% groups, , drop=FALSE]
dist_matrix_df <- dist_matrix_df[rownames(meta_df), rownames(meta_df), drop=FALSE]

# Perform Anosim
anosim_results <- anosim(dist_matrix_df, meta_df[rownames(dist_matrix_df),], permutations = opts$perms, distance="bray")

# Beta distance type
if (opts$dist == "braycurtis") {
  dist_type <- "Bray-Curtis"
} else if (opts$dist == "jaccard") {
  dist_type <- "Jaccard"
} else if (opts$dist == "weighted_unifrac") {
  dist_type <- "Weighted UniFrac"
} else if (opts$dist == "unweighted_unifrac") {
  dist_type <- "Unweighted UniFrac"
}

# Output plot and dataframe with R and P
R <- signif(anosim_results$statistic, 4)
P <- signif(anosim_results$signif, 4)

# Sort group order
anosim_results$class.vec <- factor(anosim_results$class.vec, levels=c("Between", groups))

# Color palette
color_count <- length(groups) + 1
color_pal <- brewer.pal(n = 10, "Set1")
if (color_count <= 10) {
  get_palette = colorRampPalette(color_pal[1:color_count])
} else {
  get_palette = colorRampPalette(color_pal)
}

output_png <- file.path(opts$outdir, glue("{opts$prefix}.anosim.png"))
png(output_png, height=900, width=1000, res=100)
plot(anosim_results, col = get_palette(color_count),
     main=glue("Anosim: R={R}, P={P}"), xlab="Group", ylab=glue("{dist_type} Rank"))
dev.off()

output_pdf <- file.path(opts$outdir, glue("{opts$prefix}.anosim.pdf"))
pdf(output_pdf, height=9, width=10)
plot(anosim_results, col = get_palette(color_count),
     main=glue("Anosim: R={R}, P={P}"), xlab="Group", ylab=glue("{dist_type} Rank"))
dev.off()

output_xls <- file.path(opts$outdir, glue("{opts$prefix}.anosim.xls"))
output_df <- data.frame(list("Contrast"=opts$contrast, "R"=R, "P"=P))
write.table(output_df, 
            file=output_xls,
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE,
            sep = "\t")

