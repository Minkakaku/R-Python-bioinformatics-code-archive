# R script to
#   1. Calculate spearman correlation p-value and coefficient among Genus-level taxonomy
#   2. Generate association network plot using R::igraph

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(igraph)
library(psych)
library(dplyr)
library(tidyverse)
library(RColorBrewer)

#########
# Input #
#########

# Read Genus-level relative abundance table
reabund_df <- read.table(snakemake@input[["reabundance"]], head=T, check.names=F, row.names=NULL, sep="\t") %>%
  dplyr::mutate(Phylum = unlist(lapply(Taxonomy, function(x) unlist(strsplit(x, ";"))[2]))) %>%
  dplyr::select(-Taxonomy) %>%
  dplyr::group_by(Genus, Phylum) %>%
  dplyr::summarise_all(sum) %>%
  dplyr::ungroup()

# Read metadata
meta_df <- read.table(snakemake@params[["meta"]], header=T, check.names=F, row.names=NULL, sep="\t")

# Filter samples in group
if (snakemake@wildcards[["group"]] != "all") {
  reabund_df <- reabund_df %>%
    dplyr::select(Genus, Phylum, meta_df[meta_df$Group == snakemake@wildcards, "SampleID"])
}

# Spearman correlation coefficient and p-value threshold in association network analysis
p_cutoff <- snakemake@params[["p_cutoff"]]
r_cutoff <- snakemake@params[["r_cutoff"]]

##################################
# Calculate Spearman correlation #
##################################

# Calculate spearman correlation p-value and coefficient
reabund_corr_df <- reabund_df %>%
  dplyr::select(-Phylum) %>%
  column_to_rownames(var = "Genus")
# row: Genus; column: sample
reabund_corr_df <-  as.data.frame(t(reabund_corr_df))
# remove genus with all count equals zero
reabund_corr_df <- reabund_corr_df[,colSums(reabund_corr_df) != 0]
# Select top50 Genus with total relative abundance among all samples
if (ncol(reabund_corr_df) < 50) {
  reabund_corr_df <- reabund_corr_df[,order(colSums(reabund_corr_df), decreasing = TRUE)]
} else {
  reabund_corr_df <- reabund_corr_df[,order(colSums(reabund_corr_df), decreasing = TRUE)[1:50]]
}
# Perform calculation using corr.test in R::psych
if (nrow(reabund_corr_df) <= 3) {
  corr_method <- "pearson"
} else {
  corr_method <- "spearman"
}
occor <- corr.test(reabund_corr_df, use="pairwise", method=corr_method, adjust="fdr", alpha=.05)
occor.r <- occor$r # r: spearman coefficnet
occor.p <- occor$p # p: p-value

# Output Spearman correlation p-value and coefficient matrix (top50 Genus)
write.table(occor.p,
            file = snakemake@output[["p_table"]],
            row.names = TRUE,
            col.names = NA,
            quote = FALSE,
            sep = "\t")
write.table(occor.r,
            file = snakemake@output[["r_table"]],
            row.names = TRUE,
            col.names = NA,
            quote = FALSE,
            sep = "\t")

############################
# Plot association network #
############################

# Filter out association using some threshold (not present in network plot)
occor.r[occor.p >= p_cutoff | abs(occor.r) <= r_cutoff] <- 0

# Construct igraph object
igraph <- graph_from_adjacency_matrix(occor.r, mode="undirected", weighted=TRUE, diag=FALSE)

# Remove isolated nodes
bad.vs <- V(igraph)[degree(igraph) == 0]
igraph <- delete.vertices(igraph, bad.vs)

# Save igraph$weight to a variable
igraph.weight <- E(igraph)$weight
# Remove igraph$weight before plotting
E(igraph)$weight <- NA

# Set edge color
# Postive correlation: red
# Negative correlation: blue
#c_scale <-  colorRamp(c('#3100ff','#7d5efc','#ff8080','#ff0000'))
c_scale <-  colorRamp(c('#3100ff','white','#ff0000'))
E(igraph)$color <- apply(c_scale((igraph.weight+1)/2), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )

# Set edge width by correlation coefficient
E(igraph)$width <- abs(igraph.weight)*6

# Add Genus relative abundance and Phylum info
# Prepare a dataframe for downstream analysis: Genus, Phylum, Total_reabundance, color
attr_df <- reabund_df %>%
  dplyr::select(-Genus, -Phylum)
attr_df <- as.data.frame(cbind(reabund_df$Genus, reabund_df$Phylum, rowSums(attr_df)/ncol(attr_df)))
colnames(attr_df) <- c("Genus", "Phylum", "abundance_mean")
attr_df$abundance_mean <- as.numeric(attr_df$abundance_mean)
rownames(attr_df) <- attr_df$Genus
attr_df <- attr_df[V(igraph)$name,]
attr_df <- attr_df[order(attr_df$abundance_mean, decreasing = TRUE),]
Phylum_order <- attr_df %>%
  dplyr::group_by(Phylum) %>%
  dplyr::summarise(sum = sum(abundance_mean)) %>%
  dplyr::arrange(desc(sum)) %>%
  dplyr::select(Phylum)
# Color palette
color_count <- length(unique(Phylum_order[[1]]))
color_pal <- brewer.pal(n = 12, "Set3")
if (color_count <= 12) {
  get_palette <- colorRampPalette(color_pal[1:color_count])
} else {
  get_palette <- colorRampPalette(color_pal)
}
colors <- get_palette(color_count)
Phylum_attr <- cbind(Phylum_order, colors)
attr_df <- merge(attr_df, Phylum_attr, by="Phylum")
rownames(attr_df) <- attr_df$Genus

# Set vertices size using relative abundance
rescale <- function(x,a,b,c,d){c + (x-a)/(b-a)*(d-c)}
igraph.size <- attr_df[V(igraph)$name, "abundance_mean"]
igraph.size1 <- rescale(igraph.size, min(igraph.size), max(igraph.size), 10, 20)
V(igraph)$size <- igraph.size1

# Set vertices color
igraph.col <- attr_df[V(igraph)$name, ]
levels(igraph.col$Phylum) <- Phylum_order[[1]] # 直接修改levles可以连值全部对应替换
V(igraph)$color <- attr_df[V(igraph)$name, "colors"]

# Prepare legend (Genus and Phylum info)
legend_df <- attr_df[V(igraph)$name, ]
legend_df$num <- 1:nrow(legend_df)
legend_df$label_show <- paste0(legend_df$num, ": ", legend_df$Genus)

# Output plot
pdf(file=snakemake@output[["plot_pdf"]], width=20, height=10)
par(mfrow=c(1,2))
set.seed(123)
p <- plot(igraph,
          main=NA,
          vertex.frame.color=NA,
          vertex.label=1:length(V(igraph)$name),
          vertex.label.color= "black",
          layout = layout_with_kk,
          edge.lty=1,
          edge.curved=TRUE,
          margin=c(0,0,0,0))
plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=0:1, ylim=0:1)
legend("topleft",
       title="Phylum",
       title.adj=c(0,0),
       title.cex=1.5,
       legend=Phylum_attr$Phylum,
       pt.cex=1.2,
       pch=16,
       col=Phylum_attr$colors,
       cex=0.9,
       bty="n",
       ncol=4)
legend("bottomleft",
       title="Genus",
       title.adj=c(0,0),
       title.cex=1.5,
       legend=legend_df$label_show,
       pt.cex=1.2,
       pch=16,
       col=legend_df$colors,
       text.width=0.5,
       bty="n",
       ncol=2,
       cex=0.9)
dev.off()

png(snakemake@output[["plot_png"]], width=2000, height=1000, res=100)
par(mfrow=c(1,2))
set.seed(123)
p <- plot(igraph,
          main=NA,
          vertex.frame.color=NA,
          vertex.label=1:length(V(igraph)$name),
          vertex.label.color= "black",
          layout = layout_with_kk,
          edge.lty=1,
          edge.curved=TRUE,
          margin=c(0,0,0,0))
plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=0:1, ylim=0:1)
legend("topleft",
       title="Phylum",
       title.adj=c(0,0),
       title.cex=1.5,
       legend=Phylum_attr$Phylum,
       pt.cex=1.2,
       pch=16,
       col=Phylum_attr$colors,
       cex=0.9,
       bty="n",
       ncol=4)
legend("bottomleft",
       title="Genus",
       title.adj=c(0,0),
       title.cex=1.5,
       legend=legend_df$label_show,
       pt.cex=1.2,
       pch=16,
       col=legend_df$colors,
       text.width=0.5,
       bty="n",
       ncol=2,
       cex=0.9)
dev.off()
