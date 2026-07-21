log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)

# Input
featureCounts_count <- snakemake@input[["counts"]]
sample_list <- colnames(read.table(featureCounts_count, header=TRUE, sep="\t", check.names=FALSE))[-1]
featureCounts_dir <- dirname(featureCounts_count)
metadata <- read.table(snakemake@params[["meta"]], header=TRUE, sep="\t")

# Output
fpkm_count <- snakemake@output[["fpkm"]]
boxplot_png <- snakemake@output[["boxplot"]][1]
boxplot_pdf <- snakemake@output[["boxplot"]][2]
heatmap_png <- snakemake@output[["heatmap"]][1]
heatmap_pdf <- snakemake@output[["heatmap"]][2]
pca_png <- snakemake@output[["pca"]][1]
pca_pdf <- snakemake@output[["pca"]][2]

# Analysis
# 1. Get the length of all genes
n <- 0
for (s in sample_list) {
  featureCounts <- file.path(featureCounts_dir, s, paste0(s, "_count.tsv"))
  
  featureCounts_df <- read.table(featureCounts, sep="\t", header=TRUE, comment.char = "#")
  
  if (n == 0) {
    geneLength <- featureCounts_df %>%
      dplyr::select(Geneid, Length)
    n <- 1
  } else {
    geneLength_tmp <- featureCounts_df %>%
      dplyr::select(Geneid, Length)
    geneLength <- rbind(geneLength, geneLength_tmp) %>%
      dplyr::distinct(Geneid, .keep_all = TRUE)
  }
}

# 2. Convert featureCounts to FPKM
featureCounts_merged_df <- read.table(featureCounts_count, sep="\t", header=TRUE, check.names=FALSE)
featureCounts_merged_withLength_df <- merge(x = featureCounts_merged_df, 
                                            y = geneLength, 
                                            by.x = "Gene",
                                            by.y = "Geneid",
                                            all.x = TRUE)
Gene <- featureCounts_merged_withLength_df$Gene
length <- featureCounts_merged_withLength_df$Length
kb <- length / 1000
counts <- featureCounts_merged_withLength_df[,2:(ncol(featureCounts_merged_withLength_df)-1)]
rpk <- counts / kb
fpkm <- as.data.frame(t(t(rpk)/colSums(counts) * 10^6))
fpkm_df <- cbind(Gene, fpkm)
write.table(fpkm_df, 
            file = fpkm_count,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

# 3. FPKM distribution Boxplot
plot_df <- melt(fpkm_df, 
                id.vars = "Gene") %>%
  dplyr::mutate(log_value = log10(value)) %>% 
  dplyr::rename(sample = variable)
plot_df <- merge(x = plot_df,
                 y = metadata,
                 by.x = "sample",
                 by.y = "sample",
                 all.x = TRUE)
plot_df <- plot_df
n_group <- length(unique(plot_df$group))

color_pal <- brewer.pal(n = 5, "Set1")
if (n_group <= 5) {
  get_palette = colorRampPalette(color_pal[1:n_group])
} else {
  get_palette = colorRampPalette(color_pal)
}

p_boxplot <- ggplot(plot_df, aes(x=sample, y=log_value, fill=group, color=group, alpha=group)) + 
  geom_boxplot() +
  scale_alpha_manual(values = rep(0.5, n_group)) + 
  scale_fill_manual(values = get_palette(n_group)) +
  scale_color_manual(values = get_palette(n_group)) +
  labs(x="Sample", y = "log10(FPKM)") + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, vjust = 1))

ggsave(filename = boxplot_png,
       plot = p_boxplot,
       device = "png",
       height = 8,
       width = 10)
ggsave(filename = boxplot_pdf,
       plot = p_boxplot,
       device = "pdf",
       height = 8,
       width = 10)

# 4. Correlation heatmap
fpkm_cor_df <- fpkm_df
rownames(fpkm_cor_df) <- fpkm_cor_df$Gene
fpkm_cor_df <- fpkm_cor_df %>% 
  dplyr::select(-Gene)
fpkm_cor_df <- fpkm_cor_df[rowSums(fpkm_cor_df[])>0,]
fpkm_cor_df <- log10(fpkm_cor_df+1)

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cormat <- reorder_cormat(cor(fpkm_cor_df))
melted_cormat <- melt(cormat, na.rm = TRUE) %>%
  mutate(value = value*value)
melted_cormat$value <- round(melted_cormat$value, 3)
min_cor <- min(melted_cormat$value)
max_cor <- max(melted_cormat$value)

if (length(sample_list) <= 12){
  p_heatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = (min_cor+max_cor)/2, 
                         limit = c(min_cor, max_cor), 
                         space = "Lab", 
                         name = expression(~R^2)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1),
          axis.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(face = "bold")) +
    labs(x = "Sample", 
         y = "Sample") + 
    coord_fixed()
} else {
  p_heatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    #geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = (min_cor+max_cor)/2, 
                         limit = c(min_cor, max_cor), 
                         space = "Lab", 
                         name = expression(~R^2)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1),
          axis.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(face = "bold")) +
    labs(x = "Sample", 
         y = "Sample") + 
    coord_fixed()
}


ggsave(filename = heatmap_png,
       plot = p_heatmap,
       device = "png",
       height = 8,
       width = 10)
ggsave(filename = heatmap_pdf,
       plot = p_heatmap,
       device = "pdf",
       height = 8,
       width = 10)

# 5. PCA
fpkm_pca_df <- fpkm_df
rownames(fpkm_pca_df) <- fpkm_pca_df$Gene
fpkm_pca_df <- fpkm_pca_df %>% 
  dplyr::select(-Gene)

# remove genes with all count value 0
fpkm_pca_df <- fpkm_pca_df[rowSums(fpkm_pca_df[])>0,]
fpkm_pca_df <- as.data.frame(t(fpkm_pca_df))

# remove genes with the same count value among all samples (a constant column cannot be rescaled to unit variance)
fpkm_pca_df <- fpkm_pca_df[ , which(apply(fpkm_pca_df, 2, var) != 0)]


fpkm_pca_df$sample <- rownames(fpkm_pca_df)
fpkm_pca_withGroup_df <- merge(x = fpkm_pca_df, 
                               y = metadata, 
                               by.x = "sample",
                               by.y = "sample",
                               all.x = TRUE)
rownames(fpkm_pca_withGroup_df) <- fpkm_pca_withGroup_df$sample

colnames(metadata)
fpkm.pca <- prcomp(fpkm_pca_withGroup_df[,!names(fpkm_pca_withGroup_df) %in% colnames(metadata)],
                   center = TRUE,
                   scale. = TRUE)

#autoplot(fpkm.pca, data = fpkm_pca_withGroup_df, colour = "group")

# The principal components are scaled in the autoplot method
# Manually scale them if ggplot2 is used
df <- cbind(fpkm.pca$x[,1:2], fpkm_pca_withGroup_df$group, fpkm_pca_withGroup_df$sample) %>% 
  as.data.frame()
df$PC1 <- as.numeric(df$PC1) / (fpkm.pca$sdev[1] * sqrt(nrow(fpkm_pca_withGroup_df)))
df$PC2 <- as.numeric(df$PC2) / (fpkm.pca$sdev[2] * sqrt(nrow(fpkm_pca_withGroup_df)))
df$V3 <- as.factor(df$V3)

# ggplot2
color_count <- length(unique(df$V3))
color_pal <- brewer.pal(n = 5, "Set1")
if (color_count <= 5) {
  get_palette = colorRampPalette(color_pal[1:color_count])
} else {
  get_palette = colorRampPalette(color_pal)
}

p_pca <- ggplot(df, aes(PC1, PC2, colour = V3)) +
  geom_point(aes(shape = V3), 
             size = 3) + 
  geom_text(aes(label = V4),
            hjust = 0.5,
            vjust = -1,
            size = 3) + 
  scale_color_manual(values = get_palette(color_count)) + 
  labs(x = paste0("PC1 (", round(100*summary(fpkm.pca)$importance[2, 1], 2), "%)"),
       y = paste0("PC2 (", round(100*summary(fpkm.pca)$importance[2, 2], 2), "%)"),
       colour = "Group",
       shape = "Group") + 
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        axis.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

ggsave(filename = pca_png,
       plot = p_pca,
       device = "png",
       height = 8,
       width = 10)
ggsave(filename = pca_pdf,
       plot = p_pca,
       device = "pdf",
       height = 8,
       width = 10)
