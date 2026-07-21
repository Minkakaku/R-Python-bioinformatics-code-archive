# R script to perform PCA, PCoA, NMDS analysis based on ASVs table (PCA) or beta diversity distance matrix (PCoA, NMDS)
# by Yewei on 2024.4.15
# Output: {opts$outdir}/{opts$contrast}.{opts$dist}.{png,pdf,xls}

library(optparse)
library(dplyr)
library(glue)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggforce)
library(doBy)

#############
# Arguments #
#############

# Parse arguments
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="ASVs table (PCA) or beta diversity distance matrix (PCoA, NMDS)"),
    make_option(c("-m", "--metadata"), type="character", default=NULL,
                help="Metadata file for all samples"),
    make_option(c("-c", "--contrast"), type="character", default=NULL,
                help="Contrast with groups to be included"),
    make_option(c("-d", "--dist"), type="character", default="braycurtis",
                help="Beta distance type for PCoA and NMDS (braycurtis, jaccard, weighted_unifrac, unweighted_unifrac) [default %default]"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="Group variable name in metadata [default %default]"),
    make_option(c("-s", "--sample"), type="character", default="sample",
                help="Sample variable name in metadata [default %default]"),
    make_option(c("-o", "--outdir"), type="character", default=".",
                help="Output directory [default %default]"),
    make_option(c("-t", "--type"), type="character", default="PCA",
                help="Analysis type (PCA, PCoA, NMDS) [default %default]"),
    make_option(c("-p", "--prefix"), type="character", default="output",
                help="Output file name prefix [default %default]"),
    make_option(c("-e", "--ellipse"), type="logical", action="store_true", default=FALSE,
                help="Add ellipse in PCA plot [default %default]"),
    make_option(c("-l", "--line"), type="logical", action="store_true", default=FALSE,
                help="Add line segments in PCA plot [default %default]"),
    make_option(c("-a", "--label"), type="logical", action="store_true", default=FALSE,
                help="Show sample labels in PCA plot [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

#############
# Functions #
#############

# Function to generate PCA and PCoA plot
# Input data.frame columns: "sample", "group", "PC1", "PC2"
pca_plot <- function(pca_df, group_list, x_var, y_var, x_label, y_label, p_title) {
  # ggplot2
  color_count <- length(unique(pca_df$group))
  color_pal <- brewer.pal(n = 10, "Set1")
  if (color_count <= 10) {
    get_palette = colorRampPalette(color_pal[1:color_count])
  } else {
    get_palette = colorRampPalette(color_pal)
  }

  # If any group has number of samples <= 3, use ggforce::geom_mark_ellipse() to draw ellipse otherwise use stat_ellipse()
  if (min(as.data.frame(table(pca_df$group))$Freq) > 3) {
    pca_mean <- summaryBy(PC1 + PC2 ~ group, pca_df, FUN = mean)
    pca_df <- merge(pca_df, pca_mean, by = 'group')
    pca_df$group <- factor(pca_df$group, levels = group_list)
    p_pca <- ggplot(pca_df, aes(PC1, PC2, colour=group)) +
      geom_point(aes(shape=group),
                 size=3)

    if (opts$label) {
      p_pca <- p_pca +
        geom_text(aes(label=sample),
                  hjust=0.5,
                  vjust=-1,
                  size=3,
                  show.legend=FALSE)
    }

    if (opts$ellipse) {
      p_pca <- p_pca +
        stat_ellipse(aes(fill=group, color=group),
                     type="norm",
                     geom="polygon",
                     level=0.95,
                     alpha=0.1,
                     linetype='dashed',
                     linewidth=1)
    }

    if (opts$line) {
      p_pca <- p_pca +
        geom_segment(aes(x = PC1.mean,
                         y = PC2.mean,
                         xend= PC1,
                         yend = PC2,
                         color = group),
                     alpha= 0.6,
                     show.legend = FALSE)
    }

    p_pca <- p_pca +
      scale_fill_manual(values=get_palette(color_count)) +
      scale_color_manual(values=get_palette(color_count)) +
      labs(title = p_title,
           x=paste0(x_label, " (", round(x_var, 2), "%)"),
           y=paste0(y_label, " (", round(y_var, 2), "%)"),
           colour="Group",
           shape="Group",
           fill="Group") +
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            legend.title = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5, face = "bold"))
  } else {
    pca_mean <- summaryBy(PC1 + PC2 ~ group, pca_df, FUN = mean)
    pca_df <- merge(pca_df, pca_mean, by = 'group')
    pca_df$group <- factor(pca_df$group, levels = group_list)
    p_pca <- ggplot(pca_df, aes(PC1, PC2, colour=group)) +
      geom_point(aes(shape=group),
                 size=3)

    if (opts$label) {
      p_pca <- p_pca +
        geom_text(aes(label=sample),
                  hjust=0.5,
                  vjust=-1,
                  size=3,
                  show.legend=FALSE)
    }

    if (opts$ellipse) {
      p_pca <- p_pca +
        ggforce::geom_mark_ellipse(aes(fill=group, color=group),
                                   type="norm",
                                   geom="polygon",
                                   level=0.95,
                                   alpha=0.1,
                                   linetype='dashed',
                                   linewidth=1)
    }

    if (opts$line) {
      p_pca <- p_pca +
        geom_segment(aes(x = PC1.mean,
                         y = PC2.mean,
                         xend= PC1,
                         yend = PC2,
                         color = group),
                     alpha= 0.6,
                     show.legend = FALSE)
    }

    p_pca <- p_pca +
      scale_x_continuous(expand = c(0.4, 0)) +
      scale_y_continuous(expand = c(0.4, 0)) +
      scale_fill_manual(values=get_palette(color_count)) +
      scale_color_manual(values=get_palette(color_count)) +
      labs(title = p_title,
           x=paste0(x_label, " (", round(x_var, 2), "%)"),
           y=paste0(y_label, " (", round(y_var, 2), "%)"),
           colour="Group",
           shape="Group",
           fill="Group") +
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            legend.title = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5, face = "bold"))
  }
  
  return(p_pca)
}


# Function to perform PCA analysis and generate a plot
pca_analysis <- function(abundance_df, meta_df, group_list) {
  abundance_pca_df <- abundance_df[, meta_df[meta_df$group %in% group_list,]$sample]

  # make all samples have the same total counts
  ratio_list <- min(colSums(abundance_pca_df)) / colSums(abundance_pca_df)
  for (i in 1:length(ratio_list)) {
    abundance_pca_df[, i] <- abundance_pca_df[, i] * ratio_list[i]
  }

  # select top80 ASVs with total counts
  abundance_pca_df <- abundance_pca_df[order(rowSums(abundance_pca_df),decreasing=T),][1:80,]

  # remove ASVs with all count value 0
  abundance_pca_df <- abundance_pca_df[rowSums(abundance_pca_df[])>0,]
  abundance_pca_df <- as.data.frame(t(abundance_pca_df))
  
  # remove ASVs with the same count value among all samples (a constant column cannot be rescaled to unit variance)
  abundance_pca_df <- abundance_pca_df[ , which(apply(abundance_pca_df, 2, var) != 0)]
  
  # merge group info
  abundance_pca_df$sample <- rownames(abundance_pca_df)
  abundance_pca_withGroup_df <- merge(x = abundance_pca_df, y = meta_df, by.x = "sample", by.y = "sample", all.x = TRUE)
  rownames(abundance_pca_withGroup_df) <- abundance_pca_withGroup_df$sample
  
  # perform pca analysis
  fpkm.pca <- prcomp(abundance_pca_withGroup_df[,!names(abundance_pca_withGroup_df) %in% colnames(meta_df)],
                     center = TRUE,
                     scale. = TRUE)
  
  # scale principal components
  df <- cbind(fpkm.pca$x[,1:2], abundance_pca_withGroup_df$group, abundance_pca_withGroup_df$sample) %>% 
    as.data.frame()
  df$PC1 <- as.numeric(df$PC1) / (fpkm.pca$sdev[1] * sqrt(nrow(abundance_pca_withGroup_df)))
  df$PC2 <- as.numeric(df$PC2) / (fpkm.pca$sdev[2] * sqrt(nrow(abundance_pca_withGroup_df)))
  df$V3 <- as.factor(df$V3)
  df <- df %>% 
    dplyr::rename(sample=V4,
                  group=V3)
  
  # calculate PC1 and PC2 importance
  x_val <- round(100*summary(fpkm.pca)$importance[2, 1], 2)
  y_val <- round(100*summary(fpkm.pca)$importance[2, 2], 2)
  
  # plot
  p_title <- "PCA"
  p_pca <- pca_plot(pca_df=df, 
                    group_list=group_list,
                    x_var=x_val, 
                    y_var=y_val, 
                    x_label="PC1", 
                    y_label="PC2",
                    p_title=p_title)
  
  # output
  return(list("p"=p_pca, "df"=df[,c("sample", "group", "PC1", "PC2")]))
}

# Function to perform PCoA analysis and generate a plot
pcoa_analysis <- function(distance_matrix_df, meta_df, group_list, dist) {
  library(vegan)

  # filter samples in contrast groups
  sample_list <- meta_df[meta_df$group %in% group_list,]$sample
  distance_matrix_df <- distance_matrix_df[sample_list, sample_list]
  
  # perform PCoA analysis
  pcoa <- cmdscale(distance_matrix_df, k=2, eig=TRUE)
  
  # generate data frame
  df <- as.data.frame(pcoa$points)
  colnames(df) <- c("PC1", "PC2")
  df$sample <- rownames(df)
  df <- merge(x=df, y=meta_df, by.x="sample", by.y="sample", all.x=TRUE)
  
  # calculate PC1 and PC2 importance
  var_exp <- (eigenvals(pcoa)/sum(eigenvals(pcoa)))[1:2]
  x_val <- as.numeric(var_exp[1]*100)
  y_val <- as.numeric(var_exp[2]*100)
  
  # plot
  p_title <- glue("PCoA - {dist}")
  p_pca <- pca_plot(pca_df=df, 
                    group_list=group_list,    
                    x_var=x_val, 
                    y_var=y_val, 
                    x_label="PCo1", 
                    y_label="PCo2",
                    p_title=p_title)
  
  # output
  df <- df %>% 
    rename(PCo1=PC1, PCo2=PC2)
  return(list("p"=p_pca, "df"=df[,c("sample", "group", "PCo1", "PCo2")]))
}

# perform NMDS analysis and generate a plot
nmds_analysis <- function(distance_matrix_df, meta_df, group_list, dist) {
  library(vegan)

  # filter samples in contrast groups
  sample_list <- meta_df[meta_df$group %in% group_list,]$sample
  distance_matrix_df <- distance_matrix_df[sample_list, sample_list]
  
  # perform NMDS analysis
  nmds <- metaMDS(distance_matrix_df, k=2)
  stress <- signif(nmds$stress, 4)
  
  # generate data frame
  df <- as.data.frame(nmds$points)
  df$sample <- rownames(df)
  df <- merge(x=df, y=meta_df, by.x="sample", by.y="sample", all.x=TRUE) %>%
    rename(NMDS1=MDS1, NMDS2=MDS2)
  
  # ggplot2
  color_count <- length(unique(df$group))
  color_pal <- brewer.pal(n = 10, "Set1")
  if (color_count <= 10) {
    get_palette = colorRampPalette(color_pal[1:color_count])
  } else {
    get_palette = colorRampPalette(color_pal)
  }

  p_title <- glue("NMDS - {dist}")
  
  if (stress < 0.001) {
    x_title <- "NMDS1 (Stress < 0.001)"
  } else {
    x_title <- glue("NMDS1 (Stress = {stress})")
  }

  # If any group has number of samples <= 3, use ggforce::geom_mark_ellipse() to draw ellipse otherwise use stat_ellipse()
  if (min(as.data.frame(table(df$group))$Freq) > 3) {
    nmds_mean <- summaryBy(NMDS1 + NMDS2 ~ group, df, FUN = mean)
    nmds_df <- merge(df, nmds_mean, by = 'group')
    nmds_df$group <- factor(nmds_df$group, levels = group_list)
    p_nmds <- ggplot(nmds_df, aes(NMDS1, NMDS2, colour=group)) +
        geom_point(aes(shape=group),
                   size=3)

    if (opts$label) {
      p_nmds <- p_nmds +
        geom_text(aes(label=sample),
                  hjust=0.5,
                  vjust=-1,
                  size=3,
                  show.legend=FALSE)
    }

    if (opts$ellipse) {
      p_nmds <- p_nmds +
        stat_ellipse(aes(fill=group, color=group),
                     type="norm",
                     geom="polygon",
                     level=0.95,
                     alpha=0.1,
                     linetype='dashed',
                     linewidth=1)
    }

    if (opts$line) {
      p_nmds <- p_nmds +
        geom_segment(aes(x = NMDS1.mean,
                         y = NMDS2.mean,
                         xend= NMDS1,
                         yend = NMDS2,
                         color = group),
                     alpha= 0.6,
                     show.legend = FALSE)
    }

    p_nmds <- p_nmds +
      scale_fill_manual(values=get_palette(color_count)) +
      scale_color_manual(values=get_palette(color_count)) +
      labs(title=p_title,
           x=x_title,
           y="NMDS2",
           colour="Group",
           shape="Group",
           fill="Group") +
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            legend.title = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5, face = "bold"))
  } else {
    nmds_mean <- summaryBy(NMDS1 + NMDS2 ~ group, df, FUN = mean)
    nmds_df <- merge(df, nmds_mean, by = 'group')
    nmds_df$group <- factor(nmds_df$group, levels = group_list)
    p_nmds <- ggplot(nmds_df, aes(NMDS1, NMDS2, colour=group)) +
      geom_point(aes(shape=group),
                 size=3)

    if (opts$label) {
      p_nmds <- p_nmds +
        geom_text(aes(label=sample),
                  hjust=0.5,
                  vjust=-1,
                  size=3,
                  show.legend=FALSE)
    }

    if (opts$ellipse) {
      p_nmds <- p_nmds +
        ggforce::geom_mark_ellipse(aes(fill=group, color=group),
                                   type="norm",
                                   geom="polygon",
                                   level=0.95,
                                   alpha=0.1,
                                   linetype='dashed',
                                   linewidth=1)
    }

    if (opts$line) {
      p_nmds <- p_nmds +
        geom_segment(aes(x = NMDS1.mean,
                         y = NMDS2.mean,
                         xend= NMDS1,
                         yend = NMDS2,
                         color = group),
                     alpha= 0.6,
                     show.legend = FALSE)
    }

    p_nmds <- p_nmds +
      scale_x_continuous(expand = c(0.4, 0)) +
      scale_y_continuous(expand = c(0.4, 0)) +
      scale_fill_manual(values=get_palette(color_count)) +
      scale_color_manual(values=get_palette(color_count)) +
      labs(title=p_title,
           x=x_title,
           y="NMDS2",
           colour="Group",
           shape="Group",
           fill="Group") +
      theme_bw() +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            legend.title = element_text(face = "bold"),
            plot.title = element_text(hjust = 0.5, face = "bold"))
  }
  
  # output
  return(list("p"=p_nmds, "df"=df[,c("sample", "group", "NMDS1", "NMDS2")]))
}

########
# Main #
########

# Convert contrast to group list
group_list <- unlist(strsplit(opts$contrast, "-vs-"))

# Read beta diversity distance matrix (PCoA, NMDS) or ASVs abundance table (PCA)
beta_matrix_df <- read.table(opts$input, header=T, check.names=F, row.names=1, sep="\t", quote="")

# Read metadata
meta_df <- read.table(opts$meta, header=T, check.names=F, row.names=NULL, sep="\t") %>% 
  dplyr::rename("sample" = opts$sample, "group" = opts$group) %>%
  dplyr::select(sample, group)

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
          
# Perform PCA, PCoA, NMDS analysis
if (tolower(opts$type) == "pca") {
  output_list <- pca_analysis(abundance_df=beta_matrix_df, meta_df=meta_df, group_list=group_list)
  output_df <- output_list$df
  output_p <- output_list$p
  analysis_type <- "PCA"
} else if (tolower(opts$type) == "pcoa") {
  output_list <- pcoa_analysis(distance_matrix_df=beta_matrix_df, meta_df=meta_df, group_list=group_list, dist=dist_type)
  output_df <- output_list$df
  output_p <- output_list$p
  analysis_type <- "PCoA"
} else if (tolower(opts$type) == "nmds") {
  output_list <- nmds_analysis(distance_matrix_df=beta_matrix_df, meta_df=meta_df, group_list=group_list, dist=dist_type)
  output_df <- output_list$df
  output_p <- output_list$p
  analysis_type <- "NMDS"
}

# Output
write.table(output_df, 
            file=file.path(opts$outdir, glue("{opts$prefix}.{analysis_type}.xls")), 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE,
            sep = "\t")
ggsave(filename = file.path(opts$outdir, glue("{opts$prefix}.{analysis_type}.png")), 
       device = "png",
       plot = output_p,
       width = 10,
       height = 8)
ggsave(filename = file.path(opts$outdir, glue("{opts$prefix}.{analysis_type}.pdf")), 
       device = "pdf",
       plot = output_p,
       width = 10,
       height = 8)