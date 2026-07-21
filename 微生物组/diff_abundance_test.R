# R script to perform Differential Abundance Test on each taxonomy's abundance among contrast groups
# by Yewei on 2024.5.10
# Output: {opts$outdir}/{opts$contrast}.{opts$level}.{opts$type}.{png,pdf,xls}

library(optparse)
library(dplyr)
library(glue)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggpubr)

#############
# Arguments #
#############

# Parse arguments
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Matrix file of relative taxonomic abundace for all samples"),
    make_option(c("-m", "--metadata"), type="character", default=NULL,
                help="Metadata file for all samples"),
    make_option(c("-c", "--contrast"), type="character", default=NULL,
                help="Contrast with groups to be included"),
    make_option(c("-n", "--group"), type="character", default="group",
                help="Group variable name in metadata [default %default]"),
    make_option(c("-s", "--sample"), type="character", default="sample",
                help="Sample variable name in metadata [default %default]"),
    make_option(c("-l", "--level"), type="character", default="Species",
                help="Taxonomic level (Kingdom, Phylum, Class, Order, Family, Genus, Species) [default %default]"),
    make_option(c("-o", "--outdir"), type="character", default=".",
                help="Output directory [default %default]"),
    make_option(c("-t", "--type"), type="character", default="rank_sum",
                help="Differential abundance test type (rank_sum or anova) [default %default]"),
    make_option(c("-b", "--bar_width"), type="numeric", default=0.5,
                help="Bar width barplot [default %default]")
  )
  opts = parse_args(OptionParser(option_list=option_list))
}

#############
# Functions #
#############

# Function to perform test using ggpubr::compare_means()
# Input:
#   1. abundance_df (opts$level, samples);
#   2. taxonomy_level_df (opts$level, taxonomy concatenated by " | " for duplicate opts$level)
#   3. meta_df (sample, group)
#   4. test_type ("rank_sum" or "anova")
#   5. group list
perform_test <- function(abundance_df, taxonomy_level_df, meta_df, test_type, group_list) {
  # Convert wide format to long format
  abundance_long_df <- melt(abundance_df, id.vars = "level") %>%
    dplyr::rename(sample=variable, abundance=value)

  # Add group info
  abundance_long_test_df <- merge(abundance_long_df, meta_df, by.x = "sample", by.y="sample")

  # Calculate mean and sd for each group
  mean_sd_df <- abundance_long_test_df %>%
    dplyr::select(level, group, abundance) %>%
    dplyr::group_by(level, group) %>%
    dplyr::summarise(mean = mean(abundance), sd = sd(abundance))

  # perform test
  if (test_type == "anova") {
    test <- "anova"
  } else if (test_type == "rank_sum" & length(group_list) == 2) {
    test <- "wilcox.test"
  } else if (test_type == "rank_sum" & length(group_list) > 2) {
    test <- "kruskal.test"
  }
  test_results_df <- compare_means(abundance ~ group,
                                   data = abundance_long_test_df,
                                   group.by = "level",
                                   method = test,
                                   p.adjust.method = "fdr")

  # Add mean and sd into results dataframe
  test_results_output_df <- test_results_df
  var_list <- c()
  for (g in group_list) {
    mean_sd_df_tmp <- mean_sd_df %>%
      dplyr::filter(group == g) %>%
      dplyr::select(level, mean, sd)
    colnames(mean_sd_df_tmp)[2:3] <- c(glue("mean_{g}"), glue("sd_{g}"))
    var_list <- c(var_list, c(glue("mean_{g}"), glue("sd_{g}")))
    test_results_output_df <- merge(test_results_output_df, mean_sd_df_tmp, by.x="level", by.y="level")
  }

  # Add taxonomy into results dataframe
  test_results_output_df <- merge(test_results_output_df, taxonomy_level_df, by.x="level", by.y="level")

  # Prepare results dataframe for output
  test_results_output_df <- test_results_output_df %>%
    dplyr::select(level, Taxonomy, var_list, p, p.adj, method) %>%
    dplyr::mutate(mean_max = max(select(., starts_with("mean")))) %>%
    dplyr::arrange(p, desc(mean_max)) %>%
    dplyr::select(-mean_max)

  return(list("test_results_output_df"=test_results_output_df,
              "abundance_long_test_df"=abundance_long_test_df))
}

# Function to generate barplot based on test results dataframe from perform_test()
generate_barplot<- function(test_results_output_df, abundance_long_test_df, group_list, level_name) {
  # Remove rows without p value (NaN)
  test_results_output_df <- test_results_output_df[!is.na(test_results_output_df$p),]
  # Select levels with p < 0.05 (top20 if more than 20)
  if (sum(test_results_output_df$p < 0.05) > 0) {
    if (sum(test_results_output_df$p < 0.05) > 20) {
      top_list <- test_results_output_df[test_results_output_df$p < 0.05, "level"][1:20]
    } else {
      top_list <- test_results_output_df[test_results_output_df$p < 0.05, "level"]
    }
  } else {
    # If no levels with p < 0.05, just select top5
    top_list <- test_results_output_df[1:5, "level"]
  }

  plot_df <- abundance_long_test_df[abundance_long_test_df$level %in% top_list,] %>%
  dplyr::group_by(level, group) %>%
  dplyr::summarise(n = n(), mean = mean(abundance), sd = sd(abundance)) %>%
  dplyr::mutate(se = sd/sqrt(n))

  order_df <- plot_df %>%
    dplyr::group_by(level) %>%
    dplyr::summarise(mean_max = max(mean)) %>%
    dplyr::arrange(desc(mean_max))

  plot_df$level <- factor(plot_df$level, levels=order_df$level)
  plot_df$group <- factor(plot_df$group, levels=group_list)

  # Color palette
  color_count <- length(group_list)
  color_pal <- brewer.pal(n = 10, "Set1")
  if (color_count <= 10) {
    get_palette = colorRampPalette(color_pal[1:color_count])
  } else {
    get_palette = colorRampPalette(color_pal)
  }

  p <- ggplot(plot_df, aes(level, mean, fill = group)) +
    geom_bar(position = "dodge", stat = "identity", width = opts$bar_width) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(opts$bar_width)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
          axis.text = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.margin = margin(, , , 3, "cm")) +
    labs(x = level_name,
         y = "Relative Abundance",
         fill = "Group") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    scale_fill_manual(values=get_palette(color_count))

  return(p)
}

########
# Main #
########

# Read meta data
meta_df <- read.table(opts$metadata, check.names=F, header=T, sep="\t") %>%
  dplyr::select(opts$sample, opts$group) %>%
  dplyr::rename(sample = opts$sample, group = opts$group)

# Read contrast groups
group_list <- unlist(strsplit(opts$contrast, "-vs-"))
meta_df <- meta_df %>%
  dplyr::filter(group %in% group_list)

# Read abundance matrix (sum up for duplicate taxonomy names)
abundance_df <- read.table(opts$input, header=T, sep="\t", check.names=F, quote="") %>%
  dplyr::rename(level = opts$level)
taxonomy_level_df <- abundance_df %>%
  dplyr::select(Taxonomy, level) %>%
  dplyr::group_by(level) %>%
  dplyr::summarise(Taxonomy = paste(Taxonomy, collapse=" | "))
abundance_df <- abundance_df %>%
  dplyr::select(level, meta_df$sample) %>%
  dplyr::group_by(across(all_of("level"))) %>%
  dplyr::summarise_all(sum)

# Perform test
output_list <- perform_test(abundance_df, taxonomy_level_df, meta_df, opts$type, group_list)
test_results_output_df <- output_list$test_results_output_df
abundance_long_test_df <- output_list$abundance_long_test_df

# Generate barplot
p <- generate_barplot(test_results_output_df, abundance_long_test_df, group_list, opts$level)

# Output
colnames(test_results_output_df)[1] <- opts$level
write.table(test_results_output_df,
            file=file.path(opts$outdir, glue("{opts$contrast}.{opts$level}.{opts$type}.xls")),
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE,
            sep = "\t")
ggsave(filename = file.path(opts$outdir, glue("{opts$contrast}.{opts$level}.{opts$type}.png")),
       device = "png",
       plot = p,
       width = 14,
       height = 10)
ggsave(filename = file.path(opts$outdir, glue("{opts$contrast}.{opts$level}.{opts$type}.pdf")),
       device = "pdf",
       plot = p,
       width = 14,
       height = 10)