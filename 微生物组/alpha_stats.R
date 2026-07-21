# R script to generate alpha index boxplot and perform pairwise t-test
# Input: alpha index merged table
# By Yewei on 2024.4.15

# Log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Library
library(reshape2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(glue)
library(ggpubr)
library(scales)

# Input
alpha_merged_df <- read.table(snakemake@input[["xls"]], check.names=F, header=T, sep="\t")
contrast <- snakemake@wildcards[["contrast"]]
group_list <- unlist(strsplit(contrast, "-vs-"))
meta_df <- read.table(snakemake@params[["meta"]], check.names=F, header=T, sep="\t") %>%
  dplyr::select(SampleID, Group) %>%
  dplyr::rename(sample=SampleID, group=Group) %>%
  dplyr::filter(group %in% group_list)
alpha_index_list <- c("observed_features", "ace", "chao1", "simpson", "shannon", "faith_pd")
alpha_name_list <- c("Observed Features", "ACE Index", "Chao1 Index", "Simpson Index", "Shannon Index", "Faith's PD Index")

# Main
# Sample-level data
sample_data_df <- alpha_merged_df[alpha_merged_df$sample %in% meta_df$sample,]
write.table(sample_data_df,
            file=snakemake@output[["sample_data"]],
            row.names=FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Group-level data
group_data_df <- merge(alpha_merged_df, meta_df, by="sample", all.y=TRUE) %>%
  dplyr::select(-sample) %>%
  dplyr::group_by(group) %>%
  summarise(observed_features_mean=mean(observed_features),
            observed_features_stderr=sd(observed_features)/sqrt(n()),
            ace_mean=mean(ace),
            ace_stderr=sd(ace)/sqrt(n()),
            chao1_mean=mean(chao1),
            chao1_stderr=sd(chao1)/sqrt(n()),
            simpson_mean=mean(simpson),
            simpson_stderr=sd(simpson)/sqrt(n()),
            shannon_mean=mean(shannon),
            shannon_stderr=sd(shannon)/sqrt(n()),
            goods_coverages_mean=mean(goods_coverage),
            goods_coverage_stderr=sd(goods_coverage)/sqrt(n()),
            faith_pd_mean=mean(faith_pd),
            faith_pd_stderr=sd(faith_pd)/sqrt(n()))
write.table(group_data_df,
            file=snakemake@output[["group_data"]],
            row.names=FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Perform t-test
test_data_df <- merge(alpha_merged_df, meta_df, by="sample", all.y=TRUE)
for (alpha in alpha_index_list) {
  test_result_df <- as.data.frame(compare_means(as.formula(glue("{alpha}~group")),
                                                data = test_data_df,
                                                method = snakemake@params[["table_pairwise_test_method"]],
                                                p.adjust.method = snakemake@params[["table_pairwise_test_p_adjust_method"]]))[,c("group1", "group2", "p.adj")]
  colnames(test_result_df)[3] <- alpha
  if (exists("test_output_df")) {
    test_output_df <-  merge(test_output_df, test_result_df, by=c("group1", "group2"))
  } else {
    test_output_df <- test_result_df
  }
}
write.table(test_output_df,
            file=snakemake@output[["test_p"]],
            row.names=FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Boxplot
for (alpha in alpha_index_list) {
  boxplot_df <- merge(alpha_merged_df, meta_df, by = "sample", all.y=TRUE) %>%
    dplyr::select(sample, group, alpha) %>%
    dplyr::rename(value=alpha)

  n_group <- length(unique(boxplot_df$group))
  color_pal <- brewer.pal(n = 5, "Set2")
  if (n_group <= 5) {
    get_palette <- colorRampPalette(color_pal[1:n_group])
  } else {
    get_palette <- colorRampPalette(color_pal)
  }

  is_outlier <- function(x) {
    return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
  }

  boxplot_df <- boxplot_df %>%
    mutate(outlier = sample) %>%
    group_by(group) %>%
    mutate(is_outlier=ifelse(is_outlier(value), value, as.numeric(NA)))
  boxplot_df$outlier[which(is.na(boxplot_df$is_outlier))] <- as.numeric(NA)
  boxplot_df$remove_outlier <- ifelse(is.na(boxplot_df$outlier), boxplot_df$value, NA)

  boxplot_df$group <- factor(boxplot_df$group, levels=group_list)

  alpha_name <- alpha_name_list[match(alpha, alpha_index_list)]
  p_boxplot <- ggplot(boxplot_df, aes(x=group, y=value)) +
    stat_boxplot(geom = "errorbar", width=0.3, aes(group=group)) +
    geom_boxplot(aes(fill=group), color="black", width=0.5) +
    geom_text(aes(label=outlier, color=group), na.rm=TRUE, hjust=1.2, size=3) +
    geom_jitter(aes(x=group, y=remove_outlier), height = 0, width=0.2, color="black", size=1.5) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Set2", guide=NULL) +
    scale_y_continuous(labels = label_comma()) +
    labs(x = "Group",
         y = alpha_name,
         title = glue("Boxplot of {alpha_name}: {contrast}"),
         fill = "Group") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, vjust = 1),
          legend.title = element_text(face = "bold"))

  # Whther to show pairwise test results in the boxplot
  if (snakemake@params[["plot_pairwise_test"]] == "True") {
    hide_ns_flag <- ifelse(snakemake@params[["plot_pairwise_test_hideNS"]] == "True", TRUE, FALSE)
    p_boxplot <- p_boxplot +
      geom_pwc(method = snakemake@params[["plot_pairwise_test_method"]],
               p.adjust.method = snakemake@params[["plot_pairwise_test_p_adjust_method"]],
               hide.ns = hide_ns_flag,
               label = "{p.adj.format} {p.adj.signif}",
               tip.length = 0.01)
  }

  # Whther to show global test results in the boxplot
  if (snakemake@params[["plot_global_test"]] == "True" && length(group_list) > 2) {
    p_boxplot <- p_boxplot +
      stat_compare_means(method = snakemake@params[["plot_global_test_method"]],
                         label.x = -Inf,
                         hjust = -0.2,
                         label.y = Inf,
                         vjust = 1.5,
                         size = 4)
  }

  outdir <- dirname(snakemake@output[["boxplot_png"]][1])
  ggsave(filename = glue("{outdir}/{contrast}.{alpha}.boxplot.png"),
         device = "png",
         plot = p_boxplot,
         width = 10,
         height = 9)
  ggsave(filename = glue("{outdir}/{contrast}.{alpha}.boxplot.pdf"),
         device = "pdf",
         plot = p_boxplot,
         width = 10,
         height = 9)
}

