# R script to generate stats tables and plots based on FUNGuild results

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(RColorBrewer)

#########
# Input #
#########

# Read meta data
meta_df <- read.table(snakemake@params[["meta"]], check.names=F, header=T, sep="\t") %>%
  dplyr::select("SampleID", "Group") %>%
  dplyr::rename(sample = "SampleID", group = "Group")

# Read contrast groups
group_list <- unlist(strsplit(snakemake@wildcards[["contrast"]], "-vs-"))
meta_df <- meta_df %>%
  dplyr::filter(group %in% group_list)

# Read abundance matrix (sum up for duplicate taxonomy names)
funguild_results_df <- read.table(snakemake@input[["funguild_results"]], header=T, sep="\t", check.names=F, quote="") %>%
  dplyr::filter(`Confidence Ranking` %in% c("Highly Probable", "Probable"))

###########################
# Analysis - Trophic Mode #
###########################

trophic_plot_df <- funguild_results_df %>%
  dplyr::select(`Trophic Mode`, meta_df$sample) %>%
  dplyr::rename(trophic = "Trophic Mode") %>%
  dplyr::mutate(trophic = trimws(trophic)) %>%
  dplyr::group_by(trophic) %>%
  dplyr::summarise_each(funs(sum))

trophic_plot_df <- separate_rows(trophic_plot_df, 1, sep = "-") %>%
  dplyr::group_by(trophic) %>%
  dplyr::summarise_each(funs(sum))

trophic_plot_df[,2:ncol(trophic_plot_df)] <- t(t(trophic_plot_df[,2:ncol(trophic_plot_df)]) / colSums(trophic_plot_df[,2:ncol(trophic_plot_df)]))

write.table(trophic_plot_df,
            file = snakemake@output[["trophic_table"]],
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

trophic_plot_df <- reshape2::melt(trophic_plot_df, id.vars="trophic") %>%
  rename(sample = variable, percentage = value)

# Color palette
color_count <- length(unique(trophic_plot_df$trophic))
color_pal <- brewer.pal(n = 8, "Set1")
if (color_count <= 8) {
  get_palette = colorRampPalette(color_pal[1:color_count])
} else {
  get_palette = colorRampPalette(color_pal)
}

p <- ggplot(trophic_plot_df, aes(fill=trophic, y=percentage, x=sample)) +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = "bold")) +
  labs(x = "Sample",
       y = "Relative Abundance",
       fill = "Tropihc Mode") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values=rev(get_palette(color_count)))

ggsave(filename = snakemake@output[["trophic_barplot_png"]],
       device = "png",
       plot = p,
       width = 14,
       height = 10)
ggsave(filename = snakemake@output[["trophic_barplot_pdf"]],
       device = "pdf",
       plot = p,
       width = 14,
       height = 10)

####################
# Analysis - Guild #
####################

guild_plot_df <- funguild_results_df %>%
  dplyr::select(Guild, meta_df$sample) %>%
  dplyr::mutate(Guild = trimws(Guild)) %>%
  dplyr::group_by(Guild) %>%
  dplyr::summarise_each(funs(sum))

guild_plot_df <- separate_rows(guild_plot_df, 1, sep = "-") %>%
  dplyr::group_by(Guild) %>%
  dplyr::summarise_each(funs(sum)) %>%
  dplyr::filter(!startsWith(Guild, "|"))

guild_plot_sort_df <- guild_plot_df[order(rowMeans(guild_plot_df[, 2:ncol(guild_plot_df)]), decreasing = T),]

# Select top10 Guild
if (nrow(guild_plot_sort_df) >= 10) {
  guild_plot_top_df <- guild_plot_sort_df[1:10,]
} else {
  guild_plot_top_df <- guild_plot_sort_df
}

others_line <- c("Others", as.vector(colSums(guild_plot_sort_df[, 2:ncol(guild_plot_sort_df)]) - colSums(guild_plot_top_df[, 2:ncol(guild_plot_top_df)])))
guild_plot_top_df <- rbind(guild_plot_top_df, others_line)
for (i in 2:ncol(guild_plot_top_df)) {
  guild_plot_top_df[[i]] <- as.numeric(guild_plot_top_df[[i]])
}

order_vector <- rev(guild_plot_top_df[[1]])

guild_plot_top_df[,2:ncol(guild_plot_top_df)] <- t(t(guild_plot_top_df[,2:ncol(guild_plot_top_df)]) / colSums(guild_plot_top_df[,2:ncol(guild_plot_top_df)]))

write.table(guild_plot_top_df,
            file = snakemake@output[["guild_table"]],
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

guild_plot_top_df <- reshape2::melt(guild_plot_top_df, id.vars="Guild") %>%
  rename(sample = variable, percentage = value)

guild_plot_top_df$Guild <- factor(guild_plot_top_df$Guild, levels=order_vector)

# Color palette
color_count <- length(unique(guild_plot_top_df$Guild))
color_pal <- brewer.pal(n = 12, "Paired")
if (color_count <= 12) {
  get_palette = colorRampPalette(color_pal[1:color_count])
} else {
  get_palette = colorRampPalette(color_pal)
}

p <- ggplot(guild_plot_top_df, aes(fill=Guild, y=percentage, x=sample)) +
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = "bold")) +
  labs(x = "Sample",
       y = "Relative Abundance",
       fill = "Guild") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values=rev(get_palette(color_count)))

ggsave(filename = snakemake@output[["guild_barplot_png"]],
       device = "png",
       plot = p,
       width = 14,
       height = 10)
ggsave(filename = snakemake@output[["guild_barplot_pdf"]],
       device = "pdf",
       plot = p,
       width = 14,
       height = 10)
