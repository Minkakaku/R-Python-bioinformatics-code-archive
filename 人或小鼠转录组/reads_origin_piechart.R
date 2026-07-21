log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(reshape2)
library(ggplot2)
library(dplyr)
library(scales)
library(ggforce)

qualimap_rnaseq_genome_results_file <- snakemake@input[[1]]

reads_origin_df <- read.table(qualimap_rnaseq_genome_results_file, sep="\t", header=TRUE) %>%
  select(Sample, reads_aligned_exonic, reads_aligned_intronic, reads_aligned_intergenic) %>%
  rename(Exon = reads_aligned_exonic,
         Intron = reads_aligned_intronic,
         Intergenic = reads_aligned_intergenic)

for (f in snakemake@output[["png"]]) {
   sample <- gsub("_reads_origin_piechart.png", "", basename(f))
   
   plot_df <- reads_origin_df %>%
     filter(Sample == sample)
   plot_df <- melt(plot_df, id.vars = "Sample")
   plot_df <- plot_df %>% 
     mutate(end = 2 * pi * cumsum(value)/sum(value),
            start = lag(end, default = 0),
            middle = 0.5 * (start + end),
            hjust = ifelse(middle > pi, 1, 0),
            vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))
   plot_df$label <- paste(round(((plot_df$value/sum(plot_df$value))*100), 2),"%")
   
   color <- c("Exon" = "#A6CEE3",
              "Intron" = "#B2DF8A",
              "Intergenic" = "#FDBF6F")
   
   pie <- ggplot(plot_df) + 
     geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                      start = start, end = end, fill = variable),
                  color = "transparent") +
     geom_text(aes(#x = 1.05 * sin(middle), 
                   #y = 1.05 * cos(middle), 
                   x = 1.05*sin(middle),
                   y = 1.05*cos(middle), 
                   label = label,
                   hjust = 0.5, 
                   vjust = vjust)) +
     coord_fixed() +
     scale_x_continuous(limits = c(-1.1, 1.1),  # Adjust so labels are not cut off
                        name = "", breaks = NULL, labels = NULL) +
     scale_y_continuous(limits = c(-1.1, 1.1),    # Adjust so labels are not cut off
                        name = "", breaks = NULL, labels = NULL) +
     theme_bw() + 
     theme(panel.border = element_blank(),
           legend.title = element_blank(),
           legend.text = element_text(size=9)) + 
     scale_fill_manual(values = color)
  
   ggsave(filename = f,
          plot = pie,
          device = "png",
          height = 8,
          width = 10)
   
   ggsave(filename = gsub("png", "pdf", f),
          plot = pie,
          device = "pdf",
          height = 8,
          width = 10) 

}
