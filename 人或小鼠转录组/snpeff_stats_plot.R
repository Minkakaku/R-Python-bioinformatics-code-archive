log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

for (type in c("var_region", "var_funcion", "var_impact")) {
  input_f <- snakemake@input[[type]]
  input_df <- read.table(input_f, sep="\t", header=T, check.names=FALSE)
  
  input_df2 <- input_df[, c(1, grep("_count", colnames(input_df)))]
  colnames(input_df2) <- str_replace(colnames(input_df2), "_count", "")
  
  input_df3 <- input_df2
  input_df3$var_type <- input_df3$var_type
  gathercols <- colnames(input_df3)[-1]
  input_df3 <- gather(input_df3, sample, count, gathercols, factor_key=TRUE)
  
  p <- ggplot(input_df3, aes(fill = var_type, y = count, x = sample)) + 
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_brewer(type,
                      palette = "Paired") +
    labs(x = "Sample",
         y = "Variant Count") + 
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) + 
    theme_bw() +
    theme(axis.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold", angle=45, hjust = 1, vjust = 1))
  
  ggsave(filename = snakemake@output[[type]][1],
         plot = p,
         width = 9,
         height = 8)
  
  ggsave(filename = snakemake@output[[type]][2],
         plot = p,
         width = 9,
         height = 8)
}
