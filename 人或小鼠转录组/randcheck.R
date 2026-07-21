log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)
library(ggplot2)

# input: coverage_profile_along_genes_(total).txt by qualimap
input <- read.table(snakemake@input[[1]], sep="\t", comment.char="", header=T)
colnames(input) <- c("pos", "reads")
input <- input %>%
  mutate(percentage = round(reads/sum(reads)*100, 2),
         pos = round(pos/100, 2))

# plot
p <- ggplot(input) +
  geom_line(aes(pos, percentage), 
            color = "#377EB8") + 
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold",
                                  hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = "bold")) +
  labs(x = "Relativave Position in Genes (5' to 3')", 
       y = "Percentage of Reads (%)")

ggsave(p, 
       filename = snakemake@output[["plot_pdf"]], 
       device = "pdf",
       width = 10, 
       height = 6)
ggsave(p, 
       filename = snakemake@output[["plot_png"]], 
       device = "png",
       width = 10, 
       height = 6)
