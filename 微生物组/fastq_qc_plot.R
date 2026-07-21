log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(dplyr)
library(ggplot2)
library(rjson)
library(reshape2)

# Input
## fastp results json
fastp_json <- snakemake@input[[1]]
sample <- gsub("_fastp.json", "", basename(fastp_json))
qual_plot_png <- snakemake@output[["quality_plot_png"]]
qual_plot_pdf <- snakemake@output[["quality_plot_pdf"]]
error_plot_png <- snakemake@output[["error_plot_png"]]
error_plot_pdf <- snakemake@output[["error_plot_pdf"]]
gc_plot_png <- snakemake@output[["gc_plot_png"]]
gc_plot_pdf <- snakemake@output[["gc_plot_pdf"]]

# Read json file
qc_data <- fromJSON(file=fastp_json)

## Base quality score and Base error rate
reads1_quality_df <- as.data.frame(qc_data$read1_before_filtering$quality_curves) %>% 
  dplyr::select(mean) %>%
  dplyr::rename(quality_score = mean) %>%
  mutate(read = "Read1")
reads1_quality_df$base <- as.integer(rownames(reads1_quality_df))

reads2_quality_df <- as.data.frame(qc_data$read2_before_filtering$quality_curves) %>% 
  dplyr::select(mean) %>%
  dplyr::rename(quality_score = mean) %>%
  mutate(read = "Read2")
reads2_quality_df$base <- as.integer(rownames(reads2_quality_df))

read_quality_df <- rbind(reads1_quality_df, reads2_quality_df)

# quality_score = -10*log10(error_rate) 
read_quality_df$error_rate <- 10^(read_quality_df$quality_score/(-10))

## Base GC content
reads1_gc_df <- as.data.frame(qc_data$read1_before_filtering$content_curves) %>% 
  dplyr::select(-GC)
reads1_gc_df$base <- rownames(reads1_gc_df)
reads1_gc_df_long <- melt(reads1_gc_df, id.vars=c("base")) %>%
  mutate(read = "Read1")

reads2_gc_df <- as.data.frame(qc_data$read2_before_filtering$content_curves) %>% 
  dplyr::select(-GC)
reads2_gc_df$base <- rownames(reads2_gc_df)
reads2_gc_df_long <- melt(reads2_gc_df, id.vars=c("base")) %>%
  mutate(read = "Read2")

read_gc_df <- rbind(reads1_gc_df_long, reads2_gc_df_long) %>%
  dplyr::rename(atcgn = variable, percent = value) %>%
  mutate(base = as.integer(base),
         percent = percent*100)

# Plot
base_max <- max(read_quality_df$base)
## Plot1: Base Quality Score Distribution
p_qual <- ggplot(read_quality_df, 
                 aes(x = base, y = quality_score)) +
  geom_line(color = "#377EB8") + 
  theme_bw() +
  facet_wrap(~read, scales = "free_x", strip.position = "bottom") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold",
                                  hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = "bold")) +
  labs(title = paste0("Base Quality Score Distribution", " (", sample, ")"), 
       x = NULL, 
       y = "Average Base Quality Score") +
  scale_x_continuous(expand = c(0, 5),
                     limits = c(0, base_max+1),
                     breaks = seq(0, base_max, 50)) +
  scale_y_continuous(expand = c(0, 10),
                     breaks = seq(0, 40, 10),
                     limits = c(0, 40))

ggsave(filename = qual_plot_png,
       plot = p_qual,
       width = 10,
       height = 8)
ggsave(filename = qual_plot_pdf,
       plot = p_qual,
       width = 10,
       height = 8)

## Plot2: Base Error Rate Distribution
error_rate_max <- max(read_quality_df$error_rate)
read_quality_df$error_rate[is.na(read_quality_df$error_rate)]

p_error <- ggplot(read_quality_df, 
                  aes(x = base, y = error_rate)) +
  geom_bar(stat = "identity",
           fill = "#377EB8") + 
  theme_bw() +
  facet_wrap(~read, scales = "free_x", strip.position = "bottom") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold",
                                  hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = "bold")) +
  labs(title = paste0("Base Error Rate Distribution", " (", sample, ")"), 
       x = NULL, 
       y = "Average Base Error Rate") +
  scale_x_continuous(expand = c(0, 5),
                     limits = c(0, base_max+1),
                     breaks = seq(0, base_max, 50)) +
  scale_y_continuous(expand = c(0.05, 0),
                     labels = scales::comma,
                     limits = c(0, error_rate_max*1.5))

ggsave(filename = error_plot_png, 
       plot = p_error,
       width = 10,
       height = 8)
ggsave(filename = error_plot_pdf, 
       plot = p_error,
       width = 10,
       height = 8)

## Plot3: Base GC Content Distribution
p_gc <- ggplot(read_gc_df, 
            aes(x = base, y = percent, color = atcgn)) +
  geom_line() + 
  theme_bw() +
  facet_wrap(~read, scales = "free_x", strip.position = "bottom") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold",
                                  hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = "bold")) +
  labs(title = paste0("Base GC Content Distribution", " (", sample, ")"), 
       x = NULL, 
       y = "Average Base GC Content") +
  scale_x_continuous(expand = c(0, 5),
                     limits = c(0, base_max+1),
                     breaks = seq(0, base_max, 50)) +
  scale_y_continuous(expand = c(0, 2),
                     breaks = seq(0, 50, 10),
                     limits = c(0, 50))

ggsave(filename = gc_plot_png, 
       plot = p_gc,
       width = 10,
       height = 8)
ggsave(filename = gc_plot_pdf, 
       plot = p_gc,
       width = 10,
       height = 8)
