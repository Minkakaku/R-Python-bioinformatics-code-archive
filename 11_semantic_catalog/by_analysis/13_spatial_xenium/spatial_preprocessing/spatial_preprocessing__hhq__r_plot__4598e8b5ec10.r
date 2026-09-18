setwd("/home/hhq/project/FWSP20255020-02_mouse_lung_xenium/")
getwd()
ls()

library(tidyverse)
library(ggplot2)
library(Seurat)
library(ggpubr)

plot(1)

plot_data = read.csv("/home/hhq/project/FWSP20255020-02_mouse_lung_xenium/metadata.csv")

head(plot_data)

plot_data <- plot_data  %>% select(group, sample, M129, celltype)

head(plot_data)


plot_data %>% group_by(sample) %>%
  summarise(M129_ratio = mean(!is.na(M129))) %>% head()

final_plot_data <- plot_data %>%
  group_by(sample) %>%
  summarise(M129_ratio = mean(M129 > 0, na.rm = TRUE))

head(final_plot_data)

final_plot_data$group = substr(final_plot_data$sample, 1, 3)
final_plot_data$group = as.factor(final_plot_data$group)

g <- ggplot(final_plot_data, aes(group, M129_ratio)) + geom_boxplot(aes(fill=group), show.legend = FALSE) + geom_point(aes(fill=sample, color=sample), , size=2) + theme_classic() + 
        theme(strip.text = element_text(size=15),axis.title.y = element_text(size=14),
              # text = element_text(size = 20, color='black'), 
              axis.text = element_text(color='black', size=14), 
              plot.title = element_text(face = "bold", size=21), 
              strip.background = element_blank(), strip.text.x = element_text(hjust = 0)) + 
        ggsci::scale_fill_d3() + ggsci::scale_color_d3() + xlab("") +  guides( color=guide_legend(title="Sample", override.aes = list(size=2)), size = F, fill=F)
g<- g + geom_pwc(aes(group = group), tip.length = 0,
                 method = "wilcox_test", label = "{p.adj.format}{p.adj.signif}",
                 p.adjust.method = "bonferroni", p.adjust.by = "group", 
                 # method.args = list(alternative = "greater"), 
                 hide.ns = F)

ggsave("boxplot_all_prop.png", g, dpi=300, height=5, width=4)




# ------------------------------
final_plot_data <- plot_data %>%
  group_by(sample, celltype) %>%
  summarise(M129_ratio = mean(M129 > 0, na.rm = TRUE))

head(final_plot_data)

final_plot_data$group = substr(final_plot_data$sample, 1, 3)
final_plot_data$group = as.factor(final_plot_data$group)

length(unique(final_plot_data$celltype))

g <- ggplot(final_plot_data, aes(group, M129_ratio)) + geom_boxplot(aes(fill=group)) + geom_point(size=2) + theme_classic() + ggforce::facet_wrap_paginate(~ celltype, ncol = 4, nrow = 3, scales='free') + 
        theme(strip.text = element_text(size=15),axis.title.y = element_text(size=14),
              # text = element_text(size = 20, color='black'), 
              axis.text = element_text(color='black', size=14), 
              plot.title = element_text(face = "bold", size=21), 
              strip.background = element_blank(), strip.text.x = element_text(hjust = 0)) + 
        ggsci::scale_fill_d3() + NoLegend() + xlab("")
g<- g + geom_pwc(aes(group = group), tip.length = 0,
                 method = "wilcox_test", label = "{p.adj.format}{p.adj.signif}",
                 p.adjust.method = "bonferroni", p.adjust.by = "group", 
                 # method.args = list(alternative = "greater"), 
                 hide.ns = F)


g <- ggplot(final_plot_data, aes(group, M129_ratio)) + geom_boxplot(aes(fill=group), show.legend = FALSE) + geom_point(aes(fill=sample, color=sample), , size=2) + theme_classic() + ggforce::facet_wrap_paginate(~ celltype, ncol = 4, nrow = 3, scales='free') + 
        theme(strip.text = element_text(size=15),axis.title.y = element_text(size=14),
              # text = element_text(size = 20, color='black'), 
              axis.text = element_text(color='black', size=14), 
              plot.title = element_text(face = "bold", size=21), 
              strip.background = element_blank(), strip.text.x = element_text(hjust = 0)) + 
        ggsci::scale_fill_d3() + ggsci::scale_color_d3() + xlab("") +  guides( color=guide_legend(title="Sample", override.aes = list(size=2)), size = F, fill=F)
g<- g + geom_pwc(aes(group = group), tip.length = 0,
                 method = "wilcox_test", label = "{p.adj.format}{p.adj.signif}",
                 p.adjust.method = "bonferroni", p.adjust.by = "group", 
                 # method.args = list(alternative = "greater"), 
                 hide.ns = F)


ggsave("boxplot_celltype_prop.png", g, dpi=300, height=14, width=14)


# barplot
head(final_plot_data)
