setwd('/mnt/NexData_file/home/zwl/Project/FWSC20251792/require10/')
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(ComplexHeatmap)
library(argparser)
# library(pdftools)
library(stringr)
library(jpeg)
library(cowplot)
library(Seurat)
library(RColorBrewer)

load("/mnt/NexData_file/home/zwl/Project/FWSC20251792/require10/EDA1_Control1/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/zwl/Project/FWSC20251792/require10/"

data = object.list %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        source %in%  c("Mono_Vcan","other"),
        target  %in% c("CD4_Cd5","CD4_Lef1","CD8_Ccl4","CD8_Ccl5","CD8_Ccl6","CD8_Cxcr4","CD8_Gzmf","CD8_Mki67","CD8_Xcl1"),
        # interaction_name_2 %in% pathway,
        # target  %in% c('KC Trem2+', 'HSC'),
        # pathway_name %in% c('APP','LAIR1'),
        pval < 0.05
      )  %>%
      mutate(cell_pair = str_c(source, target, sep = " -> \n")) 
  }) %>%
  enframe('group') %>%
  unnest(value)
  
write.csv(file = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require10/plot_data.csv',data)
p1 =  ggplot(data,aes(x = cell_pair , y = interaction_name_2)) +
  geom_point(aes(color = prob, size = pval)) +
  # facet_wrap(~ cell_pair, nrow = 1) +
  theme_classic() +
  theme(strip.background = element_rect(colour = NA),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) +
  scale_color_gradientn(colours  = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu"))) +
  scale_size(transform = "reverse", range = c(1, 3)) +
  Seurat::RotatedAxis()
p1
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require10/cellchat_mono_dotplot.png',p1,w=15,h=20,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require10/cellchat_mono_dotplot.pdf',p1,w=15,h=20,dpi=150)
p1 =  ggplot(data,aes(x = group , y = interaction_name_2)) +
  geom_point(aes(color = prob, size = pval)) +
  facet_wrap(~ cell_pair, nrow = 1) +
  theme_classic() +
  theme(strip.background = element_rect(colour = NA),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) +
  scale_color_gradientn(colours  = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu"))) +
  scale_size(transform = "reverse", range = c(1, 3)) +
  Seurat::RotatedAxis()
p1
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require10/cellchat_mono_dotplot_group.png',p1,w=25,h=20,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251792/require10/cellchat_mono_dotplot_group.pdf',p1,w=25,h=20,dpi=150)
