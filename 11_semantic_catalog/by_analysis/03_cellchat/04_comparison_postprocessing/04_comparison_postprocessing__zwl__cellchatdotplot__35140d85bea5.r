setwd('/mnt/NexData_file/home/zwl/Project/SC2511076348/require9')
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

load("/mnt/NexData_file/home/zwl/Project/SC2511076348/require8/6.cellchat/DM_N/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/zwl/Project/SC2511076348/require9/12.cellchat/pathway/"

p1 = object.list %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        source %in%  c("Macro1","Osteoclast","Macro2","BMSC","Osteoblast"),
        target  %in% c("Macro1","Osteoclast","Macro2","BMSC","Osteoblast"),
        # interaction_name_2 %in% pathway,
        # target  %in% c('KC Trem2+', 'HSC'),
        # pathway_name %in% c('APP','LAIR1'),
        pval < 0.05
      )  %>%
      mutate(cell_pair = str_c(source, target, sep = " -> \n")) 
  }) %>%
  enframe('group') %>%
  unnest(value) %>%
  
  
  ggplot(aes(x = group , y = interaction_name_2)) +
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
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require9/12.cellchat/pathway/cellchat_mono_dotplot.png',p1,w=30,h=30,dpi=150)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/SC2511076348/require9/12.cellchat/pathway/cellchat_mono_dotplot.pdf',p1,w=30,h=30,dpi=150)
