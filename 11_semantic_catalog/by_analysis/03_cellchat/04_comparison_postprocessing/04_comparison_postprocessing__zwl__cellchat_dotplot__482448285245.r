load("/mnt/NexData_file/home/zwl/Project/SC2509024535-2/require3/cellchat/FGF9_NC/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/zwl/Project/SC2509024535-2/require3/cellchat/"
setwd(outdir)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)
p1 = object.list %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        source %in%  c("Fibro1",'Kupffer cells', 'Macrophages'),
        target  %in% c("Fibro1",'Kupffer cells', 'Macrophages'),
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
ggsave(filename = './cellchat_dotplot.png',p1,w=15,h=15,dpi=150)
ggsave(filename = './cellchat_dotplot.pdf',p1,w=15,h=15,dpi=150)
