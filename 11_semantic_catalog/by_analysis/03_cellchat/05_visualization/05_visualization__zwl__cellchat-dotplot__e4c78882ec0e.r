library(tidyverse)
load("/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/CCl4NaCl_ctr/tmp_object.RData")
object.list.n.c <- object.list
load("/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/CCl4_FMT_ctr/tmp_object.RData")
object.list.c.f <- object.list
load("/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/CCl4_FMT_CCl4NaCl/tmp_object.RData")
object.list.n.f <- object.list
p1 = object.list.n.f %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        # source %in% 'KC Trem2+', 
        # target  %in% 'HSC', 
        
        # target  %in% c('KC Trem2+', 'HSC'),
        # pathway_name == "TGFb",

        source %in% 'Ca-LSEC',
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
  theme(strip.background = element_rect(colour = NA)) +
  scale_color_gradientn(colours  = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu"))) +
  scale_size(transform = "reverse", range = c(1, 3)) +
  Seurat::RotatedAxis()
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.KC Trem2_HSC_dotplot.png',p1,w=4,h=6)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.KC Trem2_HSC_dotplot.pdf',p1,w=4,h=6)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.TGFb_dotplot1.png',p1,w=17,h=5)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.TGFb_dotplot1.pdf',p1,w=17,h=5)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.Ca-LSEC_dotplot1.png',p1,w=17,h=15)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.Ca-LSEC_dotplot1.pdf',p1,w=17,h=15)

p2 = object.list.n.c %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        # source %in% 'KC Trem2+',
        # target  %in% 'HSC',
        
        # target  %in% c('KC Trem2+', 'HSC'),
        # pathway_name == "TGFb",
        
        source %in% 'Ca-LSEC',
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
  theme(strip.background = element_rect(colour = NA)) +
  scale_color_gradientn(colours  = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu")))+
  scale_size(transform = "reverse", range = c(1, 3)) +
  Seurat::RotatedAxis()
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.KC Trem2_HSC_dotplot2.png',p2,w=4,h=6)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.KC Trem2_HSC_dotplot2.pdf',p2,w=4,h=6)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.TGFb_dotplot2.png',p2,w=17,h=5)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.TGFb_dotplot2.pdf',p2,w=17,h=5)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.Ca-LSEC_dotplot2.png',p2,w=17,h=15)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.Ca-LSEC_dotplot2.pdf',p2,w=17,h=15)
p3 = object.list.c.f %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        # target  %in% c('KC Trem2+', 'HSC'),
        # pathway_name == "TGFb",
        
        # source %in% 'KC Trem2+',
        # target  %in% 'HSC', 

        source %in% 'Ca-LSEC',
        pval < 0.05
      )  %>%
      mutate(cell_pair = str_c(source, target, sep = " -> \n ")) 
  }) %>%
  enframe('group') %>%
  unnest(value) %>%
  
  
  ggplot(aes(x = group , y = interaction_name_2)) +
  geom_point(aes(color = prob, size = pval)) +
  facet_wrap(~ cell_pair, nrow = 1) +
  theme_classic() +
  theme(strip.background = element_rect(colour = NA)) +
  scale_color_gradientn(colours  = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu"))) +
  scale_size(transform = "reverse", range = c(1, 3)) +
  Seurat::RotatedAxis()

ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.KC Trem2_HSC_dotplot3.png',p3,w=4,h=6)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.KC Trem2_HSC_dotplot3.pdf',p3,w=4,h=6)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.TGFb_dotplot3.png',p3,w=17,h=5)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.TGFb_dotplot3.pdf',p3,w=17,h=5)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.Ca-LSEC_dotplot3.png',p3,w=17,h=15)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.Ca-LSEC_dotplot3.pdf',p3,w=17,h=15)
