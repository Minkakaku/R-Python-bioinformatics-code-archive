library(tidyverse)
setwd('/mnt/NexData_file/home/zwl/Project/SC2507233885/require4')
library(CellChat)
rm(list = ls())
gc()
load("/mnt/NexData_file/home/zwl/Project/SC2507233885/require3/cellchat/CASP3_CT/tmp_object.RData")
object.list.n.f <- object.list
p1 = object.list.n.f %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        source %in% 'L2-3 IT-MT',
        target  %in% 'L5/6 NP2',
        
        # target  %in% c('KC Trem2+', 'HSC'),
        # pathway_name == "TGFb",

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
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/SC2507233885/require4/2.L2-3 IT-MT_L5-6 NP2_dotplot.png',p1,w=6,h=12)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/SC2507233885/require4/2.L2-3 IT-MT_L5-6 NP2_dotplot.pdf',p1,w=6,h=12)

# p2 = object.list.n.c %>% 
#   map(subsetCommunication) %>%
#   map(function(x){
#     x %>% 
#       dplyr::filter(
#         # source %in% 'KC Trem2+',
#         # target  %in% 'HSC',
#         
#         # target  %in% c('KC Trem2+', 'HSC'),
#         # pathway_name == "TGFb",
#         
#         source %in% 'Ca-LSEC',
#         pval < 0.05
#       )  %>%
#       mutate(cell_pair = str_c(source, target, sep = " -> \n")) 
#   }) %>%
#   enframe('group') %>%
#   unnest(value) %>%
#   
#   
#   ggplot(aes(x = group , y = interaction_name_2)) +
#   geom_point(aes(color = prob, size = pval)) +
#   facet_wrap(~ cell_pair, nrow = 1) +
#   theme_classic() +
#   theme(strip.background = element_rect(colour = NA)) +
#   scale_color_gradientn(colours  = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu")))+
#   scale_size(transform = "reverse", range = c(1, 3)) +
#   Seurat::RotatedAxis()
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.KC Trem2_HSC_dotplot2.png',p2,w=4,h=6)
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.KC Trem2_HSC_dotplot2.pdf',p2,w=4,h=6)
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.TGFb_dotplot2.png',p2,w=17,h=5)
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.TGFb_dotplot2.pdf',p2,w=17,h=5)
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.Ca-LSEC_dotplot2.png',p2,w=17,h=15)
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.Ca-LSEC_dotplot2.pdf',p2,w=17,h=15)
# p3 = object.list.c.f %>% 
#   map(subsetCommunication) %>%
#   map(function(x){
#     x %>% 
#       dplyr::filter(
#         # target  %in% c('KC Trem2+', 'HSC'),
#         # pathway_name == "TGFb",
#         
#         # source %in% 'KC Trem2+',
#         # target  %in% 'HSC', 
#         
#         source %in% 'Ca-LSEC',
#         pval < 0.05
#       )  %>%
#       mutate(cell_pair = str_c(source, target, sep = " -> \n ")) 
#   }) %>%
#   enframe('group') %>%
#   unnest(value) %>%
#   
#   
#   ggplot(aes(x = group , y = interaction_name_2)) +
#   geom_point(aes(color = prob, size = pval)) +
#   facet_wrap(~ cell_pair, nrow = 1) +
#   theme_classic() +
#   theme(strip.background = element_rect(colour = NA)) +
#   scale_color_gradientn(colours  = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu"))) +
#   scale_size(transform = "reverse", range = c(1, 3)) +
#   Seurat::RotatedAxis()
# 
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.KC Trem2_HSC_dotplot3.png',p3,w=4,h=6)
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.KC Trem2_HSC_dotplot3.pdf',p3,w=4,h=6)
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.TGFb_dotplot3.png',p3,w=17,h=5)
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.TGFb_dotplot3.pdf',p3,w=17,h=5)
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.Ca-LSEC_dotplot3.png',p3,w=17,h=15)
# ggsave(filename = '/mnt/NexData_file/home/zwl/Project/FWSC20251383-02/require4_cellchat/6.Ca-LSEC_dotplot3.pdf',p3,w=17,h=15)


setwd('/mnt/NexData_file/home/zwl/Project/SC2507233885/require4')
CASP = read.csv('../require3/cellchat/CASP3_CT/CASP3-net_lr.txt',sep = '\t')
CT = read.csv('../require3/cellchat/CASP3_CT/CT-net_lr.txt',sep = '\t')
# 合并两个数据表
merged_data <- merge(CASP, CT, by = c("relation", "source", "target", "ligand", "receptor", "interaction_name",
                                      "interaction_name_2", "pathway_name", "annotation", "evidence"),)
colnames(merged_data) = c("relation","source","target","ligand","receptor","interaction_name","interaction_name_2",
                          "pathway_name","annotation","evidence","prob_CASP","pval_CASP","prob_CT","pval_CT" )
merged_data$CASP_CT = merged_data$prob_CASP- merged_data$prob_CT
write.csv(file = '/mnt/NexData_file/home/zwl/Project/SC2507233885/require4/2.cellchat_CASP3-CT.csv',merged_data,row.names = F)
