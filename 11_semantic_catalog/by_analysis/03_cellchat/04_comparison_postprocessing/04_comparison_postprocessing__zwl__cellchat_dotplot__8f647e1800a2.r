library(sceasy)
library(Seurat)
library(reticulate)
library(CellChat)
rm(list = ls())
gc()
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/SC2509024535/require2/cellchat/anno_raw.h5ad', from="anndata", to="seurat",
                      outFile='Mono.rds')

CellChatDB <-  CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB,search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))

setwd('/mnt/NexData_file/home/zwl/Project/SC2509024535/require2/cellchat/FGF9_NC_sub/')
outdir = '/mnt/NexData_file/home/zwl/Project/SC2509024535/require2/cellchat/FGF9_NC_sub/'
load("/mnt/NexData_file/home/zwl/Project/SC2509024535/require2/cellchat/FGF9_NC_sub/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))


vertex_cal = function (object, signaling, slot.name = "netP"){
  centr0 <- slot(object, slot.name)$centr[signaling]
  centr0 <- centr0[[1]]
  mat <- matrix(unlist(centr0), ncol = length(centr0), byrow = FALSE)
  mat <- t(mat)
  rownames(mat) <- names(centr0); colnames(mat) <- names(centr0$outdeg)
  mat <- mat[c("outdeg","indeg"),]
  mat <- as.data.frame(t(mat))
  top3_indices <- head(order(mat$outdeg, decreasing = TRUE), 3)
}

pathway_outdir<-paste0(outdir,"pathway/")
dir.create(pathway_outdir, recursive = T)
select_pathway<-c("FGF",'CSF','IL6','CCL','TNF','TGFb','VCAM')
for(k in str_c(select_pathway)){
  pairLR_top <- extractEnrichedLR(object.list[[1]], signaling = k, geneLR.return = FALSE)
  pairLR_top <- pairLR_top[1,]
  vertex_receiver <- vertex_cal(object.list[[1]],k)
  pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-centrality_heatmap.pdf"),width = 8,height= 3)
  gg <- netAnalysis_signalingRole_network(object.list[[1]], signaling = k, width = 14, height = 2.5, font.size = 10)
  dev.off()
  
  pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-contribution_barplot.pdf"),width = 8,height= 3)
  print(netAnalysis_contribution(object.list[[1]], signaling = k))
  dev.off()
  
  pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
  netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
  dev.off()
}

for(k in str_c(select_pathway)){
  pairLR_top <- extractEnrichedLR(object.list[[2]], signaling = k, geneLR.return = FALSE)
  pairLR_top <- pairLR_top[1,]
  vertex_receiver <- vertex_cal(object.list[[2]],k)
  pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-centrality_heatmap.pdf"),width = 8,height= 3)
  gg <- netAnalysis_signalingRole_network(object.list[[2]], signaling = k, width = 14, height = 2.5, font.size = 10)
  dev.off()
  
  pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-contribution_barplot.pdf"),width = 8,height= 3)
  print(netAnalysis_contribution(object.list[[2]], signaling = k))
  dev.off()
  
  pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
  netVisual_individual(object.list[[2]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
  dev.off()
}

pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
  file_name <- pdf_files[i]
  pdf_image <- pdf_render_page(file_name, dpi = 300)
  writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}
pdf_files <- list.files(path = pathway_outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
  file_name <- pdf_files[i]
  pdf_image <- pdf_render_page(file_name, dpi = 300)
  writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}


p1 = object.list %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(

        (source %in% 'Osteocytes'),#|(target  %in% 'Osteocytes'), 
        #pathway_name %in% c("FGF",'CSF','IL6','CCL','TNF','TGFb','VCAM'),
        
        # source %in% 'Ca-LSEC',
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