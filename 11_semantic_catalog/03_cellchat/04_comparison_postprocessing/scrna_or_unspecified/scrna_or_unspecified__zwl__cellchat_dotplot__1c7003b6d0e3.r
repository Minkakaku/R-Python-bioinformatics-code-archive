library(sceasy)
library(Seurat)
library(reticulate)
library(CellChat)
rm(list = ls())
gc()
sceasy::convertFormat('/mnt/NexData_file/home/zwl/Project/SC2509024535/require2/cellchat/anno_raw.h5ad', from="anndata", to="seurat",
                      outFile='Mono.rds')
marker = list(     'T cells' = c('CD3D','CD3E','CD3G'),     
                   'NK cells' = c('GNLY','NCAM1','NCR1','NKG7','KLRD1'),     
                   'B cells' = c('CD19','CD22','CD79A','CD79B','MS4A1'),     
                   'Plasma' = c('CD79A','JCHAIN','TNFRSF17'),     
                   'Neutrophil' = c('CSF3R','FCGR3B','LRG1','S100A8','S100A9'),     
                   'Mast cell' = c('CPA3','GATA2','MS4A2','TPSAB1'),     
                   'Epithelium' = c('EPCAM','KRT14','KRT15','KRT18','KRT19','MUC1'),     
                   'Fibroblast' = c('COL1A1','COL1A2','DCN','LUM'),     
                   'Endothelium' = c('PECAM1','VEGFA','CLDN5'),     
                   'Macrophage' = c('CD68','CD163','MS4A4A'),     
                   'DC' = c('CD1C','CLEC9A','FCER1A'),     
                   'Mast cell' = c('CPA3','GATA2','MS4A2','TPSAB1')   )
CellChatDB <-  CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB,search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"))

setwd('/mnt/NexData_file/home/zwl/Project/FWSC20241172/require13/cellchat/Immune_EXP_CON')
outdir = '/mnt/NexData_file/home/zwl/Project/FWSC20241172/require13/cellchat/Immune_EXP_CON/'
load("/mnt/NexData_file/home/zwl/Project/FWSC20241172/require13/cellchat/Immune_EXP_CON/tmp_object.RData")
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

pathway_outdir<-paste0(outdir,"/pathway/")
dir.create(pathway_outdir, recursive = T)
select_pathway<-c('IFN-II')
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
