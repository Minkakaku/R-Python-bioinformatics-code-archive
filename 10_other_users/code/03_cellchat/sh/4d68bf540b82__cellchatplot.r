library(cowplot)
library(Seurat)
library(RColorBrewer)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(ComplexHeatmap)
library(argparser)
library(pdftools)
library(stringr)
library(jpeg)
library(rmarkdown)

load("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/LAHS_NLAHS/tmp_object.RData")
load("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/LAHS_NC/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC20250829/19.Mission/cellchat_plot/"
# pdf(paste0(outdir,"Non-classical Monocytes_Number_netVisual_circle.pdf"),width = 12,height = 12)
# netVisual_diffInteraction(cellchat,measure = "count", comparison = c(2, 1),weight.scale = T,label.edge = T,edge.label.cex = 1,title.name = 'Non-classical Monocytes',sources.use ='Non-classical Monocytes')
# dev.off()
# pdf(paste0(outdir,"Non-classical Monocytes_weight_netVisual_circle.pdf"),width = 12,height = 12)
# netVisual_diffInteraction(cellchat,measure = "weight", comparison = c(2, 1),weight.scale = T,label.edge = T,edge.label.cex = 1,title.name = 'Non-classical Monocytes',sources.use ='Non-classical Monocytes')
# dev.off()

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

pathway_outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC2024-192/19.Mission/cellchat_plot/"
# dir.create(pathway_outdir, recursive = T)
select_pathway<-c('THBS')
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

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-",pairLR_top,"-top_hierarchy.pdf"),width = 18,height= 8)
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

pdf_files <- list.files(path = pathway_outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
    file_name <- pdf_files[i]
    pdf_image <- pdf_render_page(file_name, dpi = 300)
    writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}