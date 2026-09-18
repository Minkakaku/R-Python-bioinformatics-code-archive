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

load("/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/Fibroblast-CD8T/tmp_object.RData")
# cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/Fibroblast-CD8T/"
data<-netVisual_bubble(object.list[[1]], sources.use = c(1:5), targets.use = c(6), angle.x = 45,return.data = TRUE)
data<-data$communication[data$communication$annotation=="Secreted Signaling",]
signaling_use<-unique(data$pathway_name)
pdf(paste0(outdir,"netVisual_bubble.pdf"),width = 12,height= 6)
netVisual_bubble(object.list[[1]], sources.use = c(1:5), targets.use = c(6), angle.x = 45,signaling=signaling_use, remove.isolate = TRUE, dot.size.min = 3, dot.size.max = 12, font.size = 15, font.size.title = 15) + coord_flip()
dev.off()

# pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
# for (i in seq_along(pdf_files)) {
    # file_name <- pdf_files[i]
		file_name <- paste0(outdir,"netVisual_bubble.pdf")
    pdf_image <- pdf_render_page(file_name, dpi = 300)
    writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
# }

load("/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/Macrophage-CD8T/tmp_object.RData")
outdir<-"/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/Macrophage-CD8T/"
data<-netVisual_bubble(object.list[[1]], sources.use = c(1:6), targets.use = c(7), angle.x = 45,return.data = TRUE)
data<-data$communication[data$communication$annotation=="Secreted Signaling",]
signaling_use<-unique(data$pathway_name)
pdf(paste0(outdir,"netVisual_bubble.pdf"),width = 12,height= 6)
netVisual_bubble(object.list[[1]], sources.use = c(1:6), targets.use = c(7), angle.x = 45,signaling=signaling_use, remove.isolate = TRUE, dot.size.min = 3, dot.size.max = 12, font.size = 15, font.size.title = 15) + coord_flip()
dev.off()

file_name <- paste0(outdir,"netVisual_bubble.pdf")
pdf_image <- pdf_render_page(file_name, dpi = 300)
writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))