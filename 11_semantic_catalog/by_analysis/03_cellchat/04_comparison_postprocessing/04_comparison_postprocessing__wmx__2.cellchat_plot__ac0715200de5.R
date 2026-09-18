library(CellChat)
library(jpeg)
library(pdftools,lib.loc = "/mnt/NexData_file/share/scRNA-tools/envs/pip_cellchat/lib/R/library/")
setwd("/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou13/")
outdir = "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou13/cellchat_plot/"
dir.create(outdir)
load("/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou12/cellchat/LARC_N/tmp_object.RData")

cellchat <- mergeCellChat(object.list,add.names = names(object.list))

for(j in 1:2){
  mat_count = cellchat@net[[j]]$count
  mat_weight = cellchat@net[[j]]$weight
  # par(mfrow=n,mar=c(0.3,0.3,0.3,0.3),xpd=TRUE)
  for (i in 1:nrow(mat_count)) {
    pdf(paste0(outdir,"",names(cellchat@net)[j],"_",rownames(mat_count)[i],"_Number_netVisual_circle.pdf"),width = 10,height = 10)
    mat_count2 = matrix(0,nrow = nrow(mat_count),ncol = ncol(mat_count),dimnames = dimnames(mat_count))
    mat_count2[i,] = mat_count[i,]
    netVisual_circle(mat_count2,vertex.weight = 11,weight.scale = T,label.edge = T,edge.weight.max = max(mat_count),edge.label.cex = 1,title.name = rownames(mat_count)[i])
    dev.off()
    }
  # par(mfrow=n,mar=c(0.3,0.3,0.3,0.3),xpd=TRUE)
  for (i in 1:nrow(mat_count)) {
    pdf(paste0(outdir,"",names(cellchat@net)[j],"_",rownames(mat_count)[i],"_Weight_netVisual_circle.pdf"),width = 10,height = 10)
    mat_weight2 = matrix(0,nrow = nrow(mat_weight),ncol = ncol(mat_weight),dimnames = dimnames(mat_weight))
    mat_weight2[i,] = mat_weight[i,]
    netVisual_circle(mat_weight2,vertex.weight = 11,weight.scale = T,label.edge = T,edge.weight.max = max(mat_weight),edge.label.cex = 1,title.name = rownames(mat_weight)[i])
    dev.off()
    }
}

pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
  file_name <- pdf_files[i]
  pdf_image <- pdf_render_page(file_name, dpi = 300)
  writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}
cat("Processed", length(pdf_files), "PDF files.\n")

