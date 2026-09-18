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

load("/mnt/NexData_file/home/sh/Projects/AE2512027077-01/01.Mission/OVX_bone-CNS_sham_bone-CNS/tmp_object.RData")
outdir<-"/mnt/NexData_file/home/sh/Projects/AE2512027077-01/12.Mission/OVX_bone-CNS_sham_bone-CNS_plot/"
	#提前生成标题
# cellchat <- mergeCellChat(object.list,add.names = names(object.list))
title_Number <- c(paste0("Number of interactions - ",names(object.list)[1]),paste0("Number of interactions - ",names(object.list)[2]))
title_Weight <- c(paste0("Weight of interactions - ",names(object.list)[1]),paste0("Weight of interactions - ",names(object.list)[2]))
title_name <- c(names(object.list)[1],names(object.list)[2])

title_list <- list(title_Number,title_Weight)
labeledge<-"F"
## 细胞通讯网络可视化
#分组细胞类型相互作用数量弦图，以最大count做scale
count.max <- getMaxWeight(object.list, attribute = c("idents","count"))
weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
#### 初版
pdf(paste0(outdir,"Number_of_interactions_netplot_raw.pdf"),width = 12,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= labeledge, edge.weight.max = count.max[2], edge.width.max = 12, title.name = title_Number[i])
	}
dev.off()
#分组细胞类型相互作用权重弦图，以最大weight做scale
pdf(paste0(outdir,"Weight_of_interactions_netplot_raw.pdf"),width = 12,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= labeledge, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = title_Weight[i])
	}
dev.off()
#### 无标签版
pdf(paste0(outdir,"Number_of_interactions_netplot_unlabel.pdf"),width = 12,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= labeledge, edge.weight.max = count.max[2], edge.width.max = 12, title.name = title_Number[i])
	}
dev.off()
#分组细胞类型相互作用权重弦图，以最大weight做scale
pdf(paste0(outdir,"Weight_of_interactions_netplot_unlabel.pdf"),width = 12,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= labeledge, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = title_Weight[i])
	}
dev.off()
#### 放大版
pdf(paste0(outdir,"Number_of_interactions_netplot_p.pdf"),width = 24,height = 14)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= labeledge, edge.weight.max = count.max[2], edge.width.max = 12, title.name = title_Number[i])
	}
dev.off()
#分组细胞类型相互作用权重弦图，以最大weight做scale
pdf(paste0(outdir,"Weight_of_interactions_netplot_p.pdf"),width = 24,height = 14)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= labeledge, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = title_Weight[i])
	}
dev.off()
#### 筛选版
selected_cell_type <- c("Activate memoryB","GCB","memoryB","Microglia","Neutrophil_G1","Neutrophil_G2","Neutrophil_G3","Neutrophil_G4","Neutrophil_G5")
pdf(paste0(outdir,"Number_of_interactions_netplot_select.pdf"),width = 14,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	sub_matrix_count <- object.list[[i]]@net$count[selected_cell_type, selected_cell_type]
	netVisual_circle(sub_matrix_count, weight.scale = T, label.edge= "F", edge.width.max = 12, title.name = title_Number[i])
	}
dev.off()
#分组细胞类型相互作用权重弦图，以最大weight做scale
pdf(paste0(outdir,"Weight_of_interactions_netplot_select.pdf"),width = 14,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	sub_matrix_weight <- object.list[[i]]@net$weight[selected_cell_type, selected_cell_type]
	netVisual_circle(sub_matrix_weight, weight.scale = T, label.edge= "F", edge.width.max = 12, title.name = title_Weight[i])
	}
dev.off()

#### 筛选版贝壳图
selected_cell_type <- c("Activate memoryB","GCB","memoryB","Microglia","Neutrophil_G1","Neutrophil_G2","Neutrophil_G3","Neutrophil_G4","Neutrophil_G5")
pdf(paste0(outdir,"Number_of_interactions_circle_select.pdf"),width = 14,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	sub_matrix_count <- object.list[[i]]@net$count[selected_cell_type, selected_cell_type]
	mat_count2 = matrix(0,nrow = nrow(sub_matrix_count),ncol = ncol(sub_matrix_count),dimnames = dimnames(sub_matrix_count))
	mat_count2[c(8,9),] = sub_matrix_count[c(8,9),]
	netVisual_circle(mat_count2, weight.scale = T, label.edge= "F", edge.width.max = 12, title.name = title_Number[i])
	}
dev.off()
#分组细胞类型相互作用权重弦图，以最大weight做scale
pdf(paste0(outdir,"Weight_of_interactions_circle_select.pdf"),width = 14,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	sub_matrix_weight <- object.list[[i]]@net$weight[selected_cell_type, selected_cell_type]
	mat_weight2 = matrix(0,nrow = nrow(sub_matrix_weight),ncol = ncol(sub_matrix_weight),dimnames = dimnames(sub_matrix_weight))
	mat_weight2[c(8,9),] = sub_matrix_weight[c(8,9),]
	netVisual_circle(mat_weight2, weight.scale = T, label.edge= "F", edge.width.max = 12, title.name = title_Weight[i])
	}
dev.off()

#### 筛选版弦图
pathways.show <- c('CCL')
pdf(paste0(outdir,"Number_of_interactions_chord_select_p.pdf"),width = 18,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_aggregate(object.list[[i]],signaling = pathways.show,sources.use = c(18,19), remove.isolate = TRUE,layout = "chord", label.edge= "F", signaling.name = paste("all", names(object.list)[i]))
	}
dev.off()
pdf(paste0(outdir,"Weight_of_interactions_chord_select_p.pdf"),width = 18,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_aggregate(object.list[[i]],signaling = pathways.show,sources.use = c(18,19),remove.isolate = TRUE, layout = "chord", label.edge= "F",signaling.name = paste("all", names(object.list)[i]))
	}
dev.off()

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
pdf(paste0(outdir,"interactions_chord_gene-OVX_bone-CNS.pdf"),width = 10,height = 10)
netVisual_chord_gene(object.list[[1]], sources.use = c(18,19), signaling = c("CCL"),legend.pos.x = 8)
dev.off()
pdf(paste0(outdir,"interactions_chord_gene-sham_bone-CNS.pdf"),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], sources.use = c(18,19), signaling = c("CCL"),legend.pos.x = 8)
dev.off()

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
levels(object.list[[1]]@idents)
pdf(paste0(outdir,"interactions_chord_gene_pairL-OVX_bone-CNS.pdf"),width = 10,height = 10)
netVisual_chord_gene(object.list[[1]], sources.use = c(19,20),targets.use = 13, slot.name = "netP", legend.pos.x = 10)
dev.off()
pdf(paste0(outdir,"interactions_chord_gene_pairL-sham_bone-CNS.pdf"),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], sources.use = c(19,20),targets.use = 13, slot.name = "netP", legend.pos.x = 10)
dev.off()

pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
    file_name <- pdf_files[i]
    pdf_image <- pdf_render_page(file_name, dpi = 300)
    writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}




load("/mnt/NexData_file/home/sh/Projects/AE2512027077-01/04.Mission/OVX_CNS-defect_sham_CNS-defect/tmp_object.RData")
outdir<-"/mnt/NexData_file/home/sh/Projects/AE2512027077-01/12.Mission/OVX_CNS-defect_sham_CNS-defect_plot/"
	#提前生成标题
# cellchat <- mergeCellChat(object.list,add.names = names(object.list))
title_Number <- c(paste0("Number of interactions - ",names(object.list)[1]),paste0("Number of interactions - ",names(object.list)[2]))
title_Weight <- c(paste0("Weight of interactions - ",names(object.list)[1]),paste0("Weight of interactions - ",names(object.list)[2]))
title_name <- c(names(object.list)[1],names(object.list)[2])

title_list <- list(title_Number,title_Weight)
labeledge<-"F"
## 细胞通讯网络可视化
#分组细胞类型相互作用数量弦图，以最大count做scale
count.max <- getMaxWeight(object.list, attribute = c("idents","count"))
weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
#### 初版
pdf(paste0(outdir,"Number_of_interactions_netplot_raw.pdf"),width = 12,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= 'T', edge.weight.max = count.max[2], edge.width.max = 12, title.name = title_Number[i])
	}
dev.off()
#分组细胞类型相互作用权重弦图，以最大weight做scale
pdf(paste0(outdir,"Weight_of_interactions_netplot_raw.pdf"),width = 12,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= 'T', edge.weight.max = weight.max[2], edge.width.max = 12, title.name = title_Weight[i])
	}
dev.off()
#### 无标签版
pdf(paste0(outdir,"Number_of_interactions_netplot_unlabel.pdf"),width = 12,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= labeledge, edge.weight.max = count.max[2], edge.width.max = 12, title.name = title_Number[i])
	}
dev.off()
#分组细胞类型相互作用权重弦图，以最大weight做scale
pdf(paste0(outdir,"Weight_of_interactions_netplot_unlabel.pdf"),width = 12,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= labeledge, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = title_Weight[i])
	}
dev.off()
#### 放大版
pdf(paste0(outdir,"Number_of_interactions_netplot_p.pdf"),width = 24,height = 14)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= labeledge, edge.weight.max = count.max[2], edge.width.max = 12, title.name = title_Number[i])
	}
dev.off()
#分组细胞类型相互作用权重弦图，以最大weight做scale
pdf(paste0(outdir,"Weight_of_interactions_netplot_p.pdf"),width = 24,height = 14)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= labeledge, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = title_Weight[i])
	}
dev.off()
#### 筛选版
selected_cell_type <- c("Activate memoryB","GCB","NaiveB","Microglia","Neutrophil_G1","Neutrophil_G2","Neutrophil_G3","Neutrophil_G4","Neutrophil_G5")
pdf(paste0(outdir,"Number_of_interactions_netplot_select.pdf"),width = 14,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	sub_matrix_count <- object.list[[i]]@net$count[selected_cell_type, selected_cell_type]
	netVisual_circle(sub_matrix_count, weight.scale = T, label.edge= "F", edge.width.max = 12, title.name = title_Number[i])
	}
dev.off()
#分组细胞类型相互作用权重弦图，以最大weight做scale
pdf(paste0(outdir,"Weight_of_interactions_netplot_select.pdf"),width = 14,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	sub_matrix_weight <- object.list[[i]]@net$weight[selected_cell_type, selected_cell_type]
	netVisual_circle(sub_matrix_weight, weight.scale = T, label.edge= "F", edge.width.max = 12, title.name = title_Weight[i])
	}
dev.off()

#### 筛选版贝壳图
selected_cell_type <- c("Activate memoryB","GCB","NaiveB","Microglia","Neutrophil_G1","Neutrophil_G2","Neutrophil_G3","Neutrophil_G4","Neutrophil_G5")
pdf(paste0(outdir,"Number_of_interactions_circle_select.pdf"),width = 14,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	sub_matrix_count <- object.list[[i]]@net$count[selected_cell_type, selected_cell_type]
	mat_count2 = matrix(0,nrow = nrow(sub_matrix_count),ncol = ncol(sub_matrix_count),dimnames = dimnames(sub_matrix_count))
	mat_count2[c(3,4),] = sub_matrix_count[c(3,4),]
	netVisual_circle(mat_count2, weight.scale = T, label.edge= "F", edge.width.max = 12, title.name = title_Number[i])
	}
dev.off()
#分组细胞类型相互作用权重弦图，以最大weight做scale
pdf(paste0(outdir,"Weight_of_interactions_circle_select.pdf"),width = 14,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	sub_matrix_weight <- object.list[[i]]@net$weight[selected_cell_type, selected_cell_type]
	mat_weight2 = matrix(0,nrow = nrow(sub_matrix_weight),ncol = ncol(sub_matrix_weight),dimnames = dimnames(sub_matrix_weight))
	mat_weight2[c(3,4),] = sub_matrix_weight[c(3,4),]
	netVisual_circle(mat_weight2, weight.scale = T, label.edge= "F", edge.width.max = 12, title.name = title_Weight[i])
	}
dev.off()

#### 筛选版弦图
pathways.show <- c('APP')
pdf(paste0(outdir,"Number_of_interactions_chord_select_p.pdf"),width = 18,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_aggregate(object.list[[i]],signaling = pathways.show,targets.use = 12, remove.isolate = TRUE,layout = "chord", label.edge= "F", signaling.name = paste("all", names(object.list)[i]))
	}
dev.off()
pdf(paste0(outdir,"Weight_of_interactions_chord_select_p.pdf"),width = 18,height = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
	netVisual_aggregate(object.list[[i]],signaling = pathways.show,targets.use = 12,remove.isolate = TRUE, layout = "chord", label.edge= "F",signaling.name = paste("all", names(object.list)[i]))
	}
dev.off()

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
pdf(paste0(outdir,"interactions_chord_gene-OVX_CNS-defect.pdf"),width = 10,height = 10)
netVisual_chord_gene(object.list[[1]],targets.use = 12, signaling = c("APP"),legend.pos.x = 8)
dev.off()
pdf(paste0(outdir,"interactions_chord_gene-sham_CNS-defect.pdf"),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]],targets.use = 12, signaling = c("APP"),legend.pos.x = 8)
dev.off()

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
levels(object.list[[1]]@idents)
pdf(paste0(outdir,"interactions_chord_gene_pairL-OVX_CNS-defect.pdf"),width = 10,height = 10)
netVisual_chord_gene(object.list[[1]], sources.use = 11,targets.use = 12, slot.name = "netP", legend.pos.x = 10)
dev.off()
pdf(paste0(outdir,"interactions_chord_gene_pairL-sham_CNS-defect.pdf"),width = 10,height = 10)
netVisual_chord_gene(object.list[[2]], sources.use = 11,targets.use = 12, slot.name = "netP", legend.pos.x = 10)
dev.off()

pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
    file_name <- pdf_files[i]
    pdf_image <- pdf_render_page(file_name, dpi = 300)
    writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}
