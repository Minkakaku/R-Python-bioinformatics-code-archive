.libPaths('/mnt/NexData_file/home/wyt/software/minicoda/envs/R_envs/lib/R/library')
options(stringsAsFactors = FALSE)
library(CellChat)
library(Seurat)
library(tidyverse)
library(viridis)
library(argparser)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(NMF)
library(ggalluvial)

parser = argparse::ArgumentParser(description="celllchat contrast analysis")
parser$add_argument('--pos.dataset')
parser$add_argument('--other.dataset')
parser$add_argument('--pos_rds')
parser$add_argument('--other_rds')
parser$add_argument('--output')

args = parser$parse_args()
pos.dataset = args$pos.dataset
other.dataset = args$other.dataset
pos_rds  = args$pos_rds
other_rds  = args$other_rds
output = args$output

cols <- c(brewer.pal(9,"Set1")[c(1:5,7,6,8,9)],brewer.pal(9,"Pastel1")[c(7,2,3,5,1)],brewer.pal(12,"Set3")[c(-2,-4,-6,-9,-11)],brewer.pal(8,"Set2")[c(3,4,7)],brewer.pal(8,"Pastel2")[8],brewer.pal(12,"Paired")[c(5,6,9,10,12)])

source("/home/wyt/test/protocol/cellchat_Spatial/ce.R")
dir.create(paste0(output,"/",pos.dataset,"_",other.dataset))
setwd(paste0(output,"/",pos.dataset,"_",other.dataset))
cellchat.other <- readRDS(other_rds)
cellchat.pos <- readRDS(pos_rds)
library(purrr)

cellchat.pos_filtered <- cellchat_filtered(cellchat.1 = cellchat.pos,cellchat.2 =cellchat.other)
cellchat.other_filtered <- cellchat_filtered(cellchat.1 = cellchat.other,cellchat.2 =cellchat.pos)
object.list <- setNames(
  list(cellchat.pos_filtered, cellchat.other_filtered),
  c(pos.dataset, other.dataset)
)
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
saveRDS(cellchat,paste0(output,"/",pos.dataset,"_",other.dataset,"/merge_Cellchat.rds"))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf(paste0(pos.dataset,"_",other.dataset,"_internation_barplot.pdf"),12,6)
gg1 + gg2
dev.off()

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf(paste0(pos.dataset,"_",other.dataset,"_differential_internation_heatmap.pdf"),12,8)
gg1 
gg2
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
pdf(paste0(pos.dataset,"_",other.dataset,"_2.pathway_barplot.pdf"),12,8)
gg1 + gg2
dev.off()
library(ComplexHeatmap)
library(patchwork)

rank_results <- rankNet(cellchat, 
                        mode = "comparison", 
                        stacked = T, 
                        do.stat = TRUE)
stat_info_data <- rank_results$data %>% filter(group == pos.dataset & pvalues < 0.05) 
stat_info_data$group <- ifelse(
  is.infinite(stat_info_data$contribution.relative.1),  paste0("single_",other.dataset),
  ifelse(
    stat_info_data$contribution.relative.1 == 0,  paste0("single_",pos.dataset),
    ifelse(
      stat_info_data$contribution.relative.1 >= 1, paste0("High_",other.dataset), paste0("High_",pos.dataset)
    )
  )
)
write.csv(stat_info_data,"stat_info_data.csv",row.names = T)
####################
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 15, height = 30, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 15, height = 30, color.heatmap = "OrRd")
pdf(paste0(pos.dataset,"_",other.dataset,"_2.all_pathway.pdf"),15,15)
ht1
ht2
dev.off()

pdf(paste0(pos.dataset,"_",other.dataset,"_target_up_down_L-R.pdf"),12,12)
netVisual_bubble(cellchat,sources.use = c(1:3), targets.use = 1,  comparison = c(1, 2),angle.x = 45)
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf(paste0(pos.dataset,"_",other.dataset,"_netVisual_circle.pdf"),20,30)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

mat <- object.list[[1]]@net$count
par(mfrow = c(3,4), xpd=TRUE)
pdf(paste0(pos.dataset,"_weigth2.pdf"),12,12)
for (i in 1:nrow(mat)) {
  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(mat2, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = "Number of interactions - pos")
  
}
dev.off()

mat <- object.list[[2]]@net$count
par(mfrow = c(3,4), xpd=TRUE)
pdf(paste0(other.dataset,"_weigth2.pdf"),12,12)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = "Number of interactions - other")
}
dev.off()
par(mfrow = c(1, 1))


library(ggplot2)
all_plots <- list()
pdf("netVisual_diffInteraction.pdf",12,12)
for (i in 1:nrow(mat)) {
  p <- netVisual_diffInteraction(
    object = cellchat,
    comparison = c(2,1),
    sources.use = i,
    label.edge = TRUE,
    edge.label.cex = 1,
    edge.width.max = 2
  )
  print(p)
  all_plots[[i]] <- p
}
dev.off()

i = 1
pos_pathway <- subset(stat_info_data,group !=  paste0("single_",other.dataset))$name
pdf(paste0(pos.dataset,"_",other.dataset,pos.dataset,"_pathway_heatmap.pdf"),15,15)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pos_pathway, title = names(object.list)[i], width = 14, height =15, color.heatmap = "BuGn")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pos_pathway, title = names(object.list)[i], width = 14, height = 15, color.heatmap = "OrRd")
ht1 + ht2
dev.off()

other_pathway <- subset(stat_info_data,group != paste0("single_",pos.dataset))$name
pdf(paste0(pos.dataset,"_",other.dataset,other.dataset,"_pathway_heatmap.pdf"),15,15)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = other_pathway, title = names(object.list)[i+1], width = 14, height = 15, color.heatmap = "BuGn")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling =  other_pathway, title = names(object.list)[i+1], width = 14, height = 15, color.heatmap = "OrRd")
ht1 + ht2
dev.off()

###################plot########################
library(igraph)
pos.net.source <- subsetCommunication(cellchat.pos) %>% 
  select(c("source","target","interaction_name","ligand","receptor","prob","pathway_name")) %>% 
  filter(pathway_name %in% stat_info_data$name) %>%
  mutate(group =pos.dataset)

other.net.source <- subsetCommunication(cellchat.other) %>% 
  select(c("source","target","interaction_name","ligand","receptor","prob","pathway_name")) %>% 
  filter(pathway_name %in% stat_info_data$name) %>%
  mutate(group =other.dataset) 

merged_data <- merge(
  pos.net.source, 
  other.net.source, 
  by = c("source", "target", "interaction_name","ligand","receptor","pathway_name"),
  suffixes = c(paste0("_",pos.dataset),paste0("_",other.dataset)),
  all = TRUE 
)

merged_data$prob_diff <- ifelse(
  is.na(merged_data[[paste0("prob_",pos.dataset)]]) & !is.na(merged_data[[paste0("prob_",other.dataset)]]),
  -merged_data[[paste0("prob_",other.dataset)]],  
  ifelse(
    !is.na(merged_data[[paste0("prob_",pos.dataset)]]) & is.na(merged_data[[paste0("prob_",other.dataset)]]),
    merged_data[[paste0("prob_",pos.dataset)]],
    merged_data[[paste0("prob_",pos.dataset)]] - merged_data[[paste0("prob_",other.dataset)]]
  )
)
write.csv(merged_data,paste0(pos.dataset,"_",other.dataset,"_merged_data.csv"))
dir.create(paste0(output,"/",pos.dataset,"_",other.dataset,"/circos_celltype_plot"), recursive = TRUE, 
           showWarnings = FALSE)

for (i in unique(merged_data$source)) {
  png(paste0(output,"/",pos.dataset,"_",other.dataset,"/circos_celltype_plot/",i,"_circos_celltype_plot.png"),
      width = 8, height = 10, units = "in", res = 300)
  result <- circos_celltype_plot(data=merged_data,source.use = i,right_value2=-0.5,right_value1 =-0.15,bottom_value=0.2,title_value = -7.5)
  dev.off()
}
for (i in unique(merged_data$source)) {
  pdf(paste0(output,"/",pos.dataset,"_",other.dataset,"/circos_celltype_plot/",i,"_circos_celltype_plot.pdf"),
      width = 8, height = 10,onefile = FALSE)
  result <- circos_celltype_plot(data=merged_data,source.use = i,right_value2=-0.5,right_value1 =-0.15,bottom_value=0.2,title_value = -7.5)
  dev.off()
}

dir.create(paste0(output,"/",pos.dataset,"_",other.dataset,"/circos_gene_plot"), recursive = TRUE, 
           showWarnings = FALSE)
for (i in unique(merged_data$source)) {
  png(paste0(output,"/",pos.dataset,"_",other.dataset,"/circos_gene_plot/",i,"_circos_gene_plot.png"),
      width = 8, height = 8, units = "in", res = 300)
  circos_gene_plot(source.use = i,title_value = -1,bottom_value = 0,right_value = -0.3,num = 10,data = merged_data)
  dev.off()
}

for (i in unique(merged_data$source)) {
  pdf(paste0(output,"/",pos.dataset,"_",other.dataset,"/circos_gene_plot/",i,"_circos_gene_plot.pdf"),
      width = 8, height = 8)
  circos_gene_plot(source.use = i,title_value = -1,bottom_value = 0,right_value = -0.3,num = 10,data = merged_data)
  dev.off()
}

dir.create(paste0(output,"/",pos.dataset,"_",other.dataset,"/pathway_name"), recursive = TRUE, 
           showWarnings = FALSE)
for (j in unique(merged_data$pathway_name)) {
  png(paste0(output,"/",pos.dataset,"_",other.dataset,"/pathway_name/",j,"_circos_gene_plot.png"),
      width = 8, height = 8, units = "in", res = 300)
  circos_gene_plot(signaling =j,title_value = -1,bottom_value = 0,right_value = -0.3,num = 10,data = merged_data)
  dev.off()
}

for (j in unique(merged_data$pathway_name)) {
  pdf(paste0(output,"/",pos.dataset,"_",other.dataset,"/pathway_name/",j,"_circos_gene_plot.pdf"),
      width = 8, height = 8)
  circos_gene_plot(signaling =j,title_value = -1,bottom_value = 0,right_value = -0.3,num = 10,data = merged_data)
  dev.off()
}

dir.create(paste0(output,"/",pos.dataset,"_",other.dataset,"/barplot_celltype"), recursive = TRUE, 
           showWarnings = FALSE)
for (z in unique(c(merged_data$source,merged_data$target))) {
  p <- barplot_celltype(cellchat.pos = cellchat.pos,cellchat.other = cellchat.other,celltype_label = z)
  ggsave(p,filename = paste0(output,"/",pos.dataset,"_",other.dataset,"/barplot_celltype/",z,"_barplot_celltype_plot.pdf"),width = 12,height = 6)
}
for (z in unique(c(merged_data$source,merged_data$target))) {
  p <- barplot_celltype(cellchat.pos = cellchat.pos,cellchat.other = cellchat.other,celltype_label = z)
  ggsave(p,filename = paste0(output,"/",pos.dataset,"_",other.dataset,"/barplot_celltype/",z,"_barplot_celltype_plot.png"),width = 12,height = 6)
}
##########################################################
num <- cellchat@idents$joint %>% unique() %>% length()
pdf(paste0(output,"/",pos.dataset,"_",other.dataset,"/LBH+basal_KC_source_up_down_L-R.pdf"),20,8)
netVisual_bubble(cellchat,sources.use =c("LBH+basal_KC"),targets.use =c(1:num),comparison = c(1, 2),angle.x = 90)
dev.off()

pdf(paste0(output,"/",pos.dataset,"_",other.dataset,"/LBH+basal_KC_target_up_down_L-R.pdf"),20,8)
netVisual_bubble(cellchat,sources.use =c(1:num),targets.use =c("LBH+basal_KC"),comparison = c(1, 2),angle.x = 90)
dev.off()

pdf(paste0(output,"/",pos.dataset,"_",other.dataset,"/SFRP1+Fibroblast_source_up_down_L-R.pdf"),20,8)
netVisual_bubble(cellchat,sources.use =c("SFRP1+Fibroblast"),targets.use =c(1:num),comparison = c(1, 2),angle.x = 90)
dev.off()

pdf(paste0(output,"/",pos.dataset,"_",other.dataset,"/SFRP1+Fibroblast_target_up_down_L-R.pdf"),20,8)
netVisual_bubble(cellchat,sources.use =c(1:num),targets.use =c("SFRP1+Fibroblast"),comparison = c(1, 2),angle.x = 90)
dev.off()

pdf(paste0(output,"/",pos.dataset,"_",other.dataset,"/LECs_source_up_down_L-R.pdf"),20,8)
netVisual_bubble(cellchat,sources.use =c("LECs"),targets.use =c(1:num),comparison = c(1, 2),angle.x = 90)
dev.off()

pdf(paste0(output,"/",pos.dataset,"_",other.dataset,"/LECs_target_up_down_L-R.pdf"),20,8)
netVisual_bubble(cellchat,sources.use =c(1:num),targets.use =c("LECs"),comparison = c(1, 2),angle.x = 90)
dev.off()
