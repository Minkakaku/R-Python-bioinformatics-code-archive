setwd("./")
getwd()
ls()

ptm = Sys.time()
library(Seurat)
library(tidyverse)
library(RColorBrewer)

library(NMF)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

cellchat.Mod_96 <- readRDS("/mnt/NexData_file/home/hhq/project/Xenium/FWSP20255020-02_mouse_lung_xenium/cellchat/Mod_96_cellchat_final.rds")
cellchat.Mod_94 <- readRDS("/mnt/NexData_file/home/hhq/project/Xenium/FWSP20255020-02_mouse_lung_xenium/cellchat/Mod_94/Mod94_cellchat_final.rds")
cellchat.Mod_94s <- readRDS("/mnt/NexData_file/home/hhq/project/Xenium/FWSP20255020-02_mouse_lung_xenium/cellchat/Mod_94s/Mod_94s_cellchat_final.rds")
cellchat.Mod_95 <- readRDS("/mnt/NexData_file/home/hhq/project/Xenium/FWSP20255020-02_mouse_lung_xenium/cellchat/Mod_95/Mod95_cellchat_final.rds")

cellchat.Tcm_28 <- readRDS("/mnt/NexData_file/home/hhq/project/Xenium/FWSP20255020-02_mouse_lung_xenium/cellchat/Tcm_28/Tcm_28_cellchat_final.rds")
cellchat.Tcm_28s <- readRDS("/mnt/NexData_file/home/hhq/project/Xenium/FWSP20255020-02_mouse_lung_xenium/cellchat/Tcm_28s/Tcm_28s_cellchat_final.rds")
cellchat.Tcm_29 <- readRDS("/mnt/NexData_file/home/hhq/project/Xenium/FWSP20255020-02_mouse_lung_xenium/cellchat/Tcm_29/Tcm_29_cellchat_final.rds")
cellchat.Tcm_30 <- readRDS("/mnt/NexData_file/home/hhq/project/Xenium/FWSP20255020-02_mouse_lung_xenium/cellchat/Tcm_30/Tcm_30_cellchat_final.rds")

# group.new = levels(cellchat.LS@idents)
# cellchat.IV <- liftCellChat(cellchat.IV, group.new)

object.list <- list(Mod_96 = cellchat.Mod_96, Mod_94 = cellchat.Mod_94, Mod_94s=cellchat.Mod_94s, Mod_95 = cellchat.Mod_95, Tcm_28 = cellchat.Tcm_28, Tcm_28s = cellchat.Tcm_28s, Tcm_29 = cellchat.Tcm_29, Tcm_30 = cellchat.Tcm_30)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
#> An object of class CellChat created from a merged object with multiple datasets 
#>  593 signaling genes.
#>  7563 cells. 
#> CellChat analysis of single cell RNA-seq data!
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 2.763455

# Users can now export the merged CellChat object and the list of the two separate objects for later use
save(object.list, file = "cellchat_object.list_humanlung_I_II_IV.RData")
save(cellchat, file = "cellchat_merged_humanLung_I_II_IV.RData")

ptm = Sys.time()
pdf("Barplot_compare_weight_number.pdf")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2
dev.off()

pdf("diff_circle_tumor_source.pdf", width=20, height=16)
par(mfrow = c(1,2), xpd=TRUE)
n1 <- netVisual_diffInteraction(cellchat, weight.scale = T, sources.use='Tumor', vertex.label.cex = 0.7)
n2 <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", sources.use='Tumor', vertex.label.cex = 0.7)
dev.off()

pdf("diff_circle_tumor_source_.pdf", width=30, height=16)
par(mfrow = c(1,3), xpd=TRUE, cex.main = 2)
n1 <- netVisual_diffInteraction(cellchat, weight.scale = T, sources.use='Tumor', vertex.label.cex = 1.5, comparison = c(1,2), label.edge = T, vertex.size.max = 10, edge.label.cex = 1.5, title.name = 'II vs. I')
n2 <- netVisual_diffInteraction(cellchat, weight.scale = T, sources.use='Tumor', vertex.label.cex = 1.5, comparison = c(1,3), label.edge = T, vertex.size.max = 10, edge.label.cex = 1.5, title.name = 'IV vs. I')
n3 <- netVisual_diffInteraction(cellchat, weight.scale = T, sources.use='Tumor', vertex.label.cex = 1.5, comparison = c(2,3), label.edge = T, vertex.size.max = 10, edge.label.cex = 1.5, title.name = 'IV vs. II')
dev.off()

pdf("diff_circle_tumor_target_.pdf", width=30, height=16)
par(mfrow = c(1,3), xpd=TRUE, cex.main = 2)
n1 <- netVisual_diffInteraction(cellchat, weight.scale = T, targets.use ='Tumor', vertex.label.cex = 1.5, comparison = c(1,2), label.edge = T, vertex.size.max = 10, edge.label.cex = 1.5, title.name = 'II vs. I')
n2 <- netVisual_diffInteraction(cellchat, weight.scale = T, targets.use='Tumor', vertex.label.cex = 1.5, comparison = c(1,3), label.edge = T, vertex.size.max = 10, edge.label.cex = 1.5, title.name = 'IV vs. I')
n3 <- netVisual_diffInteraction(cellchat, weight.scale = T, targets.use='Tumor', vertex.label.cex = 1.5, comparison = c(2,3), label.edge = T, vertex.size.max = 10, edge.label.cex = 1.5, title.name = 'IV vs. II')
dev.off()

pdf("heatmap_diff_counts.pdf", onefile = T)
gg1 <- netVisual_heatmap(cellchat)
gg1
dev.off()
#> Do heatmap based on a merged object

pdf("heatmap_diff_weight.pdf", onefile = T)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg2
dev.off()

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(unlist(num.link)), max(unlist(num.link))) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

pdf("signal_change_scatter_II_vs_I.pdf", width=7, height=5, onefile = T)
for (i in as.character(unique(cellchat@meta$labels))){
    gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = i)
    print(gg1)
}
dev.off()

pdf("signal_change_scatter_IV_vs_II.pdf", width=7, height=5, onefile = T)
for (i in as.character(unique(cellchat@meta$labels))){
    gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = i, comparison = c(2, 3))
    print(gg1)
}
dev.off()

pdf("tumor_signal_change.pdf", width=15, height=15)
netAnalysis_signalingChanges_scatter(cellchat, idents.use = 'Tumor', comparison = c(1,2))
dev.off()

pdf("tumor_signal_change_IV_vs_I.pdf", width=15, height=15)
netAnalysis_signalingChanges_scatter(cellchat, idents.use = 'Tumor', comparison = c(1,3))
dev.off()

# ptm = Sys.time()

# cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
# #> Compute signaling network similarity for datasets 1 2
# cellchat <- netEmbedding(cellchat, type = "functional")
# #> Manifold learning of the signaling networks for datasets 1 2
# cellchat <- netClustering(cellchat, type = "functional")
# #> Classification learning of the signaling networks for datasets 1 2
# # Visualization in 2D-space
# netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
# #> 2D visualization of signaling networks from datasets 1 2

# cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
# cellchat <- netEmbedding(cellchat, type = "structural")
# cellchat <- netClustering(cellchat, type = "structural")
# # Visualization in 2D-space
# netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
# netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

# rankSimilarity(cellchat, type = "functional")
# #> Compute the distance of signaling networks between datasets 1 2

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = c('Tumor'), stacked = T, do.stat = TRUE, title = 'Target: Tumor')
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = c('Tumor'), stacked = F, do.stat = TRUE, title = 'Target: Tumor')
pdf('rankNet_target_tumor.pdf')
gg1 + gg2
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = c('Tumor'), stacked = T, do.stat = TRUE, title = 'Target: Tumor', comparison = c(1, 3))
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = c('Tumor'), stacked = F, do.stat = TRUE, title = 'Target: Tumor', comparison = c(1, 3))
pdf('rankNet_target_tumor_IV_I.pdf')
gg1 + gg2
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = c('Tumor'), stacked = T, do.stat = TRUE, title = 'Target: Tumor', comparison = c(2, 3))
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = c('Tumor'), stacked = F, do.stat = TRUE, title = 'Target: Tumor', comparison = c(2, 3))
pdf('rankNet_target_tumor_IV_II.pdf')
gg1 + gg2
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = 'Tumor', targets.use = NULL, stacked = T, do.stat = TRUE, title = 'Source: Tumor')
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = 'Tumor', targets.use = NULL, stacked = F, do.stat = TRUE, title = 'Source: Tumor')
pdf('rankNet_source_tumor.pdf')
gg1 + gg2
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = 'Tumor', targets.use = NULL, stacked = T, do.stat = TRUE, title = 'Source: Tumor', comparison = c(1,3))
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = 'Tumor', targets.use = NULL, stacked = F, do.stat = TRUE, title = 'Source: Tumor', comparison = c(1,3))
pdf('rankNet_source_tumor_IV_I.pdf')
gg1 + gg2
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = 'Tumor', targets.use = NULL, stacked = T, do.stat = TRUE, title = 'Source: Tumor', comparison = c(2,3))
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = 'Tumor', targets.use = NULL, stacked = F, do.stat = TRUE, title = 'Source: Tumor', comparison = c(2,3))
pdf('rankNet_source_tumor_IV_II.pdf')
gg1 + gg2
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
pdf('rankNet_total.pdf', height=14, width=8)
gg1 + gg2
dev.off()

i=1
library(ComplexHeatmap)
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1

ht2

netVisual_bubble(cellchat, sources.use = 'Tumor', targets.use = NULL,  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat, sources.use = 'Tumor', targets.use = NULL,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in II", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 'Tumor', targets.use = NULL,  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in I", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("diff_bubble_source_tumor_I_vs_II.pdf", width=15)
gg1 + gg2
dev.off()

gg1 <- netVisual_bubble(cellchat, sources.use = 'Tumor', targets.use = NULL,  comparison = c(1, 3), max.dataset = 3, title.name = "Increased signaling in IV", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 'Tumor', targets.use = NULL,  comparison = c(1, 3), max.dataset = 1, title.name = "Increased signaling in I", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("diff_bubble_source_tumor_I_vs_IV.pdf", width=15)
gg1 + gg2
dev.off()

gg1 <- netVisual_bubble(cellchat, sources.use = 'Tumor', targets.use = NULL,  comparison = c(2, 3), max.dataset = 3, title.name = "Increased signaling in IV", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 'Tumor', targets.use = NULL,  comparison = c(2, 3), max.dataset = 2, title.name = "Increased signaling in II", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("diff_bubble_source_tumor_II_vs_IV.pdf", width=22)
gg1 + gg2
dev.off()

gg1 <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = 'Tumor',  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in II", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = 'Tumor',  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in I", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("diff_bubble_target_tumor_I_vs_II.pdf", width=15, height=8)
gg1 + gg2
dev.off()

gg1 <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = 'Tumor',  comparison = c(1, 3), max.dataset = 3, title.name = "Increased signaling in IV", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = 'Tumor',  comparison = c(1, 3), max.dataset = 1, title.name = "Increased signaling in I", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("diff_bubble_target_tumor_I_vs_IV.pdf", width=15, height=8)
gg1 + gg2
dev.off()

gg1 <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = 'Tumor',  comparison = c(2, 3), max.dataset = 3, title.name = "Increased signaling in IV", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = 'Tumor',  comparison = c(2, 3), max.dataset = 2, title.name = "Increased signaling in II", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
pdf("diff_bubble_target_tumor_II_vs_IV.pdf", width=22, height=8)
gg1 + gg2
dev.off()


