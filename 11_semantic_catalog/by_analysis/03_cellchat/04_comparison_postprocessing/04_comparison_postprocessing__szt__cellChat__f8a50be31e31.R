
library(patchwork)
library(CellChat)
library(ggplotify)

my_theme2 <- theme_bw() + theme(axis.text.x = element_text(hjust = 0.5, size = 16, color = "black"),
                                axis.text.y = element_text(hjust = 0.5, size = 16, color = "black"),
                                axis.title = element_text(size = 18, color = "black"),
                                plot.title = element_text(size = 22, color = "black"),
                                axis.ticks = element_line(linewidth = 0.2, color="black"),
                                axis.ticks.length = unit(0.2, "cm"),
                                strip.text = element_text(size = 20),
                                panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank())


#加载每个数据集的cellchat对象，然后合并在一起
cellchat.activate <- readRDS(file = "activate/cellchat_activate.rds")
cellchat.non_activate <- readRDS(file = "non_activate/cellchat_non_activate.rds")

object.list <- list(activate = cellchat.activate, non_activate = cellchat.non_activate)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
################################################################################

#第一部分：预测细胞通信的一般原理

################################################################################
#为了回答细胞-细胞通信是否增强的问题，CellChat比较了来自不同生物条件的细胞通信网络的相互作用数量和强度。
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), color.use = c("#E73F74", "#219EBC"))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight", color.use = c("#E73F74", "#219EBC"))
gg <- gg1 + my_theme2 + gg2 + my_theme2
gg
ggsave(plot = gg, filename = "compareInteractions.pdf", width = 7, height = 5, device = "pdf")
################################################################################
#不同细胞群之间的相互作用数量或强度的差异
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg <- gg1 + gg2
gg

pdf(file = "netVisual_heatmap.pdf", width = 10, height = 7)
print(gg)
dev.off()
################################################################################
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf(file = "netVisual_circle.pdf", width = 9, height = 5)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()
################################################################################
#比较 2D 空间中的主要来源和目标
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
  #+geom_hline(aes(yintercept = 1))+geom_vline(aes(xintercept = 1))
}
pdf(file = "netAnalysis_signalingRole_scatter.pdf", width = 12, height = 6)
patchwork::wrap_plots(plots = gg)
dev.off()
########################################################## ######################

#第二部分：识别保守和环境特异的信号通路

################################################################################
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional", do.parallel = F)
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
#-------------------------------------------------------------------------------
#计算和可视化通路距离
pdf(file = "rankSimilarity.pdf", height = 10, width = 5)
rankSimilarity(cellchat, type = "functional")
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, #sources.use = 4, #targets.use = 1,
               color.use = c("#E73F74", "#219EBC"), cutoff.pvalue = 0.01, do.flip = F)

#gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, color.use = c("#E73F74", "#219EBC"))
gg <- gg1 #+ gg2
pdf(file = "rankNet.pdf", width = 12, height = 4)
print(gg)
dev.off()
#-------------------------------------------------------------------------------
#比较与每个细胞群相关的传出（或传入）信号
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 14)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 14)
pdf(file = "outgoing.pdf", width = 8, height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 14, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 14, color.heatmap = "GnBu")
pdf(file = "incoming.pdf", width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 14, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 14, color.heatmap = "OrRd")
pdf(file = "all.pdf", width = 8, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
################################################################################

#第三部分：识别上调和下调的信号配体对

################################################################################
netVisual_bubble(cellchat, 
                 sources.use = c(1:13), targets.use = c(1:13),  
                 comparison = c(1, 2), angle.x = 45)


gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:13),  comparison = c(1, 2), max.dataset = 2, title.name = "Iactivatereased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:13),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
################################################################################

#第四部分：使用层次结构图、圆图或和弦图可视比较细胞-细胞通信

################################################################################
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(file = paste0(pathways.show, "_netVisual_aggregate.pdf"), width = 8, height = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 5, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
################################################################################

#第五部分：比较不同数据集之间的信号基因表达分布

################################################################################
cellchat@meta$datasets <- factor(cellchat@meta$datasets, levels = c("activate", "non_activate")) # set factor level
p <- plotGeneExpression(cellchat, signaling = pathways.show, type = c("violin"), features = "",
                        split.by = "datasets", colors.ggplot = T, color.use = c("#FFE76F", "#91B822"))
ggsave(plot = p, filename = paste0(pathways.show, ".pdf"), width = 12, height = 16)

saveRDS(cellchat, file = "cellchat_comparisonAnalysis_activate_vs_non_activate.rds")
