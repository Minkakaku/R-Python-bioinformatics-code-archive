library(CellChat)



metadata = read.csv(file = '/home/zwl/Project/FWSC20241362/09.require3/3.ALL_celltype_metadata.csv',row.names = 1)
load('/home/zwl/Project/FWSC20241362/09.require3/cellchat/cellchat_object.RData')

target_bar = metadata %>% 
  dplyr::filter(cell_type %in% c("T cells","B cells","NK cells","Plasma",
                                 "Fibroblasts","Epithelial","Endothelial"))
target_cell = subsetCellChat(cellchat,cells.use = rownames(target_bar))

p1 = netVisual_diffInteraction(target_cell, weight.scale = T,comparison = c(2,1))
recor_save("sucelltype_Differential number of interactions.net",p1)

#Identify cell populations with significant changes in sending or receiving signals

num.link <- sapply(cellchat.objet.list.new, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(cellchat.objet.list.new)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cellchat.objet.list.new[[i]], title = names(cellchat.objet.list.new)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
g1 = patchwork::wrap_plots(plots = gg)
ggplot2::ggsave(filename = 'cell_populations_with_significant_changes.png',g1,w=12,h=6)

# 特异上调的下调的配受体关系对

features.name = pos.dataset
target_cell_dif = identifyOverExpressedGenes(target_cell,
                                          group.dataset = "datasets",
                                          pos.dataset = pos.dataset, 
                                          features.name = features.name,
                                          only.pos = FALSE,
                                          thresh.pc = 0.25, 
                                          thresh.fc = 1, 
                                          thresh.p = 0.05)
net = netMappingDEG(target_cell_dif, features.name = features.name)
net.up = subsetCommunication(target_cell_dif,
                             net = net,
                             datasets = pos.dataset,
                             # sources.use = 3,
                             ligand.logFC = 0.2,
                             receptor.logFC = NULL)
net.down = subsetCommunication(target_cell_dif,
                               net = net, 
                               datasets = ref.dataset,
                               # sources.use = 3,
                               ligand.logFC = -0.1,
                               receptor.logFC = -0.1)
gene.up = extractGeneSubsetFromPair(net.up, target_cell_dif)
gene.down = extractGeneSubsetFromPair(net.down, target_cell_dif)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(target_cell_dif, pairLR.use = pairLR.use.up,
                        # sources.use = 3, targets.use = c(6,11),
                        comparison = c(1, 2),  
                        angle.x = 90, 
                        remove.isolate = F,
                        title.name = paste0("Up-regulated signaling in ", pos.dataset))
gg1
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(target_cell, pairLR.use = pairLR.use.down,
                        sources.use = 3, targets.use = c(6,11),
                        comparison = c(1, 2), 
                        angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", names(cellchat.objet.list.new)[1]))
#> Comparing communications on a merged object
gg1 + gg2

#################################################single cell plot#######################################################
# mat <- cellchat@net$weight
# num_rows <- nrow(mat)
# num_cols <- 2
# num_plot_rows <- ceiling(num_rows / num_cols)
# 
# png_file_path <- "/data/work/cell_interaction/net_weight_separated.png"
# png(png_file_path, width = 800, height = num_plot_rows * 400)
# 
# par(mfrow = c(num_plot_rows, num_cols), xpd=TRUE)
# 
# for (i in 1:num_rows) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }
# 
# dev.off()


groupSize <- as.numeric(table(cellchat@idents[["EG"]]))
mat <- cellchat@net$EG$weight
num_rows <- nrow(mat)
num_cols <- 2
num_plot_rows <- ceiling(num_rows / num_cols)


par(mfrow = c(num_plot_rows, num_cols), xpd=TRUE)

mat
pdf('./Immu_T_cell.pdf', width = 10, height = 10)

for (i in c(17,25)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()



