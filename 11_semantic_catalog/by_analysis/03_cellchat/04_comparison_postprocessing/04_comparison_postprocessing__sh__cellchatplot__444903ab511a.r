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
source("/mnt/NexData_file/home/sh/Projects/SC2510235895/10.Mission/script/cellchat_f.r")

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

load("/mnt/NexData_file/home/sh/Projects/SC2510235895/10.Mission/cellchat_sub/tmp_object.RData")
# cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/sh/Projects/SC2510235895/10.Mission/cellchat_sub_plot/"
pathway_outdir<-outdir
select_pathway<-c('CCL','GALECTIN','CD137')
for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(object.list[[1]], signaling = k, geneLR.return = FALSE)
	vertex_receiver <- vertex_cal(object.list[[1]],k)
	vertex_receiver <- c(1:9)
	# pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-centrality_heatmap.pdf"),width = 8,height= 3)
	# gg <- netAnalysis_signalingRole_network(object.list[[1]], signaling = k, width = 14, height = 2.5, font.size = 10)
	# dev.off()
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-centrality_heatmap.pdf"),width = 8,height= 3)
	gg <- netAnalysis_signalingRole_network_mod(object.list[[1]], signaling = k, width = 14, height = 2.5, font.size = 10)
	dev.off()

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-contribution_barplot.pdf"),width = 5,height= 5)
	print(netAnalysis_contribution(object.list[[1]], signaling = k))
	dev.off()

	# pdf(paste0(pathway_outdir,k,"-all-",pairLR_top[9,],"-top_hierarchy.pdf"),width = 12,height= 8)
	# netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top[9,], vertex.receiver = vertex_receiver, layout = "hierarchy")
	# dev.off()
	# pdf(paste0(pathway_outdir,k,"-all-",pairLR_top[10,],"-top_hierarchy.pdf"),width = 12,height= 8)
	# netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top[10,], vertex.receiver = vertex_receiver, layout = "hierarchy")
	# dev.off()
	# pdf(paste0(pathway_outdir,k,"-all-",pairLR_top[14,],"-top_hierarchy.pdf"),width = 12,height= 8)
	# netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top[14,], vertex.receiver = vertex_receiver, layout = "hierarchy")
	# dev.off()
	pdf(paste0(pathway_outdir,k,"-all-",pairLR_top[1,],"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual_mod(object.list[[1]], signaling = k, pairLR.use = pairLR_top[1,], vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
	pdf(paste0(pathway_outdir,k,"-all-",pairLR_top[9,],"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual_mod(object.list[[1]], signaling = k, pairLR.use = pairLR_top[9,], vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
	pdf(paste0(pathway_outdir,k,"-all-",pairLR_top[10,],"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual_mod(object.list[[1]], signaling = k, pairLR.use = pairLR_top[10,], vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
	pdf(paste0(pathway_outdir,k,"-all-",pairLR_top[14,],"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual_mod(object.list[[1]], signaling = k, pairLR.use = pairLR_top[14,], vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}



#### 后面都没做
for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(cellchats, signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[1,]
	vertex_receiver <- c(6,4,12)
	pdf(paste0(pathway_outdir,k,"-all-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(cellchats, signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}
for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(cellchats, signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[1,]
	vertex_receiver <- c(6,4,12)
	pdf(paste0(pathway_outdir,k,"-all-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(cellchats, signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}

pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
	file_name <- pdf_files[i]
	pdf_image <- pdf_render_page(file_name, dpi = 300)
	writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
	}


data<-netVisual_bubble(object.list[[1]], sources.use = c(2:3), targets.use = c(1), angle.x = 45,return.data = TRUE)
data<-data$communication[data$communication$annotation=="Secreted Signaling",]
pdf(paste0(outdir,"netVisual_bubble_all.pdf"),width = 25,height= 6)
netVisual_bubble(object.list[[1]], sources.use = c(2:3), targets.use = c(1), angle.x = 45, remove.isolate = TRUE, dot.size.min = 3, dot.size.max = 12, font.size = 15, font.size.title = 15) + coord_flip()
dev.off()

signaling_use<-unique(data$pathway_name)
signaling_use<-c("CCL","CXCL","IL1","GALECTIN","SPP1","CD137")

pdf(paste0(outdir,"netVisual_bubble_sub.pdf"),width = 12,height= 4)
netVisual_bubble(object.list[[1]], sources.use = c(2:3), targets.use = c(1), angle.x = 45,signaling=signaling_use, remove.isolate = TRUE, dot.size.min = 3, dot.size.max = 12, font.size = 15, font.size.title = 15) + coord_flip()
dev.off()

pdf(paste0(outdir,"netVisual_bubble_Macrophages.pdf"),width = 8,height= 2)
netVisual_bubble(object.list[[1]], sources.use = c(3), targets.use = c(1), angle.x = 45,signaling=signaling_use, remove.isolate = TRUE, dot.size.min = 3, dot.size.max = 12, font.size = 15, font.size.title = 15) + coord_flip()
dev.off()
pdf(paste0(outdir,"netVisual_bubble_Fibroblasts.pdf"),width = 8,height= 2)
netVisual_bubble(object.list[[1]], sources.use = c(2), targets.use = c(1), angle.x = 45,signaling=signaling_use, remove.isolate = TRUE, dot.size.min = 3, dot.size.max = 12, font.size = 15, font.size.title = 15) + coord_flip()
dev.off()

pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
    file_name <- pdf_files[i]
		# file_name <- paste0(outdir,"netVisual_bubble.pdf")
    pdf_image <- pdf_render_page(file_name, dpi = 300)
    writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}



load("/mnt/NexData_file/home/sh/Projects/SC2510235895/09.Mission/cellchat_sub_group/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list),merge.data = TRUE)
outdir<-"/mnt/NexData_file/home/sh/Projects/SC2510235895/09.Mission/cellchat_sub_group/"

# 批量处理并合并
combined_data <- do.call(rbind, lapply(seq_along(object.list), function(i) {
  data <- netVisual_bubble(object.list[[i]], angle.x = 45, return.data = TRUE)
  data_filtered <- data$communication[data$communication$annotation == "Secreted Signaling", ]
  data_filtered <- data_filtered[data_filtered$pathway_name %in% signaling_use, ]
  
  if (nrow(data_filtered) > 0) {
    data_filtered$group <- names(object.list)[i]
  }
  return(data_filtered)
}))
combined_data<-combined_data[combined_data$source.target %in% c("Macrophages -> CD8T","Fibroblasts -> CD8T"),]
combined_data1<-combined_data[combined_data$source=="Macrophages",]

df<-combined_data1
df$group<-factor(df$group,levels=c("asham","aday1","aday3","aday7","aday14"))
g <- ggplot(df, aes(x = group, y = interaction_name_2,
        color = prob, size = pval)) + geom_point(pch = 16) +
        theme_linedraw() + theme(panel.grid.major = element_blank()) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1,
            vjust = 1), axis.title.x = element_blank(),
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
dot.size.max = max(df$pval)
dot.size.min = min(df$pval)

g <- g + scale_radius(range = c(dot.size.min*3, dot.size.max*3),
    breaks = sort(unique(df$pval)), labels = names(values)[values %in%
        sort(unique(df$pval))], name = "p-value")

color.use <- RColorBrewer::brewer.pal(n = 10, name = "Spectral")
color.use <- rev(color.use)

g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),
    na.value = "white", limits = c(quantile(df$prob,
        0, na.rm = T), quantile(df$prob, 1, na.rm = T)),
    breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob,
        1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5,
     title = "Commun. Prob."))

g <- g + theme(text = element_text(size = 10), plot.title = element_text(size = 10)) +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8))

if (length(unique(df$source.target)) > 1) {
    g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) -
        0.5, 1), lwd = 0.1, colour = "grey90")
}
    if (length(unique(df$interaction_name_2)) > 1) {
        g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) -
            0.5, 1), lwd = 0.1, colour = "grey90")
    }

# comparison = c(1,2,3,4,5)
# dataset.name <- names(cellchat@net)
# group.names0 <- df$source.target

# xintercept = seq(0.5 + length(dataset.name[comparison]),
#     length(group.names0) * length(dataset.name[comparison]),
#     by = length(dataset.name[comparison]))
# g <- g + geom_vline(xintercept = xintercept, linetype = "dashed",
#     color = "grey60", size = 0.2)

# group <- 1:length(comparison)
# names(group) <- dataset.name[comparison]
# color <- c('#2F6BFF','#E74C3C','green','yellow','grey')

# names(color) <- names(group[!duplicated(group)])
# color <- color[group]
# dataset.name.order <- levels(df$source.target)
# dataset.name.order <- stringr::str_match(dataset.name.order,
#     "\\(.*\\)")
# dataset.name.order <- stringr::str_sub(dataset.name.order,
#     2, stringr::str_length(dataset.name.order) -
#       1)
# xtick.color <- color[dataset.name.order]
# g <- g + theme(axis.text.x = element_text(colour = xtick.color))

# 保存
pdf(paste0(outdir, "bubble_plot_Macrophages.pdf"), width = 5, height = 10)
print(g)
dev.off()
png(paste0(outdir, "bubble_plot_Macrophages.png"), width = 500, height = 1000)
print(g)
dev.off()

df<-df[df$interaction_name %in% c("TNFSF9_TNFRSF9","LGALS9_CD45","CCL7_CCR2","CCL6_CCR2","CCL2_CCR2","CCL12_CCR2"),]
g <- ggplot(df, aes(x = group, y = interaction_name_2,
        color = prob, size = pval)) + geom_point(pch = 16) +
        theme_linedraw() + theme(panel.grid.major = element_blank()) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1,
            vjust = 1), axis.title.x = element_blank(),
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
dot.size.max = max(df$pval)
dot.size.min = min(df$pval)

g <- g + scale_radius(range = c(dot.size.min*3, dot.size.max*3),
    breaks = sort(unique(df$pval)), labels = names(values)[values %in%
        sort(unique(df$pval))], name = "p-value")

color.use <- RColorBrewer::brewer.pal(n = 10, name = "Spectral")
color.use <- rev(color.use)

g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),
    na.value = "white", limits = c(quantile(df$prob,
        0, na.rm = T), quantile(df$prob, 1, na.rm = T)),
    breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob,
        1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5,
     title = "Commun. Prob."))

g <- g + theme(text = element_text(size = 10), plot.title = element_text(size = 10)) +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8))

if (length(unique(df$source.target)) > 1) {
    g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) -
        0.5, 1), lwd = 0.1, colour = "grey90")
}
    if (length(unique(df$interaction_name_2)) > 1) {
        g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) -
            0.5, 1), lwd = 0.1, colour = "grey90")
    }
# 保存
pdf(paste0(outdir, "bubble_plot_Macrophages_select.pdf"), width = 5, height = 5)
print(g)
dev.off()
png(paste0(outdir, "bubble_plot_Macrophages_select.png"), width = 500, height = 500)
print(g)
dev.off()


combined_data2<-combined_data[combined_data$source=="Fibroblasts",]
df<-combined_data2
df$group<-factor(df$group,levels=c("asham","aday1","aday3","aday7","aday14"))
g <- ggplot(df, aes(x = group, y = interaction_name_2,
        color = prob, size = pval)) + geom_point(pch = 16) +
        theme_linedraw() + theme(panel.grid.major = element_blank()) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1,
            vjust = 1), axis.title.x = element_blank(),
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
dot.size.max = max(df$pval)
dot.size.min = min(df$pval)

g <- g + scale_radius(range = c(dot.size.min*3, dot.size.max*3),
    breaks = sort(unique(df$pval)), labels = names(values)[values %in%
        sort(unique(df$pval))], name = "p-value")

color.use <- RColorBrewer::brewer.pal(n = 10, name = "Spectral")
color.use <- rev(color.use)

g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),
    na.value = "white", limits = c(quantile(df$prob,
        0, na.rm = T), quantile(df$prob, 1, na.rm = T)),
    breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob,
        1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5,
     title = "Commun. Prob."))

g <- g + theme(text = element_text(size = 10), plot.title = element_text(size = 10)) +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8))

if (length(unique(df$source.target)) > 1) {
    g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) -
        0.5, 1), lwd = 0.1, colour = "grey90")
}
    if (length(unique(df$interaction_name_2)) > 1) {
        g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) -
            0.5, 1), lwd = 0.1, colour = "grey90")
    }

# 保存
pdf(paste0(outdir, "bubble_plot_Fibroblasts.pdf"), width = 5, height = 10)
print(g)
dev.off()
png(paste0(outdir, "bubble_plot_Fibroblasts.png"), width = 500, height = 1000)
print(g)
dev.off()

df<-df[df$interaction_name %in% c("TNFSF9_TNFRSF9","LGALS9_CD45","CCL8_CCR2","CCL7_CCR2","CCL2_CCR2"),]
g <- ggplot(df, aes(x = group, y = interaction_name_2,
        color = prob, size = pval)) + geom_point(pch = 16) +
        theme_linedraw() + theme(panel.grid.major = element_blank()) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1,
            vjust = 1), axis.title.x = element_blank(),
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
dot.size.max = max(df$pval)
dot.size.min = min(df$pval)

g <- g + scale_radius(range = c(dot.size.min*3, dot.size.max*3),
    breaks = sort(unique(df$pval)), labels = names(values)[values %in%
        sort(unique(df$pval))], name = "p-value")

color.use <- RColorBrewer::brewer.pal(n = 10, name = "Spectral")
color.use <- rev(color.use)

g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),
    na.value = "white", limits = c(quantile(df$prob,
        0, na.rm = T), quantile(df$prob, 1, na.rm = T)),
    breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob,
        1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5,
     title = "Commun. Prob."))

g <- g + theme(text = element_text(size = 10), plot.title = element_text(size = 10)) +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8))

if (length(unique(df$source.target)) > 1) {
    g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) -
        0.5, 1), lwd = 0.1, colour = "grey90")
}
    if (length(unique(df$interaction_name_2)) > 1) {
        g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) -
            0.5, 1), lwd = 0.1, colour = "grey90")
    }
# 保存
pdf(paste0(outdir, "bubble_plot_Fibroblasts_select.pdf"), width = 5, height = 5)
print(g)
dev.off()
png(paste0(outdir, "bubble_plot_Fibroblasts_select.png"), width = 500, height = 500)
print(g)
dev.off()


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


pathway_outdir<-"/mnt/NexData_file/home/sh/Projects/SC2510235895/09.Mission/cellchat_all/"
select_pathway<-c('CCL','GALECTIN','CD137')
for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(cellchats, signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[9,]
	pairLR_top <- pairLR_top[10,]
	vertex_receiver <- vertex_cal(cellchats,k)
	vertex_receiver <- c(6,4,12)
	# pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-centrality_heatmap.pdf"),width = 8,height= 3)
	# gg <- netAnalysis_signalingRole_network(object.list[[1]], signaling = k, width = 14, height = 2.5, font.size = 10)
	# dev.off()

	# pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-contribution_barplot.pdf"),width = 8,height= 3)
	# print(netAnalysis_contribution(object.list[[1]], signaling = k))
	# dev.off()

	# pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
	# netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	# dev.off()
	pdf(paste0(pathway_outdir,k,"-all-",pairLR_top[9,],"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(cellchats, signaling = k, pairLR.use = pairLR_top[9,], vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
	pdf(paste0(pathway_outdir,k,"-all-",pairLR_top[10,],"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(cellchats, signaling = k, pairLR.use = pairLR_top[10,], vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
# 	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-",pairLR_top[4,],"-top_hierarchy.pdf"),width = 12,height= 8)
# 	netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top[4,], vertex.receiver = vertex_receiver, layout = "hierarchy")
# 	dev.off()
}
for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(cellchats, signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[1,]
	vertex_receiver <- c(6,4,12)
	pdf(paste0(pathway_outdir,k,"-all-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(cellchats, signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}
for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(cellchats, signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[1,]
	vertex_receiver <- c(6,4,12)
	pdf(paste0(pathway_outdir,k,"-all-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(cellchats, signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}

	pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
	for (i in seq_along(pdf_files)) {
		file_name <- pdf_files[i]
		pdf_image <- pdf_render_page(file_name, dpi = 300)
		writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
		}