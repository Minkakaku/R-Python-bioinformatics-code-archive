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

load("/mnt/NexData_file/home/sh/Projects/AE2512027077-01/04.Mission/OVX_CNS-defect_sham_CNS-defect/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/sh/Projects/AE2512027077-01/06.Mission/OVX_CNS-defect_sham_CNS-defect_plot/"

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

pathway_outdir<-outdir
select_pathway<-c('APP')
for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(object.list[[1]], signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[1,]
	vertex_receiver <- vertex_cal(object.list[[1]],k)

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-centrality_heatmap.pdf"),width = 8,height= 3)
	gg <- netAnalysis_signalingRole_network(object.list[[1]], signaling = k, width = 16, height = 2.5, font.size = 10)
	dev.off()

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-contribution_barplot.pdf"),width = 8,height= 8)
	print(netAnalysis_contribution(object.list[[1]], signaling = k))
	dev.off()
	vertex_receiver <- c(11,12) #### 指定细胞类型
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}

for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(object.list[[2]], signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[1,]
	# vertex_receiver <- vertex_cal(object.list[[2]],k)
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-centrality_heatmap.pdf"),width = 8,height= 3)
	gg <- netAnalysis_signalingRole_network(object.list[[2]], signaling = k, width = 16, height = 2.5, font.size = 10)
	dev.off()

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-contribution_barplot.pdf"),width = 8,height= 8)
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

select<-c("Microglia")
prefix<-select
source("/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/script/cellchat_netVisual_bubble.r")
data <- netVisual_bubble_mod(cellchat, thresh = 1,
                sources.use = select,targets.use ='NaiveB',max.dataset = 2,comparison = c(2, 1), 
                remove.isolate = T,return.data = TRUE)

df<-data$communication
df_control<-df
df_control$prob<-0
df_control$pval<-1
df_control$source.target<-gsub("sham_CNS-defect","OVX_CNS-defect",df_control$source.target)
df_control$dataset<-"OVX_CNS-defect"
df<-rbind(df,df_control)
g <- ggplot(df, aes(x = source.target, y = interaction_name_2,
        color = prob, size = pval)) + geom_point(pch = 16) +
        theme_linedraw() + theme(panel.grid.major = element_blank()) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1,
            vjust = 1), axis.title.x = element_blank(),
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
dot.size.max = max(df$pval)
dot.size.min = min(df$pval)
dot.size.max = 8
g <- g + scale_radius(range = c(dot.size.min*3, dot.size.max*3),
    breaks = sort(unique(df$pval)), labels = names(values)[values %in%
        sort(unique(df$pval))], name = "p-value")

color.use <- RColorBrewer::brewer.pal(n = 10, name = "Spectral")
color.use <- rev(color.use)
if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),
        na.value = "white", limits = c(quantile(df$prob,
            0, na.rm = T), quantile(df$prob, 1, na.rm = T)),
        breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob,
            1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5,
        title = "Commun. Prob."))
}else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),
        na.value = "white") + guides(color = guide_colourbar(barwidth = 0.5,
        title = "Commun. Prob."))
}
g <- g + theme(text = element_text(size = 20), plot.title = element_text(size = 20)) +
    theme(legend.title = element_text(size = 15), legend.text = element_text(size = 12))

if (length(unique(df$source.target)) > 1) {
    g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) -
        0.5, 1), lwd = 0.1, colour = "grey90")
}
if (length(unique(df$interaction_name_2)) > 1) {
     g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) -
        0.5, 1), lwd = 0.1, colour = "grey90")
}

comparison = c(2, 1)
dataset.name <- names(cellchat@net)
group.names0 <- df$source.target
line.on<-TRUE
color.text.use<-TRUE
color.text<-NULL
# group<-NULL
    if (!is.null(comparison)) {
        if (line.on) {
            # xintercept = seq(0.5 + length(dataset.name[comparison]),
            #     length(group.names0) * length(dataset.name[comparison]),
            #     by = length(dataset.name[comparison]))
            # g <- g + geom_vline(xintercept = xintercept, linetype = "dashed",
            #     color = "grey60", size = 0.2)
        }
        if (color.text.use) {
            # if (is.null(group)) {
                group <- 1:length(comparison)
                names(group) <- dataset.name[comparison]
            # }
            # if (is.null(color.text)) {
                # color <- ggPalette(length(unique(group)))
								color <- c('#2F6BFF','#E74C3C')
            # }
            # else {
                # color <- color.text
            # }
            names(color) <- names(group[!duplicated(group)])
            # color <- color[group]
            dataset.name.order <- unique(df$source.target)
            # dataset.name.order <- stringr::str_match(dataset.name.order,
            #     "\\(.*\\)")
            # dataset.name.order <- stringr::str_sub(dataset.name.order,
            #     2, stringr::str_length(dataset.name.order) -
            #       1)
            xtick.color <- color[dataset.name.order]
            g <- g + theme(axis.text.x = element_text(colour = xtick.color))
        }
    }

pdf(paste0(outdir,'/',prefix,"_netVisual_bubble.pdf"),width = 6,height= 8)#2：1，下调
print(g)
dev.off()
png(paste0(outdir,'/',prefix,"_netVisual_bubble.png"),width = 600,height= 800)#2：1，下调
print(g)
dev.off()



load("/mnt/NexData_file/home/sh/Projects/AE2512027077-01/01.Mission/OVX_bone-CNS_sham_bone-CNS/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/sh/Projects/AE2512027077-01/06.Mission/OVX_bone-CNS_sham_bone-CNS_plot/"

pdf(paste0(outdir,"Neutrophil_G4_5-Microglia_Number_netVisual_circle.pdf"),width = 12,height = 12)
netVisual_diffInteraction(cellchat,measure = "count", comparison = c(2, 1),weight.scale = T,label.edge = T,edge.label.cex = 1,title.name = 'Neutrophil_G4_5',targets.use =c('Neutrophil_G4','Neutrophil_G5'),sources.use ='Microglia',vertex.size = 5)

pathway_outdir<-outdir
select_pathway<-c('CCL')
for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(object.list[[1]], signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[8,]
	vertex_receiver <- vertex_cal(object.list[[1]],k)

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-centrality_heatmap.pdf"),width = 8,height= 3)
	gg <- netAnalysis_signalingRole_network(object.list[[1]], signaling = k, width = 16, height = 2.5, font.size = 10)
	dev.off()

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-contribution_barplot.pdf"),width = 8,height= 8)
	print(netAnalysis_contribution(object.list[[1]], signaling = k))
	dev.off()
	vertex_receiver <- c(19,20,13) #### 指定细胞类型
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}

for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(object.list[[2]], signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[8,]
	# vertex_receiver <- vertex_cal(object.list[[2]],k)
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-centrality_heatmap.pdf"),width = 8,height= 3)
	gg <- netAnalysis_signalingRole_network(object.list[[2]], signaling = k, width = 16, height = 2.5, font.size = 10)
	dev.off()

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-contribution_barplot.pdf"),width = 8,height= 8)
	print(netAnalysis_contribution(object.list[[2]], signaling = k))
	dev.off()

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(object.list[[2]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}

for(i in 1:7){
for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(object.list[[1]], signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[i,]
	vertex_receiver <- c(19,20,13) #### 指定细胞类型
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}

for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(object.list[[2]], signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[i,]
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(object.list[[2]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}
}

pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
    file_name <- pdf_files[i]
    pdf_image <- pdf_render_page(file_name, dpi = 300)
    writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}

select<-c("Neutrophil_G4","Neutrophil_G5")
prefix<-select
# pdf(paste0(outdir,'/',prefix,"_netVisual_bubble_all.pdf"),width = 5,height= 5)#2：1，下调
# print(netVisual_bubble_mod(cellchat, thresh = 0.05,
#                 sources.use = select,targets.use ='Microglia', 
#                 max.dataset = 2,
#                 remove.isolate = T,
#                 comparison = c(2, 1), 
# 								angle.x = 45))
# dev.off()

data <- netVisual_bubble_mod(cellchat, thresh = 1,
                sources.use = select,targets.use ='Microglia',max.dataset = 2,comparison = c(2, 1), 
                remove.isolate = T,return.data = TRUE)
df<-data$communication
df_control<-df[c(11,12,13,14),]
df_control$prob<-0
df_control$pval<-1
df_control$source.target<-gsub("sham_bone-CNS","OVX_bone-CNS",df_control$source.target)
df_control$dataset<-"OVX_bone-CNS"
df_control1<-df[c(5:10),]
df_control1$prob<-0
df_control1$pval<-1
df_control1$source.target<-gsub("OVX_bone-CNS","sham_bone-CNS",df_control1$source.target)
df_control1$dataset<-"sham_bone-CNS"
df<-rbind(df,df_control,df_control1)
# df$source.target<-df$dataset
# df<-df[which(df$target %in% c('Astrocytes','Ependymal cells','Microglia','Neuronendocrine cells','Oligodendrocytes','OPC')),]
# df<-df[-c(12:13),]
g <- ggplot(df, aes(x = source.target, y = interaction_name_2,
        color = prob, size = pval)) + geom_point(pch = 16) +
        theme_linedraw() + theme(panel.grid.major = element_blank()) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1,
            vjust = 1), axis.title.x = element_blank(),
            axis.title.y = element_blank()) + scale_x_discrete(position = "bottom")
values <- c(1, 2, 3)
names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
dot.size.max = max(df$pval)
dot.size.min = min(df$pval)
dot.size.max = 8
g <- g + scale_radius(range = c(dot.size.min*3, dot.size.max*3),
    breaks = sort(unique(df$pval)), labels = names(values)[values %in%
        sort(unique(df$pval))], name = "p-value")

color.use <- RColorBrewer::brewer.pal(n = 10, name = "Spectral")
color.use <- rev(color.use)
if (min(df$prob, na.rm = T) != max(df$prob, na.rm = T)) {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),
        na.value = "white", limits = c(quantile(df$prob,
            0, na.rm = T), quantile(df$prob, 1, na.rm = T)),
        breaks = c(quantile(df$prob, 0, na.rm = T), quantile(df$prob,
            1, na.rm = T)), labels = c("min", "max")) + guides(color = guide_colourbar(barwidth = 0.5,
        title = "Commun. Prob."))
}else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),
        na.value = "white") + guides(color = guide_colourbar(barwidth = 0.5,
        title = "Commun. Prob."))
}
g <- g + theme(text = element_text(size = 20), plot.title = element_text(size = 20)) +
    theme(legend.title = element_text(size = 15), legend.text = element_text(size = 12))

if (length(unique(df$source.target)) > 1) {
    g <- g + geom_vline(xintercept = seq(1.5, length(unique(df$source.target)) -
        0.5, 1), lwd = 0.1, colour = "grey90")
}
if (length(unique(df$interaction_name_2)) > 1) {
     g <- g + geom_hline(yintercept = seq(1.5, length(unique(df$interaction_name_2)) -
        0.5, 1), lwd = 0.1, colour = "grey90")
}

comparison = c(2, 1)
dataset.name <- names(cellchat@net)
group.names0 <- df$source.target
line.on<-TRUE
color.text.use<-TRUE
color.text<-NULL
# group<-NULL
    if (!is.null(comparison)) {
        if (line.on) {
            xintercept = seq(0.5 + length(dataset.name[comparison]),
                length(group.names0) * length(dataset.name[comparison]),
                by = length(dataset.name[comparison]))
            g <- g + geom_vline(xintercept = xintercept, linetype = "dashed",
                color = "grey60", size = 0.2)
        }
        if (color.text.use) {
            # if (is.null(group)) {
                group <- 1:length(comparison)
                names(group) <- dataset.name[comparison]
            # }
            # if (is.null(color.text)) {
                # color <- ggPalette(length(unique(group)))
								color <- c('#2F6BFF','#E74C3C')
            # }
            # else {
                # color <- color.text
            # }
            names(color) <- names(group[!duplicated(group)])
            # color <- color[group]
            dataset.name.order <- levels(df$source.target)
            dataset.name.order <- stringr::str_match(dataset.name.order,
                "\\(.*\\)")
            dataset.name.order <- stringr::str_sub(dataset.name.order,
                2, stringr::str_length(dataset.name.order) -
                  1)
            xtick.color <- color[dataset.name.order]
            g <- g + theme(axis.text.x = element_text(colour = xtick.color))
        }
    }

pdf(paste0(outdir,"Neutrophil_G4_5_netVisual_bubble.pdf"),width = 6,height= 8)#2：1，下调
print(g)
dev.off()
png(paste0(outdir,"Neutrophil_G4_5_netVisual_bubble.png"),width = 600,height= 800)#2：1，下调
print(g)
dev.off()