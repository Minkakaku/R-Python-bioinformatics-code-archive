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

load("/mnt/NexData_file/home/sh/Projects/FWSC20241287/07.Mission/newcellchat/MODEL_CONTROL/tmp_object.RData")
load("/mnt/NexData_file/home/sh/Projects/FWSC20241287/07.Mission/newcellchat/PF11-H_MODEL/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC20241287/07.Mission/cellchat_plot/"

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
select_pathway<-c('LAMININ')
for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(object.list[[1]], signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[26,]
	vertex_receiver <- vertex_cal(object.list[[1]],k)
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-centrality_heatmap.pdf"),width = 8,height= 3)
	gg <- netAnalysis_signalingRole_network(object.list[[1]], signaling = k, width = 14, height = 2.5, font.size = 10)
	dev.off()

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-contribution_barplot.pdf"),width = 8,height= 10)
	print(netAnalysis_contribution(object.list[[1]], signaling = k))
	dev.off()

	vertex_receiver<-c(5,10,15)
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-",pairLR_top,"-top_hierarchy.pdf"),width = 15,height= 8)
	netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}

for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(object.list[[2]], signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[26,]
	vertex_receiver <- vertex_cal(object.list[[2]],k)
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-centrality_heatmap.pdf"),width = 8,height= 3)
	gg <- netAnalysis_signalingRole_network(object.list[[2]], signaling = k, width = 14, height = 2.5, font.size = 10)
	dev.off()

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-contribution_barplot.pdf"),width = 8,height= 10)
	print(netAnalysis_contribution(object.list[[2]], signaling = k))
	dev.off()

	vertex_receiver<-c(5,10,15)
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-",pairLR_top,"-top_hierarchy.pdf"),width = 15,height= 8)
	netVisual_individual(object.list[[2]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}

pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
    file_name <- pdf_files[i]
    pdf_image <- pdf_render_page(file_name, dpi = 300)
    writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}

cellchat <- mergeCellChat(list(object.list1$CONTROL,object.list1$MODEL,object.list[['PF11-H']]),add.names = c('CONTROL','MODEL','PF11-H'))
measure_list <- c("count","weight")
heatmap_color <- c("white","red2")
	for (i in seq_along(measure_list)) {
		gg1 = netVisual_heatmap(object.list1$CONTROL, measure = measure_list[i],color.heatmap = heatmap_color,title.name = 'CONTROL')
		gg2 = netVisual_heatmap(object.list1$MODEL, measure = measure_list[i],color.heatmap = heatmap_color,title.name = 'MODEL')
		gg3 = netVisual_heatmap(object.list[['PF11-H']], measure = measure_list[i],color.heatmap = heatmap_color,title.name = 'PF11-H')
		gg <- gg1+gg2+gg3
		if(i==1){label <- "Number"}else{label <- "Weight"}
		pdf(paste0(outdir,"group_",label,"_netVisual_heatmap.pdf"),width = 18,height = 7)
		print(gg)
		dev.off()
		}

library(ComplexHeatmap)
library(circlize)
# 提取group注释信息
group_anno <- gg1@left_annotation@anno_list$group
group_colors <- setNames(group_anno@color_mapping@colors, group_anno@color_mapping@levels)

# 创建统一的注释（注意：应该是列注释放在top/bottom，行注释放在left/right）
column_ha <- HeatmapAnnotation(
    cell_type = group_anno@color_mapping@levels,
    col = list(cell_type = group_colors),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10),
		show_annotation_name = FALSE,
    show_legend = FALSE
)

row_ha <- rowAnnotation(
    cell_type = group_anno@color_mapping@levels,
    col = list(cell_type = group_colors),
    annotation_name_side = "bottom",
    annotation_name_gp = gpar(fontsize = 10),
		show_annotation_name = FALSE,
    show_legend = TRUE
)
# 白红配色函数
col_fun <- colorRamp2(
    breaks = seq(0, 2, length.out = 100),
    colors = colorRampPalette(c("#BDA0FF", "#DDA0DD",'white',
                  "#E6B0D0", "#FAB9B6", "#FFB6C1",
                  "#FF9FAA", "#FF8994", "#FF7D89",
                  "#FF6773", "#FF515D", "#FF464F",
                  "#FF3035", "#FF1A1A", "#E61717"))(100)
)
# 创建热图列表
create_heatmap <- function(mat, top_anno, left_anno,title,label) {
    Heatmap(
        mat,
        name = "Interaction",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = col_fun,
        bottom_annotation = top_anno,
        left_annotation = left_anno,
        show_row_names = TRUE,row_names_side = c("left"),
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        column_names_rot = 45,
				column_title = title, column_title_side = c("top"),
        heatmap_legend_param = list(
            title = label,
            direction = "vertical",
            title_position = "topcenter"
        )
    )
}
gg1@matrix[is.na(gg1@matrix)] <- 0
gg2@matrix[is.na(gg2@matrix)] <- 0
gg3@matrix[is.na(gg3@matrix)] <- 0
gg1 <- create_heatmap(gg1@matrix, column_ha, row_ha, 'CONTROL', 'Weight of interaction')
gg2 <- create_heatmap(gg2@matrix, column_ha, row_ha, 'MODEL', 'Weight of interaction')
gg3 <- create_heatmap(gg3@matrix, column_ha, row_ha, 'PF11-H', 'Weight of interaction')
gg <- gg1 + gg2 + gg3
pdf(paste0(outdir, "group_Weight_netVisual_heatmap.pdf"), width = 22, height = 8)
draw(gg, 
     heatmap_legend_side = "right",      # 热图图例放在右侧
     annotation_legend_side = "right",   # 注释图例也放在右侧
     merge_legend = TRUE,                 # 合并图例
     legend_grouping = "original")        # 保持图例分组
dev.off()

col_fun <- colorRamp2(
    breaks = seq(0, 200, length.out = 100),
    colors = colorRampPalette(c("#BDA0FF", "#DDA0DD",'white',
                  "#E6B0D0", "#FAB9B6", "#FFB6C1",
                  "#FF9FAA", "#FF8994", "#FF7D89",
                  "#FF6773", "#FF515D", "#FF464F",
                  "#FF3035", "#FF1A1A", "#E61717"))(100)
)
gg1@matrix[is.na(gg1@matrix)] <- 0
gg2@matrix[is.na(gg2@matrix)] <- 0
gg3@matrix[is.na(gg3@matrix)] <- 0
gg1 <- create_heatmap(gg1@matrix, column_ha, row_ha, 'CONTROL', 'Number of interaction')
gg2 <- create_heatmap(gg2@matrix, column_ha, row_ha, 'MODEL', 'Number of interaction')
gg3 <- create_heatmap(gg3@matrix, column_ha, row_ha, 'PF11-H', 'Number of interaction')
gg <- gg1 + gg2 + gg3
pdf(paste0(outdir, "group_Number_netVisual_heatmap.pdf"), width = 22, height = 8)
draw(gg, 
     heatmap_legend_side = "right",      # 热图图例放在右侧
     annotation_legend_side = "right",   # 注释图例也放在右侧
     merge_legend = TRUE,                 # 合并图例
     legend_grouping = "original")        # 保持图例分组
dev.off()

pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
    file_name <- pdf_files[i]
    pdf_image <- pdf_render_page(file_name, dpi = 300)
    writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}


#### 画气泡图
data <- netVisual_bubble(
    cellchat,
    # sources.use = 'Remodeling Mac-Svep1',
    # targets.use = "Mesangial cells",
    comparison = c(1, 2, 3),  # 比较三组
    title.name = "CellChat Comparison",
    font.size = 8,
    angle.x = 45,
    return.data = TRUE
)
df<-data$communication[which(data$communication$source == 'Remodeling Mac-Svep1'),]
df<-df[which(df$target == 'Mesangial cells'),]
source("/mnt/NexData_file/home/sh/Projects/FWSC20250829/09.Mission/script/cellchat_netVisual_bubble.r")
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

comparison = c(1,2,3)
dataset.name <- names(cellchat@net)
group.names0 <- df$source.target

xintercept = seq(0.5 + length(dataset.name[comparison]),
    length(group.names0) * length(dataset.name[comparison]),
    by = length(dataset.name[comparison]))
g <- g + geom_vline(xintercept = xintercept, linetype = "dashed",
    color = "grey60", size = 0.2)

group <- 1:length(comparison)
names(group) <- dataset.name[comparison]
color <- c('#2F6BFF','#E74C3C','green')

names(color) <- names(group[!duplicated(group)])
color <- color[group]
dataset.name.order <- levels(df$source.target)
dataset.name.order <- stringr::str_match(dataset.name.order,
    "\\(.*\\)")
dataset.name.order <- stringr::str_sub(dataset.name.order,
    2, stringr::str_length(dataset.name.order) -
      1)
xtick.color <- color[dataset.name.order]
g <- g + theme(axis.text.x = element_text(colour = xtick.color))

# 保存
pdf(paste0(outdir, "bubble_plot_comparison.pdf"), width = 5, height = 36)
print(g)
dev.off()
pdf(paste0(outdir, "bubble_plot_comparison.pdf"), width = 5, height = 36)
print(g)
dev.off()