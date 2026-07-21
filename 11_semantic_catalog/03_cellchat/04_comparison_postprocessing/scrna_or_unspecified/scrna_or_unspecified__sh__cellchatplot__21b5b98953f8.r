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

load("/mnt/NexData_file/home/sh/Projects/AE2511106372-01/05.Mission/MI_DZ_MI-cell_type/tmp_object.RData")
#### 单独做细胞类型贝壳图 ####
for(j in 1:2){
	mat_weight = cellchat@net[[j]]$weight
	pdf(paste0(outdir,"",names(cellchat@net)[j],"_Weight_netVisual_circle.pdf"),width = 7,height = 7)
	# par(mfrow=n,mar=c(0.3,0.3,0.3,0.3),xpd=TRUE)
	# for (i in 1:nrow(mat_count)) {
	i=2
		mat_weight2 = matrix(0,nrow = nrow(mat_weight),ncol = ncol(mat_weight),dimnames = dimnames(mat_weight))
		mat_weight2[i,] = mat_weight[i,]
		netVisual_circle(mat_weight2,edge.width.max = 20,vertex.weight = 11,weight.scale = T,label.edge = T,edge.weight.max = max(mat_weight),edge.label.cex = 1,title.name = rownames(mat_weight)[i])
	# }
	dev.off()
	}


cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/sh/Projects/AE2511106372-01/06.Mission/"
# pdf(paste0(outdir,"Endothelial cells_weight_netVisual_circle.pdf"),width = 12,height = 12)
# netVisual_diffInteraction(cellchat,measure = "weight", comparison = c(2, 1),weight.scale = T,label.edge = T,edge.label.cex = 1,title.name = 'Endothelial cells',sources.use ='Endothelial cells')
# dev.off()
# pdf(paste0(outdir,"Macrophages_weight_netVisual_circle.pdf"),width = 12,height = 12)
# netVisual_diffInteraction(cellchat,measure = "weight", comparison = c(2, 1),weight.scale = T,label.edge = T,edge.label.cex = 1,title.name = 'Macrophages',sources.use ='Macrophages')
# dev.off()

pathway_outdir<-paste0(outdir,"pathway/")
dir.create(pathway_outdir, recursive = T)
select_pathway<-c('ADGRE','ADGRG','IGF','MIF','OSM','SEMA4','TNF','COLLAGEN','GALECTIN','CXCL','ApoE','LAMININ','JAM','ICAM','PECAM2','NOTCH')
for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(object.list[[1]], signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[1,]
	vertex_receiver <- vertex_cal(object.list[[1]],k)
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-centrality_heatmap.pdf"),width = 8,height= 3)
	gg <- netAnalysis_signalingRole_network(object.list[[1]], signaling = k, width = 14, height = 2.5, font.size = 10)
	dev.off()

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-contribution_barplot.pdf"),width = 8,height= 3)
	print(netAnalysis_contribution(object.list[[1]], signaling = k))
	dev.off()

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[1],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
	netVisual_individual(object.list[[1]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
	dev.off()
}

for(k in str_c(select_pathway)){
	pairLR_top <- extractEnrichedLR(object.list[[2]], signaling = k, geneLR.return = FALSE)
	pairLR_top <- pairLR_top[1,]
	vertex_receiver <- vertex_cal(object.list[[2]],k)
	pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-centrality_heatmap.pdf"),width = 8,height= 3)
	gg <- netAnalysis_signalingRole_network(object.list[[2]], signaling = k, width = 14, height = 2.5, font.size = 10)
	dev.off()

	pdf(paste0(pathway_outdir,k,"-",names(object.list)[2],"-contribution_barplot.pdf"),width = 8,height= 3)
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
pdf_files <- list.files(path = pathway_outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
    file_name <- pdf_files[i]
    pdf_image <- pdf_render_page(file_name, dpi = 300)
    writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}

pairLRuse<-data.frame(interaction_name=c("THBS1_CD47","THBS1_CD36","TGM2_ADGRG1",'SEMA4A_NRP1_PLXNA2','OSM_OSMR_IL6ST','MIF_ACKR3','ITGAV_ITGB1_ADGRE5'))

cellchat_sub<-subsetCommunication(cellchat)
pdf(paste0(outdir,"Interactions_diff_netplot.pdf"),width = 12,height = 7)
par(mfrow = c(1,2), xpd=TRUE)
for(i in c("count","weight")){
	netVisual_diffInteraction(cellchat_sub, comparison = c(2, 1), weight.scale = T,measure = i)
	}
dev.off()

signalinguse<-c("ADGRG","APP","ApoE","COLLAGEN","CXCL","GALECTIN","GAS","ICAM","IGFBP","JAM","LAMININ","MIF","NOTCH","PECAM2","PROS","SEMA6","SEMA7","TGFb")
signalinguse<-c("THBS","ADGRG","SEMA4","OSM","MIF","ADGRE")
signalinguse<-c('MHC-II','LCK','GP1BA','IL16','GDF','NPR1','CD96','NT','MHC-I','ANNEXIN','ADGRL','CD86','SLIT','FGF','SLITRK','NOTCH')
# pairLR <- searchPair(signaling = signalinguse, pairLR.use = cellchat@LR$LRsig,
#         key = "pathway_name", matching.exact = T, pair.only = T)
pdf(paste0(outdir,"Interactions_diff_netplot.pdf"),width = 12,height = 7)
# 2. 获取合并后的最大权重，用于统一两幅图的颜色/线条比例
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
# 3. 设置并排绘图
par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(object.list)) {
		mat_weight = cellchat@net[[i]]$weight
    netVisual_aggregate(
        object.list[[i]],idents.use = c("Endothelial cells","Macrophages"), 
        signaling = signalinguse,  # 只展示你感兴趣的配受体对聚合后的网络
        layout = "circle",
        weight.scale = TRUE,
        # edge.weight.max = weight.max[2],   # 统一两幅图的边权重范围，使对比更公平
        edge.width.max = 20,
        # title.name = NULL
    )
}
dev.off()

pdf(paste0(outdir,"netVisual_bubble.pdf"))
netVisual_bubble(cellchat, sources.use = "Macrophages", targets.use = c('Endothelial cells','Fibroblasts'), 
                thresh = 1, max.dataset = 2,comparison = c(2, 1),pairLR.use = pairLRuse,
								dot.size.max = 10,dot.size.min = 5,font.size = 15,
                remove.isolate = F,return.data = FALSE)
dev.off()
png(paste0(outdir,"netVisual_bubble.png"))
netVisual_bubble(cellchat, sources.use = "Macrophages", targets.use = c('Endothelial cells','Fibroblasts'), 
                thresh = 1, max.dataset = 2,comparison = c(2, 1),pairLR.use = pairLRuse,
								dot.size.max = 10,dot.size.min = 5,font.size = 15,
                remove.isolate = F,return.data = FALSE)
dev.off()



#### 后面没用 ####
select<-c("Macrophages")
prefix<-select
pairLR.use<-as.data.frame(c('MHC-II','LCK','GP1BA','IL16','GDF','NPR1','CD96','NT','MHC-I','ANNEXIN','ADGRL','CD86','SLIT','FGF','SLITRK','NOTCH'))
colnames(pairLR.use) <- 'interaction_name'
source("/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/script/cellchat_netVisual_bubble.r")
data <- netVisual_bubble_mod(cellchat, thresh = 1,targets.use = 'Endothelial cells',
                 max.dataset = 2,comparison = c(2, 1),
                remove.isolate = F,return.data = TRUE)
# sources.use = select,
pathway<-c('MHC-II','LCK','GP1BA','IL16','GDF','NPR1','CD96','NT','MHC-I','ANNEXIN','ADGRL','CD86','SLIT','FGF','SLITRK','NOTCH')
df<-data$communication[which(data$communication$pathway_name %in% pathway),]
df_control<-df
df_control$prob<-0
df_control$pval<-1
df_control$source.target<-gsub("case","control",df_control$source.target)
df_control$dataset<-"control"
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
}
else {
    g <- g + scale_colour_gradientn(colors = colorRampPalette(color.use)(99),
        na.value = "white") + guides(color = guide_colourbar(barwidth = 0.5,
        title = "Commun. Prob."))
}
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

comparison = c(2, 1)
dataset.name <- names(cellchat@net)
group.names0 <- df$source.target
line.on<-TRUE
color.text.use<-TRUE
color.text<-NULL
    if (!is.null(comparison)) {
        if (line.on) {
            xintercept = seq(0.5 + length(dataset.name[comparison]),
                length(group.names0) * length(dataset.name[comparison]),
                by = length(dataset.name[comparison]))
            g <- g + geom_vline(xintercept = xintercept, linetype = "dashed",
                color = "grey60", size = 0.2)
        }
        if (color.text.use) {
            if (is.null(group)) {
                group <- 1:length(comparison)
                names(group) <- dataset.name[comparison]
            }
            if (is.null(color.text)) {
                # color <- ggPalette(length(unique(group)))
								color <- c('#2F6BFF','#E74C3C')
            }
            else {
                color <- color.text
            }
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
        }
    }

pdf(paste0(outdir,'/',prefix,"_netVisual_bubble.pdf"),width = 18,height= 6)#2：1，下调
print(g)
dev.off()
ggsave(paste0(outdir,'/',prefix,"_netVisual_bubble.pdf"), width= 12, height= 6, units='in')
ggsave(paste0(outdir,'/',prefix,"_netVisual_bubble.png"), width= 12, height= 6, units='in')

pdf(paste0(outdir,'/',prefix,"_netVisual_bubble.pdf"),width = 12,height= 15)#2：1，下调
print(netVisual_bubble_mod(cellchat, thresh = 0.05,
                sources.use = select, 
								# targets.use = select, 
                max.dataset = 2,
                remove.isolate = F,
                comparison = c(2, 1), 
								angle.x = 45))
dev.off()