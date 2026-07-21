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
rm(list = ls())
setwd('/mnt/NexData_file/home/zwl/Project/FWSC20241362/sample12/cellchat/OSCC_OSF/')
outdir<-"/mnt/NexData_file/home/zwl/Project/FWSC20241362/sample12/cellchat/"
######################################合并三个数据######################################
load("/mnt/NexData_file/home/zwl/Project/FWSC20241362/sample12/cellchat/OSF_NL/tmp_object.RData")
object.list_NL = object.list
load("/mnt/NexData_file/home/zwl/Project/FWSC20241362/sample12/cellchat/OSCC_OSF/tmp_object.RData")
object.list_can = object.list
object.list_all = object.list_can;object.list_all[['NL']] = object.list_NL[['NL']]

idents.list = lapply(object.list_all,FUN = function(x){
  length(levels(x@idents))
})

cellchat_all <- mergeCellChat(object.list_all,add.names = names(object.list_all))
# 细胞类型向最多的看齐
cellchat_all <- liftCellChat(cellchat_all , group.new = levels( cellchat_all@idents$joint))


#############################提取数据绘制气泡图###################################################
# OSF和OSCC
load("/mnt/NexData_file/home/zwl/Project/FWSC20241362/sample12/cellchat/OSCC_OSF/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))

rownames(cellchat@net$OSF$weight)
#netVisual_bubble(cellchat, sources.use = c(1:9), targets.use = c(1:9))

object.list %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        source %in% c('myCAFs','Epi_SLIT2','BTNL9+ cECs'),
        target  %in% c('CD8 Tex'),
        # pathway_name == "TGFb",
        pval < 0.05
      )  %>%
      mutate(cell_pair = str_c(source, target, sep = " -> \n")) 
  }) %>%
  enframe('group') %>%
  unnest(value) %>%
  
  ggplot(aes(x = group , y = interaction_name_2)) +
  geom_point(aes(color = prob, size = pval)) +
  facet_wrap(~ cell_pair, nrow = 1) +
  theme_classic() +
  theme(strip.background = element_rect(colour = NA)) +
  scale_color_gradientn(colours  = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu"))) +
  scale_size(transform = "reverse", range = c(1, 3)) +
  Seurat::RotatedAxis()


##################################正常绘制#############################################
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

pathway_outdir<-paste0(outdir,"pathway/")
dir.create(pathway_outdir, recursive = T)
select_pathway<-c('MHC-I','MIF')
for (si in length(object.list_all)) {
  for(k in str_c(select_pathway)){
    pairLR_top <- extractEnrichedLR(object.list_all[[si]], signaling = k, geneLR.return = FALSE)
    pairLR_top <- pairLR_top[1,]
    vertex_receiver <- vertex_cal(object.list_all[[si]],k)
    pdf(paste0(pathway_outdir,k,"-",names(object.list_all)[si],"-centrality_heatmap.pdf"),width = 10,height= 3)
    gg <- netAnalysis_signalingRole_network(object.list_all[[si]], signaling = k, width = 17, height = 2.5, font.size = 10)
    dev.off()
    
    pdf(paste0(pathway_outdir,k,"-",names(object.list_all)[si],"-contribution_barplot.pdf"),width = 10,height= 3)
    print(netAnalysis_contribution(object.list_all[[si]], signaling = k))
    dev.off()
    
    pdf(paste0(pathway_outdir,k,"-",names(object.list_all)[si],"-",pairLR_top,"-top_hierarchy.pdf"),width = 12,height= 8)
    netVisual_individual(object.list_all[[si]], signaling = k, pairLR.use = pairLR_top, vertex.receiver = vertex_receiver, layout = "hierarchy")
    dev.off()
  }
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



select<-c("Non-classical Monocytes")
prefix<-select
pairLR.use<-as.data.frame(c('IL16_CD4','TNFSF10_TNFRSF10A','TNFSF10_TNFRSF10B'))
colnames(pairLR.use) <- 'interaction_name'
data <- netVisual_bubble_mod(cellchat, thresh = 1,
                             sources.use = select, max.dataset = 2,comparison = c(2, 1), 
                             remove.isolate = F,return.data = TRUE)
source("/mnt/NexData_file/home/sh/Projects/FWSC20250829/07.Mission/script/cellchat_netVisual_bubble.r")


pathway<-c("IL16","TRAIL")
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