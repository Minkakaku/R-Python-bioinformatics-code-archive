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

load("/mnt/NexData_file/home/sh/Projects/FWSC20250829/04.Mission/case_control/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/sh/Projects/FWSC20250829/09.Mission/cellchat_plot/"

select<-c("CD4 Treg")
prefix<-select
pairLR.use<-as.data.frame(c('LGALS9_CD45','LGALS9_HAVCR2','LGALS9_P4HB'))
colnames(pairLR.use) <- 'interaction_name'
source("/mnt/NexData_file/home/sh/Projects/FWSC20250829/09.Mission/script/cellchat_netVisual_bubble.r")

data <- netVisual_bubble_mod(cellchat, thresh = 1,
                sources.use = select, max.dataset = 2,comparison = c(2, 1), 
                remove.isolate = F,return.data = TRUE)

pathway<-c("GALECTIN")
df<-data$communication[which(data$communication$pathway_name %in% pathway),]

targets_select<-c("Immature NK","Mature NK c1","Mature NK c2","NKT")
targets_select<-c('Classical Monocytes CD302', 'Classical Monocytes IFI44', 'Classical Monocytes NFKBIA', 'Classical Monocytes RETN ','Macrophage', 'Non-classical Monocytes', 'cDCs', 'pDCs')
targets_select<-c("CD4 Cytotoxic Tcells","CD4 Tem","CD4 Treg","CD4 naive Tcells","CD8 Cytotoxic Tcells","CD8 Tcm","CD8 naive Tcells","Gamma-delta Tcells","MAIT")

pdf(paste0(outdir,'/NK_netVisual_bubble.pdf'),width = 12,height= 15)#2：1，下调
print(netVisual_bubble(cellchat, thresh = 0.05,
                sources.use = select, 
                targets.use = targets_select, 
                pairLR.use = pairLR.use, 
                max.dataset = 2,
                remove.isolate = T,
                comparison = c(2, 1), 
                angle.x = 60))
dev.off()


pdf_files <- list.files(path = outdir, pattern = "*.pdf", full.names = TRUE)
for (i in seq_along(pdf_files)) {
    file_name <- pdf_files[i]
    pdf_image <- pdf_render_page(file_name, dpi = 300)
    writeJPEG(pdf_image, str_replace(string = file_name,pattern = "pdf",replacement = "png"))
}

pairLR.use<-as.data.frame(c('LGALS9_CD45','LGALS9_HAVCR2','LGALS9_P4HB'))
colnames(pairLR.use) <- 'interaction_name'

data <- netVisual_bubble_mod(cellchat, thresh = 1,
                sources.use = select,targets.use = targets_select, max.dataset = 2,comparison = c(2, 1), 
                remove.isolate = F,return.data = TRUE)
pathway<-c("GALECTIN")
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

comparison = c(2, 1)
dataset.name <- names(cellchat@net)
group.names0 <- df$source.target

xintercept = seq(0.5 + length(dataset.name[comparison]),
    length(group.names0) * length(dataset.name[comparison]),
    by = length(dataset.name[comparison]))
g <- g + geom_vline(xintercept = xintercept, linetype = "dashed",
    color = "grey60", size = 0.2)

group <- 1:length(comparison)
names(group) <- dataset.name[comparison]
color <- c('#2F6BFF','#E74C3C')

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

ggsave(paste0(outdir,"NK_netVisual_bubble.pdf"), width= 8, height= 5, units='in')
ggsave(paste0(outdir,"NK_netVisual_bubble.png"), width= 8, height= 5, units='in')

ggsave(paste0(outdir,"Myeloid_netVisual_bubble.pdf"), width= 12, height= 5, units='in')
ggsave(paste0(outdir,"Myeloid_netVisual_bubble.png"), width= 12, height= 5, units='in')

ggsave(paste0(outdir,"T_netVisual_bubble.pdf"), width= 15, height= 5, units='in')
ggsave(paste0(outdir,"T_netVisual_bubble.png"), width= 15, height= 5, units='in')