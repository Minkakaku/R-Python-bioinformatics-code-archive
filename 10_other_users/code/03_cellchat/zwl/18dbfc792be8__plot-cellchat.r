setwd('/mnt/NexData_file/home/zwl/Project/SC2601160389')
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(ComplexHeatmap)
library(argparser)
# library(pdftools)
library(stringr)
library(jpeg)
library(cowplot)
library(Seurat)
library(RColorBrewer)
rm(list = ls())
gc()

################################################################################
##################cellchat-heatmap
################################################################################
load("/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/Stent_Con/tmp_object.RData")
load("/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/Stent2_Con2/tmp_object.RData")

library(ComplexHeatmap)
library(circlize)

# 设定保存路径
save_path <- "/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/"

# 开始循环遍历 object.list
for (i in seq_along(object.list)) {
  
  # 1. 获取当前对象的名称和矩阵数据
  obj_name <- names(object.list)[i]
  count_mat_Con <- object.list[[i]]@net$count
  
  # 2. 矩阵预处理
  count_mat <- count_mat_Con + t(count_mat_Con)
  count_mat[upper.tri(count_mat)] <- NA
  count_mat <- as.matrix(count_mat)
  
  # 3. 定义热图绘制函数 (方便重复调用生成不同格式)
  draw_ht <- function() {
    # 注意：在循环中直接调用 Heatmap 可能不显式画图，需要用 draw() 包裹
    ht <- Heatmap(
      count_mat,
      name = "Interaction",
      col = colorRampPalette(c("#4682B4", "#ffffbf", "firebrick"))(100),
      na_col = "white",
      border = FALSE, 
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      rect_gp = gpar(col = "white", lwd = 1),
      row_names_side = "left",
      column_names_side = "bottom",
      column_names_rot = 45,
      cell_fun = function(j, i, x, y, width, height, fill) {
        if(!is.na(count_mat[i, j])) {
          grid.text(sprintf("%d", count_mat[i, j]), x, y, gp = gpar(fontsize = 10))
        }
      },
      heatmap_legend_param = list(title = "Count")
    )
    draw(ht)
  }
  
  # 4. 构建文件名基础
  file_base <- paste0('13.',obj_name, "_Number_Interaction_heatmap")
  
  # 5. 保存为 PDF
  pdf(file = paste0(save_path, file_base, ".pdf"), width = 8, height = 7)
  draw_ht()
  dev.off()
  
  # 6. 保存为 PNG
  png(filename = paste0(save_path, file_base, ".png"), width = 800, height = 700, res = 100)
  draw_ht()
  dev.off()
  
  message(paste0("Done: ", obj_name))
}

################################################################################
##################cellchat-dotplot
################################################################################

load("/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/Stent_Con/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/"
print(colnames(object.list[[1]]@net$count))


p1 = object.list %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        source %in%  c("B cells","Endothelial cells","Erythroid cells","Mast cells","Myeloids",
                       "Neutrophils","Plasma cells","Schwann cells","Smooth muscle cells","T cells"),
        target  %in% c("Fib_other","sC6"),
        # interaction_name_2 %in% pathway,
        # target  %in% c('KC Trem2+', 'HSC'),
        # pathway_name %in% c('APP','LAIR1'),
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
  theme(strip.background = element_rect(colour = NA),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) +
  scale_color_gradientn(colours  = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu"))) +
  scale_size(transform = "reverse", range = c(1, 3)) +
  Seurat::RotatedAxis()
p1
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/14.cellchat_Stent_Con_dotplot.png',p1,w=30,h=15,dpi=300)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/14.cellchat_Stent_Con_dotplot.pdf',p1,w=30,h=15,dpi=300)

###################################################################
load("/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/Stent2_Con2/tmp_object.RData")
cellchat <- mergeCellChat(object.list,add.names = names(object.list))
outdir<-"/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/"
print(colnames(object.list[[1]]@net$count))

p1 = object.list %>% 
  map(subsetCommunication) %>%
  map(function(x){
    x %>% 
      dplyr::filter(
        source %in%  c("B cells","Endothelial cells","Erythroid cells","Mast cells","Myeloids",
                       "Neutrophils","Plasma cells","Schwann cells","Smooth muscle cells","T cells"),
        target  %in% c("Fib_other","sC6"),
        # interaction_name_2 %in% pathway,
        # target  %in% c('KC Trem2+', 'HSC'),
        # pathway_name %in% c('APP','LAIR1'),
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
  theme(strip.background = element_rect(colour = NA),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) +
  scale_color_gradientn(colours  = rev(RColorBrewer::brewer.pal(Inf, "RdYlBu"))) +
  scale_size(transform = "reverse", range = c(1, 3)) +
  Seurat::RotatedAxis()
p1
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/14.cellchat_Stent2_Con2_dotplot.png',p1,w=30,h=15,dpi=300)
ggsave(filename = '/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/14.cellchat_Stent2_Con2_dotplot.pdf',p1,w=30,h=15,dpi=300)

################################################################################
##################cellchat-和弦图
################################################################################

print(colnames(object.list[[1]]@net$count))

color_use = c(
  "firebrick",   # 耐火砖红
  "#1F77B4",
  "#FF7F0E",
  "#2CA02C",
  "#9467BD",
  "#8C564B",
  "#E377C2",
  "#7F7F7F",
  "#BCBD22",
  "#17BECF",
  "#AEC7E8",
  "#FFBB78"
)
outdir <- "/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/"
load("/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/Stent_Con/tmp_object.RData")
# 保存PNG  保存PDF
png(paste0(outdir, "15.cellchat_",names(object.list)[1],"_TGFb_PERIOSTIN_chord_diagram.png"), width = 2800, height = 2000, res = 300)
netVisual_chord_gene(object.list[[1]], sources.use = 1:12, targets.use = 1:12, 
                     signaling = c("TGFb","PERIOSTIN"),color.use = color_use,legend.pos.x = 8)
dev.off()

pdf(paste0(outdir, "15.cellchat_",names(object.list)[1],"_TGFb_PERIOSTIN_chord_diagram.pdf"),  width = 12, height = 8)
netVisual_chord_gene(object.list[[1]], sources.use = 1:12, targets.use = 1:12, 
                     signaling = c("TGFb","PERIOSTIN"),color.use = color_use,legend.pos.x = 8)
dev.off()

# 保存PNG  保存PDF
png(paste0(outdir, "15.cellchat_",names(object.list)[2],"_TGFb_PERIOSTIN_chord_diagram.png"), width = 2800, height = 2000, res = 300)
netVisual_chord_gene(object.list[[2]], sources.use = 1:12, targets.use = 1:12, 
                     signaling = c("TGFb","PERIOSTIN"),color.use = color_use,legend.pos.x = 8)
dev.off()

pdf(paste0(outdir, "15.cellchat_",names(object.list)[2],"_TGFb_PERIOSTIN_chord_diagram.pdf"),  width = 12, height = 8)
netVisual_chord_gene(object.list[[2]], sources.use = 1:12, targets.use = 1:12, 
                     signaling = c("TGFb","PERIOSTIN"),color.use = color_use,legend.pos.x = 8)
dev.off()



load("/mnt/NexData_file/home/zwl/Project/SC2601160389/售后3/cellchat/Stent2_Con2/tmp_object.RData")
# 保存PNG  保存PDF
png(paste0(outdir, "15.cellchat_",names(object.list)[1],"_TGFb_PERIOSTIN_chord_diagram.png"), width = 2800, height = 2000, res = 300)
netVisual_chord_gene(object.list[[1]], sources.use = 1:12, targets.use = 1:12, 
                     signaling = c("TGFb","PERIOSTIN"),color.use = color_use,legend.pos.x = 8)
dev.off()

pdf(paste0(outdir, "15.cellchat_",names(object.list)[1],"_TGFb_PERIOSTIN_chord_diagram.pdf"),  width = 12, height = 8)
netVisual_chord_gene(object.list[[1]], sources.use = 1:12, targets.use = 1:12, 
                     signaling = c("TGFb","PERIOSTIN"),color.use = color_use,legend.pos.x = 8)
dev.off()

# 保存PNG  保存PDF
png(paste0(outdir, "15.cellchat_",names(object.list)[2],"_TGFb_PERIOSTIN_chord_diagram.png"), width = 2800, height = 2000, res = 300)
netVisual_chord_gene(object.list[[2]], sources.use = 1:12, targets.use = 1:12, 
                     signaling = c("TGFb","PERIOSTIN"),color.use = color_use,legend.pos.x = 8)
dev.off()

pdf(paste0(outdir, "15.cellchat_",names(object.list)[2],"_TGFb_PERIOSTIN_chord_diagram.pdf"),  width = 12, height = 8)
netVisual_chord_gene(object.list[[2]], sources.use = 1:12, targets.use = 1:12, 
                     signaling = c("TGFb","PERIOSTIN"),color.use = color_use,legend.pos.x = 8)
dev.off()

