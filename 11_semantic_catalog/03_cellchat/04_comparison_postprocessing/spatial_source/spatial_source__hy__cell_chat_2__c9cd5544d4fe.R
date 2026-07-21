rm(list = ls())

library(cowplot)
library(Seurat)
library(RColorBrewer)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(patchwork)
library(ComplexHeatmap)

options(stringsAsFactors = FALSE)

outdir = "/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/core/make_band/cellchat/in"
inrds = "/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/core/make_band/cellchat/in/adata_celltype_anno.rds"
species = 'human'

setwd(outdir)

load("/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/core/make_band/cellchat/in/cellchat_object_list.RData")
load("/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/core/make_band/cellchat/in/cellchat_object.RData")

group = c("pCR","non-pCR")

plot_save = function(filename, plt, dpi=300, w=7, h=7){
  ggsave(paste0(filename,'.pdf'), plt, height = h, width = w)
}

fig_save = function(filename, plt, dpi=300, w=7, h=7){
  pdf(paste0(filename,'.pdf'), height = h, width = w)
  draw(plt)
  dev.off()
}

recor_save = function(filename, plt, dpi=300, w=7, h=7){
  pdf(paste0(filename,'.pdf'), height = h, width = w)
  replayPlot(plt)
  dev.off()
}

# 获取所有细胞类型（从 cellchat 对象中提取）
all_celltypes <- levels(cellchat@idents)
if (length(all_celltypes) == 0) {
  all_celltypes <- dimnames(cellchat@net$count)[[1]]
}
cat("检测到细胞类型：", paste(all_celltypes, collapse = ", "), "\n")

# 定义一个函数，将细胞类型名称转换为安全的文件名（去除或替换特殊字符）
safe_name <- function(name) {
  gsub("[^A-Za-z0-9_]", "_", name)
}

for (i in 1:length(group)) {
  for (j in (1:length(group))[(1:length(group)) != i]) {
    compare.list = c(i, j)
    
    # ========== 整体差异网络图（放大画布，避免拥挤） ==========
    # 整体图使用 15x15 英寸，标签字体 0.8
    pdf(paste0(group[i], "VS", group[j], "_Differential number of interactions.net.pdf"),
        height = 15, width = 15)
    netVisual_diffInteraction(cellchat, comparison = compare.list,
                              weight.scale = TRUE, measure = "count",
                              vertex.label.cex = 0.8)
    dev.off()
    
    pdf(paste0(group[i], "VS", group[j], "_Differential interaction strength.net.pdf"),
        height = 15, width = 15)
    netVisual_diffInteraction(cellchat, comparison = compare.list,
                              weight.scale = TRUE, measure = "weight",
                              vertex.label.cex = 0.8)
    dev.off()
    
    # ========== 为每个细胞类型作为发送方生成单独网络图 ==========
    # 画布尺寸 30x30 英寸，标签字体 0.5，节点最大尺寸 3（相对更小，间距更大）
    pdf_width <- 30
    pdf_height <- 30
    label_cex <- 0.5
    vertex_max_size <- 3
    
    for (src in all_celltypes) {
      # 生成安全的文件名（避免特殊字符导致无法创建文件）
      safe_src <- safe_name(src)
      cat("正在生成细胞类型 '", src, "' 的差异网络图...\n", sep = "")
      
      # 基于互作数量 (count)
      count_file <- paste0(group[i], "VS", group[j], "_", safe_src, "_outgoing_diff_count.pdf")
      tryCatch({
        pdf(count_file, height = pdf_height, width = pdf_width)
        netVisual_diffInteraction(cellchat,
                                  comparison = compare.list,
                                  measure = "count",
                                  sources.use = src,
                                  weight.scale = TRUE,
                                  vertex.label.cex = label_cex,
                                  vertex.size.max = vertex_max_size)
        dev.off()
        cat("  已保存：", count_file, "\n")
      }, error = function(e) {
        cat("  错误：无法生成", count_file, " - ", e$message, "\n")
        dev.off()  # 确保关闭设备
      })
      
      # 基于互作强度 (weight)
      weight_file <- paste0(group[i], "VS", group[j], "_", safe_src, "_outgoing_diff_weight.pdf")
      tryCatch({
        pdf(weight_file, height = pdf_height, width = pdf_width)
        netVisual_diffInteraction(cellchat,
                                  comparison = compare.list,
                                  measure = "weight",
                                  sources.use = src,
                                  weight.scale = TRUE,
                                  vertex.label.cex = label_cex,
                                  vertex.size.max = vertex_max_size)
        dev.off()
        cat("  已保存：", weight_file, "\n")
      }, error = function(e) {
        cat("  错误：无法生成", weight_file, " - ", e$message, "\n")
        dev.off()
      })
    }
    
    # ========== 原有热图和信息流排名图（保持原有尺寸） ==========
    gg5 = netVisual_heatmap(cellchat, comparison = compare.list)
    fig_save(paste0(group[i],"VS",group[j],"_Differential number of interactions.heatmap"), gg5)
    
    gg6 = netVisual_heatmap(cellchat, measure = "weight", comparison = compare.list)
    fig_save(paste0(group[i],"VS",group[j],"_Differential interaction strength.heatmap"), gg6)
    
    gg7 = rankNet(cellchat,
                  mode = "comparison",
                  measure = "weight",
                  sources.use = NULL,
                  targets.use = NULL,
                  comparison = compare.list,
                  stacked = TRUE,
                  do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Relative information flow.weight"), gg7)
    
    gg8 = rankNet(cellchat,
                  mode = "comparison",
                  measure = "weight",
                  sources.use = NULL,
                  targets.use = NULL,
                  comparison = compare.list,
                  stacked = FALSE,
                  do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Information flow.weight"), gg8)
    
    gg9 = rankNet(cellchat,
                  mode = "comparison",
                  measure = "count",
                  sources.use = NULL,
                  targets.use = NULL,
                  comparison = compare.list,
                  stacked = TRUE,
                  do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Relative information flow.count"), gg9)
    
    gg10 = rankNet(cellchat,
                   mode = "comparison",
                   measure = "count",
                   sources.use = NULL,
                   targets.use = NULL,
                   comparison = compare.list,
                   stacked = FALSE,
                   do.stat = TRUE)
    plot_save(paste0(group[i],"VS",group[j],"_Information flow.count"), gg10)
    
    dif_pathway = as.character(gg8$data[c(1:3, (nrow(gg8$data)-3):nrow(gg8$data)), 1])
    
    # ========== 比较群体间的传入/传出/总信号热图 ==========
    pathway.list = lapply(cellchat@netP, function(x){
      unique(x$pathways)
    }) %>% unlist() %>% unique()
    
    ht1 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "outgoing", 
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[i], 
                                            width = 10, height = length(pathway.list)/3)
    ht2 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "outgoing",
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[j],
                                            width = 10, height = length(pathway.list)/3)
    fig_save(paste0(group[i],"VS",group[j],'_Outgoing compare.heatmap'),
             plt = draw(ht1 + ht2, ht_gap = unit(0.5, "cm")),
             dpi = 300, w = 10, h = length(pathway.list)/3)
    
    ht3 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "incoming", 
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[i], 
                                            width = 10, height = length(pathway.list)/3)
    ht4 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "incoming",
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[j], 
                                            width = 10, height = length(pathway.list)/3)
    fig_save(paste0(group[i],"VS",group[j],'_Incoming compare.heatmap'),
             plt = draw(ht3 + ht4, ht_gap = unit(0.5, "cm")),
             dpi = 300, w = 10, h = length(pathway.list)/3)
    
    ht5 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[i]], pattern = "all", 
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[i], 
                                            width = 10, height = length(pathway.list)/3)
    ht6 = netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[j]], pattern = "all",
                                            signaling = pathway.list,
                                            title = names(cellchat.objet.list.new)[j], 
                                            width = 10, height = length(pathway.list)/3)
    fig_save(paste0(group[i],"VS",group[j],'_All compare.heatmap'),
             plt = draw(ht5 + ht6, ht_gap = unit(0.5, "cm")),
             dpi = 300, w = 10, h = length(pathway.list)/3)
  }
}
