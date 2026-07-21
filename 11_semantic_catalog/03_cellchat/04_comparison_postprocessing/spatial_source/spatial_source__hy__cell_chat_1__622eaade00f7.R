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

outdir <- "/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/core/make_band/cellchat/in"
setwd(outdir)

load("cellchat_object_list.RData")
load("cellchat_object.RData")

group <- c("pCR", "non-pCR")

# 辅助函数：安全文件名
safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)  # 合并多个下划线
  x <- gsub("^_|_$", "", x) # 去掉首尾下划线
  return(x)
}

# 获取所有细胞类型
all_celltypes <- levels(cellchat@idents)
if (length(all_celltypes) == 0) {
  all_celltypes <- dimnames(cellchat@net$count)[[1]]
}
cat("共有", length(all_celltypes), "个细胞类型：", paste(all_celltypes, collapse = ", "), "\n")

# 定义比较组合
comp_list <- list(c(1,2), c(2,1))

# 绘图参数
plot_params <- list(
  width = 20,           # 英寸
  height = 20,
  label_cex = 0.6,
  vertex_max_size = 4
)

# 循环处理每个比较组
for (comp in comp_list) {
  group1 <- group[comp[1]]
  group2 <- group[comp[2]]
  prefix <- paste0(group1, "VS", group2)
  cat("\n====== 正在处理比较组：", prefix, " ======\n")
  
  # ---------- 1. 整体差异网络图 ----------
  cat("  生成整体网络图...\n")
  # count
  pdf(paste0(prefix, "_Differential_number_of_interactions.net.pdf"),
      width = 15, height = 15)
  netVisual_diffInteraction(cellchat, comparison = comp,
                            measure = "count", weight.scale = TRUE,
                            vertex.label.cex = 0.8)
  dev.off()
  # weight
  pdf(paste0(prefix, "_Differential_interaction_strength.net.pdf"),
      width = 15, height = 15)
  netVisual_diffInteraction(cellchat, comparison = comp,
                            measure = "weight", weight.scale = TRUE,
                            vertex.label.cex = 0.8)
  dev.off()
  
  # ---------- 2. 为每个细胞类型生成单独图（同时保存单独文件和合并PDF） ----------
  # 对 count 和 weight 分别处理
  for (meas in c("count", "weight")) {
    cat("  正在处理 measure =", meas, "\n")
    
    # 合并 PDF 文件名
    combined_pdf <- paste0(prefix, "_all_celltypes_outgoing_diff_", meas, ".pdf")
    pdf(combined_pdf, width = plot_params$width, height = plot_params$height)
    
    for (src in all_celltypes) {
      safe_src <- safe_name(src)
      cat("    绘制发送方：", src, " ... ")
      
      # 单独文件（可选，如果不需要可注释掉）
      single_file <- paste0(prefix, "_", safe_src, "_outgoing_diff_", meas, ".pdf")
      tryCatch({
        # 先绘制到单独文件
        pdf(single_file, width = plot_params$width, height = plot_params$height)
        netVisual_diffInteraction(cellchat,
                                  comparison = comp,
                                  measure = meas,
                                  sources.use = src,
                                  weight.scale = TRUE,
                                  vertex.label.cex = plot_params$label_cex,
                                  vertex.size.max = plot_params$vertex_max_size)
        dev.off()
        
        # 再绘制到合并 PDF（当前设备）
        netVisual_diffInteraction(cellchat,
                                  comparison = comp,
                                  measure = meas,
                                  sources.use = src,
                                  weight.scale = TRUE,
                                  vertex.label.cex = plot_params$label_cex,
                                  vertex.size.max = plot_params$vertex_max_size)
        cat(" [成功]\n")
      }, error = function(e) {
        cat(" [失败：", e$message, "]\n")
        # 在合并 PDF 中写入错误提示
        plot.new()
        text(0.5, 0.5, paste("Error plotting", src, ":", e$message), cex = 1.2)
        # 确保设备关闭（如果单独文件出错，设备可能仍开着）
        try(dev.off(), silent = TRUE)
      })
    }
    
    dev.off()  # 关闭合并 PDF 设备
    cat("    合并 PDF 已保存：", combined_pdf, "\n")
  }
  
  # ---------- 3. 原有的热图和排名图（保持不变） ----------
  cat("  生成热图和排名图...\n")
  # 差异热图
  gg5 <- netVisual_heatmap(cellchat, comparison = comp)
  pdf(paste0(prefix, "_Differential_number_of_interactions.heatmap.pdf"), height = 7, width = 7)
  draw(gg5)
  dev.off()
  
  gg6 <- netVisual_heatmap(cellchat, measure = "weight", comparison = comp)
  pdf(paste0(prefix, "_Differential_interaction_strength.heatmap.pdf"), height = 7, width = 7)
  draw(gg6)
  dev.off()
  
  # 信息流排名图（使用 ggsave 保存）
  gg7 <- rankNet(cellchat, mode = "comparison", measure = "weight",
                 comparison = comp, stacked = TRUE, do.stat = TRUE)
  ggsave(paste0(prefix, "_Relative_information_flow.weight.pdf"), gg7, width = 7, height = 7)
  
  gg8 <- rankNet(cellchat, mode = "comparison", measure = "weight",
                 comparison = comp, stacked = FALSE, do.stat = TRUE)
  ggsave(paste0(prefix, "_Information_flow.weight.pdf"), gg8, width = 7, height = 7)
  
  gg9 <- rankNet(cellchat, mode = "comparison", measure = "count",
                 comparison = comp, stacked = TRUE, do.stat = TRUE)
  ggsave(paste0(prefix, "_Relative_information_flow.count.pdf"), gg9, width = 7, height = 7)
  
  gg10 <- rankNet(cellchat, mode = "comparison", measure = "count",
                  comparison = comp, stacked = FALSE, do.stat = TRUE)
  ggsave(paste0(prefix, "_Information_flow.count.pdf"), gg10, width = 7, height = 7)
  
  # 信号角色热图
  pathway.list <- unique(unlist(lapply(cellchat@netP, function(x) unique(x$pathways))))
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[comp[1]]], pattern = "outgoing",
                                           signaling = pathway.list, title = group[comp[1]],
                                           width = 10, height = length(pathway.list)/3)
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[comp[2]]], pattern = "outgoing",
                                           signaling = pathway.list, title = group[comp[2]],
                                           width = 10, height = length(pathway.list)/3)
  pdf(paste0(prefix, "_Outgoing_compare.heatmap.pdf"), width = 10, height = length(pathway.list)/3)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  ht3 <- netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[comp[1]]], pattern = "incoming",
                                           signaling = pathway.list, title = group[comp[1]],
                                           width = 10, height = length(pathway.list)/3)
  ht4 <- netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[comp[2]]], pattern = "incoming",
                                           signaling = pathway.list, title = group[comp[2]],
                                           width = 10, height = length(pathway.list)/3)
  pdf(paste0(prefix, "_Incoming_compare.heatmap.pdf"), width = 10, height = length(pathway.list)/3)
  draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  ht5 <- netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[comp[1]]], pattern = "all",
                                           signaling = pathway.list, title = group[comp[1]],
                                           width = 10, height = length(pathway.list)/3)
  ht6 <- netAnalysis_signalingRole_heatmap(cellchat.objet.list.new[[comp[2]]], pattern = "all",
                                           signaling = pathway.list, title = group[comp[2]],
                                           width = 10, height = length(pathway.list)/3)
  pdf(paste0(prefix, "_All_compare.heatmap.pdf"), width = 10, height = length(pathway.list)/3)
  draw(ht5 + ht6, ht_gap = unit(0.5, "cm"))
  dev.off()
}

cat("\n========== 所有图形生成完毕！==========\n")
cat("请检查当前目录下的 PDF 文件，包括：\n")
cat("  - 整体差异网络图：*_Differential_*.net.pdf\n")
cat("  - 每个细胞类型的单独图：*_outgoing_diff_*.pdf\n")
cat("  - 合并的多页 PDF：*_all_celltypes_outgoing_diff_*.pdf\n")
cat("  - 热图和排名图等\n")
