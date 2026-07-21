library(Seurat)
library(dplyr)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(gt)
setwd('G:\\FWSC20242121-05\\售后1\\DEG_FunctionEnrichment')
# === 读取数据 ===
adata <- readRDS("result.rds")
#adata <- subset(adata, idents = c("IgG Plasma C1", "IgG Plasma C2"))

# === 整理元数据 ===
#adata$cell_type <- droplevels(adata$cell_type)
adata$group <- as.character(adata$group)
adata$sample <- as.character(adata$sample)

adata$group <- as.factor(adata$group)
adata$sample <- as.factor(adata$sample)

# 设定身份标识为 group
Idents(adata) <- "group"

# === 设置输出路径和比较组 ===
outdir <- "G:\\FWSC20242121-05\\售后1\\DEG_FunctionEnrichment\\"
#dir.create(outdir, recursive = TRUE)

diffgroup <- c("HC/PTB")
diffgroup <- strsplit(diffgroup, "/")  # 得到 list，每个元素是 [对照组, 对比组]

# === 遍历每个 cell_type ===
for (cluster_diff_sub in unique(adata$cell_type)) {
  message("正在处理细胞类型：", cluster_diff_sub)
  tryCatch({
    # 仅保留当前细胞类型
    barcodes <- rownames(adata@meta.data[adata$cell_type == cluster_diff_sub, ])
    x <- subset(adata, cells = barcodes)
    #x$cell_type <- droplevels(x$cell_type)
    
    # 遍历每一组差异比较
    for (m in seq_along(diffgroup)) {
      diff1 <- diffgroup[[m]][1]
      diff2 <- diffgroup[[m]][2]
      
      # 差异分析
      markers <- FindMarkers(x, ident.1 = diff1, ident.2 = diff2)
      markers$gene <- rownames(markers)
      markers <- select(markers, gene, everything())
      
      # 保存全部结果
      dir.create(paste0(outdir, cluster_diff_sub, "\\"), recursive = TRUE)
      file_all <- paste0(outdir, cluster_diff_sub, "\\", diff1, "-", diff2, "_diffgenes.csv")
      write.table(markers, file = file_all, sep = ",", row.names = FALSE, quote = FALSE)
      
      # 保存显著差异基因结果
      markers_sig <- markers %>%
        filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25, pct.1 > pct.2, pct.1 > 0.1)
      file_sig <- paste0(outdir, cluster_diff_sub, "\\", diff1, "-", diff2, "_diffgenes_significant.csv")
      write.table(markers_sig, file = file_sig, sep = ",", row.names = FALSE, quote = FALSE)
      
      message("完成：", cluster_diff_sub, " - ", diff1, " vs ", diff2)
    }
  }, error = function(e) {
    message("❌ 错误（", cluster_diff_sub, "）: ", e$message)
  })
}

cluster_diff=unique(adata$cell_type)
for(n in 1:length(cluster_diff)){
  cluster_diff_sub<-cluster_diff[[n]]
  scale_row_no_extreme <- function(mat, cap = 2) {
  mat_scaled <- t(scale(t(mat)))  # 按行标准化
  mat_scaled[mat_scaled > cap] <- cap
  mat_scaled[mat_scaled < -cap] <- -cap
  return(mat_scaled)
}

# 处理表达矩阵，排除极值影响
expr_matrix_scaled <- scale_row_no_extreme(expr_matrix)

  print(cluster_diff_sub)
  for(m in 1:length(diffgroup)){
    diff1<-diffgroup[[m]][1]
    diff2<-diffgroup[[m]][2]
    print(paste0(diff1,"vs",diff2))
    data = read.table(paste0(outdir,cluster_diff_sub,"\\",diff1,"-",diff2, "_diffgenes.csv"),sep = ",",head=TRUE)
    colnames(data) <- c("gene_name", "PValue", "log2FC", "pct_nz_group" , "pct_nz_reference", "FDR")
    data$log10PValue = -log10(data$PValue)
    
    data$log10PValue[is.infinite(data$log10PValue)] = 300
    cut_off_PValue =0.05 #设置FDR的阈值
    cut_off_log2FC =0.25 #设置log2FC的阈值
    data$Sig = ifelse(data$PValue < cut_off_PValue &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                        abs(data$log2FC) >= cut_off_log2FC,  #abs绝对值
                      ifelse(data$log2FC > cut_off_log2FC ,'Up','Down'),'no')
    data = data.frame(data)
    table(data$Sig)
    # data = data %>%
    #   filter(pct_nz_group >= 0.2) %>%
    #   filter(pct_nz_reference <= 0.6) %>%
    #   arrange(desc(abs(log2FC)),FDR)
    xlim_min<-min(data$log2FC)-0.3
    xlim_max<-max(data$log2FC)+0.3
    ylim_max <- min(ceiling(max(data$log10PValue) / 50) * 50,300)
    ###绘图——基础火山图###
    p1 <- ggplot(data, aes(x =log2FC, y=log10PValue, colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
      geom_point(alpha=0.65, size=3) +  #点的透明度、大小
      scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757")) + 
      geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
      geom_hline(yintercept = -log10(cut_off_PValue), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
      labs(x="log2FC", y="-log10PValue") +  #x、y轴标签
      ggtitle("Volcano Plot") + #标题
      scale_x_continuous(limits = c(xlim_min, xlim_max)) +
      scale_y_continuous(limits = c(0, ylim_max), breaks = seq(0, ylim_max, 50)) +
      theme_bw() + # 主题，help(theme)查找其他个性化设置
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 30),
            legend.text = element_text(size = 15),
            legend.position="right", 
            legend.title = element_blank()
      ) 
    
    Up_top_10 =(     #筛选差异显著上调的前10个Gene
      data %>%
        filter(Sig == 'Up') %>%
        filter(pct_nz_group >= 0.1) %>%
        arrange(desc(abs(log2FC)),PValue) %>%
        head(10)
    )
    Up_top_10 %>% gt() #数据制成表
    Down_top_10 = (       #筛选差异显著下调的前10个Gene
      data %>%
        filter(Sig == 'Down') %>%
        filter(pct_nz_reference >= 0.1) %>%
        arrange( desc(abs(log2FC)),PValue) %>%
        head(10)
    )
    Down_top_10 %>% gt() #数据制成表
    
    p2 <- p1 +
      geom_label_repel(data = Up_top_10,
                       aes(log2FC, log10PValue, label = gene_name),
                       size = 3, fill="#CCFFFF",
                       alpha = 0.65, color = "black")+
      geom_label_repel(data = Down_top_10,
                       aes(log2FC, log10PValue, label = gene_name),
                       size = 3, fill="#FFCCCC",
                       alpha = 0.65, color = "black")
    
    p2 #出图
    ggsave(paste0(cluster_diff_sub,"\\",diff1,"-",diff2, ".Volcanoplot.pdf"),width = 10,height = 10)
    ggsave(paste0(cluster_diff_sub,"\\",diff1,"-",diff2, ".Volcanoplot.png"),width = 10,height = 10)
    
    
    # 合并上下调共20个基因
    top_genes <- union(Up_top_10$gene_name, Down_top_10$gene_name)
    
    # 提取表达矩阵并筛选基因
    expr_matrix <- as.matrix(x[["RNA"]]@data[top_genes, ])
    
    # 获取元信息
    meta_data <- x@meta.data
    meta_data <- meta_data[colnames(expr_matrix), ]
    
    # 仅保留 diff1 和 diff2 两组的细胞
    cells_use <- meta_data$group %in% c(diff1, diff2)
    expr_matrix <- expr_matrix[, cells_use]
    meta_data <- meta_data[cells_use, , drop = FALSE]
    
    # 横轴排序：diff1 在前，diff2 在后
    meta_data$group <- factor(meta_data$group, levels = c(diff1, diff2))
    meta_data <- meta_data[order(meta_data$group), ]
    expr_matrix <- expr_matrix[, rownames(meta_data)]
    scale_row_no_extreme <- function(mat, cap = 2) {
      mat_scaled <- t(scale(t(mat)))  # 按行标准化
      mat_scaled[mat_scaled > cap] <- cap
      mat_scaled[mat_scaled < -cap] <- -cap
      return(mat_scaled)
    }
    
    # 处理表达矩阵，排除极值影响
    expr_matrix_scaled <- scale_row_no_extreme(expr_matrix)
    
    
    # 列注释
    ann_col <- data.frame(Group = meta_data$group)
    rownames(ann_col) <- rownames(meta_data)
      # 保存热图
      pdf(paste0(cluster_diff_sub, "/", diff1, "-", diff2, ".Top20_Heatmap.pdf"), width = 10, height = 8)
      pheatmap(expr_matrix_scaled,
               scale = "none",
               cluster_rows = TRUE,
               cluster_cols = FALSE,  # 不对列聚类
               show_colnames = FALSE,
               annotation_col = ann_col,
               fontsize_row = 8,
               main = paste0(cluster_diff_sub, ": ", diff1, " vs ", diff2))
      dev.off()
      
      png(filename = paste0(cluster_diff_sub, "/", diff1, "-", diff2, ".Top20_Heatmap.png"),
          width = 1200, height = 1000, res = 150)
      pheatmap(expr_matrix_scaled,
               scale = "none",
               cluster_rows = TRUE,
               cluster_cols = FALSE,  # 不对列聚类
               show_colnames = FALSE,
               annotation_col = ann_col,
               fontsize_row = 8,
               main = paste0(cluster_diff_sub, ": ", diff1, " vs ", diff2))
      dev.off()
    
    
  }}

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)




