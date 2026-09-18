rm(list = ls())
gc()

library(Seurat)
library(slingshot)
library(SingleCellExperiment)
#library(qs)
library(tidyverse)
library(RColorBrewer)
library(tradeSeq)

result_path <- "/mnt/NexData_file/home/xzx/project/SC2512257794/shouhou4"
dir.create(result_path)
setwd(result_path)

sc <- readRDS("../shouhou2/anno_allsub.rds")
table(sc$cell_type)
DimPlot(sc,group.by = 'cell_type')

### 提取需要的矩阵构建SingleCellExperiment对象
scale.data <- sc@assays$RNA@counts ###提取基因 按需求 直接全部基因
scale.gene <- rownames(scale.data)

counts <- sc@assays$RNA@counts
counts <- counts[scale.gene,]

# 将表达矩阵转换为SingleCellExperiment对象
# 输入需要counts矩阵，否则影响下游分析
sim <- SingleCellExperiment(assays = List(counts = counts)) 

# umap reduction
umap = sc@reductions$umap@cell.embeddings
colnames(umap) = c('UMAP-1', 'UMAP-2')
# 将降维的结果添加到SingleCellExperiment对象中
reducedDims(sim) = SimpleList(UMAP = umap)

# metadata
meta = sc@meta.data
# colData(sim)相当于meta.data，但他不是data.frame格式
# 所以需要分步赋予
colData(sim)$sample = meta$sample
colData(sim)$cell_type = meta$cell_type
colData(sim)$group = meta$group


#plot(umap, col = rgb(0,0,0,.5), pch=16, asp = 1)
sim

sim <- slingshot(sim, 
                 clusterLabels = 'cell_type',  # 选择colData中细胞注释的列名
                 reducedDim = 'UMAP',  
                 start.clus= "Basal_cells_c3",  # 选择起点
                 end.clus = NULL     # 这里我们不指定终点
)
colnames(colData(sim))

# 基础轨迹图：点=细胞，线=推断的发育轨迹
plot(
  reducedDims(sim)$UMAP,         # UMAP坐标
  pch = 16,                 # 点形状（16=实心圆）
  asp = 1                   # 纵横比=1，避免图形拉伸
)
# 添加轨迹线（lwd=线宽，col=配色）
lines(SlingshotDataSet(sim), lwd = 2, col = brewer.pal(9, "Set1"))
# 添加图例
legend(
  "left",  legend = c("Lineage1","Lineage2","Lineage3"),      # 轨迹名称（多轨迹时需对应修改）
  col = unique(brewer.pal(6, "Set1")),
  inset = 0.8,                    # 图例位置（避免遮挡细胞）
  pch = 16
)


cell_type_colors <- c(
  'Basal_cells_c1'=	'#1f77b4',
  'Basal_cells_c2'=	'#ff7f0e',
  'Basal_cells_c3'=	'#279e68',
  'Ciliated_cells'=	'#d62728',
  'Club_cells'=	'#aa40fc',
  'Deuterosomal_cells'=	'#8c564b',
  'Goblet_cells_c1'=	'#e377c2',
  'Goblet_cells_c2'=	'#b5bd61',
  'Goblet_cells_c3'=	'#17becf',
  'Ionocyte'=	'#aec7e8',
  'Suprabasal_cells_c1'=	'#ffbb78',
  'Suprabasal_cells_c2'=	'#98df8a',
  'Tuft cells'=	'#ff9896')

cell_types <- as.vector(sim$cell_type)  # 替换成你的细胞类型

# --------------- 可视化类型 1 ---------------
cell_colors_clust <- cell_type_colors[cell_types]

pdf(file = "./slingshot/01_lineage_type1.pdf",width = 7,height = 7)
plot(reducedDims(sim)$UMAP,  pch = 16, cex = 0.5, col =cell_colors_clust )
lines(
  SlingshotDataSet(sim),
  lwd = 3,
  type = "lineages",
  col ='black'
)
dev.off()

png(file = "./slingshot/01_lineage_type1.png",width = 700,height = 700,units = 'px')
plot(reducedDims(sim)$UMAP,  pch = 16, cex = 0.5, col =cell_colors_clust )
lines(
  SlingshotDataSet(sim),
  lwd = 3,
  type = "lineages",
  col ='black'
)
dev.off()

pdf(file = "./slingshot/01_lineage_type2.pdf",width = 7,height = 7)
plot(reducedDims(sim)$UMAP,  pch = 16, cex = 0.5, col =cell_colors_clust )
lines(SlingshotDataSet(sim),lwd=2,col=brewer.pal(9,"Set1"))
legend("topleft",legend=paste0("lineage",1:3),col=unique(brewer.pal(6,"Set1")),inset=0.8,pch=16)
dev.off()

png(file = "./slingshot/01_lineage_type2.png",width = 700,height = 700,units = 'px')
plot(reducedDims(sim)$UMAP,  pch = 16, cex = 0.5, col =cell_colors_clust )
lines(SlingshotDataSet(sim),lwd=2,col=brewer.pal(9,"Set1"))
legend("topleft",legend=paste0("lineage",1:3),col=unique(brewer.pal(6,"Set1")),inset=0.8,pch=16)
dev.off()



# --------------- 可视化类型 2 按轨迹绘制伪时间 ---------------

# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) # 把这些颜色变成100个梯度，模拟渐变色
# plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)] # 这里我们用cut函数把 lineage2 分割成100个区间，同一区间的细胞给与同一个颜色
# plotcol[is.na(plotcol)] <- "lightgrey" # 不属于 lineage3 的点为NA，我们把他们变成灰色
# plotcol
# pdf(file = "./slingshot/02_lineage1.pdf",width = 7,height = 5)
# plot(reducedDims(sim)$UMAP, col = plotcol, pch=16, asp = 1)
# lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))
# legend("left",
#        legend = paste0("lineage",1:3),
#        col = unique(brewer.pal(6,"Set1")),
#        inset=0.8,
#        pch = 16)
# dev.off()
# 
# png(file = "./slingshot/02_lineage1.png",width = 700,height = 500)
# plot(reducedDims(sim)$UMAP, col = plotcol, pch=16, asp = 1)
# lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(9,"Set1"))
# legend("left",
#        legend = paste0("lineage",1:3),
#        col = unique(brewer.pal(6,"Set1")),
#        inset=0.8,
#        pch = 16)
# dev.off()


#============ 批量绘制 3 条 lineage 伪时间图 =====================
# 获取总谱系数量（自动识别 3 条）
n_lineages <- length(slingLineages(sim))

# 循环绘制每一条 lineage
for (lin in 1:n_lineages) {
  
  # 1. 获取当前 lineage 的伪时间（其他谱系=NA）
  pseudo_col <- paste0("slingPseudotime_", lin)
  pseudotime <- sim[[pseudo_col]]  # 提取对应 lineage 伪时间
  
  # 2. 生成渐变色（100 个梯度）
  colors <- colorRampPalette(brewer.pal(11, 'Spectral')[-6])(100)
  
  # 3. 给细胞上色：当前谱系渐变，其他灰色
  plotcol <- colors[cut(pseudotime, breaks = 100)]
  plotcol[is.na(plotcol)] <- "lightgrey"  # 非当前谱系 = 灰色
  
  # 4. 固定 3 条轨迹颜色（红/蓝/绿，清晰不冲突）
  lineage_cols <- c("#E41A1C", "#377EB8", "#4DAF4A")  # lineage1/2/3
  
  # ===================== 输出 PDF =====================
  pdf(file = paste0("./slingshot/02_lineage", lin, ".pdf"), width = 7, height = 5)
  plot(reducedDims(sim)$UMAP, 
       col = plotcol, 
       pch = 16, 
       asp = 1,
       main = paste0("Lineage ", lin, " Pseudotime"))  # 标题自动标注
  
  # 绘制 3 条轨迹（固定颜色）
  lines(SlingshotDataSet(sim), lwd = 2, col = lineage_cols)
  
  # 图例（正确、不遮挡）
  legend("topright",
         legend = paste0("Lineage ", 1:3),
         col = lineage_cols,
         lwd = 2,
         bty = "n",
         cex = 0.9)
  dev.off()
  
  # ===================== 输出 PNG =====================
  png(file = paste0("./slingshot/02_lineage", lin, ".png"), width = 700, height = 500)
  plot(reducedDims(sim)$UMAP, 
       col = plotcol, 
       pch = 16, 
       asp = 1,
       main = paste0("Lineage ", lin, " Pseudotime"))
  
  lines(SlingshotDataSet(sim), lwd = 2, col = lineage_cols)
  
  legend("topright",
         legend = paste0("Lineage ", 1:3),
         col = lineage_cols,
         lwd = 2,
         bty = "n",
         cex = 0.9)
  dev.off()
}


########
# 提取计数矩阵（SCE对象的assay数据）
counts <- sim@assays@data$counts
# 转换为SlingshotDataSet对象（tradeSeq的标准输入）
#SDS <- SlingshotDataSet(sim)
# 1. 确定GAM模型的最佳k值（k= knots数量，控制模型灵活性）
# 随机抽取 200~300 个细胞（足够判断k值，速度暴快）
# ==============================================================
library(dplyr)

meta$cell_id <- rownames(meta)
set.seed(123)

cells_use <- meta %>%
  group_by(cell_type) %>%
  slice_head(n = 200) %>%   # 每组最多200，不够自动取全部
  pull(cell_id)

sim_sampled <- sim[, cells_use]
counts_sampled <- counts[, cells_use]
sim_sampled <- slingshot(
  sim_sampled,
  clusterLabels = sim_sampled$cell_type,
  reducedDim = "UMAP",
  start.clus= "Basal_cells_c3",  # 选择起点
  end.clus = NULL     # 这里我们不指定终点
)
#plot(reducedDims(sim_sampled)$UMAP,  pch = 16, cex = 0.5, col =cell_colors_clust )
#lines(SlingshotDataSet(sim_sampled),lwd=2,col=brewer.pal(9,"Set1"))
SDS_sampled <- SlingshotDataSet(sim_sampled)

icMat <- evaluateK(
  counts = counts_sampled,
  sds = SDS_sampled,
  k = 5:10,
  nGenes = 300,
  verbose = TRUE
)


# 后续用icMat选择最小AIC的k值（如k=4）
# 2. 获取伪时间和细胞权重
SDS <- SlingshotDataSet(sim)
pseudotime <- slingPseudotime(SDS, na = FALSE)      # 每个细胞的伪时间
cell_Weight <- slingCurveWeights(SDS)           # 细胞在各轨迹的权重（分支时用）
# 3. 【演示加速】选择TOP200高变异基因
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc <- FindVariableFeatures(
  object = sc,
  selection.method = "vst",
  nfeatures = 500,  # 调整为目标基因数量
  verbose = FALSE
)
hvgs <- VariableFeatures(sc)
#top_genes <- head(order(rowVars(counts), decreasing = TRUE), 500)
top_counts <- counts[hvgs, ]
#top_counts <- counts

# 拟合广义加性模型（GAM）：描述基因表达随伪时间的非线性变化

system.time({
  sce <- fitGAM(
    counts = top_counts,
    pseudotime = pseudotime,
    cellWeights = cell_Weight,
    nknots = 6,  # 【关键】knots数量，基于evaluateK结果调整.如果这个点小于6，则选择对应的数字；如果这个点大于等于6，则选择6
    verbose = TRUE
  )
})
table(rowData(sce)$tradeSeq$converged)

# 提取收敛状态为 TRUE 的基因索引
converged_genes <- rowData(sce)$tradeSeq$converged == TRUE
# 对 sce 对象进行子集筛选，只保留收敛的基因
sce_converged <- sce[converged_genes, ]

# 1. 分析基因与伪时间的整体相关性（associationTest）
Relation <- associationTest(sce_converged)  # p值越小，基因与伪时间关联越强
# 2. 分析轨迹起止点的差异基因（startVsEndTest）
sET <- startVsEndTest(sce_converged,lineages = TRUE)  # 比较“起始细胞”与“终止细胞”的基因表达
# 关键指标：waldStat（值越大，起止点差异越显著）、logFClineage1（起始→终止的log2倍变化）
sET %>%
  tibble::rownames_to_column("gene") %>%
  write_csv("./slingshot/Pseudotime_ordergene.csv")
# 3. 筛选top差异基因（以waldStat排序）
# order_sET <- order(sET$waldStat, decreasing = TRUE)
# topGeneStart <- names(sce)[order_sET[1]]  # 差异最显著的第一个基因（如示例中的"FTL"）

# 1. 确保两个结果的行名（基因名）是一致的
# Relation 来自 associationTest, sET 来自 startVsEndTest



#################### 3. 可视化 ###########################
# 2. 合并两个结果
combined_res <- merge(as.data.frame(Relation), 
                      as.data.frame(sET), 
                      by = "row.names", 
                      suffixes = c("_AT", "_sET"))

# 将第一列（Row.names）重新设为行名
rownames(combined_res) <- combined_res$Row.names
combined_res <- combined_res[, -1]

# 3. 筛选 Lineage 1 中显著变化的基因（不分上下调）
target_genes_L1_all <- combined_res %>%
  filter(pvalue_AT < 0.05 & pvalue_lineage1 < 0.05) %>%
  # 计算 logFC 的绝对值，并结合 waldStat 排序
  filter(abs(logFClineage1) > 1) %>%
  arrange(desc(waldStat_lineage1)) %>%
  head(20) %>% rownames()

# # 4. 筛选 Lineage 1 中显著上升的前 5 个基因
# top_up_L1 <- combined_res %>%
#   filter(pvalue_AT < 0.05 & pvalue_lineage1 < 0.05 & logFClineage1 > 1) %>%
#   arrange(desc(waldStat_lineage1)) %>%
#   head(5) %>% rownames()
# 
# # 5. 筛选 Lineage 1 中显著下降的前 5 个基因
# top_down_L1 <- combined_res %>%
#   filter(pvalue_AT < 0.05 & pvalue_lineage1 < 0.05 & logFClineage1 < -1) %>%
#   arrange(desc(waldStat_lineage1)) %>%
#   head(5) %>% rownames()
# 
# # 3. 合并进行可视化
# genes_to_plot <- c(top_up_L1, top_down_L1)


#取celltpye和配色信息
coldata <- data.frame(celltype =sim@colData$cell_type,group=sim@colData$group)
rownames(coldata)=colnames(sim)
#把sce中的3000个细胞对应信息取出
filter_coldata<-coldata[colnames(sce),]

# 提取 Lineage 1 的权重
l1_weights <- cell_Weight[, "Lineage1"]
# 定义一个阈值（通常为 0.5，或者直接保留权重 > 0 的细胞）
# 只有权重较高的细胞才被认为属于该路径
is_l1_cell <- l1_weights > 0.5 
#添加拟时序信息
filter_coldata$Pseudotime=sce$crv$pseudotime.Lineage1
#tagetgene
top_exp=sce@assays@data$counts[target_genes_L1_all,]
#top6_exp= log2(top6_exp+1) %>% t()
#获得最终的清洁数据
plt_data=cbind(filter_coldata, t(top_exp))
colnames(plt_data)

plt_data_l1 <- plt_data[is_l1_cell, ]

# 3. 检查筛选后的数据
message("原细胞数: ", nrow(plt_data))
message("Lineage 1 相关细胞数: ", nrow(plt_data_l1))


library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
#getPalette = colorRampPalette(brewer.pal(12, "Paired"))
#mycolors = getPalette(length(unique(plt_data$celltype)))

# 为了拼图美观，我们把legend隐藏掉
plt_list = list()
for (gene in target_genes_L1_all) {
  print(gene)
  p = ggscatter(data = plt_data_l1,
                x = 'Pseudotime',
                y = gene,
                color = 'celltype',
                size = 0.1)+
    geom_smooth(se = F, color = 'black')+
    theme_bw()+
    scale_color_manual(values = cell_type_colors)+
    theme(legend.position = 'none')
  plt_list[[gene]] = p
}

library(patchwork)
pp <- wrap_plots(plt_list)
pp
ggsave(plot = pp,filename = './slingshot/03_topgenes_Lineage1.png', width = 9, height = 5)
ggsave(plot = pp,filename = './slingshot/03_topgenes_Lineage1.pdf', width = 9, height = 5)



################################################################################
lineage_names <- paste0("Lineage", 1:n_lineages)

# 基因筛选阈值
logFC_cutoff <- 1
pval_cutoff  <- 0.05
top_n_genes  <- 20  # 每个谱系取前20个基因

# 细胞权重阈值（判断属于该谱系）
weight_cutoff <- 0.5
# ========================================================

# 批量循环：每个 Lineage 自动分析 + 绘图
for (i in 1:n_lineages) {
  lin_name <- paste0("Lineage", i)
  message("正在处理：", lin_name)
  
  # ===================== 1. 筛选当前谱系显著基因 =====================
  target_genes <- combined_res %>%
    filter(pvalue_AT < pval_cutoff & 
             get(paste0("pvalue_", tolower(lin_name))) < pval_cutoff) %>%
    filter(abs(get(paste0("logFC", tolower(lin_name)))) > logFC_cutoff) %>%
    arrange(desc(get(paste0("waldStat_", tolower(lin_name))))) %>%
    head(top_n_genes) %>%
    rownames()
  
  # ===================== 2. 准备绘图元数据 =====================
  coldata <- data.frame(
    celltype = sim@colData$cell_type,
    group    = sim@colData$group,
    row.names = colnames(sim)
  )
  filter_coldata <- coldata[colnames(sce), ]
  
  # 添加当前谱系伪时间
  pseudo_col <- paste0("pseudotime.", lin_name)
  filter_coldata$Pseudotime <- sce$crv[[pseudo_col]]
  
  # 提取基因表达矩阵
  top_exp <- sce@assays@data$counts[target_genes, ]
  
  # 合并（关键：转置表达矩阵！）
  plt_data <- cbind(filter_coldata, t(top_exp))
  
  # ===================== 3. 筛选属于当前谱系的细胞 =====================
  l1_weights <- cell_Weight[, lin_name]
  is_lin_cell <- l1_weights > weight_cutoff
  plt_data_sub <- plt_data[is_lin_cell, ]
  
  message(lin_name, " 筛选后细胞：", nrow(plt_data_sub))
  
  # ===================== 4. 批量绘制基因趋势图 =====================
  plt_list <- list()
  for (gene in target_genes) {
    p <- ggscatter(
      data = plt_data_sub,
      x = "Pseudotime",
      y = gene,
      color = "celltype",
      size = 0.1
    ) +
      geom_smooth(se = FALSE, color = "black", linewidth = 0.8) +
      theme_bw() +
      scale_color_manual(values = cell_type_colors) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(gene)
    
    plt_list[[gene]] <- p
  }
  
  # ===================== 5. 拼图 + 保存 =====================
  pp <- wrap_plots(plt_list)
  pp_final <- pp +
    plot_layout(guides = "collect", ncol = 4) &
    theme(legend.position = "right")
  ggsave(plot = pp_final, 
         filename = paste0("./slingshot/03_topgenes_", lin_name, "_celltype.pdf"), 
         width = 10, height = 10)
  ggsave(plot = pp_final, 
         filename = paste0("./slingshot/03_topgenes_", lin_name, "_celltype.png"), 
         width = 10, height = 10, dpi = 300)
  
  message(lin_name, " 绘图完成！\n")
}



# 批量循环：每个 Lineage 自动分析 + 绘图
for (i in 1:n_lineages) {
  lin_name <- paste0("Lineage", i)
  message("正在处理：", lin_name)
  
  # ===================== 1. 筛选当前谱系显著基因 =====================
  target_genes <- combined_res %>%
    filter(pvalue_AT < pval_cutoff & 
             get(paste0("pvalue_", tolower(lin_name))) < pval_cutoff) %>%
    filter(abs(get(paste0("logFC", tolower(lin_name)))) > logFC_cutoff) %>%
    arrange(desc(get(paste0("waldStat_", tolower(lin_name))))) %>%
    head(top_n_genes) %>%
    rownames()
  
  # ===================== 2. 准备绘图元数据 =====================
  coldata <- data.frame(
    celltype = sim@colData$cell_type,
    group    = sim@colData$group,
    row.names = colnames(sim)
  )
  filter_coldata <- coldata[colnames(sce), ]
  
  # 添加当前谱系伪时间
  pseudo_col <- paste0("pseudotime.", lin_name)
  filter_coldata$Pseudotime <- sce$crv[[pseudo_col]]
  
  # 提取基因表达矩阵
  top_exp <- sce@assays@data$counts[target_genes, ]
  
  # 合并（关键：转置表达矩阵！）
  plt_data <- cbind(filter_coldata, t(top_exp))
  
  # ===================== 3. 筛选属于当前谱系的细胞 =====================
  l1_weights <- cell_Weight[, lin_name]
  is_lin_cell <- l1_weights > weight_cutoff
  plt_data_sub <- plt_data[is_lin_cell, ]
  
  message(lin_name, " 筛选后细胞：", nrow(plt_data_sub))
  
  # ===================== 4. 批量绘制基因趋势图 =====================
  plt_list <- list()
  for (gene in target_genes) {
    p <- ggscatter(
      data = plt_data_sub,
      x = "Pseudotime",
      y = gene,
      color = "group",
      size = 0.1
    ) +
      geom_smooth(se = FALSE, color = "black", linewidth = 0.8) +
      theme_bw() +
      scale_color_manual(values = c("#E41A1C", "#377EB8")) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(gene)
    
    plt_list[[gene]] <- p
  }
  
  # ===================== 5. 拼图 + 保存 =====================
  pp <- wrap_plots(plt_list)
  pp_final <- pp +
    plot_layout(guides = "collect", ncol = 4) &
    theme(legend.position = "right")
  #pp_final
  
  ggsave(plot = pp_final, 
         filename = paste0("./slingshot/03_topgenes_", lin_name, "_group.pdf"), 
         width = 10, height = 10)
  ggsave(plot = pp_final, 
         filename = paste0("./slingshot/03_topgenes_", lin_name, "_group.png"), 
         width = 10, height = 10, dpi = 300)
  
  message(lin_name, " 绘图完成！\n")
}


#################################################################################

save.image("./slingshotDT.Rdata")

