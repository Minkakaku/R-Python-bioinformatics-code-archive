rm(list = ls())
gc()

library(Seurat)
library(slingshot)
library(SingleCellExperiment)
#library(qs)
library(tidyverse)
library(RColorBrewer)
library(tradeSeq)

result_path <- "/mnt/NexData_file/home/xzx/project/SC2411133007-07/shouhou5"
dir.create(result_path)
setwd(result_path)

sc <- readRDS("../shouhou2/adata_epi_fin.rds")
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
colData(sim)$celltype = meta$cell_type
colData(sim)$group = meta$group


#plot(umap, col = rgb(0,0,0,.5), pch=16, asp = 1)
sim

sim <- slingshot(sim, 
                 clusterLabels = 'celltype',  # 选择colData中细胞注释的列名
                 reducedDim = 'UMAP',  
                 start.clus= "Normal_ECs",  # 选择起点
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
  "left",  legend = "Lineage1",      # 轨迹名称（多轨迹时需对应修改）
  col = unique(brewer.pal(6, "Set1")),
  inset = 0.8,                    # 图例位置（避免遮挡细胞）
  pch = 16
)


# 生成100段渐变颜色（Spectral配色，去掉中间的黄色避免刺眼）
colors <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)
# 将伪时间分段，匹配渐变颜色（slingPseudotime_1=第一条轨迹的伪时间）
plotcol <- colors[cut(sim$slingPseudotime_1, breaks = 100)]
# 伪时间为NA的细胞（不在主轨迹上）设为灰色
plotcol[is.na(plotcol)] <- "grey"
# 绘制伪时间UMAP图
plot(
  reducedDims(sim)$UMAP,
  col = plotcol,    # 颜色=伪时间（渐变）
  #pch = 12,       # 点形状（15=正方形，更清晰）
  asp = 1
)
# 叠加轨迹线
lines(SlingshotDataSet(sim), lwd = 2, col = brewer.pal(9, "Set1"))
# 图例
legend(
  "left",
  legend = "Lineage1 (Pseudotime)",
  col = unique(brewer.pal(6, "Set1")),
  inset = 0.8,
  pch = 16
)

# 1. 提取伪时间，保留第一条轨迹（多轨迹需调整）
sling_Pseudotime <- as.data.frame(slingPseudotime(sim)) #%>% select(-last_col())
colnames(sling_Pseudotime) <- "pseudotime"  # 重命名为"pseudotime"
# 2. 将伪时间添加到Seurat的meta.data中
sc <- AddMetaData(sc, metadata = sling_Pseudotime)
# 3. 可视化整合结果
# （1）细胞类型UMAP图
color <- c(pal_npg()(10), pal_jco()(10), pal_lancet()(9), pal_aaas()(10))  # 多配色组合
DimPlot(sc, group.by = "cell_type")
# （2）伪时间分布UMAP图（蓝色=起始，红色=终止）
FeaturePlot(sc, features = "pseudotime", cols = c("blue", "red"))
# （3）不同细胞类型的伪时间箱线图
VlnPlot(
  sc,
  features = "pseudotime",
  group.by = "cell_type",
  pt.size = 0#,  # 隐藏单个点，只显示箱线
  #cols = color
) + theme(axis.text.x = element_text(angle = 40, hjust = 1))  # 旋转x轴标签
# （4）不同处理组（condition）的伪时间比较（如对照组vs处理组）
VlnPlot(
  sc,
  features = "pseudotime",
  group.by = "group",
  pt.size = 0#,
  #cols = color
) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
# （5）终极整合图：UMAP+细胞类型+伪时间（点大小=伪时间）
data <- cbind(
  sc@meta.data,
  as.data.frame(sc@reductions$umap@cell.embeddings)[rownames(sling_Pseudotime), ]
)
ggplot(data, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    aes(
      x = UMAP_1,
      y = UMAP_2,
      size = pseudotime,    # 点大小=伪时间（越大越接近终止）
      colour = cell_type   # 颜色=细胞类型
    ),
    data = data
  ) +
  #scale_color_manual(values = color) +
  scale_size_continuous(range = c(0.3, 1.5)) +  # 控制点大小范围
  theme_classic() +
  theme(
    axis.line = element_line(
      size = 0.6,
      colour = "black",
      arrow = arrow(length = unit(0.1, "inches"), ends = "last", type = "open")
    ) 
  )


########
# 提取计数矩阵（SCE对象的assay数据）
counts <- sim@assays@data$counts
# 转换为SlingshotDataSet对象（tradeSeq的标准输入）
SDS <- SlingshotDataSet(sim)
# 1. 确定GAM模型的最佳k值（k= knots数量，控制模型灵活性）
set.seed(123)  # 设随机种子，确保结果可重复
icMat <- evaluateK(
  counts = counts,
  sds = SDS,
  k = 5:12,          # 测试k=3、4、5（k越大越灵活，易过拟合）
  nGenes = 500,        # 【演示用】仅用15个基因加速，实际需设为500+
  verbose = T         # 显示计算过程
)
# 后续用icMat选择最小AIC的k值（如k=4）
# 2. 获取伪时间和细胞权重
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
sET <- startVsEndTest(sce_converged)  # 比较“起始细胞”与“终止细胞”的基因表达
# 关键指标：waldStat（值越大，起止点差异越显著）、logFClineage1（起始→终止的log2倍变化）
# 3. 筛选top差异基因（以waldStat排序）
# order_sET <- order(sET$waldStat, decreasing = TRUE)
# topGeneStart <- names(sce)[order_sET[1]]  # 差异最显著的第一个基因（如示例中的"FTL"）


# 1. 提取细胞拟时间（从sce_converged中获取，按细胞排序）
#cell_order <- order(pseudotime)  # 按拟时间对细胞排序（关键：保证热图中细胞按时间从左到右排列）
pseudotime_sorted <- as.data.frame(pseudotime)
cell_order <- order(pseudotime_sorted$Lineage1)
pseudotime_sorted <- pseudotime_sorted[cell_order]  # 排序后的拟时间

# 2. 提取显著基因（可选：用之前的associationTest结果，筛选padj<0.05的基因，减少噪声）
#sig_genes <- rownames(Relation[Relation$padj < 0.05, ])  # 显著关联基因
# 或用startVsEndTest的显著差异基因
sig_genes <- rownames(sET[sET$pvalue < 0.05, ])

# 3. 提取基因表达矩阵（标准化后的表达量，按「基因×排序后细胞」整理）
expr_mat <- assay(sce_converged, "counts")  # 提取log标准化表达矩阵（tradeSeq拟合后默认存在）
expr_mat_sig <- expr_mat[sig_genes, cell_order]  # 仅保留显著基因+排序后细胞

save.image("./heatmap.Rdata")
