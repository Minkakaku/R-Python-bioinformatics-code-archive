library(tidyverse)
library(Seurat)
getwd()
library(patchwork)
options(java.parameters = '-Xmx10g')  # options('java.parameters')
library(rJava)
library(RBioFormats)
library(ggpubr)
library(SeuratDisk)
library(ggrepel)
library(loupeR)

result_dir <- 'result' 
dir.create(result_dir)

rctd_result_ls <- readRDS("rctd_result_ls.rds")
xenium_rctd <- readRDS("xenium_rctd.rds")

xenium_rctd$first_type <- 
  xenium_rctd$first_type %>% 
  str_replace_na("Others")


plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", 
          "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", 
          "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879", 
          "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", 
          "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", "#ffc300", 
          "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

xenium@meta.data$first_type %>% unique()

###把Keratinocytes和Malignant cells合并成Epithelial这一种类别----
# 创建新列，初始值与 first_type 相同
xenium@meta.data$first_type_coEpi <- xenium@meta.data$first_type
# 将 Malignant cells 和 Keratinocytes 重命名为 Epithelial
xenium@meta.data$first_type_coEpi[xenium@meta.data$first_type_coEpi %in% c("Malignant cells", "Keratinocytes")] <- "Epithelial"
# 验证合并结果
table(xenium@meta.data$first_type_coEpi,xenium@meta.data$first_type)


#1、UMAP----
# 按免疫属性重新排序（免疫细胞在前，非免疫细胞在后）
cell_types_order <- c(
  # 免疫细胞
  "Natural killer cells",  # 先天免疫核心[6](@ref)
  "CD8+ T cells",          # 细胞毒性T细胞
  "CD4+ T cells",          # 辅助T细胞
  "Tregs",                 # 调节性T细胞
  "B cells",               # B淋巴细胞
  "Plasma cells",          # B细胞终末分化态[8](@ref)
  "DC",                    # 树突状细胞（抗原提呈）
  "pDC",                   # 浆细胞样树突细胞
  "Mono-Mac",              # 单核/巨噬细胞（吞噬功能）[5](@ref)
  "Mast cells",            # 肥大细胞（过敏反应）
  "LC",                    # 朗格汉斯细胞（皮肤免疫）
  
  # 非免疫
  "Epithelial",
  #"Malignant cells",       # 肿瘤细胞
  "Proliferating cells",   # 增殖态细胞（可含免疫细胞）
  "Endothelial cells",     # 血管内皮细胞[8](@ref)
  "Fibroblasts",           # 成纤维细胞（结缔组织）
  #"Keratinocytes",         # 角质细胞（表皮）
  "Melanocytes",           # 黑色素细胞
  "Others"                 # 未定义细胞
)

# 将元数据中的细胞类型转换为有序因子
xenium@meta.data$first_type_coEpi <- factor(
  xenium@meta.data$first_type_coEpi,
  levels = cell_types_order
)

color_mapping <- setNames(
  c(
    "#FF6D31", "#FFB85F", "#FFCB18", "#FBF579", "#F0D879", "#FEEE7D", 
    "#ED5485", "#E8A0B8", "#D9534F", "#EF5285", "#F17D80",
    "#7F95D1", "#02C9C9", "#4695D6", "#73B66B", 
    #"#8BC24C", 
    "#600473", "#999999"
  ),
  cell_types_order
)

# 1. 不分组的UMAP
pdf(file.path(result_dir, 'RCTD_umap1.pdf'), height = 12, width = 14)
DimPlot(xenium, 
        group.by = 'first_type_coEpi', 
        label = FALSE, 
        cols = color_mapping)  # 应用指定颜色 [6,7](@ref)
dev.off()

# 2. 按group.y分组的UMAP
pdf(file.path(result_dir, 'RCTD_group_umap1.pdf'), height = 10, width = 18)
DimPlot(xenium, 
        group.by = 'first_type_coEpi', 
        label = FALSE, 
        split.by = 'group1',
        cols = color_mapping)  # 应用相同颜色方案 [1,5](@ref)
dev.off()

# pdf(file.path(result_dir, 'RCTD_umap.pdf'), height = 12, width = 14)
# DimPlot(xenium_rctd, group.by = 'first_type', label = F, cols = plan) 
# dev.off()
# 
# pdf(file.path(result_dir, 'RCTD_group_umap.pdf'), height = 10, width = 18)
# DimPlot(xenium_rctd, group.by = 'first_type', label = F, split.by = 'group.y',cols = plan)
# dev.off()

#寻找标记基因----

res <- "first_type_coEpi"
# Idents(All_T) <- "RNA_snn_res.0.6"
markers <-
  FindAllMarkers(xenium,
    # JoinLayers(All_T),
    group.by = res,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

top5_markers <- 
  markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 5) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(top5_markers) <- sub('.*\\.','',names(top5_markers))


#提取Top10,deepseek
top10_markers <-
  markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)


# dotplot
pdf(file.path(result_dir, str_glue("dotplot top5_marker.pdf")), width = 35, height = 6)
DotPlot(xenium, group.by = res, features = top5_markers)+RotatedAxis()
dev.off()

p <- DotPlot(
  xenium,
  features = top5_markers,  # 使用去重后的基因
  group.by = 'first_type_coEpi',
  cols = c("#ffeead", "red"), # 颜色梯度：低表达→高表达
  dot.scale = 6,                # 点大小缩放因子（默认4）
  cluster.idents = FALSE        # 禁止自动重排聚类顺序
) +
  # 移除顶部分组标签的关键设置
  #facet_null() +
  theme(# 旋转X轴标签90°，并调整垂直/水平对齐（vjust/hjust）
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        # 关键：隐藏顶部的分面标签（这里对应 group.by 的分组名称显示）
        strip.text.x = element_blank(),  
        strip.background = element_blank()
        )+ 
  # 旋转X轴标签90°，并调整垂直/水平对齐（vjust/hjust）
  labs(x = "", y = "Cell")


ggsave("result/Cell_markers2.pdf", plot = p, width = 20, height = 6, dpi = 300)


#2、cell_gene_Dotplot----

#xenium_rctd的原始数据叫count，默认是counts；重新生成了一下xenium对象，然后重新norm的数据
xenium_rctd2 <- JoinLayers(xenium_rctd,layers = 'count')

xenium_rctd2 <- CreateSeuratObject(
  LayerData(xenium_rctd, layer = "count"),
  meta.data = xenium_rctd[[]]
)

#空值替换
xenium_rctd2 <- NormalizeData(xenium_rctd2)
xenium_rctd2$first_type <- 
  xenium_rctd2$first_type %>% 
  str_replace_na("Others")

#绘图基因list
# gene_list <- c(
#   "EPCAM", "BCAM",
#   "MLANA", 
#   "CD3E", "CD2", "CD4", 
#   "IL2RA", "FOXP3", 
#   "CD8A", "GZMA", 
#   "KLRC1", "PRF1", 
#   "MKI67", 
#   "MS4A1", "CD19", 
#   "CD38", 
#   "IRF8", "CLEC4C", 
#   "LAMP3", "CCR7", 
#   "CD1A", "CD207", 
#   "CD68", "CD14", 
#   "KIT", "SOCS1"
# )

gene_list <- c(
  "KRT8",    # 角蛋白8，腺上皮骨架核心成分
  "KRT18",   # 角蛋白18，胃肠道/肝脏上皮特征
  "KRT19",   # 角蛋白19，乳腺癌转移灶标记
  "EPCAM",   # 上皮细胞黏附分子，泛上皮标志
  #"BCAM",
  "CLDN7",   # 紧密连接蛋白7，维持上皮屏障
  "MUC1",    # 黏蛋白1，促肿瘤侵袭因子
  "MUC2",    # 黏蛋白2，肠化生组织特征
  "KRT7",    # 角蛋白7，鉴别腺癌来源
  "KRT20",   # 角蛋白20，结直肠癌特异性标记
  "KRT10",    # 角蛋白10，鳞状上皮分化标记
  "PDGFRA",  # CAF核心标记，促纤维化
  "DPEP1",   # 基质重塑酶，胰腺癌CAF高表达
  "COL1A1",  # I型胶原α1链，ECM沉积标志
  "COL3A1",  # III型胶原α1链，肿瘤间质硬度关联
  "FAP",     # 成纤维细胞活化蛋白，免疫抑制性CAF
  "ACTA2",   # α-平滑肌肌动蛋白，肌成纤维细胞特征
  "THY1",    # 间充质干细胞来源CAF标记
  "DCN",     # 核心蛋白多糖，维持基质弹性
  "LUM",     # 光蛋白聚糖，胶原组装调控
  "FN1",      # 纤连蛋白1，促进癌细胞迁移
  "PECAM1",  # CD31，泛内皮标记
  "CDH5",    # VE-钙黏蛋白，内皮屏障功能
  "VWF",     # 血管性血友病因子，静脉内皮特征
  "KDR",     # VEGFR2，血管生成受体
  "MLANA","TYR","MITF","TYPR1","DCT",
  "CD3E", "CD2", "CD4",
  "IL2RA", "FOXP3",
  "CD8A", "GZMA",
  "KLRC1", "PRF1",
  "MKI67",
  "MS4A1", "CD19",
  "CD38",
  "IRF8", "CLEC4C",
  "LAMP3", "CCR7",
  "CD1A", "CD207",
  "CD68", "CD14",
  "KIT", "SOCS1"
)

valid_genes <- gene_list[gene_list %in% rownames(xenium)]
cat("有效基因:", length(valid_genes), "/", length(gene_list))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_genes <- setdiff(gene_list, valid_genes)

# xenium_rctd <- JoinLayers(xenium_rctd,layers = 'count')
# DefaultAssay(xenium_rctd) <- "RNA"
# 
# undebug(NormalizeData)
# xenium_rctd <- NormalizeData(xenium_rctd,
#                              scale.factor = median(xenium_rctd$nCount_RNA))
# print(Layers(xenium_rctd))

# # 高变基因
# xenium_rctd2 <- FindVariableFeatures(xenium_rctd2)
# #缩放
# xenium_rctd2 <- ScaleData(xenium_rctd2, features = rownames(xenium_rctd2))
# # 运行PCA
# xenium_rctd2 <- RunPCA(xenium_rctd2)
# # 运行UMAP
# xenium_rctd2 <- RunUMAP(xenium_rctd2, dims = 1:20)  # dims使用PCA的主成分数
# # 运行t-SNE
# xenium_rctd2 <- RunTSNE(xenium_rctd2, dims = 1:20)


cell_order <- c(
  # "Malignant cells",
  # "Keratinocytes",
  "Epithelial",
  "Fibroblasts",
  "Endothelial cells",
  "Melanocytes",
  "CD4+ T cells",
  "Tregs",
  "CD8+ T cells",
  "Natural killer cells",
  "Proliferating cells",
  "B cells",
  "Plasma cells",
  "pDC",
  "DC",
  "LC",
  "Mono-Mac",
  "Mast cells",
  "Others"
)

#将分组列转换为有序因子
xenium@meta.data$first_type_coEpi <- factor(
  xenium@meta.data$first_type_coEpi,
  levels = cell_order  # 强制锁定顺序
)

p <- DotPlot(
  xenium,
  features = valid_genes,  # 使用去重后的基因
  group.by = 'first_type_coEpi',
  cols = c("#ffeead", "red"), # 颜色梯度：低表达→高表达
  dot.scale = 6,                # 点大小缩放因子（默认4）
  cluster.idents = FALSE        # 禁止自动重排聚类顺序
) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+ 
  # 旋转X轴标签90°，并调整垂直/水平对齐（vjust/hjust）
  labs(x = "", y = "Cell Clusters")

ggsave("result/Cell_markers1.pdf", plot = p, width = 10, height = 6, dpi = 300)




###热图----
valid_gene_list <- gene_list[gene_list %in% rownames(xenium)]
cat("有效基因:", length(valid_gene_list), "/", length(gene_list))
# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(gene_list, valid_gene_list)
gene <- valid_gene_list

# 假设 immune_genes 是已分类的命名向量（如 'T_func' = "CD8A", ...）
ordered_gene <- names(gene)  # 直接提取基因名（忽略分类名）
# 或按分类组排序（保持组内顺序）：
ordered_gene <- unname(gene)  # 仅保留基因名，顺序与定义一致


#xenium <- ScaleData(xenium, vars.to.regress = c('nCount_RNA'))

H <- DoHeatmap(
  object = xenium, 
  features = ordered_gene,   # 关键：传入有序基因向量
  slot = "scale.data",    
  group.by = "first_type_coEpi", 
  group.colors = color_mapping, 
  disp.min = -2.5, 
  disp.max = 2.5,
  angle = -45, 
  hjust = 0.8,
  size = 3,
  raster = TRUE               # 加速渲染（大数据集必备）[1,6](@ref)
)+
  scale_fill_gradientn(
    colors = c("navy", "white", "firebrick3"),
    values = scales::rescale(c(-2.5, 0, 2.5))  # 精确匹配disp.min/max[4](@ref)
  ) 

ggsave("result/MarkerHeatmap.pdf",plot = H,width = 12,height = 12,dpi = 300)



H1 <- DoHeatmap(
  object = xenium, 
  features = top5_markers,   # 关键：传入有序基因向量
  slot = "scale.data",    
  group.by = "first_type_coEpi", 
  group.colors = color_mapping, 
  disp.min = -2.5, 
  disp.max = 2.5,
  angle = -45, 
  hjust = 0.8,
  size = 3,
  raster = TRUE               # 加速渲染（大数据集必备）[1,6](@ref)
)+
  scale_fill_gradientn(
    colors = c("navy", "white", "firebrick3"),
    values = scales::rescale(c(-2.5, 0, 2.5))  # 精确匹配disp.min/max[4](@ref)
  ) 

ggsave("result/MarkerHeatmap1.png",plot = H1,width = 12,height = 12,dpi = 300)




#3、细胞比例----

#3.1按groupx分组----
# 指定样本顺序
groupx_order <- c("FYN", "WSL", "GYX", "CYS", "ZKJ", "ZBH", "SQQ", "HYY", "CZY","ZL")
# 细胞比例统计图----
celltype_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  mutate(
    # 将group.x转换为因子并指定顺序
    group.x = factor(group.x, levels = groupx_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(first_type, group.x) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(group.x) %>%
  mutate(prop = n / sum(n))

# 基础 ggplot2 堆叠柱状图
p1 <- ggplot(celltype_counts, aes(x = .data[['group.x']], y = prop, fill = .data[['first_type']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = color_mapping) +
  # 坐标轴与图例优化
  labs(x = "Sample", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(
  filename = "result/1.1Celltype_proportions.pdf",  # 保存文件名
  plot = p1,                                # 指定要保存的图形对象
  width = 10,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
)


# 细胞样本比例统计图----
# 按样本和细胞类型分组计数
celltype_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  mutate(
    # 将group.x转换为因子并指定顺序
    group.x = factor(group.x, levels = groupx_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(first_type, group.x) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(first_type) %>%
  mutate(prop = n / sum(n))

# 基础 ggplot2 堆叠柱状图
p2 <- ggplot(celltype_counts, aes(x = .data[['first_type']], y = prop, fill = .data[['group.x']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = plan) +
  # 坐标轴与图例优化
  labs(x = "Sample", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_minimal() +
  coord_flip() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2


ggsave(
  filename = "result/1.2CelltypeSample_proportions.pdf",  # 保存文件名
  plot = p2,                                # 指定要保存的图形对象
  width = 10,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
)

#细胞数绝对值的柱状图----
# 准备配色方案
color_mapping <- setNames(
  c("#FF6D31", "#FFB85F", "#FFCB18", "#FBF579", "#F0D879", "#FEEE7D", 
    "#ED5485", "#E8A0B8", "#D9534F", "#EF5285", "#F17D80",
    "#7F95D1", "#02C9C9", "#4695D6", "#73B66B", "#8BC24C", "#600473", "#999999"),
  cell_types_order
)

# 计算各样本的细胞类型数量
celltype_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%
  mutate(
    group.x = factor(group.x, levels = groupx_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  count(group.x, first_type, name = "cell_count")

# 计算柱体位置
bar_width <- 0.7  # 总柱宽
n_types <- length(cell_types_order)
offset <- seq(-bar_width/2, bar_width/2, length.out = n_types + 2)[2:(n_types+1)]

p_group <- ggplot(celltype_counts, aes(x = group.x)) +
  # 1. 使用geom_col替代手动柱体绘制
  geom_col(
    aes(y = cell_count, fill = first_type),
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  # # 2. 修复标签层 - 关键修改点
  # geom_text(
  #   aes(
  #     y = cell_count + max(cell_count)*0.05,  # 设置y坐标
  #     label = ifelse(cell_count > 50, cell_count, ""),
  #     group = first_type
  #   ),
  #   position = position_dodge(0.8),  # 与柱体相同的dodge宽度
  #   vjust = 0,  # 垂直对齐方式
  #   size = 3.5
  # ) +
  # 配色方案保持不变
  scale_fill_manual(values = color_mapping)

ggsave("result/1.3CelltypeCounts_by_sample.pdf", p_group, 
       width = 14, height = 8, dpi = 300)


#3.2按照有效无效分组----
# 指定样本顺序
group1_order <- c("Effective", "Ineffective")

# 细胞比例统计图----
celltype_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  mutate(
    # 将group.x转换为因子并指定顺序
    group1 = factor(group1, levels = group1_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(first_type, group1) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(group1) %>%
  mutate(prop = n / sum(n))

# 基础 ggplot2 堆叠柱状图
p1 <- ggplot(celltype_counts, aes(x = .data[['group1']], y = prop, fill = .data[['first_type']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = color_mapping) +
  # 坐标轴与图例优化
  labs(x = "group", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(
  filename = "result/2.1Group_Celltype_proportions.pdf",  # 保存文件名
  plot = p1,                                # 指定要保存的图形对象
  width = 10,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
)



# 细胞样本比例统计图----
# 按样本和细胞类型分组计数
celltype_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  mutate(
    # 将group.x转换为因子并指定顺序
    group1 = factor(group1, levels = group1_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(first_type, group1) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(first_type) %>%
  mutate(prop = n / sum(n))

# 基础 ggplot2 堆叠柱状图
p2 <- ggplot(celltype_counts, aes(x = .data[['first_type']], y = prop, fill = .data[['group1']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = plan) +
  # 坐标轴与图例优化
  labs(x = "group", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_minimal() +
  coord_flip() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2


ggsave(
  filename = "result/2.2Celltypegroup_proportions.pdf",  # 保存文件名
  plot = p2,                                # 指定要保存的图形对象
  width = 10,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
)

###只画Melanocytes----
celltype_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  mutate(
    # 将group.x转换为因子并指定顺序
    group1 = factor(group1, levels = group1_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(first_type, group1) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(first_type) %>%
  mutate(prop = n / sum(n))

Melanocytes_counts <- celltype_counts %>% 
  filter(first_type == "Melanocytes")


# 基础 ggplot2 堆叠柱状图
p4 <- ggplot(Melanocytes_counts , aes(x = .data[['first_type']], y = prop, fill = .data[['group1']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = plan) +
  # 坐标轴与图例优化
  labs(x = "group", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_minimal() +
  #coord_flip() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p4


ggsave(
  filename = "result/Melanocytesgroup_proportions.pdf",  # 保存文件名
  plot = p4,                                # 指定要保存的图形对象
  width = 5,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
)


#箱线图----
###整体箱线图----
# 数据准备：计算每个样本的细胞比例
celltype_sample_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%
  mutate(
    group1 = factor(group1, levels = group1_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(sample, group1, first_type) %>% 
  tally() %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n))

# 分面箱线图（带显著性检验）
library(ggpubr)
ggplot(celltype_sample_counts, 
       aes(x = group1, y = prop, fill = group1)) +
  geom_boxplot(
    width = 0.7, 
    outlier.shape = 19,
    alpha = 0.8
  ) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +  # 显示样本点
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",        # 改为星号显示（*:p<0.05; **:p<0.01; ***:p<0.001）
    label.x = 1.5,             # 水平居中位置（两组中间）
    size = 6,                  # 增大星号字号（更醒目）
    vjust = -0.5               # 垂直微调避免重叠
  ) +
  scale_fill_manual(values = c("#4DBBD5", "#E64B35")) +  # 响应组配色
  facet_wrap(~ first_type, scales = "free_y", ncol = 3) +  # 关键：按细胞类型分面
  labs(x = "Group", y = "Proportion per Sample") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 10)  # 分面标题加粗
  )
ggsave("result/boxplot/CellType_Boxplots.pdf", width = 12, height = 20)


###分细胞类型箱线图----
box_dir <- 'result/boxplot'
dir.create(box_dir)

library(reshape2)
celltype_wide <- dcast(celltype_sample_counts, 
                       sample + group1 ~ first_type, 
                       value.var = "prop")

#画这种循环的需要把细胞名里的+和空格改掉
# 创建安全的列名映射（特殊字符转下划线）
safe_names <- gsub("[^[:alnum:]]", "_", cell_types_order) %>%
  make.unique()  # 避免重名冲突

# 重命名数据框中的列
colnames(celltype_wide) <- ifelse(
  colnames(celltype_wide) %in% cell_types_order,
  safe_names[match(colnames(celltype_wide), cell_types_order)],
  colnames(celltype_wide)  # 保留非细胞类型列
)

# 更新细胞类型顺序向量
safe_cell_types_order <- safe_names

library(ggpubr)

for (i in seq_along(cell_types_order)) {
  orig_name <- cell_types_order[i]       # 原始名称（用于标题）
  safe_name <- safe_cell_types_order[i]  # 安全名称（用于数据引用）
  
  if (safe_name %in% colnames(celltype_wide)) {
    p <- ggboxplot(
      celltype_wide, 
      x = "group1", 
      y = safe_name,  # 使用安全列名
      fill = "group1", 
      palette = c("#4DBBD5", "#E64B35"),
      add = "jitter"
    ) +
      stat_compare_means(
        method = "wilcox.test",
        label = "p.signif", # 改为星号显示（*:p<0.05; **:p<0.01; ***:p<0.001
        label.y = max(celltype_wide[[safe_name]], na.rm = TRUE) * 1.1
      ) +
      labs(
        title = paste("Proportion of", orig_name),  # 标题用原始名称
        x = "", 
        y = "Proportion"
      ) +
      theme(legend.position = "none")
    
    # 文件名也使用安全名称
    ggsave(
      filename = file.path(box_dir, paste0(safe_name, "_Boxplot.pdf")),
      plot = p,
      width = 4, 
      height = 5
    )
  }
}


#细胞数绝对值的柱状图----
# 准备配色方案
color_mapping <- setNames(
  c("#FF6D31", "#FFB85F", "#FFCB18", "#FBF579", "#F0D879", "#FEEE7D", 
    "#ED5485", "#E8A0B8", "#D9534F", "#EF5285", "#F17D80",
    "#7F95D1", "#02C9C9", "#4695D6", "#73B66B", "#8BC24C", "#600473", "#999999"),
  cell_types_order
)

# 计算各样本的细胞类型数量
celltype_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%
  mutate(
    group1 = factor(group1, levels = group1_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  count(group1, first_type, name = "cell_count")

# 计算柱体位置
bar_width <- 0.7  # 总柱宽
n_types <- length(cell_types_order)
offset <- seq(-bar_width/2, bar_width/2, length.out = n_types + 2)[2:(n_types+1)]

p_group <- ggplot(celltype_counts, aes(x = group1)) +
  # 1. 使用geom_col替代手动柱体绘制
  geom_col(
    aes(y = cell_count, fill = first_type),
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  # # 2. 修复标签层 - 关键修改点
  # geom_text(
  #   aes(
  #     y = cell_count + max(cell_count)*0.05,  # 设置y坐标
  #     label = ifelse(cell_count > 50, cell_count, ""),
  #     group = first_type
  #   ),
  #   position = position_dodge(0.8),  # 与柱体相同的dodge宽度
  #   vjust = 0,  # 垂直对齐方式
  #   size = 3.5
  # ) +
# 配色方案保持不变
scale_fill_manual(values = color_mapping)

ggsave("result/CelltypeCounts_by_group1.pdf", p_group, 
       width = 14, height = 8, dpi = 300)


#3.3免疫，非免疫----
# 定义免疫细胞类型列表
immune_cell_types <- c(
  "Natural killer cells", 
  "CD8+ T cells", 
  "CD4+ T cells", 
  "Tregs", 
  "B cells", 
  "Plasma cells", 
  "DC", 
  "pDC", 
  "Mono-Mac", 
  "Mast cells", 
  "LC"
)

xenium@meta.data <- xenium@meta.data %>%
  mutate(cell_group = ifelse(
    first_type %in% immune_cell_types, 
    "Immune_cell", 
    "Non_immune_cell"
  ))

table(xenium@meta.data$cell_group, useNA = "ifany")

# 生成列联表：统计cell_group下各first_type的细胞数量
table(xenium@meta.data$cell_group, xenium@meta.data$first_type) %>% 
  as.data.frame() %>%
  filter(Freq > 0) %>%  # 过滤0值
  dplyr::rename(cell_group = Var1, first_type = Var2) %>%
  arrange(cell_group, desc(Freq))

# 指定样本顺序
groupx_order <- c("FYN", "WSL", "GYX", "CYS", "ZKJ", "ZBH", "SQQ", "HYY", "CZY","ZL")
# 细胞比例统计图----
cellgroup_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  mutate(
    # 将group.x转换为因子并指定顺序
    group.x = factor(group.x, levels = groupx_order)
  ) %>%
  group_by(cell_group, group.x) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(group.x) %>%
  mutate(prop = n / sum(n))

# 基础 ggplot2 堆叠柱状图
p1 <- ggplot(cellgroup_counts, aes(x = .data[['group.x']], y = prop, fill = .data[['cell_group']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = plan) +
  # 坐标轴与图例优化
  labs(x = "Sample", y = "Cell_group Proportion", fill = "cell_group") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(
  filename = "result/3.1Cellgroup_Sample_proportions.pdf",  # 保存文件名
  plot = p1,                                # 指定要保存的图形对象
  width = 10,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
)


# 指定样本顺序
groupx_order <- c("FYN", "WSL", "GYX", "CYS", "ZKJ", "ZBH", "SQQ", "HYY", "CZY","ZL")
# 细胞比例统计图(免疫细胞内)----
cellgroup_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  filter(cell_group != "Non_immune_cell") %>% # 过滤非免疫细胞
  mutate(
    # 将group.x转换为因子并指定顺序
    group.x = factor(group.x, levels = groupx_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(first_type, group.x) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(group.x) %>%
  mutate(prop = n / sum(n))

# 基础 ggplot2 堆叠柱状图
p1 <- ggplot(cellgroup_counts, aes(x = .data[['group.x']], y = prop, fill = .data[['first_type']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = color_mapping) +
  # 坐标轴与图例优化
  labs(x = "Sample", y = "Cell_group Proportion", fill = "Immune_cell") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(
  filename = "result/3.2Immune_cellgroup_proportions.pdf",  # 保存文件名
  plot = p1,                                # 指定要保存的图形对象
  width = 10,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
)

# 细胞比例统计图(非免疫细胞内)----
cellgroup_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  filter(cell_group != "Immune_cell") %>% # 过滤非免疫细胞
  mutate(
    # 将group.x转换为因子并指定顺序
    group.x = factor(group.x, levels = groupx_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(first_type, group.x) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(group.x) %>%
  mutate(prop = n / sum(n))

# 基础 ggplot2 堆叠柱状图
p1 <- ggplot(cellgroup_counts, aes(x = .data[['group.x']], y = prop, fill = .data[['first_type']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = color_mapping) +
  # 坐标轴与图例优化
  labs(x = "Sample", y = "Cell_group Proportion", fill = "Non_immune_cell") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(
  filename = "result/3.3Non_immune_cell_group_proportions.pdf",  # 保存文件名
  plot = p1,                                # 指定要保存的图形对象
  width = 10,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
)

#细胞数绝对值的柱状图(免疫细胞内)----
# 计算各样本的细胞类型数量
celltype_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%
  filter(cell_group != "Non_immune_cell") %>% # 过滤非免疫细胞
  mutate(
    group.x = factor(group.x, levels = groupx_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  count(group.x, first_type, name = "cell_count")

# 计算柱体位置
bar_width <- 0.7  # 总柱宽
n_types <- length(cell_types_order)
offset <- seq(-bar_width/2, bar_width/2, length.out = n_types + 2)[2:(n_types+1)]

p_group <- ggplot(celltype_counts, aes(x = group.x)) +
  # 1. 使用geom_col替代手动柱体绘制
  geom_col(
    aes(y = cell_count, fill = first_type),
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  # # 2. 修复标签层 - 关键修改点
  # geom_text(
  #   aes(
  #     y = cell_count + max(cell_count)*0.05,  # 设置y坐标
  #     label = ifelse(cell_count > 50, cell_count, ""),
  #     group = first_type
  #   ),
  #   position = position_dodge(0.8),  # 与柱体相同的dodge宽度
  #   vjust = 0,  # 垂直对齐方式
  #   size = 3.5
  # ) +
# 配色方案保持不变
scale_fill_manual(values = color_mapping)

ggsave("result/3.4ImmunecelltypeCounts_by_Groupx.pdf", p_group, 
       width = 14, height = 8, dpi = 300)



#4、免疫基因空间featureplot----
spatialfeature_dir <- 'result/spatialfeature' 
dir.create(spatialfeature_dir)


cells_boundary_ls<-readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls<-readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))


immune_genes <- c(
  "ARG1", "CCL2", "CCL5", "CD3D", "CD3E", "CD4", "CD68", "CD69", "CD74", "CD8A", 
  "CD163", "CD274", "CIITA", "CTLA4", "CXCL9", "CXCL10", "CXCL13", "FOXP3", 
  "GZMA", "GZMB", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "ICOS", 
  "IDO1", "IFNG", "IL1B", "IL6", "IL10", "IL12B", "LAG3", "NOS2", "PDCD1", 
  "PRF1", "TGFB1", "TIGIT", "TNF", "TNFRSF9", "TOX"
)
valid_immune_genes <- immune_genes[immune_genes %in% rownames(xenium_rctd2)]
cat("有效基因:", length(valid_immune_genes), "/", length(immune_genes))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(immune_genes, valid_immune_genes)

gene_cus <- valid_immune_genes

#no FOV空间featureplot----
xenium_rctd2_meta <- xenium_rctd2[[]] %>% rownames_to_column('cell_id')
dat_expr <- xenium_rctd2 %>% 
  FetchData(vars = unlist(gene_cus), layer = "data") %>% 
  rownames_to_column('cell_id')
stopifnot(all(xenium_rctd2_meta$cell_id == dat_expr$cell_id))
dat_expr_coords_full <- left_join(xenium_rctd2_meta, dat_expr, by = 'cell_id')

###每个样本所有基因一张图----
marker_p_ls <- 
  gene_cus %>%
  setNames(., .) %>% 
  map(function(x) {
    p_combined <-
      cells_boundary_ls %>%
      imap(function(cells_boundary, sample) { # browser()
        if (!x %in% colnames(dat_expr_coords_full)) return(NULL)
        
        dat <-
          cells_boundary %>%
          right_join(dat_expr_coords_full, by = "cell_id")      
        
        if(all(is.na(dat$vertex_x))) return(NULL) 
        if(any(is.na(dat$vertex_x))) warning("vertex_x have NA values ...", immediate. = TRUE) 
        
        p <- dat %>%
          arrange(.data[[x]]) %>% 
          as.data.frame() %>%
          ggplot(aes(x = vertex_x, y = -vertex_y, color = .data[[x]])) +
          scattermore::geom_scattermore(aes(alpha = .data[[x]]), pointsize = 0.5, show.legend = TRUE) +
          guides(alpha = "none") +
          theme_void() +
          coord_equal() +
          ggtitle(sample)
        
        # all_values <- unlist(lapply(p_combined$patches$plots, function(p) p$data[[p$labels$colour]]))
        p <- p + scale_color_gradientn(
          colours = c("blue", "white", "red"),
          # limits = c(0, ceiling(max(dat_expr[-1]))),  # 统一颜色范围
          limits = c(0, 1), 
          name = "Expression"
        )
        
        p
      }) %>%
      purrr::discard(is.null) %>% 
      patchwork::wrap_plots(ncol = 3, guides = "collect") &
      theme(legend.position = "right")  # 全局设置图例位置
    
    p_combined
  })


marker_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(spatialfeature_dir,  
                  str_glue("spatial feature plot {grp}_sameLegend_red2blue.pdf")), 
        height = 24, 
        width = 20)
    print(p)
    dev.off()
  })    


###每个样本每个基因一张图----
# 创建输出目录（如果不存在）
output_dir <- file.path(spatialfeature_dir, "individual_plots")
if (!dir.exists(output_dir)) dir.create(output_dir)

# 遍历每个基因
gene_cus %>% walk(function(gene) {
  # 遍历每个样本
  cells_boundary_ls %>% iwalk(function(cells_boundary, sample) {
    # 检查基因是否存在于数据中
    if (!gene %in% colnames(dat_expr_coords_full)) return()
    
    # 合并空间坐标与表达数据
    dat <- cells_boundary %>% 
      right_join(dat_expr_coords_full, by = "cell_id")
    
    # 检查空间坐标有效性
    if (all(is.na(dat$vertex_x))) return()
    if (any(is.na(dat$vertex_x))) {
      warning(sprintf("样本 %s 的 vertex_x 包含 NA 值", sample), immediate. = TRUE)
    }
    
    # 创建动态文件名
    pdf_file <- file.path(output_dir, sprintf("spatial_%s_%s.pdf", gene, sample))
    
    # 启动PDF设备
    pdf(pdf_file, width = 8, height = 6)  # 单图尺寸
    
    # 绘制单个基因-样本空间表达图
    p <- dat %>%
      arrange(.data[[gene]]) %>% 
      ggplot(aes(x = vertex_x, y = -vertex_y, color = .data[[gene]])) +
      scattermore::geom_scattermore(
        aes(alpha = .data[[gene]]), 
        pointsize = 1,  # 适当增大点尺寸
        show.legend = TRUE
      ) +
      scale_color_gradientn(
        colours = c("blue", "white", "red"),
        limits = c(0, 1),
        name = "Expression"
      ) +
      guides(alpha = "none") +
      labs(title = sprintf("%s 在 %s 中的表达", gene, sample)) +
      theme_void() +
      coord_equal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "right"
      )
    
    print(p)
    dev.off()
  })
})



#add FOV画ImageDimPlot----
source("/mnt/NexData_file/home/hyj/h5ad2rds(1).R")

xenium_rctd3 <-  JoinLayers(xenium_rctd,layers = 'count')
rownames(xenium_rctd3@meta.data) <-  paste0(xenium_rctd3@meta.data$sample,"_",xenium_rctd3@meta.data$cell_id)
rownames(xenium_rctd3@assays$RNA@layers$count) <- rownames(xenium_rctd3)
colnames(xenium_rctd3@assays$RNA@layers$count) <- colnames(xenium_rctd3)
xenium <- CreateSeuratObject(
  xenium_rctd3@assays$RNA@layers$count,
  assay = "Xenium")

xenium@meta.data <- xenium_rctd3@meta.data

umap <- xenium_rctd3@reductions$umap@cell.embeddings
rownames(umap) <- rownames(xenium@meta.data)
metadata <- xenium@meta.data 
xenium <- Add_meta_umapData_Xenium(object = xenium,umapfile = umap,metedata = metadata)
colnames(xenium@assays$Xenium) <-  colnames(xenium)
row.names(xenium@reductions$umap@cell.embeddings) <- colnames(xenium)

xenium <- RenameCells(xenium,new.names = colnames(xenium))
label <- unique(xenium@meta.data$sample)

###加入1个FOV测试----
library(arrow)
coord1 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041445__23_03967__20241107__051512",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "CYS_23_03967",object = xenium,cov = NULL)
xenium@images[["CYS-23-03967"]] = coord1

options(future.globals.maxSize = 10 * 1024^3)
xenium <- NormalizeData(xenium)

p <- ImageFeaturePlot(xenium, fov = "CYS-23-03967",features = c("CD3E"),max.cutoff = "q90",axes = FALSE,
                      cols = c('#003366', '#FFD700'))
ggsave(p,filename = paste0(i,"_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)





###加入18个样本FOV(正式)----

coord1 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041445__23_03967__20241107__051512",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "CYS_23_03967",object = xenium,cov = NULL)
  xenium@images[["CYS_23_03967"]] = coord1
  
coord2 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041445__23_06802__20241107__051512",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "GYX_23_06802",object = xenium,cov = NULL)
  xenium@images[["GYX_23_06802"]] = coord2
  
coord3 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041445__23_10358__20241107__051512",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "CYS_23_10358",object = xenium,cov = NULL)
  xenium@images[["CYS_23_10358"]] = coord3
  
coord4 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041445__23_11342__20241107__051512",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "FYN_23_11342",object = xenium,cov = NULL)
  xenium@images[["FYN_23_11342"]] = coord4

coord5 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041445__23_18277__20241107__051512",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "ZKJ_23_18277",object = xenium,cov = NULL)
  xenium@images[["ZKJ_23_18277"]] = coord5
  
coord6 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041445__23_18305__20241107__051512",assay = 'Xenium',
                             mols.qv.threshold = 20,label = "WSL_23_18305",object = xenium,cov = NULL)
xenium@images[["WSL_23_18305"]] = coord6

coord7 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041445__23_19464_24_03997__20241107__051512",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "GYX_24_03997",object = xenium,cov = NULL)
xenium@images[["GYX_24_03997"]] = coord7

coord8 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041445__23_19464_24_03997__20241107__051512",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "ZKJ_23_19464",object = xenium,cov = NULL)
xenium@images[["ZKJ_23_19464"]] = coord8

coord9 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041445__23_20895__20241107__051512",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "WSL_23_20895",object = xenium,cov = NULL)
xenium@images[["WSL_23_20895"]] = coord9

coord10 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041447__23_23753__20241107__051512",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "ZBH_23_23753",object = xenium,cov = NULL)
xenium@images[["ZBH_23_23753"]] = coord10

coord11 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041447__23_25688__20241107__051512",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "ZBH_23_25688",object = xenium,cov = NULL)
xenium@images[["ZBH_23_25688"]] = coord11

coord12 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041447__23_28671__20241107__051512",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "ZBH_23_28671",object = xenium,cov = NULL)
xenium@images[["ZBH_23_28671"]] = coord12

coord13 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041447__23_30679__20241107__051512",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "ZL_23_30769",object = xenium,cov = NULL)
xenium@images[["ZL_23_30769"]] = coord13

coord14 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041447__23_32098__20241107__051512",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "ZL_23_32098",object = xenium,cov = NULL)
xenium@images[["ZL_23_32098"]] = coord14

coord15 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041447__23_33594_24_01346__20241107__051512",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "HYY_24_01346",object = xenium,cov = NULL)
xenium@images[["HYY_24_01346"]] = coord15

coord16 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041447__23_33594_24_01346__20241107__051512",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "SQQ_23_33594",object = xenium,cov = NULL)
xenium@images[["SQQ_23_33594"]] = coord16

coord17 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041447__24_02802__20241107__051512",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "HYY_24_02802",object = xenium,cov = NULL)
xenium@images[["HYY_24_02802"]] = coord17

coord18 <- Read_fov_spatial(data.dir = "/data/xenium/20241107__051016__20241107_FWSP20241358/output-XETG00305__0041447__24_06201__20241107__051512",assay = 'Xenium',
                            mols.qv.threshold = 20,label = "CZY_24_06201",object = xenium,cov = NULL)
xenium@images[["CZY_24_06201"]] = coord18


options(future.globals.maxSize = 10 * 1024^3)
xenium <- NormalizeData(xenium)

saveRDS(xenium,paste0('xenium_coor',"_final.rds"))
xenium <- readRDS("xenium_coor_final.rds")

p <- ImageFeaturePlot(xenium, fov = "CZY_24_06201",features = c("CD3E"),max.cutoff = "q90",axes = FALSE,
                      cols = c('#003366', "#E7475E"))
p


for( i in unique(xenium@meta.data$sample)){
  p <- ImageFeaturePlot(xenium, fov = "i",features = gene_cus, max.cutoff = "q90",axes = FALSE,
                        cols = c('#003366', '#FFD700'))
  ggsave(p,filename = paste0("i","_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)
}


#####不带segmentation----
# 确保目标目录存在
if (!dir.exists(Imagefeature_dir)) dir.create(Imagefeature_dir, recursive = TRUE)

# 遍历每个样本和基因
for (sample_id in unique(xenium@meta.data$sample)) {
  for (gene in gene_cus) {
    tryCatch({
      # 生成当前基因的空间表达图
      p <- ImageFeaturePlot(
        object = xenium,
        fov = sample_id,
        features = gene,
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('#003366', '#FFD700'),  # 深蓝到金色渐变
        size = 0.8,  # 增大点尺寸增强可见性
        dark.background = TRUE  # 黑色背景增强对比度
      )
      
      # 构造保存路径：样本_基因.png
      file_name <- paste0(sample_id, "_", gene, "_spatial.png")
      save_path <- file.path(Imagefeature_dir, file_name)
      
      # 保存图像 (300dpi高清)
      ggsave(
        filename = save_path,
        plot = p,
        width = 12,
        height = 12,
        dpi = 300
      )
      
      # 进度提示
      message("已保存: ", save_path)
      
    }, error = function(e) {
      # 错误处理
      warning(sprintf("样本 %s 的基因 %s 绘图失败: %s", 
                      sample_id, gene, e$message))
    })
  }
}

#####带segmentation----
#根据需求调整gene名
immune_genes <- c("CD86", "CD80", "IL1A", "TLR2", "TLR4")
gene_cus <- immune_genes

SegImagefeature_dir <- 'result/SegImagefeature' 
dir.create(SegImagefeature_dir)

for (sample_id in unique(xenium@meta.data$sample)) {
  # 智能匹配样本ID
  j <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  
  if (length(j) == 0) {
    warning("跳过未匹配样本: ", sample_id)
    next
  }
  
  for (gene in gene_cus) {
    tryCatch({
      # 空值检查
      if (!is.null(xenium@images[[j]])) {
        DefaultBoundary(xenium@images[[j]]) <- "segmentation"
      }
      
      p <- ImageFeaturePlot(
        object = xenium,
        fov = j,
        nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
        border.color = "white",  # 细胞边界的颜色（深蓝色）
        border.size = 0.1,    # 细胞边界的线宽
        features = gene,
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('#003366', "#E7475E"),
        size = 0.8,
        dark.background = TRUE
      )
      
      ggsave(
        filename = file.path(SegImagefeature_dir, paste0(sample_id, "_", gene, "_spatial.png")),
        plot = p,
        width = 12,
        height = 12,
        dpi = 300
      )
      
      message("成功保存: ", sample_id, "_", gene)
      
    }, error = function(e) {
      warning("失败[", sample_id, ":", gene, "]: ", e$message)
    })
  }
}



#热图----
immune_genes <- c(
  # ==== T细胞功能 (T-cell function)
  'T_function' = "CD8A",      # MHC-I限制性细胞毒性T细胞核心标志[1](@ref)
  'T_function' = "CD3D",      # TCR-CD3复合物信号转导组分，参与T细胞激活[27](@ref)
  'T_function' = "CD3E",      # TCR-CD3复合物关键链，介导抗原识别信号[27](@ref)
  'T_function' = "CD4",       # MHC-II限制性T细胞共受体，增强TCR-pMHC结合[1](@ref)
  'T_function' = "GZMB",      # 颗粒酶B诱导靶细胞caspase依赖性凋亡[2](@ref)
  'T_function' = "GZMA",      # 颗粒酶A触发线粒体凋亡通路[2](@ref)
  'T_function' = "PRF1",      # 穿孔素形成膜孔道促进颗粒酶递送[2](@ref)
  'T_function' = "IFNG",      # γ干扰素激活巨噬细胞并增强抗原呈递[24,25](@ref)
  'T_function' = "TNF",       # 肿瘤坏死因子诱导靶细胞坏死并促进炎症[24](@ref)
  
  # ==== T细胞活化/衰竭 (T-cell activation/exhaustion)
  'T_AE' = "PDCD1",      # PD-1抑制TCR信号传导，介导耗竭表型[21,22](@ref)
  'T_AE' = "CTLA4",      # 竞争性抑制CD28共刺激信号[28](@ref)
  'T_AE' = "LAG3",       # 结合MHC-II抑制T细胞增殖[7](@ref)
  'T_AE' = "TIGIT",      # 干扰CD226激活信号，抑制NK/T细胞功能[8](@ref)
  'T_AE' = "TOX",        # 转录因子驱动耗竭程序关键调控因子[6](@ref)
  'T_AE' = "CD69",       # 早期活化标志，抑制S1P介导的细胞迁出[5](@ref)
  'T_AE' = "IL2RA",      # IL-2受体α链(CD25)，维持调节性T细胞功能[2](@ref)
  'T_AE' = "ICOS",       # 诱导共刺激分子促进滤泡辅助T细胞分化[5](@ref)
  
  # ==== 抗原呈递 (Antigen presentation)
  'AP' = "HLA-DRA",       # MHC-II α链，呈递外源性抗原肽[9,10](@ref)
  'AP' = "HLA-DRB1",      # MHC-II β链多态性决定抗原结合特异性[9](@ref)
  'AP' = "HLA-A",         # MHC-I经典分子呈递内源性抗原[9](@ref)
  'AP' = "HLA-B",         # MHC-I多态性影响CTL识别效率[9](@ref)
  'AP' = "HLA-C",         # MHC-I分子同时参与NK细胞抑制信号[11](@ref)
  'AP' = "CD74",          # MHC-II恒定链(li)，引导肽段装载[10](@ref)
  'AP' = "CIITA",         # MHC-II转录激活因子主调控因子[10](@ref)
  
  # ==== 巨噬细胞极化 (Macrophage polarization)
  'Mac_pol' = "CD68",     # 泛巨噬细胞溶酶体标志物[13](@ref)
  'Mac_pol' = "CD163",    # M2型血红蛋白清道夫受体[12](@ref)
  'Mac_pol' = "IL10",     # M2型抗炎因子抑制免疫应答[13](@ref)
  'Mac_pol' = "IL12B",    # M1型促炎因子驱动Th1分化[13](@ref)
  'Mac_pol' = "ARG1",     # M2型精氨酸酶抑制NO产生[12](@ref)
  'Mac_pol' = "NOS2",     # M1型一氧化氮合酶增强杀菌活性[13](@ref)
  'Mac_pol' = "CXCL9",    # IFN-γ诱导型T细胞趋化因子[15](@ref)
  'Mac_pol' = "CXCL10",   # M1型标志性IP-10趋化因子[15](@ref)
  
  # ==== 细胞因子/趋化因子 (Cytokines/Chemokines)
  'CC' = "CCL2",          # 单核细胞趋化蛋白-1招募髓系细胞[15](@ref)
  'CC' = "CCL5",          # RANTES吸引T细胞和嗜酸性粒细胞[16](@ref)
  'CC' = "CXCL9",         # IFN-γ诱导Th1细胞趋化[15,24](@ref)
  'CC' = "CXCL10",        # IP-10增强T细胞跨内皮迁移[16](@ref)
  'CC' = "CXCL13",        # B细胞趋化因子促进生发中心形成[5](@ref)
  'CC' = "IL6",           # 促炎因子激活STAT3通路[17](@ref)
  'CC' = "IL1B",          # 炎症小体下游关键细胞因子[14](@ref)
  'CC' = "TNF",           # 肿瘤坏死因子激活内皮细胞炎症[24](@ref)
  
  # ==== 免疫调节 (Immunoregulation)
  'Immune_reg' = "FOXP3",    # 调节性T细胞核心转录因子[2](@ref)
  'Immune_reg' = "IL10",     # 抑制性细胞因子限制炎症损伤[14](@ref)
  'Immune_reg' = "TGFB1",    # 转化生长因子β抑制效应T细胞[14](@ref)
  'Immune_reg' = "IDO1",     # 色氨酸分解酶诱导T细胞凋亡[8](@ref)
  'Immune_reg' = "CD274",    # PD-L1结合PD-1抑制T细胞活性[21,22](@ref)
  'Immune_reg' = "TNFRSF9"   # 4-1BB共刺激分子增强T细胞存活[28](@ref)
)
valid_immune_genes <- immune_genes[immune_genes %in% rownames(xenium)]
cat("有效基因:", length(valid_immune_genes), "/", length(immune_genes))
# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(immune_genes, valid_immune_genes)
gene <- valid_immune_genes

# 假设 immune_genes 是已分类的命名向量（如 'T_func' = "CD8A", ...）
ordered_gene <- names(gene)  # 直接提取基因名（忽略分类名）
# 或按分类组排序（保持组内顺序）：
ordered_gene <- unname(gene)  # 仅保留基因名，顺序与定义一致


#xenium <- ScaleData(xenium, vars.to.regress = c('nCount_RNA'))

H <- DoHeatmap(
  object = xenium, 
  features = ordered_gene,   # 关键：传入有序基因向量
  slot = "scale.data",    
  group.by = "sample", 
  group.colors = plan, 
  disp.min = -2.5, 
  disp.max = 2.5,
  angle = -30, 
  hjust = 0.8,
  size = 3,
  raster = TRUE               # 加速渲染（大数据集必备）[1,6](@ref)
)+
  scale_fill_gradientn(
    colors = c("navy", "white", "firebrick3"),
    values = scales::rescale(c(-2.5, 0, 2.5))  # 精确匹配disp.min/max[4](@ref)
  ) 

ggsave("result/immunegene_Heatmap1.png",plot = H,width = 12,height = 12,dpi = 300)

pdf("result/immunegene_Heatmap.pdf", width = 12, height = 10)
H
dev.off()



###基因集----
modules <- names(valid_immune_genes) |> unique()
modules

xenium <- AddModuleScore(object = xenium,
                         features = split(valid_immune_genes, modules),
                         name = "ImmuneModule_",
                         ctrl = 50,              # 对照基因数（建议50-100）
                         seed = 2025)             # 确保可重复性

# virus_module_cols <- c("IFN_I", "IFN_II", "PRR", "ASP", "CS_D", "IBA")
# # 批量删除这些列
# for (col in virus_module_cols) {
#   xenium@meta.data[[col]] <- NULL
# }


# 定义新旧列名映射关系
new_names <- c("T_function", "T_AE", "AP", "Mac_pol", "CC", "Immune_reg")
old_names <- paste0("ImmuneModule_", 1:6)

# 直接修改metadata列名
colnames(xenium@meta.data)[match(old_names, colnames(xenium@meta.data))] <- new_names


saveRDS(xenium,paste0('xenium_coor_group_Module',"_final.rds"))
# xenium <- readRDS("xenium_coor_group_final.rds")


###免疫评分空间featureplt----
ImmuneModuleSegImagefeature_dir <- 'result/ImmuneModuleSegImagefeature_18sample' 
dir.create(ImmuneModuleSegImagefeature_dir)

module_names <- new_names

for (sample_id in unique(xenium@meta.data$sample)) {
  # 智能匹配样本ID
  j <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  
  if (length(j) == 0) {
    warning("跳过未匹配样本: ", sample_id)
    next
  }
  
  for (module in module_names) {
    tryCatch({
      # 空值检查
      if (!is.null(xenium@images[[j]])) {
        DefaultBoundary(xenium@images[[j]]) <- "segmentation"
      }
      
      p <- ImageFeaturePlot(
        object = xenium,
        fov = j,
        nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
        border.color = "white",  # 细胞边界的颜色（深蓝色）
        border.size = 0.1,    # 细胞边界的线宽
        features = module,
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('#003366', "#E7475E"),
        size = 0.8,
        dark.background = TRUE
      )+
        ggtitle(paste("Module",module, "Score")) +
        theme(
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(size = 10)
        )
      
      ggsave(
        filename = file.path(ImmuneModuleSegImagefeature_dir, paste0(sample_id, "_", module,"_spatial.png")),
        plot = p,
        width = 12,
        height = 12,
        dpi = 300
      )
      
      message("成功保存: ", sample_id, "_", module)
      
    }, error = function(e) {
      warning("失败[", sample_id, ":", "ImmuneModule_",module,"]: ", e$message)
    })
  }
}

###免疫评分箱线图（带检验）----
library(Seurat)
library(ggpubr)
library(tidyverse)

ImmuneModuleBox_dir <- 'result/ImmuneModuleBox' 
dir.create(ImmuneModuleBox_dir)

# saveRDS(xenium,paste0('xenium_coor_group',"_final.rds"))
# xenium <- readRDS("xenium_coor_group_final.rds")

###循环箱线图（按group分）----
#按照new_names <- c("IFN_I", "IFN_II", "PRR", "ASP", "CS_D", "IBA")绘制"Effective"和"Ineffective"两组的箱线图，加入统计学检验
# 1. 提取模块评分数据和分组信息
# 从Seurat对象metadata中提取指定列（模块评分和分组信息）
plot_data <- xenium@meta.data %>%
  # 选择新命名的模块列（new_names）和分组列（group1）
  select(all_of(new_names), group1) %>%
  # 将宽数据转为长数据：多个模块列合并为两列（Module列存放模块名，Score列存放评分值）
  pivot_longer(cols = all_of(new_names), 
               names_to = "Module", 
               values_to = "Score") %>%
  # 将Module列转换为因子，固定模块顺序（避免绘图时按字母排序）
  mutate(Module = factor(Module, levels = new_names))

# 循环绘制每个模块
for (module in new_names) {
  # 1. 提取当前模块数据
  module_data <- plot_data %>%
    filter(Module == module)
  
  # 2. 动态计算y轴上限（用于统计标注位置）
  y_max <- max(module_data$Score, na.rm = TRUE) * 1.15
  
  # 3. 绘制单个模块箱线图
  p <- ggplot(module_data, aes(x = group1, y = Score, fill = group1)) +
    # 箱线图层
    geom_boxplot(
      width = 0.5,
      alpha = 0.85,
      outlier.shape = 19,
      outlier.size = 0.8,
      outlier.alpha = 0.6
    ) +
    # # 样本点层（显示每个样本）
    # geom_jitter(
    #   width = 0.15,
    #   size = 2,
    #   alpha = 0.5,
    #   color = "black"
    # ) +
    # # 均值标记（菱形符号）
    # stat_summary(
    #   fun = mean,
    #   geom = "point",
  #   shape = 18,
  #   size = 5,
  #   color = "gold"
  # ) +
  # 统计检验标注
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif", # 改为星号显示（*:p<0.05; **:p<0.01; ***:p<0.001）
    label.x = 1.5,
    size = 5,
    label.y = y_max,
    vjust = -0.5
  ) +
    # 配色方案
    scale_fill_manual(
      values = c("Effective" = "#4DBBD5", "Ineffective" = "#E64B35")
    ) +
    # 标签和标题
    labs(
      title = paste("Module:", module),
      #subtitle = "Distribution by Treatment Response",
      x = "Group",
      y = "Module Score",
      fill = "Group"
    ) +
    # 主题设置
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      #plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转45度并右对齐
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey70")
    )
  
  # 4. 保存图像
  ggsave(
    filename = file.path(ImmuneModuleBox_dir, paste0('group_', module, "_boxplot.png")),
    plot = p,
    width = 6,
    height = 7,
    dpi = 300
  )
  
  message("已保存: ", module, " 的箱线图")
}

###循环箱线图（按sample分）----
#按照new_names <- c("IFN_I", "IFN_II", "PRR", "ASP", "CS_D", "IBA")绘制"Effective"和"Ineffective"两组的箱线图，加入统计学检验

# 1. 提取模块评分数据和分组信息
# 从Seurat对象metadata中提取指定列（模块评分和分组信息）
plot_data <- xenium@meta.data %>%
  # 选择新命名的模块列（new_names）和分组列（group1）
  select(all_of(new_names), sample) %>%
  # 将宽数据转为长数据：多个模块列合并为两列（Module列存放模块名，Score列存放评分值）
  pivot_longer(cols = all_of(new_names), 
               names_to = "Module", 
               values_to = "Score") %>%
  # 将Module列转换为因子，固定模块顺序（避免绘图时按字母排序）
  mutate(Module = factor(Module, levels = new_names))

# 循环绘制每个模块
for (module in new_names) {
  # 1. 提取当前模块数据
  module_data <- plot_data %>%
    filter(Module == module)
  
  # 2. 动态计算y轴上限（用于统计标注位置）
  y_max <- max(module_data$Score, na.rm = TRUE) * 1.15
  
  # 3. 绘制单个模块箱线图
  p <- ggplot(module_data, aes(x = sample, y = Score, fill = sample)) +
    # 箱线图层
    geom_boxplot(
      width = 0.5,
      alpha = 0.85,
      outlier.shape = 19,
      outlier.size = 0.8,
      outlier.alpha = 0.6
    ) +
    # # 样本点层（显示每个样本）
    # geom_jitter(
    #   width = 0.15,
    #   size = 2,
    #   alpha = 0.5,
    #   color = "black"
    # ) +
    # # 均值标记（菱形符号）
    # stat_summary(
    #   fun = mean,
    #   geom = "point",
  #   shape = 18,
  #   size = 5,
  #   color = "gold"
  # ) +
  # 统计检验标注
  # stat_compare_means(
  #   method = "wilcox.test",
  #   label = "p.signif", # 改为星号显示（*:p<0.05; **:p<0.01; ***:p<0.001）
  #   label.x = 1.5,
  #   size = 5,
  #   label.y = y_max,
  #   vjust = -0.5
  # ) +
    # 配色方案
    scale_fill_manual(
      values = plan
    ) +
    # 标签和标题
    labs(
      title = paste("Module:", module),
      #subtitle = "Distribution by Treatment Response",
      x = "sample",
      y = "Module Score",
      #fill = "sample"
    ) +
    # theme(legend.position = "none")
    # 主题设置
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      #plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转45度并右对齐
      axis.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey70")
    )
  
  # 4. 保存图像
  ggsave(
    filename = file.path(ImmuneModuleBox_dir, paste0('sample_',module, "_boxplot.png")),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  message("已保存: ", module, " 的箱线图")
}



#5、治疗有效组治疗前后的病毒反应基因空间表达图和病毒反应评分比较----
#病毒反应基因空间表达图----
VriusSegImagefeature_dir <- 'result/VriusSegImagefeature' 
dir.create(VriusSegImagefeature_dir)

sample_5 <- c("GYX_23_06802", "FYN_23_11342", "WSL_23_18305", "GYX_24_03997", "WSL_23_20895")

virus_genes <- c(
  # ==== Type I Interferon Response (I型干扰素响应)
  'IFN_I' = "IFIT1",     # 阻断病毒核糖体组装
  'IFN_I' = "IFIT2",     # 抑制eIF4E翻译起始
  'IFN_I' = "IFIT3",     # 稳定IFIT1抗病毒复合物
  'IFN_I' = "OAS1",      # 激活RNase L降解病毒RNA
  'IFN_I' = "OAS2",      # 增强病毒RNA降解效率
  'IFN_I' = "OAS3",      # 协同激活RNase L通路
  'IFN_I' = "MX1",       # GTP酶阻断病毒核转运
  'IFN_I' = "MX2",       # 抑制HIV病毒核输入
  'IFN_I' = "ISG15",     # 泛素化修饰病毒蛋白
  'IFN_I' = "ISG20",     # 外切酶降解病毒ssRNA
  'IFN_I' = "IFI27",     # 促进线粒体凋亡增强清除
  'IFN_I' = "IFI44",     # 抑制病毒mRNA出核
  'IFN_I' = "IFI44L",    # 干扰病毒复制复合体
  'IFN_I' = "IFI6",      # 拮抗凋亡维持细胞存活
  'IFN_I' = "IFITM1",    # 阻断病毒膜融合过程
  'IFN_I' = "RSAD2",     # 抑制病毒聚合酶活性
  
  # ==== Type II Interferon Response (II型干扰素响应)
  'IFN_II' = "GBP1",     # 破坏病原体膜结构
  'IFN_II' = "GBP2",     # 激活NLRP3炎症小体
  'IFN_II' = "GBP4",     # 抑制病毒复制周期
  'IFN_II' = "GBP5",     # 增强抗原呈递效率
  'IFN_II' = "CXCL9", # 招募CD8+ T细胞至感染部位
  'IFN_II' = "CXCL10",# 增强T细胞跨内皮迁移
  'IFN_II' = "CXCL11",# 抑制血管生成限制扩散
  'IFN_II' = "STAT1",    # 介导IFN-γ信号传导
  'IFN_II' = "IRF1",     # 上调MHC-I表达
  'IFN_II' = "IRF8",     # 促进树突细胞成熟
  
  # ==== Pattern Recognition Receptors (模式识别受体) 
  'PRR' = "DDX58",       # RIG-I识别5'-三磷酸RNA
  'PRR' = "IFIH1",       # MDA5感知长链dsRNA
  'PRR' = "TLR3",        # 内体区识别病毒dsRNA
  'PRR' = "TLR7",        # 感知ssRNA病毒
  'PRR' = "TLR8",        # 识别病毒ssRNA激活NF-κB
  'PRR' = "TLR9",        # 结合未甲基化CpG DNA
  'PRR' = "TMEM173",     # STING激活cGAS通路
  
  # ==== Antiviral Signaling pathway (抗病毒信号通路)
  'ASP' = "IRF3",# 启动IFN-β转录
  'ASP' = "IRF7",# 扩增I型干扰素表达
  'ASP' = "IRF9",# ISGF3复合体组分
  'ASP' = "STAT1",# 介导ISGF3/GAS形成
  'ASP' = "STAT2",# ISGF3转录复合体核心
  'ASP' = "JAK1", # 磷酸化STAT蛋白
  'ASP' = "TYK2", # 介导IFNAR下游信号
  
  # ==== Cell stress/death (细胞应激/死亡)
  'CS_D' = "CASP1",# 炎症性caspase切割GSDMD
  'CS_D' = "CASP4",# 非经典炎症小体效应器
  'CS_D' = "GSDMD", # 焦亡孔道形成蛋白
  'CS_D' = "GSDME", # 次级焦亡执行蛋白
  'CS_D' = "MLKL",  # 坏死性凋亡终末执行者
  'CS_D' = "BCL2",   # 抗凋亡维持稳态
  'CS_D' = "BAX",    # 促凋亡诱导线粒体通透
  'CS_D' = "BAK1",   # 与BAX协同形成膜孔
  'CS_D' = "RIPK1",  # 坏死性凋亡调控节点
  'CS_D' = "RIPK3",  # 磷酸化MLKL诱导坏死
  'CS_D' = "NLRP3", # 炎症小体传感器*
  
  # ==== Inflammatory body activation (炎症小体)
  'IBA' = "NLRP3", # 识别PAMPs/DAMPs*
  'IBA' = "NLRC4", # 感知胞内鞭毛蛋白
  'IBA' = "PYCARD",# ASC衔接蛋白桥接caspase
  'IBA' = "IL1B",  # 促炎细胞因子
  'IBA' = "IL18"   # 诱导IFN-γ增强免疫
)



valid_Vriusgenes <- virus_genes[virus_genes %in% rownames(xenium)]
cat("有效基因:", length(valid_Vriusgenes), "/", length(virus_genes))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_genes <- setdiff(virus_genes, valid_Vriusgenes)

gene_list <- valid_Vriusgenes

for (sample_id in sample_5) {
  # 智能匹配样本ID
  j <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  
  if (length(j) == 0) {
    warning("跳过未匹配样本: ", sample_id)
    next
  }
  
  for (gene in gene_list) {
    tryCatch({
      # 空值检查
      if (!is.null(xenium@images[[j]])) {
        DefaultBoundary(xenium@images[[j]]) <- "segmentation"
      }
      
      p <- ImageFeaturePlot(
        object = xenium,
        fov = j,
        nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
        border.color = "white",  # 细胞边界的颜色（深蓝色）
        border.size = 0.1,    # 细胞边界的线宽
        features = gene,
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('#003366', "#E7475E"),
        size = 0.8,
        dark.background = TRUE
      )
      
      ggsave(
        filename = file.path(VriusSegImagefeature_dir, paste0(sample_id, "_", gene, "_spatial.png")),
        plot = p,
        width = 12,
        height = 12,
        dpi = 300
      )
      
      message("成功保存: ", sample_id, "_", gene)
      
    }, error = function(e) {
      warning("失败[", sample_id, ":", gene, "]: ", e$message)
    })
  }
}

#病毒反应评分比较----
###新增分组：Treatment is Effective or Ineffective----
group_map <- c(
  "CYS_23_03967" = "Ineffective",
  "GYX_23_06802" = "Effective",
  "CYS_23_10358" = "Ineffective",
  "FYN_23_11342" = "Effective",
  "ZKJ_23_18277" = "Ineffective",
  "WSL_23_18305" = "Effective",
  "GYX_24_03997" = "Effective",
  "ZKJ_23_19464" = "Ineffective",
  "WSL_23_20895" = "Effective",
  "ZBH_23_23753" = "Ineffective",
  "ZBH_23_25688" = "Ineffective",
  "ZBH_23_28671" = "Ineffective",
  "ZL_23_30769" = "Ineffective",
  "ZL_23_32098" = "Ineffective",
  "HYY_24_01346" = "Ineffective",
  "SQQ_23_33594" = "Ineffective",
  "HYY_24_02802" = "Ineffective",
  "CZY_24_06201" = "Ineffective"
)

# 1. 确认样本标识列名（假设为 'orig.ident'）
sample_col <- "sample"  # 若列名为 'sample' 则修改为 "sample"

# 2. 提取当前样本标识向量
samples <- as.character(xenium@meta.data[[sample_col]])

# 3. 创建分组向量（按样本名精确匹配）
group_vector <- group_map[samples]
names(group_vector) <- colnames(xenium)  # 确保分组向量与细胞barcode一一对应

# 4. 使用 AddMetaData 添加分组列
xenium <- AddMetaData(
  object = xenium,
  metadata = group_vector,
  col.name = "group1"
)

# #删除
# xenium@meta.data["group1"] <- NULL

# # 统计每个 group 下有多少个样本及细胞数
# group_sample_summary <- xenium@meta.data %>%
#   group_by(group1, sample) %>%
#   summarise(cell_count = n()) %>%  # 统计每个组合的细胞数
#   group_by(group1) %>%
#   summarise(
#     sample_count = n_distinct(sample),  # 每个 group.y 下的样本数
#     total_cells = sum(cell_count),      # 每个 group.y 下的总细胞数
#     samples = paste(unique(sample), collapse = ", ")  # 列出样本名称
#   )
# 
# print(group_sample_summary)

###基因集----
modules <- names(valid_Vriusgenes) |> unique()
modules

xenium <- AddModuleScore(object = xenium,
                         features = split(valid_Vriusgenes, modules),
                         name = "VirusModule_",
                         ctrl = 50,              # 对照基因数（建议50-100）
                         seed = 2025)             # 确保可重复性

# virus_module_cols <- c("IFN_I", "IFN_II", "PRR", "ASP", "CS_D", "IBA")
# # 批量删除这些列
# for (col in virus_module_cols) {
#   xenium@meta.data[[col]] <- NULL
# }


# 定义新旧列名映射关系
new_names <- c("IFN_I", "IFN_II", "PRR", "ASP", "CS_D", "IBA")
old_names <- paste0("VirusModule_", 1:6)

# 直接修改metadata列名
colnames(xenium@meta.data)[match(old_names, colnames(xenium@meta.data))] <- new_names

###病毒评分空间featureplt----
VriusModuleSegImagefeature_dir <- 'result/VriusModuleSegImagefeature_18sample' 
dir.create(VriusModuleSegImagefeature_dir)

module_names <- new_names

for (sample_id in unique(xenium@meta.data$sample)) {
  # 智能匹配样本ID
  j <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  
  if (length(j) == 0) {
    warning("跳过未匹配样本: ", sample_id)
    next
  }
  
  for (module in module_names) {
    tryCatch({
      # 空值检查
      if (!is.null(xenium@images[[j]])) {
        DefaultBoundary(xenium@images[[j]]) <- "segmentation"
      }
      
      p <- ImageFeaturePlot(
        object = xenium,
        fov = j,
        nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
        border.color = "white",  # 细胞边界的颜色（深蓝色）
        border.size = 0.1,    # 细胞边界的线宽
        features = module,
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('#003366', "#E7475E"),
        size = 0.8,
        dark.background = TRUE
      )+
        ggtitle(paste(module, "Score")) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right",
          legend.title = element_text(size = 10)
        )
      
      ggsave(
        filename = file.path(VriusModuleSegImagefeature_dir, paste0(sample_id, "_", module,"_spatial.png")),
        plot = p,
        width = 12,
        height = 12,
        dpi = 300
      )
      
      message("成功保存: ", sample_id, "_", module)
      
    }, error = function(e) {
      warning("失败[", sample_id, ":", "VirusModule_",module,"]: ", e$message)
    })
  }
}

###病毒反应评分箱线图（带检验）----
library(Seurat)
library(ggpubr)
library(tidyverse)

VriusModuleBox_dir <- 'result/VriusModuleBox' 
dir.create(VriusModuleBox_dir)

# saveRDS(xenium,paste0('xenium_coor_group',"_final.rds"))
# xenium <- readRDS("xenium_coor_group_final.rds")

###循环箱线图（按group分）----
#按照new_names <- c("IFN_I", "IFN_II", "PRR", "ASP", "CS_D", "IBA")绘制"Effective"和"Ineffective"两组的箱线图，加入统计学检验
# 1. 提取模块评分数据和分组信息
# 从Seurat对象metadata中提取指定列（模块评分和分组信息）
plot_data <- xenium@meta.data %>%
  # 选择新命名的模块列（new_names）和分组列（group1）
  select(all_of(new_names), group1) %>%
  # 将宽数据转为长数据：多个模块列合并为两列（Module列存放模块名，Score列存放评分值）
  pivot_longer(cols = all_of(new_names), 
               names_to = "Module", 
               values_to = "Score") %>%
  # 将Module列转换为因子，固定模块顺序（避免绘图时按字母排序）
  mutate(Module = factor(Module, levels = new_names))

# 循环绘制每个模块
for (module in new_names) {
  # 1. 提取当前模块数据
  module_data <- plot_data %>%
    filter(Module == module)
  
  # 2. 动态计算y轴上限（用于统计标注位置）
  y_max <- max(module_data$Score, na.rm = TRUE) * 1.15
  
  # 3. 绘制单个模块箱线图
  p <- ggplot(module_data, aes(x = group1, y = Score, fill = group1)) +
    # 箱线图层
    geom_boxplot(
      width = 0.5,
      alpha = 0.85,
      outlier.shape = 19,
      outlier.size = 0.8,
      outlier.alpha = 0.6
    ) +
    # # 样本点层（显示每个样本）
    # geom_jitter(
    #   width = 0.15,
    #   size = 2,
    #   alpha = 0.5,
    #   color = "black"
    # ) +
    # # 均值标记（菱形符号）
    # stat_summary(
    #   fun = mean,
    #   geom = "point",
  #   shape = 18,
  #   size = 5,
  #   color = "gold"
  # ) +
  # 统计检验标注
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif", # 改为星号显示（*:p<0.05; **:p<0.01; ***:p<0.001）
    label.x = 1.5,
    size = 5,
    label.y = y_max,
    vjust = -0.5
  ) +
    # 配色方案
    scale_fill_manual(
      values = c("Effective" = "#4DBBD5", "Ineffective" = "#E64B35")
    ) +
    # 标签和标题
    labs(
      title = paste("Module:", module),
      #subtitle = "Distribution by Treatment Response",
      x = "Group",
      y = "Module Score",
      fill = "Group"
    ) +
    # 主题设置
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      #plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转45度并右对齐
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey70")
    )
  
  # 4. 保存图像
  ggsave(
    filename = file.path(VriusModuleBox_dir, paste0('group_', module, "_boxplot.png")),
    plot = p,
    width = 6,
    height = 7,
    dpi = 300
  )
  
  message("已保存: ", module, " 的箱线图")
}

###循环箱线图（按sample分）----
#按照new_names <- c("IFN_I", "IFN_II", "PRR", "ASP", "CS_D", "IBA")绘制"Effective"和"Ineffective"两组的箱线图，加入统计学检验

# 1. 提取模块评分数据和分组信息
# 从Seurat对象metadata中提取指定列（模块评分和分组信息）
plot_data <- xenium@meta.data %>%
  # 选择新命名的模块列（new_names）和分组列（group1）
  select(all_of(new_names), sample) %>%
  # 将宽数据转为长数据：多个模块列合并为两列（Module列存放模块名，Score列存放评分值）
  pivot_longer(cols = all_of(new_names), 
               names_to = "Module", 
               values_to = "Score") %>%
  # 将Module列转换为因子，固定模块顺序（避免绘图时按字母排序）
  mutate(Module = factor(Module, levels = new_names))

# 循环绘制每个模块
for (module in new_names) {
  # 1. 提取当前模块数据
  module_data <- plot_data %>%
    filter(Module == module)
  
  # 2. 动态计算y轴上限（用于统计标注位置）
  y_max <- max(module_data$Score, na.rm = TRUE) * 1.15
  
  # 3. 绘制单个模块箱线图
  p <- ggplot(module_data, aes(x = sample, y = Score, fill = sample)) +
    # 箱线图层
    geom_boxplot(
      width = 0.5,
      alpha = 0.85,
      outlier.shape = 19,
      outlier.size = 0.8,
      outlier.alpha = 0.6
    ) +
    # # 样本点层（显示每个样本）
    # geom_jitter(
    #   width = 0.15,
    #   size = 2,
    #   alpha = 0.5,
    #   color = "black"
    # ) +
    # # 均值标记（菱形符号）
    # stat_summary(
    #   fun = mean,
    #   geom = "point",
  #   shape = 18,
  #   size = 5,
  #   color = "gold"
  # ) +
  # 统计检验标注
  # stat_compare_means(
  #   method = "wilcox.test",
  #   label = "p.signif", # 改为星号显示（*:p<0.05; **:p<0.01; ***:p<0.001）
  #   label.x = 1.5,
  #   size = 5,
  #   label.y = y_max,
  #   vjust = -0.5
  # ) +
  # 配色方案
  scale_fill_manual(
    values = plan
  ) +
    # 标签和标题
    labs(
      title = paste("Module:", module),
      #subtitle = "Distribution by Treatment Response",
      x = "sample",
      y = "Module Score",
      #fill = "sample"
    ) +
    # theme(legend.position = "none")
    # 主题设置
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      #plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转45度并右对齐
      axis.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey70")
    )
  
  # 4. 保存图像
  ggsave(
    filename = file.path(VriusModuleBox_dir, paste0('sample_',module, "_boxplot.png")),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  message("已保存: ", module, " 的箱线图")
}



#6、WSL/GYX治疗前后相互间不同细胞柱状图、箱线图----
#   治疗前         治疗后
#"WSL_23_18305","WSL_23_20895"
#"GYX_24_03997","GYX_23_06802"

xenium@meta.data$sample %>% unique()

# 指定样本顺序
sample_order <- c("GYX_24_03997", "GYX_23_06802")

# 细胞比例统计图----
celltype_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  filter(sample %in% sample_order) %>%
  mutate(
    # 将group.x转换为因子并指定顺序
    sample = factor(sample, levels = sample_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(first_type, sample) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(sample) %>%
  mutate(prop = n / sum(n))

# 基础 ggplot2 堆叠柱状图
p1 <- ggplot(celltype_counts, aes(x = .data[['sample']], y = prop, fill = .data[['first_type']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = color_mapping) +
  # 坐标轴与图例优化
  labs(x = "sample", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(
  filename = "result/Cell_prop/6.2GYX_Celltype_proportions.pdf",  # 保存文件名
  plot = p1,                                # 指定要保存的图形对象
  width = 5,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
)


# "WSL_23_18305","WSL_23_20895"
# "GYX_24_03997","GYX_23_06802"

# 指定样本顺序
sample_order <- c("WSL_23_18305","WSL_23_20895")

#各细胞箱线图----
# 数据准备：计算每个样本的细胞比例
celltype_sample_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%
  filter(sample %in% sample_order) %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(sample, first_type) %>% 
  tally() %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n))

box_dir <- 'result/boxplot/6.1WSL'
dir.create(box_dir)

library(reshape2)
celltype_wide <- dcast(celltype_sample_counts, 
                       sample ~ first_type, 
                       value.var = "prop")

#画这种循环的需要把细胞名里的+和空格改掉
# 创建安全的列名映射（特殊字符转下划线）
safe_names <- gsub("[^[:alnum:]]", "_", cell_types_order) %>%
  make.unique()  # 避免重名冲突

# 重命名数据框中的列
colnames(celltype_wide) <- ifelse(
  colnames(celltype_wide) %in% cell_types_order,
  safe_names[match(colnames(celltype_wide), cell_types_order)],
  colnames(celltype_wide)  # 保留非细胞类型列
)

# 更新细胞类型顺序向量
safe_cell_types_order <- safe_names

library(ggpubr)

for (i in seq_along(cell_types_order)) {
  orig_name <- cell_types_order[i]       # 原始名称（用于标题）
  safe_name <- safe_cell_types_order[i]  # 安全名称（用于数据引用）
  
  if (safe_name %in% colnames(celltype_wide)) {
    p <- ggbarplot(
      celltype_wide, 
      x = "sample", 
      y = safe_name,  # 使用安全列名
      fill = "sample", 
      palette = c("#4DBBD5", "#E64B35"),
    ) +
      # stat_compare_means(
      #   method = "wilcox.test",
      #   label = "p.signif", # 改为星号显示（*:p<0.05; **:p<0.01; ***:p<0.001
      #   label.y = max(celltype_wide[[safe_name]], na.rm = TRUE) * 1.1
      # ) +
      labs(
        title = paste("Proportion of", orig_name),  # 标题用原始名称
        x = "", 
        y = "Proportion"
      ) +
      theme(legend.position = "none")
    
    # 文件名也使用安全名称
    ggsave(
      filename = file.path(box_dir, paste0("WSL_", safe_name, "_Boxplot.pdf")),
      plot = p,
      width = 4, 
      height = 5
    )
  }
}


# 细胞比例统计图(免疫细胞内)----
#   治疗前         治疗后
# "WSL_23_18305","WSL_23_20895"
# "GYX_24_03997","GYX_23_06802"
xenium@meta.data$sample %>% unique()
# 指定样本顺序
sample_order <- c("GYX_24_03997","GYX_23_06802")

cellgroup_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  filter(cell_group != "Non_immune_cell") %>% # 过滤非免疫细胞
  filter(sample %in% sample_order) %>%
  mutate(
    # 将group.x转换为因子并指定顺序
    sample = factor(sample, levels = sample_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(first_type, sample) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(sample) %>%
  mutate(prop = n / sum(n))

# 基础 ggplot2 堆叠柱状图
p1 <- ggplot(cellgroup_counts, aes(x = .data[['sample']], y = prop, fill = .data[['first_type']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = color_mapping) +
  # 坐标轴与图例优化
  labs(x = "Sample", y = "Cell_group Proportion", fill = "Immune_cell") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(
  filename = "result/Cell_prop/6.2GYX_Immune_cellgroup_proportions.pdf",  # 保存文件名
  plot = p1,                                # 指定要保存的图形对象
  width = 5,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
  )





#7、除WSL/GYX,其他患者治疗前后细胞组成变化柱状图，箱线图+空转映射----
sss <- xenium@meta.data$sample %>% unique()
library(stringr)
sorted_sss <- sss[order(
  str_extract(sss, "^[A-Z]+") # 提取前缀
)]
sorted_sss

#   治疗前         治疗后
#"WSL_23_18305","WSL_23_20895"###
#"GYX_24_03997","GYX_23_06802"###

#"CYS_23_03967","CYS_23_10358"
#"ZKJ_23_18277","ZKJ_23_19464"
#"ZBH_23_23753","ZBH_23_28671","ZBH_23_25688"
#"ZL_23_30769","ZL_23_32098"
#"HYY_24_01346","HYY_24_02802"

# 指定样本顺序
sample_order <- c("HYY_24_01346","HYY_24_02802")

# 细胞比例统计图----
celltype_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  filter(sample %in% sample_order) %>%
  mutate(
    # 将group.x转换为因子并指定顺序
    sample = factor(sample, levels = sample_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(first_type, sample) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(sample) %>%
  mutate(prop = n / sum(n))

# 基础 ggplot2 堆叠柱状图
p1 <- ggplot(celltype_counts, aes(x = .data[['sample']], y = prop, fill = .data[['first_type']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = color_mapping) +
  # 坐标轴与图例优化
  labs(x = "sample", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(
  filename = "result/Cell_prop/7.5HYY_Celltype_proportions.pdf",  # 保存文件名
  plot = p1,                                # 指定要保存的图形对象
  width = 5,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
)



#各细胞箱线图----
# 数据准备：计算每个样本的细胞比例
celltype_sample_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%
  filter(sample %in% sample_order) %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(sample, first_type) %>% 
  tally() %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n))

box_dir <- 'result/boxplot/7.5HYY'
dir.create(box_dir)

library(reshape2)
celltype_wide <- dcast(celltype_sample_counts, 
                       sample ~ first_type, 
                       value.var = "prop")

#画这种循环的需要把细胞名里的+和空格改掉
# 创建安全的列名映射（特殊字符转下划线）
safe_names <- gsub("[^[:alnum:]]", "_", cell_types_order) %>%
  make.unique()  # 避免重名冲突

# 重命名数据框中的列
colnames(celltype_wide) <- ifelse(
  colnames(celltype_wide) %in% cell_types_order,
  safe_names[match(colnames(celltype_wide), cell_types_order)],
  colnames(celltype_wide)  # 保留非细胞类型列
)

# 更新细胞类型顺序向量
safe_cell_types_order <- safe_names

library(ggpubr)

for (i in seq_along(cell_types_order)) {
  orig_name <- cell_types_order[i]       # 原始名称（用于标题）
  safe_name <- safe_cell_types_order[i]  # 安全名称（用于数据引用）
  
  if (safe_name %in% colnames(celltype_wide)) {
    p <- ggbarplot(
      celltype_wide, 
      x = "sample", 
      y = safe_name,  # 使用安全列名
      fill = "sample", 
      palette = c("#4DBBD5", "#E64B35"
                  #,"#FFB85F"
                  ),
    ) +
      # stat_compare_means(
      #   method = "wilcox.test",
      #   label = "p.signif", # 改为星号显示（*:p<0.05; **:p<0.01; ***:p<0.001
      #   label.y = max(celltype_wide[[safe_name]], na.rm = TRUE) * 1.1
      # ) +
      labs(
        title = paste("Proportion of", orig_name),  # 标题用原始名称
        x = "", 
        y = "Proportion"
      ) +
      theme(legend.position = "none")
    
    # 文件名也使用安全名称
    ggsave(
      filename = file.path(box_dir, paste0("HYY_", safe_name, "_Boxplot.pdf")),
      plot = p,
      width = 4, 
      height = 5
    )
  }
}


# 细胞比例统计图(免疫细胞内)----

cellgroup_counts <- xenium@meta.data %>%
  filter(first_type != "Others") %>%  # 过滤掉"Others"细胞
  filter(cell_group != "Non_immune_cell") %>% # 过滤非免疫细胞
  filter(sample %in% sample_order) %>%
  mutate(
    # 将group.x转换为因子并指定顺序
    sample = factor(sample, levels = sample_order),
    first_type = factor(first_type, levels = cell_types_order)
  ) %>%
  group_by(first_type, sample) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(sample) %>%
  mutate(prop = n / sum(n))

# 基础 ggplot2 堆叠柱状图
p1 <- ggplot(cellgroup_counts, aes(x = .data[['sample']], y = prop, fill = .data[['first_type']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = color_mapping) +
  # 坐标轴与图例优化
  labs(x = "Sample", y = "Cell_group Proportion", fill = "Immune_cell") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(
  filename = "result/Cell_prop/7.5HYY_Immune_cellgroup_proportions.pdf",  # 保存文件名
  plot = p1,                                # 指定要保存的图形对象
  width = 5,                              # 宽度（默认单位为英寸）
  height = 8,                              # 高度
  dpi = 300,                               # 分辨率（dpi）
  bg = "white"                             # 背景颜色
)

###Imagedimplot----
# library(RColorBrewer)
# n_groups <- length(unique(xenium$SubtypeT_Fib))  # 获取分组数量
# colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_groups)

cell_types_order <- c(
  # 免疫细胞
  "Natural killer cells",  # 先天免疫核心[6](@ref)
  "CD8+ T cells",          # 细胞毒性T细胞
  "CD4+ T cells",          # 辅助T细胞
  "Tregs",                 # 调节性T细胞
  "B cells",               # B淋巴细胞
  "Plasma cells",          # B细胞终末分化态[8](@ref)
  "DC",                    # 树突状细胞（抗原提呈）
  "pDC",                   # 浆细胞样树突细胞
  "Mono-Mac",              # 单核/巨噬细胞（吞噬功能）[5](@ref)
  "Mast cells",            # 肥大细胞（过敏反应）
  "LC",                    # 朗格汉斯细胞（皮肤免疫）
  
  # 非免疫
  "Epithelial",
  #"Malignant cells",       # 肿瘤细胞
  "Proliferating cells",   # 增殖态细胞（可含免疫细胞）
  "Endothelial cells",     # 血管内皮细胞[8](@ref)
  "Fibroblasts",           # 成纤维细胞（结缔组织）
  #"Keratinocytes",         # 角质细胞（表皮）
  "Melanocytes",           # 黑色素细胞
  "Others"                 # 未定义细胞
)

# 将元数据中的细胞类型转换为有序因子
xenium@meta.data$first_type_coEpi <- factor(
  xenium@meta.data$first_type_coEpi,
  levels = cell_types_order
)

color_mapping <- setNames(
  c(
    "#FF6D31", "#FFB85F", "#FFCB18", "#FBF579", "#F0D879", "#FEEE7D", 
    "#ED5485", "#E8A0B8", "#D9534F", "#EF5285", "#F17D80",
    "#7F95D1", "#02C9C9", "#4695D6", "#73B66B", 
    #"#8BC24C", 
    "#600473", "#999999"
  ),
  cell_types_order
)

sample_label <- unique(xenium@meta.data$sample)

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = T,
                    border.color = "white", 
                    border.size = 0.1,
                    cols = color_mapping,
                    #shuffle.cols = T, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "first_type_coEpi",
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0("RCTD/result/Imagedimplottest/",i,"_spatial.png"),width = 12 ,height = 12,dpi = 300)
}



#8、关键免疫细胞和肿瘤细胞的空间关系分析----
# mapping <- data.frame(
#   orig.ident = c("23_03967", "23_06802", "23_10358", "23_11342", 
#                  "23_18277", "23_18305", "23_19464_24_03997", 
#                  "23_20895", "23_23753", "23_25688", "23_28671", 
#                  "23_30679", "23_32098", "23_33594_24_01346", 
#                  "24_02802", "24_06201"),
#   Slide_ID = c(rep("ID0041445", 8), rep("ID0041447", 8))
# )

# 在 all_data@meta.data 中添加新列
#xenium@meta.data$Slide_ID <- mapping$Slide_ID[match(xenium@meta.data$orig.ident, mapping$orig.ident)]

# 验证结果
#table(xenium@meta.data$Slide_ID,xenium@meta.data$orig.ident)


#meta.data$cellid需要和cells_position_ls里的cell_id一致
xenium@meta.data$cell_id <- paste(
  xenium$project,
  xenium$cell_id2,
  sep = "_"
  )

xenium@meta.data %>% head()

rownames(xenium@meta.data) %>% head()

cells_boundary_ls<-readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
cells_position_ls<-readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))

library(Seurat)
library(tidyverse)

source("spatial analysis.R")

# proximity_out_dir <- "proximity_out_dir"
# dir.create(proximity_out_dir)

###sample名称和slide名称需一致，slide名称是cells_position_ls各FOV的名字
# Proximity analysis -------
# 1. input data ------------
cells_position_ls %>% head() 

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- xenium[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(first_type_coEpi, sample, "cell_id", "group1") %>% 
      #rownames_to_column("cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(first_type_coEpi = factor(first_type_coEpi)) %>%
      dplyr::rename(Classification1 = first_type_coEpi,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample,
                    Group = group1)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos %>% head()
all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
all_data_meta_pos <- all_data_meta_pos %>%
  as.data.frame() %>%       # 确保为数据框类型
  remove_rownames()     # 强制清除现有行名


# ########### only select 100,000 cells ##############
# all_data_meta_pos_subset <- 
#   all_data_meta_pos %>% 
#   rownames_to_column("cell_id") %>%
#   slice_sample(n = 100000)


# select ROI cells --------------------------------------------------------------
# all_tumor_roi_selected <- readRDS('tumor_roi_selected_result3/all_tumor_roi_selected3.rds')
# all_data_meta_pos_subset <- 
#   all_data_meta_pos %>% 
#   rownames_to_column("cell_id") %>% 
#   dplyr::filter(cell_id %in% Cells(all_tumor_roi_selected))

# all_data_meta_pos_subset %>% head

# cell_id        X        Y Classification1 group       Slide_ID
# 1 1_aaaaaagg-1 1670.389 185.1358      Macrophage  CD-F X_202413076_32
# 2 1_aaadlghn-1 1681.172 176.4662          Plasma  CD-F X_202413076_32
# 3 1_aaadmplc-1 1677.561 182.6699               T  CD-F X_202413076_32
# 4 1_aaadoogj-1 1677.109 165.6637          Plasma  CD-F X_202413076_32
# 5 1_aaadpljf-1 1671.887 176.8148          Plasma  CD-F X_202413076_32
# 6 1_aaafaiko-1 1719.466 192.8071          Plasma  CD-F X_202413076_32
# ......


# 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

## 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )



# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}


# res_pa <- readRDS('proximity/res_pa.rds')


# summarized by group and celltype
pa_groups <-
  with(
    all_data_meta_pos,
    str_c(Classification1, Group, Slide_ID, sep = '..')
  )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y")+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path("result/proximity", "proximity analysis.pdf"),
  height = 14,
  width  = 18
)
p_res_pa %>% print()
dev.off()



## 2.2 proximity analysis - series analysis -----
init_cells <- c(
  '10'  = 750,
  '50'  = 800,
  '100' = 1000,
  '150' = 1200,
  '200' = 1500
)

res_pa_ls <-
  c(10, 50, 100, 150, 200) %>%
  set_names(., .) %>%
  map(function(x){
    k = init_cells[as.character(x)]
    res_pa <-
      ProximutyAnalysis_NN_multi(
        all_data_meta_pos,
        query = all_data_meta_pos,
        x = "X",
        y = "Y",
        group.by = "Slide_ID",
        classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
        k = k, # 'auto',
        searchtype = "radius",
        radius = x, # 邻域半径
        mpp = 1, # !important, 1 for xenium
        returntype = "count",
        cell_enough_ratio_cutoff = 90
      )
    res_pa
  })

# res_pa_ls <- readRDS("data/res_pa_ls.rds")

# 将同心圆统计 转换成 环状同心的统计
# 注意一定要倒序操作
# 200 = 200 - 150
# 150 = 150 - 100 
# ...
# 
res_pa_ls2 <- res_pa_ls
diff_idx_ls <-
  c(10, 50, 100, 150, 200) %>% 
  as.character() %>%
  rev() %>%
  {matrix(c(.[-length(.)], .[-1]), ncol = 2)} %>%
  apply(1, c, simplify = FALSE)

for(i in diff_idx_ls){
  first  <- i[1]
  second <- i[2]
  res_pa_ls2[[ first ]] <- res_pa_ls2[[ first ]] - res_pa_ls2[[ second ]] 
}

# count to freq
res_pa_ls2 <- 
  res_pa_ls2 %>%
  map(function(x){# browser()
    x_colsum <- colSums(x)
    xx <- sweep(x, 2, x_colsum, FUN = "/")
    if(any(x_colsum == 0)){
      xx[, x_colsum == 0] <- 0
    }
    xx
  })

# 每个样本中的每种Target的其他细胞的计数取中位数
res_pa_count_tidy <- 
  res_pa_ls2 %>%
  imap(function(x, distance){# browser()
    # 分组
    pa_groups <-
      with(
        all_data_meta_pos,
        str_c(Classification1, Group, Slide_ID, sep = '..')
      )
    # 分组内取中位数
    res <- apply(x, 1, function(y){tapply(y, INDEX = pa_groups, FUN = "median")})
    
    # data tidy
    res_df  <-
      res %>%
      as.data.frame() %>%
      rownames_to_column('group') %>%
      pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
      separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 
    res_df$distance <- distance
    
    # return
    res_df
  }) %>%
  bind_rows()


# plot 
# only plot "Epithelial cells" ....
res_pa_count_tidy_epi <-
  res_pa_count_tidy %>%
  mutate(distance = factor(distance, levels = sort(as.numeric(unique(distance))))) %>%
  dplyr::filter(Target == "Epithelial")

res_pa_count_tidy_epi_median <-
  res_pa_count_tidy_epi %>%
  group_by(Others, Group, distance) %>% 
  summarise(cell_props = median(cell_props))

p_proximity_series_epi <-
  res_pa_count_tidy_epi %>%
  ggplot(aes(x = distance, y = cell_props, color = Group)) +
  geom_boxplot( ) +
  ggsci::scale_color_aaas() +
  # geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = 0.5)) +
  geom_point(data = res_pa_count_tidy_epi_median, position = position_dodge(width = 0.75)) +
  geom_smooth(data = res_pa_count_tidy_epi_median, aes(group = Group), position = position_dodge(width = 0.75)) +
  facet_wrap(~ Others, scales = 'free_y') +
  ggpubr::theme_pubclean()


# export plot
pdf(
  file.path("result/proximity", "proximity analysis_series distance.pdf"),
  height = 14,
  width  = 18
)
p_proximity_series_epi %>% print()
dev.off()


# 3. Nearest neighbor analysis -------
# 由于最近邻只需要最近的那个细胞的信息
# 因此不同于proximity analysis，有以下变化：
# classification调为NULL
# searchtype 邻近度使用radius
# k不需要太大，10即可，最近邻分析只需要第二个细胞的信息
# radius适当调大100um，防止部分细胞比较离散，距离其他细胞较远
# returntype与cell_enough_ratio_cutoff也没有意义了

# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification NULL即可，返回每一个query细胞的细胞编号、距离、query和total细胞的id
# k 计算的最近邻的细胞数
# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325

# all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
# 
# all_data_meta_pos <- all_data_meta_pos %>%
#   as.data.frame() %>%       # 确保为数据框类型
#   remove_rownames()     # 强制清除现有行名
#   

res_nn_ls <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = NULL, # 若为NULL: 提前返回window的细胞近邻的list
    k = 10, # 'auto',
    searchtype = "radius",
    radius = 100, 
    mpp = 1, # !important, 1 for xenium
    cell_enough_ratio_cutoff = 20
  )

# res_nn_ls <- readRDS("data/res_nn_ls.rds")

res_nn_df <-
  res_nn_ls %>%
  map(function(x){ # browser()
    cell_idx  <- x$nn.idx[, 2]
    cell_idx[cell_idx == 0] <- NA # some cells have no nn cells
    cell_id   <- x$total_cell_id[cell_idx]
    cell_dists <- x$nn.dists[, 2]
    
    all_data_used <- all_data_meta_pos %>% column_to_rownames("cell_id")
    cell_anno <- all_data_used[cell_id, 'Classification1']
    
    res <-
      data.frame(
        Others = cell_anno,
        distance = cell_dists
      ) %>%
      cbind(all_data_used[x$query_cell_id, c('Classification1', 'Group', 'Slide_ID')]) %>%
      rownames_to_column("cell_id") %>%
      dplyr::rename(Target = Classification1) 
    
    res
  }) %>%
  bind_rows()

if(any(is.na(res_nn_df$Others) )){
  message(str_glue("{sum(is.na(res_nn_df$Others))} cells have no nearest neighbor cell"))
  res_nn_df <- drop_na(res_nn_df)
}

# plot
p_nn_ls <-
  res_nn_df %>%
  group_by(Target) %>%
  nest() %>%
  deframe() %>%
  imap(function(x, target){ # browser()
    x %>%
      ggplot(aes(x = distance, color = Group, fill = Group)) + 
      geom_density(alpha = 0.5) +
      facet_wrap(~ Others) +
      ggpubr::theme_pubclean() +
      labs(
        title = str_glue("Target cell: {target}"),
        x = "distance(um)"
      )
    
  })

# export plot
pdf(
  file.path("result/nearest neighbor", "nearest neighbor analysis.pdf"),
  height = 14,
  width  = 18
)
p_nn_ls %>% print()
dev.off()



#9、FYN、WSL、GYX治疗前后免疫-溶瘤病毒相关通路空间表达变化----
library(Seurat)
library(clusterProfiler)
library(ggplot2)
library(pheatmap)
#   治疗前         治疗后
#"WSL_23_18305","WSL_23_20895"###
#"GYX_24_03997","GYX_23_06802"###
#"CYS_23_03967","CYS_23_10358"
#"ZKJ_23_18277","ZKJ_23_19464"
#"ZBH_23_23753","ZBH_23_28671","ZBH_23_25688"
#"ZL_23_30769","ZL_23_32098"
#"HYY_24_01346","HYY_24_02802"
sample_group_list <- list(
  "CYS_23_03967" = "before",
  "GYX_23_06802" = "after",
  "CYS_23_10358" = "after",
  "FYN_23_11342" = "before",
  "ZKJ_23_18277" = "before",
  "WSL_23_18305" = "before",
  "GYX_24_03997" = "before",
  "ZKJ_23_19464" = "after",
  "WSL_23_20895" = "after",
  "ZBH_23_23753" = "before",
  "ZBH_23_25688" = "after",
  "ZBH_23_28671" = "after",
  "ZL_23_30769" = "before",
  "ZL_23_32098" = "after",
  "HYY_24_01346" = "before",
  "SQQ_23_33594" = "before",
  "HYY_24_02802" = "after",
  "CZY_24_06201" = "before"
)
# 提取 xenium 的 meta.data
meta_data <- xenium@meta.data
# 新增 group2 列，根据 sample 列匹配对应的值
meta_data$group2 <- sapply(meta_data$sample, function(x) sample_group_list[[x]])
# 将修改后的 meta.data 重新赋值回 xenium 对象
xenium@meta.data <- meta_data

xenium_plot <- xenium

#通路基因集----
pathway_gene_sets <- list(
  # 免疫反应类型相关通路
  "Type1_Immune_Response" = c("TBX21", "IFNG", "IL12B", "IL12RB1", "IL12RB2", "IL2", "TNF", 
                              "CXCR3", "CXCL9", "CXCL10", "CXCL11", "STAT1", "STAT4", "GZMB", 
                              "PRF1", "CD8A", "CD8B"),
  "Type2_Immune_Response" = c("GATA3", "IL4", "IL5", "IL13", "IL33", "IL1RL1", "CCR4", "CCR8", 
                              "STAT6", "PTGDR2", "IL10", "ARG1", "CCL17", "CCL22"),
  "Type17_Immune_Response" = c("RORC", "IL17A", "IL17F", "IL22", "IL23", "IL23R", "CCR6", 
                               "STAT3", "IL6", "IL1B", "CCL20", "IL21", "IL6ST"),
  "Regulatory_Immune_Response" = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IL10", "TGFB1", 
                                   "IKZF2", "ENTPD1", "ICOS", "TIGIT", "LAG3", "PDCD1"),
  "Tfh_Cell_Response" = c("BCL6", "CXCR5", "PDCD1", "ICOS", "SH2D1A", "IL21", "CD40LG", 
                          "BTLA", "MAF", "BATF"),
  "B_Cell_Mediated_Immunity" = c("CD19", "MS4A1", "CD79A", "CD79B", "IGHM", "IGHG1", "IGHG4", 
                                 "IGHA1", "IGHA2", "CD27", "TNFRSF13B", "TNFRSF13C", "TNFRSF17", 
                                 "SDC1"),
  "ILC1_Mediated_Immunity" = c("TBX21", "IFNG", "TNF", "ID2", "NFIL3", "EOMES", "KLRB1", 
                               "IL12RB1", "IL12RB2", "NCR1", "NCR2", "NCR3"),
  "ILC2_Mediated_Immunity" = c("GATA3", "IL5", "IL13", "RORA", "ID2", "NFIL3", "IL1RL1", 
                               "IL17RB", "KLRB1", "PTGDR2", "KIT", "IL4", "AREG"),
  "ILC3_Mediated_Immunity" = c("RORC", "IL17A", "IL17F", "IL22", "ID2", "NFIL3", "NCR1", 
                               "NCR2", "KIT", "IL23R", "IL1R1", "AHR", "LTB"),
  "Myeloid_Mediated_Immunity" = c("CD68", "CD14", "CD16", "CD11b", "CD11c", "TLR2", "TLR4", 
                                  "TLR7", "TLR8", "TLR9", "NLRP3", "PYCARD", "IL1B", "IL18", "MYD88"),
  "NK_Cell_Mediated_Cytotoxicity" = c("NCAM1", "FCGR3A", "KLRC1", "KLRD1", "NCR1", "GZMB", 
                                      "GNLY", "PRF1", "IFNG", "TNF", "KIR family genes"),
  "Complement_Mediated_Immunity" = c("C1QA", "C1QB", "C1QC", "C2", "C3", "C4A", "C4B", "C5", 
                                     "C6", "C7", "C8A", "C8B", "C8G", "C9", "CFB", "CFD", "CFH", 
                                     "CFI", "CR1", "CR2"),
  "Antiviral_Interferon_Response" = c("IFNA1", "IFNA2", "IFNA3", "IFNB1", "IFNAR1", "IFNAR2", 
                                      "IRF3", "IRF7", "IRF9", "STAT1", "STAT2", "ISG15", "MX1", 
                                      "MX2", "OAS1", "OAS2", "OAS3", "RSAD2", "IFIT1", "IFIT2", 
                                      "IFIT3", "DDX58", "IFIH1"),
  "Inflammatory_Response" = c("IL1A", "IL1B", "IL6", "TNF", "CXCL8", "CCL2", "ICAM1", "VCAM1", 
                              "SELE", "PTGS2", "IL1RN", "NFKB1", "RELA", "NFKBIA", "NLRP3", "CASP1"),
  "Antigen_Presentation" = c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", 
                             "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA", "TAP1", "TAP2", 
                             "PSMB8", "PSMB9", "B2M"),
  
  # 表观遗传调控相关通路
  "DNA_Methylation" = c("DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3", "MBD1", "MBD2", 
                        "MBD3", "MBD4"),
  "Histone_Modification" = c("EZH2", "SUZ12", "KMT2A", "KMT2B", "KMT2C", "KMT2D", "HDAC1", 
                             "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7", "HDAC8", 
                             "HDAC9", "HDAC10", "HDAC11", "KDM1A", "KDM5A", "KDM5B", "KDM5C", 
                             "KDM5D", "KDM6A", "KDM6B"),
  "Chromatin_Remodeling" = c("SMARCA4", "SMARCB1", "ARID1A", "ARID2", "CHD4", "CHD7", "BRD4"),
  "Non-coding_RNA_Regulation" = c("HOTAIR", "MALAT1", "NEAT1", "XIST", "H19", "PVT1"),
  
  # 肿瘤干细胞特征相关通路
  "Melanoma_Stem_Cell_Markers" = c("ABCB5", "CD271", "CD133", "ALDH1A1", "SOX2", "SOX10"),
  "Stemness_Signaling_Pathways" = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "JAG1", "JAG2", 
                                    "DLL1", "DLL3", "DLL4", "HES1", "HEY1"),
  "Pluripotency_Markers" = c("NANOG", "POU5F1", "SOX2", "KLF4", "MYC"),
  "Quiescence_Activation_Markers" = c("CXCR6", "LGR5", "ZEB1", "TWIST1", "SNAI1", "SNAI2"),
  
  # 恶性进展特征相关通路
  "EMT_Related" = c("CDH1", "CDH2", "VIM", "TWIST1", "SNAI1", "SNAI2", "ZEB1", "ZEB2"),
  "Invasion_Metastasis" = c("MMP2", "MMP9", "MMP14", "TIMP1", "TIMP2", "TIMP3", "TIMP4", 
                            "PLAU", "PLAUR", "LOX"),
  "Angiogenesis" = c("VEGFA", "VEGFB", "VEGFC", "FLT1", "KDR", "FLT4", "ANGPT1", "ANGPT2", "TEK"),
  "Neuroinvasion" = c("GAP43", "NCAM1", "BDNF", "NGF", "NTRK1", "NTRK2", "NTRK3", "SERPINE2", "L1CAM"),
  
  # 病毒-肿瘤相互作用相关通路
  "Viral_Recognition_Receptors" = c("TLR3", "TLR7", "TLR8", "TLR9", "DDX58", "IFIH1", "TMEM173"),
  "Intracellular_Viral_Sensing" = c("IRF3", "IRF7", "MAVS", "TBKBP1", "TANK", "TBK1", "MYD88", "TRIF"),
  "Oncolytic_Viral_Replication" = c("EIF2AK2", "OAS1", "OAS2", "OAS3", "RNASEL", "ADAR", "APOBEC3A", 
                                    "APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H"),
  "Viral_Immune_Evasion" = c("USP18", "SOCS1", "SOCS3", "ISG15", "TRIM25", "HERC5", "UBE2L6"),
  
  # 治疗耐药机制相关通路
  "Immune_Checkpoints" = c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "HAVCR2", "LAG3", "TIGIT"),
  "Immunosuppressive_Cytokines" = c("TGFB1", "TGFB2", "TGFB3", "IL10", "VEGFA", "CSF1", "IDO1", "PTGS2"),
  "Drug_Resistance_TFs" = c("STAT3", "NFKB1", "RELA", "FOXO3", "HIF1A", "MYC", "ETV1"),
  "Drug_Efflux_Transporters" = c("ABCB1", "ABCC1", "ABCG2", "ABCB5"),
  
  # 组织微环境相互作用相关通路
  "Fibroblast_Activation" = c("ACTA2", "FAP", "PDGFRB", "TGFB1", "POSTN", "SPARC", "WNT5A"),
  "Immune_Cell_Recruitment" = c("CCL2", "CCL5", "CXCL9", "CXCL10", "CXCL11", "CXCL12", "CCR2", 
                                "CCR5", "CXCR3", "CXCR4"),
  "ECM_Remodeling" = c("COL1A1", "COL4A1", "LAMB1", "LAMC1", "FN1", "VCAN", "DCN", "LUM"),
  "Tumor_Stroma_Signaling" = c("WNT", "SHH", "IHH", "DHH", "PTCH1", "SMO", "GLI1", "GLI2", "GLI3", 
                               "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"),
  
  # 细胞衰老特征相关通路
  "Cell_Cycle_Regulation" = c("CDKN1A", "CDKN2A", "RB1", "TP53", "MDM2", "E2F1", "CCND1", "CDK4", "CDK6"),
  "SASP" = c("IL1A", "IL1B", "IL6", "IL8", "MMP1", "MMP3", "MMP10", "CXCL1", "CXCL2"),
  "DNA_Damage_Response" = c("H2AFX", "ATM", "ATR", "CHEK1", "CHEK2", "BRCA1", "BRCA2", "RAD51"),
  "Telomerase_Activity" = c("TERT", "TERC", "TERF1", "TERF2", "POT1", "TINF2", "ACD"),
  
  # 细胞应激反应相关通路
  "ER_Stress" = c("XBP1", "ERN1", "ATF6", "EIF2AK3", "DDIT3", "HSPA5", "ATF4"),
  "Heat_Shock_Response" = c("HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1", "HSPB1", "HSF1"),
  "Oxidative_Stress" = c("NFE2L2", "KEAP1", "HMOX1", "NQO1", "GPX1", "GPX2", "GPX3", "GPX4", 
                         "GCLC", "GCLM", "GSR"),
  "Autophagy_Stress" = c("BECN1", "MAP1LC3B", "ATG5", "ATG7", "SQSTM1", "ULK1", "AMBRA1", "GABARAP")
)

# library(GSVA)
# # 提取表达矩阵，假设是 counts 矩阵（也可用归一化后的矩阵，根据需求选择）
# expr_matrix <- GetAssayData(xenium_plot, assay = "Xenium", slot = "counts")  
# # 进行 GSVA 分析
# gsva_scores <- gsva(expr_matrix, pathway_gene_sets, method = "gsva")  
# # 将 GSVA 结果转换为数据框并添加到 Seurat 对象元数据中（或者直接用于热图）
# gsva_scores_df <- as.data.frame(t(gsva_scores))
# xenium_plot@meta.data <- cbind(xenium_plot@meta.data, gsva_scores_df)
# # 后续整理成矩阵绘制热图的步骤类似，把提取的列换成 GSVA 打分对应的列即可

all_genes <- rownames(xenium)
# 筛选存在于Xenium数据的通路基因集
filtered_pathways <- lapply(pathway_gene_sets, function(genes) {
  intersect(genes, all_genes)
})
# 移除空基因集
filtered_pathways <- filtered_pathways[sapply(filtered_pathways, length) > 0]
pathway_names <- names(filtered_pathways)

xenium_plot <- AddModuleScore(
  object = xenium_plot,
  features = filtered_pathways,
  name = "Pathway_",
  ctrl = 50,  # 使用100个随机基因作为背景校正[1](@ref)
  seed = 42
)

# 重命名打分列（修复自动编号问题）
colnames(xenium_plot@meta.data)[grep("Pathway_", colnames(xenium_plot@meta.data))] <- pathway_names

# 提取metadata并筛选before组----
library(dplyr)
meta_before <- xenium_plot@meta.data %>% 
  filter(group2 == "before") %>%
  dplyr::select(sample, all_of(pathway_names))

# 计算每个样本的通路平均活性
pathway_means <- meta_before %>%
  group_by(sample) %>%
  summarise(across(all_of(pathway_names), mean, na.rm = TRUE)) %>%
  tibble::column_to_rownames("sample") %>%
  as.matrix()

# 2. 创建颜色映射（蓝-白-红表示低-中-高）
my_palette <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)

# 3. 绘制热图
library(pheatmap)
H <- pheatmap(
  mat = t(as.matrix(pathway_means)),
  color = my_palette,
  border_color = NA,
  cluster_rows = F,  # 是否对行聚类
  cluster_cols = TRUE,  # 是否对列聚类
  clustering_method = "ward.D2",  # 聚类方法
  show_rownames = TRUE,  # 显示行名
  show_colnames = TRUE,  # 显示列名
  fontsize_row = 9,      # 行名字体大小
  fontsize_col = 10,     # 列名字体大小
  angle_col = 45,        # 列名角度
  main = "Pathway Activity",  # 标题
  legend = TRUE,         # 显示图例
  scale = "column",       # 列标准化（使不同通路可比）
  #height = 15,
  #width = 5,
  cellwidth = 18,         # 单元格宽度
  cellheight = 10         # 单元格高度
)

ggsave(H, filename = "before_PathwayActivity_heatmap.pdf", height = 10, width = 18,)



# 提取metadata并筛选after组----
library(dplyr)
meta_after <- xenium_plot@meta.data %>% 
  filter(group2 == "after") %>%
  dplyr::select(sample, all_of(pathway_names))

# 计算每个样本的通路平均活性
pathway_means <- meta_after %>%
  group_by(sample) %>%
  summarise(across(all_of(pathway_names), mean, na.rm = TRUE)) %>%
  tibble::column_to_rownames("sample") %>%
  as.matrix()

# 2. 创建颜色映射（蓝-白-红表示低-中-高）
my_palette <- colorRampPalette(c("#4575b4", "white", "#d73027"))(100)

# 3. 绘制热图
library(pheatmap)
H <- pheatmap(
  mat = t(as.matrix(pathway_means)),
  color = my_palette,
  border_color = NA,
  cluster_rows = F,  # 是否对行聚类
  cluster_cols = TRUE,  # 是否对列聚类
  clustering_method = "ward.D2",  # 聚类方法
  show_rownames = TRUE,  # 显示行名
  show_colnames = TRUE,  # 显示列名
  fontsize_row = 9,      # 行名字体大小
  fontsize_col = 10,     # 列名字体大小
  angle_col = 45,        # 列名角度
  main = "Pathway Activity",  # 标题
  legend = TRUE,         # 显示图例
  scale = "column",       # 列标准化（使不同通路可比）
  #height = 15,
  #width = 5,
  cellwidth = 18,         # 单元格宽度
  cellheight = 10         # 单元格高度
)

ggsave(H, filename = "after_PathwayActivity_heatmap.pdf", height = 10, width = 18,)






#10、治疗有效患者(WSL/GYX)治疗前（基线）和治疗无效患者(基线)细胞比例变化柱状图、箱线图、关键细胞类型映射----



#11、巨噬评分（不是整张片子）----
Macrophage_out_dir <- "result/Macrophage" 
dir.create(Macrophage_out_dir)

Macrophage <- subset(xenium, subset = first_type %in% c("Mono-Mac"))

## PCA ---
Macrophage <- FindVariableFeatures(Macrophage)
Macrophage <- ScaleData(Macrophage, vars.to.regress = 'nCount_RNA')
Macrophage <- RunPCA(Macrophage)

### pca inspector ---
pdf(file.path(Macrophage_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage, ndims = 50)
dev.off()

pdf(file.path(Macrophage_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 1:5)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 6:10)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 11:15)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage, dims = 16:20)
dev.off()

## dimension reduction ---
dims <- 20 # 18

library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)

Macrophage <- RunUMAP(Macrophage, dims = 1:dims,return.model = TRUE)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} sample.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE) + NoLegend()
dev.off()

Macrophage <- Macrophage %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    # theta = 4, # default 2, bigger value -- more diversity per cluster
    # sigma = 0.5,
    # lambda = 0.5, # default, smaller value -- more integrated
    reduction.save = "harmony"
  )

Macrophage <- RunUMAP(
  Macrophage, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster ---
Macrophage  <-
  FindNeighbors(
    Macrophage ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Macrophage <-
  FindClusters(
    Macrophage,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.1.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

Idents(Macrophage) <- "RNA_snn_res.0.6"

Macrophage <- readRDS("RCTD/result/Macrophage/Macrophage.rds")

#巨噬基因集----
macrophage_gene_sets <- c(
  # ==== M1型极化 (Classical activation)
  'M1' = "CD86",      # M1共刺激分子，促进T细胞活化[1,6](@ref)
  'M1' = "TNF",       # 促炎因子，诱导细胞坏死并激活内皮细胞[7,8](@ref)
  'M1' = "IL1B",      # 炎症小体下游因子，驱动急性期反应[1,7](@ref)
  'M1' = "NOS2",      # 一氧化氮合酶，增强杀菌活性（M1核心标志）[1,6](@ref)
  'M1' = "CXCL9",     # IFN-γ诱导型T细胞趋化因子[1,5](@ref)
  'M1' = "CXCL10",    # IP-10趋化因子，招募Th1细胞[1,5](@ref)
  
  # ==== M2型极化 (Alternative activation)
  'M2' = "ARG1",      # 精氨酸酶，抑制NO产生并促进组织修复[1,6](@ref)
  'M2' = "IL10",      # 抗炎因子，抑制巨噬细胞活化[1,7](@ref)
  'M2' = "TGFB1",     # 免疫抑制因子，促进纤维化及Treg分化[1,8](@ref)
  'M2' = "IL4R",      # IL-4受体，介导M2极化信号转导[1,6](@ref)
  'M2' = "MRC1",      # 甘露糖受体(CD206)，内吞病原体（M2典型标记）[1,5](@ref)
  
  # ==== 肿瘤相关巨噬细胞 (TAM signature)
  'TAM' = "TREM2",    # 脂质感应受体，促进免疫抑制性微环境[3,5](@ref)
  'TAM' = "SPP1",     # 骨桥蛋白，介导肿瘤转移与血管生成[3,5](@ref)
  'TAM' = "C1QA",     # 补体成分，参与肿瘤免疫逃逸[3,6](@ref)
  'TAM' = "APOE",     # 载脂蛋白，调控胆固醇代谢（TAM代谢重编程）[3,5](@ref)
  'TAM' = "FOLR2",    # 叶酸受体β，组织驻留巨噬细胞标记[5,6](@ref)
  
  # ==== 吞噬功能 (Phagocytosis)
  'Pha' = "MARCO",      # 清道夫受体，结合细菌表面多糖[6,8](@ref)
  'Pha' = "MSR1",       # 巨噬细胞清道夫受体1，内化氧化LDL[6,8](@ref)
  'Pha' = "CD36",       # 脂肪酸转运蛋白，介导凋亡细胞识别[4,6](@ref)
  'Pha' = "MERTK",      # 酪氨酸激酶受体，调控吞噬体成熟[4,6](@ref)
  
  # ==== 代谢重编程 (Metabolism)
  'Met' = "IDO1",       # 色氨酸分解酶，抑制T细胞功能[4,7](@ref)
  'Met' = "LDHA",       # 乳酸脱氢酶，促进糖酵解（Warburg效应）[4,9](@ref)
  'Met' = "HK2",        # 己糖激酶2，糖酵解限速酶[4,9](@ref)
  'Met' = "SLC2A1"      # 葡萄糖转运蛋白GLUT1（标准命名），增强葡萄糖摄取[4,9](@ref)
)

valid_macrophage_gene_sets <- macrophage_gene_sets[macrophage_gene_sets %in% rownames(xenium)]
cat("有效基因:", length(valid_macrophage_gene_sets), "/", length(macrophage_gene_sets))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_genes <- setdiff(macrophage_gene_sets, valid_macrophage_gene_sets)

###基因集---
modules <- names(valid_macrophage_gene_sets) |> unique()
modules

Macrophage <- AddModuleScore(object = Macrophage,
                         features = split(valid_macrophage_gene_sets, modules),
                         name = "VirusModule_",
                         ctrl = 50,              # 对照基因数（建议50-100）
                         seed = 2025)             # 确保可重复性

# virus_module_cols <- c("IFN_I", "IFN_II", "PRR", "ASP", "CS_D", "IBA")
# # 批量删除这些列
# for (col in virus_module_cols) {
#   xenium@meta.data[[col]] <- NULL
# }


# 定义新旧列名映射关系
new_names <- c("M1","M2","TAM","Pha","Met")
old_names <- paste0("VirusModule_", 1:5)

# 直接修改metadata列名
colnames(Macrophage@meta.data)[match(old_names, colnames(Macrophage@meta.data))] <- new_names

xenium@meta.data %>% head()

# 提取巨噬细胞评分数据（确保细胞ID对齐）
macrophage_scores <- Macrophage@meta.data[, c("cell_id", "M1", "M2", "TAM", "Pha", "Met")]  # 假设细胞ID列名为"cell_id"

# 添加到Xenium对象的metadata
xenium <- AddMetaData(
  object = xenium,
  metadata = macrophage_scores,
  col.name = c("M1", "M2", "TAM", "Pha", "Met")  # 指定新列名[6](@ref)
)

测试
Idents(xenium) <- "first_type"
ImageDimPlot(xenium, fov = "CZY_24_06201",nmols = 20000,cols ='#FFD700',cells = WhichCells(xenium, idents = "Mono-Mac" ))


###病毒评分空间featureplt----
MacModuleSegImagefeature_dir <- 'result/MacModuleSegImagefeature_18sample' 
dir.create(MacModuleSegImagefeature_dir)

module_names <- new_names

for (sample_id in unique(xenium@meta.data$sample)) {
  # 智能匹配样本ID
  j <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  
  if (length(j) == 0) {
    warning("跳过未匹配样本: ", sample_id)
    next
  }
  
  for (module in module_names) {
    tryCatch({
      # 空值检查
      if (!is.null(xenium@images[[j]])) {
        DefaultBoundary(xenium@images[[j]]) <- "segmentation"
      }
      
      p <- ImageFeaturePlot(
        object = xenium,
        fov = j,
        nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
        border.color = "white",  # 细胞边界的颜色（深蓝色）
        border.size = 0.01,    # 细胞边界的线宽
        features = module,
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('blue', "red"),
        size = 0.8,
        dark.background = TRUE
      )+
        ggtitle(paste(module, "Score")) +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right",
          legend.title = element_text(size = 10)
        )
      
      ggsave(
        filename = file.path(MacModuleSegImagefeature_dir, paste0(sample_id, "_", module,"_spatial.png")),
        plot = p,
        width = 12,
        height = 12,
        dpi = 300
      )
      
      message("成功保存: ", sample_id, "_", module)
      
    }, error = function(e) {
      warning("失败[", sample_id, ":", "MacModule_",module,"]: ", e$message)
    })
  }
}

###病毒反应评分箱线图（带检验）----
library(Seurat)
library(ggpubr)
library(tidyverse)

MacModuleBox_dir <- 'result/MacModuleBox' 
dir.create(MacModuleBox_dir)

# saveRDS(xenium,paste0('xenium_coor_group',"_final.rds"))
# xenium <- readRDS("xenium_coor_group_final.rds")

###循环箱线图（按group分）----
#按照new_names <- c("IFN_I", "IFN_II", "PRR", "ASP", "CS_D", "IBA")绘制"Effective"和"Ineffective"两组的箱线图，加入统计学检验

# 定义新旧列名映射关系
new_names <- c("M1","M2","TAM","Pha","Met")

# 1. 提取模块评分数据和分组信息
# 从Seurat对象metadata中提取指定列（模块评分和分组信息）
plot_data <- xenium@meta.data %>%
  # 选择新命名的模块列（new_names）和分组列（group1）
  select(all_of(new_names), group1) %>%
  # 将宽数据转为长数据：多个模块列合并为两列（Module列存放模块名，Score列存放评分值）
  pivot_longer(cols = all_of(new_names), 
               names_to = "Module", 
               values_to = "Score") %>%
  # 将Module列转换为因子，固定模块顺序（避免绘图时按字母排序）
  mutate(Module = factor(Module, levels = new_names))

# 循环绘制每个模块
for (module in new_names) {
  # 1. 提取当前模块数据
  module_data <- plot_data %>%
    filter(Module == module)
  
  # 2. 动态计算y轴上限（用于统计标注位置）
  y_max <- max(module_data$Score, na.rm = TRUE) * 1.15
  
  # 3. 绘制单个模块箱线图
  p <- ggplot(module_data, aes(x = group1, y = Score, fill = group1)) +
    # 箱线图层
    geom_boxplot(
      width = 0.5,
      alpha = 0.85,
      outlier.shape = 19,
      outlier.size = 0.8,
      outlier.alpha = 0.6
    ) +
    # # 样本点层（显示每个样本）
    # geom_jitter(
    #   width = 0.15,
    #   size = 2,
    #   alpha = 0.5,
    #   color = "black"
    # ) +
    # # 均值标记（菱形符号）
    # stat_summary(
    #   fun = mean,
    #   geom = "point",
  #   shape = 18,
  #   size = 5,
  #   color = "gold"
  # ) +
  # 统计检验标注
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif", # 改为星号显示（*:p<0.05; **:p<0.01; ***:p<0.001）
    label.x = 1.5,
    size = 5,
    label.y = y_max,
    vjust = -0.5
  ) +
    # 配色方案
    scale_fill_manual(
      values = c("Effective" = "#4DBBD5", "Ineffective" = "#E64B35")
    ) +
    # 标签和标题
    labs(
      title = paste("Module:", module),
      #subtitle = "Distribution by Treatment Response",
      x = "Group",
      y = "Module Score",
      fill = "Group"
    ) +
    # 主题设置
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      #plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转45度并右对齐
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey70")
    )
  
  # 4. 保存图像
  ggsave(
    filename = file.path(MacModuleBox_dir, paste0('group_', module, "_boxplot.png")),
    plot = p,
    width = 6,
    height = 7,
    dpi = 300
  )
  
  message("已保存: ", module, " 的箱线图")
}

###循环箱线图（按sample分）----
#按照new_names <- c("IFN_I", "IFN_II", "PRR", "ASP", "CS_D", "IBA")绘制"Effective"和"Ineffective"两组的箱线图，加入统计学检验

# 1. 提取模块评分数据和分组信息
# 从Seurat对象metadata中提取指定列（模块评分和分组信息）
plot_data <- xenium@meta.data %>%
  # 选择新命名的模块列（new_names）和分组列（group1）
  select(all_of(new_names), sample) %>%
  # 将宽数据转为长数据：多个模块列合并为两列（Module列存放模块名，Score列存放评分值）
  pivot_longer(cols = all_of(new_names), 
               names_to = "Module", 
               values_to = "Score") %>%
  # 将Module列转换为因子，固定模块顺序（避免绘图时按字母排序）
  mutate(Module = factor(Module, levels = new_names))

# 循环绘制每个模块
for (module in new_names) {
  # 1. 提取当前模块数据
  module_data <- plot_data %>%
    filter(Module == module)
  
  # 2. 动态计算y轴上限（用于统计标注位置）
  y_max <- max(module_data$Score, na.rm = TRUE) * 1.15
  
  # 3. 绘制单个模块箱线图
  p <- ggplot(module_data, aes(x = sample, y = Score, fill = sample)) +
    # 箱线图层
    geom_boxplot(
      width = 0.5,
      alpha = 0.85,
      outlier.shape = 19,
      outlier.size = 0.8,
      outlier.alpha = 0.6
    ) +
    # # 样本点层（显示每个样本）
    # geom_jitter(
    #   width = 0.15,
    #   size = 2,
    #   alpha = 0.5,
    #   color = "black"
    # ) +
    # # 均值标记（菱形符号）
    # stat_summary(
    #   fun = mean,
    #   geom = "point",
  #   shape = 18,
  #   size = 5,
  #   color = "gold"
  # ) +
  # 统计检验标注
  # stat_compare_means(
  #   method = "wilcox.test",
  #   label = "p.signif", # 改为星号显示（*:p<0.05; **:p<0.01; ***:p<0.001）
  #   label.x = 1.5,
  #   size = 5,
  #   label.y = y_max,
  #   vjust = -0.5
  # ) +
  # 配色方案
  scale_fill_manual(
    values = plan
  ) +
    # 标签和标题
    labs(
      title = paste("Module:", module),
      #subtitle = "Distribution by Treatment Response",
      x = "sample",
      y = "Module Score",
      #fill = "sample"
    ) +
    # theme(legend.position = "none")
    # 主题设置
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      #plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转45度并右对齐
      axis.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey70")
    )
  
  # 4. 保存图像
  ggsave(
    filename = file.path(MacModuleBox_dir, paste0('sample_',module, "_boxplot.png")),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  message("已保存: ", module, " 的箱线图")
}




#保存数据----
#saveRDS(xenium,paste0('xenium_coor',"_final.rds"))
xenium <- readRDS("RCTD/xenium_coor_final.rds")
