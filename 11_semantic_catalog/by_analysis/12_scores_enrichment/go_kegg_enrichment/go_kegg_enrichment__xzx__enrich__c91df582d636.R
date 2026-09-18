library(grid)
library(ggplot2)
library(gridExtra)
library(clusterProfiler)
#library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(org.Mm.eg.db)
library(GOplot)

setwd("/mnt/NexData_file/home/xzx/project/FWSC20252792/con_tst/")

sc <- readRDS('/mnt/NexData_file/home/xzx/project/FWSC20252792/con_tst/con_TST.rds')
sc <- NormalizeData(sc)
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

Idents(sc) <- sc$group
markers <- FindAllMarkers(object = sc, only.pos = TRUE)
table(markers$cluster)

library(dplyr)

gene = markers[markers$p_val <0.05,]$gene

entrezIDs = bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb= "org.Mm.eg.db", drop = TRUE)
#> Warning in bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb =
#> "org.Hs.eg.db", : 7.1% of input gene IDs are fail to map...
gene<- entrezIDs$ENTREZID

marker_new = markers[markers$gene %in% entrezIDs$SYMBOL,]
marker_new = marker_new[!duplicated(marker_new$gene),]
identical(marker_new$gene, entrezIDs$SYMBOL)
#> [1] FALSE
p = identical(marker_new$gene, entrezIDs$SYMBOL);p
#> [1] FALSE
if(!p) entrezIDs = entrezIDs[match(marker_new$gene,entrezIDs$SYMBOL),]
marker_new$ID = entrezIDs$ENTREZID


## GO
gcSample=split(marker_new$ID, marker_new$cluster) 
###参数可以更改，看看?compareCluster
#One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
xx <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Mm.eg.db",
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05
)

p <- dotplot(xx)
p

# 提取富集结果信息，这里以获取富集条目的描述、p值、各分组基因比率等为例
enrich_res <- as.data.frame(xx@compareClusterResult)

# 先把 GeneRatio 从 "5/100" 转换为数值（5/100 -> 0.05）
enrich_df <- enrich_res %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    if (grepl("/", x)) {
      eval(parse(text = x))   # 把 "5/100" 变成 0.05
    } else {
      as.numeric(x)           # 如果已经是数值
    }
  }))

enrich_df <- enrich_df[enrich_df$p.adjust < 0.01,]
# 宽格式
enrich_wide <- pivot_wider(enrich_df,
                           id_cols = "Description",
                           names_from = "Cluster",
                           values_from = "GeneRatio")

# 把 NA 替换成 0（只对数值列）
enrich_wide <- enrich_wide %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

library(pheatmap)
# 5. 设置 Description 为行名，准备矩阵
mat <- as.data.frame(enrich_wide)
rownames(mat) <- mat$Description
mat <- mat[ , -1, drop = FALSE]   # 去掉 Description 列，保留数值

pheatmap(mat,
         scale = "none",              # 行标准化，方便比较
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         fontsize_row = 6,           # 通路名字体调小
         fontsize_col = 8)


################################################################################
# 1. GeneRatio 转数值
enrich_res <- enrich_res %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    if (grepl("/", x)) eval(parse(text = x)) else as.numeric(x)
  }))

# 2. 按 Cluster 分组，取每组 p.adjust 最小的前20条
top_terms <- enrich_res %>%
  group_by(Cluster) %>%
  arrange(p.adjust) %>%
  slice_head(n = 20) %>%
  ungroup()

# 3. 得到最终需要保留的通路集合
keep_desc <- unique(top_terms$Description)

# 4. 构建宽表（GeneRatio 或 -log10(p.adjust) 都可以，这里用 GeneRatio）
enrich_wide <- enrich_res %>%
  filter(Description %in% keep_desc) %>%
  dplyr::select(Description, Cluster, GeneRatio) %>%
  pivot_wider(id_cols = "Description",
              names_from = "Cluster",
              values_from = "GeneRatio") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

# 5. 准备矩阵
mat <- as.data.frame(enrich_wide)
rownames(mat) <- mat$Description
mat <- mat[ , -1, drop = FALSE]

# 6. 画热图
pheatmap(mat,
         scale = "none",
         color = colorRampPalette(c("white","orange","red"))(50),
         fontsize_row = 6,
         fontsize_col = 8,
         clustering_method = "complete")

#################################################################################
plot_enrich_heatmap <- function(enrich_res, top_n = 20, value_col = c("GeneRatio", "logFDR"),legend_title = NULL) {
  value_col <- match.arg(value_col)
  
  # 转换 GeneRatio
  enrich_res <- enrich_res %>%
    mutate(GeneRatio = sapply(GeneRatio, function(x) {
      if (grepl("/", x)) eval(parse(text = x)) else as.numeric(x)
    }),
    logFDR = -log10(p.adjust))
  
  # 每个 Cluster 取前 top_n
  top_terms <- enrich_res %>%
    group_by(Cluster) %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  keep_desc <- unique(top_terms$Description)
  
  # 构建宽表
  enrich_wide <- enrich_res %>%
    filter(Description %in% keep_desc) %>%
    dplyr::select(Description, Cluster, all_of(value_col)) %>%
    pivot_wider(id_cols = "Description",
                names_from = "Cluster",
                values_from = all_of(value_col)) %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
  
  # 准备矩阵
  mat <- as.data.frame(enrich_wide)
  rownames(mat) <- mat$Description
  mat <- mat[ , -1, drop = FALSE]
  # 默认 legend 标题
  if (is.null(legend_title)) {
    legend_title <- ifelse(value_col == "GeneRatio", "Gene Ratio", "-log10(FDR)")
  }
  # 绘制热图
  pheatmap(mat,
           scale = "none",
           color = colorRampPalette(c("white","orange","red"))(50),
           fontsize_row = 14,
           fontsize_col = 8,
           cluster_cols = F,
           legend = TRUE,
           legend_labels = legend_title,
           clustering_method = "complete")
}

# 
# 用 GeneRatio
plot_enrich_heatmap(as.data.frame(xx@compareClusterResult), top_n = 20, value_col = "GeneRatio",
                    legend_title = "")

# 用 -log10(p.adjust)
p <- plot_enrich_heatmap(as.data.frame(xx@compareClusterResult), 
                    top_n = 20, 
                    value_col = "logFDR")
p
ggsave(plot = p,filename = "enrichGO_heatmap.png",
       width = 9,height = 12,dpi = 300,bg = "white")
ggsave(plot = p,filename = "enrichGO_heatmap.pdf",
       width = 9,height = 12,dpi = 300,bg = "white")

write.csv(enrich_res,"enrich_GO.csv")
write.csv(markers,"DEGs_group.csv")
