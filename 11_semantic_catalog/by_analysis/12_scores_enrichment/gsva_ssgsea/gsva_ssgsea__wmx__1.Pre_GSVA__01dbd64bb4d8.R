# 提取rds中的癌症样品
# 进行Psedobulk化，计算各个MP的特征的GSVA评分
# 进行高低组样品特征的GSVA评分
# 计算MP与高低组的GSVA评分的相关性

library(Seurat)
library(dplyr)
library(Matrix)
library(edgeR)


# 1.读入的数据，需要按照group提取LARC组
# 2.按照sample列分组，可命名为Immune_group：
# HIGH（7人）:DYD0703C（2545441），GYH0708C（2546436），LYF0725C（2550992），JL0728（2551432），QHX0729C（2551787），GC0812C（2545456）， WHL0717（2548977）
# LOW（9人）: YBN0710C（2547142），TXL0725C（2551168），BCR0807C（2554268），HJM0808C（2554603），WZQ0811C（2554936），ZYH0808C（2554439），ZBZ0808C（2554443），NX0708C（2546273），WHT0811C（2554933）
# 3.进行Pseudobulk话，后续使用GSVA分析

# =========================
# 输入输出
# =========================
rds_file <- "/mnt/NexData_file/home/sh/Projects/AE2511046228-02/注释数据/Main_mapping.rds"

outdir <- "/mnt/NexData_file/home/wmx/Project/AE2511046228-02/shouhou15/Pseudobulk_GSVA/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 1. 读入 Seurat 对象
# =========================
obj <- readRDS(rds_file)

DefaultAssay(obj) <- "RNA"

# 如果是 Seurat v5，建议合并 layers
if ("Assay5" %in% class(obj[["RNA"]])) {
  obj <- JoinLayers(obj, assay = "RNA")
}

# =========================
# 2. 提取 LARC 组
# =========================
if (!"group" %in% colnames(obj@meta.data)) {
  stop("metadata 中没有 group 列")
}

obj <- subset(obj, subset = group == "LARC")

cat("LARC cells:", ncol(obj), "\n")

# =========================
# 3. 定义 HIGH / LOW 样本
# =========================

HIGH_samples <- c(
  "DYD0703C",
  "GYH0708C",
  "LYF0725C",
  "JL0728C", # JL0728
  "QHX0729C",
  "GC0812C",
  "WHL0717C" # WHL0717
)

LOW_samples <- c(
  "YBJ0710C", # YBN0710C
  "TXL0725C",
  "BCR0807C",
  "HJM0808C",
  "WZQ0811C",
  "ZYH0808C",
  "ZBZ0808C",
  "NX0708C",
  "WHT0811C"
)



if (!"sample" %in% colnames(obj@meta.data)) {
  stop("metadata 中没有 sample 列")
}

obj$sample <- as.character(obj$sample)

obj$Immune_group <- dplyr::case_when(
  obj$sample %in% HIGH_samples ~ "HIGH",
  obj$sample %in% LOW_samples  ~ "LOW",
  TRUE ~ NA_character_
)

cat("Immune_group 分布:\n")
print(table(obj$Immune_group, useNA = "ifany"))

# 只保留成功匹配 HIGH / LOW 的样本
# obj <- subset(obj, subset = !is.na(Immune_group))

obj$Immune_group <- factor(obj$Immune_group, levels = c("LOW", "HIGH"))

cat("过滤后 cells:", ncol(obj), "\n")
print(table(obj$sample, obj$Immune_group))

# =========================
# 4. 保存整理后的 metadata
# =========================
write.table(
  obj@meta.data,
  file = file.path(outdir, "LARC_metadata_with_Immune_group.xls"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# =========================
# 5. Pseudobulk by sample
# =========================

data <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
)

sample_vec <- obj$sample

# gene x sample 的 pseudobulk counts
pb_data <- rowsum(
  t(as.matrix(data)),
  group = sample_vec
)

pb_data <- t(pb_data)

cat("Pseudobulk matrix:", nrow(pb_data), "genes x", ncol(pb_data), "samples\n")

# =========================
# 6. 构建 pseudobulk metadata
# =========================
pb_meta <- obj@meta.data %>%
  dplyr::select(sample, group, Immune_group) %>%
  distinct() %>%
  arrange(sample)

pb_meta <- pb_meta[match(colnames(pb_data), pb_meta$sample), ]

rownames(pb_meta) <- pb_meta$sample

# 检查顺序
stopifnot(all(colnames(pb_data) == rownames(pb_meta)))

# =========================
# 7. 输出 pseudobulk counts 和 metadata
# =========================
write.csv(
  pb_data,
  file = file.path(outdir, "LARC_pseudobulk_data_by_sample.csv")
)

write.table(
  pb_meta,
  file = file.path(outdir, "LARC_pseudobulk_metadata.xls"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# =========================
# 9. 保存对象
# =========================
saveRDS(
  obj,
  file = file.path(outdir, "LARC_seurat_with_Immune_group.rds")
)

cat("Done! Results saved to:", outdir, "\n")
