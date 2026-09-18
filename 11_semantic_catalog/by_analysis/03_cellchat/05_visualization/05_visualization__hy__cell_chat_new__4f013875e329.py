# cellchat_celltype_specific_diff_net_flexible.R
# 功能：为每个细胞类型单独绘制差异通讯网络（数量/强度），画布随细胞类型数量自动扩大
# 特点：健壮的细胞类型提取，支持手动回退

rm(list = ls())

library(CellChat)
library(ggplot2)
library(scales)

# ------------------------------
# 参数设置（与原脚本一致）
# ------------------------------
options(stringsAsFactors = FALSE)

outdir <- "/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/core/make_band/cellchat/in"
setwd(outdir)

load("/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/core/make_band/cellchat/in/cellchat_object.RData")

group <- c("pCR", "non-pCR")
compare.list <- c(1, 2)

# ------------------------------
# 细胞类型名称提取（多级回退）
# ------------------------------
idents <- NULL

# 回退1：直接从 net$count 矩阵行名获取（最可靠）
if (!is.null(cellchat@net$count) && length(cellchat@net$count) >= compare.list[2]) {
  tmp_mat <- cellchat@net$count[[compare.list[1]]]
  if (!is.null(tmp_mat) && nrow(tmp_mat) > 0) {
    idents <- rownames(tmp_mat)
  }
}

# 回退2：尝试 cellchat@idents（常见于旧版 CellChat 或单个对象）
if (is.null(idents)) {
  if (exists("idents", envir = as.environment(cellchat))) {
    idents <- levels(cellchat@idents)
  }
}

# 回退3：尝试 meta 数据框中的 "ident" 或 "labels" 列
if (is.null(idents)) {
  meta <- cellchat@meta
  if ("ident" %in% colnames(meta)) {
    idents <- unique(as.character(meta$ident))
  } else if ("labels" %in% colnames(meta)) {
    idents <- unique(as.character(meta$labels))
  }
}

# 回退4：尝试 cellchat@active.ident
if (is.null(idents)) {
  try({
    active.ident <- cellchat@active.ident
    if (is.factor(active.ident)) idents <- levels(active.ident)
    else if (is.character(active.ident)) idents <- unique(active.ident)
  }, silent = TRUE)
}

# 如果仍为空，输出对象结构并提示手动输入
if (is.null(idents) || length(idents) == 0) {
  message("Error: Unable to automatically extract cell types.")
  message("Here is a summary of the cellchat object structure:")
  print(str(cellchat, max.level = 2))
  
  # 手动输入提示（取消注释并填入你的细胞类型向量）
  # idents <- c("CD8+ T", "CD4+ T", "B cell", "Macrophage", ...)  
  stop("Please provide cell type names manually by uncommenting the 'idents <- ...' line above.")
}

message("Found cell types (", length(idents), "): ", paste(idents, collapse = ", "))

# ------------------------------
# 获取颜色（自动或手动）
# ------------------------------
color.use <- NULL

# 尝试从 meta$color 获取
if (!is.null(cellchat@meta$color) && length(cellchat@meta$color) == length(idents)) {
  color.use <- setNames(as.character(cellchat@meta$color), idents)
}

# 尝试从 cellchat@color 获取
if (is.null(color.use)) {
  try({
    tmp <- cellchat@color
    if (length(tmp) == length(idents)) color.use <- setNames(as.character(tmp), idents)
  }, silent = TRUE)
}

# 若仍无，自动生成
if (is.null(color.use)) {
  color.use <- hue_pal()(length(idents))
  names(color.use) <- idents
}

message("Color assignment complete.")

# ------------------------------
# 计算差异矩阵
# ------------------------------
net1_count <- cellchat@net$count[[compare.list[1]]]
net2_count <- cellchat@net$count[[compare.list[2]]]
net1_weight <- cellchat@net$weight[[compare.list[1]]]
net2_weight <- cellchat@net$weight[[compare.list[2]]]

diff_count <- net1_count - net2_count
diff_weight <- net1_weight - net2_weight

# ------------------------------
# 动态画布大小
# ------------------------------
n_ct <- length(idents)
canvas_width <- max(10, ceiling(n_ct / 3) * 2)
canvas_height <- canvas_width
message("Canvas size set to ", canvas_width, " x ", canvas_height, " inches")

# ------------------------------
# 循环绘制
# ------------------------------
for (ct in idents) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  
  # --- 数量网络 ---
  mat_num <- matrix(0, nrow = n_ct, ncol = n_ct, dimnames = list(idents, idents))
  mat_num[ct, ] <- diff_count[ct, ]
  mat_num[, ct] <- diff_count[, ct]
  
  vw_num <- rowSums(abs(mat_num))
  if (sum(vw_num) == 0) vw_num <- rep(1, n_ct)
  
  pdf(paste0(group[1], "VS", group[2], "_", ct_safe, "_Diff_number_net.pdf"),
      width = canvas_width, height = canvas_height)
  netVisual_circle(mat_num,
                   weight.scale = TRUE,
                   vertex.weight = vw_num,
                   vertex.color = color.use[rownames(mat_num)],
                   vertex.label.cex = 0.8,
                   edge.weight.max = max(abs(mat_num)),
                   title = paste0("Diff number of interactions (", ct, ")"))
  dev.off()
  
  # --- 强度网络 ---
  mat_wt <- matrix(0, nrow = n_ct, ncol = n_ct, dimnames = list(idents, idents))
  mat_wt[ct, ] <- diff_weight[ct, ]
  mat_wt[, ct] <- diff_weight[, ct]
  
  vw_wt <- rowSums(abs(mat_wt))
  if (sum(vw_wt) == 0) vw_wt <- rep(1, n_ct)
  
  pdf(paste0(group[1], "VS", group[2], "_", ct_safe, "_Diff_weight_net.pdf"),
      width = canvas_width, height = canvas_height)
  netVisual_circle(mat_wt,
                   weight.scale = TRUE,
                   vertex.weight = vw_wt,
                   vertex.color = color.use[rownames(mat_wt)],
                   vertex.label.cex = 0.8,
                   edge.weight.max = max(abs(mat_wt)),
                   title = paste0("Diff interaction strength (", ct, ")"))
  dev.off()
  
  message("Completed: ", ct)
}

message("All cell-type specific networks saved.")
