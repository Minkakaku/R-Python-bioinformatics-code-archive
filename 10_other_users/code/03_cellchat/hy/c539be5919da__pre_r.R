# net_source_only_each_celltype_enhanced.R
# 功能：从 net_lr 文件读取，为每个细胞类型单独绘制以该细胞为 source 的网络图
#       分别绘制“数量”和“强度”两个版本，节点、字体、箭头、线宽均显著放大
# 修改：中心细胞与周围细胞节点大小统一为15，曲线曲度加大至0.5
# 输入：当前目录下的 "PRE-R.net_lr.txt"

rm(list = ls())
library(igraph)
library(scales)

# ------------------------------
# 1. 读取数据
# ------------------------------
net_file <- "PRE-R.net_lr.txt"
if (!file.exists(net_file)) stop("File not found: ", net_file)

net <- read.table(net_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                  comment.char = "", quote = "")
stopifnot(all(c("source", "target", "prob") %in% colnames(net)))

# ------------------------------
# 2. 获取所有细胞类型
# ------------------------------
celltypes <- sort(union(unique(net$source), unique(net$target)))
n_ct <- length(celltypes)
message("Total unique cell types: ", n_ct)

# ------------------------------
# 3. 构建强度矩阵和数量矩阵
# ------------------------------
agg_strength <- aggregate(prob ~ source + target, data = net, FUN = sum)
agg_count <- aggregate(prob ~ source + target, data = net, FUN = length)
colnames(agg_count)[3] <- "count"

mat_strength <- matrix(0, nrow = n_ct, ncol = n_ct, dimnames = list(celltypes, celltypes))
mat_count    <- matrix(0, nrow = n_ct, ncol = n_ct, dimnames = list(celltypes, celltypes))

for (i in seq_len(nrow(agg_strength))) {
  s <- agg_strength$source[i]; t <- agg_strength$target[i]
  mat_strength[s, t] <- agg_strength$prob[i]
  mat_count[s, t] <- agg_count$count[i]   # 注意 agg_count 与 agg_strength 行顺序相同
}

# ------------------------------
# 4. 颜色分配
# ------------------------------
color_pal <- hue_pal()(n_ct)
names(color_pal) <- celltypes

# ------------------------------
# 5. 动态画布（增大基础尺寸）
# ------------------------------
canvas_size <- max(16, ceiling(n_ct / 3) * 2.5)  # 基础14英寸，每3个细胞类型加2.5英寸
message("Canvas size set to ", canvas_size, " x ", canvas_size, " inches")

# ------------------------------
# 6. 绘图函数（统一节点大小，加大曲度）
# ------------------------------
draw_net <- function(mat_vals, ct, ct_safe, suffix) {
  # 只保留中心细胞作为 source 的行
  mat_ct <- matrix(0, nrow = n_ct, ncol = n_ct, dimnames = dimnames(mat_vals))
  mat_ct[ct, ] <- mat_vals[ct, ]
  
  g <- graph_from_adjacency_matrix(mat_ct, mode = "directed", weighted = TRUE)
  lay <- layout_in_circle(g)
  
  # 节点大小：中心细胞与周围细胞统一为15
  v_size <- rep(15, n_ct)            # 所有节点大小一致
  
  # 字体大小
  v_label_cex <- 2
  
  # 边宽度缩放因子
  edge_scale <- 2.5   # 可调，越大线越粗
  e_weight <- E(g)$weight
  if (max(e_weight) > 0) {
    e_width <- edge_scale * (1.25 + 5.0 * (e_weight / max(e_weight)))
  } else {
    e_width <- rep(0, ecount(g))
  }
  
  # 边颜色：深灰半透明
  e_col <- rgb(0.2, 0.2, 0.2, alpha = 0.8)
  
  pdf_file <- paste0("PRE-R_", ct_safe, "_", suffix, "_network.pdf")
  pdf(pdf_file, width = canvas_size, height = canvas_size)
  
  par(mar = c(1, 1, 5, 1))
  
  plot(g,
       layout = lay,
       vertex.size = v_size,
       vertex.color = color_pal[V(g)$name],
       vertex.label = V(g)$name,
       vertex.label.cex = v_label_cex,
       vertex.label.color = "black",
       vertex.frame.color = NA,
       edge.width = e_width,
       edge.arrow.size = 0.8,
       edge.color = e_col,
       edge.curved = 0.5,            # 曲度加大，避免与节点圆重叠
       main = paste0("PRE-R ", suffix, " (source only): ", ct),
       cex.main = 1.8)
  
  dev.off()
  message("Saved: ", pdf_file)
}

# ------------------------------
# 7. 循环绘制每个细胞类型的强度图和数量图
# ------------------------------
for (ct in celltypes) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  draw_net(mat_strength, ct, ct_safe, "strength")
  draw_net(mat_count,    ct, ct_safe, "count")
}

message("All enhanced networks generated.")
