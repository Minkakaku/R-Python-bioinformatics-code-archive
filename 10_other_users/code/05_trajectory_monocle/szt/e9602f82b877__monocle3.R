suppressPackageStartupMessages({
  library(optparse); library(Seurat); library(monocle3)
  library(tidyverse); library(patchwork); library(ggridges)
  library(viridis); library(pheatmap); library(anndataR); library(CytoTRACE2)
  set.seed(1234)
})

# ===== 命令行参数 =====
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "输入文件(.rds/.h5ad) [必需]"),
  make_option(c("-o", "--outdir"), type = "character", default = "./monocle3_results", help = "输出目录"),
  make_option("--celltype_col", type = "character", default = "celltype", help = "细胞类型列名"),
  make_option("--group_col", type = "character", default = "group", help = "分组列名(留空跳过)"),
  make_option("--sample_col", type = "character", default = "sample", help = "样本列名(下采样分层用)"),
  make_option("--subsample_prop", type = "double", default = 1.0, help = "下采样比例(0-1)"),
  make_option("--species", type = "character", default = "human", help = "物种(CytoTRACE2用)"),
  # 轨迹关键参数
  make_option("--root_type", type = "character", default = "", help = "起始细胞类型(留空=CytoTRACE2 highest potency自动选根)"),
  make_option("--close_loop", type = "character", default = "false", help = "是否闭合轨迹环路"),
  make_option("--num_dim", type = "integer", default = 30, help = "PCA维度数"),
  make_option("--cluster_resolution", type = "double", default = 1e-4, help = "轨迹聚类分辨率(值越低partition越少)"),
  make_option("--top_genes", type = "integer", default = 50, help = "可视化top基因数"),
  # 并行与输出
  make_option("--workers", type = "integer", default = 6, help = "并行核心数"),
  make_option("--width", type = "double", default = 10, help = "图片宽度(英寸)"),
  make_option("--height", type = "double", default = 8, help = "图片高度(英寸)"),
  make_option("--dpi", type = "integer", default = 300, help = "图片分辨率")
)
opt_parser <- OptionParser(option_list = option_list, usage = "%prog -i <input.rds|.h5ad> -o <output_dir> [options]")
opt <- parse_args(opt_parser)

# ---- DEBUG（正式运行时注释掉） ----
setwd('E:/Projects/scRNA/Human_PBMC')
opt <- list(
  input="./output/adata_celltype_neutr.rds", celltype_col="celltype",
  sample_col="sample", group_col="group", outdir="./output/monocle3",
  root_type="", close_loop="false", num_dim=30L,
  cluster_resolution=1e-4, subsample_prop=1, species="human",
  top_genes=50L, workers=6L, width=10, height=8, dpi=300L
)

# 参数后处理 & 校验
opt$close_loop <- tolower(as.character(opt$close_loop)) %in% c("true", "1", "yes")
if (opt$group_col == "") opt$group_col <- NULL
if (opt$root_type == "") opt$root_type <- NULL
if (is.null(opt$input)) { print_help(opt_parser); stop("错误: 必须指定 --input", call. = FALSE) }
if (!file.exists(opt$input)) stop(paste("错误: 文件不存在:", opt$input), call. = FALSE)
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

cat("=== Monocle3 拟时序分析 ===\n")
cat("输入:", opt$input, "| 输出:", opt$outdir, "| 细胞类型列:", opt$celltype_col, "\n")

# ===== 辅助函数 =====

# 将一组细胞 ID 映射到 principal graph 中出现最多的节点
cells_to_root <- function(cds, cell_ids) {
  cv <- as.matrix(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[colnames(cds), ])
  igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(cv[cell_ids, ]))))
  ]
}

# 根据细胞类型选根节点
get_root_node <- function(cds, celltype_col, target_type) {
  ids <- which(colData(cds)[[celltype_col]] == target_type)
  if (length(ids) == 0) { warning(paste("未找到细胞类型:", target_type)); return(NULL) }
  cells_to_root(cds, ids)
}

# 基于 CytoTRACE2 highest potency 自动选根（取 top_n 高分化潜能细胞）
auto_root_node <- function(cds, ct2_scores, top_n = 200) {
  common <- intersect(names(ct2_scores), colnames(cds))
  top_cells <- names(sort(ct2_scores[common], decreasing = TRUE))[1:min(top_n, length(common))]
  cells_to_root(cds, top_cells)
}

# 安全保存 ggplot / pheatmap
save_plot <- function(filename, plot, w = opt$width, h = opt$height, ...) {
  tryCatch(ggsave(filename, plot, width = w, height = h, dpi = opt$dpi, ...),
           error = function(e) warning(paste("保存失败:", basename(filename), "-", e$message)))
}
save_pheatmap <- function(filename, mat, title, w = opt$width, h = opt$height, fcol = 9) {
  tryCatch({
    png(filename, width = w, height = h, units = "in", res = opt$dpi)
    pheatmap(mat, scale = "row", clustering_method = "ward.D2",
             fontsize_row = 10, fontsize_col = fcol,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100), main = title)
    dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE); warning(paste("热图失败:", e$message)) })
}

# 构建模块聚合热图的 cell_group tibble
make_cell_group <- function(cds, ...) {
  cols <- list(...)
  tibble(cell = rownames(colData(cds)),
         cell_group = Reduce(function(a, b) paste(a, b, sep = " | "),
                             lapply(cols, function(col) as.character(colData(cds)[[col]]))))
}

# ===== Step 1: 数据加载 =====
cat("\n[Step 1] 加载数据...\n")
input_ext <- tolower(tools::file_ext(opt$input))
seurat_obj <- switch(input_ext,
  "rds" = readRDS(opt$input),
  "h5ad" = { cat("  h5ad → Seurat (anndataR)...\n"); read_h5ad(opt$input, as = "Seurat") },
  stop(paste("错误: 不支持的格式:", input_ext), call. = FALSE))
cat("细胞:", ncol(seurat_obj), "| 基因:", nrow(seurat_obj), "\n")

# 校验
stopifnot(opt$celltype_col %in% colnames(seurat_obj@meta.data))
if (!is.null(opt$group_col)) stopifnot(opt$group_col %in% colnames(seurat_obj@meta.data))
umap_key <- if ("umap" %in% Reductions(seurat_obj)) "umap" else
            if ("X_umap" %in% Reductions(seurat_obj)) "X_umap" else
            stop("错误: 未找到 UMAP 降维", call. = FALSE)
cat("UMAP key:", umap_key, "\n\n细胞类型分布:\n")
print(table(seurat_obj@meta.data[[opt$celltype_col]]))

# ===== Step 2: 下采样（可选） =====
if (opt$subsample_prop < 1.0) {
  cat("\n[Step 2] 分层下采样, 比例:", opt$subsample_prop, "\n")
  meta <- seurat_obj@meta.data %>% mutate(..cell_id = rownames(seurat_obj@meta.data))
  gv <- opt$celltype_col
  if (opt$sample_col %in% colnames(meta)) gv <- c(gv, opt$sample_col)
  ids <- meta %>% group_by(across(all_of(gv))) %>%
    slice_sample(prop = opt$subsample_prop) %>% pull(..cell_id)
  seurat_obj <- seurat_obj[, ids]
  cat("下采样后:", ncol(seurat_obj), "细胞\n")
} else cat("\n[Step 2] 不下采样\n")

# ===== Step 3: CytoTRACE2 分化潜能预测 =====
cat("\n[Step 3] CytoTRACE2 分化潜能预测...\n")
cytotrace2_result <- cytotrace2(
  seurat_obj, 
  is_seurat = TRUE, 
  species = opt$species, 
  slot_type = "counts",
  batch_size = 20000,   # 大数据集分批处理
  ncores = 4,
  seed = 1234
)
annotation <- data.frame(phenotype = seurat_obj@meta.data[[opt$celltype_col]],
                          row.names = colnames(seurat_obj))
plots <- plotData(cytotrace2_result, annotation, is_seurat = TRUE)
p <- (plots$CytoTRACE2_UMAP + plots$CytoTRACE2_Potency_UMAP +
      plots$CytoTRACE2_Relative_UMAP + plots$CytoTRACE2_Boxplot_byPheno) + plot_layout(ncol = 2)
save_plot(file.path(opt$outdir, "cytotrace2.png"), p)

# ===== Step 4: Seurat → CDS + 预处理 + 注入 UMAP =====
cat("\n[Step 4] 构建 CDS + 预处理 + 注入 UMAP...\n")
data <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
cds <- new_cell_data_set(data, seurat_obj@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(data), row.names = rownames(data)))
# preprocess: size factor 归一化 + PCA（graph_test 需要 size factor）
cds <- preprocess_cds(cds, num_dim = opt$num_dim, method = "PCA")
# 直接注入 Seurat UMAP（learn_graph/cluster_cells 只在 UMAP 空间操作）
reducedDims(cds)[["UMAP"]] <- Embeddings(seurat_obj, reduction = umap_key)[colnames(cds), ]
cat("CDS:", ncol(cds), "cells |", nrow(cds), "genes | UMAP:", umap_key, "\n")

# ===== Step 5: 聚类 + 学习轨迹 =====
cat("\n[Step 5] 聚类 & 学习轨迹...\n")
cds <- cluster_cells(cds, resolution = opt$cluster_resolution, k = 20)
cds <- learn_graph(cds, use_partition = FALSE, close_loop = opt$close_loop)
cat("Partitions:", length(unique(partitions(cds))),
    "| Clusters:", length(unique(clusters(cds))), "\n")

# ===== Step 6: 选根 + Pseudotime =====
cat("\n[Step 6] 选择根节点 & 计算 Pseudotime...\n")
root_node <- if (!is.null(opt$root_type)) get_root_node(cds, opt$celltype_col, opt$root_type) else NULL
if (is.null(root_node)) {
  ct2_scores <- setNames(cytotrace2_result@meta.data$CytoTRACE2_Score, colnames(cytotrace2_result))
  root_node <- auto_root_node(cds, ct2_scores)
  cat("自动选根 (CytoTRACE2 highest potency):", root_node, "\n")
} else cat("根节点:", root_node, "(", opt$root_type, ")\n")

cds <- order_cells(cds, root_pr_nodes = root_node)
pt <- pseudotime(cds)
cat("Pseudotime: [", round(min(pt[is.finite(pt)]), 2), ",",
    round(max(pt[is.finite(pt)]), 2), "] | Inf:", sum(is.infinite(pt)), "/", length(pt), "\n")

# ===== Step 7: 轨迹可视化 =====
cat("\n[Step 7] 轨迹可视化...\n")
p_ct <- plot_cells(cds, color_cells_by = opt$celltype_col, label_groups_by_cluster = FALSE,
  label_leaves = TRUE, label_branch_points = TRUE, group_label_size = 5, graph_label_size = 3) +
  ggtitle("Trajectory — Cell Type")
p_pt <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
  label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 3) +
  scale_color_viridis_c() + ggtitle("Trajectory — Pseudotime")
p_par <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) +
  ggtitle("Partitions")
save_plot(file.path(opt$outdir, "trajectory_celltype.png"), p_ct)
save_plot(file.path(opt$outdir, "trajectory_pseudotime.png"), p_pt)
save_plot(file.path(opt$outdir, "trajectory_partition.png"), p_par)

p_grp <- NULL
if (!is.null(opt$group_col)) {
  p_grp <- plot_cells(cds, color_cells_by = opt$group_col, label_groups_by_cluster = FALSE,
    label_leaves = FALSE, label_branch_points = FALSE) + ggtitle(paste("Trajectory —", opt$group_col))
  save_plot(file.path(opt$outdir, "trajectory_group.png"), p_grp)
}
save_plot(file.path(opt$outdir, "trajectory_combined.png"),
          (p_ct + p_pt) / (p_par + (p_grp %||% plot_spacer())),
          w = opt$width * 2, h = opt$height * 2, limitsize = FALSE)

# ===== Step 8: Pseudotime 分布 =====
cat("\n[Step 8] Pseudotime 分布可视化...\n")
pt_df <- data.frame(pseudotime = pseudotime(cds), cell_type = colData(cds)[[opt$celltype_col]])
if (!is.null(opt$group_col)) pt_df$group <- colData(cds)[[opt$group_col]]
pt_df <- pt_df[is.finite(pt_df$pseudotime), ]

save_plot(file.path(opt$outdir, "pseudotime_ridge.png"),
  ggplot(pt_df, aes(x = pseudotime, y = cell_type, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "Pseudotime", option = "C") +
    labs(title = "Pseudotime by Cell Type", x = "Pseudotime", y = NULL) +
    theme_classic())

if (!is.null(opt$group_col)) {
  save_plot(file.path(opt$outdir, "pseudotime_density.png"),
    ggplot(pt_df, aes(x = pseudotime, fill = group)) + geom_density(alpha = 0.5) +
      labs(title = "Pseudotime by Group", fill = opt$group_col) + theme_classic(),
    h = opt$height * 0.7)
  save_plot(file.path(opt$outdir, "pseudotime_density_facet.png"),
    ggplot(pt_df, aes(x = pseudotime, fill = group)) + geom_density(alpha = 0.5) +
      facet_wrap(~cell_type, scales = "free_y") +
      labs(title = "Pseudotime by Cell Type & Group", fill = opt$group_col) + theme_classic(),
    w = opt$width * 1.5, h = opt$height * 1.2, limitsize = FALSE)
}

# ===== Step 9: graph_test — 轨迹 DEG =====
cat("\n[Step 9] graph_test (Moran's I)...\n")
Track_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = opt$workers)
Track_genes_sig <- Track_genes %>% filter(q_value < 0.05, morans_I > 0.15) %>% arrange(desc(morans_I))
n_sig <- nrow(Track_genes_sig)
cat("  显著基因:", n_sig, "(q<0.05, I>0.15)\n")
write.csv(Track_genes, file.path(opt$outdir, "graph_test_all.csv"), row.names = TRUE)
write.csv(Track_genes_sig, file.path(opt$outdir, "graph_test_significant.csv"), row.names = TRUE)

# Top 基因 UMAP
if (n_sig >= 6)
  save_plot(file.path(opt$outdir, "top_genes_umap.png"),
    plot_cells(cds, genes = head(Track_genes_sig$gene_short_name, 6),
      show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE),
    w = opt$width * 1.5, limitsize = FALSE)

# Top 基因沿 Pseudotime 趋势
top_n <- min(opt$top_genes, n_sig)
if (top_n > 0) {
  vg <- intersect(as.character(head(Track_genes_sig$gene_short_name, top_n)), rownames(cds))
  if (length(vg) > 0) tryCatch(
    save_plot(file.path(opt$outdir, "genes_in_pseudotime.png"),
      plot_genes_in_pseudotime(cds[vg, ], color_cells_by = opt$celltype_col, min_expr = 0.5, ncol = 2),
      h = opt$height * 1.5, limitsize = FALSE),
    error = function(e) warning(paste("genes_in_pseudotime 失败:", e$message)))
}

# ===== Step 10: 共表达模块分析 =====
cat("\n[Step 10] 共表达模块分析...\n")
module_genes <- head(Track_genes_sig, min(500, n_sig)) %>% pull(gene_short_name) %>% as.character()

if (length(module_genes) >= 20) {
  tryCatch({
    gene_module_df <- find_gene_modules(cds[module_genes, ], resolution = 1e-2, cores = opt$workers)
    cat("  模块数:", length(unique(gene_module_df$module)), "\n")
    write.csv(gene_module_df, file.path(opt$outdir, "gene_modules.csv"), row.names = FALSE)

    # 按细胞类型聚合热图
    agg <- aggregate_gene_expression(cds, gene_module_df, make_cell_group(cds, opt$celltype_col))
    rownames(agg) <- paste0("Module ", rownames(agg))
    save_pheatmap(file.path(opt$outdir, "gene_modules_heatmap.png"), agg, "Gene Modules × Cell Types")

    # 按分组聚合热图（可选）
    if (!is.null(opt$group_col)) {
      agg2 <- aggregate_gene_expression(cds, gene_module_df,
                                         make_cell_group(cds, opt$celltype_col, opt$group_col))
      rownames(agg2) <- paste0("Module ", rownames(agg2))
      save_pheatmap(file.path(opt$outdir, "gene_modules_heatmap_by_group.png"),
                    agg2, "Gene Modules × Cell Type | Group", w = opt$width * 1.5, fcol = 7)
    }
  }, error = function(e) warning(paste("模块分析失败:", e$message)))
} else cat("  显著基因不足20个，跳过\n")

# ===== Step 11: 保存 =====
cat("\n[Step 11] 保存结果...\n")
saveRDS(cds, file.path(opt$outdir, "monocle3_cds.rds"))
cat("\n=== 完成 ===\n输出目录:", opt$outdir, "\n")
for (f in list.files(opt$outdir)) cat("  ", f, "\n")
