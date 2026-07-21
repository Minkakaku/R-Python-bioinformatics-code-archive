suppressPackageStartupMessages({
  library(optparse)
  library(CellChat)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  set.seed(1234)
})

# =============================================================================
# 命令行参数
# =============================================================================
option_list <- list(
  # 输入输出
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "输入Seurat RDS文件 [必需]"),
  make_option(c("-o", "--outdir"), type = "character", default = "./cellchat_results", help = "输出目录"),
  # 核心参数
  make_option(c("-c", "--celltype_col"), type = "character", default = "celltype", help = "细胞类型列名"),
  make_option(c("--species"), type = "character", default = "human", help = "物种: human 或 mouse"),
  make_option(c("--db_category"), type = "character", default = "all", help = "数据库类型: all/Secreted Signaling/ECM-Receptor/Cell-Cell Contact"),
  make_option(c("--min_cells"), type = "integer", default = 10, help = "最少细胞数过滤"),
  make_option(c("--compute_type"), type = "character", default = "triMean", help = "通讯概率计算方法: triMean/truncatedMean"),
  make_option(c("--workers"), type = "integer", default = 4, help = "并行核心数"),
  # 分组比较参数
  make_option(c("--group_col"), type = "character", default = NULL, help = "分组列名(指定后启用比较分析模式)"),
  # 输出设置
  make_option(c("--width"), type = "double", default = 10, help = "图片宽度(英寸)"),
  make_option(c("--height"), type = "double", default = 8, help = "图片高度(英寸)"),
  make_option(c("--dpi"), type = "integer", default = 300, help = "图片分辨率")
)

opt_parser <- OptionParser(option_list = option_list, usage = "%prog -i <input.rds> -o <output_dir> [options]")
opt <- parse_args(opt_parser)

# 参数检查
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("错误: 必须指定 --input", call. = FALSE)
}
if (!file.exists(opt$input)) stop(paste("错误: 文件不存在:", opt$input), call. = FALSE)
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

cat("=== CellChat 细胞通讯分析 ===\n")
cat("输入:", opt$input, "| 输出:", opt$outdir, "\n")
cat("物种:", opt$species, "| 细胞类型列:", opt$celltype_col, "\n")
if (!is.null(opt$group_col)) cat("分组比较模式: 按", opt$group_col, "列分组\n")

# =============================================================================
# 公共函数: 对单组数据运行 CellChat 分析
# =============================================================================
run_cellchat <- function(seurat_sub, celltype_col, species, db_category, min_cells, compute_type, workers) {
  cellchat <- createCellChat(object = seurat_sub, group.by = celltype_col, assay = "RNA")

  # 设置数据库
  if (species == "mouse") {
    CellChatDB <- CellChatDB.mouse
  } else {
    CellChatDB <- CellChatDB.human
  }
  if (db_category != "all") {
    CellChatDB.use <- subsetDB(CellChatDB, search = db_category, key = "annotation")
  } else {
    CellChatDB.use <- CellChatDB
  }
  cellchat@DB <- CellChatDB.use

  # 预处理 + 推断
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = workers)
  cellchat <- tryCatch(
    identifyOverExpressedGenes(cellchat, do.fast = FALSE),
    error = function(e) identifyOverExpressedGenes(cellchat)
  )
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = compute_type)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat)

  return(cellchat)
}

# =============================================================================
# Step 1: 数据加载
# =============================================================================
cat("\n[Step 1] 加载数据...\n")
seurat_obj <- readRDS(opt$input)
cat("细胞数:", ncol(seurat_obj), "| 基因数:", nrow(seurat_obj), "\n")

if (!opt$celltype_col %in% colnames(seurat_obj@meta.data)) {
  stop(paste("错误: metadata中不存在列:", opt$celltype_col), call. = FALSE)
}

# =============================================================================
# 路由: 单组分析 vs 分组比较
# =============================================================================
if (is.null(opt$group_col)) {
  # ======================= 单组分析模式 =======================
  cat("\n========== 单组分析模式 ==========\n")
  cat("\n细胞类型分布:\n")
  print(table(seurat_obj@meta.data[[opt$celltype_col]]))

  cat("\n[Step 2] 运行CellChat...\n")
  cellchat <- run_cellchat(seurat_obj, opt$celltype_col, opt$species, opt$db_category, opt$min_cells, opt$compute_type, opt$workers)

  # 导出通讯数据
  df_net_lr <- subsetCommunication(cellchat)
  df_net_pathway <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df_net_lr, file.path(opt$outdir, "communications_LR.csv"), row.names = FALSE)
  write.csv(df_net_pathway, file.path(opt$outdir, "communications_pathway.csv"), row.names = FALSE)
  cat("L-R对数:", nrow(df_net_lr), "| 通路数:", nrow(df_net_pathway), "\n")

  # 可视化
  cat("\n[Step 3] 生成可视化...\n")
  groupSize <- as.numeric(table(cellchat@idents))

  # 聚合网络
  png(file.path(opt$outdir, "network_circle.png"), width = opt$width * 2, height = opt$height, units = "in", res = opt$dpi)
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
  dev.off()

  # 各细胞类型信号发送
  mat <- cellchat@net$weight
  n_ct <- nrow(mat)
  n_col <- min(4, n_ct)
  n_row <- ceiling(n_ct / n_col)
  png(file.path(opt$outdir, "network_circle_per_celltype.png"), width = n_col * 4, height = n_row * 4, units = "in", res = opt$dpi)
  par(mfrow = c(n_row, n_col), xpd = TRUE)
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  dev.off()

  # 信号通路热图
  tryCatch({
    png(file.path(opt$outdir, "heatmap_pathway.png"), width = opt$width, height = opt$height, units = "in", res = opt$dpi)
    netVisual_heatmap(cellchat)
    dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE); warning(e$message) })

  # Bubble Plot
  tryCatch({
    p_bubble <- netVisual_bubble(cellchat, remove.isolate = FALSE)
    ggsave(file.path(opt$outdir, "bubble_plot.png"), p_bubble, width = opt$width, height = opt$height * 1.5, dpi = opt$dpi, limitsize = FALSE)
  }, error = function(e) warning(e$message))

  # 信号角色热图
  tryCatch({
    png(file.path(opt$outdir, "signaling_role_heatmap.png"), width = opt$width, height = opt$height, units = "in", res = opt$dpi)
    netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 12)
    dev.off()
    png(file.path(opt$outdir, "signaling_role_heatmap_incoming.png"), width = opt$width, height = opt$height, units = "in", res = opt$dpi)
    netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 12)
    dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE); warning(e$message) })

  # 保存
  saveRDS(cellchat, file.path(opt$outdir, "cellchat_object.rds"))
  cat("\n=== 单组分析完成 ===\n")

} else {
  # ======================= 分组比较模式 =======================
  cat("\n========== 分组比较模式 ==========\n")

  if (!opt$group_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("错误: metadata中不存在分组列:", opt$group_col), call. = FALSE)
  }

  groups <- unique(seurat_obj@meta.data[[opt$group_col]])
  if (length(groups) != 2) {
    stop(paste("错误: 分组列必须恰好有2个分组，当前有:", length(groups), "个 ->", paste(groups, collapse = ", ")), call. = FALSE)
  }
  cat("分组:", paste(groups, collapse = " vs "), "\n")

  # Step 2: 分别对每组运行 CellChat
  object.list <- list()
  for (g in groups) {
    cat("\n[Step 2] 运行CellChat:", g, "...\n")
    cells <- which(seurat_obj@meta.data[[opt$group_col]] == g)
    seurat_sub <- seurat_obj[, cells]
    cat("  细胞数:", ncol(seurat_sub), "\n")
    print(table(seurat_sub@meta.data[[opt$celltype_col]]))
    object.list[[g]] <- run_cellchat(seurat_sub, opt$celltype_col, opt$species, opt$db_category, opt$min_cells, opt$compute_type, opt$workers)
  }

  # Step 3: 合并
  cat("\n[Step 3] 合并CellChat对象...\n")
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))

  # 导出各组通讯数据
  for (g in names(object.list)) {
    df_lr <- subsetCommunication(object.list[[g]])
    df_pw <- subsetCommunication(object.list[[g]], slot.name = "netP")
    write.csv(df_lr, file.path(opt$outdir, paste0("communications_LR_", g, ".csv")), row.names = FALSE)
    write.csv(df_pw, file.path(opt$outdir, paste0("communications_pathway_", g, ".csv")), row.names = FALSE)
  }

  # Step 4: 比较可视化
  cat("\n[Step 4] 生成比较可视化...\n")

  # 4.1 总交互数/强度比较条形图
  cat("  [4.1] 交互数量/强度比较...\n")
  tryCatch({
    gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
    gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
    ggsave(file.path(opt$outdir, "compare_interactions.png"), gg1 + gg2, width = opt$width, height = opt$height * 0.5, dpi = opt$dpi)
  }, error = function(e) warning(paste("交互比较失败:", e$message)))

  # 4.2 差异网络 Circle Plot
  cat("  [4.2] 差异网络Circle Plot...\n")
  tryCatch({
    png(file.path(opt$outdir, "diff_network_circle.png"), width = opt$width * 2, height = opt$height, units = "in", res = opt$dpi)
    par(mfrow = c(1, 2), xpd = TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
    dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE); warning(e$message) })

  # 4.3 差异热图
  cat("  [4.3] 差异热图...\n")
  tryCatch({
    gg1 <- netVisual_heatmap(cellchat)
    gg2 <- netVisual_heatmap(cellchat, measure = "weight")
    png(file.path(opt$outdir, "diff_heatmap.png"), width = opt$width * 2, height = opt$height, units = "in", res = opt$dpi)
    print(gg1 + gg2)
    dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE); warning(e$message) })

  # 4.4 各组聚合网络对比
  cat("  [4.4] 各组网络对比...\n")
  tryCatch({
    weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
    png(file.path(opt$outdir, "compare_network_count.png"), width = opt$width * 2, height = opt$height, units = "in", res = opt$dpi)
    par(mfrow = c(1, 2), xpd = TRUE)
    for (i in 1:length(object.list)) {
      netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE, label.edge = FALSE,
                       edge.weight.max = weight.max[2], edge.width.max = 12,
                       title.name = paste0("Number of interactions - ", names(object.list)[i]))
    }
    dev.off()

    weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))
    png(file.path(opt$outdir, "compare_network_weight.png"), width = opt$width * 2, height = opt$height, units = "in", res = opt$dpi)
    par(mfrow = c(1, 2), xpd = TRUE)
    for (i in 1:length(object.list)) {
      netVisual_circle(object.list[[i]]@net$weight, weight.scale = TRUE, label.edge = FALSE,
                       edge.weight.max = weight.max[2], edge.width.max = 12,
                       title.name = paste0("Interaction strength - ", names(object.list)[i]))
    }
    dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE); warning(e$message) })

  # 4.5 信号通路信息流排名（Information Flow）
  cat("  [4.5] 信号通路信息流排名...\n")
  tryCatch({
    gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE)
    gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", stacked = FALSE, do.stat = TRUE)
    ggsave(file.path(opt$outdir, "information_flow_compare.png"), gg1 + gg2, width = opt$width * 1.5, height = opt$height, dpi = opt$dpi, limitsize = FALSE)
  }, error = function(e) warning(paste("信息流排名失败:", e$message)))

  # 4.6 各组信号角色热图对比
  cat("  [4.6] 信号角色对比热图...\n")
  tryCatch({
    library(ComplexHeatmap)
    pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)

    png(file.path(opt$outdir, "compare_signaling_outgoing.png"), width = opt$width * 1.5, height = opt$height * 1.2, units = "in", res = opt$dpi)
    ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[1], width = 5, height = 6)
    ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[2], width = 5, height = 6)
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    dev.off()

    png(file.path(opt$outdir, "compare_signaling_incoming.png"), width = opt$width * 1.5, height = opt$height * 1.2, units = "in", res = opt$dpi)
    ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[1], width = 5, height = 6, color.heatmap = "GnBu")
    ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[2], width = 5, height = 6, color.heatmap = "GnBu")
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE); warning(paste("信号角色热图失败:", e$message)) })

  # Step 5: 保存
  cat("\n[Step 5] 保存结果...\n")
  saveRDS(cellchat, file.path(opt$outdir, "cellchat_merged.rds"))
  saveRDS(object.list, file.path(opt$outdir, "cellchat_object_list.rds"))

  cat("\n=== 分组比较分析完成 ===\n")
  cat("比较可视化:\n")
  cat("  compare_interactions.png - 交互数量/强度条形图对比\n")
  cat("  diff_network_circle.png - 差异网络图（红=增强，蓝=减弱）\n")
  cat("  diff_heatmap.png - 差异热图\n")
  cat("  compare_network_count/weight.png - 各组网络对比\n")
  cat("  information_flow_compare.png - 信号通路信息流排名\n")
  cat("  compare_signaling_outgoing/incoming.png - 信号角色对比热图\n")
  cat("  communications_LR/pathway_*.csv - 各组通讯数据表\n")
}
