## ----------------------------
## 2. Helper functions
## ----------------------------

normalize_mode <- function(mode) {
  mode <- tolower(as.character(mode))
  if (mode == "dadius") {
    mode <- "radius"
  }
  if (!mode %in% c("knn", "radius")) {
    stop("neighbor_mode must be either 'knn' or 'radius'.")
  }
  mode
}

format_mode_value <- function(x) {
  format(x, trim = TRUE, scientific = FALSE)
}

make_one_hot_matrix <- function(x) {
  levels_sorted <- sort(unique(as.character(x)))
  factor_x <- factor(as.character(x), levels = levels_sorted)
  one_hot <- stats::model.matrix(~ factor_x - 1)
  colnames(one_hot) <- levels_sorted
  one_hot
}

get_knn_neighbors <- function(coords, k) {
  n_cells <- nrow(coords)
  k_use <- min(as.integer(k), n_cells)
  if (k_use < 1) {
    stop("k_neighbors must be >= 1.")
  }

  knn_res <- FNN::get.knnx(data = coords, query = coords, k = k_use)
  neighbors <- split(knn_res$nn.index, row(knn_res$nn.index))
  lapply(seq_len(n_cells), function(i) unique(c(i, as.integer(neighbors[[i]]))))
}

get_radius_neighbors <- function(coords, radius) {
  if (is.null(radius) || radius <= 0) {
    stop("radius_value must be > 0 when neighbor_mode = 'radius'.")
  }

  frnn_res <- dbscan::frNN(coords, eps = radius, sort = FALSE)
  lapply(seq_len(nrow(coords)), function(i) unique(c(i, as.integer(frnn_res$id[[i]]))))
}

build_neighborhood_matrix <- function(df,
                                      region_col,
                                      cluster_col,
                                      x_col,
                                      y_col,
                                      neighbor_mode,
                                      k_neighbors,
                                      radius_value) {
  required_cols <- c(region_col, cluster_col, x_col, y_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  na_counts <- vapply(required_cols, function(col) sum(is.na(df[[col]])), numeric(1))
  if (any(na_counts > 0)) {
    bad_cols <- paste(names(na_counts)[na_counts > 0], na_counts[na_counts > 0], sep = "=", collapse = ", ")
    stop("Required columns contain NA values: ", bad_cols)
  }

  df[[region_col]] <- as.character(df[[region_col]])
  df[[cluster_col]] <- as.character(df[[cluster_col]])

  one_hot <- make_one_hot_matrix(df[[cluster_col]])

  neighborhood_counts <- matrix(
    0,
    nrow = nrow(df),
    ncol = ncol(one_hot),
    dimnames = list(rownames(df), colnames(one_hot))
  )

  region_ids <- unique(df[[region_col]])

  for (region_id in region_ids) {
    idx <- which(df[[region_col]] == region_id)
    coords_region <- as.matrix(df[idx, c(x_col, y_col), drop = FALSE])
    storage.mode(coords_region) <- "double"

    region_one_hot <- one_hot[idx, , drop = FALSE]

    if (neighbor_mode == "knn") {
      neighbor_list <- get_knn_neighbors(coords_region, k_neighbors)
    } else {
      neighbor_list <- get_radius_neighbors(coords_region, radius_value)
    }

    for (local_i in seq_along(neighbor_list)) {
      local_neighbors <- neighbor_list[[local_i]]
      neighborhood_counts[idx[local_i], ] <- colSums(region_one_hot[local_neighbors, , drop = FALSE])
    }
  }

  window_df <- as.data.frame(neighborhood_counts, check.names = FALSE, stringsAsFactors = FALSE)
  window_df[[x_col]] <- df[[x_col]]
  window_df[[y_col]] <- df[[y_col]]
  window_df[[region_col]] <- df[[region_col]]
  window_df[[cluster_col]] <- df[[cluster_col]]

  list(
    window_df = window_df,
    feature_matrix = neighborhood_counts
  )
}

run_cn_analysis <- function(df,
                            region_col,
                            cluster_col,
                            x_col,
                            y_col,
                            neighbor_mode,
                            k_neighbors,
                            radius_value,
                            n_clusters,
                            random_seed,
                            kmeans_nstart,
                            kmeans_iter_max) {
  neighbor_mode <- normalize_mode(neighbor_mode)

  if (!is.numeric(n_clusters) || length(n_clusters) != 1 || is.na(n_clusters) || n_clusters <= 0) {
    stop("n_clusters must be a positive integer.")
  }
  n_clusters <- as.integer(n_clusters)

  neighborhood_res <- build_neighborhood_matrix(
    df = df,
    region_col = region_col,
    cluster_col = cluster_col,
    x_col = x_col,
    y_col = y_col,
    neighbor_mode = neighbor_mode,
    k_neighbors = k_neighbors,
    radius_value = radius_value
  )

  feature_matrix <- neighborhood_res$feature_matrix
  if (nrow(feature_matrix) < n_clusters) {
    stop(
      "n_clusters (", n_clusters, ") is larger than the number of available neighborhoods (",
      nrow(feature_matrix), ")."
    )
  }

  set.seed(random_seed)
  kmeans_res <- stats::kmeans(
    x = feature_matrix,
    centers = n_clusters,
    nstart = kmeans_nstart,
    iter.max = kmeans_iter_max
  )

  cn_labels <- as.integer(kmeans_res$cluster) - 1L

  if (neighbor_mode == "knn") {
    cn_col <- paste0("CN_k", k_neighbors, "_n", n_clusters)
  } else {
    cn_col <- paste0("CN_r", format_mode_value(radius_value), "_n", n_clusters)
  }

  result_df <- df
  result_df[[cn_col]] <- cn_labels

  list(
    result_df = result_df,
    window_df = neighborhood_res$window_df,
    feature_matrix = feature_matrix,
    kmeans_res = kmeans_res,
    cn_col = cn_col
  )
}
#
cn_exp_heatmap <- function(df, 
                             cluster_col, 
                             cn_col, 
                             row_clus = TRUE,  # 是否对行（邻域）进行聚类
                             col_clus = TRUE,  # 是否对列（细胞类型）进行聚类
                             vmin = NULL, 
                             vmax = NULL) {
  
  # ===========================
  # 1. 数据准备与矩阵计算
  # ===========================
  
  # 计算列联表
  counts <- table(df[[cn_col]], df[[cluster_col]])
  # 强制转换为基础矩阵
  counts_mat <- as.matrix(counts)
  
  # 归一化行 (Neighborhoods)
  niche_mat <- prop.table(counts_mat, 1) * 10
  
  # 计算 Tissue Averages (全局背景)
  tissue_avgs <- table(df[[cluster_col]]) / nrow(df)
  tissue_avgs <- as.numeric(tissue_avgs[colnames(niche_mat)]) 
  
  # 构造背景矩阵
  tissue_mat <- matrix(rep(tissue_avgs, each = nrow(niche_mat)), 
                       nrow = nrow(niche_mat), 
                       dimnames = dimnames(niche_mat))
  
  # 计算 Log2 Fold Change
  num <- niche_mat + tissue_mat
  denom <- rowSums(num)
  normalized_num <- sweep(num, 1, denom, "/")
  fc_2_mat <- log2(normalized_num / tissue_mat)
  
  # ===========================
  # 2. 行列分别无监督聚类
  # ===========================
  
  # --- 行聚类 (CN/Neighborhoods) ---
  if(row_clus && nrow(fc_2_mat) > 2) {
    # 替换 Inf/-Inf/NA 为 0 以便计算距离
    fc_clean <- fc_2_mat
    fc_clean[!is.finite(fc_clean)] <- 0 
    
    # 欧氏距离 + Complete Linkage 聚类
    dist_row <- dist(fc_clean, method = "euclidean")
    hclust_row <- hclust(dist_row, method = "complete")
    row_order <- rownames(fc_2_mat)[hclust_row$order]
  } else {
    # 不聚类则按字母顺序或原始顺序
    row_order <- sort(rownames(fc_2_mat))
  }
  
  # --- 列聚类 (Cell Type) ---
  if(col_clus && ncol(fc_2_mat) > 2) {
    fc_clean <- fc_2_mat
    fc_clean[!is.finite(fc_clean)] <- 0
    
    # 对转置矩阵计算距离
    dist_col <- dist(t(fc_clean), method = "euclidean")
    hclust_col <- hclust(dist_col, method = "complete")
    col_order <- colnames(fc_2_mat)[hclust_col$order]
  } else {
    col_order <- sort(colnames(fc_2_mat))
  }
  
  # ===========================
  # 3. 绘图数据转换
  # ===========================
  
  # 使用 as.data.frame.matrix 避免 "table" 类导致的 pivot_longer 报错
  plot_df <- as.data.frame.matrix(fc_2_mat) %>%
    rownames_to_column("CN") %>%
    pivot_longer(cols = -CN, names_to = "CellType", values_to = "Log2FC")
  
  # 找出每行最大值打星号
  plot_df <- plot_df %>%
    group_by(CN) %>%
    mutate(is_max = Log2FC == max(Log2FC, na.rm = TRUE)) %>% 
    ungroup() %>%
    mutate(label = ifelse(is_max, "*", ""))
  
  # === 关键：应用聚类顺序 ===
  plot_df$CN <- factor(plot_df$CN, levels = row_order)
  plot_df$CellType <- factor(plot_df$CellType, levels = col_order)
  
  # ===========================
  # 4. 动态颜色范围
  # ===========================
  
  # 排除 Inf 计算最大值
  valid_vals <- plot_df$Log2FC[is.finite(plot_df$Log2FC)]
  max_abs <- if(length(valid_vals) > 0) max(abs(valid_vals), na.rm = TRUE) else 1
  if(max_abs == 0) max_abs <- 1
  
  use_vmin <- if(is.null(vmin)) -max_abs else vmin
  use_vmax <- if(is.null(vmax)) max_abs else vmax
  
  # ===========================
  # 5. 绘图 (ggplot2)
  # ===========================
  
  p <- ggplot(plot_df, aes(x = CellType, y = CN, fill = Log2FC)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(aes(label = label), vjust = 0.75, size = 5, color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                         midpoint = 0, limit = c(use_vmin, use_vmax),
                         name = "Log2 FC\n(vs Tissue Avg)",
                         na.value = "grey90") +
    theme_minimal() +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid = element_blank(),
      legend.position = "right"
    ) 
  
  return(p)
}