library(dplyr)
library(tibble)
library(tidyr) 

#' 计算组内细胞间的最小距离 (V7 - 移除 rownames_to_column)
#'
#' @param input_df 包含所有样本数据的完整数据框。
#' @param source_type 源细胞类型。如果为 NULL，则触发所有配对模式。
#' @param target_type 目标细胞类型。如果为 NULL，则触发所有配对模式。
#' @param group_col 用于分组的列名 (e.g., "sample")
#' @param celltype_col 包含细胞类型注释的列名 (e.g., "sub")
#' @param x_col 空间 X 坐标的列名 (e.g., "X")
#' @param y_col 空间 Y 坐标的列名 (e.g., "Y")
#'
#' @return 一个数据框，包含计算的最小距离。
#'
calculate_distances <- function(input_df, 
                                         source_type = NULL, 
                                         target_type = NULL, 
                                         group_col = "sample",
                                         celltype_col = "celltype",
                                         x_col = "X",
                                         y_col = "Y") {
  
  # --- 内部辅助函数 1: 距离计算 ---
  .calculate_min_distance <- function(cell_data, src_type, tgt_type, ctype_col, xc, yc) {
    source <- cell_data[cell_data[[ctype_col]] == src_type, c(xc, yc)]
    target <- cell_data[cell_data[[ctype_col]] == tgt_type, c(xc, yc)]
    if (nrow(target) == 0 || nrow(source) == 0) return(rep(NA, nrow(source)))
    dist_matrix <- sqrt(
      outer(source[[xc]], target[[xc]], "-")^2 + 
      outer(source[[yc]], target[[yc]], "-")^2
    )
    apply(dist_matrix, 1, min)
  }
  
  # --- 内部辅助函数 2: 处理单个组 ---
  .process_single_group <- function(group_df, src, tgt, ctype_col, xc, yc) {
    source_indices <- which(group_df[[ctype_col]] == src)
    group_df$min_dist <- NA
    if (length(source_indices) > 0) {
      group_df$min_dist[source_indices] <- .calculate_min_distance(
        group_df, src, tgt, ctype_col, xc, yc
      )
    }
    return(group_df)
  }
  
  # --- 内部辅助函数 3: 计算单个配对 (已修改) ---
  .compute_single_pair <- function(input_df, src, tgt, grp_col, ctype_col, xc, yc) {
    
    result_df <- input_df %>%
      dplyr::filter(!!sym(ctype_col) == src | !!sym(ctype_col) == tgt) %>%
      #
      # !! 这一行已被移除: tibble::rownames_to_column("barcode") !!
      #
      dplyr::group_by(!!sym(grp_col)) %>%
      dplyr::group_modify(~ .process_single_group(.x, 
                                                  src = src, 
                                                  tgt = tgt, 
                                                  ctype_col = ctype_col, 
                                                  xc = xc, 
                                                  yc = yc)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(min_dist))
      
    if(nrow(result_df) > 0) {
      result_df$source_type <- src
      result_df$target_type <- tgt
    }
    return(result_df)
  }
  
  # --- 函数主逻辑：检查模式 ---
  
  # 检查列名时，我们不关心 'barcode'，只关心用于计算的列
  req_cols <- c(group_col, celltype_col, x_col, y_col)
  missing_cols <- setdiff(req_cols, colnames(input_df))
  if(length(missing_cols) > 0) {
    stop("Input data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (!is.null(source_type) && !is.null(target_type)) {
    # 模式 1: 指定了源和目标
    message(paste("Calculating distances for single pair:", source_type, "->", target_type))
    return(.compute_single_pair(input_df, source_type, target_type, group_col, celltype_col, x_col, y_col))
    
  } else if (is.null(source_type) && is.null(target_type)) {
    # 模式 2: 未指定源和目标 -> 运行所有配对
    message("Calculating distances for all pairwise combinations...")
    
    cell_types <- unique(input_df[[celltype_col]])
    if (length(cell_types) < 2) {
      warning("Less than 2 cell types found. No pairs to calculate.")
      return(data.frame())
    }
    
    pairs_df <- tidyr::expand_grid(source = cell_types, target = cell_types) %>%
      dplyr::filter(source != target)
      
    total_pairs <- nrow(pairs_df)
    message(paste("Found", total_pairs, "pairs to process."))
    
    # --- 使用 For 循环的逻辑 ---
    
    results_list <- list()
    message("Running in sequential 'for' loop mode...")
    
    for (i in 1:total_pairs) {
      current_source <- pairs_df$source[i]
      current_target <- pairs_df$target[i]
      
      message(paste0("Processing pair ", i, " of ", total_pairs, ": ", 
                     current_source, " -> ", current_target))
      
      pair_result <- .compute_single_pair(
        input_df = input_df, 
        src = current_source, 
        tgt = current_target, 
        grp_col = group_col, 
        ctype_col = celltype_col, 
        xc = x_col, 
        yc = y_col
      )
      
      results_list[[i]] <- pair_result
    }
    
    message("All pairs processed. Combining results...")
    all_results <- dplyr::bind_rows(results_list)
    
    return(all_results)
    
  } else {
    stop("Both 'source_type' and 'target_type' must be provided, or both must be NULL.")
  }
}