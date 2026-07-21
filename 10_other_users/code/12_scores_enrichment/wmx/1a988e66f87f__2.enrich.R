# 这个是为了对上次富集的结果进行挑选通路的绘制
# 改法：
# 1. 从 matchedfile 读入通路表
# 2. 用 matchedfile$ID 和 enrich 结果的 ID 严格匹配
# 3. 每个 PT + group comparison 的 GO + KEGG 合并后单独绘图

library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(readr)
library(forcats)

rm(list = ls())

# =========================
# 1. 输入输出路径
# =========================
enrich_dir <- "/mnt/NexData_file/home/wmx/Project/AE2512227671-01/shouhou8/Diff_and_enrich/Enrich"

out_base <- "/mnt/NexData_file/home/wmx/Project/AE2512227671-01/shouhou8/Diff_and_enrich/Enrich/enrich_plot/"
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

# matchedfile：这里换成你自己保存的 matched_pathways_from_gmt_with_ID.csv
matchedfile <- "/mnt/NexData_file/home/wmx/Project/AE2512227671-01/shouhou8/matched_pathways_from_gmt.csv"

# enrich 文件
efiles <- list.files(
  enrich_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = FALSE
)

# =========================
# 2. 读取 matchedfile
# =========================
matched_df_all <- read.csv(
  matchedfile,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# 必须包含 main_name 和 ID
if (!all(c("main_name", "ID") %in% colnames(matched_df_all))) {
  stop("matchedfile 必须包含 main_name 和 ID 两列")
}

matched_df_all <- matched_df_all %>%
  mutate(
    main_name = as.character(main_name),
    ID = trimws(as.character(ID))
  ) %>%
  filter(
    !is.na(main_name),
    main_name != "",
    !is.na(ID),
    ID != ""
  )

# 去重
matched_df_all <- matched_df_all %>%
  distinct(main_name, ID, .keep_all = TRUE)

print(table(matched_df_all$main_name))
print(head(matched_df_all))

# =========================
# 3. 工具函数
# =========================
safe_name <- function(x) {
  x <- gsub("[/ ]", "_", x)
  x <- gsub("α", "alpha", x)
  x <- gsub("γ", "gamma", x)
  x <- gsub("[()]", "", x)
  x <- gsub("[,:;|]", "_", x)
  x <- gsub("__+", "_", x)
  x
}

# 从文件名提取 PT + group comparison
# 例如：
# PT1_Slc5a12_C-K_enrich.go.csv      -> PT1_Slc5a12_C-K
# PT1_Slc5a12_C-K_enrich.KEGG.csv    -> PT1_Slc5a12_C-K
# PT7_Atp11a_KFG-K_enrich.go.csv     -> PT7_Atp11a_KFG-K
get_comparison_name <- function(file) {
  x <- basename(file)
  
  x <- sub("\\.csv$", "", x, ignore.case = TRUE)
  x <- sub("\\.xls$", "", x, ignore.case = TRUE)
  x <- sub("\\.txt$", "", x, ignore.case = TRUE)
  
  x <- sub("_enrich\\.go.*$", "", x, ignore.case = TRUE)
  x <- sub("_enrich\\.kegg.*$", "", x, ignore.case = TRUE)
  x <- sub("_enrichGO.*$", "", x, ignore.case = TRUE)
  x <- sub("_enrichKEGG.*$", "", x, ignore.case = TRUE)
  
  return(x)
}

# =========================
# 4. 读取单个 enrichGO / enrichKEGG 文件
# =========================
read_one_enrich_file <- function(ff) {
  
  message("[READ] ", ff)
  
  if (!file.exists(ff)) {
    message("[SKIP] file not exists: ", ff)
    return(NULL)
  }
  
  df <- tryCatch({
    if (grepl("\\.csv$", ff, ignore.case = TRUE)) {
      read.csv(ff, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      read.delim(ff, check.names = FALSE, stringsAsFactors = FALSE)
    }
  }, error = function(e) {
    message("[SKIP] read failed: ", ff)
    return(NULL)
  })
  
  if (is.null(df)) return(NULL)
  
  if (nrow(df) == 0) {
    message("[SKIP] empty enrichment result: ", ff)
    return(NULL)
  }
  
  # 所有列先转字符，避免 GO/KEGG bind_rows 类型冲突
  df <- df %>%
    mutate(across(everything(), as.character))
  
  # 兼容 Description / Term / Pathway
  if (!"Description" %in% colnames(df)) {
    if ("Term" %in% colnames(df)) {
      df$Description <- df$Term
    } else if ("Pathway" %in% colnames(df)) {
      df$Description <- df$Pathway
    } else {
      message("[SKIP] no Description/Term/Pathway: ", ff)
      return(NULL)
    }
  }
  
  df <- df[!is.na(df$Description) & df$Description != "", ]
  
  if (nrow(df) == 0) {
    message("[SKIP] no valid Description rows: ", ff)
    return(NULL)
  }
  
  # 清理 KEGG 物种尾巴
  df$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", df$Description)
  
  # 必须有 ID，用于 matchedfile 匹配
  if (!"ID" %in% colnames(df)) {
    message("[SKIP] no ID column: ", ff)
    return(NULL)
  }
  
  df$ID <- trimws(as.character(df$ID))
  df <- df[!is.na(df$ID) & df$ID != "", ]
  
  if (nrow(df) == 0) {
    message("[SKIP] no valid ID rows: ", ff)
    return(NULL)
  }
  
  # 兼容 p.adjust
  if (!"p.adjust" %in% colnames(df)) {
    if ("Adjusted.P.value" %in% colnames(df)) {
      df$`p.adjust` <- df$Adjusted.P.value
    } else if ("FDR" %in% colnames(df)) {
      df$`p.adjust` <- df$FDR
    } else if ("qvalue" %in% colnames(df)) {
      df$`p.adjust` <- df$qvalue
    } else if ("pvalue" %in% colnames(df)) {
      df$`p.adjust` <- df$pvalue
    } else {
      message("[SKIP] no p.adjust/FDR/qvalue/pvalue: ", ff)
      return(NULL)
    }
  }
  
  df$`p.adjust` <- suppressWarnings(as.numeric(df$`p.adjust`))
  
  if ("pvalue" %in% colnames(df)) {
    df$pvalue <- suppressWarnings(as.numeric(df$pvalue))
  }
  
  if ("qvalue" %in% colnames(df)) {
    df$qvalue <- suppressWarnings(as.numeric(df$qvalue))
  }
  
  # geneID 转字符
  if ("geneID" %in% colnames(df)) {
    df$geneID <- as.character(df$geneID)
  }
  
  # 兼容 Count
  if (!"Count" %in% colnames(df)) {
    if ("geneID" %in% colnames(df)) {
      df$Count <- sapply(strsplit(as.character(df$geneID), "/"), length)
    } else if ("Overlap" %in% colnames(df)) {
      df$Count <- as.numeric(sub("/.*", "", df$Overlap))
    } else {
      df$Count <- 1
    }
  }
  
  df$Count <- suppressWarnings(as.numeric(df$Count))
  
  # 判断数据库来源
  b <- basename(ff)
  df$database <- dplyr::case_when(
    grepl("kegg", b, ignore.case = TRUE) ~ "KEGG",
    grepl("go", b, ignore.case = TRUE) ~ "GO",
    TRUE ~ "Enrich"
  )
  
  return(df)
}

# =========================
# 5. 用 matchedfile$ID 严格匹配
# =========================
match_one_pathway_type <- function(enrich_df, matched_df_all, main_name) {
  
  if (!"ID" %in% colnames(enrich_df)) {
    message("[SKIP] enrich_df has no ID column")
    return(data.frame())
  }
  
  target_df <- matched_df_all %>%
    filter(main_name == !!main_name) %>%
    mutate(ID_match = trimws(as.character(ID))) %>%
    filter(!is.na(ID_match), ID_match != "")
  
  if (nrow(target_df) == 0) {
    message("[SKIP] no target ID for ", main_name)
    return(data.frame())
  }
  
  enrich_df <- enrich_df %>%
    mutate(ID_match = trimws(as.character(ID)))
  
  # 为了避免 join 后列名冲突，只保留 matchedfile 中需要附加的列
  add_cols <- c("ID_match", "main_name")
  if ("keyword" %in% colnames(target_df)) add_cols <- c(add_cols, "keyword")
  if ("pathway" %in% colnames(target_df)) add_cols <- c(add_cols, "pathway")
  if ("gmt_info" %in% colnames(target_df)) add_cols <- c(add_cols, "gmt_info")
  
  target_df2 <- target_df %>%
    select(all_of(add_cols)) %>%
    distinct()
  
  hit <- enrich_df %>%
    inner_join(target_df2, by = "ID_match")
  
  if (nrow(hit) == 0) {
    return(data.frame())
  }
  
  # 避免重复
  hit <- hit %>%
    distinct(database, ID, Description, main_name, .keep_all = TRUE)
  
  return(hit)
}

# =========================
# 6. barplot
# =========================
plot_enrich_bar <- function(df, title, out_file, showCategory = 10) {
  
  if (nrow(df) == 0) return(NULL)
  
  plot_df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = showCategory) %>%
    mutate(
      Description = fct_reorder(Description, Count)
    )
  
  p <- ggplot(plot_df, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_col(width = 0.75) +
    scale_fill_gradient(low = "#ff4757", high = "#546de5") +
    labs(
      title = title,
      x = "Count",
      y = NULL,
      fill = "p.adjust"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 12)
    )
  
  ggsave(out_file, p, width = 8, height = max(4, 0.38 * nrow(plot_df) + 2))
  
  out_png <- gsub("\\.pdf$", ".png", out_file)
  ggsave(out_png, p, width = 8, height = max(4, 0.38 * nrow(plot_df) + 2), dpi = 300)
  
  message("[SAVE] ", out_file)
}

# =========================
# 7. dotplot
# =========================
plot_enrich_dot <- function(df, title, out_file, showCategory = 10) {
  
  if (nrow(df) == 0) return(NULL)
  
  plot_df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = showCategory) %>%
    mutate(
      Description = fct_reorder(Description, Count)
    )
  
  p <- ggplot(plot_df, aes(x = Count, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#ff4757", high = "#546de5") +
    labs(
      title = title,
      x = "Count",
      y = NULL,
      color = "p.adjust",
      size = "Count"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 12)
    )
  
  ggsave(out_file, p, width = 8, height = max(4, 0.38 * nrow(plot_df) + 2))
  
  out_png <- gsub("\\.pdf$", ".png", out_file)
  ggsave(out_png, p, width = 8, height = max(4, 0.38 * nrow(plot_df) + 2), dpi = 300)
  
  message("[SAVE] ", out_file)
}

# =========================
# 8. 一个 PT + group comparison 的 GO + KEGG 合并后画图
# =========================
process_one_comparison <- function(
    comp,
    comp_files,
    matched_df_all,
    out_base,
    showCategory = 10
) {
  
  message("========== Processing comparison: ", comp, " ==========")
  
  enrich_list <- list()
  
  for (ff in comp_files) {
    tmp <- read_one_enrich_file(ff)
    if (!is.null(tmp) && nrow(tmp) > 0) {
      enrich_list[[ff]] <- tmp
    }
  }
  
  if (length(enrich_list) == 0) {
    message("[SKIP] no valid enrich files for ", comp)
    return(NULL)
  }
  
  # 再保险：合并前统一类型
  enrich_list <- lapply(enrich_list, function(x) {
    x <- x %>%
      mutate(across(everything(), as.character))
    
    x$`p.adjust` <- suppressWarnings(as.numeric(x$`p.adjust`))
    x$Count <- suppressWarnings(as.numeric(x$Count))
    
    return(x)
  })
  
  enrich_df <- bind_rows(enrich_list)
  
  if (nrow(enrich_df) == 0) {
    message("[SKIP] merged enrichment result is empty for ", comp)
    return(NULL)
  }
  
  enrich_df$ID <- trimws(as.character(enrich_df$ID))
  
  comp_dir <- file.path(out_base, safe_name(comp))
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 保存当前组别 GO + KEGG 合并表
  write.csv(
    enrich_df,
    file.path(comp_dir, paste0(safe_name(comp), "_GO_KEGG_merged.csv")),
    row.names = FALSE
  )
  
  # 按 matchedfile 中的 main_name 循环
  for (main_name in unique(matched_df_all$main_name)) {
    
    matched_df <- match_one_pathway_type(
      enrich_df = enrich_df,
      matched_df_all = matched_df_all,
      main_name = main_name
    )
    
    if (nrow(matched_df) == 0) {
      message("[SKIP] ", comp, " - ", main_name, ": no matched pathway by ID")
      next
    }
    
    # 同一个 ID 可能重复，保留 p.adjust 最小的一条
    matched_df <- matched_df %>%
      group_by(ID) %>%
      slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(p.adjust)
    
    # 保存匹配结果
    write.csv(
      matched_df,
      file.path(comp_dir, paste0(safe_name(comp), "_", safe_name(main_name), "_matched_by_ID.csv")),
      row.names = FALSE
    )
    
    bar_file <- file.path(
      comp_dir,
      paste0(safe_name(comp), "_", safe_name(main_name), "_barplot.pdf")
    )
    
    dot_file <- file.path(
      comp_dir,
      paste0(safe_name(comp), "_", safe_name(main_name), "_dotplot.pdf")
    )
    
    plot_enrich_bar(
      df = matched_df,
      title = paste0(comp, " | ", main_name, "_barplot"),
      out_file = bar_file,
      showCategory = showCategory
    )
    
    plot_enrich_dot(
      df = matched_df,
      title = paste0(comp, " | ", main_name, "_dotplot"),
      out_file = dot_file,
      showCategory = showCategory
    )
  }
}

# =========================
# 9. 按 PT + group comparison 分组
# =========================
file_info <- data.frame(
  file = efiles,
  comparison = sapply(efiles, get_comparison_name),
  stringsAsFactors = FALSE
)

print(table(file_info$comparison))

# 保存 file_info，方便检查 GO/KEGG 是否配对成功
write.csv(
  file_info,
  file.path(out_base, "file_info_comparison_grouping.csv"),
  row.names = FALSE
)

# =========================
# 10. 批量运行
# =========================
for (comp in unique(file_info$comparison)) {
  
  comp_files <- file_info$file[file_info$comparison == comp]
  
  message("[FILES for ", comp, "]")
  print(comp_files)
  
  process_one_comparison(
    comp = comp,
    comp_files = comp_files,
    matched_df_all = matched_df_all,
    out_base = out_base,
    showCategory = 10
  )
}