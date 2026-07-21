# log_con <- file(snakemake@log[[1]], "w")
# sink(log_con, type = "message")

library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(tidyr)
library(pheatmap)
library(ggpubr)
library(plotly)
library(RColorBrewer)
library(jsonlite)

library(ape)

##################### parameter #########################################
enrich.out.dir <- snakemake@params[["input_dir"]]
output.dir <- snakemake@params[["out_dir"]]

##################### function defined #########################################
create.json <- function(mat_for_heatmap, annot, pheatmap_obj, title = "") {
  ## Row tree from pheatmap clustering
  row_order <- pheatmap_obj$tree_row$order
  mat_ord <- mat_for_heatmap[row_order, , drop = FALSE]
  mat_num <- round(mat_ord, 4)
  mat_num[is.na(mat_num)] <- 0

  row_labels <- rownames(mat_ord)
  col_labels <- colnames(mat_ord)

  ## [0] Row tree dataset & newick
  row_source <- rbind(c("name"), matrix(row_labels, ncol = 1))
  row_newick <- ape::write.tree(ape::as.phylo(pheatmap_obj$tree_row))

  ## [1] Heatmap dataset
  heatmap_source <- rbind(
    c("", col_labels),
    cbind(row_labels, mat_num)
  )

  ## [2] Compare annotation bar (key-value source object)
  compare_source <- as.list(setNames(
    as.character(annot[col_labels, "Compare"]),
    col_labels
  ))

  ## [3] Regulation annotation bar (key-value source object)
  regulation_source <- as.list(setNames(
    as.character(annot[col_labels, "Regulation"]),
    col_labels
  ))

  ## Generate annotation colors (matching pheatmap defaults)
  compare_vals <- unique(as.character(annot[col_labels, "Compare"]))
  n_compare <- length(compare_vals)
  if (n_compare <= 9) {
    compare_colors <- brewer.pal(max(3, n_compare), "Set1")[1:n_compare]
  } else {
    compare_colors <- colorRampPalette(brewer.pal(9, "Set1"))(n_compare)
  }
  compare_color_map <- as.list(setNames(compare_colors, compare_vals))

  regulation_vals <- unique(as.character(annot[col_labels, "Regulation"]))
  reg_colors <- sapply(regulation_vals, function(v) {
    if (grepl("up", v, ignore.case = TRUE)) {
      return("#E41A1C")
    }
    if (grepl("down", v, ignore.case = TRUE)) {
      return("#377EB8")
    }
    return("#999999")
  })
  reg_color_map <- as.list(setNames(reg_colors, regulation_vals))

  ## Build JSON object
  json_obj <- list(
    width = 1000,
    height = 1000,
    margin = list(top = 80, right = 60, bottom = 60, left = 60),
    dataset = list(
      list(source = row_source, categories = row_newick),
      list(source = heatmap_source),
      list(source = compare_source),
      list(source = regulation_source)
    ),
    series = list(
      list(
        type = "tree",
        datasetIndex = 0,
        gridIndex = 0,
        direction = "horizontal",
        encode = list(tree = "categories"),
        label = list(show = FALSE)
      ),
      list(
        type = "heatmap",
        xAxisIndex = 1,
        yAxisIndex = 1,
        datasetIndex = 1,
        gridIndex = 1,
        legendIndex = 1,
        seriesLayoutBy = "row",
        label = list(show = FALSE),
        annotationBars = list(
          list(datasetIndex = 2, legendIndex = 2),
          list(datasetIndex = 3, legendIndex = 3)
        )
      )
    ),
    legend = list(
      list(show = FALSE),
      list(right = "15%", title = "", colors = list("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")),
      list(right = "10%", title = "Compare", colors = compare_color_map),
      list(top = "20%", right = "10%", title = "Regulation", colors = reg_color_map)
    ),
    grid = list(
      list(width = "8%"),
      list(left = "9%", width = "45%")
    ),
    xAxis = list(
      list(type = "scaleLinear", silent = TRUE, show = FALSE),
      list(type = "scaleBand", show = FALSE)
    ),
    yAxis = list(
      list(type = "scaleLinear", silent = TRUE, show = FALSE),
      list(type = "scaleBand", position = "right", show = TRUE)
    )
  )
  toJSON(json_obj, auto_unbox = TRUE, pretty = TRUE)
}
enrich.cluster <- function(enrich.data, topn = 40, th.pvalue = 0.05) {
  data <- subset(enrich.data, select = c("Compare", "Terms", "Regulation", "p.value"))
  if (nrow(data) < 3) {
    return(NULL)
  }
  make_significance_stars <- function(pvalue_score) {
    significance_stars <- matrix("", nrow = nrow(pvalue_score), ncol = ncol(pvalue_score), dimnames = dimnames(pvalue_score))
    significance_stars[pvalue_score > -log10(0.05)] <- "*"
    significance_stars[pvalue_score > -log10(0.01)] <- "**"
    significance_stars[pvalue_score > -log10(0.001)] <- "***"
    significance_stars
  }
  data$p.value <- -log10(data$p.value + 1e-276)
  # 去重：同一 Terms+Regulation+Compare 保留 p.value 最大的（即最显著的原始p值）
  data <- data.table(data)
  data <- data[, .SD[which.max(p.value)], by = .(Terms, Regulation, Compare)]
  data[, maxv := max(p.value), by = "Terms"]
  data.wid <- data.table::dcast(data, Terms + maxv ~ Regulation + Compare, value.var = "p.value", fill = 0, sep = ":") %>%
    arrange(-maxv) %>%
    filter(nchar(as.vector(Terms)) < 50) %>%
    filter(maxv > -log10(th.pvalue)) %>%
    head(topn)
  if (ncol(data.wid) < 4) {
    return(NULL)
  }
  if (nrow(data.wid) < 3) {
    return(NULL)
  }
  annot <- data.frame(ID = colnames(data.wid)[c(-1, -2)])
  rownames(annot) <- annot$ID
  annot <- separate(annot, col = "ID", sep = ":", into = c("Regulation", "Compare"))
  data.wid <- as.data.frame(data.wid)
  rownames(data.wid) <- data.wid$Terms
  scalev <- "row"
  if (ncol(data.wid) == 4) {
    scalev <- "none"
  }
  heatmap_color <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
  pvalue_score <- as.matrix(data.wid[, c(-1, -2), drop = FALSE])
  rownames(pvalue_score) <- data.wid$Terms
  significance_stars <- make_significance_stars(pvalue_score)
  p <- pheatmap(pvalue_score,
    scale = scalev, annotation = annot, silent = T, cluster_cols = F, border_color = "white",
    main = enrich.data$Category[1], show_colnames = F, color = heatmap_color,
    display_numbers = significance_stars, number_color = "black", fontsize_number = 12
  )
  mat_final <- pvalue_score
  if (scalev == "row") {
    mat_final <- t(scale(t(mat_final)))
    mat_final[is.na(mat_final)] <- 0
  }
  p.json <- create.json(mat_final, annot, p, title = enrich.data$Category[1])
  return(list(plot = p, json = p.json, data.wid = data.wid, annot = annot))
}


##################### Execute the main program #########################################
# 空数据占位图辅助函数
save_empty_plot <- function(png_path, msg = "No data available") {
  png(png_path, width = 800, height = 400)
  plot.new()
  text(0.5, 0.5, msg, cex = 2, col = "grey50")
  dev.off()
  pdf_path <- sub(png_path, pattern = "\\.png$", replacement = ".pdf")
  pdf(pdf_path, width = 8, height = 4)
  plot.new()
  text(0.5, 0.5, msg, cex = 2, col = "grey50")
  dev.off()
  json_path <- sub(png_path, pattern = "\\.png$", replacement = ".json")
  writeLines("{}", json_path)
}

files <- dir(enrich.out.dir, pattern = "Table 4. Functional enrichment of ", full.names = T)
file.list <- setNames(files, sub(sub(files, pattern = ".*Table 4. Functional enrichment of ", replacement = ""), pattern = ".xls$", replacement = ""))

# 安全读取：跳过空文件
data.list <- lapply(file.list, function(f) {
  tryCatch(
    {
      d <- fread(f)
      if (is.null(d) || nrow(d) == 0) {
        return(NULL)
      }
      d
    },
    error = function(e) NULL
  )
})
data.list <- Filter(Negate(is.null), data.list)

if (!dir.exists(output.dir)) dir.create(output.dir, recursive = TRUE)

# 所有可能的分类
all_categories <- c("GOBP", "GOCC", "GOMF", "KEGG", "PFAM")

if (length(data.list) > 0) {
  data <- rbindlist(data.list, idcol = "Compare", fill = TRUE)

  categories_in_data <- unique(data$Category)

  lapply(split(data, data$Category), function(x) {
    res <- enrich.cluster(x, topn = 40, th.pvalue = 0.05)
    i <- x$Category[1]
    if (!is.null(res)) {
      data.wid <- res$data.wid
      n_rows <- nrow(data.wid)
      n_cols <- ncol(data.wid) - 2
      plot_height <- 2 + n_rows * 0.2
      max_label_len <- max(nchar(as.character(rownames(data.wid))), na.rm = TRUE)
      if (is.infinite(max_label_len)) max_label_len <- 10
      plot_width <- 5 + n_cols * 0.3 + max_label_len * 0.05
      ggsave(res$plot$gtable, filename = file.path(output.dir, paste0("Figure 5. Enrichment cluster of ", i, ".pdf")), width = plot_width, height = plot_height, create.dir = T)
      ggsave(res$plot$gtable, filename = file.path(output.dir, paste0("Figure 5. Enrichment cluster of ", i, ".png")), width = plot_width, height = plot_height)
      writeLines(res$json, file.path(output.dir, paste0("Figure 5. Enrichment cluster of ", i, ".json")))
    } else {
      save_empty_plot(file.path(output.dir, paste0("Figure 5. Enrichment cluster of ", i, ".png")),
        msg = paste0("Insufficient data for ", i, " enrichment cluster")
      )
    }
  })

  # 为缺失的分类生成占位图
  missing_categories <- setdiff(all_categories, categories_in_data)
} else {
  missing_categories <- all_categories
}

for (cat in missing_categories) {
  save_empty_plot(file.path(output.dir, paste0("Figure 5. Enrichment cluster of ", cat, ".png")),
    msg = paste0("No enrichment data for ", cat, " cluster")
  )
}

write.csv("", file.path(output.dir, "done"))
