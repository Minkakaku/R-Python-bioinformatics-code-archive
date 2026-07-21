#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(stringr)
library(optparse)
library(readr)
library(readxl)

# 读取兼容函数
read_table_compat <- function(file_path, delim = "\t", col_names = TRUE) {
  ext <- tolower(tools::file_ext(file_path))
  if (ext == "csv") {
    data <- readr::read_csv(file_path, show_col_types = FALSE, col_names = col_names)
  } else if (ext == "tsv") {
    data <- readr::read_tsv(file_path, show_col_types = FALSE, col_names = col_names)
  } else if (ext == "xlsx" || ext == "xls") {
    data <- tryCatch(
      readxl::read_excel(file_path, col_names = col_names),
      error = function(e) {
        readr::read_delim(file_path, delim = delim, show_col_types = FALSE, col_names = col_names)
      }
    )
  } else {
    data <- readr::read_delim(file_path, delim = delim, show_col_types = FALSE, col_names = col_names)
  }
  as.data.frame(data)
}

# 参数解析
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "输入文件路径 (xls)"),
  make_option(c("--width"), type = "double", default = 10, help = "图片宽度"),
  make_option(c("--height"), type = "double", default = 8, help = "图片高度"),
  make_option(c("--top_n"), type = "integer", default = 20, help = "整体展示的Term数量"),
  make_option(c("--color_var"), type = "character", default = "p.adjust", help = "用于颜色的变量 (p.adjust 或 pvalue)"),
  make_option(c("--x_title"), type = "character", default = "Count", help = "X轴标题内容"),
  make_option(c("--x_title_size"), type = "integer", default = 14, help = "X轴标题字体大小"),
  make_option(c("--x_text_size"), type = "integer", default = 10, help = "X轴刻度字体大小"),
  make_option(c("--y_text_size"), type = "integer", default = 10, help = "Y轴刻度字体大小"),
  make_option(c("--font_family"), type = "character", default = "Arial", help = "字体 (Arial 或 Times New Roman)"),
  make_option(c("--legend_title_size"), type = "integer", default = 12, help = "图例标题字体大小"),
  make_option(c("--legend_text_size"), type = "integer", default = 10, help = "图例文字字体大小"),
  make_option(c("--low_color"), type = "character", default = "#E06663", help = "低值颜色 (显著)"),
  make_option(c("--high_color"), type = "character", default = "#327EBA", help = "高值颜色 (不显著)"),
  make_option(c("--plot_name"), type = "character", default = "enrich_plot", help = "生成图名称(包含路径，无后缀)"),
  make_option(c("--outdir"), type = "character", default = "./", help = "输出目录"),
  make_option(c("--coord_flip"), type = "character", default = "true", help = "是否旋转坐标轴 (true/false)"),
  make_option(c("--y_label_wrap"), type = "integer", default = 80, help = "Y轴标签换行长度"),
  make_option(c("--bar_width"), type = "double", default = 0.8, help = "柱子宽度"),
  make_option(c("--term_list"), type = "character", default = NULL, help = "分号分隔的Description列表")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("必须提供输入文件路径", call. = FALSE)
}

# 读取数据
data <- read_table_compat(opt$input)

# 数据预处理
if (!is.null(opt$term_list) && opt$term_list != "-") {
  target_terms <- unlist(strsplit(opt$term_list, ";"))
  target_terms <- trimws(target_terms)
  
  existing_terms <- unique(data$Description)
  missing_terms <- setdiff(target_terms, existing_terms)
  
  if (length(missing_terms) > 0) {
    warning(paste("以下Terms在数据中未找到:", paste(missing_terms, collapse = ", ")))
  }
  
  plot_data <- data %>%
    filter(Description %in% target_terms) %>%
    arrange(!!sym(opt$color_var))
} else {
  # 按照p.adjust排序，取整体前top_n
  plot_data <- data %>%
    arrange(!!sym(opt$color_var)) %>%
    slice_head(n = opt$top_n)
}

# 设置Factor以保持顺序
# 如果是 coord_flip=true (默认), Description在Y轴, ggplot从下到上. rev后, 显著的在上.
# 如果是 coord_flip=false, Description在X轴, ggplot从左到右. rev后, 显著的在左.
plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))

# 字体映射
font_map <- list(
  "Arial" = "sans",
  "Times New Roman" = "serif"
)
plot_font <- ifelse(opt$font_family %in% names(font_map), font_map[[opt$font_family]], opt$font_family)

# 绘图
if (opt$coord_flip == "true" || opt$coord_flip == "yes") {
  # 水平条形图 (Description 在 Y 轴)
  # 使用 coord_flip() 翻转
  p <- ggplot(plot_data, aes(x = Description, y = Count, fill = !!sym(opt$color_var))) +
    geom_bar(stat = "identity", width = opt$bar_width) +
    coord_flip() +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = opt$y_label_wrap)) +
    labs(y = opt$x_title, x = NULL) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = opt$x_text_size, color = "black"), # 翻转后 X 变为 Y (Description)
      axis.text.x = element_text(size = opt$y_text_size, color = "black"), # 翻转后 Y 变为 X (Count)
      axis.title.x = element_text(size = opt$x_title_size)
    )
} else {
  # 垂直柱状图 (Description 在 X 轴)
  p <- ggplot(plot_data, aes(x = Description, y = Count, fill = !!sym(opt$color_var))) +
    geom_bar(stat = "identity", width = opt$bar_width) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = opt$y_label_wrap)) +
    labs(y = opt$x_title, x = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = opt$x_text_size, color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = opt$y_text_size, color = "black"),
      axis.title.y = element_text(size = opt$x_title_size)
    )
}

# 公共主题
p <- p + 
  scale_fill_gradient(low = opt$low_color, high = opt$high_color) +
  theme(
    text = element_text(family = plot_font),
    legend.title = element_text(size = opt$legend_title_size),
    legend.text = element_text(size = opt$legend_text_size)
  )

# 保存图片
plot_dir <- dirname(opt$plot_name)
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
output_png <- paste0(opt$plot_name, ".png")
output_pdf <- paste0(opt$plot_name, ".pdf")
ggsave(output_png, plot = p, width = opt$width, height = opt$height, limitsize = FALSE)
ggsave(output_pdf, plot = p, width = opt$width, height = opt$height, limitsize = FALSE)

if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}
file.create(file.path(opt$outdir, "done"))
