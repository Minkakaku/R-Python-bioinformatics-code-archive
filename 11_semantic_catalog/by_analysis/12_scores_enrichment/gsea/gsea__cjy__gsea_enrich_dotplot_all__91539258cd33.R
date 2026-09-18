
library(optparse)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)
library(dplyr)

# 1. 解析命令行参数
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "输入GSEA结果RData文件路径"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "输出目录 [默认: .]"),
  make_option(c("-n", "--plot_name"), type = "character", default = "gsea_dotplot",
              help = "输出文件名(不含后缀) [默认: gsea_dotplot]"),
  make_option(c("-w", "--width"), type = "double", default = 13,
              help = "图片宽度 [默认: 13]"),
  make_option(c("-H", "--height"), type = "double", default = 10,
              help = "图片高度 [默认: 10]"),
  make_option(c("--font_family"), type = "character", default = "Arial",
              help = "字体族 (Arial 或 Times New Roman) [默认: Arial]"),
  make_option(c("-t", "--terms"), type = "character", default = NULL,
              help = "指定要展示的通路Description，分号分隔。若未指定，则使用top_n [可选]"),
  make_option(c("--top_n"), type = "integer", default = 30,
              help = "展示Top N通路 [默认: 30]"),
  make_option(c("--color_var"), type = "character", default = "p.adjust",
              help = "颜色映射变量 (p.adjust 或 pvalue) [默认: p.adjust]"),
  make_option(c("--low_color"), type = "character", default = "#E06663",
              help = "低值颜色 [默认: #E06663]"),
  make_option(c("--high_color"), type = "character", default = "#327EBA",
              help = "高值颜色 [默认: #327EBA]"),
  make_option(c("--x_title"), type = "character", default = "GeneRatio",
              help = "X轴标题 [默认: GeneRatio]"),
  make_option(c("--x_title_size"), type = "double", default = 16,
              help = "X轴标题字体大小 [默认: 16]"),
  make_option(c("--x_text_size"), type = "double", default = 14,
              help = "X轴刻度字体大小 [默认: 14]"),
  make_option(c("--y_text_size"), type = "double", default = 14,
              help = "Y轴刻度字体大小 [默认: 14]"),
  make_option(c("--y_label_wrap"), type = "integer", default = 80,
              help = "Y轴标签换行字符数 [默认: 80]"),
  make_option(c("--legend_title_size"), type = "double", default = 12,
              help = "图例标题字体大小 [默认: 12]"),
  make_option(c("--legend_text_size"), type = "double", default = 10,
              help = "图例文本字体大小 [默认: 10]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("必须提供输入文件 (-i/--input)", call. = FALSE)
}

# 创建输出目录
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# 2. 加载数据函数 (兼容多种对象)
load_gsea_data <- function(file_path) {
  cat("加载GSEA结果数据:", file_path, "\n")
  env <- new.env()
  load(file_path, envir = env)
  
  # 优先查找标准命名
  candidate_names <- c("GO_gseresult", "KEGG_gseresult", "Reactome_gseresult")
  found_name <- NULL
  
  for (name in candidate_names) {
    if (exists(name, envir = env)) {
      found_name <- name
      break
    }
  }
  
  # 如果未找到标准命名，扫描所有对象查找 gseaResult 类型
  if (is.null(found_name)) {
    all_objs <- ls(envir = env)
    for (obj_name in all_objs) {
      obj <- get(obj_name, envir = env)
      if (inherits(obj, "gseaResult")) {
        found_name <- obj_name
        cat("提示: 未找到标准命名对象，使用发现的 gseaResult 对象:", found_name, "\n")
        break
      }
    }
  }
  
  if (is.null(found_name)) {
    stop("错误: RData文件中未找到 gseaResult 类型的对象。请检查输入文件。", call. = FALSE)
  }
  
  gsea_result <- get(found_name, envir = env)
  return(gsea_result)
}

# 3. 主逻辑
gsea_obj <- load_gsea_data(opt$input)

# 确定要展示的通路
show_category <- opt$top_n # 默认情况

if (!is.null(opt$terms) && opt$terms != "-" && opt$terms != "") {
  # 解析terms
  terms_vec <- unlist(strsplit(opt$terms, ";"))
  terms_vec <- trimws(terms_vec)
  terms_vec <- terms_vec[terms_vec != ""]
  
  if (length(terms_vec) > 0) {
    # 匹配Description
    # gseaResult@result 包含 ID, Description 等列
    # 我们需要传递 Description 或 ID 给 dotplot 的 showCategory 参数
    # enrichplot::dotplot 这里的 showCategory 如果是数字则是top N，如果是向量则是名字
    # 对于 gseaResult, dotplot 似乎更倾向于处理 top N 或 自动。
    # 如果要指定特定的 Description，最好先过滤。
    
    # 检查 terms 是否在结果中
    matched_idx <- which(gsea_obj@result$Description %in% terms_vec)
    if (length(matched_idx) == 0) {
      warning("警告: 指定的 Terms 未在 GSEA 结果中找到。将回退到使用 Top N。", call. = FALSE)
    } else {
      show_category <- terms_vec
      cat("将展示指定的", length(matched_idx), "个 Terms。\n")
    }
  }
}

# 字体映射
font_map <- list(
  "Arial" = "sans",
  "Times New Roman" = "serif"
)
plot_font <- ifelse(opt$font_family %in% names(font_map), font_map[[opt$font_family]], opt$font_family)

# 绘图
cat("开始绘制 Dotplot...\n")

# 处理颜色变量
color_var <- opt$color_var
if (!color_var %in% c("p.adjust", "pvalue", "qvalue")) {
  warning("未知的 color_var: ", color_var, "。重置为 p.adjust")
  color_var <- "p.adjust"
}

p <- dotplot(gsea_obj, 
             showCategory = show_category,
             color = color_var,
             font.size = opt$y_text_size, # dotplot font.size 控制的是Y轴标签大小
             title = "") # 不设标题

# 自定义样式
p <- p + 
  scale_color_gradient(low = opt$low_color, high = opt$high_color) +
  xlab(opt$x_title) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = opt$y_label_wrap)) +
  theme(
    text = element_text(family = plot_font),
    axis.title.x = element_text(size = opt$x_title_size),
    axis.text.x = element_text(size = opt$x_text_size),
    axis.text.y = element_text(size = opt$y_text_size),
    legend.title = element_text(size = opt$legend_title_size),
    legend.text = element_text(size = opt$legend_text_size),
    plot.margin = unit(c(1, 1, 1, 1), 'lines')
  )

# 保存输出
plot_dir <- dirname(opt$plot_name)
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
out_pdf <- paste0(opt$plot_name, ".pdf")
out_png <- paste0(opt$plot_name, ".png")

ggsave(filename = out_pdf, plot = p, device = "pdf", width = opt$width, height = opt$height)
ggsave(filename = out_png, plot = p, device = "png", width = opt$width, height = opt$height)

cat("绘图完成。\nPDF: ", out_pdf, "\nPNG: ", out_png, "\n")

file.create(file.path(opt$outdir, "done"))
