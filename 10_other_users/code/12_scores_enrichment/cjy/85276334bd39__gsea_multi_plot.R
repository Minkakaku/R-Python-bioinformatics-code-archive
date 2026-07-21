#!/usr/bin/env Rscript

# GSEA多通路绘图脚本 (支持Description列表)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(optparse)
})

# 命令行参数定义
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "输入GSEA结果RData文件路径 [必需]"),
  make_option(c("-t", "--terms"), type = "character", default = NULL,
              help = "多个通路Description，使用英文分号分隔，例如: 'cell cycle;DNA replication'。如果不传或传'-'，则默认展示Top 5 [可选]"),
  make_option(c("--plot_name"), type = "character", default = "gsea_multi_plot",
              help = "生成图名称(包含路径，无后缀) [默认: %default]"),
  make_option(c("--outdir"), type = "character", default = "./",
              help = "输出目录 [默认: %default]"),
  make_option(c("--width"), type = "numeric", default = 16,
              help = "图形宽度（英寸）[默认: %default]"),
  make_option(c("--height"), type = "numeric", default = 15,
              help = "图形高度（英寸）[默认: %default]"),
  make_option(c("--font"), type = "character", default = "Arial",
              help = "字体族名称 (Times New Roman, Arial) [默认: %default]"),
  make_option(c("--x_axis_title_size"), type = "numeric", default = 12,
              help = "X轴标题字体大小 [默认: %default]"),
  make_option(c("--y_axis_title_size"), type = "numeric", default = 12,
              help = "Y轴标题字体大小 [默认: %default]"),
  make_option(c("--x_axis_text_size"), type = "numeric", default = 10,
              help = "X轴刻度字体大小 [默认: %default]"),
  make_option(c("--y_axis_text_size"), type = "numeric", default = 10,
              help = "Y轴刻度字体大小 [默认: %default]"),
  make_option(c("--show_pvalue"), type = "character", default = "yes",
              help = "是否显示P值表格 (yes=显示, no=隐藏) [默认: %default]"),
  make_option(c("--resolution"), type = "numeric", default = 300,
              help = "图片分辨率（DPI）[默认: %default]")
)

# 解析命令行
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 将输入字符串解析为Term向量（支持空格/分号）
parse_terms <- function(x) {
  terms <- unlist(strsplit(x, ";"))
  terms <- trimws(terms)
  terms <- terms[terms != ""]
  unique(terms)
}

# 参数验证
validate_params <- function(opt) {
  if (is.null(opt$input)) {
    stop("错误: 必须提供输入文件路径 (-i 或 --input)")
  }
  if (!file.exists(opt$input)) {
    stop("错误: 输入文件 '", opt$input, "' 不存在")
  }

  # 解析 Terms
  # 如果 opt$terms 为 NULL 或 "-"，返回 NULL 表示使用默认 Top 5
  terms <- NULL
  if (!is.null(opt$terms) && opt$terms != "-") {
    terms <- parse_terms(opt$terms)
    if (length(terms) == 0) {
      stop("错误: 提供的 terms 参数无法解析出有效内容")
    }
  }

  # 字体与数值参数校验
  valid_fonts <- c("Times New Roman", "Arial")
  if (!(opt$font %in% valid_fonts)) {
    stop("错误: 字体必须是以下之一: ", paste(valid_fonts, collapse = ", "))
  }
  numeric_params <- c("width", "height", "x_axis_title_size", "y_axis_title_size",
                      "x_axis_text_size", "y_axis_text_size", "resolution")
  for (param in numeric_params) {
    if (opt[[param]] <= 0) {
      stop("错误: ", param, " 必须为正数")
    }
  }
  valid_pvalue <- c("yes", "no")
  if (!(opt$show_pvalue %in% valid_pvalue)) {
    stop("错误: show_pvalue必须是以下之一: ", paste(valid_pvalue, collapse = ", "))
  }

  # 返回解析后的 Terms 向量（或 NULL）
  terms
}

# 加载GSEA结果对象
load_gsea_data <- function(file_path) {
  cat("加载GSEA结果数据...", "\n")
  
  # 创建临时环境以加载数据，避免污染当前环境
  env <- new.env()
  load(file_path, envir = env)
  
  # 定义支持的对象名称列表
  candidate_names <- c("GO_gseresult", "KEGG_gseresult", "Reactome_gseresult")
  
  # 查找存在的对象
  found_name <- NULL
  for (name in candidate_names) {
    if (exists(name, envir = env)) {
      found_name <- name
      break
    }
  }
  
  # 如果未找到标准命名对象，尝试搜索任何 gseaResult 类型的对象
  if (is.null(found_name)) {
    all_objs <- ls(envir = env)
    for (obj_name in all_objs) {
      obj <- get(obj_name, envir = env)
      if (inherits(obj, "gseaResult")) {
        found_name <- obj_name
        cat("提示: 未找到标准命名对象(GO/KEGG/Reactome_gseresult)，但发现 gseaResult 对象:", found_name, "\n")
        break
      }
    }
  }
  
  if (is.null(found_name)) {
    stop("错误: RData文件中未找到支持的 GSEA 结果对象 (GO_gseresult, KEGG_gseresult, Reactome_gseresult)")
  }
  
  gsea_result <- get(found_name, envir = env)
  
  if (!inherits(gsea_result, "gseaResult")) {
    stop("错误: 对象 '", found_name, "' 不是有效的 gseaResult 对象")
  }
  
  cat("成功加载 GSEA 结果对象:", found_name, "，包含", nrow(gsea_result@result), "个基因集", "\n")
  gsea_result
}

# 查找 Terms 在结果中的位置
find_term_positions <- function(gsea_result, terms) {
  # 如果 terms 为 NULL，则默认取 Top 5
  if (is.null(terms)) {
    cat("未指定具体 Terms，默认选择 Top 5 通路...\n")
    n_total <- nrow(gsea_result@result)
    n_select <- min(5, n_total)
    positions <- seq_len(n_select)
    
    cat("选定的 Top", n_select, "通路:\n")
    for (i in positions) {
      cat(sprintf("[%d] %s (ID: %s)\n", i, gsea_result@result$Description[i], gsea_result@result$ID[i]))
    }
    return(positions)
  }

  cat("查找 Term 列表:", paste(terms, collapse = " ; "), "\n")
  all_descriptions <- gsea_result@result$Description
  
  # 使用 exact match
  positions <- match(terms, all_descriptions)

  missing_idx <- which(is.na(positions))
  if (length(missing_idx) > 0) {
    missing <- terms[missing_idx]
    # 尝试模糊匹配提示（可选，这里先报错）
    stop("错误: 在GSEA结果中未找到以下 Description: ", paste(missing, collapse = " ; "))
  }

  # 输出每个Term的位置，便于核对
  for (i in seq_along(terms)) {
    cat("找到 Term '", terms[i], "' 位置:", positions[i], "; ID:",
        gsea_result@result$ID[positions[i]], "\n")
  }
  positions
}

# 生成叠加GSEA图
generate_gsea_plot <- function(gsea_result, positions, opt) {
  cat("生成GSEA叠加图...", "\n")
  
  # 注意：gseaplot2 支持传入 numeric index 作为 geneSetID
  p <- gseaplot2(
    gsea_result,
    geneSetID = positions,
    color = "red", # 如果用户没指定特定颜色列表，这里使用参数中的颜色（gseaplot2会自动处理单色或多色）
    pvalue_table = (opt$show_pvalue == "yes")
  )
  
  # 如果 curve_color 是单个颜色且有多条曲线，gseaplot2 可能只会用这个颜色画所有线
  # 如果需要区分颜色，通常 gseaplot2 会默认分配。
  # 这里的 opt$curve_color 默认是 "red"。如果用户想保持原样则不动。
  # 如果 positions 数量 > 1 且 color 只给了一个值，所有线都会是那个颜色。
  # 但通常多条线我们希望不同颜色。
  # 原始脚本逻辑是传入 "red"，这意味着所有线都是红色？
  # 查看原脚本： make_option(c("--curve_color"), ..., default = "red")
  # gseaplot2(..., color = "red")
  # 确实如此。如果不改这部分逻辑，就保持原样。

  # 关键：创建可调字体与大小的主题
  custom_theme <- theme(
    text = element_text(family = opt$font),
    axis.title.x = element_text(size = opt$x_axis_title_size, family = opt$font),
    axis.title.y = element_text(size = opt$y_axis_title_size, family = opt$font),
    axis.text.x = element_text(size = opt$x_axis_text_size, family = opt$font),
    axis.text.y = element_text(size = opt$y_axis_text_size, family = opt$font)
  )

  # 关键：根据返回对象类型应用主题
  if (inherits(p, "ggplot")) {
    cat("检测到单个ggplot对象，应用自定义主题", "\n")
    p <- p + custom_theme
  } else if (inherits(p, "gglist") || is.list(p)) {
    cat("检测到gglist复合图形对象，对子图应用自定义主题", "\n")
    for (i in c(1, 3)) {
      if (inherits(p[[i]], "ggplot")) {
        if (i == 1) {
          # 关键：子图1为ES曲线，移除X轴刻度文本以贴合示例风格
          subplot_theme <- theme(
            text = element_text(family = opt$font),
            axis.title.x = element_text(size = opt$x_axis_title_size, family = opt$font),
            axis.title.y = element_text(size = opt$y_axis_title_size, family = opt$font),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = opt$y_axis_text_size, family = opt$font)
          )
        } else {
          subplot_theme <- custom_theme
        }
        p[[i]] <- p[[i]] + subplot_theme
      }
    }
  } else {
    cat("检测到其他类型的图形对象:", class(p), "，无法应用自定义主题", "\n")
  }

  p
}

# 保存为多种格式
save_plots <- function(plot, plot_name, opt) {
  cat("保存图片...", "\n")
  
  # 创建图片目录
  plot_dir <- dirname(plot_name)
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  pdf_file <- paste0(plot_name, ".pdf")
  pdf(file = pdf_file, width = opt$width, height = opt$height)
  print(plot)
  dev.off()
  cat("已保存PDF文件:", pdf_file, "\n")

  png_file <- paste0(plot_name, ".png")
  png(filename = png_file, width = opt$width, height = opt$height,
      units = "in", res = opt$resolution)
  print(plot)
  dev.off()
  cat("已保存PNG文件:", png_file, "\n")
}

# 主流程
main <- function() {
  cat("=== GSEA多通路绘图工具 ===", "\n")

  terms <- validate_params(opt)  # 关键：返回已验证的 Terms 向量 或 NULL
  gsea_result <- load_gsea_data(opt$input)
  positions <- find_term_positions(gsea_result, terms)  # 获得位置列表 (Indices)
  gsea_plot <- generate_gsea_plot(gsea_result, positions, opt)
  save_plots(gsea_plot, opt$plot_name, opt)
  
  if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE)
  }
  file.create(file.path(opt$outdir, "done"))

  cat("=== 绘图完成 ===", "\n")
}

if (!interactive()) {
  main()
}
