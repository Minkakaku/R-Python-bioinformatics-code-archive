library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(ggsci)
library(plotly)
library(readxl)
library(optparse)


# 1. 参数设置 ==================================================================
option_list <- list(
  make_option(c("-i", "--data_input"), type="character", default="input_data.txt", help="输入数据文件 (txt/xlsx)"),
  make_option(c('--enrich_type'), type="character", default="KEGG", help="绘图类型: 'Biological process', 'Cellular component', 'Molecular function', 'KEGG', 'Pfam' 等"),
  make_option(c("--p_cutoff"), type="numeric", default=0.05, help="P值过滤阈值"),
  make_option(c("--top_n"), type="numeric", default=10, help="每个分组显示前N个条目 (0表示显示全部)"),
  make_option(c("--label_size"), type="numeric", default=12, help="坐标轴/标签字体大小"),
  make_option(c("--palette"), type="character", default="Reds", help="配色方案: 
              连续渐变(P值): RdYlBu,Spectral,Blues,Reds 等 (RColorBrewer);
              离散颜色(分组): npg,jco,lancet (ggsci)"),
  make_option(c("--legend_pos"), type="character", default="right", help="图例位置 (top,bottom,right,left,none)"),
  make_option(c("--str_wrap"), type="numeric", default=45, help="Y轴文字换行长度"),
  make_option(c("--fig_width"), type="numeric", default=12, help="图片宽度"),
  make_option(c("--fig_height"), type="numeric", default=10, help="图片高度"),
  make_option(c("--plot_name"), type="character", help="生成图名称(包含路径和文件名，无后缀)"),
  make_option(c("--outdir"), type="character", default="./", help="输出目录")
)
opt <- parse_args(OptionParser(option_list = option_list))

## 检查并创建输出目录 --------------------------------------------------------
dir.create(file.path(opt$outdir, "5-Functional_enrichment"), showWarnings=FALSE)

# 2. 检查并创建输出目录 ========================================================
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# 3. 读入数据 ==================================================================
inFile <- opt$data_input
ext <- tolower(tools::file_ext(inFile))
if (ext == "xlsx") {
  enrich.data = as.data.frame(read_xlsx(inFile))    
} else {
  enrich.data = fread(inFile)
}

colnames(enrich.data)[colnames(enrich.data) == "Fold enrich"] <- "Fold_Enrich"
colnames(enrich.data)[colnames(enrich.data) == "p.value"] <- "P_Value"
colnames(enrich.data)[colnames(enrich.data) == "mapping"] <- "Count" 

# 检查必须的列是否存在
required_cols <- c("Category", "Regulation", "Terms", "Fold_Enrich", "P_Value", "Count")
missing_cols <- setdiff(required_cols, colnames(enrich.data))
if(length(missing_cols) > 0) {
  stop(paste("错误: 数据中缺少以下列 (请检查拼写/大小写):", paste(missing_cols, collapse=", ")))
}


# 4. 根据 enrich_type 筛选数据 =================================================

target_type_input <- toupper(opt$enrich_type) 
target_cats <- c()

if (target_type_input == "GO") {
  target_cats <- c("Biological process", "Cellular component", "Molecular function")
  cat("[INFO] 模式: GO (合并展示 BP, CC, MF)\n")
} else if (target_type_input == "Biological process") {
  target_cats <- c("Biological process")
  cat("[INFO] 模式: BP (仅展示 Biological process)\n")
} else if (target_type_input == "Cellular component") {
  target_cats <- c("Cellular component")
  cat("[INFO] 模式: CC (仅展示 Cellular component)\n")
} else if (target_type_input == "Molecular function") {
  target_cats <- c("Molecular function")
  cat("[INFO] 模式: MF (仅展示 Molecular function)\n")
} else {
  target_cats <- c(opt$enrich_type)
  cat(sprintf("[INFO] 模式: 自定义分类 (%s)\n", opt$enrich_type))
}

# 过滤数据
plot_data <- enrich.data %>%
  filter(Category %in% target_cats) %>%
  filter(P_Value < opt$p_cutoff)

if (nrow(plot_data) == 0) {
  stop(sprintf("错误: 在筛选条件 (Category=%s, P<%s) 下无数据。\n请检查 'Category' 列内容或放宽 P 值。",
               paste(target_cats, collapse="/"), opt$p_cutoff))
}

# 5. 数据处理 (Top N & 排序) ===================================================

# 处理 Top N
# 逻辑：如果 top_n > 0，则截取；如果为 0，则保留所有数据
if (opt$top_n > 0) {
  cat(sprintf("[INFO] 筛选前 %d 个显著条目\n", opt$top_n))
  plot_data <- plot_data %>%
    group_by(Category, Regulation) %>%
    arrange(P_Value) %>%
    slice_head(n = opt$top_n) %>%
    ungroup()
} else {
  cat("[INFO] Top_N 设置为 0，显示所有符合 P 值的条目\n")
}

# 文本换行处理
plot_data$Terms <- str_wrap(plot_data$Terms, width = opt$str_wrap)

# 因子排序 (使气泡在图中按 Fold_Enrich 大小排列)
plot_data <- plot_data %>%
  arrange(Category, Regulation, Fold_Enrich) %>%
  mutate(Terms = factor(Terms, levels = unique(Terms)))

# 6. 绘图配置与生成 ============================================================

cat("[INFO] 正在绘图...\n")

# 判断配色模式：是离散(ggsci) 还是 连续(Brewer)
# 常用 ggsci 包列表
ggsci_palettes <- c("npg", "aaas", "nejm", "lancet", "jama", "jco", "ucscgb", "d3", "locuszoom", "igv", "ucheck", "startrek", "tron", "futurama", "rickandmorty")

is_discrete_palette <- opt$palette %in% ggsci_palettes

# 基础绘图对象
p <- ggplot(plot_data, aes(y = Terms, x = log2(Fold_Enrich)))

# 根据配色方案决定颜色映射逻辑
if (is_discrete_palette) {
  # === 模式 A: 离散配色 (如 npg, jco) ===
  # 颜色映射到 Regulation (Up/Down)，P值仅用于气泡大小或其他(此处P值暂不映射颜色，或者可以映射透明度)
  # 这里为了美观，让颜色代表分组(Up/Down)，大小代表Count
  cat(sprintf("[INFO] 使用离散配色方案: %s (颜色映射到 Regulation)\n", opt$palette))
  
  p <- p + geom_point(aes(size = Count, color = Regulation))
  
  # 动态调用 ggsci 的配色函数
  scale_func_name <- paste0("scale_color_", opt$palette)
  if (exists(scale_func_name)) {
    p <- p + match.fun(scale_func_name)()
  } else {
    warning("未找到指定的 ggsci 配色函数，使用默认 npg")
    p <- p + scale_color_npg()
  }
  
} else {
  # === 模式 B: 连续配色 (默认, 如 RdYlBu) ===
  # 颜色映射到 P_Value (显著性)
  cat(sprintf("[INFO] 使用连续配色方案: %s (颜色映射到 P_Value)\n", opt$palette))
  
  p <- p + geom_point(aes(size = Count, color = P_Value)) +
    # direction = -1 通常让 P值越小(越显著)颜色越深/红
    scale_color_distiller(palette = opt$palette, direction = -1)
}

# 继续添加通用图层
p <- p +
  # 分面设置
  facet_grid(~ Regulation, scales = "free", space = "free_y") +
  
  # 气泡大小范围
  scale_size(range = c(2, 6)) +
  
  # 主题设置
  theme_test(base_size = opt$label_size) + # 使用 label_size 控制基准字体
  theme(
    axis.text.y = element_text(color = "black", size = opt$label_size),
    axis.text.x = element_text(color = "black", size = opt$label_size),
    strip.text = element_text(face = "bold", size = opt$label_size + 1),
    legend.position = opt$legend_pos, # 使用参数控制图例位置
    legend.title = element_text(face = "bold")
  ) +
  
  # 标签
  labs(
    x = "log2(Fold Enrichment)",
    y = NULL,
    color = ifelse(is_discrete_palette, "Regulation", "p.value"),
    size = "mapping"
  )

ggsave(plot=p,filename=paste0(opt$plot_name, '.png'), width = opt$fig_width,height = opt$fig_height)
ggsave(plot=p,filename=paste0(opt$plot_name, '.pdf'), width = opt$fig_width,height = opt$fig_height)

## 统一结束标识符
done_file <- file.path(opt$outdir, "done")
file.create(done_file)
