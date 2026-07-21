.libPaths("/mnt/NexData_file/share/scRNA-rawdata/envs/scRNA_p/lib/R/library")
library(optparse)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggpubr)

# 定义命令行参数
option_list = list(
  make_option(c("-d", "--dat"), type = "character", help = "Path to the input data file (tab-separated)"),
  make_option(c("-o", "--output"), type = "character", help = "Path to the output directory")
)

# 解析命令行参数
opt = parse_args(OptionParser(option_list = option_list))

# 检查参数
if (is.null(opt$dat) | is.null(opt$output)) {
  stop("Both --dat and --output must be provided.")
}

# 读取数据
dat <- read.csv(opt$dat, sep = "\t", row.names = 1)

# 将除 cell_type 和 group 外的列转换为数字类型
dat[, !(names(dat) %in% c("cell_type", "group"))] <- 
  lapply(dat[, !(names(dat) %in% c("cell_type", "group"))], as.numeric)

# 获取所有的组别
group <- unique(dat$group)

# 获取所有的细胞类型
cell_type <- unique(dat$cell_type)

# 根据组别长度自动定义颜色
custom_colors <- brewer.pal(length(group), "Set3")
names(custom_colors) <- group

# 创建输出目录（如果不存在）
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# 获取细胞类型的数量
cell_type_count <- length(cell_type)

# 计算图形宽度和高度
# 假设每一列 facet_wrap 会占据固定的宽度与高度
facet_col <- 3  # 每行2列显示
facet_row <- ceiling(cell_type_count / facet_col)  # 根据细胞类型数量计算行数

# 如果只有一个 group，只画 cell_type 间比较的图
if (length(group) == 1) {
  
  # 不分组，cell_type 间比较
  for (term in setdiff(names(dat), c("cell_type", "group"))) {
    gg2 <- ggplot(dat, aes(x = cell_type, y = .data[[term]])) +  # 使用 .data[[term]] 获取动态列
      geom_violin(aes(fill = cell_type), trim = FALSE) +  # 使用 geom_violin 绘制小提琴图
      scale_fill_manual(values = brewer.pal(length(unique(dat$cell_type)), "Set3")) +
      theme_classic() +
      stat_compare_means(method = "wilcox.test", comparisons = combn(unique(dat$cell_type), 2, simplify = FALSE), label = "p.format") +
      labs(title = paste(term, "Scores - Cell Type Comparison"))
    
    # 保存不分组的 `cell_type` 小提琴图
    ggsave(file.path(opt$output, paste0(term, "_CellType_Comparison_ViolinPlot.pdf")), gg2, width = 10, height = 8+ 0.15* length(cell_type))
    ggsave(file.path(opt$output, paste0(term, "_CellType_Comparison_ViolinPlot.png")), gg2, width = 10, height = 8+ 0.15* length(cell_type))
  }
  
} else if (length(cell_type) == 1) {
  
  # 如果只有一个 cell_type，只画不分 cell_type, group 组间比较的图
  for (term in setdiff(names(dat), c("cell_type", "group"))) {
    gg1 <- ggplot(dat, aes(x = group, y = .data[[term]])) +  # 使用 .data[[term]] 获取动态列
      geom_violin(aes(fill = group), trim = FALSE) +  # 使用 geom_violin 绘制小提琴图
      scale_fill_manual(values = custom_colors) +
      theme_classic() +
      stat_compare_means(method = "wilcox.test", comparisons = combn(group, 2, simplify = FALSE), label = "p.format") +
      labs(title = paste(term, "Scores - Group Comparison"))
    
    # 保存不分 `cell_type` 的小提琴图
    ggsave(file.path(opt$output, paste0(term, "_Group_Comparison_ViolinPlot.pdf")), gg1, width = 10, height = 8+ 0.15* length(group))
    ggsave(file.path(opt$output, paste0(term, "_Group_Comparison_ViolinPlot.png")), gg1, width = 10, height = 8+ 0.15* length(group))
  }
  
} else {
  
  # 正常情况下，执行三种情况的图绘制
  for (term in setdiff(names(dat), c("cell_type", "group"))) {
    
    # 1. 不分 cell_type, group 组间比较
    gg1 <- ggplot(dat, aes(x = group, y = .data[[term]])) +  # 使用 .data[[term]] 获取动态列
      geom_violin(aes(fill = group), trim = FALSE) +  # 使用 geom_violin 绘制小提琴图
      scale_fill_manual(values = custom_colors) +
      theme_classic() +
      stat_compare_means(method = "wilcox.test", comparisons = combn(group, 2, simplify = FALSE), label = "p.format") +
      labs(title = paste(term, "Scores - Group Comparison"))
    
    # 保存不分 `cell_type` 的小提琴图
    ggsave(file.path(opt$output, paste0(term, "_Group_Comparison_ViolinPlot.pdf")), gg1, width = 10, height = 8+ 0.15* length(group))
    ggsave(file.path(opt$output, paste0(term, "_Group_Comparison_ViolinPlot.png")), gg1, width = 10, height = 8+ 0.15* length(group))

    
    # 2. 不分组，cell_type 间比较
    gg2 <- ggplot(dat, aes(x = cell_type, y = .data[[term]])) +  # 使用 .data[[term]] 获取动态列
      geom_violin(aes(fill = cell_type), trim = FALSE) +  # 使用 geom_violin 绘制小提琴图
      scale_fill_manual(values = brewer.pal(length(unique(dat$cell_type)), "Set3")) +
      theme_classic() +
      stat_compare_means(method = "wilcox.test", comparisons = combn(unique(dat$cell_type), 2, simplify = FALSE), label = "p.format") +
      labs(title = paste(term, "Scores - Cell Type Comparison"))
    
    # 保存不分组的 `cell_type` 小提琴图
    ggsave(file.path(opt$output, paste0(term, "_CellType_Comparison_ViolinPlot.pdf")), gg2, width = 10, height = 8+ 0.15* length(cell_type))
    ggsave(file.path(opt$output, paste0(term, "_CellType_Comparison_ViolinPlot.png")), gg2, width = 10, height = 8+ 0.15* length(cell_type))

    
    # 3. 对每个组，单独画 cell_type 间比较的小提琴图
    for (grp in group) {
      subset_data <- dat[dat$group == grp, ]
      
      gg3 <- ggplot(subset_data, aes(x = cell_type, y = .data[[term]])) +  # 使用 .data[[term]] 获取动态列
        geom_violin(aes(fill = cell_type), trim = FALSE) +  # 使用 geom_violin 绘制小提琴图
        scale_fill_manual(values = brewer.pal(length(unique(subset_data$cell_type)), "Set3")) +
        theme_classic() +
        stat_compare_means(method = "wilcox.test", comparisons = combn(unique(subset_data$cell_type), 2, simplify = FALSE), label = "p.format") +
        labs(title = paste(term, "Scores - Cell Type Comparison within", grp))
      
      # 保存每个组的 `cell_type` 小提琴图
      ggsave(file.path(opt$output, paste0(term, "_", grp, "_CellType_Comparison_ViolinPlot.pdf")), gg3, width = 10, height = 8+ 0.15* length(cell_type))
      ggsave(file.path(opt$output, paste0(term, "_", grp, "_CellType_Comparison_ViolinPlot.png")), gg3, width = 10, height = 8+ 0.15* length(cell_type))

    # 4. 每个cell_type,组间比较的小提琴图
    for (term in setdiff(names(dat), c("cell_type", "group"))) {
    gg1 <- ggplot(dat, aes(x = group, y = .data[[term]])) +  # 使用 .data[[term]] 获取动态列
    geom_violin(aes(fill = group), trim = FALSE) +  # 使用 geom_violin 绘制小提琴图
    scale_fill_manual(values = custom_colors) +
    theme_classic() +
    stat_compare_means(method = "wilcox.test", comparisons = combn(group, 2, simplify = FALSE), label = "p.format") +
    facet_wrap(~ cell_type, scales = 'free', ncol = facet_col) +  # 动态调整列数
    labs(title = paste(term, "Scores"))
  
  # 保存图形，文件名为列名 + '_Scores_ViolinPlot'
  ggsave(file.path(opt$output, paste0(term, "_Scores_ViolinPlot.pdf")), gg1, width = 15, height = 5 * facet_row)  # 高度根据行数调整
  ggsave(file.path(opt$output, paste0(term, "_Scores_ViolinPlot.png")), gg1, width = 15, height = 5 * facet_row)  # 高度根据行数调整
}
    }
  }
}

cat("Violin plots saved to ", opt$output, "\n")
