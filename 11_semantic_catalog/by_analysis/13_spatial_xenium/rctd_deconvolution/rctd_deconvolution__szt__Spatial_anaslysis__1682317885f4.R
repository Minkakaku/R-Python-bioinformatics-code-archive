# spatial analysis of HD ------

library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggplot2)
source("spatial analysis.R")

## read data
all_data <- readRDS("data_for_analysis/all_data_projected_rctd_20250902.rds")
head(all_data@meta.data)


# rctd remove other --------
# Date <- Sys.Date()
Date <- "2025-09-09"
spatial_out_dir <- str_glue("spatial_analysis-{Date}")
dir.create(spatial_out_dir)

all_data@meta.data$annotation <- all_data@meta.data$rctd_20250902
all_data_temp <- all_data
all_data <- subset(all_data_temp, annotation != 'Other')
Idents(all_data) <- "annotation"

sample_mapping <- c(
  "X24M.16.008um" = "24M-1",
  "X24M.16.008um.2" = "24M-2",
  "X4M.7.008um" = "4M-1",
  "X4M.7.008um.4" = "4M-2"
)

cells_position_ls <- 
  all_data@images %>% 
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) %>%
  # setNames(., str_remove(names(.), "\\..*$"))
  setNames(., sample_mapping[names(.)])

cells_position_df <- cells_position_ls %>% bind_rows() 



# CN analysis ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

CN_out_dir <- file.path(spatial_out_dir, "CN analysis")
dir.create(CN_out_dir, recursive = T)


all_data_meta_for_cn_analysis <- 
  all_data[[c("project", "group", "annotation")]] %>%
  dplyr::mutate(annotation_full = annotation ) %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x', 'y')], by = "cell_id") %>%
  rename(X = x, Y = y) %>%
  mutate(
    annotation_full = factor(annotation_full), # import!
    obj_id = cell_id # import!
  )


### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn_10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn_10_heat$gtable)
dev.off()


#### CN annotaton --------------------
# cn_annotation <- c(
#   "CN1"  = "Immune core",
#   "CN2"  = "Plasma core",
#   "CN3"  = "Tregs immune",
#   "CN4"  = "Mac_Fibro_Ec",
#   "CN5"  = "Myo Stroma",
#   "CN6"  = "Fibro Stroma",
#   "CN7"  = "CD163 Macophage Stroma",
#   "CN8"  = "Enterocytes",
#   "CN9"  = "Myeloid enterocytes",
#   "CN10" = "Amplifying area"
# )

cn_annotation <- c(
  "CN1"  = "CN1",
  "CN2"  = "CN2",
  "CN3"  = "CN3",
  "CN4"  = "CN4",
  "CN5"  = "CN5",
  "CN6"  = "CN6",
  "CN7"  = "CN7",
  "CN8"  = "CN8",
  "CN9"  = "CN9",
  "CN10" = "CN10"
)

p_cn_10_heat2 <-
  cn_cluster_FC %>%
  pheatmap::pheatmap(
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = cn_annotation,
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap_annotation.pdf"),
  height = 6,
  width  = 11
)
cowplot::plot_grid(p_cn_10_heat2$gtable)
dev.off()




### CN k15 -----
res_cn_15 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 15
)

saveRDS(res_cn_15, file.path(CN_out_dir, "res_cn_15.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn15_cluster_mean <- res_cn_15$cluster_centers
cn15_baseMean     <- res_cn_15$base_mean

# log2(FC)
cn15_cluster_FC <- log2(1 + sweep(cn15_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn15_baseMean))

# heatmap
p_cn15_heat <-
  pheatmap::pheatmap(
    cn15_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn15_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_15_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn15_heat$gtable)
dev.off()




### spatial plot -----------------
sp_out_dir <- file.path(CN_out_dir, 'spatial plot') 
dir.create(sp_out_dir)


p_sp_plot_anno <-
  # Images(all_data) %>%
  # setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  # imap(function(img, img_name){
  sample_mapping %>% 
  imap(function(img_name, img){
    # print(img_name)
    # print(img)
    # return()
    p <- SpatialDimPlot(
      all_data,
      images = img, 
      group.by = "annotation",
      pt.size.factor	= 24
    ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

# export
pdf(
  file.path(sp_out_dir, "spatial plot annotation.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_anno %>% walk(print)
dev.off()

### CN spatial plot -----------------
cn_sp_out_dir <- file.path(CN_out_dir, 'CN spatial plot') 
dir.create(cn_sp_out_dir)

# cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

# cn10 
all_data_meta_cn10 <- 
  all_data[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn10, file.path(CN_out_dir, 'all_data_meta_cn10.rds') )
# all_data_meta_cn10 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn10.rds'))

cn10_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


# cn15 
all_data_meta_cn15 <- 
  all_data[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn15, file.path(CN_out_dir, 'all_data_meta_cn15.rds') )
# all_data_meta_cn15 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn15.rds'))

cn15_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn15, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn15_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN15 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

# ### CN2 spatial plot -----------------
# cn_sp_out_dir <- file.path(CN_out_dir, 'cn spatial plot') 
# dir.create(cn_sp_out_dir)
# 
# # cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# # image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')
# 
# all_data_meta <- 
#   all_data[[]] %>%
#   rownames_to_column('cell_id') %>%
#   left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
#   mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
# 
# saveRDS(all_data_meta, file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# 
# cn_sp_p_ls <-
#   cells_position_ls %>% 
#   imap(function(cells_position, project){ # browser()
#     
#     # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
#     # group_info  <- anno_data$group[anno_data$project == project]
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = y, y = x)) +
#       geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#       scale_y_reverse() +
#       coord_equal() +
#       ggsci::scale_fill_igv() +
#       theme_void() +
#       labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
#     p
#   })
# 
# # export image
# cn_sp_p_ls %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir, str_glue("CN spatial plot {name}.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })


### CN spatial plot2 -----------------
# cn_sp_out_dir2 <- file.path(CN_out_dir, 'cn spatial plot2') 
# dir.create(cn_sp_out_dir2)
# 
# all_data_meta <- readRDS( file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# cn_sp_p_ls2 <-
#   cells_position_ls %>% # .[1] %>%
#   imap(function(cells_position, project){ # browser()
#     cn_color <- structure(
#       names = paste0("CN", 1:10),
#       ggsci::pal_igv()(10)
#     ) %>%
#       c(
#         Others = "gray90"
#       )
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p_ls <-
#       paste0("CN", 1:10) %>% # .[1] %>%
#       set_names(., .) %>%
#       map(function(cn){
#         dat2 <- dat
#         
#         dat2$CN2 <- as.character(dat2$CN2)
#         dat2$CN2[dat2$CN2 != cn] <- "Others"
#         
#         p <- dat2 %>%
#           ggplot(aes(x = y, y = x)) +
#           geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#           scale_y_reverse() +
#           coord_equal() +
#           # ggsci::scale_fill_igv() +
#           scale_fill_manual(values = cn_color) +
#           theme_void() +
#           labs(fill = NULL) 
#         p
#       })
#     
#     p_ls
#   })

# # export image
# cn_sp_p_ls2 %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir2, str_glue("CN spatial plot {name} - split by CN.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })

### cn10 cell props ---------
all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn10_props <- 
  all_cn10_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

# all_cn_count_all_zeros_df <-
#   expand.grid(
#     CN2     = all_cn_count_df$CN2 %>% unique(),
#     project = all_cn_count_df$project %>% unique(),
#     group   = all_cn_count_df$group %>% unique()
#   )

p_cn10_props_group <- 
  all_cn10_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+ 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn10_props %>% print
p_cn10_props_group %>% print
dev.off()


### cn15 cell props ---------
all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn15_props <- 
  all_cn15_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

p_cn15_props_group <- 
  all_cn15_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+ 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN15_boxplot_and_barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn15_props %>% print
p_cn15_props_group %>% print
dev.off()










## 修改图片CN10 --------
all_data_meta_cn10 <- readRDS('spatial_analysis-2025-09-09/CN analysis/all_data_meta_cn10.rds')

all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))


all_cn10_count_df %>% head()


# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn10_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis-2025-09-09/CN analysis/shapiro_test_results_CN10.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn10_count_df, CN == "CN1")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.12, 0.16),
    breaks = seq(0.12, 0.16, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn10_count_df, CN == "CN2")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn10_count_df, CN == "CN3")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.12, 0.2),
    breaks = seq(0.12, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn10_count_df, CN == "CN4")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.02, 0.12),
    breaks = seq(0.02, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN4.pdf", p_detailed, width = 10, height = 8)





plot_data_CN5  <- subset(all_cn10_count_df, CN == "CN5")
print(plot_data_CN5)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.08, 0.16),
    breaks = seq(0.08, 0.16, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN5",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN5.pdf", p_detailed, width = 10, height = 8)





plot_data_CN6  <- subset(all_cn10_count_df, CN == "CN6")
print(plot_data_CN6)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.06, 0.1),
    breaks = seq(0.06, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN6",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn10_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.06, 0.12),
    breaks = seq(0.06, 0.12, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn10_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.025, 0.085),
    breaks = seq(0.025, 0.085, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn10_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.12, 0.18),
    breaks = seq(0.12, 0.18, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn10_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN10-CN10.pdf", p_detailed, width = 10, height = 8)











## 修改图片CN15 --------
all_data_meta_cn15 <- readRDS('spatial_analysis-2025-09-09/CN analysis/all_data_meta_cn15.rds')

all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn15_count_df %>% head()


# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn15_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis-2025-09-09/CN analysis/shapiro_test_results_CN15.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn15_count_df, CN == "CN1")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.07, 0.09),
    breaks = seq(0.07, 0.09, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn15_count_df, CN == "CN2")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.12),
    breaks = seq(0, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn15_count_df, CN == "CN3")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.02, 0.06),
    breaks = seq(0.02, 0.06, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn15_count_df, CN == "CN4")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.02, 0.12),
    breaks = seq(0.02, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN4.pdf", p_detailed, width = 10, height = 8)





plot_data_CN5  <- subset(all_cn15_count_df, CN == "CN5")
print(plot_data_CN5)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.11),
    breaks = seq(0.05, 0.11, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN5",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN5.pdf", p_detailed, width = 10, height = 8)





plot_data_CN6  <- subset(all_cn15_count_df, CN == "CN6")
print(plot_data_CN6)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.1),
    breaks = seq(0.04, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN6",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn15_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.1),
    breaks = seq(0.04, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn15_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.075, 0.12),
    breaks = seq(0.075, 0.12, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn15_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.1),
    breaks = seq(0.05, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn15_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.12),
    breaks = seq(0.04, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN10.pdf", p_detailed, width = 10, height = 8)






plot_data_CN11  <- subset(all_cn15_count_df, CN == "CN11")
print(plot_data_CN11)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN11, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.08),
    breaks = seq(0.05, 0.08, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN11",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN11.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN11.pdf", p_detailed, width = 10, height = 8)






plot_data_CN12  <- subset(all_cn15_count_df, CN == "CN12")
print(plot_data_CN12)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN12, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.1),
    breaks = seq(0.05, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN12",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN12.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN12.pdf", p_detailed, width = 10, height = 8)






plot_data_CN13  <- subset(all_cn15_count_df, CN == "CN13")
print(plot_data_CN13)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN13, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.07),
    breaks = seq(0.04, 0.07, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN13",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN13.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN13.pdf", p_detailed, width = 10, height = 8)






plot_data_CN14  <- subset(all_cn15_count_df, CN == "CN14")
print(plot_data_CN14)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN14, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.05),
    breaks = seq(0, 0.05, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN14",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN14.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN14.pdf", p_detailed, width = 10, height = 8)






plot_data_CN15  <- subset(all_cn15_count_df, CN == "CN15")
print(plot_data_CN15)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN15, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.09, 0.14),
    breaks = seq(0.09, 0.14, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN15",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN15.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis-2025-09-09/CN analysis/CN15-CN15.pdf", p_detailed, width = 10, height = 8)











# Proximity analysis -------
### 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

All_type_Proximity_out_dir <- file.path(spatial_out_dir, "Proximity analysis")
dir.create(All_type_Proximity_out_dir)

cells_position_ls <- 
  all_data@images %>%
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) %>% 
  # setNames(., str_remove(names(.), "\\..*$"))
  setNames(., sample_mapping[names(.)])



# all_data_cn_meta
all_data_meta_pos <- 
  cells_position_ls %>% 
  imap(function(cells_position, slide){# browser()
    all_data_meta <- 
      all_data[[]] %>% 
      dplyr::mutate(annotation_full = annotation)
    proj_meta_sub <- 
      all_data_meta %>% 
      dplyr::select(annotation_full, group, project) %>% 
      rownames_to_column("cell_id") %>% 
      dplyr::filter(project == slide)
    cells_position[, 1:3] %>% 
      inner_join(proj_meta_sub, by = "cell_id") %>% 
      mutate(annotation_full = factor(annotation_full)) %>% 
      dplyr::rename(Classification1 = annotation_full, 
                    X = x, 
                    Y = y, 
                    Slide_ID = project)
  }) %>% 
  enframe() %>% 
  unnest(cols = value) %>% 
  dplyr::select(-name) %>% 
  column_to_rownames("cell_id")

all_data_meta_pos$Slide_ID <- str_c(all_data_meta_pos$group, all_data_meta_pos$Slide_ID, sep = "_")

all_data_meta_pos_subset <- all_data_meta_pos %>% rownames_to_column("cell_id")
# readRDS('data/all_data_meta_pos_subset.rds')

all_data_meta_pos_subset %>% head

# cell_id        X        Y Classification1 group       Slide_ID
# 1 1_aaaaaagg-1 1670.389 185.1358      Macrophage  CD-F X_202413076_32
# 2 1_aaadlghn-1 1681.172 176.4662          Plasma  CD-F X_202413076_32
# 3 1_aaadmplc-1 1677.561 182.6699               T  CD-F X_202413076_32
# 4 1_aaadoogj-1 1677.109 165.6637          Plasma  CD-F X_202413076_32
# 5 1_aaadpljf-1 1671.887 176.8148          Plasma  CD-F X_202413076_32
# 6 1_aaafaiko-1 1719.466 192.8071          Plasma  CD-F X_202413076_32
# ......


### 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

#### 2.1 proximity analysis --------
res_pa <- 
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset,
    query = all_data_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 500, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 8 / all_data@images$X24M.16.008um@scale.factors$spot, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
saveRDS(res_pa, file = file.path(All_type_Proximity_out_dir, "res_pa.rds"))
res_pa <- readRDS(file.path(All_type_Proximity_out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')

# summarized by group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos_subset,
#     str_c(Classification1, group, Slide_ID, sep = '..')
#   )
pa_groups <- 
  all_data_meta_pos_subset %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE) 

pa_groups2 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn10[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

pa_groups3 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn15[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
res_pa_df2 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups2, FUN = "median")})
res_pa_df3 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups3, FUN = "median")})


res_pa_df  <- 
  res_pa_df %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project"), sep = "\\.\\.") 

res_pa_df2  <- 
  res_pa_df2 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.")
res_pa_df2$CN <- factor(res_pa_df2$CN, levels = paste0("CN", 1:10))

res_pa_df3  <- 
  res_pa_df3 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df3$CN <- factor(res_pa_df3$CN, levels = paste0("CN", 1:15))


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_res_pa_ls_cn10 <-
  res_pa_df2 %>%
  with(split(res_pa_df2, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

p_res_pa_ls_cn15 <-
  res_pa_df3 %>%
  with(split(res_pa_df3, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

# p_res_pa_ls <-
#   res_pa_df2 %>%
#   with(split(res_pa_df2, CN2)) %>%
#   imap(function(dat, name){
#     dat %>%
#   # dplyr::filter(Target == "Epithelial cells") %>%
#   ggplot(aes(x = Others, y = cell_props, fill = Group)) +
#   geom_boxplot(width = 0.45, outlier.colour = NA) +
#   # geom_jitter(width = 0.25) +
#   facet_wrap(~ Target, scales = "free_y", ncol = 2)+
#   ggpubr::theme_pubclean() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   ggtitle(name)
#   })

# export plot
pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis.pdf"),
  height = 50,
  width  = 25
)
p_res_pa%>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN10.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn10 %>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN15.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn15 %>% print()
dev.off()











# 修改图片
res_pa_df

shapiro_results <- res_pa_df %>%
  group_by(Target) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(cell_props)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(cell_props)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis-2025-09-09/Proximity analysis/shapiro_test_results.csv", row.names = FALSE)







#### 2.2 proximity analysis 2 -------
# ##total data  clean
# all_data_meta_pa <-
#   all_data[[]] %>%
#   mutate(annotation4 = ifelse(
#     as.character(annotation_full) == "Fibroblast",
#     as.character(annotation3),
#     as.character(annotation_full)
#   )) %>%
#   dplyr::select(-annotation_full) %>%
#   dplyr::rename(annotation_full = annotation4) %>%
#   mutate(annotation_final = annotation_full %>%
#            as.character() %>%
#            factor(., levels = unique(.))
#   ) %>%
#   mutate(
#     Group = group
#   )
# 
# 
# 
# all_data_meta_pos <-
#   cells_position_ls  %>%
#   imap(function(cells_position, slide){# browser()
#     
#     if('cell_id' %in% colnames(all_data_meta_pa) ){
#       all_data_meta_pa <- all_data_meta_pa %>%
#         remove_rownames() %>%
#         column_to_rownames("cell_id")
#     }
#     
#     proj_meta_sub <-
#       all_data_meta_pa %>%
#       dplyr::select(annotation_final, Group, project) %>%
#       rownames_to_column("cell_id") %>%
#       dplyr::filter(project == slide)
#     
#     if(! 'cell_id' %in% colnames(proj_meta_sub) ){
#       proj_meta_sub <- rownames_to_column(proj_meta_sub, "cell_id")
#     }
#     
#     cells_position[, 1:3] %>%
#       inner_join(proj_meta_sub, by = "cell_id") %>%
#       dplyr::rename(Classification1 = annotation_final,
#                     X = x_centroid,
#                     Y = y_centroid,
#                     Slide_ID = project)
#   }) %>%
#   enframe() %>%
#   unnest(cols = value) %>%
#   dplyr::select(-name) %>%
#   column_to_rownames("cell_id")

mpp <- 8 / all_data@images$TSA_233841.008um@scale.factors$spot

res_pa_ls <- list()

res_pa_ls[['10']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 10, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
res_pa_ls[['50']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 50, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['100']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 100, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['150']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 150, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

res_pa_ls[['200']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", #NULL: 提前返回window的细胞近邻的list
    k = 1200, # 'auto',
    searchtype = "radius",
    radius = 200, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

# saveRDS(res_pa_ls, file.path(CN_out_dir, "res_pa_ls.rds"))


# epithelial proximity
# proximity_cells <- c(
#   "THBS1 Fibroblast",
#   "TNXB Fibroblast",
#   "COL5A1 Fibroblast",
#   "EEF1G Fibroblast",
#   "ADAMDEC1 Fibroblast"
# )
proximity_cells <- c(
  "CEACAM6 Goblet"
)

proximity_all_data_ls <- 
  proximity_cells %>%
  setNames(., .) %>%
  map(function(cell){
    cell_idx <- all_data_meta_pos$Classification1 == cell
    
    stopifnot(sum(cell_idx) > 0)
    
    # 1 all_data
    all_data_sub  <- all_data_meta_pos[cell_idx, ]
    
    # 2 res_pa
    res_pa_sub_ls <- res_pa_ls %>% map(~ .x[, cell_idx])
    
    # distance diff
    diff_idx_ls <- seq_along(res_pa_sub_ls) %>% rev() %>% { matrix(c(.[-length(.)], .[-1]), ncol=2)} %>% apply(1, c, simplify = FALSE)
    
    
    # diff_idx_ls
    res_pa_sub_diff_ls <- vector(length = 5, mode = "list")
    
    for(idx in diff_idx_ls){
      res_pa_sub_diff_ls[[ idx[1] ]] <- res_pa_sub_ls[[ idx[1] ]] - res_pa_sub_ls[[ idx[2] ]]
    }
    res_pa_sub_diff_ls[[1]] <- res_pa_sub_ls[[1]]
    names(res_pa_sub_diff_ls) <- names(res_pa_sub_ls)
    
    res_pa_sub_count <-
      res_pa_sub_diff_ls %>% 
      map(~ apply(.x, 1, function(x) tapply(x, all_data_sub$Slide_ID, median)))  
    
    res_pa_sub_count <-
      res_pa_sub_count %>% 
      imap(function(x, y){
        x <- x %>% as.data.frame(check.names = FALSE) %>% rownames_to_column('Slide_ID')
        x$distance <- y
        x
      }) %>% 
      bind_rows() 
    
    
    # add Group 
    all_data_group_sub <- all_data_sub[c('group', 'Slide_ID')] %>% distinct()
    
    # return
    list(
      all_data = all_data_group_sub,
      res_pa   = res_pa_sub_count
    )
  })


# proximity_all_data_ls$Fibroblast$res_pa
# # res_pa_epi_count %>%
#   pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
#   left_join(all_data_group_epi, by = "Slide_ID") 


p_proximity_cell_count_distribution_ls <-
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      mutate(distance = as.numeric(distance)) %>%
      
      ggplot(aes(x = Cell, y = count, fill = group)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ distance, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # scale_fill_manual(values = grp_cols_sp) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p 
  })



res_pa_count_tidy_ls <- 
  proximity_all_data_ls %>%
  map(function(dat_ls){ # browser()
    
    all_data_group_sub <- dat_ls$all_data %>% map_dfc(as.character)
    
    res_pa_count_raw <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID") %>%
      mutate(distance = as.numeric(distance)) 
    
    
    res_pa_count_tidy <-
      1:nrow(res_pa_count_raw) %>%  
      map(function(x){
        count <- res_pa_count_raw$count[x] 
        if(count < 1) return(NULL) # <- count是median count，有可能是0.5
        xx <- res_pa_count_raw[x, ] %>% unlist(use.names = F) %>% rep(count) %>% matrix(nrow = count, byrow = T) %>% as.data.frame()
        xx
      }) %>%
      purrr::discard(is.null) %>%
      bind_rows()
    colnames(res_pa_count_tidy) <- colnames(res_pa_count_raw)
    
    res_pa_count_tidy$count <- as.numeric(res_pa_count_tidy$count)
    
    res_pa_count_tidy
  })


# only selecet two macrohpages
p_proximity_cell_count_distribution_ls2 <-
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      # mutate(distance = as.numeric(distance)) %>%
      mutate(distance = factor(distance, levels = c("10", "50", "100", "150", "200"))) %>%
      
      dplyr::filter(Cell %in% c("C1QC Macrophage", "SELENOP Macrophage")) %>%
      
      ggplot(aes(x = Cell, y = count, fill = distance)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ group, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # scale_fill_manual(values = grp_cols_sp) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p 
  })

# # show
# head(res_pa_epi_count_tidy, 10)

p_proximity_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    # p <-
    #   res_pa_count_tidy %>%
    #   mutate(distance = as.numeric(distance)) %>%
    #   
    #   ggplot(aes(x = distance, fill = Group)) + #, color = Group)) +
    #   geom_density(alpha = 0.6) +
    #   # geom_smooth(aes(y = count), method = "loess") +
    #   facet_wrap(~ Cell, scales = "free_y") +
    #   theme_bw() +
    #   scale_fill_manual(values = grp_cols_sp) +
    #   labs(
    #     x = 'Distance(um)',
    #     y = 'Median count'
    #   )
    
    p <-
      res_pa_count_tidy %>%
      group_by(Cell, distance, group) %>%
      summarise(count = median(count)) %>%
      
      mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
      
      ggplot(aes(x = distance, y = count, color = group)) +
      # geom_density(alpha = 0.6) +
      geom_smooth(method = 'loess') +
      facet_wrap(~ Cell, scales = "free_y") +
      ggpubr::theme_pubr() +
      theme(
        strip.background = element_blank(),
        panel.grid = element_blank()
      ) +
      # scale_color_manual(values = grp_cols_sp) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })




# export pdf
pdf(file.path(CN_out_dir, "proximity analysis - barplot.pdf"), width = 14, height = 8)
p_proximity_cell_count_distribution_ls %>% walk(print)
dev.off()

pdf(file.path(CN_out_dir, "proximity analysis - barplot2.pdf"), width = 14, height = 6)
p_proximity_cell_count_distribution_ls2 %>% walk(print)
dev.off()

# # export pdf
# pdf("proximity analysis - bxolplot - combined.pdf", width = 18, height = 8)
# print(p_proximity_cell_count_distribution)
# dev.off()


# export pdf
pdf(file.path(CN_out_dir, "proximity analysis - smooth_lines.pdf"), width = 14, height = 8)
p_proximity_ls %>% walk(print)
dev.off()


# radar plot 

p_proximity_radar_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    p <-
      res_pa_count_tidy %>%
      
      # mutate(distance = factor(distance, levels = c('10', '50', '100', '150', '200' ))) %>%
      mutate(distance = as.numeric(distance)) %>%
      group_by(group, Cell, distance) %>%
      summarise(count = median(count)) %>%
      ggplot(aes(x = Cell, y = distance, size = count, color = group)) +
      geom_vline(xintercept = 0:12 + 0.5, color = "gray") +
      geom_hline(yintercept = c(10, 50, 100, 150, 200), color = "gray") +
      geom_point(position = position_jitterdodge(
        jitter.width = 0,
        jitter.height = 10,
        dodge.width = 0.75
      )) +
      geom_text(
        data = data.frame(x = 0, y = c(10, 50, 100, 150, 200)),
        aes(x = x, y = y, label = paste0(y, "um")),
        inherit.aes = FALSE,
        vjust = -1.5
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(1,1,1,6, 'mm')
      ) +
      scale_size(range = c(3, 12)) +
      scale_y_continuous(limits = c(-100, 220)) +
      # scale_color_manual(values = grp_cols_sp)+
      coord_polar(theta = 'x', clip = 'off') +
      labs(x = NULL, y = NULL, title = name)
    p 
  })





# export pdf
pdf(file.path(CN_out_dir, "proximity analysis - radarplot.pdf"), width = 10, height = 8)
p_proximity_radar_ls %>% walk(print)
dev.off()















# rctd remove erythrocyte --------
# Date <- Sys.Date()
Date <- "2025-09-09"
spatial_out_dir <- str_glue("spatial_analysis_remove_erythrocyte-{Date}")
dir.create(spatial_out_dir)

all_data@meta.data$annotation <- all_data@meta.data$rctd_remove_erythrocyte_20250902
all_data <- subset(all_data_temp, annotation != 'Other')
Idents(all_data) <- "annotation"

sample_mapping <- c(
  "X24M.16.008um" = "24M-1",
  "X24M.16.008um.2" = "24M-2",
  "X4M.7.008um" = "4M-1",
  "X4M.7.008um.4" = "4M-2"
)

cells_position_ls <- 
  all_data@images %>% 
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) %>%
  # setNames(., str_remove(names(.), "\\..*$"))
  setNames(., sample_mapping[names(.)])

cells_position_df <- cells_position_ls %>% bind_rows() 



# CN analysis ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

CN_out_dir <- file.path(spatial_out_dir, "CN analysis")
dir.create(CN_out_dir, recursive = T)


all_data_meta_for_cn_analysis <- 
  all_data[[c("project", "group", "annotation")]] %>%
  dplyr::mutate(annotation_full = annotation ) %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x', 'y')], by = "cell_id") %>%
  rename(X = x, Y = y) %>%
  mutate(
    annotation_full = factor(annotation_full), # import!
    obj_id = cell_id # import!
  )


### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))
res_cn_10 <- readRDS(file.path(CN_out_dir, "res_cn_10.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn_10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn_10_heat$gtable)
dev.off()


#### CN annotaton --------------------
# cn_annotation <- c(
#   "CN1"  = "Immune core",
#   "CN2"  = "Plasma core",
#   "CN3"  = "Tregs immune",
#   "CN4"  = "Mac_Fibro_Ec",
#   "CN5"  = "Myo Stroma",
#   "CN6"  = "Fibro Stroma",
#   "CN7"  = "CD163 Macophage Stroma",
#   "CN8"  = "Enterocytes",
#   "CN9"  = "Myeloid enterocytes",
#   "CN10" = "Amplifying area"
# )

cn_annotation <- c(
  "CN1"  = "CN1",
  "CN2"  = "CN2",
  "CN3"  = "CN3",
  "CN4"  = "CN4",
  "CN5"  = "CN5",
  "CN6"  = "CN6",
  "CN7"  = "CN7",
  "CN8"  = "CN8",
  "CN9"  = "CN9",
  "CN10" = "CN10"
)

p_cn_10_heat2 <-
  cn_cluster_FC %>%
  pheatmap::pheatmap(
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = cn_annotation,
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap_annotation.pdf"),
  height = 6,
  width  = 11
)
cowplot::plot_grid(p_cn_10_heat2$gtable)
dev.off()




### CN k15 -----
res_cn_15 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 15
)

saveRDS(res_cn_15, file.path(CN_out_dir, "res_cn_15.rds"))
res_cn_15 <- readRDS(file.path(CN_out_dir, "res_cn_15.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn15_cluster_mean <- res_cn_15$cluster_centers
cn15_baseMean     <- res_cn_15$base_mean

# log2(FC)
cn15_cluster_FC <- log2(1 + sweep(cn15_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn15_baseMean))

# heatmap
p_cn15_heat <-
  pheatmap::pheatmap(
    cn15_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn15_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_15_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn15_heat$gtable)
dev.off()




### spatial plot -----------------
sp_out_dir <- file.path(CN_out_dir, 'spatial plot') 
dir.create(sp_out_dir)


p_sp_plot_anno <-
  # Images(all_data) %>%
  # setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  # imap(function(img, img_name){
  sample_mapping %>% 
  imap(function(img_name, img){
    # print(img_name)
    # print(img)
    # return()
    p <- SpatialDimPlot(
      all_data,
      images = img, 
      group.by = "annotation",
      pt.size.factor	= 24
    ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

# export
pdf(
  file.path(sp_out_dir, "spatial plot annotation.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_anno %>% walk(print)
dev.off()

### CN spatial plot -----------------
cn_sp_out_dir <- file.path(CN_out_dir, 'CN spatial plot') 
dir.create(cn_sp_out_dir)

# cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

# cn10 
all_data_meta_cn10 <- 
  all_data[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn10, file.path(CN_out_dir, 'all_data_meta_cn10.rds') )
all_data_meta_cn10 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn10.rds'))

cn10_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


# cn15 
all_data_meta_cn15 <- 
  all_data[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn15, file.path(CN_out_dir, 'all_data_meta_cn15.rds') )
all_data_meta_cn15 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn15.rds'))

cn15_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn15, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn15_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN15 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

# ### CN2 spatial plot -----------------
# cn_sp_out_dir <- file.path(CN_out_dir, 'cn spatial plot') 
# dir.create(cn_sp_out_dir)
# 
# # cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# # image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')
# 
# all_data_meta <- 
#   all_data[[]] %>%
#   rownames_to_column('cell_id') %>%
#   left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
#   mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
# 
# saveRDS(all_data_meta, file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# 
# cn_sp_p_ls <-
#   cells_position_ls %>% 
#   imap(function(cells_position, project){ # browser()
#     
#     # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
#     # group_info  <- anno_data$group[anno_data$project == project]
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = y, y = x)) +
#       geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#       scale_y_reverse() +
#       coord_equal() +
#       ggsci::scale_fill_igv() +
#       theme_void() +
#       labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
#     p
#   })
# 
# # export image
# cn_sp_p_ls %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir, str_glue("CN spatial plot {name}.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })


### CN spatial plot2 -----------------
# cn_sp_out_dir2 <- file.path(CN_out_dir, 'cn spatial plot2') 
# dir.create(cn_sp_out_dir2)
# 
# all_data_meta <- readRDS( file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# cn_sp_p_ls2 <-
#   cells_position_ls %>% # .[1] %>%
#   imap(function(cells_position, project){ # browser()
#     cn_color <- structure(
#       names = paste0("CN", 1:10),
#       ggsci::pal_igv()(10)
#     ) %>%
#       c(
#         Others = "gray90"
#       )
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p_ls <-
#       paste0("CN", 1:10) %>% # .[1] %>%
#       set_names(., .) %>%
#       map(function(cn){
#         dat2 <- dat
#         
#         dat2$CN2 <- as.character(dat2$CN2)
#         dat2$CN2[dat2$CN2 != cn] <- "Others"
#         
#         p <- dat2 %>%
#           ggplot(aes(x = y, y = x)) +
#           geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#           scale_y_reverse() +
#           coord_equal() +
#           # ggsci::scale_fill_igv() +
#           scale_fill_manual(values = cn_color) +
#           theme_void() +
#           labs(fill = NULL) 
#         p
#       })
#     
#     p_ls
#   })

# # export image
# cn_sp_p_ls2 %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir2, str_glue("CN spatial plot {name} - split by CN.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })

### cn10 cell props ---------
all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn10_props <- 
  all_cn10_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

# all_cn_count_all_zeros_df <-
#   expand.grid(
#     CN2     = all_cn_count_df$CN2 %>% unique(),
#     project = all_cn_count_df$project %>% unique(),
#     group   = all_cn_count_df$group %>% unique()
#   )

p_cn10_props_group <- 
  all_cn10_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  # stat_compare_means(label = "p.format", size=3, aes(group = group))+ 
  stat_compare_means(method = "t.test",label = "p.format", size=3)+ 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn10_props %>% print
p_cn10_props_group %>% print
dev.off()


### cn15 cell props ---------
all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn15_props <- 
  all_cn15_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

p_cn15_props_group <- 
  all_cn15_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  # stat_compare_means(label = "p.format", size=3, aes(group = group))+ 
  stat_compare_means(label = "p.format", size=3) + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN15_boxplot_and_barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn15_props %>% print
p_cn15_props_group %>% print
dev.off()













## 修改图片CN10 --------
all_data_meta_cn10 <- readRDS('spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/all_data_meta_cn10.rds')

all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn10_count_df %>% head()


# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn10_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/shapiro_test_results_CN10.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn10_count_df, CN == "CN1")

# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.07, 0.13),
    breaks = seq(0.07, 0.13, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn10_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.11),
    breaks = seq(0.01, 0.11, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn10_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.1, 0.2),
    breaks = seq(0.1, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn10_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.03, 0.13),
    breaks = seq(0.03, 0.13, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN4.pdf", p_detailed, width = 10, height = 8)





plot_data_CN5  <- subset(all_cn10_count_df, CN == "CN5")
print(plot_data_CN5)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.12, 0.2),
    breaks = seq(0.12, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN5",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN5.pdf", p_detailed, width = 10, height = 8)





plot_data_CN6  <- subset(all_cn10_count_df, CN == "CN6")
print(plot_data_CN6)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.06, 0.12),
    breaks = seq(0.06, 0.12, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN6",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn10_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.12, 0.16),
    breaks = seq(0.12, 0.16, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn10_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.06),
    breaks = seq(0.01, 0.06, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn10_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.08, 0.13),
    breaks = seq(0.08, 0.13, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn10_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN10-CN10.pdf", p_detailed, width = 10, height = 8)








## 修改图片CN15 --------
all_data_meta_cn15 <- readRDS('spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/all_data_meta_cn15.rds')

all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn15_count_df %>% head()


# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn15_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/shapiro_test_results_CN15.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn15_count_df, CN == "CN1")
print(plot_data_CN1)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.14),
    breaks = seq(0.04, 0.14, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn15_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn15_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.09, 0.12),
    breaks = seq(0.09, 0.12, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn15_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.02, 0.14),
    breaks = seq(0.02, 0.14, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN4.pdf", p_detailed, width = 10, height = 8)





plot_data_CN5  <- subset(all_cn15_count_df, CN == "CN5")
print(plot_data_CN5)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.08, 0.15),
    breaks = seq(0.08, 0.15, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN5",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN5.pdf", p_detailed, width = 10, height = 8)





plot_data_CN6  <- subset(all_cn15_count_df, CN == "CN6")
print(plot_data_CN6)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.1),
    breaks = seq(0.04, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN6",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn15_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.08),
    breaks = seq(0.04, 0.08, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN7.pdf", p_detailed, width = 10, height = 8)





# plot_data_CN8  <- subset(all_cn15_count_df, CN == "CN8")
# print(plot_data_CN8)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0.075, 0.12),
#     breaks = seq(0.075, 0.12, by = 0.01)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN8",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn15_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.1),
    breaks = seq(0.04, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn15_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.03, 0.1),
    breaks = seq(0.03, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN10.pdf", p_detailed, width = 10, height = 8)






plot_data_CN11  <- subset(all_cn15_count_df, CN == "CN11")
print(plot_data_CN11)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN11, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.07),
    breaks = seq(0.01, 0.07, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN11",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN11.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN11.pdf", p_detailed, width = 10, height = 8)






plot_data_CN12  <- subset(all_cn15_count_df, CN == "CN12")
print(plot_data_CN12)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN12, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN12",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN12.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN12.pdf", p_detailed, width = 10, height = 8)






plot_data_CN13  <- subset(all_cn15_count_df, CN == "CN13")
print(plot_data_CN13)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN13, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.06),
    breaks = seq(0, 0.06, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN13",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN13.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN13.pdf", p_detailed, width = 10, height = 8)






plot_data_CN14  <- subset(all_cn15_count_df, CN == "CN14")
print(plot_data_CN14)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN14, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.1),
    breaks = seq(0.05, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN14",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN14.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN14.pdf", p_detailed, width = 10, height = 8)






# plot_data_CN15  <- subset(all_cn15_count_df, CN == "CN15")
# print(plot_data_CN15)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN15, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0.09, 0.14),
#     breaks = seq(0.09, 0.14, by = 0.01)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN15",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN15.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2025-09-09/CN analysis/CN15-CN15.pdf", p_detailed, width = 10, height = 8)









# Proximity analysis -------
### 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

All_type_Proximity_out_dir <- file.path(spatial_out_dir, "Proximity analysis")
dir.create(All_type_Proximity_out_dir)

cells_position_ls <- 
  all_data@images %>%
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) %>% 
  # setNames(., str_remove(names(.), "\\..*$"))
  setNames(., sample_mapping[names(.)])



# all_data_cn_meta
all_data_meta_pos <- 
  cells_position_ls %>% 
  imap(function(cells_position, slide){# browser()
    all_data_meta <- 
      all_data[[]] %>% 
      dplyr::mutate(annotation_full = annotation)
    proj_meta_sub <- 
      all_data_meta %>% 
      dplyr::select(annotation_full, group, project) %>% 
      rownames_to_column("cell_id") %>% 
      dplyr::filter(project == slide)
    cells_position[, 1:3] %>% 
      inner_join(proj_meta_sub, by = "cell_id") %>% 
      mutate(annotation_full = factor(annotation_full)) %>% 
      dplyr::rename(Classification1 = annotation_full, 
                    X = x, 
                    Y = y, 
                    Slide_ID = project)
  }) %>% 
  enframe() %>% 
  unnest(cols = value) %>% 
  dplyr::select(-name) %>% 
  column_to_rownames("cell_id")

all_data_meta_pos$Slide_ID <- str_c(all_data_meta_pos$group, all_data_meta_pos$Slide_ID, sep = "_")

all_data_meta_pos_subset <- all_data_meta_pos %>% rownames_to_column("cell_id")
# readRDS('data/all_data_meta_pos_subset.rds')

all_data_meta_pos_subset %>% head

# cell_id        X        Y Classification1 group       Slide_ID
# 1 1_aaaaaagg-1 1670.389 185.1358      Macrophage  CD-F X_202413076_32
# 2 1_aaadlghn-1 1681.172 176.4662          Plasma  CD-F X_202413076_32
# 3 1_aaadmplc-1 1677.561 182.6699               T  CD-F X_202413076_32
# 4 1_aaadoogj-1 1677.109 165.6637          Plasma  CD-F X_202413076_32
# 5 1_aaadpljf-1 1671.887 176.8148          Plasma  CD-F X_202413076_32
# 6 1_aaafaiko-1 1719.466 192.8071          Plasma  CD-F X_202413076_32
# ......


### 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

#### 2.1 proximity analysis --------
res_pa <- 
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset,
    query = all_data_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 500, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 8 / all_data@images$X24M.16.008um@scale.factors$spot, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
saveRDS(res_pa, file = file.path(All_type_Proximity_out_dir, "res_pa.rds"))
# res_pa <- readRDS(file.path(All_type_Proximity_out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')

# summarized by group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos_subset,
#     str_c(Classification1, group, Slide_ID, sep = '..')
#   )
pa_groups <- 
  all_data_meta_pos_subset %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE) 

pa_groups2 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn10[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

pa_groups3 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn15[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
res_pa_df2 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups2, FUN = "median")})
res_pa_df3 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups3, FUN = "median")})


res_pa_df  <- 
  res_pa_df %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project"), sep = "\\.\\.") 

res_pa_df2  <- 
  res_pa_df2 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.")
res_pa_df2$CN <- factor(res_pa_df2$CN, levels = paste0("CN", 1:10))

res_pa_df3  <- 
  res_pa_df3 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df3$CN <- factor(res_pa_df3$CN, levels = paste0("CN", 1:15))


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_res_pa_ls_cn10 <-
  res_pa_df2 %>%
  with(split(res_pa_df2, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

p_res_pa_ls_cn15 <-
  res_pa_df3 %>%
  with(split(res_pa_df3, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

# p_res_pa_ls <-
#   res_pa_df2 %>%
#   with(split(res_pa_df2, CN2)) %>%
#   imap(function(dat, name){
#     dat %>%
#   # dplyr::filter(Target == "Epithelial cells") %>%
#   ggplot(aes(x = Others, y = cell_props, fill = Group)) +
#   geom_boxplot(width = 0.45, outlier.colour = NA) +
#   # geom_jitter(width = 0.25) +
#   facet_wrap(~ Target, scales = "free_y", ncol = 2)+
#   ggpubr::theme_pubclean() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   ggtitle(name)
#   })

# export plot
pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis.pdf"),
  height = 50,
  width  = 25
)
p_res_pa%>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN10.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn10 %>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN15.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn15 %>% print()
dev.off()










# 修改图片
# res_pa_df

shapiro_results <- res_pa_df %>%
  group_by(Target) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(cell_props)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(cell_props)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2025-09-09/Proximity analysis/shapiro_test_results.csv", row.names = FALSE)












# <------------
# <------------
# <------------
# <------------
# rctd remove erythrocyte --------
# Date <- Sys.Date()
Date <- "2025-11-21"
spatial_out_dir <- str_glue("spatial_analysis_remove_erythrocyte-{Date}")
dir.create(spatial_out_dir)

all_data <- readRDS("data_for_analysis/all_data_projected_rctd_20250902.rds")
all_data@meta.data$annotation <- all_data@meta.data$rctd_remove_erythrocyte_20250902
unique(all_data$annotation)
# all_data <- subset(all_data, annotation != 'Other')
all_data <- all_data[, !all_data@meta.data$annotation %in% c("Other", "Granulocytes-sub4", "Synovial fibroblasts-sub2", "Synovial intimal fibroblasts-sub10")]

Idents(all_data) <- "annotation"

sample_mapping <- c(
  "X24M.16.008um" = "24M-1",
  "X24M.16.008um.2" = "24M-2",
  "X4M.7.008um" = "4M-1",
  "X4M.7.008um.4" = "4M-2"
)

cells_position_ls <- 
  all_data@images %>% 
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) %>%
  # setNames(., str_remove(names(.), "\\..*$"))
  setNames(., sample_mapping[names(.)])

cells_position_df <- cells_position_ls %>% bind_rows() 



# CN analysis ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

CN_out_dir <- file.path(spatial_out_dir, "CN analysis")
dir.create(CN_out_dir, recursive = T)


all_data_meta_for_cn_analysis <- 
  all_data[[c("project", "group", "annotation")]] %>%
  dplyr::mutate(annotation_full = annotation ) %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x', 'y')], by = "cell_id") %>%
  rename(X = x, Y = y) %>%
  mutate(
    annotation_full = factor(annotation_full), # import!
    obj_id = cell_id # import!
  )


### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))
# res_cn_10 <- readRDS(file.path(CN_out_dir, "res_cn_10.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn_10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn_10_heat$gtable)
dev.off()


#### CN annotaton --------------------
# cn_annotation <- c(
#   "CN1"  = "Immune core",
#   "CN2"  = "Plasma core",
#   "CN3"  = "Tregs immune",
#   "CN4"  = "Mac_Fibro_Ec",
#   "CN5"  = "Myo Stroma",
#   "CN6"  = "Fibro Stroma",
#   "CN7"  = "CD163 Macophage Stroma",
#   "CN8"  = "Enterocytes",
#   "CN9"  = "Myeloid enterocytes",
#   "CN10" = "Amplifying area"
# )

cn_annotation <- c(
  "CN1"  = "CN1",
  "CN2"  = "CN2",
  "CN3"  = "CN3",
  "CN4"  = "CN4",
  "CN5"  = "CN5",
  "CN6"  = "CN6",
  "CN7"  = "CN7",
  "CN8"  = "CN8",
  "CN9"  = "CN9",
  "CN10" = "CN10"
)

p_cn_10_heat2 <-
  cn_cluster_FC %>%
  pheatmap::pheatmap(
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = cn_annotation,
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap_annotation.pdf"),
  height = 6,
  width  = 11
)
cowplot::plot_grid(p_cn_10_heat2$gtable)
dev.off()




### CN k15 -----
res_cn_15 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 15
)

saveRDS(res_cn_15, file.path(CN_out_dir, "res_cn_15.rds"))
# res_cn_15 <- readRDS(file.path(CN_out_dir, "res_cn_15.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn15_cluster_mean <- res_cn_15$cluster_centers
cn15_baseMean     <- res_cn_15$base_mean

# log2(FC)
cn15_cluster_FC <- log2(1 + sweep(cn15_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn15_baseMean))

# heatmap
p_cn15_heat <-
  pheatmap::pheatmap(
    cn15_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn15_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_15_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn15_heat$gtable)
dev.off()




### spatial plot -----------------
sp_out_dir <- file.path(CN_out_dir, 'spatial plot') 
dir.create(sp_out_dir)


p_sp_plot_anno <-
  # Images(all_data) %>%
  # setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  # imap(function(img, img_name){
  sample_mapping %>% 
  imap(function(img_name, img){
    # print(img_name)
    # print(img)
    # return()
    p <- SpatialDimPlot(
      all_data,
      images = img, 
      group.by = "annotation",
      pt.size.factor	= 24
    ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

# export
pdf(
  file.path(sp_out_dir, "spatial plot annotation.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_anno %>% walk(print)
dev.off()

### CN spatial plot -----------------
cn_sp_out_dir <- file.path(CN_out_dir, 'CN spatial plot') 
dir.create(cn_sp_out_dir)

# cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

# cn10 
all_data_meta_cn10 <- 
  all_data[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn10, file.path(CN_out_dir, 'all_data_meta_cn10.rds') )
# all_data_meta_cn10 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn10.rds'))

cn10_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


# cn15 
all_data_meta_cn15 <- 
  all_data[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn15, file.path(CN_out_dir, 'all_data_meta_cn15.rds') )
# all_data_meta_cn15 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn15.rds'))

cn15_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn15, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn15_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN15 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

# ### CN2 spatial plot -----------------
# cn_sp_out_dir <- file.path(CN_out_dir, 'cn spatial plot') 
# dir.create(cn_sp_out_dir)
# 
# # cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# # image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')
# 
# all_data_meta <- 
#   all_data[[]] %>%
#   rownames_to_column('cell_id') %>%
#   left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
#   mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
# 
# saveRDS(all_data_meta, file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# 
# cn_sp_p_ls <-
#   cells_position_ls %>% 
#   imap(function(cells_position, project){ # browser()
#     
#     # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
#     # group_info  <- anno_data$group[anno_data$project == project]
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = y, y = x)) +
#       geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#       scale_y_reverse() +
#       coord_equal() +
#       ggsci::scale_fill_igv() +
#       theme_void() +
#       labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
#     p
#   })
# 
# # export image
# cn_sp_p_ls %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir, str_glue("CN spatial plot {name}.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })


### CN spatial plot2 -----------------
# cn_sp_out_dir2 <- file.path(CN_out_dir, 'cn spatial plot2') 
# dir.create(cn_sp_out_dir2)
# 
# all_data_meta <- readRDS( file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# cn_sp_p_ls2 <-
#   cells_position_ls %>% # .[1] %>%
#   imap(function(cells_position, project){ # browser()
#     cn_color <- structure(
#       names = paste0("CN", 1:10),
#       ggsci::pal_igv()(10)
#     ) %>%
#       c(
#         Others = "gray90"
#       )
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p_ls <-
#       paste0("CN", 1:10) %>% # .[1] %>%
#       set_names(., .) %>%
#       map(function(cn){
#         dat2 <- dat
#         
#         dat2$CN2 <- as.character(dat2$CN2)
#         dat2$CN2[dat2$CN2 != cn] <- "Others"
#         
#         p <- dat2 %>%
#           ggplot(aes(x = y, y = x)) +
#           geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#           scale_y_reverse() +
#           coord_equal() +
#           # ggsci::scale_fill_igv() +
#           scale_fill_manual(values = cn_color) +
#           theme_void() +
#           labs(fill = NULL) 
#         p
#       })
#     
#     p_ls
#   })

# # export image
# cn_sp_p_ls2 %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir2, str_glue("CN spatial plot {name} - split by CN.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })

### cn10 cell props ---------
all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn10_props <- 
  all_cn10_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

# all_cn_count_all_zeros_df <-
#   expand.grid(
#     CN2     = all_cn_count_df$CN2 %>% unique(),
#     project = all_cn_count_df$project %>% unique(),
#     group   = all_cn_count_df$group %>% unique()
#   )

p_cn10_props_group <- 
  all_cn10_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(method = "t.test",label = "p.format", size=3)+ 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn10_props %>% print
p_cn10_props_group %>% print
dev.off()


### cn15 cell props ---------
all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn15_props <- 
  all_cn15_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

p_cn15_props_group <- 
  all_cn15_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(label = "p.format", size=3) + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN15_boxplot_and_barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn15_props %>% print
p_cn15_props_group %>% print
dev.off()













## 修改图片CN10 --------
all_data_meta_cn10 <- readRDS('spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/all_data_meta_cn10.rds')

all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn10_count_df %>% head()

write_csv(all_cn10_count_df, file = "spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/all_cn10_count_df.csv", col_names = T, )

# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn10_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/shapiro_test_results_CN10.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn10_count_df, CN == "CN1")
print(plot_data_CN1)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.03, 0.12),
    breaks = seq(0.03, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn10_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.1, 0.17),
    breaks = seq(0.1, 0.17, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn10_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn10_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.03, 0.13),
    breaks = seq(0.03, 0.13, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN4.pdf", p_detailed, width = 10, height = 8)





plot_data_CN5  <- subset(all_cn10_count_df, CN == "CN5")
print(plot_data_CN5)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.12, 0.2),
    breaks = seq(0.12, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN5",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN5.pdf", p_detailed, width = 10, height = 8)





plot_data_CN6  <- subset(all_cn10_count_df, CN == "CN6")
print(plot_data_CN6)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.06, 0.12),
    breaks = seq(0.06, 0.12, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN6",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn10_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.12, 0.16),
    breaks = seq(0.12, 0.16, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn10_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.06),
    breaks = seq(0.01, 0.06, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn10_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.08, 0.13),
    breaks = seq(0.08, 0.13, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn10_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN10-CN10.pdf", p_detailed, width = 10, height = 8)








## 修改图片CN15 --------
all_data_meta_cn15 <- readRDS('spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/all_data_meta_cn15.rds')

all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn15_count_df %>% head()


# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn15_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/shapiro_test_results_CN15.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn15_count_df, CN == "CN1")
print(plot_data_CN1)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.14),
    breaks = seq(0.04, 0.14, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn15_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn15_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.09, 0.12),
    breaks = seq(0.09, 0.12, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn15_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.02, 0.14),
    breaks = seq(0.02, 0.14, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN4.pdf", p_detailed, width = 10, height = 8)





plot_data_CN5  <- subset(all_cn15_count_df, CN == "CN5")
print(plot_data_CN5)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.08, 0.15),
    breaks = seq(0.08, 0.15, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN5",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN5.pdf", p_detailed, width = 10, height = 8)





plot_data_CN6  <- subset(all_cn15_count_df, CN == "CN6")
print(plot_data_CN6)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.1),
    breaks = seq(0.04, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN6",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn15_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.08),
    breaks = seq(0.04, 0.08, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN7.pdf", p_detailed, width = 10, height = 8)





# plot_data_CN8  <- subset(all_cn15_count_df, CN == "CN8")
# print(plot_data_CN8)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0.075, 0.12),
#     breaks = seq(0.075, 0.12, by = 0.01)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN8",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn15_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.1),
    breaks = seq(0.04, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn15_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.03, 0.1),
    breaks = seq(0.03, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN10.pdf", p_detailed, width = 10, height = 8)






plot_data_CN11  <- subset(all_cn15_count_df, CN == "CN11")
print(plot_data_CN11)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN11, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.07),
    breaks = seq(0.01, 0.07, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN11",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN11.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN11.pdf", p_detailed, width = 10, height = 8)






plot_data_CN12  <- subset(all_cn15_count_df, CN == "CN12")
print(plot_data_CN12)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN12, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN12",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN12.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN12.pdf", p_detailed, width = 10, height = 8)






plot_data_CN13  <- subset(all_cn15_count_df, CN == "CN13")
print(plot_data_CN13)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN13, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.06),
    breaks = seq(0, 0.06, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN13",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN13.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN13.pdf", p_detailed, width = 10, height = 8)






plot_data_CN14  <- subset(all_cn15_count_df, CN == "CN14")
print(plot_data_CN14)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN14, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.1),
    breaks = seq(0.05, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN14",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN14.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN14.pdf", p_detailed, width = 10, height = 8)






# plot_data_CN15  <- subset(all_cn15_count_df, CN == "CN15")
# print(plot_data_CN15)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN15, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0.09, 0.14),
#     breaks = seq(0.09, 0.14, by = 0.01)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN15",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN15.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2025-11-21/CN analysis/CN15-CN15.pdf", p_detailed, width = 10, height = 8)









# Proximity analysis -------
### 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

All_type_Proximity_out_dir <- file.path(spatial_out_dir, "Proximity analysis")
dir.create(All_type_Proximity_out_dir)

cells_position_ls <- 
  all_data@images %>%
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) %>% 
  # setNames(., str_remove(names(.), "\\..*$"))
  setNames(., sample_mapping[names(.)])



# all_data_cn_meta
all_data_meta_pos <- 
  cells_position_ls %>% 
  imap(function(cells_position, slide){# browser()
    all_data_meta <- 
      all_data[[]] %>% 
      dplyr::mutate(annotation_full = annotation)
    proj_meta_sub <- 
      all_data_meta %>% 
      dplyr::select(annotation_full, group, project) %>% 
      rownames_to_column("cell_id") %>% 
      dplyr::filter(project == slide)
    cells_position[, 1:3] %>% 
      inner_join(proj_meta_sub, by = "cell_id") %>% 
      mutate(annotation_full = factor(annotation_full)) %>% 
      dplyr::rename(Classification1 = annotation_full, 
                    X = x, 
                    Y = y, 
                    Slide_ID = project)
  }) %>% 
  enframe() %>% 
  unnest(cols = value) %>% 
  dplyr::select(-name) %>% 
  column_to_rownames("cell_id")

all_data_meta_pos$Slide_ID <- str_c(all_data_meta_pos$group, all_data_meta_pos$Slide_ID, sep = "_")

all_data_meta_pos_subset <- all_data_meta_pos %>% rownames_to_column("cell_id")
# readRDS('data/all_data_meta_pos_subset.rds')

all_data_meta_pos_subset %>% head

# cell_id        X        Y Classification1 group       Slide_ID
# 1 1_aaaaaagg-1 1670.389 185.1358      Macrophage  CD-F X_202413076_32
# 2 1_aaadlghn-1 1681.172 176.4662          Plasma  CD-F X_202413076_32
# 3 1_aaadmplc-1 1677.561 182.6699               T  CD-F X_202413076_32
# 4 1_aaadoogj-1 1677.109 165.6637          Plasma  CD-F X_202413076_32
# 5 1_aaadpljf-1 1671.887 176.8148          Plasma  CD-F X_202413076_32
# 6 1_aaafaiko-1 1719.466 192.8071          Plasma  CD-F X_202413076_32
# ......


### 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

#### 2.1 proximity analysis --------
res_pa <- 
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset,
    query = all_data_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 500, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 8 / all_data@images$X24M.16.008um@scale.factors$spot, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
saveRDS(res_pa, file = file.path(All_type_Proximity_out_dir, "res_pa.rds"))
# res_pa <- readRDS(file.path(All_type_Proximity_out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')

# summarized by group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos_subset,
#     str_c(Classification1, group, Slide_ID, sep = '..')
#   )
pa_groups <- 
  all_data_meta_pos_subset %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE) 

pa_groups2 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn10[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

pa_groups3 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn15[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
res_pa_df2 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups2, FUN = "median")})
res_pa_df3 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups3, FUN = "median")})


res_pa_df  <- 
  res_pa_df %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project"), sep = "\\.\\.") 

res_pa_df2  <- 
  res_pa_df2 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.")
res_pa_df2$CN <- factor(res_pa_df2$CN, levels = paste0("CN", 1:10))

res_pa_df3  <- 
  res_pa_df3 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df3$CN <- factor(res_pa_df3$CN, levels = paste0("CN", 1:15))


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_res_pa_ls_cn10 <-
  res_pa_df2 %>%
  with(split(res_pa_df2, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

p_res_pa_ls_cn15 <-
  res_pa_df3 %>%
  with(split(res_pa_df3, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

# p_res_pa_ls <-
#   res_pa_df2 %>%
#   with(split(res_pa_df2, CN2)) %>%
#   imap(function(dat, name){
#     dat %>%
#   # dplyr::filter(Target == "Epithelial cells") %>%
#   ggplot(aes(x = Others, y = cell_props, fill = Group)) +
#   geom_boxplot(width = 0.45, outlier.colour = NA) +
#   # geom_jitter(width = 0.25) +
#   facet_wrap(~ Target, scales = "free_y", ncol = 2)+
#   ggpubr::theme_pubclean() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   ggtitle(name)
#   })

# export plot
pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis.pdf"),
  height = 50,
  width  = 25
)
p_res_pa%>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN10.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn10 %>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN15.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn15 %>% print()
dev.off()










# 修改图片
# res_pa_df

shapiro_results <- res_pa_df %>%
  group_by(Target) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(cell_props)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(cell_props)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2025-09-09/Proximity analysis/shapiro_test_results.csv", row.names = FALSE)

















# <------------
# <------------
# <------------
# <------------
# rctd remove erythrocyte --------
# Date <- Sys.Date()
Date <- "2026-01-06"
spatial_out_dir <- str_glue("spatial_analysis_remove_erythrocyte-{Date}")
dir.create(spatial_out_dir)

all_data <- readRDS("data_for_analysis/all_data_filt.rds")
all_data@meta.data$annotation <- all_data@meta.data$rctd_20251201
unique(all_data$annotation)
# all_data_filt <- subset(all_data, annotation != 'Other')
all_data_filt <- all_data[, !all_data@meta.data$annotation %in% c("Other")]

# 只保留原始图像，移除自动生成的重复项
original_images <- c("24M-1", "24M-2", "4M-1", "4M-2")
images_to_keep <- names(all_data_filt@images)[names(all_data_filt@images) %in% original_images]

# 如果原始图像都存在，直接子集化
if(length(images_to_keep) == 4) {
  all_data_filt@images <- all_data_filt@images[images_to_keep]
  cat("修复完成！当前图像名称:", paste(names(all_data_filt@images), collapse = ", "), "\n")
} else {
  # 如果原始图像不全，使用更安全的方法
  cat("采取备用方案...\n")
  # 方法：基于图像数据内容去重（更安全）
  image_data_list <- list()
  for(img_name in names(all_data_filt@images)) {
    img <- all_data_filt@images[[img_name]]
    # 使用关键标识来去重
    img_key <- paste(img@assay, img@key, sep = "_")
    if(!img_key %in% names(image_data_list)) {
      image_data_list[[img_key]] <- img
    }
  }
  all_data_filt@images <- image_data_list
}

Idents(all_data_filt) <- "annotation"

sample_mapping <- c(
  "X24M.16.008um" = "24M-1",
  "X24M.16.008um.2" = "24M-2",
  "X4M.7.008um" = "4M-1",
  "X4M.7.008um.4" = "4M-2"
)

cells_position_ls <- 
  all_data_filt@images %>% 
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) # %>%
  # setNames(., str_remove(names(.), "\\..*$"))
  # setNames(., sample_mapping[names(.)])

cells_position_df <- cells_position_ls %>% bind_rows() 






# CN analysis knn=25---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

CN_out_dir <- file.path(spatial_out_dir, "CN analysis - knn25")
dir.create(CN_out_dir, recursive = T)


all_data_meta_for_cn_analysis <- 
  all_data_filt[[c("project", "group", "annotation")]] %>%
  dplyr::mutate(annotation_full = annotation ) %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x', 'y')], by = "cell_id") %>%
  rename(X = x, Y = y) %>%
  mutate(
    annotation_full = factor(annotation_full), # import!
    obj_id = cell_id # import!
  )


### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 25,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))
# res_cn_10 <- readRDS(file.path(CN_out_dir, "res_cn_10.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn_10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn_10_heat$gtable)
dev.off()


#### CN annotaton --------------------
# cn_annotation <- c(
#   "CN1"  = "Immune core",
#   "CN2"  = "Plasma core",
#   "CN3"  = "Tregs immune",
#   "CN4"  = "Mac_Fibro_Ec",
#   "CN5"  = "Myo Stroma",
#   "CN6"  = "Fibro Stroma",
#   "CN7"  = "CD163 Macophage Stroma",
#   "CN8"  = "Enterocytes",
#   "CN9"  = "Myeloid enterocytes",
#   "CN10" = "Amplifying area"
# )

# cn_annotation <- c(
#   "CN1"  = "CN1",
#   "CN2"  = "CN2",
#   "CN3"  = "CN3",
#   "CN4"  = "CN4",
#   "CN5"  = "CN5",
#   "CN6"  = "CN6",
#   "CN7"  = "CN7",
#   "CN8"  = "CN8",
#   "CN9"  = "CN9",
#   "CN10" = "CN10"
# )
# 
# p_cn_10_heat2 <-
#   cn_cluster_FC %>%
#   pheatmap::pheatmap(
#     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#     border_color = NA,
#     treeheight_row = 0,
#     labels_row = cn_annotation,
#     # clustering_method = 'ward.D2',
#     cluster_rows = F,
#     scale = "column",
#     silent = TRUE
#   )
# 
# pdf(
#   file.path(CN_out_dir, "CN_10_heatmap_annotation.pdf"),
#   height = 6,
#   width  = 11
# )
# cowplot::plot_grid(p_cn_10_heat2$gtable)
# dev.off()
# 



### CN k15 -----
res_cn_15 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 25,
  k.cluster      = 15
)

saveRDS(res_cn_15, file.path(CN_out_dir, "res_cn_15.rds"))
# res_cn_15 <- readRDS(file.path(CN_out_dir, "res_cn_15.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn15_cluster_mean <- res_cn_15$cluster_centers
cn15_baseMean     <- res_cn_15$base_mean

# log2(FC)
cn15_cluster_FC <- log2(1 + sweep(cn15_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn15_baseMean))

# heatmap
p_cn15_heat <-
  pheatmap::pheatmap(
    cn15_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn15_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_15_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn15_heat$gtable)
dev.off()




### spatial plot -----------------
sp_out_dir <- file.path(CN_out_dir, 'spatial plot') 
dir.create(sp_out_dir)


p_sp_plot_anno <-
  # Images(all_data) %>%
  # setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  # imap(function(img, img_name){
  sample_mapping %>% 
  imap(function(img_name, img){
    # print(img_name)
    # print(img)
    # return()
    p <- SpatialDimPlot(
      all_data_filt,
      images = img_name, 
      group.by = "annotation",
      pt.size.factor	= 24
    ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

# export
pdf(
  file.path(sp_out_dir, "spatial plot annotation.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_anno %>% walk(print)
dev.off()

### CN spatial plot -----------------
cn_sp_out_dir <- file.path(CN_out_dir, 'CN spatial plot') 
dir.create(cn_sp_out_dir)

# cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

# cn10 
all_data_meta_cn10 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn10, file.path(CN_out_dir, 'all_data_meta_cn10.rds') )
# all_data_meta_cn10 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn10.rds'))

cn10_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


# cn15 
all_data_meta_cn15 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn15, file.path(CN_out_dir, 'all_data_meta_cn15.rds') )
# all_data_meta_cn15 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn15.rds'))

cn15_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn15, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn15_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN15 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

# ### CN2 spatial plot -----------------
# cn_sp_out_dir <- file.path(CN_out_dir, 'cn spatial plot') 
# dir.create(cn_sp_out_dir)
# 
# # cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# # image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')
# 
# all_data_meta <- 
#   all_data[[]] %>%
#   rownames_to_column('cell_id') %>%
#   left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
#   mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
# 
# saveRDS(all_data_meta, file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# 
# cn_sp_p_ls <-
#   cells_position_ls %>% 
#   imap(function(cells_position, project){ # browser()
#     
#     # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
#     # group_info  <- anno_data$group[anno_data$project == project]
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = y, y = x)) +
#       geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#       scale_y_reverse() +
#       coord_equal() +
#       ggsci::scale_fill_igv() +
#       theme_void() +
#       labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
#     p
#   })
# 
# # export image
# cn_sp_p_ls %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir, str_glue("CN spatial plot {name}.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })


### CN spatial plot2 -----------------
# cn_sp_out_dir2 <- file.path(CN_out_dir, 'cn spatial plot2') 
# dir.create(cn_sp_out_dir2)
# 
# all_data_meta <- readRDS( file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# cn_sp_p_ls2 <-
#   cells_position_ls %>% # .[1] %>%
#   imap(function(cells_position, project){ # browser()
#     cn_color <- structure(
#       names = paste0("CN", 1:10),
#       ggsci::pal_igv()(10)
#     ) %>%
#       c(
#         Others = "gray90"
#       )
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p_ls <-
#       paste0("CN", 1:10) %>% # .[1] %>%
#       set_names(., .) %>%
#       map(function(cn){
#         dat2 <- dat
#         
#         dat2$CN2 <- as.character(dat2$CN2)
#         dat2$CN2[dat2$CN2 != cn] <- "Others"
#         
#         p <- dat2 %>%
#           ggplot(aes(x = y, y = x)) +
#           geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#           scale_y_reverse() +
#           coord_equal() +
#           # ggsci::scale_fill_igv() +
#           scale_fill_manual(values = cn_color) +
#           theme_void() +
#           labs(fill = NULL) 
#         p
#       })
#     
#     p_ls
#   })

# # export image
# cn_sp_p_ls2 %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir2, str_glue("CN spatial plot {name} - split by CN.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })

### cn10 cell props ---------
all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn10_props <- 
  all_cn10_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

# all_cn_count_all_zeros_df <-
#   expand.grid(
#     CN2     = all_cn_count_df$CN2 %>% unique(),
#     project = all_cn_count_df$project %>% unique(),
#     group   = all_cn_count_df$group %>% unique()
#   )

p_cn10_props_group <- 
  all_cn10_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(method = "t.test",label = "p.format", size=3)+ 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn10_props %>% print
p_cn10_props_group %>% print
dev.off()


### cn15 cell props ---------
all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn15_props <- 
  all_cn15_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

p_cn15_props_group <- 
  all_cn15_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(label = "p.format", size=3) + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN15_boxplot_and_barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn15_props %>% print
p_cn15_props_group %>% print
dev.off()













## 修改图片CN10 --------
all_data_meta_cn10 <- readRDS('spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/all_data_meta_cn10.rds')

all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn10_count_df %>% head()

write_csv(all_cn10_count_df, file = "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/all_cn10_count_df.csv", col_names = T, )

# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn10_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/shapiro_test_results_CN10.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn10_count_df, CN == "CN1")
print(plot_data_CN1)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.15),
    breaks = seq(0, 0.15, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn10_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.2),
    breaks = seq(0.01, 0.2, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn10_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.15),
    breaks = seq(0.05, 0.15, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn10_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.12, 0.25),
    breaks = seq(0.12, 0.25, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN4.pdf", p_detailed, width = 10, height = 8)





plot_data_CN5  <- subset(all_cn10_count_df, CN == "CN5")
print(plot_data_CN5)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.1),
    breaks = seq(0.05, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN5",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN5.pdf", p_detailed, width = 10, height = 8)





plot_data_CN6  <- subset(all_cn10_count_df, CN == "CN6")
print(plot_data_CN6)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.1, 0.2),
    breaks = seq(0.1, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN6",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn10_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.15),
    breaks = seq(0.04, 0.15, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn10_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.12, 0.15),
    breaks = seq(0.12, 0.15, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn10_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.15),
    breaks = seq(0.05, 0.15, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn10_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.08),
    breaks = seq(0, 0.08, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN10.pdf", p_detailed, width = 10, height = 8)








## 修改图片CN15 --------
all_data_meta_cn15 <- readRDS('spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/all_data_meta_cn15.rds')

all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn15_count_df %>% head()


# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn15_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/shapiro_test_results_CN15.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn15_count_df, CN == "CN1")
print(plot_data_CN1)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.1),
    breaks = seq(0.05, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn15_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.12),
    breaks = seq(0, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn15_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn15_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.08, 0.15),
    breaks = seq(0.08, 0.15, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN4.pdf", p_detailed, width = 10, height = 8)





# plot_data_CN5  <- subset(all_cn15_count_df, CN == "CN5")
# print(plot_data_CN5)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0, 0.1),
#     breaks = seq(0, 0.1, by = 0.02)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN5",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN5.pdf", p_detailed, width = 10, height = 8)
# 
# 
# 
# 
# 
# plot_data_CN6  <- subset(all_cn15_count_df, CN == "CN6")
# print(plot_data_CN6)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0.04, 0.12),
#     breaks = seq(0.04, 0.12, by = 0.02)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN6",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn15_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn15_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6,
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15,
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.06, 0.1),
    breaks = seq(0.06, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn15_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.06),
    breaks = seq(0.01, 0.06, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn15_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.07),
    breaks = seq(0, 0.07, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN10.pdf", p_detailed, width = 10, height = 8)






plot_data_CN11  <- subset(all_cn15_count_df, CN == "CN11")
print(plot_data_CN11)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN11, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.1),
    breaks = seq(0.04, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN11",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN11.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN11.pdf", p_detailed, width = 10, height = 8)






plot_data_CN12  <- subset(all_cn15_count_df, CN == "CN12")
print(plot_data_CN12)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN12, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.03, 0.14),
    breaks = seq(0.03, 0.14, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN12",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN12.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN12.pdf", p_detailed, width = 10, height = 8)






plot_data_CN13  <- subset(all_cn15_count_df, CN == "CN13")
print(plot_data_CN13)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN13, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.08),
    breaks = seq(0, 0.08, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN13",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN13.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN13.pdf", p_detailed, width = 10, height = 8)






plot_data_CN14  <- subset(all_cn15_count_df, CN == "CN14")
print(plot_data_CN14)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN14, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN14",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN14.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN14.pdf", p_detailed, width = 10, height = 8)






plot_data_CN15  <- subset(all_cn15_count_df, CN == "CN15")
print(plot_data_CN15)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN15, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6,
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15,
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.08),
    breaks = seq(0.04, 0.08, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN15",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN15.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN15.pdf", p_detailed, width = 10, height = 8)










# CN analysis knn=50---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

CN_out_dir <- file.path(spatial_out_dir, "CN analysis - knn50")
dir.create(CN_out_dir, recursive = T)


all_data_meta_for_cn_analysis <- 
  all_data_filt[[c("project", "group", "annotation")]] %>%
  dplyr::mutate(annotation_full = annotation ) %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x', 'y')], by = "cell_id") %>%
  rename(X = x, Y = y) %>%
  mutate(
    annotation_full = factor(annotation_full), # import!
    obj_id = cell_id # import!
  )


### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))
# res_cn_10 <- readRDS(file.path(CN_out_dir, "res_cn_10.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn_10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn_10_heat$gtable)
dev.off()


#### CN annotaton --------------------
# cn_annotation <- c(
#   "CN1"  = "Immune core",
#   "CN2"  = "Plasma core",
#   "CN3"  = "Tregs immune",
#   "CN4"  = "Mac_Fibro_Ec",
#   "CN5"  = "Myo Stroma",
#   "CN6"  = "Fibro Stroma",
#   "CN7"  = "CD163 Macophage Stroma",
#   "CN8"  = "Enterocytes",
#   "CN9"  = "Myeloid enterocytes",
#   "CN10" = "Amplifying area"
# )

# cn_annotation <- c(
#   "CN1"  = "CN1",
#   "CN2"  = "CN2",
#   "CN3"  = "CN3",
#   "CN4"  = "CN4",
#   "CN5"  = "CN5",
#   "CN6"  = "CN6",
#   "CN7"  = "CN7",
#   "CN8"  = "CN8",
#   "CN9"  = "CN9",
#   "CN10" = "CN10"
# )
# 
# p_cn_10_heat2 <-
#   cn_cluster_FC %>%
#   pheatmap::pheatmap(
#     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#     border_color = NA,
#     treeheight_row = 0,
#     labels_row = cn_annotation,
#     # clustering_method = 'ward.D2',
#     cluster_rows = F,
#     scale = "column",
#     silent = TRUE
#   )
# 
# pdf(
#   file.path(CN_out_dir, "CN_10_heatmap_annotation.pdf"),
#   height = 6,
#   width  = 11
# )
# cowplot::plot_grid(p_cn_10_heat2$gtable)
# dev.off()
# 



### CN k15 -----
res_cn_15 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 25,
  k.cluster      = 15
)

saveRDS(res_cn_15, file.path(CN_out_dir, "res_cn_15.rds"))
# res_cn_15 <- readRDS(file.path(CN_out_dir, "res_cn_15.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn15_cluster_mean <- res_cn_15$cluster_centers
cn15_baseMean     <- res_cn_15$base_mean

# log2(FC)
cn15_cluster_FC <- log2(1 + sweep(cn15_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn15_baseMean))

# heatmap
p_cn15_heat <-
  pheatmap::pheatmap(
    cn15_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn15_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_15_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn15_heat$gtable)
dev.off()




### spatial plot -----------------
sp_out_dir <- file.path(CN_out_dir, 'spatial plot') 
dir.create(sp_out_dir)


p_sp_plot_anno <-
  # Images(all_data) %>%
  # setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  # imap(function(img, img_name){
  sample_mapping %>% 
  imap(function(img_name, img){
    # print(img_name)
    # print(img)
    # return()
    p <- SpatialDimPlot(
      all_data_filt,
      images = img_name, 
      group.by = "annotation",
      pt.size.factor	= 24
    ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

# export
pdf(
  file.path(sp_out_dir, "spatial plot annotation.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_anno %>% walk(print)
dev.off()

### CN spatial plot -----------------
cn_sp_out_dir <- file.path(CN_out_dir, 'CN spatial plot') 
dir.create(cn_sp_out_dir)

# cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

# cn10 
all_data_meta_cn10 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn10, file.path(CN_out_dir, 'all_data_meta_cn10.rds') )
# all_data_meta_cn10 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn10.rds'))

cn10_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


# cn15 
all_data_meta_cn15 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn15, file.path(CN_out_dir, 'all_data_meta_cn15.rds') )
# all_data_meta_cn15 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn15.rds'))

cn15_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn15, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn15_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN15 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

# ### CN2 spatial plot -----------------
# cn_sp_out_dir <- file.path(CN_out_dir, 'cn spatial plot') 
# dir.create(cn_sp_out_dir)
# 
# # cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# # image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')
# 
# all_data_meta <- 
#   all_data[[]] %>%
#   rownames_to_column('cell_id') %>%
#   left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
#   mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
# 
# saveRDS(all_data_meta, file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# 
# cn_sp_p_ls <-
#   cells_position_ls %>% 
#   imap(function(cells_position, project){ # browser()
#     
#     # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
#     # group_info  <- anno_data$group[anno_data$project == project]
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = y, y = x)) +
#       geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#       scale_y_reverse() +
#       coord_equal() +
#       ggsci::scale_fill_igv() +
#       theme_void() +
#       labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
#     p
#   })
# 
# # export image
# cn_sp_p_ls %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir, str_glue("CN spatial plot {name}.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })


### CN spatial plot2 -----------------
# cn_sp_out_dir2 <- file.path(CN_out_dir, 'cn spatial plot2') 
# dir.create(cn_sp_out_dir2)
# 
# all_data_meta <- readRDS( file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# cn_sp_p_ls2 <-
#   cells_position_ls %>% # .[1] %>%
#   imap(function(cells_position, project){ # browser()
#     cn_color <- structure(
#       names = paste0("CN", 1:10),
#       ggsci::pal_igv()(10)
#     ) %>%
#       c(
#         Others = "gray90"
#       )
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p_ls <-
#       paste0("CN", 1:10) %>% # .[1] %>%
#       set_names(., .) %>%
#       map(function(cn){
#         dat2 <- dat
#         
#         dat2$CN2 <- as.character(dat2$CN2)
#         dat2$CN2[dat2$CN2 != cn] <- "Others"
#         
#         p <- dat2 %>%
#           ggplot(aes(x = y, y = x)) +
#           geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#           scale_y_reverse() +
#           coord_equal() +
#           # ggsci::scale_fill_igv() +
#           scale_fill_manual(values = cn_color) +
#           theme_void() +
#           labs(fill = NULL) 
#         p
#       })
#     
#     p_ls
#   })

# # export image
# cn_sp_p_ls2 %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir2, str_glue("CN spatial plot {name} - split by CN.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })

### cn10 cell props ---------
all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn10_props <- 
  all_cn10_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

# all_cn_count_all_zeros_df <-
#   expand.grid(
#     CN2     = all_cn_count_df$CN2 %>% unique(),
#     project = all_cn_count_df$project %>% unique(),
#     group   = all_cn_count_df$group %>% unique()
#   )

p_cn10_props_group <- 
  all_cn10_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(method = "t.test",label = "p.format", size=3)+ 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn10_props %>% print
p_cn10_props_group %>% print
dev.off()


### cn15 cell props ---------
all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn15_props <- 
  all_cn15_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

p_cn15_props_group <- 
  all_cn15_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(label = "p.format", size=3) + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN15_boxplot_and_barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn15_props %>% print
p_cn15_props_group %>% print
dev.off()













## 修改图片CN10 --------
all_data_meta_cn10 <- readRDS('spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/all_data_meta_cn10.rds')

all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn10_count_df %>% head()

write_csv(all_cn10_count_df, file = "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/all_cn10_count_df.csv", col_names = T, )

# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn10_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/shapiro_test_results_CN10.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn10_count_df, CN == "CN1")
print(plot_data_CN1)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.06, 0.14),
    breaks = seq(0.06, 0.14, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn10_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.12),
    breaks = seq(0, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn10_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.2),
    breaks = seq(0, 0.2, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn10_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.1, 0.2),
    breaks = seq(0.1, 0.2, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN4.pdf", p_detailed, width = 10, height = 8)





plot_data_CN5  <- subset(all_cn10_count_df, CN == "CN5")
print(plot_data_CN5)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.06, 0.2),
    breaks = seq(0.06, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN5",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN5.pdf", p_detailed, width = 10, height = 8)





plot_data_CN6  <- subset(all_cn10_count_df, CN == "CN6")
print(plot_data_CN6)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.1, 0.2),
    breaks = seq(0.1, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN6",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn10_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.18),
    breaks = seq(0.05, 0.18, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn10_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.1, 0.2),
    breaks = seq(0.1, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn10_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.07, 0.15),
    breaks = seq(0.07, 0.15, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn10_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.08),
    breaks = seq(0, 0.08, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN10.pdf", p_detailed, width = 10, height = 8)








## 修改图片CN15 --------
all_data_meta_cn15 <- readRDS('spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/all_data_meta_cn15.rds')

all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn15_count_df %>% head()


# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn15_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/shapiro_test_results_CN15.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn15_count_df, CN == "CN1")
print(plot_data_CN1)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.1),
    breaks = seq(0.05, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn15_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn15_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn15_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.09, 0.18),
    breaks = seq(0.09, 0.18, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN4.pdf", p_detailed, width = 10, height = 8)





# plot_data_CN5  <- subset(all_cn15_count_df, CN == "CN5")
# print(plot_data_CN5)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0, 0.1),
#     breaks = seq(0, 0.1, by = 0.02)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN5",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN5.pdf", p_detailed, width = 10, height = 8)
# 
# 
# 
# 
# 
# plot_data_CN6  <- subset(all_cn15_count_df, CN == "CN6")
# print(plot_data_CN6)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0.04, 0.12),
#     breaks = seq(0.04, 0.12, by = 0.02)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN6",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn15_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn15_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6,
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15,
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.07, 0.09),
    breaks = seq(0.07, 0.09, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn15_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.06),
    breaks = seq(0.01, 0.06, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn15_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.07),
    breaks = seq(0.01, 0.07, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN10.pdf", p_detailed, width = 10, height = 8)






plot_data_CN11  <- subset(all_cn15_count_df, CN == "CN11")
print(plot_data_CN11)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN11, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.1),
    breaks = seq(0.04, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN11",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN11.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN11.pdf", p_detailed, width = 10, height = 8)






plot_data_CN12  <- subset(all_cn15_count_df, CN == "CN12")
print(plot_data_CN12)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN12, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.14),
    breaks = seq(0.04, 0.14, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN12",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN12.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN12.pdf", p_detailed, width = 10, height = 8)






plot_data_CN13  <- subset(all_cn15_count_df, CN == "CN13")
print(plot_data_CN13)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN13, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.08),
    breaks = seq(0.01, 0.08, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN13",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN13.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN13.pdf", p_detailed, width = 10, height = 8)






plot_data_CN14  <- subset(all_cn15_count_df, CN == "CN14")
print(plot_data_CN14)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN14, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN14",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN14.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN14.pdf", p_detailed, width = 10, height = 8)






plot_data_CN15  <- subset(all_cn15_count_df, CN == "CN15")
print(plot_data_CN15)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN15, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6,
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15,
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.08),
    breaks = seq(0.04, 0.08, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN15",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN15.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN15.pdf", p_detailed, width = 10, height = 8)













# Proximity analysis -------
### 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

All_type_Proximity_out_dir <- file.path(spatial_out_dir, "Proximity analysis")
dir.create(All_type_Proximity_out_dir)

cells_position_ls <- 
  all_data_filt@images %>%
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) # %>% 
  # setNames(., str_remove(names(.), "\\..*$"))
  # setNames(., sample_mapping[names(.)])



# all_data_cn_meta
all_data_meta_pos <- 
  cells_position_ls %>% 
  imap(function(cells_position, slide){# browser()
    all_data_meta <- 
      all_data_filt[[]] %>% 
      dplyr::mutate(annotation_full = annotation)
    proj_meta_sub <- 
      all_data_meta %>% 
      dplyr::select(annotation_full, group, project) %>% 
      rownames_to_column("cell_id") %>% 
      dplyr::filter(project == slide)
    cells_position[, 1:3] %>% 
      inner_join(proj_meta_sub, by = "cell_id") %>% 
      mutate(annotation_full = factor(annotation_full)) %>% 
      dplyr::rename(Classification1 = annotation_full, 
                    X = x, 
                    Y = y, 
                    Slide_ID = project)
  }) %>% 
  enframe() %>% 
  unnest(cols = value) %>% 
  dplyr::select(-name) %>% 
  column_to_rownames("cell_id")

all_data_meta_pos$Slide_ID <- str_c(all_data_meta_pos$group, all_data_meta_pos$Slide_ID, sep = "_")

all_data_meta_pos_subset <- all_data_meta_pos %>% rownames_to_column("cell_id")
# readRDS('data/all_data_meta_pos_subset.rds')

all_data_meta_pos_subset %>% head

# cell_id        X        Y Classification1 group       Slide_ID
# 1 1_aaaaaagg-1 1670.389 185.1358      Macrophage  CD-F X_202413076_32
# 2 1_aaadlghn-1 1681.172 176.4662          Plasma  CD-F X_202413076_32
# 3 1_aaadmplc-1 1677.561 182.6699               T  CD-F X_202413076_32
# 4 1_aaadoogj-1 1677.109 165.6637          Plasma  CD-F X_202413076_32
# 5 1_aaadpljf-1 1671.887 176.8148          Plasma  CD-F X_202413076_32
# 6 1_aaafaiko-1 1719.466 192.8071          Plasma  CD-F X_202413076_32
# ......


### 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

#### 2.1 proximity analysis --------
res_pa <- 
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset,
    query = all_data_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 'auto', # 'auto',500
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 8 / all_data_filt@images$`24M-1`@scale.factors$spot, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
saveRDS(res_pa, file = file.path(All_type_Proximity_out_dir, "res_pa.rds"))
# res_pa <- readRDS(file.path(All_type_Proximity_out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')

# summarized by group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos_subset,
#     str_c(Classification1, group, Slide_ID, sep = '..')
#   )
pa_groups <- 
  all_data_meta_pos_subset %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE) 

pa_groups2 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn10[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

pa_groups3 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn15[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
res_pa_df2 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups2, FUN = "median")})
res_pa_df3 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups3, FUN = "median")})


res_pa_df  <- 
  res_pa_df %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project"), sep = "\\.\\.") 

res_pa_df2  <- 
  res_pa_df2 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.")
res_pa_df2$CN <- factor(res_pa_df2$CN, levels = paste0("CN", 1:10))

res_pa_df3  <- 
  res_pa_df3 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df3$CN <- factor(res_pa_df3$CN, levels = paste0("CN", 1:15))


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_res_pa_ls_cn10 <-
  res_pa_df2 %>%
  with(split(res_pa_df2, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

p_res_pa_ls_cn15 <-
  res_pa_df3 %>%
  with(split(res_pa_df3, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

# p_res_pa_ls <-
#   res_pa_df2 %>%
#   with(split(res_pa_df2, CN2)) %>%
#   imap(function(dat, name){
#     dat %>%
#   # dplyr::filter(Target == "Epithelial cells") %>%
#   ggplot(aes(x = Others, y = cell_props, fill = Group)) +
#   geom_boxplot(width = 0.45, outlier.colour = NA) +
#   # geom_jitter(width = 0.25) +
#   facet_wrap(~ Target, scales = "free_y", ncol = 2)+
#   ggpubr::theme_pubclean() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   ggtitle(name)
#   })

# export plot
pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis.pdf"),
  height = 50,
  width  = 25
)
p_res_pa%>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN10.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn10 %>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN15.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn15 %>% print()
dev.off()










# 修改图片
# res_pa_df

shapiro_results <- res_pa_df %>%
  group_by(Target) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(cell_props)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(cell_props)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, file.path(All_type_Proximity_out_dir, "shapiro_test_results.csv"), row.names = FALSE)







#### 2.2 proximity analysis 2 -------
mpp <- 8 / all_data_filt@images$`24M-1`@scale.factors$spot

res_pa_ls <- list()

res_pa_ls[['10']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 10, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
res_pa_ls[['50']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 50, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['100']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 100, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['150']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 150, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

res_pa_ls[['200']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", #NULL: 提前返回window的细胞近邻的list
    k = 1200, # 'auto',
    searchtype = "radius",
    radius = 200, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

saveRDS(res_pa_ls, file.path(All_type_Proximity_out_dir, "res_pa_ls.rds"))


# epithelial proximity
unique(all_data_meta_pos$Classification1)

proximity_cells <- c(
  "Osteoblasts", 
  "Chondrocytes", 
  "Osteoclasts"
)

proximity_all_data_ls <- 
  proximity_cells %>%
  setNames(., .) %>%
  map(function(cell){
    cell_idx <- all_data_meta_pos$Classification1 == cell
    
    stopifnot(sum(cell_idx) > 0)
    
    # 1 all_data
    all_data_sub  <- all_data_meta_pos[cell_idx, ]
    
    # 2 res_pa
    res_pa_sub_ls <- res_pa_ls %>% map(~ .x[, cell_idx])
    
    # distance diff
    diff_idx_ls <- seq_along(res_pa_sub_ls) %>% rev() %>% { matrix(c(.[-length(.)], .[-1]), ncol=2)} %>% apply(1, c, simplify = FALSE)
    
    
    # diff_idx_ls
    res_pa_sub_diff_ls <- vector(length = 5, mode = "list")
    
    for(idx in diff_idx_ls){
      res_pa_sub_diff_ls[[ idx[1] ]] <- res_pa_sub_ls[[ idx[1] ]] - res_pa_sub_ls[[ idx[2] ]]
    }
    res_pa_sub_diff_ls[[1]] <- res_pa_sub_ls[[1]]
    names(res_pa_sub_diff_ls) <- names(res_pa_sub_ls)
    
    res_pa_sub_count <-
      res_pa_sub_diff_ls %>% 
      map(~ apply(.x, 1, function(x) tapply(x, all_data_sub$Slide_ID, median)))  
    
    res_pa_sub_count <-
      res_pa_sub_count %>% 
      imap(function(x, y){
        x <- x %>% as.data.frame(check.names = FALSE) %>% rownames_to_column('Slide_ID')
        x$distance <- y
        x
      }) %>% 
      bind_rows() 
    
    
    # add Group 
    all_data_group_sub <- all_data_sub[c('group', 'Slide_ID')] %>% distinct()
    
    # return
    list(
      all_data = all_data_group_sub,
      res_pa   = res_pa_sub_count
    )
  })


# proximity_all_data_ls$Fibroblast$res_pa
# # res_pa_epi_count %>%
#   pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
#   left_join(all_data_group_epi, by = "Slide_ID") 


p_proximity_cell_count_distribution_ls <-
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      mutate(distance = as.numeric(distance)) %>%
      
      ggplot(aes(x = Cell, y = count, fill = group)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ distance, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # scale_fill_manual(values = grp_cols_sp) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p 
  })



res_pa_count_tidy_ls <- 
  proximity_all_data_ls %>%
  map(function(dat_ls){ # browser()
    
    all_data_group_sub <- dat_ls$all_data %>% map_dfc(as.character)
    
    res_pa_count_raw <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID") %>%
      mutate(distance = as.numeric(distance)) 
    
    
    res_pa_count_tidy <-
      1:nrow(res_pa_count_raw) %>%  
      map(function(x){
        count <- res_pa_count_raw$count[x] 
        if(count < 1) return(NULL) # <- count是median count，有可能是0.5
        xx <- res_pa_count_raw[x, ] %>% unlist(use.names = F) %>% rep(count) %>% matrix(nrow = count, byrow = T) %>% as.data.frame()
        xx
      }) %>%
      purrr::discard(is.null) %>%
      bind_rows()
    colnames(res_pa_count_tidy) <- colnames(res_pa_count_raw)
    
    res_pa_count_tidy$count <- as.numeric(res_pa_count_tidy$count)
    
    res_pa_count_tidy
  })


# only selecet two macrohpages
p_proximity_cell_count_distribution_ls2 <-
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      # mutate(distance = as.numeric(distance)) %>%
      mutate(distance = factor(distance, levels = c("10", "50", "100", "150", "200"))) %>%
      
      # dplyr::filter(Cell %in% c("C1QC Macrophage", "SELENOP Macrophage")) %>%
      
      ggplot(aes(x = Cell, y = count, fill = distance)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ group, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # scale_fill_manual(values = grp_cols_sp) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p 
  })

# # show
# head(res_pa_epi_count_tidy, 10)

p_proximity_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    # p <-
    #   res_pa_count_tidy %>%
    #   mutate(distance = as.numeric(distance)) %>%
    #   
    #   ggplot(aes(x = distance, fill = Group)) + #, color = Group)) +
    #   geom_density(alpha = 0.6) +
    #   # geom_smooth(aes(y = count), method = "loess") +
    #   facet_wrap(~ Cell, scales = "free_y") +
    #   theme_bw() +
    #   scale_fill_manual(values = grp_cols_sp) +
    #   labs(
    #     x = 'Distance(um)',
    #     y = 'Median count'
    #   )
    
    p <-
      res_pa_count_tidy %>%
      group_by(Cell, distance, group) %>%
      summarise(count = median(count)) %>%
      
      mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
      
      ggplot(aes(x = distance, y = count, color = group)) +
      # geom_density(alpha = 0.6) +
      geom_smooth(method = 'loess') +
      facet_wrap(~ Cell, scales = "free_y") +
      ggpubr::theme_pubr() +
      theme(
        strip.background = element_blank(),
        panel.grid = element_blank()
      ) +
      # scale_color_manual(values = grp_cols_sp) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })




# export pdf
pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - barplot-1.pdf"), width = 18, height = 8)
p_proximity_cell_count_distribution_ls %>% walk(print)
dev.off()

pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - barplot-2.pdf"), width = 18, height = 6)
p_proximity_cell_count_distribution_ls2 %>% walk(print)
dev.off()

# # export pdf
# pdf("proximity analysis - bxolplot - combined.pdf", width = 18, height = 8)
# print(p_proximity_cell_count_distribution)
# dev.off()


# export pdf
pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - smooth_lines.pdf"), width = 14, height = 8)
p_proximity_ls %>% walk(print)
dev.off()


# radar plot 
p_proximity_radar_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    p <-
      res_pa_count_tidy %>%
      
      # mutate(distance = factor(distance, levels = c('10', '50', '100', '150', '200' ))) %>%
      mutate(distance = as.numeric(distance)) %>%
      group_by(group, Cell, distance) %>%
      summarise(count = median(count)) %>%
      ggplot(aes(x = Cell, y = distance, size = count, color = group)) +
      geom_vline(xintercept = 0:12 + 0.5, color = "gray") +
      geom_hline(yintercept = c(10, 50, 100, 150, 200), color = "gray") +
      geom_point(position = position_jitterdodge(
        jitter.width = 0,
        jitter.height = 10,
        dodge.width = 0.75
      )) +
      geom_text(
        data = data.frame(x = 0, y = c(10, 50, 100, 150, 200)),
        aes(x = x, y = y, label = paste0(y, "um")),
        inherit.aes = FALSE,
        vjust = -1.5
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(1,1,1,6, 'mm')
      ) +
      scale_size(range = c(3, 12)) +
      scale_y_continuous(limits = c(-100, 220)) +
      # scale_color_manual(values = grp_cols_sp)+
      coord_polar(theta = 'x', clip = 'off') +
      labs(x = NULL, y = NULL, title = name)
    p 
  })





# export pdf
pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - radarplot.pdf"), width = 14, height = 12)
p_proximity_radar_ls %>% walk(print)
dev.off()
















# <------------
# <------------
# <------------
# <------------
# rctd remove erythrocyte --------
# Date <- Sys.Date()
Date <- "2026-01-28"
spatial_out_dir <- str_glue("spatial_analysis_remove_erythrocyte-{Date}")
dir.create(spatial_out_dir)

all_data <- readRDS("data_for_analysis/all_data_filt.rds")
all_data@meta.data$annotation <- all_data@meta.data$rctd_20251201
unique(all_data$annotation)
# all_data_filt <- subset(all_data, annotation != 'Other')
all_data_filt <- all_data[, !all_data@meta.data$annotation %in% c("Other")]

# 只保留原始图像，移除自动生成的重复项
original_images <- c("24M-1", "24M-2", "4M-1", "4M-2")
images_to_keep <- names(all_data_filt@images)[names(all_data_filt@images) %in% original_images]

# 如果原始图像都存在，直接子集化
if(length(images_to_keep) == 4) {
  all_data_filt@images <- all_data_filt@images[images_to_keep]
  cat("修复完成！当前图像名称:", paste(names(all_data_filt@images), collapse = ", "), "\n")
} else {
  # 如果原始图像不全，使用更安全的方法
  cat("采取备用方案...\n")
  # 方法：基于图像数据内容去重（更安全）
  image_data_list <- list()
  for(img_name in names(all_data_filt@images)) {
    img <- all_data_filt@images[[img_name]]
    # 使用关键标识来去重
    img_key <- paste(img@assay, img@key, sep = "_")
    if(!img_key %in% names(image_data_list)) {
      image_data_list[[img_key]] <- img
    }
  }
  all_data_filt@images <- image_data_list
}

Idents(all_data_filt) <- "annotation"

sample_mapping <- c(
  "X24M.16.008um" = "24M-1",
  "X24M.16.008um.2" = "24M-2",
  "X4M.7.008um" = "4M-1",
  "X4M.7.008um.4" = "4M-2"
)

cells_position_ls <- 
  all_data_filt@images %>% 
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) # %>%
# setNames(., str_remove(names(.), "\\..*$"))
# setNames(., sample_mapping[names(.)])

cells_position_df <- cells_position_ls %>% bind_rows() 






# CN analysis knn=25---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

CN_out_dir <- file.path(spatial_out_dir, "CN analysis - knn25")
dir.create(CN_out_dir, recursive = T)


all_data_meta_for_cn_analysis <- 
  all_data_filt[[c("project", "group", "annotation")]] %>%
  dplyr::mutate(annotation_full = annotation ) %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x', 'y')], by = "cell_id") %>%
  rename(X = x, Y = y) %>%
  mutate(
    annotation_full = factor(annotation_full), # import!
    obj_id = cell_id # import!
  )


### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 25,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))
# res_cn_10 <- readRDS(file.path(CN_out_dir, "res_cn_10.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn_10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn_10_heat$gtable)
dev.off()


#### CN annotaton --------------------
# cn_annotation <- c(
#   "CN1"  = "Immune core",
#   "CN2"  = "Plasma core",
#   "CN3"  = "Tregs immune",
#   "CN4"  = "Mac_Fibro_Ec",
#   "CN5"  = "Myo Stroma",
#   "CN6"  = "Fibro Stroma",
#   "CN7"  = "CD163 Macophage Stroma",
#   "CN8"  = "Enterocytes",
#   "CN9"  = "Myeloid enterocytes",
#   "CN10" = "Amplifying area"
# )

# cn_annotation <- c(
#   "CN1"  = "CN1",
#   "CN2"  = "CN2",
#   "CN3"  = "CN3",
#   "CN4"  = "CN4",
#   "CN5"  = "CN5",
#   "CN6"  = "CN6",
#   "CN7"  = "CN7",
#   "CN8"  = "CN8",
#   "CN9"  = "CN9",
#   "CN10" = "CN10"
# )
# 
# p_cn_10_heat2 <-
#   cn_cluster_FC %>%
#   pheatmap::pheatmap(
#     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#     border_color = NA,
#     treeheight_row = 0,
#     labels_row = cn_annotation,
#     # clustering_method = 'ward.D2',
#     cluster_rows = F,
#     scale = "column",
#     silent = TRUE
#   )
# 
# pdf(
#   file.path(CN_out_dir, "CN_10_heatmap_annotation.pdf"),
#   height = 6,
#   width  = 11
# )
# cowplot::plot_grid(p_cn_10_heat2$gtable)
# dev.off()
# 



### CN k15 -----
res_cn_15 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 25,
  k.cluster      = 15
)

saveRDS(res_cn_15, file.path(CN_out_dir, "res_cn_15.rds"))
# res_cn_15 <- readRDS(file.path(CN_out_dir, "res_cn_15.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn15_cluster_mean <- res_cn_15$cluster_centers
cn15_baseMean     <- res_cn_15$base_mean

# log2(FC)
cn15_cluster_FC <- log2(1 + sweep(cn15_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn15_baseMean))

# heatmap
p_cn15_heat <-
  pheatmap::pheatmap(
    cn15_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn15_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_15_heatmap.pdf"),
  height = 8,
  width  = 10
)
cowplot::plot_grid(p_cn15_heat$gtable)
dev.off()




### CN k20 -----
res_cn_20 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 25,
  k.cluster      = 20
)

saveRDS(res_cn_20, file.path(CN_out_dir, "res_cn_20.rds"))
# res_cn_15 <- readRDS(file.path(CN_out_dir, "res_cn_15.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn20_cluster_mean <- res_cn_20$cluster_centers
cn20_baseMean     <- res_cn_20$base_mean

# log2(FC)
cn20_cluster_FC <- log2(1 + sweep(cn20_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn20_baseMean))

# heatmap
p_cn20_heat <-
  pheatmap::pheatmap(
    cn20_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn20_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_20_heatmap.pdf"),
  height = 10,
  width  = 10
)
cowplot::plot_grid(p_cn20_heat$gtable)
dev.off()




### CN k25 -----
res_cn_25 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 25,
  k.cluster      = 25
)

saveRDS(res_cn_25, file.path(CN_out_dir, "res_cn_25.rds"))
# res_cn_25 <- readRDS(file.path(CN_out_dir, "res_cn_25.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn25_cluster_mean <- res_cn_25$cluster_centers
cn25_baseMean     <- res_cn_25$base_mean

# log2(FC)
cn25_cluster_FC <- log2(1 + sweep(cn25_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn25_baseMean))

# heatmap
p_cn25_heat <-
  pheatmap::pheatmap(
    cn25_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn25_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_25_heatmap.pdf"),
  height = 12,
  width  = 10
)
cowplot::plot_grid(p_cn25_heat$gtable)
dev.off()





### spatial plot -----------------
sp_out_dir <- file.path(CN_out_dir, 'spatial plot') 
dir.create(sp_out_dir)


p_sp_plot_anno <-
  # Images(all_data) %>%
  # setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  # imap(function(img, img_name){
  sample_mapping %>% 
  imap(function(img_name, img){
    # print(img_name)
    # print(img)
    # return()
    p <- SpatialDimPlot(
      all_data_filt,
      images = img_name, 
      group.by = "annotation",
      pt.size.factor	= 24
    ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

# export
pdf(
  file.path(sp_out_dir, "spatial plot annotation.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_anno %>% walk(print)
dev.off()

### CN spatial plot -----------------
cn_sp_out_dir <- file.path(CN_out_dir, 'CN spatial plot') 
dir.create(cn_sp_out_dir)

# cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

# cn10 
all_data_meta_cn10 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn10, file.path(CN_out_dir, 'all_data_meta_cn10.rds') )
# all_data_meta_cn10 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn10.rds'))

cn10_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


# cn15 
all_data_meta_cn15 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn15, file.path(CN_out_dir, 'all_data_meta_cn15.rds') )
# all_data_meta_cn15 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn15.rds'))

cn15_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn15, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn15_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN15 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


# cn20 
all_data_meta_cn20 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_20$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn20, file.path(CN_out_dir, 'all_data_meta_cn20.rds') )
# all_data_meta_cn20 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn20.rds'))

cn20_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn20, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn20_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN20 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


# cn25 
all_data_meta_cn25 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_25$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn25, file.path(CN_out_dir, 'all_data_meta_cn25.rds') )
# all_data_meta_cn25 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn25.rds'))

cn25_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn25, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn25_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN25 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })



# ### CN2 spatial plot -----------------
# cn_sp_out_dir <- file.path(CN_out_dir, 'cn spatial plot') 
# dir.create(cn_sp_out_dir)
# 
# # cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# # image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')
# 
# all_data_meta <- 
#   all_data[[]] %>%
#   rownames_to_column('cell_id') %>%
#   left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
#   mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
# 
# saveRDS(all_data_meta, file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# 
# cn_sp_p_ls <-
#   cells_position_ls %>% 
#   imap(function(cells_position, project){ # browser()
#     
#     # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
#     # group_info  <- anno_data$group[anno_data$project == project]
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = y, y = x)) +
#       geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#       scale_y_reverse() +
#       coord_equal() +
#       ggsci::scale_fill_igv() +
#       theme_void() +
#       labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
#     p
#   })
# 
# # export image
# cn_sp_p_ls %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir, str_glue("CN spatial plot {name}.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })


### CN spatial plot2 -----------------
# cn_sp_out_dir2 <- file.path(CN_out_dir, 'cn spatial plot2') 
# dir.create(cn_sp_out_dir2)
# 
# all_data_meta <- readRDS( file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# cn_sp_p_ls2 <-
#   cells_position_ls %>% # .[1] %>%
#   imap(function(cells_position, project){ # browser()
#     cn_color <- structure(
#       names = paste0("CN", 1:10),
#       ggsci::pal_igv()(10)
#     ) %>%
#       c(
#         Others = "gray90"
#       )
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p_ls <-
#       paste0("CN", 1:10) %>% # .[1] %>%
#       set_names(., .) %>%
#       map(function(cn){
#         dat2 <- dat
#         
#         dat2$CN2 <- as.character(dat2$CN2)
#         dat2$CN2[dat2$CN2 != cn] <- "Others"
#         
#         p <- dat2 %>%
#           ggplot(aes(x = y, y = x)) +
#           geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#           scale_y_reverse() +
#           coord_equal() +
#           # ggsci::scale_fill_igv() +
#           scale_fill_manual(values = cn_color) +
#           theme_void() +
#           labs(fill = NULL) 
#         p
#       })
#     
#     p_ls
#   })

# # export image
# cn_sp_p_ls2 %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir2, str_glue("CN spatial plot {name} - split by CN.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })

### cn10 cell props ---------
all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn10_props <- 
  all_cn10_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

# all_cn_count_all_zeros_df <-
#   expand.grid(
#     CN2     = all_cn_count_df$CN2 %>% unique(),
#     project = all_cn_count_df$project %>% unique(),
#     group   = all_cn_count_df$group %>% unique()
#   )

p_cn10_props_group <- 
  all_cn10_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(method = "t.test",label = "p.format", size=3)+ 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10 barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn10_props %>% print
dev.off()

pdf(
  file.path(CN_out_dir, "CN10 boxplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn10_props_group %>% print
dev.off()


### cn15 cell props ---------
all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn15_props <- 
  all_cn15_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

p_cn15_props_group <- 
  all_cn15_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(label = "p.format", size=3) + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN15 barplot.pdf"),
  height = 8, 
  width  = 10
)
p_cn15_props %>% print
dev.off()

pdf(
  file.path(CN_out_dir, "CN15 boxplot.pdf"),
  height = 6, 
  width  = 12
)
p_cn15_props_group %>% print
dev.off()



### cn20 cell props ---------
all_cn20_count_df <- 
  all_data_meta_cn20 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn20_props <- 
  all_cn20_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

p_cn20_props_group <- 
  all_cn20_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(label = "p.format", size=3) + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN20 barplot.pdf"),
  height = 10, 
  width  = 10
)
p_cn20_props %>% print
dev.off()

pdf(
  file.path(CN_out_dir, "CN20 boxplot.pdf"),
  height = 6, 
  width  = 14
)
p_cn20_props_group %>% print
dev.off()



### cn25 cell props ---------
all_cn25_count_df <- 
  all_data_meta_cn25 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn25_props <- 
  all_cn25_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

p_cn25_props_group <- 
  all_cn25_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(label = "p.format", size=3) + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN25 barplot.pdf"),
  height = 12, 
  width  = 10
)
p_cn25_props %>% print
dev.off()

pdf(
  file.path(CN_out_dir, "CN25 boxplot.pdf"),
  height = 6, 
  width  = 16
)
p_cn25_props_group %>% print
dev.off()












## 修改图片CN10 --------
all_data_meta_cn10 <- readRDS('spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/all_data_meta_cn10.rds')

all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn10_count_df %>% head()

write_csv(all_cn10_count_df, file = "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/all_cn10_count_df.csv", col_names = T, )

# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn10_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/shapiro_test_results_CN10.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn10_count_df, CN == "CN1")
print(plot_data_CN1)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.15),
    breaks = seq(0, 0.15, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn10_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.2),
    breaks = seq(0.01, 0.2, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn10_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.15),
    breaks = seq(0.05, 0.15, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn10_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.12, 0.25),
    breaks = seq(0.12, 0.25, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN4.pdf", p_detailed, width = 10, height = 8)





plot_data_CN5  <- subset(all_cn10_count_df, CN == "CN5")
print(plot_data_CN5)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.1),
    breaks = seq(0.05, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN5",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN5.pdf", p_detailed, width = 10, height = 8)





plot_data_CN6  <- subset(all_cn10_count_df, CN == "CN6")
print(plot_data_CN6)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.1, 0.2),
    breaks = seq(0.1, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN6",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn10_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.15),
    breaks = seq(0.04, 0.15, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn10_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.12, 0.15),
    breaks = seq(0.12, 0.15, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn10_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.15),
    breaks = seq(0.05, 0.15, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn10_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.08),
    breaks = seq(0, 0.08, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN10-CN10.pdf", p_detailed, width = 10, height = 8)








## 修改图片CN15 --------
all_data_meta_cn15 <- readRDS('spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/all_data_meta_cn15.rds')

all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn15_count_df %>% head()


# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn15_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/shapiro_test_results_CN15.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn15_count_df, CN == "CN1")
print(plot_data_CN1)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.1),
    breaks = seq(0.05, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn15_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.12),
    breaks = seq(0, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn15_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn15_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.08, 0.15),
    breaks = seq(0.08, 0.15, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN4.pdf", p_detailed, width = 10, height = 8)





# plot_data_CN5  <- subset(all_cn15_count_df, CN == "CN5")
# print(plot_data_CN5)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0, 0.1),
#     breaks = seq(0, 0.1, by = 0.02)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN5",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN5.pdf", p_detailed, width = 10, height = 8)
# 
# 
# 
# 
# 
# plot_data_CN6  <- subset(all_cn15_count_df, CN == "CN6")
# print(plot_data_CN6)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0.04, 0.12),
#     breaks = seq(0.04, 0.12, by = 0.02)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN6",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn15_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn15_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6,
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15,
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.06, 0.1),
    breaks = seq(0.06, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn15_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.06),
    breaks = seq(0.01, 0.06, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn15_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.07),
    breaks = seq(0, 0.07, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN10.pdf", p_detailed, width = 10, height = 8)






plot_data_CN11  <- subset(all_cn15_count_df, CN == "CN11")
print(plot_data_CN11)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN11, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.1),
    breaks = seq(0.04, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN11",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN11.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN11.pdf", p_detailed, width = 10, height = 8)






plot_data_CN12  <- subset(all_cn15_count_df, CN == "CN12")
print(plot_data_CN12)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN12, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.03, 0.14),
    breaks = seq(0.03, 0.14, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN12",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN12.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN12.pdf", p_detailed, width = 10, height = 8)






plot_data_CN13  <- subset(all_cn15_count_df, CN == "CN13")
print(plot_data_CN13)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN13, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.08),
    breaks = seq(0, 0.08, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN13",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN13.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN13.pdf", p_detailed, width = 10, height = 8)






plot_data_CN14  <- subset(all_cn15_count_df, CN == "CN14")
print(plot_data_CN14)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN14, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN14",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN14.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN14.pdf", p_detailed, width = 10, height = 8)






plot_data_CN15  <- subset(all_cn15_count_df, CN == "CN15")
print(plot_data_CN15)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN15, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6,
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15,
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.08),
    breaks = seq(0.04, 0.08, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN15",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN15.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn25/CN15-CN15.pdf", p_detailed, width = 10, height = 8)










# CN analysis knn=50---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

CN_out_dir <- file.path(spatial_out_dir, "CN analysis - knn50")
dir.create(CN_out_dir, recursive = T)


all_data_meta_for_cn_analysis <- 
  all_data_filt[[c("project", "group", "annotation")]] %>%
  dplyr::mutate(annotation_full = annotation ) %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x', 'y')], by = "cell_id") %>%
  rename(X = x, Y = y) %>%
  mutate(
    annotation_full = factor(annotation_full), # import!
    obj_id = cell_id # import!
  )


### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))
# res_cn_10 <- readRDS(file.path(CN_out_dir, "res_cn_10.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn_10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn_10_heat$gtable)
dev.off()


#### CN annotaton --------------------
# cn_annotation <- c(
#   "CN1"  = "Immune core",
#   "CN2"  = "Plasma core",
#   "CN3"  = "Tregs immune",
#   "CN4"  = "Mac_Fibro_Ec",
#   "CN5"  = "Myo Stroma",
#   "CN6"  = "Fibro Stroma",
#   "CN7"  = "CD163 Macophage Stroma",
#   "CN8"  = "Enterocytes",
#   "CN9"  = "Myeloid enterocytes",
#   "CN10" = "Amplifying area"
# )

# cn_annotation <- c(
#   "CN1"  = "CN1",
#   "CN2"  = "CN2",
#   "CN3"  = "CN3",
#   "CN4"  = "CN4",
#   "CN5"  = "CN5",
#   "CN6"  = "CN6",
#   "CN7"  = "CN7",
#   "CN8"  = "CN8",
#   "CN9"  = "CN9",
#   "CN10" = "CN10"
# )
# 
# p_cn_10_heat2 <-
#   cn_cluster_FC %>%
#   pheatmap::pheatmap(
#     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#     border_color = NA,
#     treeheight_row = 0,
#     labels_row = cn_annotation,
#     # clustering_method = 'ward.D2',
#     cluster_rows = F,
#     scale = "column",
#     silent = TRUE
#   )
# 
# pdf(
#   file.path(CN_out_dir, "CN_10_heatmap_annotation.pdf"),
#   height = 6,
#   width  = 11
# )
# cowplot::plot_grid(p_cn_10_heat2$gtable)
# dev.off()
# 



### CN k15 -----
res_cn_15 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 15
)

saveRDS(res_cn_15, file.path(CN_out_dir, "res_cn_15.rds"))
# res_cn_15 <- readRDS(file.path(CN_out_dir, "res_cn_15.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn15_cluster_mean <- res_cn_15$cluster_centers
cn15_baseMean     <- res_cn_15$base_mean

# log2(FC)
cn15_cluster_FC <- log2(1 + sweep(cn15_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn15_baseMean))

# heatmap
p_cn15_heat <-
  pheatmap::pheatmap(
    cn15_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn15_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_15_heatmap.pdf"),
  height = 8,
  width  = 10
)
cowplot::plot_grid(p_cn15_heat$gtable)
dev.off()




### CN k20 -----
res_cn_20 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 20
)

saveRDS(res_cn_20, file.path(CN_out_dir, "res_cn_20.rds"))
# res_cn_15 <- readRDS(file.path(CN_out_dir, "res_cn_15.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn20_cluster_mean <- res_cn_20$cluster_centers
cn20_baseMean     <- res_cn_20$base_mean

# log2(FC)
cn20_cluster_FC <- log2(1 + sweep(cn20_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn20_baseMean))

# heatmap
p_cn20_heat <-
  pheatmap::pheatmap(
    cn20_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn20_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_20_heatmap.pdf"),
  height = 10,
  width  = 10
)
cowplot::plot_grid(p_cn20_heat$gtable)
dev.off()




### CN k25 -----
res_cn_25 <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 25
)

saveRDS(res_cn_25, file.path(CN_out_dir, "res_cn_25.rds"))
# res_cn_25 <- readRDS(file.path(CN_out_dir, "res_cn_25.rds"))


#### heatmap ----------------------
# FoldChange cluster/allCell
cn25_cluster_mean <- res_cn_25$cluster_centers
cn25_baseMean     <- res_cn_25$base_mean

# log2(FC)
cn25_cluster_FC <- log2(1 + sweep(cn25_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn25_baseMean))

# heatmap
p_cn25_heat <-
  pheatmap::pheatmap(
    cn25_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn25_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_25_heatmap.pdf"),
  height = 12,
  width  = 10
)
cowplot::plot_grid(p_cn25_heat$gtable)
dev.off()





### spatial plot -----------------
sp_out_dir <- file.path(CN_out_dir, 'spatial plot') 
dir.create(sp_out_dir)


p_sp_plot_anno <-
  # Images(all_data) %>%
  # setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  # imap(function(img, img_name){
  sample_mapping %>% 
  imap(function(img_name, img){
    # print(img_name)
    # print(img)
    # return()
    p <- SpatialDimPlot(
      all_data_filt,
      images = img_name, 
      group.by = "annotation",
      pt.size.factor	= 24
    ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

# export
pdf(
  file.path(sp_out_dir, "spatial plot annotation.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_anno %>% walk(print)
dev.off()


### CN spatial plot -----------------
cn_sp_out_dir <- file.path(CN_out_dir, 'CN spatial plot') 
dir.create(cn_sp_out_dir)

# cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

# cn10 
all_data_meta_cn10 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn10, file.path(CN_out_dir, 'all_data_meta_cn10.rds') )
# all_data_meta_cn10 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn10.rds'))

cn10_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


# cn15 
all_data_meta_cn15 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn15, file.path(CN_out_dir, 'all_data_meta_cn15.rds') )
# all_data_meta_cn15 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn15.rds'))

cn15_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn15, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn15_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN15 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


# cn20 
all_data_meta_cn20 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_20$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn20, file.path(CN_out_dir, 'all_data_meta_cn20.rds') )
# all_data_meta_cn20 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn20.rds'))

cn20_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn20, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn20_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN20 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


# cn25 
all_data_meta_cn25 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn_25$cn_cluster, "cell_id"), by = "cell_id") # %>%
# mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
saveRDS(all_data_meta_cn25, file.path(CN_out_dir, 'all_data_meta_cn25.rds') )
# all_data_meta_cn25 <- readRDS(file.path(CN_out_dir, 'all_data_meta_cn25.rds'))

cn25_sp_p_ls <- 
  cells_position_ls %>% 
  imap(function(cells_position, project){ # browser()
    # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    # group_info  <- anno_data$group[anno_data$project == project]
    dat <-
      cells_position %>%
      right_join(all_data_meta_cn25, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = y, y = x)) +
      geom_tile(aes(fill = CN), height = 10, width = 10)  +
      scale_y_reverse() +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
    p
  })

# export image
cn25_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN25 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })



# ### CN2 spatial plot -----------------
# cn_sp_out_dir <- file.path(CN_out_dir, 'cn spatial plot') 
# dir.create(cn_sp_out_dir)
# 
# # cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# # image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')
# 
# all_data_meta <- 
#   all_data[[]] %>%
#   rownames_to_column('cell_id') %>%
#   left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
#   mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
# 
# saveRDS(all_data_meta, file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# 
# cn_sp_p_ls <-
#   cells_position_ls %>% 
#   imap(function(cells_position, project){ # browser()
#     
#     # anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
#     # group_info  <- anno_data$group[anno_data$project == project]
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = y, y = x)) +
#       geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#       scale_y_reverse() +
#       coord_equal() +
#       ggsci::scale_fill_igv() +
#       theme_void() +
#       labs(fill = NULL) #, title = paste(project, group_info, sep = "_") )
#     p
#   })
# 
# # export image
# cn_sp_p_ls %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir, str_glue("CN spatial plot {name}.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })


### CN spatial plot2 -----------------
# cn_sp_out_dir2 <- file.path(CN_out_dir, 'cn spatial plot2') 
# dir.create(cn_sp_out_dir2)
# 
# all_data_meta <- readRDS( file.path(CN_out_dir, 'all_data_meta_cn.rds') )
# cn_sp_p_ls2 <-
#   cells_position_ls %>% # .[1] %>%
#   imap(function(cells_position, project){ # browser()
#     cn_color <- structure(
#       names = paste0("CN", 1:10),
#       ggsci::pal_igv()(10)
#     ) %>%
#       c(
#         Others = "gray90"
#       )
#     
#     dat <-
#       cells_position %>%
#       right_join(all_data_meta, by = "cell_id")
#     
#     p_ls <-
#       paste0("CN", 1:10) %>% # .[1] %>%
#       set_names(., .) %>%
#       map(function(cn){
#         dat2 <- dat
#         
#         dat2$CN2 <- as.character(dat2$CN2)
#         dat2$CN2[dat2$CN2 != cn] <- "Others"
#         
#         p <- dat2 %>%
#           ggplot(aes(x = y, y = x)) +
#           geom_tile(aes(fill = CN2), height = 10, width = 10)  +
#           scale_y_reverse() +
#           coord_equal() +
#           # ggsci::scale_fill_igv() +
#           scale_fill_manual(values = cn_color) +
#           theme_void() +
#           labs(fill = NULL) 
#         p
#       })
#     
#     p_ls
#   })

# # export image
# cn_sp_p_ls2 %>% 
#   iwalk(function(p, name){
#     pdf(
#       file.path(cn_sp_out_dir2, str_glue("CN spatial plot {name} - split by CN.pdf")),
#       width  = 20,
#       height = 20
#     )
#     print(p)
#     dev.off()
#   })

### cn10 cell props ---------
all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn10_props <- 
  all_cn10_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

# all_cn_count_all_zeros_df <-
#   expand.grid(
#     CN2     = all_cn_count_df$CN2 %>% unique(),
#     project = all_cn_count_df$project %>% unique(),
#     group   = all_cn_count_df$group %>% unique()
#   )

p_cn10_props_group <- 
  all_cn10_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(method = "t.test",label = "p.format", size=3)+ 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10 barplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn10_props %>% print
dev.off()

pdf(
  file.path(CN_out_dir, "CN10 boxplot.pdf"),
  height = 6, 
  width  = 10
)
p_cn10_props_group %>% print
dev.off()


### cn15 cell props ---------
all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn15_props <- 
  all_cn15_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

p_cn15_props_group <- 
  all_cn15_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(label = "p.format", size=3) + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN15 barplot.pdf"),
  height = 8, 
  width  = 10
)
p_cn15_props %>% print
dev.off()

pdf(
  file.path(CN_out_dir, "CN15 boxplot.pdf"),
  height = 6, 
  width  = 12
)
p_cn15_props_group %>% print
dev.off()



### cn20 cell props ---------
all_cn20_count_df <- 
  all_data_meta_cn20 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn20_props <- 
  all_cn20_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

p_cn20_props_group <- 
  all_cn20_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(label = "p.format", size=3) + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN20 barplot.pdf"),
  height = 10, 
  width  = 10
)
p_cn20_props %>% print
dev.off()

pdf(
  file.path(CN_out_dir, "CN20 boxplot.pdf"),
  height = 6, 
  width  = 14
)
p_cn20_props_group %>% print
dev.off()



### cn25 cell props ---------
all_cn25_count_df <- 
  all_data_meta_cn25 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

p_cn25_props <- 
  all_cn25_count_df %>% 
  ggplot(aes(x = project, y = freq, fill = CN)) +
  geom_col() + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

p_cn25_props_group <- 
  all_cn25_count_df %>% 
  # right_join(all_cn_count_all_zeros_df, by = c('CN2', 'project', 'group')) %>% 
  # mutate(freq = freq %>% str_replace_na("0") %>% as.numeric()) %>% 
  ggplot(aes(x = CN, y = freq, fill = group)) + 
  geom_boxplot() + 
  stat_compare_means(label = "p.format", size=3, aes(group = group))+
  # stat_compare_means(label = "p.format", size=3) + 
  # ggsci::scale_fill_aaas() + 
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) + 
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN25 barplot.pdf"),
  height = 12, 
  width  = 10
)
p_cn25_props %>% print
dev.off()

pdf(
  file.path(CN_out_dir, "CN25 boxplot.pdf"),
  height = 6, 
  width  = 16
)
p_cn25_props_group %>% print
dev.off()











## 修改图片CN10 --------
all_data_meta_cn10 <- readRDS('spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/all_data_meta_cn10.rds')

all_cn10_count_df <- 
  all_data_meta_cn10 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn10_count_df %>% head()

write_csv(all_cn10_count_df, file = "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/all_cn10_count_df.csv", col_names = T, )

# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn10_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/shapiro_test_results_CN10.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn10_count_df, CN == "CN1")
print(plot_data_CN1)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.06, 0.14),
    breaks = seq(0.06, 0.14, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn10_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.12),
    breaks = seq(0, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn10_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.2),
    breaks = seq(0, 0.2, by = 0.05)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn10_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.1, 0.2),
    breaks = seq(0.1, 0.2, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN4.pdf", p_detailed, width = 10, height = 8)





plot_data_CN5  <- subset(all_cn10_count_df, CN == "CN5")
print(plot_data_CN5)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.06, 0.2),
    breaks = seq(0.06, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN5",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN5.pdf", p_detailed, width = 10, height = 8)





plot_data_CN6  <- subset(all_cn10_count_df, CN == "CN6")
print(plot_data_CN6)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.1, 0.2),
    breaks = seq(0.1, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN6",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn10_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.18),
    breaks = seq(0.05, 0.18, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn10_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.1, 0.2),
    breaks = seq(0.1, 0.2, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn10_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.07, 0.15),
    breaks = seq(0.07, 0.15, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn10_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.08),
    breaks = seq(0, 0.08, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN10-CN10.pdf", p_detailed, width = 10, height = 8)








## 修改图片CN15 --------
all_data_meta_cn15 <- readRDS('spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/all_data_meta_cn15.rds')

all_cn15_count_df <- 
  all_data_meta_cn15 %>% 
  group_by(project, group, CN) %>% 
  summarise(count = n()) %>% 
  mutate(freq = count/sum(count))

all_cn15_count_df %>% head()


# 对每种细胞类型进行正态性检验
shapiro_results <- all_cn15_count_df %>%
  group_by(CN) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(freq)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(freq)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, "spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/shapiro_test_results_CN15.csv", row.names = FALSE)







plot_data_CN1  <- subset(all_cn15_count_df, CN == "CN1")
print(plot_data_CN1)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN1, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.1),
    breaks = seq(0.05, 0.1, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN1",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN1.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN1.pdf", p_detailed, width = 10, height = 8)





plot_data_CN2  <- subset(all_cn15_count_df, CN == "CN2")
print(plot_data_CN2)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN2, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN2",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN2.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN2.pdf", p_detailed, width = 10, height = 8)





plot_data_CN3  <- subset(all_cn15_count_df, CN == "CN3")
print(plot_data_CN3)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN3, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN3",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN3.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN3.pdf", p_detailed, width = 10, height = 8)





plot_data_CN4  <- subset(all_cn15_count_df, CN == "CN4")
print(plot_data_CN4)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN4, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.09, 0.18),
    breaks = seq(0.09, 0.18, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN4",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN4.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN4.pdf", p_detailed, width = 10, height = 8)





# plot_data_CN5  <- subset(all_cn15_count_df, CN == "CN5")
# print(plot_data_CN5)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN5, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0, 0.1),
#     breaks = seq(0, 0.1, by = 0.02)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN5",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN5.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN5.pdf", p_detailed, width = 10, height = 8)
# 
# 
# 
# 
# 
# plot_data_CN6  <- subset(all_cn15_count_df, CN == "CN6")
# print(plot_data_CN6)
# # 更详细的版本，包含更多自定义选项
# p_detailed <- ggplot(plot_data_CN6, aes(x = group, y = freq, fill = group)) +
#   # 小提琴图（带透明度）
#   geom_violin(
#     alpha = 0.6, 
#     width = 0.9,
#     trim = TRUE,
#     scale = "width"
#   ) +
#   # 箱线图（调整位置和样式）
#   geom_boxplot(
#     width = 0.15, 
#     alpha = 0.9,
#     outlier.color = "red",
#     outlier.shape = 16,
#     outlier.size = 2,
#     position = position_dodge(0.9)
#   ) +
#   # 添加均值点
#   stat_summary(
#     fun = mean, 
#     geom = "point", 
#     shape = 18, 
#     size = 3, 
#     color = "white"
#   ) +
#   # 组间统计比较
#   stat_compare_means(
#     method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
#     label = "p.format",
#     tip.length = 0.02,
#     bracket.size = 0.6,
#     step.increase = 0.1
#   ) +
#   # 颜色方案
#   scale_fill_manual(
#     values = c(
#       "24M" = "#fdc4b6",    # 青色
#       "4M" = "#e8a0b8"      # 深蓝色
#     )
#   ) +
#   # y轴范围（根据您的数据调整）
#   scale_y_continuous(
#     limits = c(0.04, 0.12),
#     breaks = seq(0.04, 0.12, by = 0.02)
#   ) +
#   # 标签和标题
#   labs(
#     title = "CN6",
#     x = "",
#     y = "Expression Level",
#     fill = "group"
#   ) +
#   # 主题美化
#   theme_classic() +
#   theme(
#     plot.title = element_text(
#       hjust = 0.5, 
#       size = 18, 
#       face = "bold",
#       margin = margin(b = 15)
#     ),
#     axis.title.y = element_text(
#       size = 14,
#       margin = margin(r = 10)
#     ),
#     axis.text.x = element_text(
#       size = 12,
#       angle = 0,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(size = 11),
#     legend.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 11),
#     panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(20, 20, 20, 20)
#   )
# 
# # 显示图形
# print(p_detailed)
# 
# # 保存图形
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN6.png", p_detailed, width = 10, height = 8, dpi = 300)
# ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN6.pdf", p_detailed, width = 10, height = 8)





plot_data_CN7  <- subset(all_cn15_count_df, CN == "CN7")
print(plot_data_CN7)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN7, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0, 0.1),
    breaks = seq(0, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN7",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN7.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN7.pdf", p_detailed, width = 10, height = 8)





plot_data_CN8  <- subset(all_cn15_count_df, CN == "CN8")
print(plot_data_CN8)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN8, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6,
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15,
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.07, 0.09),
    breaks = seq(0.07, 0.09, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN8",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN8.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN8.pdf", p_detailed, width = 10, height = 8)





plot_data_CN9  <- subset(all_cn15_count_df, CN == "CN9")
print(plot_data_CN9)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN9, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.06),
    breaks = seq(0.01, 0.06, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN9",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN9.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN9.pdf", p_detailed, width = 10, height = 8)





plot_data_CN10  <- subset(all_cn15_count_df, CN == "CN10")
print(plot_data_CN10)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN10, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.07),
    breaks = seq(0.01, 0.07, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN10",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN10.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN10.pdf", p_detailed, width = 10, height = 8)






plot_data_CN11  <- subset(all_cn15_count_df, CN == "CN11")
print(plot_data_CN11)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN11, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.1),
    breaks = seq(0.04, 0.1, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN11",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN11.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN11.pdf", p_detailed, width = 10, height = 8)






plot_data_CN12  <- subset(all_cn15_count_df, CN == "CN12")
print(plot_data_CN12)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN12, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.14),
    breaks = seq(0.04, 0.14, by = 0.03)
  ) +
  # 标签和标题
  labs(
    title = "CN12",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN12.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN12.pdf", p_detailed, width = 10, height = 8)






plot_data_CN13  <- subset(all_cn15_count_df, CN == "CN13")
print(plot_data_CN13)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN13, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.01, 0.08),
    breaks = seq(0.01, 0.08, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN13",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN13.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN13.pdf", p_detailed, width = 10, height = 8)






plot_data_CN14  <- subset(all_cn15_count_df, CN == "CN14")
print(plot_data_CN14)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN14, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6, 
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 3, 
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.05, 0.12),
    breaks = seq(0.05, 0.12, by = 0.02)
  ) +
  # 标签和标题
  labs(
    title = "CN14",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5, 
      size = 18, 
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN14.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN14.pdf", p_detailed, width = 10, height = 8)






plot_data_CN15  <- subset(all_cn15_count_df, CN == "CN15")
print(plot_data_CN15)
# 更详细的版本，包含更多自定义选项
p_detailed <- ggplot(plot_data_CN15, aes(x = group, y = freq, fill = group)) +
  # 小提琴图（带透明度）
  geom_violin(
    alpha = 0.6,
    width = 0.9,
    trim = TRUE,
    scale = "width"
  ) +
  # 箱线图（调整位置和样式）
  geom_boxplot(
    width = 0.15,
    alpha = 0.9,
    outlier.color = "red",
    outlier.shape = 16,
    outlier.size = 2,
    position = position_dodge(0.9)
  ) +
  # 添加均值点
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    color = "white"
  ) +
  # 组间统计比较
  stat_compare_means(
    method = "t.test",  # 根据数据正态性选择 t.test 或 wilcox.test
    label = "p.format",
    tip.length = 0.02,
    bracket.size = 0.6,
    step.increase = 0.1
  ) +
  # 颜色方案
  scale_fill_manual(
    values = c(
      "24M" = "#fdc4b6",    # 青色
      "4M" = "#e8a0b8"      # 深蓝色
    )
  ) +
  # y轴范围（根据您的数据调整）
  scale_y_continuous(
    limits = c(0.04, 0.08),
    breaks = seq(0.04, 0.08, by = 0.01)
  ) +
  # 标签和标题
  labs(
    title = "CN15",
    x = "",
    y = "Expression Level",
    fill = "group"
  ) +
  # 主题美化
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold",
      margin = margin(b = 15)
    ),
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      size = 12,
      angle = 0,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# 显示图形
print(p_detailed)

# 保存图形
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN15.png", p_detailed, width = 10, height = 8, dpi = 300)
ggsave("spatial_analysis_remove_erythrocyte-2026-01-06/CN analysis - knn50/CN15-CN15.pdf", p_detailed, width = 10, height = 8)













# Proximity analysis -------
### 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

All_type_Proximity_out_dir <- file.path(spatial_out_dir, "Proximity analysis")
dir.create(All_type_Proximity_out_dir)

cells_position_ls <- 
  all_data_filt@images %>%
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) # %>% 
# setNames(., str_remove(names(.), "\\..*$"))
# setNames(., sample_mapping[names(.)])



# all_data_cn_meta
all_data_meta_pos <- 
  cells_position_ls %>% 
  imap(function(cells_position, slide){# browser()
    all_data_meta <- 
      all_data_filt[[]] %>% 
      dplyr::mutate(annotation_full = annotation)
    proj_meta_sub <- 
      all_data_meta %>% 
      dplyr::select(annotation_full, group, project) %>% 
      rownames_to_column("cell_id") %>% 
      dplyr::filter(project == slide)
    cells_position[, 1:3] %>% 
      inner_join(proj_meta_sub, by = "cell_id") %>% 
      mutate(annotation_full = factor(annotation_full)) %>% 
      dplyr::rename(Classification1 = annotation_full, 
                    X = x, 
                    Y = y, 
                    Slide_ID = project)
  }) %>% 
  enframe() %>% 
  unnest(cols = value) %>% 
  dplyr::select(-name) %>% 
  column_to_rownames("cell_id")

all_data_meta_pos$Slide_ID <- str_c(all_data_meta_pos$group, all_data_meta_pos$Slide_ID, sep = "_")

all_data_meta_pos_subset <- all_data_meta_pos %>% rownames_to_column("cell_id")
# readRDS('data/all_data_meta_pos_subset.rds')

all_data_meta_pos_subset %>% head

# cell_id        X        Y Classification1 group       Slide_ID
# 1 1_aaaaaagg-1 1670.389 185.1358      Macrophage  CD-F X_202413076_32
# 2 1_aaadlghn-1 1681.172 176.4662          Plasma  CD-F X_202413076_32
# 3 1_aaadmplc-1 1677.561 182.6699               T  CD-F X_202413076_32
# 4 1_aaadoogj-1 1677.109 165.6637          Plasma  CD-F X_202413076_32
# 5 1_aaadpljf-1 1671.887 176.8148          Plasma  CD-F X_202413076_32
# 6 1_aaafaiko-1 1719.466 192.8071          Plasma  CD-F X_202413076_32
# ......


### 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

#### 2.1 proximity analysis --------
res_pa <- 
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset,
    query = all_data_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 'auto', # 'auto',500
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 8 / all_data_filt@images$`24M-1`@scale.factors$spot, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
saveRDS(res_pa, file = file.path(All_type_Proximity_out_dir, "res_pa.rds"))
# res_pa <- readRDS(file.path(All_type_Proximity_out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')

# summarized by group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos_subset,
#     str_c(Classification1, group, Slide_ID, sep = '..')
#   )
pa_groups <- 
  all_data_meta_pos_subset %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE) 

pa_groups2 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn10[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

pa_groups3 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn15[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
res_pa_df2 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups2, FUN = "median")})
res_pa_df3 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups3, FUN = "median")})


res_pa_df  <- 
  res_pa_df %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project"), sep = "\\.\\.") 

res_pa_df2  <- 
  res_pa_df2 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.")
res_pa_df2$CN <- factor(res_pa_df2$CN, levels = paste0("CN", 1:10))

res_pa_df3  <- 
  res_pa_df3 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df3$CN <- factor(res_pa_df3$CN, levels = paste0("CN", 1:15))


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_res_pa_ls_cn10 <-
  res_pa_df2 %>%
  with(split(res_pa_df2, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

p_res_pa_ls_cn15 <-
  res_pa_df3 %>%
  with(split(res_pa_df3, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

# p_res_pa_ls <-
#   res_pa_df2 %>%
#   with(split(res_pa_df2, CN2)) %>%
#   imap(function(dat, name){
#     dat %>%
#   # dplyr::filter(Target == "Epithelial cells") %>%
#   ggplot(aes(x = Others, y = cell_props, fill = Group)) +
#   geom_boxplot(width = 0.45, outlier.colour = NA) +
#   # geom_jitter(width = 0.25) +
#   facet_wrap(~ Target, scales = "free_y", ncol = 2)+
#   ggpubr::theme_pubclean() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   ggtitle(name)
#   })

# export plot
pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis.pdf"),
  height = 50,
  width  = 25
)
p_res_pa%>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN10.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn10 %>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN15.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn15 %>% print()
dev.off()










# 修改图片
# res_pa_df

shapiro_results <- res_pa_df %>%
  group_by(Target) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(cell_props)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(cell_props)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, file.path(All_type_Proximity_out_dir, "shapiro_test_results.csv"), row.names = FALSE)







#### 2.2 proximity analysis 2 -------
mpp <- 8 / all_data_filt@images$`24M-1`@scale.factors$spot

res_pa_ls <- list()

res_pa_ls[['10']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 10, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

res_pa_ls[['20']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 20, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

res_pa_ls[['50']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 50, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['100']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 100, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['150']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 150, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

res_pa_ls[['200']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", #NULL: 提前返回window的细胞近邻的list
    k = 1200, # 'auto',
    searchtype = "radius",
    radius = 200, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

saveRDS(res_pa_ls, file.path(All_type_Proximity_out_dir, "res_pa_ls.rds"))


# epithelial proximity
unique(all_data_meta_pos$Classification1)

proximity_cells <- c(
  "Osteoblasts", 
  "Chondrocytes", 
  "Osteoclasts"
)

proximity_all_data_ls <- 
  proximity_cells %>%
  setNames(., .) %>%
  map(function(cell){
    cell_idx <- all_data_meta_pos$Classification1 == cell
    
    stopifnot(sum(cell_idx) > 0)
    
    # 1 all_data
    all_data_sub  <- all_data_meta_pos[cell_idx, ]
    
    # 2 res_pa
    res_pa_sub_ls <- res_pa_ls %>% map(~ .x[, cell_idx])
    
    # distance diff
    diff_idx_ls <- seq_along(res_pa_sub_ls) %>% rev() %>% { matrix(c(.[-length(.)], .[-1]), ncol=2)} %>% apply(1, c, simplify = FALSE)
    
    
    # diff_idx_ls
    res_pa_sub_diff_ls <- vector(length = 5, mode = "list")
    
    for(idx in diff_idx_ls){
      res_pa_sub_diff_ls[[ idx[1] ]] <- res_pa_sub_ls[[ idx[1] ]] - res_pa_sub_ls[[ idx[2] ]]
    }
    res_pa_sub_diff_ls[[1]] <- res_pa_sub_ls[[1]]
    names(res_pa_sub_diff_ls) <- names(res_pa_sub_ls)
    
    res_pa_sub_count <-
      res_pa_sub_diff_ls %>% 
      map(~ apply(.x, 1, function(x) tapply(x, all_data_sub$Slide_ID, median)))  
    
    res_pa_sub_count <-
      res_pa_sub_count %>% 
      imap(function(x, y){
        x <- x %>% as.data.frame(check.names = FALSE) %>% rownames_to_column('Slide_ID')
        x$distance <- y
        x
      }) %>% 
      bind_rows() 
    
    
    # add Group 
    all_data_group_sub <- all_data_sub[c('group', 'Slide_ID')] %>% distinct()
    
    # return
    list(
      all_data = all_data_group_sub,
      res_pa   = res_pa_sub_count
    )
  })


# proximity_all_data_ls$Fibroblast$res_pa
# # res_pa_epi_count %>%
#   pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
#   left_join(all_data_group_epi, by = "Slide_ID") 


p_proximity_cell_count_distribution_ls <-
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      mutate(distance = as.numeric(distance)) %>%
      
      ggplot(aes(x = Cell, y = count, fill = group)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ distance, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # scale_fill_manual(values = grp_cols_sp) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p 
  })



res_pa_count_tidy_ls <- 
  proximity_all_data_ls %>%
  map(function(dat_ls){ # browser()
    
    all_data_group_sub <- dat_ls$all_data %>% map_dfc(as.character)
    
    res_pa_count_raw <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID") %>%
      mutate(distance = as.numeric(distance)) 
    
    
    res_pa_count_tidy <-
      1:nrow(res_pa_count_raw) %>%  
      map(function(x){
        count <- res_pa_count_raw$count[x] 
        if(count < 1) return(NULL) # <- count是median count，有可能是0.5
        xx <- res_pa_count_raw[x, ] %>% unlist(use.names = F) %>% rep(count) %>% matrix(nrow = count, byrow = T) %>% as.data.frame()
        xx
      }) %>%
      purrr::discard(is.null) %>%
      bind_rows()
    colnames(res_pa_count_tidy) <- colnames(res_pa_count_raw)
    
    res_pa_count_tidy$count <- as.numeric(res_pa_count_tidy$count)
    
    res_pa_count_tidy
  })


# only selecet two macrohpages
p_proximity_cell_count_distribution_ls2 <-
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      # mutate(distance = as.numeric(distance)) %>%
      mutate(distance = factor(distance, levels = c("10", "20", "50", "100", "150", "200"))) %>%
      
      # dplyr::filter(Cell %in% c("C1QC Macrophage", "SELENOP Macrophage")) %>%
      
      ggplot(aes(x = Cell, y = count, fill = distance)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ group, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # scale_fill_manual(values = grp_cols_sp) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p 
  })

# # show
# head(res_pa_epi_count_tidy, 10)

p_proximity_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){

    p <-
      res_pa_count_tidy %>%
      group_by(Cell, distance, group) %>%
      summarise(count = median(count)) %>%
      
      mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
      
      ggplot(aes(x = distance, y = count, color = group)) +
      # geom_density(alpha = 0.6) +
      geom_smooth(method = 'loess') +
      facet_wrap(~ Cell, scales = "free_y") +
      ggpubr::theme_pubr() +
      theme(
        strip.background = element_blank(),
        panel.grid = element_blank()
      ) +
      # scale_color_manual(values = grp_cols_sp) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })


# 修改后的代码，强制包含所有刻度
p_proximity_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    
    # 先进行数据转换
    plot_data <- res_pa_count_tidy %>%
      group_by(Cell, distance, group) %>%
      summarise(count = median(count), .groups = 'drop') %>%
      mutate(distance = as.numeric(str_remove(distance, "um$")))
    
    # 检查数据
    print(paste("数据集:", name))
    print("转换后的distance值:")
    print(unique(plot_data$distance))
    
    p <-
      plot_data %>%
      ggplot(aes(x = distance, y = count, color = group)) +
      geom_smooth(method = 'loess') +
      facet_wrap(~ Cell, scales = "free_y") +
      ggpubr::theme_pubr() +
      theme(
        strip.background = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      # 关键：手动设置x轴刻度
      scale_x_continuous(
        breaks = unique(plot_data$distance),  # 使用所有唯一的distance值
        labels = paste0(unique(plot_data$distance), "um")
      ) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })



# export pdf
pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - barplot-1.pdf"), width = 18, height = 8)
p_proximity_cell_count_distribution_ls %>% walk(print)
dev.off()

pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - barplot-2.pdf"), width = 18, height = 6)
p_proximity_cell_count_distribution_ls2 %>% walk(print)
dev.off()

# # export pdf
# pdf("proximity analysis - bxolplot - combined.pdf", width = 18, height = 8)
# print(p_proximity_cell_count_distribution)
# dev.off()


# export pdf
pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - smooth_lines.pdf"), width = 20, height = 12)
p_proximity_ls %>% walk(print)
dev.off()


# radar plot 
p_proximity_radar_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    p <-
      res_pa_count_tidy %>%
      
      # mutate(distance = factor(distance, levels = c('10', '50', '100', '150', '200' ))) %>%
      mutate(distance = as.numeric(distance)) %>%
      group_by(group, Cell, distance) %>%
      summarise(count = median(count)) %>%
      ggplot(aes(x = Cell, y = distance, size = count, color = group)) +
      geom_vline(xintercept = 0:12 + 0.5, color = "gray") +
      geom_hline(yintercept = c(10, 20, 50, 100, 150, 200), color = "gray") +
      geom_point(position = position_jitterdodge(
        jitter.width = 0,
        jitter.height = 10,
        dodge.width = 0.75
      )) +
      geom_text(
        data = data.frame(x = 0, y = c(10, 20, 50, 100, 150, 200)),
        aes(x = x, y = y, label = paste0(y, "um")),
        inherit.aes = FALSE,
        vjust = -1.5
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(1,1,1,6, 'mm')
      ) +
      scale_size(range = c(3, 12)) +
      scale_y_continuous(limits = c(-100, 220)) +
      # scale_color_manual(values = grp_cols_sp)+
      coord_polar(theta = 'x', clip = 'off') +
      labs(x = NULL, y = NULL, title = name)
    p 
  })





# export pdf
pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - radarplot.pdf"), width = 14, height = 12)
p_proximity_radar_ls %>% walk(print)
dev.off()

















# 20260128 analysis -------
Date <- "2026-01-28"
spatial_out_dir <- str_glue("spatial_analysis_remove_erythrocyte-{Date}")
dir.create(spatial_out_dir)

all_data <- readRDS("data_for_analysis/all_data_filt.rds")

all_data@meta.data$annotation <- all_data@meta.data$rctd_20251201
all_data_filt <- all_data[, !all_data@meta.data$annotation %in% c("Other")]
unique(all_data_filt$annotation)

# Proximity analysis -------
### 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

All_type_Proximity_out_dir <- file.path(spatial_out_dir, "Proximity analysis")
dir.create(All_type_Proximity_out_dir)

cells_position_ls <- 
  all_data_filt@images %>%
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) # %>% 
# setNames(., str_remove(names(.), "\\..*$"))
# setNames(., sample_mapping[names(.)])



# all_data_cn_meta
all_data_meta_pos <- 
  cells_position_ls %>% 
  imap(function(cells_position, slide){# browser()
    all_data_meta <- 
      all_data_filt[[]] %>% 
      dplyr::mutate(annotation_full = annotation)
    proj_meta_sub <- 
      all_data_meta %>% 
      dplyr::select(annotation_full, group, project) %>% 
      rownames_to_column("cell_id") %>% 
      dplyr::filter(project == slide)
    cells_position[, 1:3] %>% 
      inner_join(proj_meta_sub, by = "cell_id") %>% 
      mutate(annotation_full = factor(annotation_full)) %>% 
      dplyr::rename(Classification1 = annotation_full, 
                    X = x, 
                    Y = y, 
                    Slide_ID = project)
  }) %>% 
  enframe() %>% 
  unnest(cols = value) %>% 
  dplyr::select(-name) %>% 
  column_to_rownames("cell_id")

all_data_meta_pos$Slide_ID <- str_c(all_data_meta_pos$group, all_data_meta_pos$Slide_ID, sep = "_")

all_data_meta_pos_subset <- all_data_meta_pos %>% rownames_to_column("cell_id")
# readRDS('data/all_data_meta_pos_subset.rds')

all_data_meta_pos_subset %>% head

# cell_id        X        Y Classification1 group       Slide_ID
# 1 1_aaaaaagg-1 1670.389 185.1358      Macrophage  CD-F X_202413076_32
# 2 1_aaadlghn-1 1681.172 176.4662          Plasma  CD-F X_202413076_32
# 3 1_aaadmplc-1 1677.561 182.6699               T  CD-F X_202413076_32
# 4 1_aaadoogj-1 1677.109 165.6637          Plasma  CD-F X_202413076_32
# 5 1_aaadpljf-1 1671.887 176.8148          Plasma  CD-F X_202413076_32
# 6 1_aaafaiko-1 1719.466 192.8071          Plasma  CD-F X_202413076_32
# ......


### 2. proximity analysis -------
# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification 一般是细胞的注释
# k 计算的最近邻的细胞数
# searchtype 邻近度使用radius

# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325
# returntype 可以每个细胞多少radius内的细胞计数count，也可以返回细胞频率freq
# cell_enough_ratio_cutoff 为了降低计算时间，k值不建议设置太高，
# 如果k值未能完全覆盖全部的邻近细胞，则cell_enough_ratio_cutoff会自动加k值，
# 直到默认90%的细胞均已经计算了全部邻近细胞

#### 2.1 proximity analysis --------
res_pa <- 
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset,
    query = all_data_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 'auto', # 'auto',500
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 8 / all_data_filt@images$`24M-1`@scale.factors$spot, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
saveRDS(res_pa, file = file.path(All_type_Proximity_out_dir, "res_pa.rds"))
# res_pa <- readRDS(file.path(All_type_Proximity_out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS('data/res_pa.rds')

# summarized by group and celltype
# pa_groups <-
#   with(
#     all_data_meta_pos_subset,
#     str_c(Classification1, group, Slide_ID, sep = '..')
#   )
pa_groups <- 
  all_data_meta_pos_subset %>%
  group_by(Slide_ID) %>% # compute pa_groups per Slide_ID
  nest() %>%
  deframe() %>%
  imap(function(dat, Slide_ID){
    with(
      dat,
      str_c(Classification1, group, Slide_ID, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE) 

pa_groups2 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn10[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

pa_groups3 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn15[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
res_pa_df2 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups2, FUN = "median")})
res_pa_df3 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups3, FUN = "median")})


res_pa_df  <- 
  res_pa_df %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project"), sep = "\\.\\.") 

res_pa_df2  <- 
  res_pa_df2 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.")
res_pa_df2$CN <- factor(res_pa_df2$CN, levels = paste0("CN", 1:10))

res_pa_df3  <- 
  res_pa_df3 %>% 
  as.data.frame() %>% 
  rownames_to_column('group') %>% 
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>% 
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df3$CN <- factor(res_pa_df3$CN, levels = paste0("CN", 1:15))


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_res_pa_ls_cn10 <-
  res_pa_df2 %>%
  with(split(res_pa_df2, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

p_res_pa_ls_cn15 <-
  res_pa_df3 %>%
  with(split(res_pa_df3, CN)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      ggplot(aes(x = Others, y = cell_props, fill = group)) +
      geom_boxplot(width = 0.45, outlier.colour = NA) +
      # geom_jitter(width = 0.25) +
      facet_wrap(~ Target, scales = "free_y", ncol = 2)+
      ggpubr::theme_pubclean() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(name)
  })

# p_res_pa_ls <-
#   res_pa_df2 %>%
#   with(split(res_pa_df2, CN2)) %>%
#   imap(function(dat, name){
#     dat %>%
#   # dplyr::filter(Target == "Epithelial cells") %>%
#   ggplot(aes(x = Others, y = cell_props, fill = Group)) +
#   geom_boxplot(width = 0.45, outlier.colour = NA) +
#   # geom_jitter(width = 0.25) +
#   facet_wrap(~ Target, scales = "free_y", ncol = 2)+
#   ggpubr::theme_pubclean() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   ggtitle(name)
#   })

# export plot
pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - 50um.pdf"),
  height = 50,
  width  = 25
)
p_res_pa%>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN10.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn10 %>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "proximity analysis - CN15.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_ls_cn15 %>% print()
dev.off()










# 修改图片
# res_pa_df

shapiro_results <- res_pa_df %>%
  group_by(Target) %>%
  summarise(
    n = n(),  # 样本数量
    shapiro_statistic = ifelse(n() >= 3, shapiro.test(cell_props)$statistic, NA),  # W统计量
    shapiro_p_value = ifelse(n() >= 3, shapiro.test(cell_props)$p.value, NA),     # p值
    .groups = "drop"
  ) %>%
  mutate(
    normality = case_when(
      is.na(shapiro_p_value) ~ "样本量不足",
      shapiro_p_value > 0.05 ~ "符合正态分布",
      TRUE ~ "不符合正态分布"
    )
  )

# 查看结果
print(shapiro_results, n = Inf)

# 保存结果到CSV文件
write.csv(shapiro_results, file.path(All_type_Proximity_out_dir, "shapiro_test_results.csv"), row.names = FALSE)







#### 2.2 proximity analysis 2 -------
mpp <- 8 / all_data_filt@images$`24M-1`@scale.factors$spot

res_pa_ls <- list()

res_pa_ls[['10']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 10, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

res_pa_ls[['20']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 20, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

res_pa_ls[['50']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 50, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

res_pa_ls[['100']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 100, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

res_pa_ls[['150']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 150, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

res_pa_ls[['200']] <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", #NULL: 提前返回window的细胞近邻的list
    k = 1200, # 'auto',
    searchtype = "radius",
    radius = 200, # <--------
    mpp = mpp, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

saveRDS(res_pa_ls, file.path(All_type_Proximity_out_dir, "res_pa_ls.rds"))


# epithelial proximity
unique(all_data_meta_pos$Classification1)

proximity_cells <- c(
  "Osteoblasts", 
  "Chondrocytes", 
  "Osteoclasts"
)

proximity_all_data_ls <- 
  proximity_cells %>%
  setNames(., .) %>%
  map(function(cell){
    cell_idx <- all_data_meta_pos$Classification1 == cell
    
    stopifnot(sum(cell_idx) > 0)
    
    # 1 all_data
    all_data_sub  <- all_data_meta_pos[cell_idx, ]
    
    # 2 res_pa
    res_pa_sub_ls <- res_pa_ls %>% map(~ .x[, cell_idx])
    
    # distance diff
    diff_idx_ls <- seq_along(res_pa_sub_ls) %>% rev() %>% { matrix(c(.[-length(.)], .[-1]), ncol=2)} %>% apply(1, c, simplify = FALSE)
    
    
    # diff_idx_ls
    res_pa_sub_diff_ls <- vector(length = 5, mode = "list")
    
    for(idx in diff_idx_ls){
      res_pa_sub_diff_ls[[ idx[1] ]] <- res_pa_sub_ls[[ idx[1] ]] - res_pa_sub_ls[[ idx[2] ]]
    }
    res_pa_sub_diff_ls[[1]] <- res_pa_sub_ls[[1]]
    names(res_pa_sub_diff_ls) <- names(res_pa_sub_ls)
    
    res_pa_sub_count <-
      res_pa_sub_diff_ls %>% 
      map(~ apply(.x, 1, function(x) tapply(x, all_data_sub$Slide_ID, median)))  
    
    res_pa_sub_count <-
      res_pa_sub_count %>% 
      imap(function(x, y){
        x <- x %>% as.data.frame(check.names = FALSE) %>% rownames_to_column('Slide_ID')
        x$distance <- y
        x
      }) %>% 
      bind_rows() 
    
    
    # add Group 
    all_data_group_sub <- all_data_sub[c('group', 'Slide_ID')] %>% distinct()
    
    # return
    list(
      all_data = all_data_group_sub,
      res_pa   = res_pa_sub_count
    )
  })


# proximity_all_data_ls$Fibroblast$res_pa
# # res_pa_epi_count %>%
#   pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
#   left_join(all_data_group_epi, by = "Slide_ID") 


p_proximity_cell_count_distribution_ls <-
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      mutate(distance = as.numeric(distance)) %>%
      
      ggplot(aes(x = Cell, y = count, fill = group)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ distance, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # scale_fill_manual(values = grp_cols_sp) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p 
  })



res_pa_count_tidy_ls <- 
  proximity_all_data_ls %>%
  map(function(dat_ls){ # browser()
    
    all_data_group_sub <- dat_ls$all_data %>% map_dfc(as.character)
    
    res_pa_count_raw <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID") %>%
      mutate(distance = as.numeric(distance)) 
    
    
    res_pa_count_tidy <-
      1:nrow(res_pa_count_raw) %>%  
      map(function(x){
        count <- res_pa_count_raw$count[x] 
        if(count < 1) return(NULL) # <- count是median count，有可能是0.5
        xx <- res_pa_count_raw[x, ] %>% unlist(use.names = F) %>% rep(count) %>% matrix(nrow = count, byrow = T) %>% as.data.frame()
        xx
      }) %>%
      purrr::discard(is.null) %>%
      bind_rows()
    colnames(res_pa_count_tidy) <- colnames(res_pa_count_raw)
    
    res_pa_count_tidy$count <- as.numeric(res_pa_count_tidy$count)
    
    res_pa_count_tidy
  })


# only selecet two macrohpages
p_proximity_cell_count_distribution_ls2 <-
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      # mutate(distance = as.numeric(distance)) %>%
      mutate(distance = factor(distance, levels = c("10", "20", "50", "100", "150", "200"))) %>%
      
      # dplyr::filter(Cell %in% c("C1QC Macrophage", "SELENOP Macrophage")) %>%
      
      ggplot(aes(x = Cell, y = count, fill = distance)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ group, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      # scale_fill_manual(values = grp_cols_sp) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p 
  })

# # show
# head(res_pa_epi_count_tidy, 10)

p_proximity_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    # p <-
    #   res_pa_count_tidy %>%
    #   mutate(distance = as.numeric(distance)) %>%
    #   
    #   ggplot(aes(x = distance, fill = Group)) + #, color = Group)) +
    #   geom_density(alpha = 0.6) +
    #   # geom_smooth(aes(y = count), method = "loess") +
    #   facet_wrap(~ Cell, scales = "free_y") +
    #   theme_bw() +
    #   scale_fill_manual(values = grp_cols_sp) +
    #   labs(
    #     x = 'Distance(um)',
    #     y = 'Median count'
    #   )
    
    p <-
      res_pa_count_tidy %>%
      group_by(Cell, distance, group) %>%
      summarise(count = median(count)) %>%
      
      mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
      
      ggplot(aes(x = distance, y = count, color = group)) +
      # geom_density(alpha = 0.6) +
      geom_smooth(method = 'loess') +
      facet_wrap(~ Cell, scales = "free_y") +
      ggpubr::theme_pubr() +
      theme(
        strip.background = element_blank(),
        panel.grid = element_blank()
      ) +
      # scale_color_manual(values = grp_cols_sp) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })




# export pdf
pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - barplot-1.pdf"), width = 18, height = 8)
p_proximity_cell_count_distribution_ls %>% walk(print)
dev.off()

pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - barplot-2.pdf"), width = 18, height = 6)
p_proximity_cell_count_distribution_ls2 %>% walk(print)
dev.off()

# # export pdf
# pdf("proximity analysis - bxolplot - combined.pdf", width = 18, height = 8)
# print(p_proximity_cell_count_distribution)
# dev.off()


# export pdf
pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - smooth_lines.pdf"), width = 14, height = 8)
p_proximity_ls %>% walk(print)
dev.off()


# radar plot 
p_proximity_radar_ls <- 
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    p <-
      res_pa_count_tidy %>%
      
      # mutate(distance = factor(distance, levels = c('10', '50', '100', '150', '200' ))) %>%
      mutate(distance = as.numeric(distance)) %>%
      group_by(group, Cell, distance) %>%
      summarise(count = median(count)) %>%
      ggplot(aes(x = Cell, y = distance, size = count, color = group)) +
      geom_vline(xintercept = 0:12 + 0.5, color = "gray") +
      geom_hline(yintercept = c(10, 20, 50, 100, 150, 200), color = "gray") +
      geom_point(position = position_jitterdodge(
        jitter.width = 0,
        jitter.height = 10,
        dodge.width = 0.75
      )) +
      geom_text(
        data = data.frame(x = 0, y = c(10, 20, 50, 100, 150, 200)),
        aes(x = x, y = y, label = paste0(y, "um")),
        inherit.aes = FALSE,
        vjust = -1.5
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(1,1,1,6, 'mm')
      ) +
      scale_size(range = c(3, 12)) +
      scale_y_continuous(limits = c(-100, 220)) +
      # scale_color_manual(values = grp_cols_sp)+
      coord_polar(theta = 'x', clip = 'off') +
      labs(x = NULL, y = NULL, title = name)
    p 
  })





# export pdf
pdf(file.path(All_type_Proximity_out_dir, "proximity analysis - radarplot.pdf"), width = 14, height = 12)
p_proximity_radar_ls %>% walk(print)
dev.off()











