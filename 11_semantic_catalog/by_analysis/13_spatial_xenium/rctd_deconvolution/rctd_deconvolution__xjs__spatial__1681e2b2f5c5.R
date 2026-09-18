##邻域分析------
#筛选掉verification组样本，T细胞按照CD4_T、CD8_T、Cyclin_T、Treg，其他按照大类分析
# CN_data <- subset(
#   xenium,
#   subset = group != "verification_verification")

CN_data <- all_data_temp
CN_data = CN_data[,!CN_data@meta.data$first_type %in% c("Unknown")]
CN_data@meta.data$first_type %>% table()
CN_data@meta.data$first_type <- factor(CN_data@meta.data$first_type)
CN_data@meta.data$first_type %>% unique()


library(Seurat)
library(tidyverse)
source("spatial analysis.R")
##data ------
CN_out_dir <- "RCTD_new/CN analysis"
dir.create(CN_out_dir)

#all_data <- readRDS("data_for_analysis/all_data_annotation_rctd.rData")
Idents(CN_data) <- "first_type"

cells_position_ls<-readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_position_ls.rds"))
cells_position_df <- cells_position_ls %>% bind_rows()

#meta.data$cell_id需要和cells_position_ls里的cell_id一致

# rownames(CN_data@meta.data) <- CN_data@meta.data$cell_id1
table1 <- data_frame(rownames(CN_data@meta.data),CN_data@meta.data$cell_id)


## CN analysis ---------
# options(warn = 0)
# Classification1 must be a factor
# add a column: obj_id

all_data_submeta_for_cn_analysis <- 
  CN_data[[c("sample", "group", "first_type")]] %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id")

colnames(all_data_submeta_for_cn_analysis)[colnames(all_data_submeta_for_cn_analysis) %in% c("x_centroid", "y_centroid")] <- c("X", "Y")

all_data_submeta_for_cn_analysis <-all_data_submeta_for_cn_analysis %>% 
  mutate(
    first_type = factor(first_type), # import!
    obj_id = cell_id # import!
  )

### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "first_type",
  k.nn           = 10,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))

### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn10_heat <-
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
cowplot::plot_grid(p_cn10_heat$gtable)
dev.off()


### CN k15 ---------
res_cn_15 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "first_type",
  k.nn           = 10,
  k.cluster      = 15
)

saveRDS(res_cn_15, file.path(CN_out_dir, "res_cn_15.rds"))

### heatmap ----------------------
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
  height = 9,
  width  = 10
)
cowplot::plot_grid(p_cn15_heat$gtable)
dev.off()



### CN annotaton --------------------
# CN k10
cn_annotation <- c(
  "CN1"  = "CAF_Myeloid",
  "CN2"  = "CAF_Mast",
  "CN3"  = "Stroma",
  "CN4"  = "immune",
  "CN5"  = "immune",
  "CN6"  = "immune",
  "CN7"  = "Basal",
  "CN8"  = "Tumor_MoMac",
  "CN9"  = "Tumor_MoMac",
  "CN10" = "Tumor"
)

p_cn_heat2 <-
  cn10_cluster_FC %>%
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
  file.path(CN_out_dir, "CN_heatmap_annotation.pdf"),
  height = 6,
  width  = 11
)
cowplot::plot_grid(p_cn_heat2$gtable)
dev.off()


### CN10 spatial plot -----------------
cn10_sp_out_dir <- file.path(CN_out_dir, 'CN10 spatial plot') 
dir.create(cn10_sp_out_dir)

cells_boundary_ls<-readRDS(file.path("~/R_OUTPUT/data_for_analysis", "cells_boundary_ls.rds"))
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

# CN_data@meta.data$cell_id3 <- CN_data@meta.data$cell_id
# CN_data@meta.data$cell_id <- NULL

All_type_meta_cn10 <- 
  CN_data[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn10, 'RCTD_new/CN analysis/all_data_meta_cn10.rds')
All_type_meta_cn10 <- readRDS("RCTD_new/CN analysis/all_data_meta_cn10.rds")

cn10_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn10[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$group1[anno_data$sample == project]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })


# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn10_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


###筛选掉的细胞NA的坐标绘图(未经过筛选用上一步绘图)-----
cn10_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project) {
    # 获取组信息
    anno_data   <- All_type_meta_cn10[c("sample", "group1")] %>% distinct()
    group_info  <- anno_data$group1[anno_data$sample == project]
    
    # 合并数据并清理无效坐标
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn10, by = "cell_id") %>%
      # 清理无效坐标值
      filter(
        !is.na(vertex_x), !is.na(vertex_y),
        is.finite(vertex_x), is.finite(vertex_y),
        !is.nan(vertex_x), !is.nan(vertex_y)
      )
    
    # 检查是否有有效数据
    if (nrow(dat) == 0) {
      warning("项目 ", project, " 无有效坐标数据")
      return(ggplot() + 
               geom_text(aes(0,0), label = "无有效数据", size = 12) +
               theme_void() +
               labs(title = paste(project, group_info, sep = "_")))
    }
    
    # 创建图表
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id)) +
      scale_y_reverse() +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_"))
    
    return(p)
  })

# 安全导出图像
cn10_sp_p_ls %>% 
  iwalk(function(p, name) {
    pdf_file <- file.path(cn10_sp_out_dir, paste0("CN10_spatial_plot_", name, ".pdf"))
    pdf(pdf_file, width = 20, height = 20)
    
    tryCatch({
      print(p)
    }, error = function(e) {
      message("绘图失败: ", name, " - ", e$message)
      # 创建错误提示图
      error_plot <- ggplot() + 
        geom_text(aes(0,0), label = "绘图错误", size = 12) +
        theme_void()
      print(error_plot)
    })
    
    dev.off()
  })



All_type_meta_cn10$group %>% unique()
All_type_meta_cn10$CN %>% unique()


### cell props ------
All_type_cn10_count_df <-
  All_type_meta_cn10 %>% 
  group_by(sample, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

# All_type_cn10_count_df <- factor

p_cn10_props <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

p_cn10_props_group <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  stat_compare_means(label = "p.format",size=3, aes(group = group)) +
  ggsci::scale_fill_aaas(
    name = "Group"#,  # 图例标题
    #breaks = group1_level  # 强制图例顺序与group1_level一致
  )+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_group_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn10_props %>% print
p_cn10_props_group %>% print
dev.off()

Table1 <- data_frame(All_type_meta_cn10$CN,All_type_meta_cn10$sample,All_type_meta_cn10$group4)

# all_colors <- c(
#   "red", "#A6CEE3", "#5B9EC9", "#FF6B3A", "#2D82AF", #0-4
#   "#98D277","#52AF43", "#6F9E4C", "#FDBB69", "#F78620",#5-9
#   "#00F5A7", "#F15BB5", "#7D54A5", "#02c9c9", "#600473",#10-14
#   "#B294C7", "pink", "#9E8099", "blue", "#DBB466",#15-19
#   "#B15928", "#4ECDC4", "#F0EB99"          #20-22,other
# )
# 
# cols_brewer = c("#f94144","#f25c27","#eb770a","#f89f12","#f9c74f","#90be6d",
#                 "#43aa8b","#4d908e","#577590","#277da1","#967aa1",
#                 "#F2698D","#F391AE","#F3B9CF","#AF8EB8","#7B84B7","#3F6FAF",
#                 "#035AA6","#0473B4","#8491B7","#737AAB","#5E5E97","gray90")
# 
# group4_colors <- c(
#   "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", 
#   "#9467BD", "#8C564B", "#E377C2", "#7F7F7F"
# )
# 
# p_cn10_props_group <-
#   All_type_cn10_count_df %>%
#   
#   ggplot(aes(x = group4, y = freq, fill = group4)) +
#   geom_boxplot() +
#   ggsci::scale_fill_aaas(
#     name = "Group"#,  # 图例标题
#   )+
#   theme(axis.text.x = element_blank()
#   )+
#   labs(
#     x = NULL
#   )+
#   facet_wrap(~CN, nrow = 1)
# 
# 
# pdf(
#   file.path(CN_out_dir, "CN10_boxplot_group4.pdf"),
#   height = 6,
#   width  = 15
# )
# 
# p_cn10_props_group %>% print
# dev.off()


### CN15 spatial plot -----------------
cn15_sp_out_dir <- file.path(CN_out_dir, 'CN15 spatial plot') 
dir.create(cn15_sp_out_dir)

cells_boundary_ls <- readRDS("~/data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn15 <- 
  CN_data[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn15, file.path(CN_out_dir, "all_data_meta_cn15.rds"))

cn15_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn15[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$sample == project]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn15, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })


# export image
cn15_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn15_sp_out_dir, str_glue("CN15 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

### cell props ------
All_type_cn15_count_df <-
  All_type_meta_cn15 %>% 
  group_by(sample, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn15_props <-
  All_type_cn15_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN)) +
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
  All_type_cn15_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  stat_compare_means(label = "p.format",size=3,aes(group = group))+
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN15_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn15_props %>% print
p_cn15_props_group %>% print
dev.off()



#邻近度和最近邻分析----
# cells_boundary_ls<-readRDS(file.path("~/data_for_analysis", "cells_boundary_ls.rds"))
# cells_position_ls<-readRDS(file.path("~/data_for_analysis", "cells_position_ls.rds"))

# library(tidyverse)
# 
# # 从cells_position_ls中提取cell_id → 样本名的映射
# cell_name_mapping <- imap_dfr(
#   cells_position_ls,
#   ~ tibble(cell_id = .x$cell_id, name = .y)  # .y是列表元素名（样本名），.x是tibble
# )
# 
# # 合并映射表到元数据
# CN_data@meta.data <- CN_data@meta.data %>%
#   left_join(cell_name_mapping, by = "cell_id")
# 
# rownames(CN_data@meta.data) <- CN_data@meta.data$cell_id

###meta.data$cell_id需要和cells_position_ls里的cell_id一致
CN_data@meta.data %>% head()
CN_data@meta.data$cell_id <- rownames(CN_data@meta.data)

Table1 <- data_frame(rownames(CN_data@meta.data),CN_data@meta.data$cell_id,CN_data@meta.data$sample,CN_data@meta.data$name)

table(CN_data@meta.data$sample,CN_data@meta.data$name)

library(Seurat)
library(tidyverse)

source("spatial analysis.R")

proximity_out_dir <- "RCTD_new/proximity"
dir.create(proximity_out_dir)

###sample名称和slide名称需一致，slide名称是cells_position_ls各FOV的名字
###name列即为对应cells_position_ls的样本名sample列

### Proximity analysis -------
#### 1. input data ------------
cells_position_ls %>% head() 

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){  #browser()
    
    all_data_meta <- CN_data[[]]  
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(first_type, sample, "cell_id", "group") %>% 
      #rownames_to_column("cell_id") %>%
      dplyr::filter(sample == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(first_type = factor(first_type)) %>%
      dplyr::rename(Classification1 = first_type,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = sample,
                    Group = group)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos %>% head()
all_data_meta_pos$cell_id <- rownames(all_data_meta_pos)
all_data_meta_pos <- all_data_meta_pos %>%
  as.data.frame() %>%       # 确保为数据框类型
  remove_rownames()     # 强制清除现有行名


#### ########### only select 100,000 cells ##############
# all_data_meta_pos_subset <- 
#   all_data_meta_pos %>% 
#   rownames_to_column("cell_id") %>%
#   slice_sample(n = 100000)


#### select ROI cells --------------------------------------------------------------
# all_tumor_roi_selected <- readRDS('tumor_roi_selected_result3/all_tumor_roi_selected3.rds')
# all_data_meta_pos_subset <- 
#   all_data_meta_pos %>% 
#   rownames_to_column("cell_id") %>% 
#   dplyr::filter(cell_id %in% Cells(all_tumor_roi_selected))

# all_data_meta_pos_subset %>% head

# cell_id        X        Y Classification1 group       Slide_ID
# 1 1_aaaaaagg-1 1670.389 185.1358      Macrophage  CD-F X_202413076_32
# 2 1_aaadlghn-1 1681.172 176.4662          Plasma  CD-F X_202413076_32
# 3 1_aaadmplc-1 1677.561 182.6699               T  CD-F X_202413076_32
# 4 1_aaadoogj-1 1677.109 165.6637          Plasma  CD-F X_202413076_32
# 5 1_aaadpljf-1 1671.887 176.8148          Plasma  CD-F X_202413076_32
# 6 1_aaafaiko-1 1719.466 192.8071          Plasma  CD-F X_202413076_32
# ......



## 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 'auto', #750 
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )


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
pa_groups <-
  with(
    all_data_meta_pos,
    str_c(Classification1, Group, Slide_ID, sep = '..')
  )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.35, outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y")+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# export plot
pdf(
  file.path("RCTD_new/proximity", "proximity analysis.pdf"),
  height = 14,
  width  = 18
)
p_res_pa %>% print()
dev.off()



## 2.2 proximity analysis - series analysis -----

# # 删除 Slide_ID 为 "P5_POST2" 的行
# all_data_meta_pos <- all_data_meta_pos[all_data_meta_pos$Slide_ID != "P5_POST2", ]

# # 删除 Slide_ID 为 "P5_POST2" 的行
# all_data_meta_pos <- all_data_meta_pos %>%
#   filter(Slide_ID != "P5_POST2")

init_cells <- c(
  '10'  = 750,
  '50'  = 800,
  '100' = 1000,
  '150' = 1200,
  '200' = 1500
)

res_pa_ls <-
  c(10, 50, 100, 150, 200) %>%
  set_names(., .) %>%
  map(function(x){
    k = init_cells[as.character(x)]
    res_pa <-
      ProximutyAnalysis_NN_multi(
        all_data_meta_pos,
        query = all_data_meta_pos,
        x = "X",
        y = "Y",
        group.by = "Slide_ID",
        classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
        k = 'auto', #k,
        searchtype = "radius",
        radius = x, # 邻域半径
        mpp = 1, # !important, 1 for xenium
        returntype = "count",
        cell_enough_ratio_cutoff = 90
      )
    res_pa
  })


# res_pa_ls <- readRDS("data/res_pa_ls.rds")
# 将同心圆统计 转换成 环状同心的统计
# 注意一定要倒序操作
# 200 = 200 - 150
# 150 = 150 - 100 
# ...
# 
res_pa_ls2 <- res_pa_ls
diff_idx_ls <-
  c(10, 50, 100, 150, 200) %>% 
  as.character() %>%
  rev() %>%
  {matrix(c(.[-length(.)], .[-1]), ncol = 2)} %>%
  apply(1, c, simplify = FALSE)

for(i in diff_idx_ls){
  first  <- i[1]
  second <- i[2]
  res_pa_ls2[[ first ]] <- res_pa_ls2[[ first ]] - res_pa_ls2[[ second ]] 
}

# count to freq
res_pa_ls2 <- 
  res_pa_ls2 %>%
  map(function(x){# browser()
    x_colsum <- colSums(x)
    xx <- sweep(x, 2, x_colsum, FUN = "/")
    if(any(x_colsum == 0)){
      xx[, x_colsum == 0] <- 0
    }
    xx
  })

# 每个样本中的每种Target的其他细胞的计数取中位数
res_pa_count_tidy <- 
  res_pa_ls2 %>%
  imap(function(x, distance){# browser()
    # 分组
    pa_groups <-
      with(
        all_data_meta_pos,
        str_c(Classification1, Group, Slide_ID, sep = '..')
      )
    # 分组内取中位数
    res <- apply(x, 1, function(y){tapply(y, INDEX = pa_groups, FUN = "median")})
    
    # data tidy
    res_df  <-
      res %>%
      as.data.frame() %>%
      rownames_to_column('group') %>%
      pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
      separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 
    res_df$distance <- distance
    
    # return
    res_df
  }) %>%
  bind_rows()


# plot 
# only plot "Epithelial cells" ....
res_pa_count_tidy_epi <-
  res_pa_count_tidy %>%
  mutate(distance = factor(distance, levels = sort(as.numeric(unique(distance))))) %>%
  dplyr::filter(Target == "Normal")

res_pa_count_tidy_epi_median <-
  res_pa_count_tidy_epi %>%
  group_by(Others, Group, distance) %>% 
  summarise(cell_props = median(cell_props))

p_proximity_series_epi <-
  res_pa_count_tidy_epi %>%
  ggplot(aes(x = distance, y = cell_props, color = Group)) +
  geom_boxplot( ) +
  ggsci::scale_color_aaas() +
  # geom_point(position = position_jitterdodge(jitter.width = 0.01, dodge.width = 0.5)) +
  geom_point(data = res_pa_count_tidy_epi_median, position = position_dodge(width = 0.75)) +
  geom_smooth(data = res_pa_count_tidy_epi_median, aes(group = Group), position = position_dodge(width = 0.75)) +
  facet_wrap(~ Others, scales = 'free_y') +
  ggpubr::theme_pubclean()


# export plot
pdf(
  file.path("RCTD_new/proximity", "proximity analysis_series distance.pdf"),
  height = 14,
  width  = 18
)
p_proximity_series_epi %>% print()
dev.off()


## 3. Nearest neighbor analysis -------
# 由于最近邻只需要最近的那个细胞的信息
# 因此不同于proximity analysis，有以下变化：
# classification调为NULL
# searchtype 邻近度使用radius
# k不需要太大，10即可，最近邻分析只需要第二个细胞的信息
# radius适当调大100um，防止部分细胞比较离散，距离其他细胞较远
# returntype与cell_enough_ratio_cutoff也没有意义了

# 参数介绍
# x,y 指定细胞的坐标
# group.by 一般是切片或者样本编号，代表分别计算最近邻关系
# classification NULL即可，返回每一个query细胞的细胞编号、距离、query和total细胞的id
# k 计算的最近邻的细胞数
# radius是只查看当前细胞的多少um内的细胞
# mpp 坐标与um的转换比例，xenium是1，celldive是0.325


res_nn_ls <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos,
    query = all_data_meta_pos,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = NULL, # 若为NULL: 提前返回window的细胞近邻的list
    k = 10, # 'auto',
    searchtype = "radius",
    radius = 100, 
    mpp = 1, # !important, 1 for xenium
    cell_enough_ratio_cutoff = 20
  )

# res_nn_ls <- readRDS("data/res_nn_ls.rds")

res_nn_df <-
  res_nn_ls %>%
  map(function(x){ # browser()
    cell_idx  <- x$nn.idx[, 2]
    cell_idx[cell_idx == 0] <- NA # some cells have no nn cells
    cell_id   <- x$total_cell_id[cell_idx]
    cell_dists <- x$nn.dists[, 2]
    
    all_data_used <- all_data_meta_pos %>% column_to_rownames("cell_id")
    cell_anno <- all_data_used[cell_id, 'Classification1']
    
    res <-
      data.frame(
        Others = cell_anno,
        distance = cell_dists
      ) %>%
      cbind(all_data_used[x$query_cell_id, c('Classification1', 'Group', 'Slide_ID')]) %>%
      rownames_to_column("cell_id") %>%
      dplyr::rename(Target = Classification1) 
    
    res
  }) %>%
  bind_rows()

if(any(is.na(res_nn_df$Others) )){
  message(str_glue("{sum(is.na(res_nn_df$Others))} cells have no nearest neighbor cell"))
  res_nn_df <- drop_na(res_nn_df)
}

# plot
p_nn_ls <-
  res_nn_df %>%
  group_by(Target) %>%
  nest() %>%
  deframe() %>%
  imap(function(x, target){ # browser()
    x %>%
      ggplot(aes(x = distance, color = Group, fill = Group)) + 
      geom_density(alpha = 0.5) +
      facet_wrap(~ Others) +
      ggpubr::theme_pubclean() +
      labs(
        title = str_glue("Target cell: {target}"),
        x = "distance(um)"
      )
    
  })

# export plot
pdf(
  file.path("RCTD_new/nearest neighbor", "nearest neighbor analysis.pdf"),
  height = 14,
  width  = 18
)
p_nn_ls %>% print()
dev.off()

