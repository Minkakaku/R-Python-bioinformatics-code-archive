
library(tidyverse)
library(Seurat)
library(patchwork)

# setting dir -------
img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)
mk_output <- "4. marker gene" # dir.create(mk_output)
er_output <- "5. enriment analysis" # dir.create(er_output)

dir.create(img_output)
dir.create(qc_outdir)
dir.create(dr_outdir)
dir.create(sp_output)
dir.create(mk_output)
dir.create(er_output)

# parse ROI----------------
sample_df <-
  data.frame(
    sample = c("IR7_1","IR7_2","IR3_1","IR3_2","IR1_1","IR1_2","S_1","S_2"),
    # sample = c("D534","D536","D537","D538","D546","D547","s1","S2"),
    # sample = c("5_","6_","3_","4_","1_","2_","7_","8_"),
    group  = c("IR7_1","IR7_2","IR3_1","IR3_2","IR1_1","IR1_2","S_1","S_2")
  )

roi_ls <- list.files(
  'ROI',
  full.names = TRUE
)

names(roi_ls) <-
  roi_ls %>%
  basename() %>%
  str_remove("\\..*")


roi_coords_ls <-
  roi_ls %>%
  map(function(path){
    # 读取文件内容，跳过前两行
    dat <- read.csv(path,skip = 2)
    # 将列名转换为小写
    colnames(dat) <- dat %>% colnames() %>% str_to_lower()
    # 添加一个新的列 'cell'，其值与 'cell.id' 列相同
    dat<-dat %>% mutate(cell=cell.id)%>%
      # 将 'cell.id' 列转换为行名
      column_to_rownames('cell.id')
    dat
  })

### 根据每个样本所在芯片的顺序编号，来对应all_data_ls中的cell id
.id<-c("5_","6_","3_","4_","1_","2_","7_","8_")
# .id<- str_c(names(roi_coords_ls), "_")
# .id
names(.id)<-names(roi_coords_ls)
roi_coords_ls <- roi_coords_ls %>%
  map2(.id, ~ {
    .x$cell <- str_c(.y, .x$cell)
    return(.x)
  })




library(dplyr)
library(purrr)

roi_df <- roi_coords_ls %>% 
  imap(~ .x %>% mutate(sample = .y)) %>% 
  bind_rows() %>%
  remove_rownames() %>%  # 确保转换为标准数据框并重置行名
  column_to_rownames(var = "cell")
roi_df$ROI <- paste0("roi_", roi_df$sample)

head(roi_df)



# 将roi_df添加到all_data中 -------
roi_metadata <- data.frame(
  cell = colnames(all_data), 
  ROI = "Other", 
  stringsAsFactors = FALSE
) %>% 
  left_join(
    roi_df %>% rownames_to_column("cell") %>% 
      select(cell, ROI), 
    by = "cell"
  ) %>% 
  mutate(ROI = ifelse(is.na(ROI.y), ROI.x, ROI.y)) %>% 
  select(cell, ROI) %>% 
  column_to_rownames("cell")


all_data@meta.data$ROI <- roi_metadata[colnames(all_data), "ROI"]

table(all_data$ROI)

all_data_filt <- subset(all_data, subset = ROI != "Other")


dim(all_data)
dim(all_data_filt)
table(all_data_filt$ROI)



# spatial plot ------- 
sp_output_filt <- "ROI_analysis/spatial plot"
dir.create(sp_output_filt, recursive = T)

all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

# plan <- c("#00B4F0","#7CAE00","#FF64B0","#619CFF","#F8766D","#00BA38","#00C08B","#DE8C00","#C77CFF")
# plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
#           "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
#           "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
#           "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
#           "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
#           "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
#           "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
#           "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
#           "#005995", "#fa625f", "#f17d80")
plan1 <- c("gray90", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = ROI, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan1)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - ROI.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })





# 添加一列新注释，区分ROI -------
all_data@meta.data <- all_data@meta.data %>% 
  mutate(cell_type_res05_ROI = ifelse(ROI == "Other", "Other", cell_type_res05))

table(all_data$cell_type_res05_ROI)


all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = cell_type_res05_ROI, group = cell_id))  +
      scale_y_reverse() +
      coord_equal() +
      scale_fill_manual(values = plan)+
      # coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) +
      theme_void()
    p
  })


# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output_filt, str_glue("spatial plot {grp} - celltype_ROI.pdf")),
      height = 12,
      width  = 12
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })




saveRDS(all_data, "data_for_analysis/all_data_roi_1208.rds")
saveRDS(all_data_filt, "data_for_analysis/all_data_filt_roi_1208.rds")



# spatial analysis roi cell_type_res05--------
library(Seurat)
library(tidyverse)

source("spatial analysis.R")


##  all_data ------
CN_out_dir <- "ROI_analysis/spatial_analysis/CN analysis"
dir.create(CN_out_dir, recursive = T)

all_data_filt <- readRDS("data_for_analysis/all_data_filt_roi_1208.rds")
CN_data <- all_data_filt

Idents(CN_data) <- "cell_type_res05"

cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()

## CN analysis ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

all_data_submeta_for_cn_analysis <- 
  CN_data[[c("project", "group", "cell_type_res05")]] %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    cell_type_res05 = factor(cell_type_res05), # import!
    obj_id = cell_id # import!
  )

### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "cell_type_res05",
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
  group.by       = "project",
  classification = "cell_type_res05",
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
# # CN k10
# cn_annotation <- c(
#   "CN1"  = "Macrophage",
#   "CN2"  = "T",
#   "CN3"  = "B",
#   "CN4"  = "Fibro_Stroma",
#   "CN5"  = "Myo_Stroma",
#   "CN6"  = "Pan immune",
#   "CN7"  = "NK_Endothelial_Stroma",
#   "CN8"  = "Epithelial",
#   "CN9"  = "T_Epithelial",
#   "CN10" = "T_mast_Epithelial"
# )
# 
# p_cn_heat2 <-
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
#   file.path(CN_out_dir, "CN_heatmap_annotation.pdf"),
#   height = 6,
#   width  = 11
# )
# cowplot::plot_grid(p_cn_heat2$gtable)
# dev.off()


### CN10 spatial plot -----------------
cn10_sp_out_dir <- file.path(CN_out_dir, 'CN10 spatial plot') 
dir.create(cn10_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn10 <- 
  CN_data[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn10, file.path(cn10_sp_out_dir, "all_data_meta_cn10.rds"))

cn10_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn10[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
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

### cell props ------
All_type_cn10_count_df <-
  All_type_meta_cn10 %>% 
  group_by(project, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

All_type_cn10_count_df$project <- factor(All_type_cn10_count_df$project, levels = sample_map_group[,1])
All_type_cn10_count_df$group <- factor(All_type_cn10_count_df$group, levels = unique(sample_map_group[,2]))

p_cn10_props <-
  All_type_cn10_count_df %>%
  
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

library(ggpubr)
p_cn10_props_group <-
  All_type_cn10_count_df %>%
  
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

p_cn10_props_group_2 <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif",size=3,aes(group = group))+
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn10_props %>% print
p_cn10_props_group %>% print
p_cn10_props_group_2 %>% print
dev.off()



### CN15 spatial plot -----------------
cn15_sp_out_dir <- file.path(CN_out_dir, 'CN15 spatial plot') 
dir.create(cn15_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn15 <- 
  CN_data[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn15, file.path(cn15_sp_out_dir, "all_data_meta_cn15.rds"))

cn15_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn15[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
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
  group_by(project, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

All_type_cn15_count_df$project <- factor(All_type_cn15_count_df$project, levels = sample_map_group[,1])
All_type_cn15_count_df$group <- factor(All_type_cn15_count_df$group, levels = unique(sample_map_group[,2]))

p_cn15_props <-
  All_type_cn15_count_df %>%
  
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

p_cn15_props_group_2 <-
  All_type_cn15_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif",size=3,aes(group = group))+
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
p_cn15_props_group_2 %>% print
dev.off()








## Proximity analysis -------
### 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

All_type_Proximity_out_dir <- "ROI_analysis/spatial_analysis/Proximity analysis"
dir.create(All_type_Proximity_out_dir)

# all_data_cn_meta
All_type_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){#browser()
    
    proj_meta_sub <-
      CN_data[[]] %>%  
      dplyr::select(cell_type_res05, group, project) %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%
      dplyr::rename(Classification1 = cell_type_res05,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = project)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

All_type_meta_pos_subset <- All_type_meta_pos %>% rownames_to_column("cell_id")

saveRDS(All_type_meta_pos_subset, file = file.path(All_type_Proximity_out_dir, "Subtype_meta_pos_subset.rds"))

All_type_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
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
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = "auto", # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

saveRDS(res_pa, file = file.path(All_type_Proximity_out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS(file.path(All_type_Proximity_out_dir, "res_pa.rds"))


# summarized by group and celltype
pa_groups <-
  with(
    All_type_meta_pos_subset,
    str_c(Classification1, group, Slide_ID, sep = '..')
  )

pa_groups2 <-
  All_type_meta_pos_subset %>%
  left_join(All_type_meta_cn10[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

pa_groups3 <-
  All_type_meta_pos_subset %>%
  left_join(All_type_meta_cn15[c('cell_id', 'CN')], "cell_id") %>%
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
res_pa_df$group <- factor(res_pa_df$group, levels = unique(sample_map_group[,2]))

res_pa_df2  <-
  res_pa_df2 %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df2$CN <- factor(res_pa_df2$CN, levels = paste0("CN", 1:10))
res_pa_df2$group <- factor(res_pa_df2$group, levels = unique(sample_map_group[,2]))


res_pa_df3  <-
  res_pa_df3 %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df3$CN <- factor(res_pa_df3$CN, levels = paste0("CN", 1:15))
res_pa_df3$group <- factor(res_pa_df3$group, levels = unique(sample_map_group[,2]))


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

# export plot
pdf(
  file.path(All_type_Proximity_out_dir, "Proximity analysis.pdf"),
  height = 25,
  width  = 12
)
p_res_pa%>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "Proximity analysis - CN10.pdf"),
  height = 25,
  width  = 12
)
p_res_pa_ls_cn10 %>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "Proximity analysis - CN15.pdf"),
  height = 25,
  width  = 12
)
p_res_pa_ls_cn15 %>% print()
dev.off()








### 10 50 100 150 200 -------
res_pa_ls <- list()

res_pa_ls[['10']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 10, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
res_pa_ls[['50']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 50, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['100']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 100, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['150']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 150, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

res_pa_ls[['200']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", #NULL: 提前返回window的细胞近邻的list
    k = 1200, # 'auto',
    searchtype = "radius",
    radius = 200, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

saveRDS(res_pa_ls, file.path(All_type_Proximity_out_dir, "res_pa_ls.rds"))




proximity_cells <- c(
  "Cardiomyocytes", "FB", "MAC", "EC"
  # "GME"
)

proximity_all_data_ls <- 
  proximity_cells %>%
  setNames(., .) %>%
  map(function(cell){
    cell_idx <- All_type_meta_pos_subset$Classification1 == cell
    stopifnot(sum(cell_idx) > 0)
    
    # 1 all_data
    all_data_sub  <- All_type_meta_pos_subset[cell_idx, ]
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


p_proximity_cell_count_distribution_ls <- 
  proximity_all_data_ls %>%
  imap(function(dat_ls, name){
    
    all_data_group_sub <- dat_ls$all_data
    
    p <-
      dat_ls$res_pa %>%
      pivot_longer(-c(Slide_ID, distance), names_to = "Cell", values_to = "count") %>%
      left_join(all_data_group_sub, by = "Slide_ID")  %>%
      mutate(distance = as.numeric(distance)) %>%
      
      ggplot(aes(x = Cell, y = count, fill = Slide_ID)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ distance, scales = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = plan) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    p 
  })


pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - barplot -1.pdf"), width = 12, height = 30)
p_proximity_cell_count_distribution_ls %>% walk(print)
dev.off()


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
      facet_wrap(~ Slide_ID, scales = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = plan) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    p 
  })

pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - barplot -2.pdf"), width = 12, height = 30)
p_proximity_cell_count_distribution_ls2 %>% walk(print)
dev.off()



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


# # show
# head(res_pa_epi_count_tidy, 10)

p_proximity_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    # p <-
    #   res_pa_count_tidy %>%
    #   mutate(distance = as.numeric(distance)) %>%
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
      group_by(Cell, distance, Slide_ID) %>%
      summarise(count = median(count)) %>%
      
      mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
      
      ggplot(aes(x = distance, y = count, color = Slide_ID)) +
      # geom_density(alpha = 0.6) +
      geom_smooth(method = 'loess') +
      facet_wrap(~ Cell, scales = "free_y") +
      ggpubr::theme_pubr() +
      theme(
        strip.background = element_blank(),
        panel.grid = element_blank()
      ) +
      scale_color_manual(values = plan) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })

# export pdf
pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - smooth_lines.pdf"), width = 18, height = 12)
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
      group_by(Slide_ID, Cell, distance) %>%
      summarise(count = median(count)) %>%
      ggplot(aes(x = Cell, y = distance, size = count, color = Slide_ID)) +
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
      scale_color_manual(values = plan)+
      coord_polar(theta = 'x', clip = 'off') +
      labs(x = NULL, y = NULL, title = name)
    p 
  })


# export pdf
pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - radarplot.pdf"), width = 16, height = 12)
p_proximity_radar_ls %>% walk(print)
dev.off()











# spatial analysis cell_type_sub--------
library(Seurat)
library(tidyverse)

source("spatial analysis.R")


##  all_data ------
CN_out_dir <- "ROI_analysis/spatial_analysis_all_cell/CN analysis"
dir.create(CN_out_dir, recursive = T)

all_data <- readRDS("data_for_analysis/all_data_roi_1208.rds")
sc_fibro <- readRDS("fibroblast_subcluster/modify_0918/sc_fibro_0918.rds")


library(tibble)
library(dplyr)

meta_data_fibro <- 
  sc_fibro[["annotation_fb_0918"]] %>% 
  rownames_to_column('cells') %>% 
  right_join(
    all_data[["cell_type_res05"]] %>% rownames_to_column("cells"), 
    by = "cells"
  ) %>% 
  mutate(
    cell_type_res05_sub = ifelse((as.character(cell_type_res05) == "FB"), as.character(annotation_fb_0918), 
                                 as.character(cell_type_res05))
  ) %>% 
  dplyr::select(cells, cell_type_res05_sub) %>% 
  column_to_rownames('cells')

all_data <- AddMetaData(all_data, meta_data_fibro['cell_type_res05_sub'])
saveRDS(all_data, "data_for_analysis/all_data_roi_fbsub_1218.rds")


CN_data <- all_data

Idents(CN_data) <- "cell_type_res05_sub"

cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()

## CN analysis ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

all_data_submeta_for_cn_analysis <- 
  CN_data[[c("project", "group", "cell_type_res05_sub")]] %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    cell_type_res05_sub = factor(cell_type_res05_sub), # import!
    obj_id = cell_id # import!
  )

### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "cell_type_res05_sub",
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
  group.by       = "project",
  classification = "cell_type_res05_sub",
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
# # CN k10
# cn_annotation <- c(
#   "CN1"  = "Macrophage",
#   "CN2"  = "T",
#   "CN3"  = "B",
#   "CN4"  = "Fibro_Stroma",
#   "CN5"  = "Myo_Stroma",
#   "CN6"  = "Pan immune",
#   "CN7"  = "NK_Endothelial_Stroma",
#   "CN8"  = "Epithelial",
#   "CN9"  = "T_Epithelial",
#   "CN10" = "T_mast_Epithelial"
# )
# 
# p_cn_heat2 <-
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
#   file.path(CN_out_dir, "CN_heatmap_annotation.pdf"),
#   height = 6,
#   width  = 11
# )
# cowplot::plot_grid(p_cn_heat2$gtable)
# dev.off()


### CN10 spatial plot -----------------
cn10_sp_out_dir <- file.path(CN_out_dir, 'CN10 spatial plot') 
dir.create(cn10_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn10 <- 
  CN_data[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn10, file.path(cn10_sp_out_dir, "all_data_meta_cn10.rds"))

cn10_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn10[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
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

### cell props ------
All_type_cn10_count_df <-
  All_type_meta_cn10 %>% 
  group_by(project, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

All_type_cn10_count_df$project <- factor(All_type_cn10_count_df$project, levels = sample_map_group[,1])
All_type_cn10_count_df$group <- factor(All_type_cn10_count_df$group, levels = unique(sample_map_group[,2]))

p_cn10_props <-
  All_type_cn10_count_df %>%
  
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

library(ggpubr)
p_cn10_props_group <-
  All_type_cn10_count_df %>%
  
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

p_cn10_props_group_2 <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif",size=3,aes(group = group))+
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn10_props %>% print
p_cn10_props_group %>% print
p_cn10_props_group_2 %>% print
dev.off()



### CN15 spatial plot -----------------
cn15_sp_out_dir <- file.path(CN_out_dir, 'CN15 spatial plot') 
dir.create(cn15_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn15 <- 
  CN_data[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn15, file.path(cn15_sp_out_dir, "all_data_meta_cn15.rds"))

cn15_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn15[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
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
  group_by(project, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

All_type_cn15_count_df$project <- factor(All_type_cn15_count_df$project, levels = sample_map_group[,1])
All_type_cn15_count_df$group <- factor(All_type_cn15_count_df$group, levels = unique(sample_map_group[,2]))

p_cn15_props <-
  All_type_cn15_count_df %>%
  
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

p_cn15_props_group_2 <-
  All_type_cn15_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif",size=3,aes(group = group))+
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
p_cn15_props_group_2 %>% print
dev.off()








## Proximity analysis -------
### 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

All_type_Proximity_out_dir <- "ROI_analysis/spatial_analysis_all_cell/Proximity analysis"
dir.create(All_type_Proximity_out_dir)

# all_data_cn_meta
All_type_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){#browser()
    
    proj_meta_sub <-
      CN_data[[]] %>%  
      dplyr::select(cell_type_res05_sub, group, project) %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%
      dplyr::rename(Classification1 = cell_type_res05_sub,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = project)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

All_type_meta_pos_subset <- All_type_meta_pos %>% rownames_to_column("cell_id")

saveRDS(All_type_meta_pos_subset, file = file.path(All_type_Proximity_out_dir, "Subtype_meta_pos_subset.rds"))

All_type_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
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
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = "auto", # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

saveRDS(res_pa, file = file.path(All_type_Proximity_out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS(file.path(All_type_Proximity_out_dir, "res_pa.rds"))


# summarized by group and celltype
pa_groups <-
  with(
    All_type_meta_pos_subset,
    str_c(Classification1, group, Slide_ID, sep = '..')
  )

pa_groups2 <-
  All_type_meta_pos_subset %>%
  left_join(All_type_meta_cn10[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

pa_groups3 <-
  All_type_meta_pos_subset %>%
  left_join(All_type_meta_cn15[c('cell_id', 'CN')], "cell_id") %>%
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
res_pa_df$group <- factor(res_pa_df$group, levels = unique(sample_map_group[,2]))

res_pa_df2  <-
  res_pa_df2 %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df2$CN <- factor(res_pa_df2$CN, levels = paste0("CN", 1:10))
res_pa_df2$group <- factor(res_pa_df2$group, levels = unique(sample_map_group[,2]))


res_pa_df3  <-
  res_pa_df3 %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df3$CN <- factor(res_pa_df3$CN, levels = paste0("CN", 1:15))
res_pa_df3$group <- factor(res_pa_df3$group, levels = unique(sample_map_group[,2]))


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

# export plot
pdf(
  file.path(All_type_Proximity_out_dir, "Proximity analysis.pdf"),
  height = 25,
  width  = 12
)
p_res_pa%>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "Proximity analysis - CN10.pdf"),
  height = 25,
  width  = 12
)
p_res_pa_ls_cn10 %>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "Proximity analysis - CN15.pdf"),
  height = 25,
  width  = 12
)
p_res_pa_ls_cn15 %>% print()
dev.off()








### 10 50 100 150 200 -------
res_pa_ls <- list()

res_pa_ls[['10']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 10, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
res_pa_ls[['50']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 50, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['100']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 100, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['150']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 150, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

res_pa_ls[['200']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", #NULL: 提前返回window的细胞近邻的list
    k = 1200, # 'auto',
    searchtype = "radius",
    radius = 200, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

saveRDS(res_pa_ls, file.path(All_type_Proximity_out_dir, "res_pa_ls.rds"))




proximity_cells <- c(
  "Cardiomyocytes", "f1", "MAC", "EC"
  # "GME"
)

proximity_all_data_ls <- 
  proximity_cells %>%
  setNames(., .) %>%
  map(function(cell){
    cell_idx <- All_type_meta_pos_subset$Classification1 == cell
    stopifnot(sum(cell_idx) > 0)
    
    # 1 all_data
    all_data_sub  <- All_type_meta_pos_subset[cell_idx, ]
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
      facet_wrap(~ distance, scales = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = plan) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    p 
  })


pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - barplot -1.pdf"), width = 12, height = 30)
p_proximity_cell_count_distribution_ls %>% walk(print)
dev.off()


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
      facet_wrap(~ Slide_ID, scales = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = plan) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    p 
  })

pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - barplot -2.pdf"), width = 12, height = 30)
p_proximity_cell_count_distribution_ls2 %>% walk(print)
dev.off()



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


# # show
# head(res_pa_epi_count_tidy, 10)

p_proximity_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    # p <-
    #   res_pa_count_tidy %>%
    #   mutate(distance = as.numeric(distance)) %>%
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
      scale_color_manual(values = plan) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })

# export pdf
pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - smooth_lines.pdf"), width = 24, height = 16)
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
      scale_color_manual(values = plan)+
      coord_polar(theta = 'x', clip = 'off') +
      labs(x = NULL, y = NULL, title = name)
    p 
  })


# export pdf
pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - radarplot.pdf"), width = 16, height = 12)
p_proximity_radar_ls %>% walk(print)
dev.off()













# spatial analysis cell_type_sub--------
library(Seurat)
library(tidyverse)

source("spatial analysis.R")


##  all_data ------
CN_out_dir <- "ROI_analysis/spatial_analysis_roi/CN analysis"
dir.create(CN_out_dir, recursive = T)

all_data <- readRDS("data_for_analysis/all_data_roi_fbsub_1218.rds")
all_data_filt <- subset(all_data, subset = ROI != "Other")

dim(all_data)
dim(all_data_filt)

CN_data <- all_data_filt

Idents(CN_data) <- "cell_type_res05_sub"

cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()

## CN analysis ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

all_data_submeta_for_cn_analysis <- 
  CN_data[[c("project", "group", "cell_type_res05_sub")]] %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    cell_type_res05_sub = factor(cell_type_res05_sub), # import!
    obj_id = cell_id # import!
  )

### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "cell_type_res05_sub",
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
  group.by       = "project",
  classification = "cell_type_res05_sub",
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
# # CN k10
# cn_annotation <- c(
#   "CN1"  = "Macrophage",
#   "CN2"  = "T",
#   "CN3"  = "B",
#   "CN4"  = "Fibro_Stroma",
#   "CN5"  = "Myo_Stroma",
#   "CN6"  = "Pan immune",
#   "CN7"  = "NK_Endothelial_Stroma",
#   "CN8"  = "Epithelial",
#   "CN9"  = "T_Epithelial",
#   "CN10" = "T_mast_Epithelial"
# )
# 
# p_cn_heat2 <-
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
#   file.path(CN_out_dir, "CN_heatmap_annotation.pdf"),
#   height = 6,
#   width  = 11
# )
# cowplot::plot_grid(p_cn_heat2$gtable)
# dev.off()


### CN10 spatial plot -----------------
cn10_sp_out_dir <- file.path(CN_out_dir, 'CN10 spatial plot') 
dir.create(cn10_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn10 <- 
  CN_data[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn10, file.path(cn10_sp_out_dir, "all_data_meta_cn10.rds"))

cn10_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn10[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
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

### cell props ------
All_type_cn10_count_df <-
  All_type_meta_cn10 %>% 
  group_by(project, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

All_type_cn10_count_df$project <- factor(All_type_cn10_count_df$project, levels = sample_map_group[,1])
All_type_cn10_count_df$group <- factor(All_type_cn10_count_df$group, levels = unique(sample_map_group[,2]))

p_cn10_props <-
  All_type_cn10_count_df %>%
  
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

library(ggpubr)
p_cn10_props_group <-
  All_type_cn10_count_df %>%
  
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

p_cn10_props_group_2 <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif",size=3,aes(group = group))+
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn10_props %>% print
p_cn10_props_group %>% print
p_cn10_props_group_2 %>% print
dev.off()



### CN15 spatial plot -----------------
cn15_sp_out_dir <- file.path(CN_out_dir, 'CN15 spatial plot') 
dir.create(cn15_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn15 <- 
  CN_data[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn15, file.path(cn15_sp_out_dir, "all_data_meta_cn15.rds"))

cn15_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn15[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
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
  group_by(project, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

All_type_cn15_count_df$project <- factor(All_type_cn15_count_df$project, levels = sample_map_group[,1])
All_type_cn15_count_df$group <- factor(All_type_cn15_count_df$group, levels = unique(sample_map_group[,2]))

p_cn15_props <-
  All_type_cn15_count_df %>%
  
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

p_cn15_props_group_2 <-
  All_type_cn15_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  stat_compare_means(label = "p.signif",size=3,aes(group = group))+
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
p_cn15_props_group_2 %>% print
dev.off()








## Proximity analysis -------
### 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

All_type_Proximity_out_dir <- "ROI_analysis/spatial_analysis_roi/Proximity analysis"
dir.create(All_type_Proximity_out_dir)

# all_data_cn_meta
All_type_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){#browser()
    
    proj_meta_sub <-
      CN_data[[]] %>%  
      dplyr::select(cell_type_res05_sub, group, project) %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%
      dplyr::rename(Classification1 = cell_type_res05_sub,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = project)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

All_type_meta_pos_subset <- All_type_meta_pos %>% rownames_to_column("cell_id")

saveRDS(All_type_meta_pos_subset, file = file.path(All_type_Proximity_out_dir, "Subtype_meta_pos_subset.rds"))

All_type_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
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
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = "auto", # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

saveRDS(res_pa, file = file.path(All_type_Proximity_out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS(file.path(All_type_Proximity_out_dir, "res_pa.rds"))


# summarized by group and celltype
pa_groups <-
  with(
    All_type_meta_pos_subset,
    str_c(Classification1, group, Slide_ID, sep = '..')
  )

pa_groups2 <-
  All_type_meta_pos_subset %>%
  left_join(All_type_meta_cn10[c('cell_id', 'CN')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN, sep = '..')
  )

pa_groups3 <-
  All_type_meta_pos_subset %>%
  left_join(All_type_meta_cn15[c('cell_id', 'CN')], "cell_id") %>%
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
res_pa_df$group <- factor(res_pa_df$group, levels = unique(sample_map_group[,2]))

res_pa_df2  <-
  res_pa_df2 %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df2$CN <- factor(res_pa_df2$CN, levels = paste0("CN", 1:10))
res_pa_df2$group <- factor(res_pa_df2$group, levels = unique(sample_map_group[,2]))


res_pa_df3  <-
  res_pa_df3 %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "group", "project", "CN"), sep = "\\.\\.") 
res_pa_df3$CN <- factor(res_pa_df3$CN, levels = paste0("CN", 1:15))
res_pa_df3$group <- factor(res_pa_df3$group, levels = unique(sample_map_group[,2]))


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

# export plot
pdf(
  file.path(All_type_Proximity_out_dir, "Proximity analysis.pdf"),
  height = 25,
  width  = 12
)
p_res_pa%>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "Proximity analysis - CN10.pdf"),
  height = 25,
  width  = 12
)
p_res_pa_ls_cn10 %>% print()
dev.off()


pdf(
  file.path(All_type_Proximity_out_dir, "Proximity analysis - CN15.pdf"),
  height = 25,
  width  = 12
)
p_res_pa_ls_cn15 %>% print()
dev.off()








### 10 50 100 150 200 -------
res_pa_ls <- list()

res_pa_ls[['10']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 10, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
res_pa_ls[['50']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 50, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['100']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 100, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )
res_pa_ls[['150']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1" , #NULL: 提前返回window的细胞近邻的list
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 150, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

res_pa_ls[['200']] <-
  ProximutyAnalysis_NN_multi(
    All_type_meta_pos_subset,
    query = All_type_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", #NULL: 提前返回window的细胞近邻的list
    k = 1200, # 'auto',
    searchtype = "radius",
    radius = 200, # <--------
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

saveRDS(res_pa_ls, file.path(All_type_Proximity_out_dir, "res_pa_ls.rds"))




proximity_cells <- c(
  "Cardiomyocytes", "f1", "MAC", "EC"
  # "GME"
)

proximity_all_data_ls <- 
  proximity_cells %>%
  setNames(., .) %>%
  map(function(cell){
    cell_idx <- All_type_meta_pos_subset$Classification1 == cell
    stopifnot(sum(cell_idx) > 0)
    
    # 1 all_data
    all_data_sub  <- All_type_meta_pos_subset[cell_idx, ]
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
      facet_wrap(~ distance, scales = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = plan) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    p 
  })


pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - barplot -1.pdf"), width = 12, height = 30)
p_proximity_cell_count_distribution_ls %>% walk(print)
dev.off()


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
      facet_wrap(~ Slide_ID, scales = "free_y", ncol = 1) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = plan) +
      labs(
        # x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    p 
  })

pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - barplot -2.pdf"), width = 12, height = 30)
p_proximity_cell_count_distribution_ls2 %>% walk(print)
dev.off()



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


# # show
# head(res_pa_epi_count_tidy, 10)

p_proximity_ls <-
  res_pa_count_tidy_ls %>%
  imap(function(res_pa_count_tidy, name){
    # p <-
    #   res_pa_count_tidy %>%
    #   mutate(distance = as.numeric(distance)) %>%
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
      scale_color_manual(values = plan) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })

# export pdf
pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - smooth_lines.pdf"), width = 24, height = 16)
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
      scale_color_manual(values = plan)+
      coord_polar(theta = 'x', clip = 'off') +
      labs(x = NULL, y = NULL, title = name)
    p 
  })


# export pdf
pdf(file.path(All_type_Proximity_out_dir, "Proximity analysis - radarplot.pdf"), width = 16, height = 12)
p_proximity_radar_ls %>% walk(print)
dev.off()







