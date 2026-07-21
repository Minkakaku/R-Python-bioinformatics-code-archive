library(Seurat)
library(tidyverse)

source("spatial analysis.R")

out_dir <- "CN_analysis"
dir.create(out_dir)

all_data <- readRDS("all_data.rds")
Idents(all_data) <- "annotation_full"

cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()



# implement CN ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

all_data_meta_for_cn_analysis <- 
  all_data[[c("project", "group", "annotation_full")]] %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    annotation_full = factor(annotation_full), # import!
    obj_id = cell_id # import!
  )

res_cn <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation_full",
  k.nn           = 50,
  k.cluster      = 10
)


# saveRDS(res_cn, file.path(out_dir, "res_cn.rds"))

# heatmap ----------------------
# FoldChange cluster/allCell
cn_cluster_mean <- res_cn$cluster_centers
cn_baseMean     <- res_cn$base_mean

# log2(FC)
cn_cluster_FC <- log2(1 + sweep(cn_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn_baseMean))

# heatmap
p_cn_heat <-
  pheatmap::pheatmap(
    cn_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(out_dir, "CN_heatmap_CN.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn_heat$gtable)
dev.off()

# CN annotaton --------------------
# CN k10


cn_annotation <- c(
  "CN1"  = "Macrophage",
  "CN2"  = "T",
  "CN3"  = "B",
  "CN4"  = "Fibro_Stroma",
  "CN5"  = "Myo_Stroma",
  "CN6"  = "Pan immune",
  "CN7"  = "NK_Endothelial_Stroma",
  "CN8"  = "Epithelial",
  "CN9"  = "T_Epithelial",
  "CN10" = "T_mast_Epithelial"
)

p_cn_heat2 <-
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
  file.path(out_dir, "CN_heatmap_annotation.pdf"),
  height = 6,
  width  = 11
)
cowplot::plot_grid(p_cn_heat2$gtable)
dev.off()


# CN spatial plot -----------------
cn_sp_out_dir <- file.path(out_dir, 'cn spatial plot') 
dir.create(cn_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

all_data_meta <- 
  all_data[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))

saveRDS(all_data_meta, 'all_data_meta_cn.rds')

cn_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]

    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN2, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })


# export image
cn_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

# cell props ------
all_cn_count_df <-
  all_data_meta %>% 
  group_by(project, group, CN2) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn_props <-
  all_cn_count_df %>%
  
  ggplot(aes(x = project, y = freq, fill = CN2)) +
  geom_col() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn_props_group <-
  all_cn_count_df %>%
  
  ggplot(aes(x = CN2, y = freq, fill = group)) +
  geom_boxplot() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(out_dir, "CN_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn_props %>% print
p_cn_props_group %>% print
dev.off()



# Proximity analysis -------
# 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){#browser()
    
    proj_meta_sub <-
      all_data[[]] %>%  
      dplyr::select(annotation_full, group, project) %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%
      dplyr::rename(Classification1 = annotation_full,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = project)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos_subset <- all_data_meta_pos %>% rownames_to_column("cell_id")
# saveRDS(all_data_meta_pos_subset,'./data/all_data_meta_pos_subset.rds')

all_data_meta_pos_subset %>% head

#        cell_id Classification1   Group           Slide_ID        X        Y
# 1 3_aaaacdgo-1     Fibroblasts   MUT    20231128__R01-001 2228.621 860.3357
# 2 3_aaaanalm-1   Myeloid cells   MUT    20231128__R01-001 2248.922 850.1516
# 3 3_aaaaonpl-1   Myeloid cells   MUT    20231128__R01-001 2270.447 856.8892
# 4 3_aaabciio-1         B cells   MUT    20231128__R01-001 2279.562 838.6568
# 5 3_aaabdkea-1     Fibroblasts   MUT    20231128__R01-001 2292.319 858.3909
# 6 3_aaacdgip-1         B cells   MUT    20231128__R01-001 2288.924 868.3226
# ......


# 2. proximity analysis -------
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

## 2.1 proximity analysis --------
res_pa <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset,
    query = all_data_meta_pos_subset,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    k = 750, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )

# saveRDS(res_pa, file = file.path(out_dir, "res_pa.rds"))

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
    all_data_meta_pos_subset,
    str_c(Classification1, group, Slide_ID, sep = '..')
  )

pa_groups2 <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta[c('cell_id', 'CN2')], "cell_id") %>%
  with(
    str_c(Classification1, group, Slide_ID, CN2, sep = '..')
  )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
res_pa_df2 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups2, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 


res_pa_df2  <-
  res_pa_df2 %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project", "CN2"), sep = "\\.\\.") 


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
    # dplyr::filter(Target == "Epithelial cells") %>%
    ggplot(aes(x = Others, y = cell_props, fill = Group)) +
    geom_boxplot(width = 0.45, outlier.colour = NA) +
    # geom_jitter(width = 0.25) +
    facet_wrap(~ Target, scales = "free_y", ncol = 2)+
    ggpubr::theme_pubclean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 

p_res_pa_ls <-
  res_pa_df2 %>%
  with(split(res_pa_df2, CN2)) %>%
  imap(function(dat, name){
    dat %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle(name)
  })

# export plot
pdf(
  file.path(out_dir, "proximity analysis.pdf"),
  height = 14,
  width  = 12
)
p_res_pa%>% print()
dev.off()


pdf(
  file.path(out_dir, "proximity analysis - CN.pdf"),
  height = 14,
  width  = 12
)
p_res_pa_ls %>% print()
dev.off()


