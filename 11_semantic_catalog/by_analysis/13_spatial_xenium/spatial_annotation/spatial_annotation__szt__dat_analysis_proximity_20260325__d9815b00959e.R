library(Seurat)
library(tidyverse)

source("spatial analysis.R")

out_dir_260325 <- "20260325 dat_analysis"
dir.create(out_dir_260325)


all_data_full <- readRDS("all_data_full_v6.rds")

# many cells are droped by QC
# use all_data instead
# all_data_full <- merge(visim_ls[[1]], visim_ls[-1])

cell_relevel <- c(
  'B cells'          = 'B',
  'IgG plasma cells' = 'Plasma',
  'IgA plasma cells' = 'Plasma',
  
  'tumor_cells'           = 'Epithelial',
  'enterocytes'           = 'Epithelial',
  'enteroendocrine_cells' = 'Epithelial',
  'pmc_infla'             = 'Epithelial',
  'gmc'                   = 'Epithelial',
  'goblet_cells'          = 'Epithelial',
  'pmc'                   = 'Epithelial',
  'endocrine_cells'       = 'Epithelial',
  'parietal_cells'        = 'Epithelial',
  'chief_cells'           = 'Epithelial',
  
  'EC'        = 'EC',
  'Pericyte'  = 'Pericyte',
  'CAF'       = 'CAF',
  'SMC'       = 'SMC',
  
  'Neutrophils' = 'Neutrophil',
  'Macrophage'  = 'Macrophage',
  'Monocytes'   = 'Monocyte',
  'Mast Cells'  = 'Mast Cell',
  'DC'          = 'DC',
  
  'CD4T'                = 'T',
  'CD8T'                = 'T',
  'NK'                  = 'T',
  'gamma-delta T cells' = 'T',
  
  'gial_cells'          = 'gial_cell'
)
all_data_full$celltype  <- 
  all_data_full$sub_celltype2 %>%
    fct_recode(!!!structure(
      names(cell_relevel),
      names = cell_relevel
    ))


table(is.na(all_data_full$celltype))


Idents(all_data_full) <- "celltype"

cells_position_ls <- 
  all_data_full@images %>%
  map(GetTissueCoordinates) %>%
  map(dplyr::rename, cell_id = cell) %>% 
  setNames(., str_remove(names(.), "\\..*$"))  
# map(function(x){
#   x[x$cell_id %in% all_cellid_not_others, ]
# })
cells_position_df <- cells_position_ls %>% bind_rows() 


all_data_meta_for_cn_analysis <- 
  all_data_full[[]] %>%
  dplyr::filter(celltype != "Others") %>% 
  # dplyr::filter(group_stage != "late_cancer") %>% # 20250923 remove late_cancer
  dplyr::filter(!is.na(celltype)) %>% 
  dplyr::mutate(annotation_full = celltype ) %>% 
  
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x', 'y')], by = "cell_id") %>%
  rename(X = x, Y = y) %>%
  mutate(
    annotation_full = factor(annotation_full), # important!
    obj_id = cell_id # important!
  )

stopifnot(all(!is.na(all_data_meta_for_cn_analysis$celltype)))



# Proximity analysis -------
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



# 1. input data ------------
# cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

all_data_meta <- 
  # readRDS(file.path(cn_out_dir, 'all_data_meta_cn2.rds')) %>% 
  # dplyr::mutate(annotation_full = celltype )
  all_data_meta_for_cn_analysis

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    # all_data_meta <- 
    #   all_data_full[[]] %>%
    #   dplyr::filter(celltype != "Others") %>%
    #   dplyr::mutate(annotation_full = celltype )
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(annotation_full, cell_id, group, project) %>%
      # rownames_to_column("cell_id") %>%
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


# cell_id         X         Y Classification1       group            Slide_ID
# 1 A210239_s_008um_00260_00253-1  7446.846  7798.865    stromal_cell      others      others_A210239
# 2 A210239_s_008um_00266_00304-1  8193.100  7878.610    stromal_cell      others      others_A210239
# 3 A210239_s_008um_00154_00449-1 10294.610  6219.288    stromal_cell      others      others_A210239
# 4 A210239_s_008um_00226_00670-1 13535.346  7237.088    stromal_cell      others      others_A210239
# 5 A210239_s_008um_00446_00748-1 14709.489 10439.890    plasma cells      others      others_A210239
# 6 A210239_s_008um_00235_00397-1  9547.348  7411.098     tumor_cells late_cancer late_cancer_A210239
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
    k = 300, # 'auto',
    searchtype = "radius",
    radius = 50, #  50um半径
    mpp = 8 / all_data_full@images[[1]]@scale.factors$spot, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )
# saveRDS(res_pa, file = file.path(out_dir_260325, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}

# res_pa <- readRDS(file.path(out_dir_260325, 'res_pa.rds'))


# summarized by group and celltype
pa_groups <-
  all_data_meta_pos_subset %>%
  group_by(Slide_ID) %>%
  nest() %>%
  deframe() %>%
  imap(function(dat, name){
    with(
      dat,
      str_c(Classification1, group, name, sep = '..')
    )
  }) %>%
  unlist(use.names = FALSE)


# pa_groups2 <-
#   all_data_meta_pos_subset %>%
#   left_join(all_data_meta[c('cell_id', 'CN2')], "cell_id") %>%
#   with(
#     str_c(Classification1, group, Slide_ID, CN2, sep = '..')
#   )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups, FUN = "median")})
# res_pa_df2 <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups2, FUN = "median")})

res_pa_df  <-
  res_pa_df %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 

# adjust gorup levels
res_pa_df$Group <-
  res_pa_df$Group %>% 
  fct_recode(
    NAG = "1NAG",
    # "CAG"
    # "early_cancer"
    # "HGD"
    # "IM"
    # "late_cancer"
    # "LGD"         
    # "others"
    others = "SMC",
    others = "TLS"  
  ) %>%
  fct_relevel(
    'NAG', 'CAG', 'IM', 'LGD', 'HGD',
    'early_cancer', 'late_cancer', 'others'
  ) 

# proximity props / slide props
all_props_df <-
  all_data_meta_pos_subset %>%
  group_by(Slide_ID, Classification1) %>%
  summarise(count = n()) %>%
  group_by(Slide_ID) %>%
  mutate(freq = count/sum(count))


res_pa_df2 <-
  res_pa_df %>%
  left_join(
    all_props_df[c('Slide_ID', 'Classification1', 'freq')],
    by = c(Others = 'Classification1', project = 'Slide_ID')
  ) %>%
  mutate(cell_props2 = cell_props / freq)



# res_pa_df2  <-
#   res_pa_df2 %>%
#   as.data.frame() %>%
#   rownames_to_column('group') %>%
#   pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
#   separate(group, into = c("Target", "Group", "project", "CN2"), sep = "\\.\\.") 


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa <-
  res_pa_df %>%
  # mutate(Group = factor(Group, levels = c(
  #   "NORMAL", "IM",  "LGD", "HGD",  "early_cancer", "late_cancer", "TLS", "other"        
  # ))) %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA, linewidth = 0.1) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_res_pa2 <-
  res_pa_df2 %>%
  # mutate(Group = factor(Group, levels = c(
  #   "NORMAL", "IM",  "LGD", "HGD",  "early_cancer", "late_cancer", "TLS", "other"        
  # ))) %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props2, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA, linewidth = 0.1) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  scale_y_log10() +
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


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
  file.path(out_dir_260325, "proximity analysis.pdf"),
  height = 15,
  width  = 15
)
p_res_pa %>% print()
# p_res_pa2 %>% print()
dev.off()

pdf(
  file.path(out_dir_260325, "proximity analysis2.pdf"),
  height = 15,
  width  = 15
)
p_res_pa2 %>% print()
dev.off()
# 
# pdf(
#   file.path(out_dir_260325, "proximity analysis - CN.pdf"),
#   height = 50,
#   width  = 25
# )
# p_res_pa_ls %>% print()
# dev.off()
res_pa_df %>%
  dplyr::filter(
    Target == "Epithelial",
    Others %in% c("CAF", "Macrophage", "T"),
    Group != "others"
  )  %>%
  write.csv(file.path(out_dir_260325, 'res_pa.tidy.csv'))


res_pa_df %>%
  dplyr::filter(
    Target == "Epithelial",
    Others %in% c("CAF", "Macrophage", "T"),
    Group != "others"
  ) %>%
  ggplot(aes(x = Group, y = cell_props, fill = Others, color = Others)) +
  geom_boxplot(outliers = FALSE, alpha = 0.5,  width = 0.5, color = "black") +
  geom_jitter(width = 0.15, height = 0, size = 1) +
  facet_wrap(~Others) +
  theme_classic() +
  Seurat::RotatedAxis() +
  theme(
    plot.margin = margin(1,1,1,8, "mm"),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12)
  ) +
  ggsci::scale_color_igv() +
  ggsci::scale_fill_igv() +
  labs(
    x = NULL,
    y = NULL,
    # title = name
  )

cell_cols_for_3_cells <- 
  structure(
    names = c("CAF", "Macrophage", "T"),
    ggsci::pal_igv()(3)
  )

p_epi_boxplot_ls <-
  res_pa_df %>%
  dplyr::filter(
    Target == "Epithelial",
    Others %in% c("CAF", "Macrophage", "T"),
    Group != "others"
  ) %>%
  # mutate(Group = factor(Group)) %>%
  group_by(Others) %>%
  nest() %>%
  deframe() %>%
  imap(function(dat, name){ # browser()
    
    signif_width <- 0.25
    signif_text_size <- 2.5
    signif_y_fold <- 1.01
    signif_line_size <- 0.3
    
    # beta
    plot_data_beta <- dat
    
    # !!!important  
    plot_data_beta$Group <- plot_data_beta$Group %>% factor(., levels = sort(unique(.)))
    
    signif_df_1 <-
      ggpubr::compare_means(
        as.formula(str_glue("cell_props ~ Group")),
        data = plot_data_beta,
        method = "wilcox"
      ) %>% # print()
      dplyr::filter(p < 0.05) %>% 
      
      mutate(
        group1 = factor(group1, levels = levels(plot_data_beta$Group))%>% as.numeric(),
        group2 = factor(group2, levels = levels(plot_data_beta$Group))%>% as.numeric(),
        xmin = group1,
        xmax = group2,
        ymax = 1:n()
      )
    ymax1 <- plot_data_beta$cell_props %>% max %>% ceiling() 
    p1 <-
      plot_data_beta %>%
      ggplot(aes(x = Group, y = cell_props)) +
      geom_boxplot(outliers = FALSE, alpha = 0.5,  width = 0.5, color = "black", fill = cell_cols_for_3_cells[name]) +
      geom_jitter(width = 0.15, height = 0, size = 1, color = cell_cols_for_3_cells[name]) +
      ggsignif::geom_signif(
        y_position = signif_y_fold * ymax1 + ymax1 * signif_df_1$ymax * 0.1,
        xmin = signif_df_1$xmin,
        xmax = signif_df_1$xmax,
        annotations = signif_df_1$p.format,
        tip_length = 0,
        size = signif_line_size,
        textsize = signif_text_size
      ) +
      theme_classic() +
      Seurat::RotatedAxis() +
      theme(
        plot.margin = margin(1,1,1,8, "mm"),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12)
      ) +
      labs(
        x = NULL,
        y = NULL,
        title = name
      )
    p1
  })
  


pdf(
  file.path(out_dir_260325, "proximity analysis - epithleial - v2.pdf"),
  height = 4,
  width  = 13
)
p_epi_boxplot_ls %>% patchwork::wrap_plots(nrow = 1) %>% print()
dev.off()

## 2.2 proximity analysis - 200um --------------- 
all_data_meta_pos_subset_epi <- all_data_meta_pos_subset %>% dplyr::filter(Classification1 == "Epithelial")
res_pa_200um <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset,
    query = all_data_meta_pos_subset_epi,
    x = "X",
    y = "Y",
    group.by = "Slide_ID",
    # classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
    classification = NULL,
    k = 1000, # 'auto',
    searchtype = "radius",
    radius = 200, # 邻域半径
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 90
  )


saveRDS(res_pa_200um, file = file.path(out_dir_260325, str_glue("epi_res_pa_200um.rds")))


cell_match_var <- structure(
  as.character(all_data_meta_pos_subset$Classification1),
  names = all_data_meta_pos_subset$cell_id
)


res_nn_df_ls <-
  res_pa_200um %>%
  # .[1] %>%
  imap(function(nn, grp){ # browser()
    print(grp)
    print(Sys.time())
    res_ls <- list()
    
    nn.idx  <- nn$nn.idx[, -1] + 1
    nn.dist <- nn$nn.dists[, -1]
    
    
    cell_match_var <- factor( c("Others", cell_match_var) )
    
    for(i in 1:nrow(nn$nn.idx)){
      if(i %% 5000 == 0) print(i)
      
      # head 10 cells, mean dist
      res <- tapply(nn.dist[i, ], cell_match_var[nn.idx[i, ]], function(x){mean(head(x, 10))})
      res_ls[[i]] <- res
    }
    
    print(Sys.time())
    
    res_df <- as.data.frame(bind_rows(res_ls))
    rownames(res_df) <- nn$query_cell_id
    res_df
    # res_ls
    
  })


# cell count prop 
res_nn_count_df_ls <-
  res_nn_df_ls %>%
  map(function(x){ # browser()
    xx <- x[, !colnames(x) %in% c("Others")]   # , "Epithelial")]
    xx[is.na(xx)] <- Inf
    
    dat_ls <-
      apply(xx, 2, function(y){
        x_cut <-
          cut(
            y, 
            breaks = c(0, 25, 50, 100, 150, 200, Inf),
            include.lowest = TRUE
          )
        levs <- levels(x_cut)
        new_levs <- 
          levs %>%
          str_remove_all("\\[|\\]|\\(|\\)") %>%
          str_replace(",", "-")
        
        x_cut <- fct_recode(x_cut, !!!structure(levs, names = new_levs))
        res <- table(x_cut) %>% as.data.frame()
        colnames(res) <- c("radius", "count")
        res
      })
    
    dat <-
      dat_ls %>%
      enframe(name = "cell") %>% 
      bind_rows() %>% 
      unnest(value)
    dat
  })

res_nn_count_df <-
  res_nn_count_df_ls %>%
  enframe(name = "project") %>%
  unnest(value)  %>%
  # left_join(distinct(all_data_meta_cn_subset[c("project", "group")]), by = "project") %>%
  dplyr::filter(radius != "200-Inf")


# plot  - group
p_proximity_count <-
  res_nn_count_df %>%
  ggplot(aes(x = radius, y = count, fill = cell) ) +
  geom_col() +
  # facet_wrap(~ group) +
  ggsci::scale_fill_igv() +
  theme_bw() +
  labs(
    x = NULL,
    y = NULL
  )


p_proximity_freq <-
  res_nn_count_df %>%
  group_by(radius, cell) %>%
  summarise(count = sum(count)) %>%
  
  group_by(radius) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = radius, y = freq, fill = cell)) +
  geom_col() +
  ggsci::scale_fill_igv() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(
    x = NULL,
    y = NULL
  )

# export 
ggsave(
  p_proximity_count,
  filename = file.path(out_dir_260325, str_glue("2 proximity analysis barplot - count.pdf")),
  width  = 8,
  height = 6
)

ggsave(
  p_proximity_freq,
  filename = file.path(out_dir_260325, str_glue("2 proximity analysis barplot - freq.pdf")),
  width  = 8,
  height = 6
)



# plot - project
df_proj_grp <- 
  all_data_meta_pos_subset[c('group', 'Slide_ID')] %>%
  distinct()
  
p_proximity_grp_freq <-
  res_nn_count_df %>%
  group_by(project, radius, cell) %>%
  summarise(count = sum(count)) %>%
  
  left_join(df_proj_grp, by = c("project" = "Slide_ID")) %>% 
  mutate(group = fct_recode(
      group,
      NAG = "1NAG",
      others = "SMC",
      others = "TLS"  
    ) %>%
    fct_relevel(
      'NAG', 'CAG', 'IM', 'LGD', 'HGD',
      'early_cancer', 'late_cancer', 'others'
    )
  )  %>%
  
  group_by(group, radius) %>%
  mutate(freq = count/sum(count)) %>%
  # mutate(project = str_c(project, sep = "_")) %>%
  # mutate(project = factor(project, levels = sort(unique(project)))) %>%
  
  
  ggplot(aes(x = radius, y = freq, fill = cell)) +
  geom_col(width = 0.75) +
  facet_wrap(~ group, ncol = 3) +
  ggsci::scale_fill_igv() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(
    x = NULL,
    y = NULL
  )

ggsave(
  p_proximity_grp_freq,
  filename = file.path(out_dir_260325, str_glue("2 proximity analysis barplot - group - freq.pdf")),
  width  = 12,
  height = 24
)

