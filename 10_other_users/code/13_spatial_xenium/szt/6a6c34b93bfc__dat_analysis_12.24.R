library(Seurat)
library(tidyverse)


# 0. ----------------
out_dir <- "dat_analysis_cycle3"
cn_out_dir <- "dat_analysis_cycle3/cn_analysis" 
dir.create(cn_out_dir)


all_data <- readRDS(file.path(out_dir, "subtype", "all_data_cycle3_1205.rds"))
all_data_meta_cn <- readRDS("~/dat_analysis_cycle3/cn_analysis/meta_data_cn.rds")

# subset samples
group_A <- c(
  '13638',
  '10858',
  '47983A3'
)

group_B <- c(
  '45827',
  '72696',
  '58198',
  '71402A6',
  '06938F1',
  '15201H1',
  '28519G1',
  '44994A1',
  '43158G1',
  '67718G1'
)

stopifnot( all(group_A %in% all_data$project) )
stopifnot( all(group_B %in% all_data$project) )

all_data_subset <- all_data[, all_data$project %in% c(group_A, group_B)]
all_data_meta_cn_subset <- all_data_meta_cn[all_data$project %in% c(group_A, group_B), ]

# 2. proximity analysis - only CN9 -----------

## 2.1. input data ------------
cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- 
      all_data_meta_cn_subset %>%
      # rownames_to_column("cell_id") %>% 
      dplyr::filter(CN == "CN9") %>% 
      dplyr::select(annotation_full = annotation4, group, project, cell_id) %>%
      column_to_rownames("cell_id")
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(annotation_full, group, project) %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(annotation_full = factor(annotation_full)) %>%
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


## 2.2. proximity analysis -------
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

### 2.2.1 proximity analysis --------
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

# saveRDS(res_pa, file = file.path(cn_out_dir, "res_pa.rds"))

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

# 
# res_pa_df2  <-
#   res_pa_df2 %>%
#   as.data.frame() %>%
#   rownames_to_column('group') %>%
#   pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
#   separate(group, into = c("Target", "Group", "project", "CN2"), sep = "\\.\\.") 


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa  <-
  res_pa_df %>%
  # # remove Others cells, if cell_props all equals to 0
  # group_by(Target, Others) %>%
  # dplyr::filter(rep(! all(cell_props == 0), n())) %>%
  
  # dplyr::filter(Target == "Epithelial cells") %>% 
  # dplyr::filter(Target == "apCAF") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() + 
  ggpubr::geom_pwc(
    aes(group = Group),
    tip.length = 0,
    hide.ns  = F,
    method = "t_test",
    # label = "p.format"
  ) +
  scale_y_continuous(limits = function(x) c(x[1], 1.1 * x[2])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 


  # res_pa_df %>%
  #   # dplyr::filter(Target == "Epithelial cells") %>%
  #   dplyr::filter(Target == "apCAF")  %>%
  #   ggpubr::compare_means(
  #     data = ., 
  #     cell_props ~ Group,
  #     method = "t.test",
  #     group.by = "Others"
  #   )
  
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
  file.path(cn_out_dir, "proximity analysis.pdf"),
  height = 50,
  width  = 25
)
p_res_pa %>% print()
dev.off()

# 
# pdf(
#   file.path(cn_out_dir, "proximity analysis - CN.pdf"),
#   height = 50,
#   width  = 25
# )
# p_res_pa_ls %>% print()
# dev.off()

### 2.2.2 proximity analysis 2 -------
##total data  clean
all_data_meta_pa <-
  all_data[[]] %>%
  mutate(annotation4 = ifelse(
    as.character(annotation_full) == "Fibroblast",
    as.character(annotation3),
    as.character(annotation_full)
  )) %>%
  dplyr::select(-annotation_full) %>%
  dplyr::rename(annotation_full = annotation4) %>%
  mutate(annotation_final = annotation_full %>%
           as.character() %>%
           factor(., levels = unique(.))
  ) %>%
  mutate(
    Group = group
  )



all_data_meta_pos <-
  cells_position_ls  %>%
  imap(function(cells_position, slide){# browser()
    
    if('cell_id' %in% colnames(all_data_meta_pa) ){
      all_data_meta_pa <- all_data_meta_pa %>%
        remove_rownames() %>%
        column_to_rownames("cell_id")
    }
    
    proj_meta_sub <-
      all_data_meta_pa %>%
      dplyr::select(annotation_final, Group, project) %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    if(! 'cell_id' %in% colnames(proj_meta_sub) ){
      proj_meta_sub <- rownames_to_column(proj_meta_sub, "cell_id")
    }
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%
      dplyr::rename(Classification1 = annotation_final,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = project)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")


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
    mpp = 1, # !important, 1 for xenium
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
    mpp = 1, # !important, 1 for xenium
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
    mpp = 1, # !important, 1 for xenium
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
    mpp = 1, # !important, 1 for xenium
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
    mpp = 1, # !important, 1 for xenium
    returntype = "count",
    cell_enough_ratio_cutoff = 99
  )

# saveRDS(res_pa_ls, file.path(cn_out_dir, "res_pa_ls.rds"))


# epithelial proximity
# proximity_cells <- c(
#   "THBS1 Fibroblast",
#   "TNXB Fibroblast",
#   "COL5A1 Fibroblast",
#   "EEF1G Fibroblast",
#   "ADAMDEC1 Fibroblast"
# )
proximity_cells <- c(
  'Immediate early-response (IER) fibroblasts',
  'S3 fibroblasts',
  'S1 fibroblasts',
  'Peri-cryptal S2 fibroblasts' 
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
    all_data_group_sub <- all_data_sub[c('Group', 'Slide_ID')] %>% distinct()
    
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
      
      ggplot(aes(x = Cell, y = count, fill = Group)) +
      # geom_violin(trim = FALSE) +
      geom_boxplot() +
      # geom_point() +
      facet_wrap(~ distance, scales = "free_y") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = grp_cols_sp) +
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
      group_by(Cell, distance, Group) %>%
      summarise(count = median(count)) %>%
      
      mutate(distance = as.numeric(str_remove(distance, "um$"))) %>%
      
      ggplot(aes(x = distance, y = count, color = Group)) +
      # geom_density(alpha = 0.6) +
      geom_smooth(method = 'loess') +
      facet_wrap(~ Cell, scales = "free_y") +
      ggpubr::theme_pubr() +
      theme(
        strip.background = element_blank(),
        panel.grid = element_blank()
      ) +
      scale_color_manual(values = grp_cols_sp) +
      labs(
        x = 'Distance(um)',
        y = 'Median count',
        title = name
      )
    
    p
  })




# export pdf
pdf(file.path(cn_out_dir, "proximity analysis - barplot.pdf"), width = 14, height = 8)
p_proximity_cell_count_distribution_ls %>% walk(print)
dev.off()

# # export pdf
# pdf("proximity analysis - bxolplot - combined.pdf", width = 18, height = 8)
# print(p_proximity_cell_count_distribution)
# dev.off()


# export pdf
pdf(file.path(cn_out_dir, "proximity analysis - smooth_lines.pdf"), width = 14, height = 8)
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
      group_by(Group, Cell, distance) %>%
      summarise(count = median(count)) %>%
      ggplot(aes(x = Cell, y = distance, size = count, color = Group)) +
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
      scale_color_manual(values = grp_cols_sp)+
      coord_polar(theta = 'x', clip = 'off') +
      labs(x = NULL, y = NULL, title = name)
    p 
  })





# export pdf
pdf(file.path(cn_out_dir, "proximity analysis - radarplot.pdf"), width = 10, height = 8)
p_proximity_radar_ls %>% walk(print)
dev.off()





# 3. 肿瘤细胞各种通路打分 --------------

# all_data_subset[[c(
#   'HLA', 'tumor_cell_LN_metastasis', 
#   'IFN', 'cancer_stem_cell_and_EMT', 
#   'Inflammation', 'project', 'group'
# )]]



path_de_df <- 
  c(
    'HLA', 'tumor_cell_LN_metastasis', 
    'IFN', 'cancer_stem_cell_and_EMT', 
    'Inflammation'
  ) %>%
  set_names(., .) %>%
  map(function(pathway){
    dat <- 
      all_data_meta_cn_subset %>%
      dplyr::filter(annotation4 == 'Epithelial') %>%
      .[c(pathway, 'project', 'group')] 
    
    dat <- dat %>%
      group_by(project, group) %>%
      summarise( {{pathway}} := mean(.data[[pathway]]) ) %>%
      as.data.frame()
    
    pathway_a <- dat[[pathway]][dat$group == "A"] 
    pathway_b <- dat[[pathway]][dat$group == "B"] 
    
    delta <- mean(pathway_a) - mean(pathway_b)
    wtest <- wilcox.test(pathway_a, pathway_b)
    
    w_stat <- wtest$statistic
    pvalue <- wtest$p.value
    
    data.frame(
      statistic = w_stat,
      pvalue = pvalue,
      delta = delta
    )
  }) %>%
  enframe("pathway") %>%
  unnest(value)

# plot
p_path_de_df <-
  path_de_df %>%
  # mutate(pvalue = ifelse(pvalue < 1e-30, 0, pvalue)) %>%
  mutate(pathway = reorder(pathway, delta)) %>%
  ggplot(aes(y = pathway, x = delta, fill = statistic, size = pvalue)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(Inf, "Reds")) +
  ggpubr::theme_pubclean() +
  guides(
    color = guide_colorbar(theme = theme(
      legend.key.width  = unit(10, "lines") 
    )),
    size = guide_legend(reverse = TRUE)
  ) +
  scale_size(range = c(2, 6), transform = "reverse") + 
  theme(
    # legend.spacing = unit(2, "lines")
    # legend.position = "right"
  ) +
  labs(
    fill = "wilcox statistic",
    x = "delta (pathway score A - pathway score  B)",
    y = NULL
  )

ggsave(
  p_path_de_df,
  filename = file.path(cn_out_dir, '3. dotplot plot of pathway score - A vs B - epithelial.pdf'),
  height = 4,
  width  = 8
)

## 3.2 A vs B - one CN ---------


path_de_df_cn <- 
  c(
    'HLA', 'tumor_cell_LN_metastasis', 
    'IFN', 'cancer_stem_cell_and_EMT', 
    'Inflammation'
  ) %>%
  set_names(., .) %>%
  map(function(pathway){
    dat <- 
      all_data_meta_cn_subset %>%
      dplyr::filter(annotation4 == 'Epithelial') %>%
      .[c(pathway, 'project', 'group', 'CN2')]

    dat2 <-
      dat %>%
      group_by(project, group, CN2) %>%
      summarise( {{pathway}} := mean(.data[[pathway]]) ) %>%
      as.data.frame()
    
    # print(head(dat2))
    
    dat2 %>%
      group_by(CN2) %>%
      nest() %>%
      deframe() %>%
      map(function(x){
        pathway_a <- x[[pathway]][x$group == "A"] 
        pathway_b <- x[[pathway]][x$group == "B"] 
        
        delta <- mean(pathway_a) - mean(pathway_b)
        ttest <- wilcox.test(pathway_a, pathway_b)
        
        t_stat <- ttest$statistic
        pvalue <- ttest$p.value
        
        data.frame(
          statistic = t_stat,
          pvalue = pvalue,
          delta = delta
        )
      }) %>%
      enframe("CN2") %>%
      unnest(value)
    
  }) %>%
  enframe("pathway") %>%
  unnest(value)


# plot
p_path_de_df_cn <-
  path_de_df_cn %>%
  mutate(pathway = factor(pathway, levels = path_de_df_cn %>% group_by(pathway) %>% summarise(s = sum(delta)) %>% arrange(s) %>% .[[1]])) %>%
  ggplot(aes(y = pathway, x = delta, fill = statistic, size = pvalue)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(Inf, "Reds")) + 
  facet_wrap(~ CN2) +
  ggpubr::theme_pubclean() +
  guides(
    fill = guide_colorbar(theme = theme(
      legend.key.width  = unit(10, "lines") 
    )),
    size = guide_legend(reverse = TRUE)
  ) +
  scale_size(range = c(2, 6), transform = "reverse") + 
  theme(
    # legend.spacing = unit(2, "lines")
  ) +
  labs(
    fill = "wilcox statistic",
    x = "delta (pathway score A - pathway score  B)",
    y = NULL
  )

ggsave(
  p_path_de_df_cn,
  filename = file.path(cn_out_dir, '3. dotplot plot of pathway score - A vs B - split by CN - epithelial.pdf'),
  height = 8,
  width  = 15
)


## 3.3 pairwise CN - A and B ---------

path_de_df_group_ls <- 
  c(
    'HLA', 'tumor_cell_LN_metastasis', 
    'IFN', 'cancer_stem_cell_and_EMT', 
    'Inflammation'
  ) %>%
  set_names(., .) %>%
  map(function(pathway){
    dat <- 
      all_data_meta_cn_subset %>%
      dplyr::filter(annotation4 == 'Epithelial') %>%
      .[c(pathway, 'project', 'group', 'CN2')]
    
    dat2 <-
      dat %>%
      group_by(project, group, CN2) %>%
      summarise( {{pathway}} := mean(.data[[pathway]]) ) %>%
      as.data.frame()
    
    dat2 %>%
      group_by(group) %>%
      nest() %>%
      deframe() %>%
      map(function(x){ # browser()
        x$CN2 %>%
          unique() %>%
          set_names(., .) %>% 
          map(function(cn){
            pathway_a <- x[[pathway]][x$CN2 == cn] 
            pathway_b <- x[[pathway]][x$CN2 != cn] 
            
            delta <- mean(pathway_a) - mean(pathway_b)
            ttest <- wilcox.test(pathway_a, pathway_b)
            
            t_stat <- ttest$statistic
            pvalue <- ttest$p.value
            
            data.frame(
              statistic = t_stat,
              pvalue = pvalue,
              delta = delta
            )
          }) %>%
          enframe("CN2") %>%
          unnest(value)
        
      }) %>%
      enframe("group") %>%
      unnest(value)
  })  


# plot
p_path_de_df_group_ls <-
  path_de_df_group_ls %>%
  imap(function(dat, name){
    dat %>%
      ggplot(aes(y = CN2, x = delta, fill = statistic, size = pvalue)) + 
      geom_point(shape = 21, color = "black") +
      scale_fill_gradientn(colours = RColorBrewer::brewer.pal(Inf, "Reds")) + 
      facet_wrap(~ group) +
      ggpubr::theme_pubclean() +
      guides(
        fill = guide_colorbar(theme = theme(
          legend.key.width  = unit(10, "lines") 
        )),
        size = guide_legend(reverse = TRUE)
      ) +
      scale_size(range = c(2, 6), transform = "reverse") + 
      theme(
        # legend.spacing = unit(2, "lines")
      ) +
      labs(
        title = name, 
        fill  = "wilcox statistic",
        x = "delta (pathway score of CN - average pathway score of Others CNs)",
        y = NULL
      )
  })


pdf(
  file = file.path(cn_out_dir, '3. dotplot plot of pathway score - CN - split by group - epithelial.pdf'),
  height = 6,
  width  = 10
)
p_path_de_df_group_ls %>% walk(print)
dev.off()

# 4. proximity analysis of HLA-high cells -------

# HLA score distribution
p_hla_distribution <-
  all_data_meta_cn_subset %>% 
  ggplot(aes(x = project, y = HLA)) +
  geom_boxplot() +
  theme_bw() +
  Seurat::RotatedAxis()

ggsave(
  p_hla_distribution,
  filename = file.path(cn_out_dir, "HLA score distribution.pdf"),
  height = 6, 
  width = 8
)

# proximity
all_data_meta_cn_subset_hla <-
  all_data_meta_pos_subset %>%
  left_join(all_data_meta_cn_subset[c('HLA', 'project', 'cell_id')], by = 'cell_id') %>%
  group_by(project) %>%
  mutate(HLA_group = ifelse(HLA > median(HLA), "High", "Low"))  


# summarized by group and celltype
pa_groups_hla <-
  with(
    all_data_meta_cn_subset_hla,
    str_c(Classification1, group, Slide_ID, HLA_group, sep = '..')
  )


# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df_hla <-  apply(res_pa, 1, function(x){tapply(x, INDEX = pa_groups_hla, FUN = "median")})

res_pa_df_hla  <-
  res_pa_df_hla %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project", "HLA_group"), sep = "\\.\\.") 



# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa_hla  <-
  res_pa_df_hla %>%
  # remove Others cells, if cell_props all equals to 0
  group_by(Target, Others) %>%
  dplyr::filter(rep(! all(cell_props == 0), n())) %>%
  
  dplyr::filter(Target == "Epithelial") %>%
  dplyr::filter(HLA_group == "High") %>%
  
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
  ggpubr::theme_pubclean() + 
  ggpubr::geom_pwc(
    aes(group = Group),
    tip.length = 0,
    hide.ns  = F,
    method = "t_test",
    # label = "p.format"
  ) +
  scale_y_continuous(limits = function(x) c(x[1], 1.1 * x[2])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 



# export plot
pdf(
  file.path(cn_out_dir, "proximity analysis_hla.pdf"),
  height = 8,
  width  = 15
)
p_res_pa_hla %>% print()
dev.off()


# 5. dotplot of DCs pDCs Macrophage Mono-like Macrophage --------
mye_cells_to_check <- c(
  'DCs',
  'pDCs',
  'Macrophage',
  'Mono-like Macrophage'
)
all_data_mye_check <- all_data[, all_data$annotation4 %in% mye_cells_to_check]

DotPlot(
  all_data_mye_check,
  group.by = "annotation4",
  features = "CD68"
)

all_data_mye_check[['sketch']] <- NULL
Idents(all_data_mye_check) <- "annotation4"


p_dot_check <-
  DotPlot(
    all_data_mye_check,
    group.by = "annotation4",
    features = c(
      "Myeloid" = "ITGAM",
      "Macrophage" = "CD68",
      "Macrophage" = "CD163",
      "Macrophage" = "C1QA",
      "Macrophage" = "C1QC",
      
      "Mast cell" = "CPA3",
      "Mast cell" = "CTSG",
      "Mast cell" = "TPSAB1",
      
      "Monocytes" = "CD14", 
      "Monocytes" = "FCGR3A",
      
      "Dendritic cell" = "CD83",
      'Dendritic cell' = "LILRA4",  
      
      'Dendritic cell' = "CLEC9A",
      'Dendritic cell' = "IDO1",
      
      'Dendritic cell' = "CD1C",
      'Dendritic cell' = "FCER1A"
    )
  )

ggsave(
  p_dot_check + RotatedAxis(),
  filename = file.path(cn_out_dir, "marker dotplot of 4 myeloid cells.pdf"),
  height = 4,
  width  = 12
)

# 4. 废弃 proximity analysis of HLA-high cells ------
all_data_meta_cn_subset_hla <-
  all_data_meta_cn_subset %>%
  mutate(HLA_group = ifelse(HLA > median(HLA), "High", "Low")) %>%
  dplyr::filter(HLA_group == "High")



# all_data_cn_meta
all_data_meta_pos_hla <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- 
      all_data_meta_cn_subset_hla %>%
      # rownames_to_column("cell_id") %>% 
      dplyr::filter(CN == "CN9") %>% 
      dplyr::select(annotation_full = annotation4, group, project, cell_id) %>%
      column_to_rownames("cell_id")
    
    proj_meta_sub <-
      all_data_meta %>%  
      dplyr::select(annotation_full, group, project) %>%
      rownames_to_column("cell_id") %>%
      dplyr::filter(project == slide)
    
    cells_position[, 1:3] %>%
      inner_join(proj_meta_sub, by = "cell_id") %>%  
      mutate(annotation_full = factor(annotation_full)) %>%
      dplyr::rename(Classification1 = annotation_full,
                    X = x_centroid,
                    Y = y_centroid,
                    Slide_ID = project)
  }) %>%
  enframe() %>%
  unnest(cols = value) %>%
  dplyr::select(-name) %>%
  column_to_rownames("cell_id")

all_data_meta_pos_subset_hla <- all_data_meta_pos_hla %>% rownames_to_column("cell_id")
# readRDS('data/all_data_meta_pos_subset.rds')

all_data_meta_pos_subset_hla %>% head

# cell_id        X        Y Classification1 group       Slide_ID
# 1 1_aaaaaagg-1 1670.389 185.1358      Macrophage  CD-F X_202413076_32
# 2 1_aaadlghn-1 1681.172 176.4662          Plasma  CD-F X_202413076_32
# 3 1_aaadmplc-1 1677.561 182.6699               T  CD-F X_202413076_32
# 4 1_aaadoogj-1 1677.109 165.6637          Plasma  CD-F X_202413076_32
# 5 1_aaadpljf-1 1671.887 176.8148          Plasma  CD-F X_202413076_32
# 6 1_aaafaiko-1 1719.466 192.8071          Plasma  CD-F X_202413076_32
# ......


## 4.2. proximity analysis -------
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

## 4.2.1 proximity analysis --------
res_pa_hla <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset_hla,
    query = all_data_meta_pos_subset_hla,
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

# saveRDS(res_pa_hla, file = file.path(cn_out_dir, "res_pa_hla.rds"))

# return res_pa_hla
# rows: all celltypes
# cols: query cell
res_pa_colsum_hla <- colSums(res_pa_hla)
res_pa_hla <- sweep(res_pa_hla, 2, res_pa_colsum_hla, FUN = "/")
if(any(res_pa_colsum_hla == 0)){
  res_pa_hla[, res_pa_colsum_hla == 0] <- 0
}

# res_pa_hla <- readRDS('data/res_pa_hla.rds')


# summarized by group and celltype
pa_groups_hla <- 
  with(
    all_data_meta_pos_subset_hla,
    str_c(Classification1, group, Slide_ID, sep = '..')
  )

# pa_groups2 <-
#   all_data_meta_pos_subset_hla %>%
#   left_join(all_data_meta[c('cell_id', 'CN2')], "cell_id") %>%
#   with(
#     str_c(Classification1, group, Slide_ID, CN2, sep = '..')
#   )

# 分组内取细胞比例的中位数
# apply对每行操作时会转置结果
res_pa_df_hla <-  apply(res_pa_hla, 1, function(x){tapply(x, INDEX = pa_groups_hla, FUN = "median")})


res_pa_df_hla  <-
  res_pa_df_hla %>%
  as.data.frame() %>%
  rownames_to_column('group') %>%
  pivot_longer(-group, names_to = "Others", values_to = "cell_props") %>%
  separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 


# plot boxplot
# Target:中心细胞
# Others: 其他细胞占比

p_res_pa_hla <-
  res_pa_df_hla %>%
  # dplyr::filter(Target == "Epithelial cells") %>%
  ggplot(aes(x = Others, y = cell_props, fill = Group)) +
  geom_boxplot(width = 0.45, outlier.colour = NA) +
  # geom_jitter(width = 0.25) +
  facet_wrap(~ Target, scales = "free_y", ncol = 2)+
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
  file.path(cn_out_dir, "proximity analysis_hla.pdf"),
  height = 50,
  width  = 25
)
p_res_pa_hla %>% print()
dev.off()

# 
# pdf(
#   file.path(cn_out_dir, "proximity analysis - CN.pdf"),
#   height = 50,
#   width  = 25
# )
# p_res_pa_ls %>% print()
# dev.off()












