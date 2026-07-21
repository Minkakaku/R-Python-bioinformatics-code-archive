library(Seurat)
library(tidyverse)



# 0. ----------------
out_dir <- "dat_analysis_cycle3/result_20250220" 
dir.create(out_dir, recursive = TRUE)

cn_out_dir <- "dat_analysis_cycle3/cn_analysis_250219"  
dir.create(cn_out_dir) 


all_data <- readRDS(file.path(out_dir, "subtype", "all_data_cycle3_250219.rds"))
all_data_meta_cn <- readRDS(file.path(cn_out_dir, "meta_data_cn.rds"))


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
      dplyr::filter(CN == "CN9") %>% # <-------- only CN9 ----------- 
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

saveRDS(res_pa, file = file.path(out_dir, "res_pa.rds"))

# return res_pa
# rows: all celltypes
# cols: query cell
res_pa_colsum <- colSums(res_pa)
res_pa <- sweep(res_pa, 2, res_pa_colsum, FUN = "/")
if(any(res_pa_colsum == 0)){
  res_pa[, res_pa_colsum == 0] <- 0
}


# summarized by group and celltype
pa_groups <-
  with(
    all_data_meta_pos_subset,
    str_c(Classification1, group, Slide_ID, sep = '..')
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

# 来源于dat_analysis_2024.9.5.R
# cells_to_analysis <- c(
#   "Epithelial",
#   "CD4T",
#   "CD8T",
#   "Epi0_CD4T",
#   "Epi0_CD83_T",
#   "Epi0_T"
# )

cells_to_analysis <- c(
  "Epithelial",
  "CD4T",
  "CD8T",
  "CD4T-Macrophage",
  "Epi0_CD4T",
  "Epi0_CD83_T",
  "Epi0_T",
  "CD4 Endothelial"
)


p_res_pa_ls <-
  res_pa_df %>%
  with(split(res_pa_df, CN2)) %>%
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(Target == "Epithelial cells") %>%
      dplyr::filter(Target %in% cells_to_analysis) %>%
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
  file.path(out_dir, "proximity analysis - 2.pdf"),
  height = 18, # 14,
  width  = 18
)
p_res_pa_ls %>% print()
dev.off()



## 2.2 proximity analysis - series analysis -----
init_cells <- c(
  '25'  = 750,
  '50'  = 800,
  '100' = 1000,
  '150' = 1200,
  '200' = 1500
)

proximity_series_ls <- c(25, 50, 100, 150, 200)

res_pa_ls <-
  proximity_series_ls %>%
  set_names(., .) %>%
  map(function(x){
    k = init_cells[as.character(x)]
    res_pa <-
      ProximutyAnalysis_NN_multi(
        all_data_meta_pos_subset,
        query = all_data_meta_pos_subset,
        x = "X",
        y = "Y",
        group.by = "Slide_ID",
        classification = "Classification1", # 若为NULL: 提前返回window的细胞近邻的list
        k = k, # 'auto',
        searchtype = "radius",
        radius = x, # 邻域半径
        mpp = 1, # !important, 1 for xenium
        returntype = "count",
        cell_enough_ratio_cutoff = 90
      )
    res_pa
  })


saveRDS(res_pa_ls, file.path(out_dir, "res_pa_ls.rds"))

# 将同心圆统计 转换成 环状同心的统计
# 注意一定要倒序操作
# 200 = 200 - 150
# 150 = 150 - 100 
# ...
# 
res_pa_ls2 <- res_pa_ls
diff_idx_ls <-
  proximity_series_ls %>% 
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
        all_data_meta_pos_subset,
        str_c(Classification1, group, Slide_ID, sep = '..')
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
  dplyr::filter(Target == "Epithelial")

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
  file.path(out_dir, "proximity analysis_series distance.pdf"),
  height = 14,
  width  = 18
)
p_proximity_series_epi %>% print()
dev.off()


# write out cell count of others cells of epithelial 


# 每个样本中的每种Target的其他细胞的计数取中位数
res_pa_count_tidy_sum <- 
  res_pa_ls %>%
  imap(function(x, distance){# browser()
    # 分组
    pa_groups <-
      with(
        all_data_meta_pos_subset,
        str_c(Classification1, group, Slide_ID, sep = '..')
      )
    # 分组内取加和
    res <- apply(x, 1, function(y){tapply(y, INDEX = pa_groups, FUN = "sum")})
    
    # data tidy
    res_df  <-
      res %>%
      as.data.frame() %>%
      rownames_to_column('group') %>%
      pivot_longer(-group, names_to = "Others", values_to = "cell_counts") %>%
      separate(group, into = c("Target", "Group", "project"), sep = "\\.\\.") 
    res_df$distance <- distance
    
    # return
    res_df
  }) %>%
  bind_rows()


res_pa_count_tidy_epi_sum <- res_pa_count_tidy_sum %>% dplyr::filter(Target == "Epithelial")


res_pa_count_tidy_epi_sum2 <-
  res_pa_count_tidy_epi_sum %>%
  pivot_wider(names_from = "distance", values_from = "cell_counts") %>%
  mutate(
    Total = `200`,
    `200` = `200` - `150`,
    `150` = `150` - `100`,
    `100` = `100` - `50`,
    `50`  = `50`  - `25`,
    `25`  = `25`
  ) %>%
  pivot_longer(
    c(`25`, `50`, `100`, `150`, `200`, Total),
    names_to = "distance",
    values_to = "cell_counts"
  ) %>%
  mutate(
    Others = str_c(Others, distance, sep = "_"),
    project = str_c(Group, project, sep = "_")
  ) %>%
  dplyr::select(
    project, Others, cell_counts
  ) %>%
  pivot_wider(names_from = "project", values_from = "cell_counts")
  

# export 
openxlsx::write.xlsx( 
  res_pa_count_tidy_epi_sum2, 
  file = file.path(out_dir, "proximity analysis of epithelial in CN9 - cell counts table.xlsx"), 
  asTable = TRUE 
) 









