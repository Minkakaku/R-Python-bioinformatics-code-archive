library(Seurat)
library(tidyverse)
source("spatial analysis.R")

# copy from dat_analysis_250220.R
# 更新细胞注释为最新的结果

# 0. ----------------
out_dir <- "dat_analysis_cycle3/result_20250422"
dir.create(out_dir, recursive = TRUE)

cn_out_dir <- "dat_analysis_cycle3/cn_analysis_250219"  
dir.create(cn_out_dir) 


# all_data <- readRDS(file.path(out_dir, "subtype", "all_data_cycle3_250219.rds"))
all_data <- readRDS(file.path("dat_analysis_cycle3", "subtype", "all_data_cycle3_250219.rds"))
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

# compute_cells <- "Immune_cells" 
# compute_cells <- "CD8_T_cells" 
compute_cells <- "CD4_T_cells"

# 免疫细胞大群包括：CD4T细胞、CD8T细胞、B细胞、cDC1、Macrophage、Mast、Monocyte、NK、pDCs、proliferating-T共10类
Immune_cells <- c(
  'B_cells' = 'ASC',
  'B_cells' = 'B-GC', 
  'B_cells' = 'B-memory', 
  'B_cells' = 'B-naive',
  
  'CD8_T_cells' = 'CD8-Tnaive', 
  'CD8_T_cells' = 'CD8-Tem', 
  'CD8_T_cells' = 'CD8-Tex',
  'CD8_T_cells' = 'CD8-Trm',
  
  'CD4_T_cells' = 'Treg',
  'CD4_T_cells' = 'Tfh',
  'CD4_T_cells' = 'CD4-Tnaive',
  'CD4_T_cells' = 'CD4-Tmemory',
  'CD4_T_cells' = 'CD4T-Tcm', 
  
  'proliferating-T' = 'proliferating-T', 
  'NK'   = 'NK', 
  'Mast' = 'Mast', 
  'cDC1' = 'cDC1', 
  'pDCs' = 'pDCs',
  'Monocyte'   = 'Monocyte', 
  'Macrophage' = 'Macrophage'
)

CD4_T_cells <- c(
  'Treg',
  'Tfh',
  'CD4-Tnaive',
  'CD4-Tmemory',
  'CD4T-Tcm'
)


CD8_T_cells <- c(
  'CD8-Tnaive', 
  'CD8-Tem', 
  'CD8-Tex',
  'CD8-Trm'
) 

 # "Immune_cells" "CD8_T_cells" "CD4_T_cells"

if(compute_cells == "Immune_cells"){
  
  all_data_meta_pos_subset2 <- 
    all_data_meta_pos_subset %>%
    dplyr::filter(Classification1 %in% c('Epithelial', Immune_cells)) %>%
    mutate(Classification1 = fct_recode(Classification1, !!!Immune_cells)) %>%
    mutate(Classification1 = fct_drop(Classification1))
  
}else if(compute_cells == "CD8_T_cells"){
  
  all_data_meta_pos_subset2 <- 
    all_data_meta_pos_subset %>%
    dplyr::filter(Classification1 %in% c('Epithelial', CD8_T_cells)) %>%
    mutate(Classification1 = fct_drop(Classification1))
  
}else if(compute_cells == "CD4_T_cells"){
  
  all_data_meta_pos_subset2 <- 
    all_data_meta_pos_subset %>%
    dplyr::filter(Classification1 %in% c('Epithelial', CD4_T_cells)) %>%
    mutate(Classification1 = fct_drop(Classification1))
  
}else{
  stop("compute_cells can only be Immune_cells, CD8_T_cells or CD4_T_cells")
}


all_data_meta_pos_subset2 %>% head

# cell_id        X         Y Classification1 group Slide_ID
# 1 5_aaaaadgh-1 2615.709 346.7925      Epithelial     B  06938F1
# 2 5_aaaaamag-1 2617.601 380.8063      Epithelial     B  06938F1
# 3 5_aaaaijpf-1 2603.295 377.5136      Macrophage     B  06938F1
# 4 5_aaaaoega-1 2607.563 349.8770      Epithelial     B  06938F1
# 5 5_aaajpabc-1 2552.105 349.1224         B_cells     B  06938F1
# 6 5_aadgplbj-1 2628.076 485.5199      Epithelial     B  06938F1
# ......

all_data_meta_pos_subset2$Classification1 %>% unique
# [1] Epithelial      Macrophage      B_cells         CD8_T_cells     CD4_T_cells     proliferating-T NK              Mast            cDC1           
# [10] Monocyte        pDCs           
# Levels: B_cells CD4_T_cells CD8_T_cells cDC1 Epithelial Macrophage Mast Monocyte NK pDCs proliferating-T

all_data_meta_pos_subset2_epi <- all_data_meta_pos_subset2[all_data_meta_pos_subset2$Classification1 == "Epithelial", ]

## 2.2 proximity analysis --------------- 
res_pa_200um <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset2,
    query = all_data_meta_pos_subset2_epi,
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


saveRDS(res_pa_200um, file = file.path(out_dir, str_glue("{compute_cells}_res_pa_200um.rds")))


cell_match_var <- structure(
  as.character(all_data_meta_pos_subset2$Classification1),
  names = all_data_meta_pos_subset2$cell_id
)

# res_nn_df_ls <-
#   res_pa_200um %>%
#   # .[1] %>%
#   map(function(nn){ # browser()
#     res_ls <- list()
#     for(i in 1:nrow(nn$nn.idx)){
#       if(i %% 1000 == 0) print(i)
#         
#       cell_current <- nn$total_cell_id[i]
#       
#       nn.idx  <- nn$nn.idx[i, ]
#       nn.dist <- nn$nn.dists[i, ]
#       
#       idx <- nn.idx > 1 # remove 0 1
#       nn.idx  <- nn.idx[idx]  
#       nn.dist <- nn.dist[idx]
#       
#       
#       cell_id <- nn$query_cell_id[nn.idx]
#       cells   <- cell_match_var[cell_id]
#       
#       # head 10 cells, mean
#       res <- tapply(nn.dist, cells, function(x){mean(head(x, 10))})
#       res_ls[[cell_current]] <- res
#     }
#     res_df <- as.data.frame(bind_rows(res_ls))
#     rownames(res_df) <- nn$total_cell_id[i]
#     res_df
#   })

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
    xx <- x[, !colnames(x) %in% c("Others", "Epithelial")]
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
  unnest(value) %>% 
  left_join(distinct(all_data_meta_cn_subset[c("project", "group")]), by = "project") %>%
  dplyr::filter(radius != "200-Inf")

# plot  - group
p_proximity_count <-
  res_nn_count_df %>%
  ggplot(aes(x = radius, y = count, fill = cell) ) +
  geom_col() +
  facet_wrap(~ group) +
  ggsci::scale_fill_igv() +
  theme_bw() +
  labs(
    x = NULL,
    y = NULL
  )


p_proximity_freq <-
  res_nn_count_df %>%
  group_by(group, radius, cell) %>%
  summarise(count = sum(count)) %>%
  
  group_by(group, radius) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = radius, y = freq, fill = cell)) +
  geom_col() +
  facet_wrap(~ group) +
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
  filename = file.path(out_dir, str_glue("{compute_cells} - proximity analysis barplot - count.pdf")),
  width  = 8,
  height = 6
)

ggsave(
  p_proximity_freq,
  filename = file.path(out_dir, str_glue("{compute_cells} - proximity analysis barplot - freq.pdf")),
  width  = 8,
  height = 6
)

# plot - project
p_proximity_proj_freq <-
  res_nn_count_df %>%
  group_by(group, project, radius, cell) %>%
  summarise(count = sum(count)) %>%
  
  group_by(project, radius) %>%
  mutate(freq = count/sum(count)) %>%
  mutate(project = str_c(group, project, sep = "_")) %>%
  mutate(project = factor(project, levels = sort(unique(project)))) %>%
  
  ggplot(aes(x = radius, y = freq, fill = cell)) +
  geom_col(width = 0.75) +
  facet_wrap(~ project, ncol = 3) +
  ggsci::scale_fill_igv() +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(
    x = NULL,
    y = NULL
  )

ggsave(
  p_proximity_proj_freq,
  filename = file.path(out_dir, str_glue("{compute_cells} - proximity analysis barplot - project - freq.pdf")),
  width  = 12,
  height = 24
)

# write out cell count of others cells of epithelial
res_pa_count_tidy_epi_all <-
  res_nn_count_df %>%
  mutate(
    project = str_c(group, project, sep = "__"),
    Others  = str_c(cell, radius, sep = "__")
  ) %>%
  arrange(project) %>%
  dplyr::select(project, Others, count) %>%
  pivot_wider(names_from = "project", values_from = "count")


# export 
openxlsx::write.xlsx( 
  res_pa_count_tidy_epi_all, 
  file = file.path(out_dir, str_glue("{compute_cells} - proximity analysis of epithelial in CN9 - all cell counts table.xlsx")), 
  asTable = TRUE 
) 




