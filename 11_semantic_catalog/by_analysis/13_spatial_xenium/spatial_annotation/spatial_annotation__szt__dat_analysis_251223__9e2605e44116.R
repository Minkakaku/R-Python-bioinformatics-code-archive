library(Seurat)
library(tidyverse)
source("spatial analysis.R")

# copy from dat_analysis_250220.R
# 更新细胞注释为最新的结果

# 0. ----------------
out_dir <- "dat_analysis_cycle3/result_20251225"
dir.create(out_dir, recursive = TRUE)

# cn_out_dir <- "dat_analysis_cycle3/cn_analysis_250219"  
# dir.create(cn_out_dir) 


all_data <- readRDS(file.path("dat_analysis_cycle3", "subtype", "all_data_cycle3_250219.rds"))
# all_data_meta_cn <- readRDS(file.path(cn_out_dir, "meta_data_cn.rds"))


# subset samples
group_lym_meta <- c(
  "00037G2",
  "06685N3",
  "06938F1",
  # "09925E1", # only small number of cells, removed
  "09975J1",
  "10858",
  "12876F7",
  "13638",
  "15201H1",
  "28519G1",
  "43158G1",
  "44994A1",
  "45827",
  "47983A3",
  "57986L1",
  "58198",
  "67718G1",
  "71402A6",
  "72696",
  "73479"
  # "77054", # only small number of cells, removed
)

group_non_lym_meta <- c(
  '09791H1',
  '47634I1', #（细胞注释包含大片肿瘤细胞）
  
  '09975N5', #（可能写成09975N5了，基本没切到细胞）
  '06685N2',
  '09925H4', 
  '00037G3',
  '06938I6',
  '28519M4',
  '67718F4',
  '57986I1', # （细胞注释包含大片肿瘤细胞）
  '12876F9'
)

stopifnot( all(group_lym_meta %in% all_data$project) )
stopifnot( all(group_non_lym_meta %in% all_data$project) )


# all_data_subset <- all_data[, all_data$project %in% c(group_lym_meta, group_non_lym_meta)]

# all_data_cn_meta
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- 
      all_data[[]] %>%
      rownames_to_column("cell_id") %>%
      dplyr::select(annotation_full = annotation5, group, project, cell_id) %>%
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
all_data_meta_pos_subset %>% head
# cell_id        X        Y Classification1 group Slide_ID
# 1 1_aaaaeffo-1 1431.054 378.9300            mCAF     A  00037G2
# 2 1_aaaahfio-1 1426.949 384.5132             ASC     A  00037G2
# 3 1_aaaajbca-1 1430.740 389.0237           apCAF     A  00037G2
# 4 1_aaaaohfn-1 1421.059 379.4140             ASC     A  00037G2
# 5 1_aaaaplka-1 1438.761 389.0710            iCAF     A  00037G2
# 6 1_aaabchhk-1 1440.542 395.7985      Macrophage     A  00037G2


# 1. proximity analysis - only group_lym_meta -----------
# epi, rCAF, mCAF

## 1.1 input data ------------
cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

all_data_meta_pos_subset2 <- 
  all_data_meta_pos_subset %>%
  dplyr::filter(Slide_ID %in% group_lym_meta) %>% #  ###
  dplyr::filter(Classification1 %in% c('Epithelial', 'mCAF', 'rCAF')) %>%
  mutate(Classification1 = fct_drop(Classification1))

all_data_meta_pos_subset2_epi <- all_data_meta_pos_subset2[all_data_meta_pos_subset2$Classification1 == "Epithelial", ]

## 1.2 proximity analysis --------------- 
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


saveRDS(res_pa_200um, file = file.path(out_dir, str_glue("epi_caf_res_pa_200um.rds")))


cell_match_var <- structure(
  as.character(all_data_meta_pos_subset2$Classification1),
  names = all_data_meta_pos_subset2$cell_id
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
  unnest(value)  %>%
  # left_join(distinct(all_data_meta_cn_subset[c("project", "group")]), by = "project") %>%
  dplyr::filter(radius != "200-Inf")

## 1.3 plot-----------
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
  filename = file.path(out_dir, str_glue("1 proximity analysis barplot - count.pdf")),
  width  = 8,
  height = 6
)

ggsave(
  p_proximity_freq,
  filename = file.path(out_dir, str_glue("1 proximity analysis barplot - freq.pdf")),
  width  = 8,
  height = 6
)

# plot - project
p_proximity_proj_freq <-
  res_nn_count_df %>%
  group_by(project, radius, cell) %>%
  summarise(count = sum(count)) %>%
  
  group_by(project, radius) %>%
  mutate(freq = count/sum(count)) %>%
  mutate(project = str_c(project, sep = "_")) %>%
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
  filename = file.path(out_dir, str_glue("1 proximity analysis barplot - project - freq.pdf")),
  width  = 12,
  height = 24
)

# write out cell count of others cells of epithelial
res_pa_count_tidy_epi_all <-
  res_nn_count_df %>%
  mutate(
    # project = str_c(group, project, sep = "__"),
    Others  = str_c(cell, radius, sep = "__")
  ) %>%
  arrange(project) %>%
  dplyr::select(project, Others, count) %>%
  pivot_wider(names_from = "project", values_from = "count")


# export 
openxlsx::write.xlsx( 
  res_pa_count_tidy_epi_all, 
  file = file.path(out_dir, str_glue("1 proximity analysis of epithelial - all cell counts table.xlsx")), 
  asTable = TRUE 
) 


## density plot ----------------
p_den_ls <-
  res_nn_df_ls  %>% 
  imap(function(dat, name){
    dat %>%
      pivot_longer(
        c(mCAF, rCAF),
        names_to = 'CAF',
        values_to = "dist"
      ) %>%
      ggplot(aes(x = dist, fill = CAF)) +
      # geom_density(position = "fill")
      geom_density(alpha = 0.5) +
      theme_classic() +
      ggtitle(name)
  })
  

p_den_total <-
  res_nn_df_ls  %>% 
  bind_rows() %>%
  pivot_longer(
    c(mCAF, rCAF),
    names_to = 'CAF',
    values_to = "dist"
  ) %>%
  ggplot(aes(x = dist, fill = CAF)) +
  # geom_density(position = "fill")
  geom_density(alpha = 0.5) +
  theme_classic()

# export
ggsave(
  p_den_total,
  filename = file.path(out_dir, str_glue("1 proximity analysis denstiy plot - all.pdf")),
  width  = 8,
  height = 6
)


pdf(
  file = file.path(out_dir, str_glue("1 proximity analysis denstiy plot - project.pdf")),
  width  = 8,
  height = 6
)
p_den_ls %>% walk(print)
dev.off()


# 2. proximity analysis - all slide -----------
# mCAF

## 2.1 input data ------------
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

mye_cells <- c(
  "Macrophage", 
  "Mast", 
  "Monocyte", 
  "pDCs", 
  "cDC1"
)


b_cells <- c(
  "B-GC",
  "B-memory",
  "B-naive"
)


cell_annotaion_ls <- list(
  CD4_T = CD4_T_cells,
  CD8_T = CD8_T_cells,
  mye   = mye_cells,
  b     = b_cells
)

all_data_meta_pos_subset3 <- 
  all_data_meta_pos_subset %>%
  dplyr::filter(Slide_ID %in% c(group_lym_meta, group_non_lym_meta)) %>% #  ###
  dplyr::filter(Classification1 %in% c('mCAF', CD4_T_cells, CD8_T_cells, mye_cells, b_cells)) %>%
  mutate(Classification1 = fct_drop(Classification1))

all_data_meta_pos_subset3_mCAF <- all_data_meta_pos_subset3[all_data_meta_pos_subset3$Classification1 == "mCAF", ]

## 2.2 proximity analysis --------------- 
res_pa_200um_2 <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset3,
    query = all_data_meta_pos_subset3_mCAF,
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


# saveRDS(res_pa_200um, file = file.path(out_dir, str_glue("mcaf_res_pa_200um.rds")))


cell_match_var_mcaf <- structure(
  as.character(all_data_meta_pos_subset3$Classification1),
  names = all_data_meta_pos_subset3$cell_id
)


res_nn_df_ls_mcaf <-
  res_pa_200um_2 %>%
  # .[1] %>%
  imap(function(nn, grp){ # browser()
    print(grp)
    print(Sys.time())
    res_ls <- list()
    
    nn.idx  <- nn$nn.idx[, -1] + 1
    nn.dist <- nn$nn.dists[, -1]
    
    
    cell_match_var <- factor( c("Others", cell_match_var_mcaf) )
    
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
res_nn_count_df_ls_mcaf <-
  res_nn_df_ls_mcaf %>%
  map(function(x){ # browser()
    xx <- x[, !colnames(x) %in% c("Others", "mCAF")]
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

res_nn_count_df_mcaf <-
  res_nn_count_df_ls_mcaf %>%
  enframe(name = "project") %>%
  unnest(value)  %>%
  # left_join(distinct(all_data_meta_cn_subset[c("project", "group")]), by = "project") %>%
  dplyr::filter(radius != "200-Inf")

## 2.3 plot-----------
# plot  - group
p_proximity_count_mcaf_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_count_df_mcaf %>%
      dplyr::filter(cell %in% x) %>%
      ggplot(aes(x = radius, y = count, fill = cell) ) +
      geom_col() +
      # facet_wrap(~ group) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      labs(
        title = name,
        x = NULL,
        y = NULL
      )
  })



p_proximity_freq_mcaf_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_count_df_mcaf %>%
      dplyr::filter(cell %in% x) %>%
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
        title = name,
        x = NULL,
        y = NULL
      )
  })


# export 
pdf(
  file = file.path(out_dir, str_glue("2 proximity analysis barplot - count.pdf")),
  width  = 8,
  height = 6
)
p_proximity_count_mcaf_ls %>% walk(print)
dev.off()

pdf(
  file = file.path(out_dir, str_glue("2 proximity analysis barplot - freq.pdf")),
  width  = 8,
  height = 6
)
p_proximity_freq_mcaf_ls %>% walk(print)
dev.off()


# plot - project
p_proximity_proj_freq_mcaf_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_count_df_mcaf %>%
      dplyr::filter(cell %in% x) %>%
      group_by(project, radius, cell) %>%
      summarise(count = sum(count)) %>%
      
      group_by(project, radius) %>%
      mutate(freq = count/sum(count)) %>%
      mutate(project = str_c(project, sep = "_")) %>%
      mutate(project = factor(project, levels = sort(unique(project)))) %>%
      
      ggplot(aes(x = radius, y = freq, fill = cell)) +
      geom_col(width = 0.75) +
      facet_wrap(~ project, ncol = 3) +
      ggsci::scale_fill_igv() +
      scale_y_continuous(labels = scales::percent) +
      theme_bw() +
      labs(
        title = name,
        x = NULL,
        y = NULL
      )
  })


pdf(
  file = file.path(out_dir, str_glue("2 proximity analysis barplot - project - freq.pdf")),
  width  = 12,
  height = 24
)
p_proximity_proj_freq_mcaf_ls %>% walk(print)
dev.off()


# write out cell count of others cells of epithelial
res_pa_count_tidy_mcaf_all <-
  res_nn_count_df_mcaf %>%
  mutate(
    # project = str_c(group, project, sep = "__"),
    Others  = str_c(cell, radius, sep = "__")
  ) %>%
  arrange(project) %>%
  dplyr::select(project, Others, count) %>%
  pivot_wider(names_from = "project", values_from = "count")


# export 
openxlsx::write.xlsx( 
  res_pa_count_tidy_mcaf_all, 
  file = file.path(out_dir, str_glue("2 proximity analysis of mCAF - all cell counts table.xlsx")), 
  asTable = TRUE 
) 


## density plot ----------------
p_den_ls_mcaf_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_df_ls_mcaf  %>% 
      imap(function(dat, name){
        dat2 <- dat %>%
          pivot_longer(
            all_of(x),
            names_to = 'CAF',
            values_to = "dist"
          )
        dat2$dist[is.na(dat2$dist)] <- 1000
          
        dat2 %>%
          ggplot(aes(x = dist, fill = CAF)) +
          # geom_density(position = "fill")
          geom_density(alpha = 0.5) +
          theme_classic() +
          ggtitle(name)
      })
  })
  


p_den_total_mcaf_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    dat <- res_nn_df_ls_mcaf  %>% 
      bind_rows() %>%
      pivot_longer(
        all_of(x),
        names_to = 'CAF',
        values_to = "dist"
      )  
    dat$dist[is.na(dat$dist)] <- 1000
    
    dat %>%
      ggplot(aes(x = dist, fill = CAF)) +
      # geom_density(position = "fill")
      geom_density(alpha = 0.5) +
      theme_classic() +
      ggtitle(name)
  })
  

# export
pdf(
  file = file.path(out_dir, str_glue("2 proximity analysis denstiy plot - all.pdf")),
  width  = 8,
  height = 6
)
p_den_total_mcaf_ls %>% walk(print)
dev.off()

p_den_ls_mcaf_ls %>%
  iwalk(function(p_ls, name){
    pdf(
      file = file.path(out_dir, str_glue("2 proximity analysis denstiy plot - project - {name}.pdf")),
      width  = 8,
      height = 6
    )
    p_ls %>% walk(print)
    dev.off()
    
  })


# 3. proximity analysis - all slide -----------
# rCAF

## 3.1 input data ------------
# cell_annotaion_ls <- list(
#   CD4_T = CD4_T_cells,
#   CD8_T = CD8_T_cells,
#   mye   = mye_cells,
#   b     = b_cells
# )

all_data_meta_pos_subset4 <- 
  all_data_meta_pos_subset %>%
  dplyr::filter(Slide_ID %in% c(group_lym_meta, group_non_lym_meta)) %>% #  ###
  dplyr::filter(Classification1 %in% c('rCAF', CD4_T_cells, CD8_T_cells, mye_cells, b_cells)) %>%
  mutate(Classification1 = fct_drop(Classification1))

all_data_meta_pos_subset4_rCAF <- all_data_meta_pos_subset4[all_data_meta_pos_subset4$Classification1 == "rCAF", ]

## 3.2 proximity analysis --------------- 
res_pa_200um_3 <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset4,
    query = all_data_meta_pos_subset4_rCAF,
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


# saveRDS(res_pa_200um, file = file.path(out_dir, str_glue("mcaf_res_pa_200um.rds")))


cell_match_var_rcaf <- structure(
  as.character(all_data_meta_pos_subset4$Classification1),
  names = all_data_meta_pos_subset4$cell_id
)


res_nn_df_ls_rcaf <-
  res_pa_200um_3 %>%
  # .[1] %>%
  imap(function(nn, grp){ # browser()
    print(grp)
    print(Sys.time())
    res_ls <- list()
    
    nn.idx  <- nn$nn.idx[, -1] + 1
    nn.dist <- nn$nn.dists[, -1]
    
    
    cell_match_var <- factor( c("Others", cell_match_var_rcaf) )
    
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
res_nn_count_df_ls_rcaf <-
  res_nn_df_ls_rcaf %>%
  map(function(x){ # browser()
    xx <- x[, !colnames(x) %in% c("Others", "rCAF")]
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

res_nn_count_df_rcaf <-
  res_nn_count_df_ls_rcaf %>%
  enframe(name = "project") %>%
  unnest(value)  %>%
  # left_join(distinct(all_data_meta_cn_subset[c("project", "group")]), by = "project") %>%
  dplyr::filter(radius != "200-Inf")

## 3.3 plot-----------
# plot  - group
p_proximity_count_rcaf_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_count_df_rcaf %>%
      dplyr::filter(cell %in% x) %>%
      ggplot(aes(x = radius, y = count, fill = cell) ) +
      geom_col() +
      # facet_wrap(~ group) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      labs(
        title = name,
        x = NULL,
        y = NULL
      )
  })



p_proximity_freq_rcaf_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_count_df_rcaf %>%
      dplyr::filter(cell %in% x) %>%
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
        title = name,
        x = NULL,
        y = NULL
      )
  })


# export 
pdf(
  file = file.path(out_dir, str_glue("3 proximity analysis barplot - count.pdf")),
  width  = 8,
  height = 6
)
p_proximity_count_rcaf_ls %>% walk(print)
dev.off()

pdf(
  file = file.path(out_dir, str_glue("3 proximity analysis barplot - freq.pdf")),
  width  = 8,
  height = 6
)
p_proximity_freq_rcaf_ls %>% walk(print)
dev.off()


# plot - project
p_proximity_proj_freq_rcaf_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_count_df_rcaf %>%
      dplyr::filter(cell %in% x) %>%
      group_by(project, radius, cell) %>%
      summarise(count = sum(count)) %>%
      
      group_by(project, radius) %>%
      mutate(freq = count/sum(count)) %>%
      mutate(project = str_c(project, sep = "_")) %>%
      mutate(project = factor(project, levels = sort(unique(project)))) %>%
      
      ggplot(aes(x = radius, y = freq, fill = cell)) +
      geom_col(width = 0.75) +
      facet_wrap(~ project, ncol = 3) +
      ggsci::scale_fill_igv() +
      scale_y_continuous(labels = scales::percent) +
      theme_bw() +
      labs(
        title = name,
        x = NULL,
        y = NULL
      )
  })


pdf(
  file = file.path(out_dir, str_glue("3 proximity analysis barplot - project - freq.pdf")),
  width  = 12,
  height = 24
)
p_proximity_proj_freq_rcaf_ls %>% walk(print)
dev.off()


# write out cell count of others cells of epithelial
res_pa_count_tidy_rcaf_all <-
  res_nn_count_df_rcaf %>%
  mutate(
    # project = str_c(group, project, sep = "__"),
    Others  = str_c(cell, radius, sep = "__")
  ) %>%
  arrange(project) %>%
  dplyr::select(project, Others, count) %>%
  pivot_wider(names_from = "project", values_from = "count")


# export 
openxlsx::write.xlsx( 
  res_pa_count_tidy_rcaf_all, 
  file = file.path(out_dir, str_glue("3 proximity analysis of rCAF - all cell counts table.xlsx")), 
  asTable = TRUE 
) 


## density plot ----------------
p_den_ls_rcaf_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_df_ls_rcaf  %>% 
      imap(function(dat, name){
        dat2 <- dat %>%
          pivot_longer(
            all_of(x),
            names_to = 'CAF',
            values_to = "dist"
          )
        dat2$dist[is.na(dat2$dist)] <- 1000
        
        dat2 %>%
          ggplot(aes(x = dist, fill = CAF)) +
          # geom_density(position = "fill")
          geom_density(alpha = 0.5) +
          theme_classic() +
          ggtitle(name)
      })
  })



p_den_total_rcaf_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    dat <- res_nn_df_ls_rcaf  %>% 
      bind_rows() %>%
      pivot_longer(
        all_of(x),
        names_to = 'CAF',
        values_to = "dist"
      )  
    dat$dist[is.na(dat$dist)] <- 1000
    
    dat %>%
      ggplot(aes(x = dist, fill = CAF)) +
      # geom_density(position = "fill")
      geom_density(alpha = 0.5) +
      theme_classic() +
      ggtitle(name)
  })


# export
pdf(
  file = file.path(out_dir, str_glue("3 proximity analysis denstiy plot - all.pdf")),
  width  = 8,
  height = 6
)
p_den_total_rcaf_ls %>% walk(print)
dev.off()

p_den_ls_rcaf_ls %>%
  iwalk(function(p_ls, name){
    pdf(
      file = file.path(out_dir, str_glue("3 proximity analysis denstiy plot - project - {name}.pdf")),
      width  = 8,
      height = 6
    )
    p_ls %>% walk(print)
    dev.off()
    
  })










