library(Seurat)
library(tidyverse)
source("spatial analysis.R")

# copy from dat_analysis_251223.R
# 更新细胞注释为最新的结果

# 0. ----------------
out_dir <- "dat_analysis_cycle3/result_20260308"
dir.create(out_dir, recursive = TRUE)

all_data <- readRDS(file.path("dat_analysis_cycle3", "subtype", "all_data_cycle3_250219.rds"))

cn_out_dir <- "dat_analysis_cycle3/cn_analysis_250219"  
all_data_meta_cn <- readRDS(file.path(cn_out_dir, "meta_data_cn.rds"))


# subset samples
removed_sample <- c(
  '09925E1',
  '77054',
  '09975N5'
)

group_primary <- c(
  '09791H1',
  '47634I1' 
)

group_lym_meta <- c(
  '00037G2',
  '06685N2',
  '06938F1',
  '09975J1',
  '10858',
  '12876F7',
  '13638',
  '15201H1',
  '28519G1',
  '43158G1',
  '44994A1',
  '45827',
  '47983A3',
  '57986L1',
  '58198',
  '67718G1',
  '71402A6',
  '72696',
  '73479'
)

group_non_lym_meta <- c(
  '06685N3',
  '09925H4',
  '00037G3',
  '06938I6',
  '28519M4',
  '67718F4',
  '57986I1',
  '12876F9'
)

stopifnot( all(removed_sample %in% all_data$project) )
stopifnot( all(group_primary %in% all_data$project) )
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


# 1. proximity analysis - all slide -----------
# CD8-Tex : 'all cells' or 'all CAFs'
# only Metastatic samples

## 1.1 input data ------------
caf_cells <- c(
  "apCAF",
  "iCAF",
  "mCAF",
  "rCAF",
  "vCAF"
)

cell_annotaion_ls <- list(
  CAF = caf_cells,
  all_cells = as.character(unique(all_data_meta_pos_subset$Classification1))
)


all_data_meta_pos_subset5 <- 
  all_data_meta_pos_subset %>%
  dplyr::filter(Slide_ID %in% c( group_lym_meta, group_non_lym_meta)) %>% #  ###
  # mutate(slide_group = case_when(
  #   Slide_ID %in% group_primary ~ "Primary",
  #   Slide_ID %in% group_lym_meta ~ "Metastatic ",
  #   Slide_ID %in% group_non_lym_meta ~ "Metastatic ",
  #   TRUE ~ "Others"
  # )) %>%
  
  # dplyr::filter(Classification1 %in% c('CD8-Tex')) %>%
  mutate(Classification1 = fct_drop(Classification1)) 

all_data_meta_pos_subset5_cd8tex <- all_data_meta_pos_subset5[all_data_meta_pos_subset5$Classification1 == "CD8-Tex", ]

## 1.2 proximity analysis --------------- 
res_pa_200um_5 <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset5,
    query = all_data_meta_pos_subset5_cd8tex,
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
    cell_enough_ratio_cutoff = 99
  )


# saveRDS(res_pa_200um, file = file.path(out_dir, str_glue("mcaf_res_pa_200um.rds")))


cell_match_var_cd8tex <- structure(
  as.character(all_data_meta_pos_subset5$Classification1),
  names = all_data_meta_pos_subset5$cell_id
)


res_nn_df_ls_cd8tex <-
  res_pa_200um_5 %>%
  # .[1] %>%
  imap(function(nn, grp){ # browser()
    print(grp)
    print(Sys.time())
    res_ls <- list()
    
    nn.idx  <- nn$nn.idx[, -1] + 1
    nn.dist <- nn$nn.dists[, -1]
    
    
    cell_match_var <- factor( c("Others", cell_match_var_cd8tex) ) # <------
    
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
res_nn_count_df_ls_cd8tex <-
  res_nn_df_ls_cd8tex %>%
  map(function(x){ # browser()
    xx <- x[, !colnames(x) %in% c("Others", "CD8-Tex")]
    xx[is.na(xx)] <- Inf
    
    dat_ls <-
      apply(xx, 2, function(y){
        x_cut <-
          cut(
            y, 
            breaks = c(0, 50, 100, 150, 200, Inf),
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

res_nn_count_df_cd8tex <-
  res_nn_count_df_ls_cd8tex %>%
  enframe(name = "project") %>%
  unnest(value)  %>%
  # left_join(distinct(all_data_meta_cn_subset[c("project", "group")]), by = "project") %>%
  dplyr::filter(radius != "200-Inf")

## 1.3 plot -----------
# plot  - group
p_proximity_count_cd8tex_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_count_df_cd8tex %>%
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


p_proximity_freq_cd8tex_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_count_df_cd8tex %>%
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
  file = file.path(out_dir, str_glue("1 proximity analysis barplot of CD8-Tex - count - all metastatic slide.pdf")),
  width  = 6,
  height = 8
)
p_proximity_count_cd8tex_ls %>% walk(print)
dev.off()

pdf(
  file = file.path(out_dir, str_glue("1 proximity analysis barplot of CD8-Tex - freq - all metastatic slide.pdf")),
  width  = 6,
  height = 8
)
p_proximity_freq_cd8tex_ls %>% walk(print)
dev.off()



# plot - project
p_proximity_proj_count_cd8tex_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_count_df_cd8tex %>%
      dplyr::filter(cell %in% x) %>%
      ggplot(aes(x = radius, y = count, fill = cell) ) +
      geom_col() +
      facet_wrap(~ project, ncol = 3) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      labs(
        title = name,
        x = NULL,
        y = NULL
      )
  })

p_proximity_proj_freq_cd8tex_ls <-
  cell_annotaion_ls %>%
  imap(function(x, name){
    res_nn_count_df_cd8tex %>%
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
  file = file.path(out_dir, str_glue("2 proximity analysis barplot of CD8-Tex - project - count.pdf")),
  width  = 12,
  height = 24
)
p_proximity_proj_count_cd8tex_ls %>% walk(print)
dev.off()

pdf(
  file = file.path(out_dir, str_glue("2 proximity analysis barplot of CD8-Tex - project - freq.pdf")),
  width  = 12,
  height = 24
)
p_proximity_proj_freq_cd8tex_ls %>% walk(print)
dev.off()


# write out cell count of others cells of epithelial
res_pa_count_tidy_cd8tex_all <-
  res_nn_count_df_cd8tex %>%
  mutate(
    # project = str_c(group, project, sep = "__"),
    Others  = str_c(cell, radius, sep = "__")
  ) %>%
  arrange(project) %>%
  dplyr::select(project, Others, count) %>%
  pivot_wider(names_from = "project", values_from = "count")


# export 
openxlsx::write.xlsx( 
  res_pa_count_tidy_cd8tex_all, 
  file = file.path(out_dir, str_glue("2 proximity analysis of CD8-Tex - all cell counts table.xlsx")), 
  asTable = TRUE 
) 


# 2. promixity analysis - per CN -----------

all_data_meta_cn %>%
  group_by(CN) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  iwalk(function(dat, cn_name){
    dir.create(file.path(out_dir, 'CN', cn_name), recursive = TRUE)
    
    
    ## 2.1 data input -----------
    all_data_meta_pos_subset5_cn <- 
      all_data_meta_pos_subset %>%
      dplyr::filter(Slide_ID %in% c( group_lym_meta, group_non_lym_meta)) %>% #  ###
      mutate(Classification1 = fct_drop(Classification1)) %>%
      dplyr::filter(cell_id %in% dat$cell_id) # per CN
    
    all_data_meta_pos_subset5_cd8tex_cn <- all_data_meta_pos_subset5_cn[all_data_meta_pos_subset5_cn$Classification1 == "CD8-Tex", ]
    
    
    ## 2.2 proximity analysis --------------- 
    res_pa_200um_5_cn <-
      ProximutyAnalysis_NN_multi(
        all_data_meta_pos_subset5_cn,
        query = all_data_meta_pos_subset5_cd8tex_cn,
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
        cell_enough_ratio_cutoff = 99
      )
    
    
    # saveRDS(res_pa_200um, file = file.path(out_dir, str_glue("mcaf_res_pa_200um.rds")))
    
    
    cell_match_var_cd8tex_cn <- structure(
      as.character(all_data_meta_pos_subset5_cn$Classification1),
      names = all_data_meta_pos_subset5_cn$cell_id
    )
    
    
    res_nn_df_ls_cd8tex_cn <-
      res_pa_200um_5_cn %>%
      # .[1] %>%
      imap(function(nn, grp){ # browser()
        print(grp)
        print(Sys.time())
        res_ls <- list()
        
        if(length(dim(nn$nn.idx)) < 2) return(NULL)
        
        nn.idx  <- nn$nn.idx[, -1] + 1
        nn.dist <- nn$nn.dists[, -1]
        
        
        cell_match_var <- factor( c("Others", cell_match_var_cd8tex_cn) ) # <------
        
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
        
      }) %>%
      discard(is.null)
    
    
    # cell count prop 
    res_nn_count_df_ls_cd8tex_cn <-
      res_nn_df_ls_cd8tex_cn %>%
      map(function(x){ # browser()
        xx <- x[, !colnames(x) %in% c("Others", "CD8-Tex")]
        xx[is.na(xx)] <- Inf
        
        dat_ls <-
          apply(xx, 2, function(y){
            x_cut <-
              cut(
                y, 
                breaks = c(0, 50, 100, 150, 200, Inf),
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
    
    res_nn_count_df_cd8tex_cn <-
      res_nn_count_df_ls_cd8tex_cn %>%
      enframe(name = "project") %>%
      unnest(value)  %>%
      # left_join(distinct(all_data_meta_cn_subset[c("project", "group")]), by = "project") %>%
      dplyr::filter(radius != "200-Inf")
    
    ## 2.3 plot -----------
    # plot  - group
    p_proximity_count_cd8tex_ls_cn <-
      cell_annotaion_ls %>%
      imap(function(x, name){
        res_nn_count_df_cd8tex_cn %>%
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
    
    
    
    p_proximity_freq_cd8tex_ls_cn <-
      cell_annotaion_ls %>%
      imap(function(x, name){
        res_nn_count_df_cd8tex_cn %>%
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
      file = file.path(out_dir, 'CN', cn_name, str_glue("1 proximity analysis barplot of CD8-Tex - count - all metastatic slide.pdf")),
      width  = 6,
      height = 8
    )
    p_proximity_count_cd8tex_ls_cn %>% walk(print)
    dev.off()
    
    pdf(
      file = file.path(out_dir, 'CN', cn_name, str_glue("1 proximity analysis barplot of CD8-Tex - freq - all metastatic slide.pdf")),
      width  = 6,
      height = 8
    )
    p_proximity_freq_cd8tex_ls_cn %>% walk(print)
    dev.off()
    
    
    
    # plot - project
    p_proximity_proj_count_cd8tex_ls_cn <-
      cell_annotaion_ls %>%
      imap(function(x, name){
        res_nn_count_df_cd8tex_cn %>%
          dplyr::filter(cell %in% x) %>%
          ggplot(aes(x = radius, y = count, fill = cell) ) +
          geom_col() +
          facet_wrap(~ project, ncol = 3) +
          ggsci::scale_fill_igv() +
          theme_bw() +
          labs(
            title = name,
            x = NULL,
            y = NULL
          )
      })
    
    p_proximity_proj_freq_cd8tex_ls_cn <-
      cell_annotaion_ls %>%
      imap(function(x, name){
        res_nn_count_df_cd8tex_cn %>%
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
      file = file.path(out_dir, 'CN', cn_name, str_glue("2 proximity analysis barplot of CD8-Tex - project - count.pdf")),
      width  = 12,
      height = 24
    )
    p_proximity_proj_count_cd8tex_ls_cn %>% walk(print)
    dev.off()
    
    pdf(
      file = file.path(out_dir, 'CN', cn_name, str_glue("2 proximity analysis barplot of CD8-Tex - project - freq.pdf")),
      width  = 12,
      height = 24
    )
    p_proximity_proj_freq_cd8tex_ls_cn %>% walk(print)
    dev.off()
    
    
    # write out cell count of others cells of epithelial
    res_pa_count_tidy_cd8tex_all_cn <-
      res_nn_count_df_cd8tex_cn %>%
      mutate(
        # project = str_c(group, project, sep = "__"),
        Others  = str_c(cell, radius, sep = "__")
      ) %>%
      arrange(project) %>%
      dplyr::select(project, Others, count) %>%
      pivot_wider(names_from = "project", values_from = "count")
    
    
    # export 
    openxlsx::write.xlsx( 
      res_pa_count_tidy_cd8tex_all_cn, 
      file = file.path(out_dir, 'CN', cn_name, str_glue("2 proximity analysis of CD8-Tex - all cell counts table.xlsx")), 
      asTable = TRUE 
    ) 
    
  })



