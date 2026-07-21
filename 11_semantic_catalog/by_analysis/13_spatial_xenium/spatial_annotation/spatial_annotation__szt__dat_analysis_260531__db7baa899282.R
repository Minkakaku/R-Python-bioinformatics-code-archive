library(tidyverse)
library(Seurat)
source('spatial analysis.R')

out_dir_260531 <- "dat_analysis_cycle3/result_2026531"
dir.create(out_dir_260531)

# read data --------------
cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")
all_data <- readRDS("dat_analysis_cycle3/subtype/all_data_cycle3_260417.rds")
CD8T <- read("dat_analysis_cycle3/subtype/CD8T/CD8_subtype2.rds")

# subset samples ------------------
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

length(c(group_primary, group_lym_meta, group_non_lym_meta) )

# cold and hot tumor ------------------------
group_non_meta <- c(
  '09791H1', '47634I1', '06685N3', '09925H4', '00037G3',
  '06938I6', '28519M4', '67718F4', '57986I1', '12876F9'
)
group_lym_meta_cold <- c(
  '73479', '13638', '09975J1', '06685N2', 
  # '09925E1', '77054', # removed
  '00037G2',    '72696', '45827', '28519G1',
  '71402A6', '43158G1', '67718G1'
)
group_lym_meta_hot <- c(
  '47983A3', '10858', '58198', '44994A1', '06938F1', '15201H1', '57986L1', '12876F7'
)

# 	无转移淋巴结：09791H1、47634I1（细胞注释包含大片肿瘤细胞）；09775N5（可能写成09975N5了，基本没切到细胞）、06685N3、00925H4、00037G3、06938I6、28519M4、67718F4、57986I1（细胞注释包含大片肿瘤细胞）、12876F9
# 	冷肿瘤转移淋巴结（luminal）：73479、13638、09775J1、06685N2、00925E1、00037G2、77054、72696、45827、28519G1、71402A6、43158G1、67718G1
# 	热肿瘤转移淋巴结（TNBC、Her2）：47983A3、10858、58198、44994A1、06938F1、15201H1、57986L1、12876F7

stopifnot( all(group_non_meta %in% all_data$project) )
stopifnot( all(group_lym_meta_cold %in% all_data$project) )
stopifnot( all(group_lym_meta_hot %in% all_data$project) )

# length(c(group_non_meta, group_lym_meta_cold, group_lym_meta_hot))

# cells ------------------
caf_cells <- c(
  "apCAF",
  "iCAF",
  "mCAF",
  "rCAF",
  "vCAF"
)

CD4_T_cells <- c(
  'Treg',
  'Tfh',
  'CD4-Tnaive',
  'CD4-Tmemory',
  'CD4T-Tcm'
)

CD8_T_cells <- c(
  'CD8-Tem',
  'CD8-Tpex',
  'CD8-Tex',
  'CD8-Tnaive',
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
  CD4_T     = CD4_T_cells,
  CD8_T     = CD8_T_cells,
  mye       = mye_cells,
  b         = b_cells,
  caf_cells = caf_cells
)

# all_data_meta_pos -----------
all_data_meta_pos <-
  cells_position_ls %>%
  imap(function(cells_position, slide){ # browser()
    
    all_data_meta <- 
      all_data[[]] %>%
      rownames_to_column("cell_id") %>%
      dplyr::select(annotation_full = annotation6, group, project, cell_id) %>%
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



# 1. proximity analysis of each CD8T subtype ---------------
all_data_meta_pos_subset1 <- 
  all_data_meta_pos_subset %>%
  dplyr::filter(! Slide_ID %in% removed_sample) %>% #  ###
  dplyr::filter(Classification1 %in% purrr::list_c(cell_annotaion_ls)) %>%
  mutate(Classification1 = fct_drop(Classification1)) 

res_pa_200um_1 <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset1,
    query = all_data_meta_pos_subset1,
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

# saveRDS(res_pa_200um_1, file = file.path(out_dir_260531, str_glue("res_pa_200um_1.rds")))

cell_match_var_ls <- 
  all_data_meta_pos_subset1 %>% 
  with(split(Classification1, Slide_ID)) %>%
  map(as.character)

# CD8_T_cells
center_cell_id <- all_data_meta_pos_subset1$cell_id[all_data_meta_pos_subset1$Classification1 %in% CD8_T_cells]
  
## compute mean distance  ----------------
res_nn_df_ls_cd8t <-
  res_pa_200um_1 %>%
  # .[1] %>%
  imap(function(nn, grp){ # browser()
    print(grp)
    print(Sys.time())
    res_ls <- list()
    
    cell_match_var <- cell_match_var_ls[[grp]]
    
    nn.idx  <- nn$nn.idx[, -1] + 1
    nn.dist <- nn$nn.dists[, -1]
    total.cell <- nn$total_cell_id
    
    cell_match_var <- factor( c("Others", cell_match_var) ) #### ###
    
    # only compute cencter cells
    center_idx <- nn$query_cell_id %in% center_cell_id
    
    for(i in 1:nrow(nn$nn.idx[center_idx, ])){ 
      if(i %% 5000 == 0) print(i)
      
      # head 10 cells, mean dist
      res <- tapply(nn.dist[i, ], cell_match_var[nn.idx[i, ]], function(x){mean(head(x, 10))})
      res_ls[[i]] <- res
    }
    
    print(Sys.time())
    
    res_df <- as.data.frame(bind_rows(res_ls))
    rownames(res_df) <- nn$query_cell_id[center_idx]
    res_df
    # res_ls
    
  })


# cell counts  
res_nn_count_df_ls_cd8t_center_ls <-
  CD8_T_cells %>%
  setNames(., .) %>%
  map(function(cell){
    
    cell_annotaion_ls %>%
      .[! names(.) %in% 'CD8_T'] %>%
      imap(function(cells, name){
        
        res_nn_count_df_ls <-
          res_nn_df_ls_cd8t %>%
          imap(function(x, proj){ # browser()
            cell_match_var <- cell_match_var_ls[[proj]]
            center_cell_idx <- cell_match_var == cell
            if(!any(center_cell_idx)){
              print(str_glue("cencter {cell} in {proj} have no cells and is droped."))
              return(NULL)
            }
            
            xx <- x[center_cell_idx, cells]
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
          }) %>%
          purrr::discard(is.null)
        
        res_nn_count_df <-
          res_nn_count_df_ls %>%
          enframe(name = "project") %>%
          unnest(value)  %>%
          dplyr::filter(radius != "200-Inf")
        
        res_nn_count_df
      })
  })

## plot -----------
### 1 plot  - group ------------
p_proximity_count_cd8t_center_ls_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
          ggplot(aes(x = radius, y = count, fill = cell) ) +
          geom_col() +
          ggsci::scale_fill_igv() +
          theme_bw() +
          labs(
            title = name,
            x = NULL,
            y = NULL
          )
      })
  })


p_proximity_freq_cd8t_center_ls_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
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
  })


# export 
p_proximity_count_cd8t_center_ls_ls %>%
  iwalk(function(p_ls, cd8_cell){
    pdf(
      file = file.path(
        out_dir_260531, 
        str_glue("1.1 proximity analysis barplot of {cd8_cell} - count - all slide.pdf")
      ),
      width  = 6,
      height = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })


p_proximity_freq_cd8t_center_ls_ls %>%
  iwalk(function(p_ls, cd8_cell){
    pdf(
      file = file.path(
        out_dir_260531,
        str_glue("1.1 proximity analysis barplot of {cd8_cell} - freq - all slide.pdf")),
      width  = 6,
      height = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })


### 2 plot  - group - lym_meta ------------
lym_meta_ls <- list(
  group_primary      = group_primary,
  group_lym_meta     = group_lym_meta,
  group_non_lym_meta = group_non_lym_meta
)

p_proximity_count_cd8t_center_ls_subgrp_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    lym_meta_ls %>%
      imap(function(proj_ls, grp_name){
        ###
        nn_df_ls %>%
          imap(function(nn_df, name){
            nn_df %>%
              
              dplyr::filter(project %in% proj_ls) %>% # filter by sub group
              
              ggplot(aes(x = radius, y = count, fill = cell) ) +
              geom_col() +
              ggsci::scale_fill_igv() +
              theme_bw() +
              labs(
                title = str_glue("{grp_name} - {name}"),
                x = NULL,
                y = NULL
              )
          })
        ###
      })
  })


p_proximity_freq_cd8t_center_ls_subgrp_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    lym_meta_ls %>%
      imap(function(proj_ls, grp_name){
        ###
        nn_df_ls %>%
          imap(function(nn_df, name){
            nn_df %>%
              
              dplyr::filter(project %in% proj_ls) %>% # filter by sub group
              
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
                title = str_glue("{grp_name} - {name}"),
                x = NULL,
                y = NULL
              )
          })
        ###
      })
  })


# export 
p_proximity_count_cd8t_center_ls_subgrp_ls %>%
  iwalk(function(p_ls_ls, cd8_cell){
    p_ls_ls %>%
      iwalk(function(p_ls, sub_grp){
        pdf(
          file = file.path(
            out_dir_260531, 
            str_glue("1.2 proximity analysis barplot of {cd8_cell} - count - {sub_grp}.pdf")
          ),
          width  = 6,
          height = 8
        )
        p_ls %>% walk(print)
        dev.off()
      })
  })


p_proximity_freq_cd8t_center_ls_subgrp_ls %>%
  iwalk(function(p_ls_ls, cd8_cell){
    p_ls_ls %>%
      iwalk(function(p_ls, sub_grp){
        pdf(
          file = file.path(
            out_dir_260531,
            str_glue("1.2 proximity analysis barplot of {cd8_cell} - freq - {sub_grp}.pdf")),
          width  = 6,
          height = 8
        )
        p_ls %>% walk(print)
        dev.off()
      })
  })

### 3 plot  - group - cold and hot ------------
lym_cold_hot_ls <- list(
  group_non_meta      = group_non_meta,
  group_lym_meta_cold = group_lym_meta_cold,
  group_lym_meta_hot  = group_lym_meta_hot
)

p_proximity_count_cd8t_center_ls_subgrp2_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    lym_cold_hot_ls %>%
      imap(function(proj_ls, grp_name){
        ###
        nn_df_ls %>%
          imap(function(nn_df, name){
            nn_df %>%
              
              dplyr::filter(project %in% proj_ls) %>% # filter by sub group
              
              ggplot(aes(x = radius, y = count, fill = cell) ) +
              geom_col() +
              ggsci::scale_fill_igv() +
              theme_bw() +
              labs(
                title = str_glue("{grp_name} - {name}"),
                x = NULL,
                y = NULL
              )
          })
        ###
      })
  })


p_proximity_freq_cd8t_center_ls_subgrp2_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    lym_cold_hot_ls %>%
      imap(function(proj_ls, grp_name){
        ###
        nn_df_ls %>%
          imap(function(nn_df, name){
            nn_df %>%
              
              dplyr::filter(project %in% proj_ls) %>% # filter by sub group
              
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
                title = str_glue("{grp_name} - {name}"),
                x = NULL,
                y = NULL
              )
          })
        ###
      })
  })


# export 
p_proximity_count_cd8t_center_ls_subgrp2_ls %>%
  iwalk(function(p_ls_ls, cd8_cell){
    p_ls_ls %>%
      iwalk(function(p_ls, sub_grp){
        pdf(
          file = file.path(
            out_dir_260531, 
            str_glue("1.3 proximity analysis barplot of {cd8_cell} - count - {sub_grp}.pdf")
          ),
          width  = 6,
          height = 8
        )
        p_ls %>% walk(print)
        dev.off()
      })
  })


p_proximity_freq_cd8t_center_ls_subgrp2_ls %>%
  iwalk(function(p_ls_ls, cd8_cell){
    p_ls_ls %>%
      iwalk(function(p_ls, sub_grp){
        pdf(
          file = file.path(
            out_dir_260531,
            str_glue("1.3 proximity analysis barplot of {cd8_cell} - freq - {sub_grp}.pdf")),
          width  = 6,
          height = 8
        )
        p_ls %>% walk(print)
        dev.off()
      })
  })


### 4 plot - project --------------
p_proximity_count_cd8t_center_ls_proj_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
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
  })

p_proximity_proj_freq_cd8t_center_ls_proj_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
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
  })


# export
p_proximity_count_cd8t_center_ls_proj_ls %>%
  iwalk(function(p_ls, cd8_cell){
    pdf(
      file = file.path(
        out_dir_260531,
        str_glue("1.4 proximity analysis barplot of {cd8_cell} - project - count.pdf")),
      width  = 12,
      height = 24
    )
    p_ls %>% walk(print)
    dev.off()
  })


p_proximity_proj_freq_cd8t_center_ls_proj_ls %>%
  iwalk(function(p_ls, cd8_cell){
    pdf(
      file = file.path(
        out_dir_260531,
        str_glue("1.4 proximity analysis barplot of {cd8_cell} - project - freq.pdf")),
      width  = 12,
      height = 24
    )
    p_ls %>% walk(print)
    dev.off()
  })

# write out cell count of others cells of epithelial
res_pa_count_tidy_cd8t_center_all <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  map(function(nn_df_ls){
    nn_df_ls %>%
      imap(function(nn_df, major){
        nn_df %>%
        mutate(
          Others = str_c(.env[['major']], cell, radius, sep = "__")
        ) %>%
        arrange(project) %>%
        dplyr::select(project, Others, count) %>%
        pivot_wider(names_from = "project", values_from = "count")
    }) %>%
    bind_rows()
  })
  


# export 
openxlsx::write.xlsx( 
  res_pa_count_tidy_cd8t_center_all, 
  file = file.path(out_dir_260531, str_glue("1 proximity analysis of each cluster of CD8-T - all cell counts table.xlsx")), 
  asTable = TRUE 
) 



# 2. proximity analysis of each CAF subtype ---------------
# caf_cells
center_cell_id_caf <- 
  all_data_meta_pos_subset1$cell_id[
    all_data_meta_pos_subset1$Classification1 %in% caf_cells
  ]

## compute mean distance -------------
res_nn_df_ls_caf <-
  res_pa_200um_1 %>%
  # .[1] %>%
  imap(function(nn, grp){ # browser()
    print(grp)
    print(Sys.time())
    res_ls <- list()
    
    cell_match_var <- cell_match_var_ls[[grp]]
    
    nn.idx  <- nn$nn.idx[, -1] + 1
    nn.dist <- nn$nn.dists[, -1]
    total.cell <- nn$total_cell_id
    
    cell_match_var <- factor( c("Others", cell_match_var) ) #### ###
    
    # only compute cencter cells
    center_idx <- nn$query_cell_id %in% center_cell_id_caf
    
    for(i in 1:nrow(nn$nn.idx[center_idx, ])){ 
      if(i %% 5000 == 0) print(i)
      
      # head 10 cells, mean dist
      res <- tapply(nn.dist[i, ], cell_match_var[nn.idx[i, ]], function(x){mean(head(x, 10))})
      res_ls[[i]] <- res
    }
    
    print(Sys.time())
    
    res_df <- as.data.frame(bind_rows(res_ls))
    rownames(res_df) <- nn$query_cell_id[center_idx]
    res_df
    # res_ls
    
  })


# cell counts  
res_nn_count_df_ls_caf_center_ls <-
  caf_cells %>%
  setNames(., .) %>%
  map(function(cell){
    
    cell_annotaion_ls %>%
      .[! names(.) %in% 'caf_cells'] %>%
      imap(function(cells, name){
        
        res_nn_count_df_ls <-
          res_nn_df_ls_caf %>%
          imap(function(x, proj){ # browser()
            
            cell_match_var <- cell_match_var_ls[[proj]]
            center_cell_idx <- cell_match_var == cell
            if(!any(center_cell_idx)){
              print(str_glue("cencter {cell} in {proj} have no cells and is droped."))
              return(NULL)
            }

            xx <- x[center_cell_idx, intersect(cells, colnames(x))] # CD8-Trm may not exist
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
          }) %>%
          purrr::discard(is.null)
        
        res_nn_count_df <-
          res_nn_count_df_ls %>%
          enframe(name = "project") %>%
          unnest(value)  %>%
          dplyr::filter(radius != "200-Inf")
        
        res_nn_count_df
      })
  })

## plot -----------
### 1 plot  - group ------------
p_proximity_count_caf_center_ls_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(nn_df_ls, caf_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
          ggplot(aes(x = radius, y = count, fill = cell) ) +
          geom_col() +
          ggsci::scale_fill_igv() +
          theme_bw() +
          labs(
            title = name,
            x = NULL,
            y = NULL
          )
      })
  })


p_proximity_freq_caf_center_ls_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(nn_df_ls, caf_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
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
  })


# export 
p_proximity_count_caf_center_ls_ls %>%
  iwalk(function(p_ls, caf_cell){
    pdf(
      file = file.path(
        out_dir_260531, 
        str_glue("2.1 proximity analysis barplot of {caf_cell} - count - all slide.pdf")
      ),
      width  = 6,
      height = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })


p_proximity_freq_caf_center_ls_ls %>%
  iwalk(function(p_ls, caf_cell){
    pdf(
      file = file.path(
        out_dir_260531,
        str_glue("2.1 proximity analysis barplot of {caf_cell} - freq - all slide.pdf")),
      width  = 6,
      height = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })

### 2 plot  - group - lym_meta ------------
# lym_meta_ls <- list(
#   group_primary      = group_primary,
#   group_lym_meta     = group_lym_meta,
#   group_non_lym_meta = group_non_lym_meta
# )

p_proximity_count_caf_center_ls_subgrp_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(nn_df_ls, caf_cell){
    lym_meta_ls %>%
      imap(function(proj_ls, grp_name){
        ###
        nn_df_ls %>%
          imap(function(nn_df, name){
            nn_df %>%
              
              dplyr::filter(project %in% proj_ls) %>% # filter by sub group
              
              ggplot(aes(x = radius, y = count, fill = cell) ) +
              geom_col() +
              ggsci::scale_fill_igv() +
              theme_bw() +
              labs(
                title = str_glue("{grp_name} - {name}"),
                x = NULL,
                y = NULL
              )
          })
        ###
      })
  })


p_proximity_freq_caf_center_ls_subgrp_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(nn_df_ls, caf_cell){
    lym_meta_ls %>%
      imap(function(proj_ls, grp_name){
        ###
        nn_df_ls %>%
          imap(function(nn_df, name){
            nn_df %>%
              
              dplyr::filter(project %in% proj_ls) %>% # filter by sub group
              
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
                title = str_glue("{grp_name} - {name}"),
                x = NULL,
                y = NULL
              )
          })
        ###
      })
  })


# export 
p_proximity_count_caf_center_ls_subgrp_ls %>%
  iwalk(function(p_ls_ls, caf_cell){
    p_ls_ls %>%
      iwalk(function(p_ls, sub_grp){
        pdf(
          file = file.path(
            out_dir_260531, 
            str_glue("2.2 proximity analysis barplot of {caf_cell} - count - {sub_grp}.pdf")
          ),
          width  = 6,
          height = 8
        )
        p_ls %>% walk(print)
        dev.off()
      })
  })


p_proximity_freq_caf_center_ls_subgrp_ls %>%
  iwalk(function(p_ls_ls, caf_cell){
    p_ls_ls %>%
      iwalk(function(p_ls, sub_grp){
        pdf(
          file = file.path(
            out_dir_260531,
            str_glue("2.2 proximity analysis barplot of {caf_cell} - freq - {sub_grp}.pdf")),
          width  = 6,
          height = 8
        )
        p_ls %>% walk(print)
        dev.off()
      })
  })

### 3 plot  - group - cold and hot ------------
# lym_cold_hot_ls <- list(
#   group_non_meta      = group_non_meta,
#   group_lym_meta_cold = group_lym_meta_cold,
#   group_lym_meta_hot  = group_lym_meta_hot
# )

p_proximity_count_caf_center_ls_subgrp2_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(nn_df_ls, caf_cell){
    lym_cold_hot_ls %>%
      imap(function(proj_ls, grp_name){
        ###
        nn_df_ls %>%
          imap(function(nn_df, name){
            nn_df %>%
              
              dplyr::filter(project %in% proj_ls) %>% # filter by sub group
              
              ggplot(aes(x = radius, y = count, fill = cell) ) +
              geom_col() +
              ggsci::scale_fill_igv() +
              theme_bw() +
              labs(
                title = str_glue("{grp_name} - {name}"),
                x = NULL,
                y = NULL
              )
          })
        ###
      })
  })


p_proximity_freq_caf_center_ls_subgrp2_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(nn_df_ls, caf_cell){
    lym_cold_hot_ls %>%
      imap(function(proj_ls, grp_name){
        ###
        nn_df_ls %>%
          imap(function(nn_df, name){
            nn_df %>%
              
              dplyr::filter(project %in% proj_ls) %>% # filter by sub group
              
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
                title = str_glue("{grp_name} - {name}"),
                x = NULL,
                y = NULL
              )
          })
        ###
      })
  })


# export 
p_proximity_count_caf_center_ls_subgrp2_ls %>%
  iwalk(function(p_ls_ls, caf_cell){
    p_ls_ls %>%
      iwalk(function(p_ls, sub_grp){
        pdf(
          file = file.path(
            out_dir_260531, 
            str_glue("2.3 proximity analysis barplot of {caf_cell} - count - {sub_grp}.pdf")
          ),
          width  = 6,
          height = 8
        )
        p_ls %>% walk(print)
        dev.off()
      })
  })


p_proximity_freq_caf_center_ls_subgrp2_ls %>%
  iwalk(function(p_ls_ls, caf_cell){
    p_ls_ls %>%
      iwalk(function(p_ls, sub_grp){
        pdf(
          file = file.path(
            out_dir_260531,
            str_glue("2.3 proximity analysis barplot of {caf_cell} - freq - {sub_grp}.pdf")),
          width  = 6,
          height = 8
        )
        p_ls %>% walk(print)
        dev.off()
      })
  })

### 4 plot - project --------------
p_proximity_count_caf_center_ls_proj_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(nn_df_ls, caf_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
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
  })

p_proximity_proj_freq_caf_center_ls_proj_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(nn_df_ls, caf_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
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
  })


# export
p_proximity_count_caf_center_ls_proj_ls %>%
  iwalk(function(p_ls, caf_cell){
    pdf(
      file = file.path(
        out_dir_260531,
        str_glue("2.4 proximity analysis barplot of {caf_cell} - project - count.pdf")),
      width  = 12,
      height = 24
    )
    p_ls %>% walk(print)
    dev.off()
  })


p_proximity_proj_freq_caf_center_ls_proj_ls %>%
  iwalk(function(p_ls, caf_cell){
    pdf(
      file = file.path(
        out_dir_260531,
        str_glue("2.4 proximity analysis barplot of {caf_cell} - project - freq.pdf")),
      width  = 12,
      height = 24
    )
    p_ls %>% walk(print)
    dev.off()
  })

# write out cell count of others cells of epithelial
res_pa_count_tidy_caf_center_all <-
  res_nn_count_df_ls_caf_center_ls %>%
  map(function(nn_df_ls){
    nn_df_ls %>%
      imap(function(nn_df, major){
        nn_df %>%
          mutate(
            Others = str_c(.env[['major']], cell, radius, sep = "__")
          ) %>%
          arrange(project) %>%
          dplyr::select(project, Others, count) %>%
          pivot_wider(names_from = "project", values_from = "count")
      }) %>%
      bind_rows()
  })



# export 
openxlsx::write.xlsx( 
  res_pa_count_tidy_caf_center_all, 
  file = file.path(out_dir_260531, str_glue("2 proximity analysis of each cluster of CAF - all cell counts table.xlsx")), 
  asTable = TRUE 
) 

# 3. proximity analysis of tumor cellls  ---------------
# 肿瘤细胞与各CAF亚群, 转移淋巴结
all_data_meta_pos_subset2 <- 
  all_data_meta_pos_subset %>%
  dplyr::filter(Slide_ID %in% group_lym_meta) %>% #  only lym_meta
  dplyr::filter(Classification1 %in% c('Epithelial', caf_cells)) %>%
  mutate(Classification1 = fct_drop(Classification1)) 

res_pa_200um_2 <-
  ProximutyAnalysis_NN_multi(
    all_data_meta_pos_subset2,
    query = all_data_meta_pos_subset2,
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

# saveRDS(res_pa_200um_1, file = file.path(out_dir_260531, str_glue("res_pa_200um_2.rds")))


cell_match_var_ls_epi <- 
  all_data_meta_pos_subset2 %>% 
  with(split(Classification1, Slide_ID)) %>%
  map(as.character)

# epi
center_cell_id_epi <-
  all_data_meta_pos_subset2$cell_id[
    all_data_meta_pos_subset2$Classification1 %in% "Epithelial"
  ]

## compute mean distance  ----------------
res_nn_df_ls_epi <-
  res_pa_200um_2 %>%
  # .[1] %>%
  imap(function(nn, grp){ # browser()
    print(grp)
    print(Sys.time())
    res_ls <- list()
    
    cell_match_var <- cell_match_var_ls_epi[[grp]]
    
    nn.idx  <- nn$nn.idx[, -1] + 1
    nn.dist <- nn$nn.dists[, -1]
    total.cell <- nn$total_cell_id
    
    cell_match_var <- factor( c("Others", cell_match_var) ) #### ###
    
    # only compute cencter cells
    center_idx <- nn$query_cell_id %in% center_cell_id_epi
    if(sum(center_idx) < 1) return(NULL)
    
    for(i in 1:nrow(nn$nn.idx[center_idx, ])){ 
      if(i %% 5000 == 0) print(i)
      
      # head 10 cells, mean dist
      res <- tapply(nn.dist[i, ], cell_match_var[nn.idx[i, ]], function(x){mean(head(x, 10))})
      res_ls[[i]] <- res
    }
    
    print(Sys.time())
    
    res_df <- as.data.frame(bind_rows(res_ls))
    rownames(res_df) <- nn$query_cell_id[center_idx]
    res_df
    # res_ls
    
  })


# cell counts  
res_nn_count_df_ls_epi_center_ls <-
  c("Epithelial") %>%
  setNames(., .) %>%
  map(function(cell){
    
    cell_annotaion_ls['caf_cells'] %>%
      imap(function(cells, name){
        
        res_nn_count_df_ls <-
          res_nn_df_ls_epi %>%
          imap(function(x, proj){ # browser()
            cell_match_var <- cell_match_var_ls_epi[[proj]]
            center_cell_idx <- cell_match_var == cell
            if(!any(center_cell_idx)){
              print(str_glue("center {cell} in {proj} have no cells and is droped."))
              return(NULL)
            }
            
            xx <- x[center_cell_idx, cells]
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
          }) %>%
          purrr::discard(is.null)
        
        res_nn_count_df <-
          res_nn_count_df_ls %>%
          enframe(name = "project") %>%
          unnest(value)  %>%
          dplyr::filter(radius != "200-Inf")
        
        res_nn_count_df
      })
  })


## plot -----------
### 1 plot  - group ------------
p_proximity_count_epi_center_ls_ls <-
  res_nn_count_df_ls_epi_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
          ggplot(aes(x = radius, y = count, fill = cell) ) +
          geom_col() +
          ggsci::scale_fill_igv() +
          theme_bw() +
          labs(
            title = name,
            x = NULL,
            y = NULL
          )
      })
  })


p_proximity_freq_epi_center_ls_ls <-
  res_nn_count_df_ls_epi_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
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
  })


# export 
p_proximity_count_epi_center_ls_ls %>%
  iwalk(function(p_ls, epi_cell){
    pdf(
      file = file.path(
        out_dir_260531, 
        str_glue("3.1 proximity analysis barplot of {epi_cell} - count - only lym meta slide.pdf")
      ),
      width  = 6,
      height = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })


p_proximity_freq_epi_center_ls_ls %>%
  iwalk(function(p_ls, epi_cell){
    pdf(
      file = file.path(
        out_dir_260531,
        str_glue("3.1 proximity analysis barplot of {epi_cell} - freq - only lym meta slide.pdf")),
      width  = 6,
      height = 8
    )
    p_ls %>% walk(print)
    dev.off()
  })

### 2 plot - project --------------
p_proximity_count_epi_center_ls_proj_ls <-
  res_nn_count_df_ls_epi_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
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
  })

p_proximity_proj_freq_epi_center_ls_proj_ls <-
  res_nn_count_df_ls_epi_center_ls %>%
  imap(function(nn_df_ls, cd8_cell){
    nn_df_ls %>%
      imap(function(nn_df, name){
        nn_df %>%
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
  })


# export
p_proximity_count_epi_center_ls_proj_ls %>%
  iwalk(function(p_ls, epi_cell){
    pdf(
      file = file.path(
        out_dir_260531,
        str_glue("3.2 proximity analysis barplot of {epi_cell} - project - count.pdf")),
      width  = 12,
      height = 24
    )
    p_ls %>% walk(print)
    dev.off()
  })


p_proximity_proj_freq_epi_center_ls_proj_ls %>%
  iwalk(function(p_ls, epi_cell){
    pdf(
      file = file.path(
        out_dir_260531,
        str_glue("3.2 proximity analysis barplot of {epi_cell} - project - freq.pdf")),
      width  = 12,
      height = 24
    )
    p_ls %>% walk(print)
    dev.off()
  })

# write out cell count of others cells of epithelial
res_pa_count_tidy_epi_center_all <-
  res_nn_count_df_ls_epi_center_ls %>%
  map(function(nn_df_ls){
    nn_df_ls %>%
      imap(function(nn_df, major){
        nn_df %>%
          mutate(
            Others = str_c(.env[['major']], cell, radius, sep = "__")
          ) %>%
          arrange(project) %>%
          dplyr::select(project, Others, count) %>%
          pivot_wider(names_from = "project", values_from = "count")
      }) %>%
      bind_rows()
  })



# export 
openxlsx::write.xlsx( 
  res_pa_count_tidy_epi_center_all, 
  file = file.path(out_dir_260531, str_glue("3 proximity analysis of each cluster of epithelial - all cell counts table.xlsx")), 
  asTable = TRUE 
) 

# 4. CD8 props  ---------------
## group -----------
plot_dat_cd8t_grp <-
  CD8T[[]] %>%
  dplyr::filter(!is.na(subtype2)) %>%
  group_by(project, group, subtype2) %>%
  summarise(count = n()) %>%
  group_by(project, group) %>%
  mutate(freq = count/sum(count))

p_bar_cd8t_group <- 
  plot_dat_cd8t_grp %>%
  ggplot(aes(x = group, y = freq, color = group)) +
  geom_boxplot(width = 0.75) +
  geom_jitter(width = 0.15) +
  ggpubr::geom_pwc(tip.length = 0) +
  facet_wrap(~ subtype2, scales = "free", nrow = 1) +
  theme_classic() +
  theme(strip.background = element_rect(fill = NA, colour = NA)) +
  ggsci::scale_color_d3() +
  labs(
    x = NULL,
    y = "% of all CD8T cells"
  )

# export 
ggsave(
  p_bar_cd8t_group,
  filename = file.path(out_dir_260531, "4.1 barplot of each cluster of CD8T cells - group.pdf"),
  height = 4,
  width = 16
)

## group - hot cold -------------
plot_dat_cd8t_hotcold <-
  CD8T[[]] %>%
  mutate(group2 = case_when(
    project %in% group_non_meta ~ "non_meta",
    project %in% group_lym_meta_cold ~ "cold",
    project %in% group_lym_meta_hot  ~ "hot",
    TRUE ~ 'other'
  )) %>%
  dplyr::filter(group2 != "other") %>% 
  dplyr::filter(!is.na(subtype2)) %>%
  group_by(project, group2, subtype2) %>%
  summarise(count = n()) %>%
  group_by(project, group2) %>%
  mutate(freq = count/sum(count))

p_bar_cd8t_hotcold <- 
  plot_dat_cd8t_hotcold %>%
  ggplot(aes(x = group2, y = freq, color = group2)) +
  geom_boxplot(width = 0.75) +
  geom_jitter(width = 0.15) +
  ggpubr::geom_pwc(tip.length = 0) +
  facet_wrap(~ subtype2, scales = "free", nrow = 1) +
  theme_classic() +
  theme(strip.background = element_rect(fill = NA, colour = NA)) +
  ggsci::scale_color_d3() +
  labs(
    x = NULL,
    y = "% of all CD8T cells"
  )

# export 
ggsave(
  p_bar_cd8t_hotcold,
  filename = file.path(out_dir_260531, "4.2 barplot of each cluster of CD8T cells - cold and hot.pdf"),
  height = 4,
  width = 16
)


# 5. supp ----------

# lym_meta_ls <- list(
#   group_primary      = group_primary,
#   group_lym_meta     = group_lym_meta,
#   group_non_lym_meta = group_non_lym_meta
# )
# 
# lym_cold_hot_ls <- list(
#   group_non_meta      = group_non_meta,
#   group_lym_meta_cold = group_lym_meta_cold,
#   group_lym_meta_hot  = group_lym_meta_hot
# )

list(
  lym_meta = lym_meta_ls,
  lym_cold_hot = lym_cold_hot_ls
) %>%
  map(enframe, "group", "project") %>%
  map(unnest, "project") %>% 
  append(
    list(group = distinct(all_data[[c('group', 'project')]]) )
  ) %>%
  map(as.data.frame) %>%
  
  openxlsx::write.xlsx(
    file = file.path(out_dir_260531, 'supp. sample name in each group.xlsx'),
    asTable = TRUE
  )
