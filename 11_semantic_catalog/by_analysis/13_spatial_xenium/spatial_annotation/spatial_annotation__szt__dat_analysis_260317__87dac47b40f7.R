library(Seurat)
library(tidyverse)
source("spatial analysis.R")

# copy from dat_analysis_251223.R
# 更新细胞注释为最新的结果

# 0. ----------------
out_dir <- "dat_analysis_cycle3/result_20260317"
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

group_non_meta <- c(
  group_primary,
  group_non_lym_meta
)

stopifnot( all(removed_sample %in% all_data$project) )
stopifnot( all(group_primary %in% all_data$project) )
stopifnot( all(group_lym_meta %in% all_data$project) )
stopifnot( all(group_non_lym_meta %in% all_data$project) )

# all_data_subset <- all_data[, all_data$project %in% c(group_lym_meta, group_non_lym_meta)]

cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")

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

# 20260317 change CD8Tex to CD8Tpex
all_data_meta_pos_subset$Classification1 <- 
  all_data_meta_pos_subset$Classification1 %>%
  fct_recode("CD8-Tpex" = "CD8-Tex")


# 1. proximity analysis - all slide -----------

## 1.1 input data ------------
CD8_T_cells <- c(
  'CD8-Tnaive', 
  'CD8-Tem', 
  'CD8-Tpex',
  'CD8-Trm',
  'CD8-T'
) 

caf_cells <- c(
  "apCAF",
  "iCAF",
  "mCAF",
  "rCAF",
  "vCAF"
)

cell_annotaion_ls <-
  list(
    CAF = caf_cells
  )

# RUN all analysis ...
CD8_T_cells %>%
  setdiff("CD8-Trm") %>% 
  setNames(., .) %>%
  walk(function(core_cell){
    # core_cell <- 'CD8-Tpex'
    
    if(core_cell == "CD8-T"){
      all_data_meta_pos_subset$Classification1 <-
        all_data_meta_pos_subset$Classification1 %>%
        as.character() %>%
        fct_recode(!!!structure(
          CD8_T_cells,
          names = rep("CD8-T", length = length(CD8_T_cells))
        )) 
    }
    
    dir.create(file.path(out_dir, core_cell))
    
    all_data_meta_pos_subset3 <- 
      all_data_meta_pos_subset %>%
      dplyr::filter(Slide_ID %in% c(group_primary, group_lym_meta, group_non_lym_meta)) %>% #  ###
      mutate(slide_group = case_when(
        Slide_ID %in% group_primary ~ "Primary",
        Slide_ID %in% group_lym_meta ~ "Lym_meta",
        Slide_ID %in% group_non_lym_meta ~ "Lym_non_meta",
        TRUE ~ "Others"
      )) %>% 
      mutate(slide_group2 = case_when(
        Slide_ID %in% group_non_meta ~ "non-metastasis",
        Slide_ID %in% group_lym_meta ~ "metastasis ",
        TRUE ~ "Others"
      )) %>%
      
      dplyr::filter(Classification1 %in% c(core_cell, caf_cells)) %>%
      mutate(Classification1 = fct_drop(Classification1))
    
    all_data_meta_pos_subset3_core <- all_data_meta_pos_subset3[all_data_meta_pos_subset3$Classification1 == core_cell, ]
    
    ## 1.2 proximity analysis --------------- 
    res_pa_200um_2 <-
      ProximutyAnalysis_NN_multi(
        all_data_meta_pos_subset3,
        query = all_data_meta_pos_subset3_core,
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
    
    
    cell_match_var_core <- structure(
      as.character(all_data_meta_pos_subset3$Classification1),
      names = all_data_meta_pos_subset3$cell_id
    )
    
    
    res_nn_df_ls_core <-
      res_pa_200um_2 %>%
      # .[1] %>%
      imap(function(nn, grp){ # browser()
        print(grp)
        print(Sys.time())
        res_ls <- list()
        
        nn.idx  <- nn$nn.idx[, -1] + 1
        nn.dist <- nn$nn.dists[, -1]
        
        
        cell_match_var <- factor( c("Others", cell_match_var_core) )
        
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
    res_nn_count_df_ls_core <-
      res_nn_df_ls_core %>%
      map(function(x){ # browser()
        xx <- x[, !colnames(x) %in% c("Others", core_cell)]
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
    
    res_nn_count_df_core <-
      res_nn_count_df_ls_core %>%
      enframe(name = "project") %>%
      unnest(value)  %>%
      # left_join(distinct(all_data_meta_cn_subset[c("project", "group")]), by = "project") %>%
      dplyr::filter(radius != "200-Inf")
    
    ## 1.3 plot-----------
    # plot  - group
    p_proximity_count_core_ls <-
      cell_annotaion_ls %>%
      imap(function(x, name){
        res_nn_count_df_core %>%
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
    
    
    
    p_proximity_freq_core_ls <-
      cell_annotaion_ls %>%
      imap(function(x, name){
        res_nn_count_df_core %>%
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
      file = file.path(out_dir, core_cell, str_glue("1 proximity analysis barplot - count.pdf")),
      width  = 8,
      height = 6
    )
    p_proximity_count_core_ls %>% walk(print)
    dev.off()
    
    pdf(
      file = file.path(out_dir, core_cell, str_glue("1 proximity analysis barplot - freq.pdf")),
      width  = 8,
      height = 6
    )
    p_proximity_freq_core_ls %>% walk(print)
    dev.off()
    
    
    # plot - project
    p_proximity_proj_freq_core_ls <-
      cell_annotaion_ls %>%
      imap(function(x, name){
        res_nn_count_df_core %>%
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
      file = file.path(out_dir, core_cell, str_glue("1 proximity analysis barplot - project - freq.pdf")),
      width  = 12,
      height = 24
    )
    p_proximity_proj_freq_core_ls %>% walk(print)
    dev.off()
    
    
    # primary lym_meta lym_non_meta
    p_proximity_freq_core_ls_slide_grp <-
      cell_annotaion_ls %>% 
      imap(function(x, name){
        res_nn_count_df_core %>%
          
          mutate(slide_group = case_when(
            project %in% group_primary ~ "Primary",
            project %in% group_lym_meta ~ "Lym_meta",
            project %in% group_non_lym_meta ~ "Lym_non_meta",
            TRUE ~ "Others"
          )) %>%
          
          group_by(slide_group) %>%
          nest() %>%
          deframe() %>%
          imap(function(dat, name2){
            dat %>%
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
                title = str_glue("{name} - {name2}"),
                x = NULL,
                y = NULL
              )
          })
      })
    
    
    p_proximity_freq_core_ls_slide_grp2 <-
      cell_annotaion_ls %>% 
      imap(function(x, name){
        res_nn_count_df_core %>%
          
          mutate(slide_group = case_when(
            project %in% group_non_meta ~ "non-metastasis",
            project %in% group_lym_meta ~ "metastasis",
            TRUE ~ "Others"
          )) %>%
          
          group_by(slide_group) %>%
          nest() %>%
          deframe() %>%
          imap(function(dat, name2){
            dat %>%
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
                title = str_glue("{name} - {name2}"),
                x = NULL,
                y = NULL
              )
          })
      })
    
    
    # export 
    pdf(
      file = file.path(out_dir, core_cell, str_glue("proximity analysis barplot - freq - slide group.pdf")),
      width  = 8,
      height = 6
    )
    p_proximity_freq_core_ls_slide_grp %>% walk(print)
    dev.off()
    
    pdf(
      file = file.path(out_dir, core_cell, str_glue("proximity analysis barplot - freq - slide group2.pdf")),
      width  = 8,
      height = 6
    )
    p_proximity_freq_core_ls_slide_grp2 %>% walk(print)
    dev.off()
    
    
    # write out cell count of others cells of epithelial
    res_pa_count_tidy_core_all <-
      res_nn_count_df_core %>%
      mutate(
        # project = str_c(group, project, sep = "__"),
        Others  = str_c(cell, radius, sep = "__")
      ) %>%
      arrange(project) %>%
      dplyr::select(project, Others, count) %>%
      pivot_wider(names_from = "project", values_from = "count")
    
    
    # export 
    openxlsx::write.xlsx( 
      res_pa_count_tidy_core_all, 
      file = file.path(out_dir, core_cell, str_glue("proximity analysis of {core_cell} - all cell counts table.xlsx")), 
      asTable = TRUE 
    ) 
    
    
    ## density plot ----------------
    p_den_ls_core_ls <-
      cell_annotaion_ls %>%
      imap(function(x, name){
        res_nn_df_ls_core  %>% 
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
    
    
    
    p_den_total_core_ls <-
      cell_annotaion_ls %>%
      imap(function(x, name){
        dat <- res_nn_df_ls_core  %>% 
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
    
    
    p_den_total_core_ls_slide_grp <-
      cell_annotaion_ls %>%
      imap(function(x, name){
        slide_ls <-
          list(
            "Primary" = group_primary,
            "Lym_meta"= group_lym_meta,
            "Lym_non_meta" = group_non_lym_meta
          )
        
        slide_ls %>%
          imap(function(slide, name2){
            dat <- res_nn_df_ls_core[slide]  %>% 
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
              ggtitle(str_glue("{name} - {name2}"))
          })
      })
    
    
    p_den_total_core_ls_slide_grp2 <-
      cell_annotaion_ls %>%
      imap(function(x, name){
        slide_ls2 <-
          list(
            "non-metastasis" = group_non_meta,
            "metastasis"= group_lym_meta
          )
        
        slide_ls2 %>%
          imap(function(slide, name2){
            dat <- res_nn_df_ls_core[slide]  %>% 
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
              ggtitle(str_glue("{name} - {name2}"))
          })
      })
    
    # export
    pdf(
      file = file.path(out_dir, core_cell, str_glue("proximity analysis denstiy plot - all.pdf")),
      width  = 8,
      height = 6
    )
    p_den_total_core_ls %>% walk(print)
    dev.off()
    
    pdf(
      file = file.path(out_dir, core_cell, str_glue("proximity analysis denstiy plot - slide group.pdf")),
      width  = 8,
      height = 6
    )
    p_den_total_core_ls_slide_grp %>% walk(print)
    dev.off()
    
    pdf(
      file = file.path(out_dir, core_cell, str_glue("proximity analysis denstiy plot - slide group2.pdf")),
      width  = 8,
      height = 6
    )
    p_den_total_core_ls_slide_grp2 %>% walk(print)
    dev.off()
    
    
    p_den_ls_core_ls %>%
      iwalk(function(p_ls, name){
        pdf(
          file = file.path(out_dir, core_cell, str_glue("proximity analysis denstiy plot - project - {name}.pdf")),
          width  = 8,
          height = 6
        )
        p_ls %>% walk(print)
        dev.off()
        
      })
  })


# 2. correlation analysis ----------------
slide_ls <-
  list(
    all_slide = TRUE,
    # primary = group_primary,
    lym_meta = group_lym_meta,
    non_lym_meta = group_non_lym_meta,
    non_meta = group_non_meta
  )


# run
slide_ls %>%
  iwalk(function(slide, name){
    out_dir_cor <- file.path(out_dir, "correlation", name)
    dir.create(out_dir_cor, recursive = TRUE)
    
    if(name == "all_slide"){
      all_dat_used <- all_data_meta_pos_subset
    }else{
      all_dat_used <-
        all_data_meta_pos_subset %>%
        dplyr::filter(Slide_ID %in% slide)
    }
    
    
    
    ## 2.1 CD8Tpex and others - counts ---------
    all_cells <- all_dat_used$Classification1 %>% unique()
    
    
    plot_dat1_count <-
      all_dat_used %>%
      group_by(Slide_ID, Classification1) %>%
      summarise(count = n()) %>%
      pivot_wider(
        names_from = "Classification1",
        values_from = "count",
        values_fill = 0
      ) %>%
      column_to_rownames('Slide_ID')
    
    
    p_cor_cd8tpex_count_ls <-
      all_cells %>%
      setdiff("CD8-Tpex") %>% # .[1] %>%
      setNames(., .) %>%
      map(function(cell){ # browser()
        dat <- plot_dat1_count
        
        cor_test <- cor.test(dat[['CD8-Tpex']], dat[[cell]])
        cor_val <- cor_test$estimate %>% signif(2)
        cor_p <- cor_test$p.value %>% signif(2)
        cor_label <- str_glue("R = {cor_val} p = {cor_p}")
        
        p1 <- 
          dat %>%
          ggplot(aes(x = `CD8-Tpex`, y = .data[[cell]])) +
          geom_point() +
          geom_smooth(method = "lm", color = "black") +
          
          theme_classic() +
          annotate(
            'text',
            x = Inf,
            y = ifelse(cor_val > 0, 0, Inf),
            label = cor_label,
            vjust = ifelse(cor_val > 0, -2, 2),
            hjust = 1.2
          )
        p2 <- 
          p1 +
          ggtitle('log transfrom') +
          scale_x_continuous(trans =  "log1p") +
          scale_y_continuous(trans =  "log1p") 
        
        p1 + p2
      })
    
    # export
    pdf(
      file.path(out_dir_cor, '1. correlation of CD8Tpex and others cells - cell count.pdf'),
      height = 4,
      width = 9
    )
    p_cor_cd8tpex_count_ls %>% walk(print)
    dev.off()
    
    
    ## 2.2 CD8Tpex and others - freq ---------
    plot_dat2_freq <-
      all_dat_used %>%
      group_by(Slide_ID, Classification1) %>%
      summarise(count = n()) %>%
      mutate(count = count/sum(count)) %>%
      pivot_wider(
        names_from = "Classification1",
        values_from = "count",
        values_fill = 0
      ) %>%
      column_to_rownames('Slide_ID')
    
    
    p_cor_cd8tpex_freq_ls <-
      all_cells %>%
      setdiff("CD8-Tpex") %>% # .[1] %>%
      setNames(., .) %>%
      map(function(cell){ # browser()
        dat <- plot_dat2_freq
        
        cor_test <- cor.test(dat[['CD8-Tpex']], dat[[cell]])
        cor_val <- cor_test$estimate %>% signif(2)
        cor_p <- cor_test$p.value %>% signif(2)
        cor_label <- str_glue("R = {cor_val} p = {cor_p}")
        
        p1 <- 
          dat %>%
          ggplot(aes(x = `CD8-Tpex`, y = .data[[cell]])) +
          geom_point() +
          geom_smooth(method = "lm", color = "black") +
          
          theme_classic() +
          annotate(
            'text',
            x = Inf,
            y = ifelse(cor_val > 0, 0, Inf),
            label = cor_label,
            vjust = ifelse(cor_val > 0, -2, 2),
            hjust = 1.2
          )
        
        p1
      })
    
    # export
    pdf(
      file.path(out_dir_cor, '2. correlation of CD8Tpex and others cells - cell freq.pdf'),
      height = 4,
      width = 5
    )
    p_cor_cd8tpex_freq_ls %>% walk(print)
    dev.off()
    
    ## 2.3 CD8Tpex/immune cells, CAF/CAF ----------------
    immune_cells <- c(
      'Macrophage',
      'B-naive',
      'CD8-Tem',
      'CD4-Tmemory',
      'CD4-Tnaive',
      'CD8-Tnaive',
      'CD4T-Tcm',
      'Tfh',
      'Treg',
      'proliferating-T',
      'Mast',
      'pDCs',
      'NK',
      'cDC1',
      'B-memory',
      'CD8-Tpex',
      'B-GC',
      'Monocyte',
      'CD8-Trm'
    )
    
    caf_cells
    
    plot_dat3_freq_immune <-
      all_dat_used %>%
      dplyr::filter(Classification1 %in% immune_cells) %>%
      group_by(Slide_ID, Classification1) %>%
      summarise(count = n()) %>%
      mutate(count = count/sum(count)) %>%
      pivot_wider(
        names_from = "Classification1",
        values_from = "count",
        values_fill = 0
      ) %>%
      column_to_rownames('Slide_ID')
    
    plot_dat3_freq_caf <-
      all_dat_used %>%
      dplyr::filter(Classification1 %in% caf_cells) %>%
      group_by(Slide_ID, Classification1) %>%
      summarise(count = n()) %>%
      mutate(count = count/sum(count)) %>%
      pivot_wider(
        names_from = "Classification1",
        values_from = "count",
        values_fill = 0
      ) %>%
      column_to_rownames('Slide_ID')
    
    stopifnot(all(rownames(plot_dat3_freq_immune) == rownames(plot_dat3_freq_caf)))
    
    plot_dat3_freq <- cbind(
      plot_dat3_freq_caf,
      plot_dat3_freq_immune['CD8-Tpex']
    )
    
    p_cor_cd8tpex_caf_ls <-
      caf_cells %>%
      setdiff("CD8-Tpex") %>% # .[1] %>%
      setNames(., .) %>%
      map(function(cell){ # browser()
        dat <- plot_dat3_freq
        
        cor_test <- cor.test(dat[['CD8-Tpex']], dat[[cell]])
        cor_val <- cor_test$estimate %>% signif(2)
        cor_p <- cor_test$p.value %>% signif(2)
        cor_label <- str_glue("R = {cor_val} p = {cor_p}")
        
        p1 <- 
          dat %>%
          ggplot(aes(x = `CD8-Tpex`, y = .data[[cell]])) +
          geom_point() +
          geom_smooth(method = "lm", color = "black") +
          
          theme_classic() +
          annotate(
            'text',
            x = Inf,
            y = ifelse(cor_val > 0, 0, Inf),
            label = cor_label,
            vjust = ifelse(cor_val > 0, -2, 2),
            hjust = 1.2
          )
        
        p1
      })
    
    # export
    pdf(
      file.path(out_dir_cor, '3. correlation of CD8Tpex in immunes cells and caf cells - cell freq.pdf'),
      height = 4,
      width = 5
    )
    p_cor_cd8tpex_caf_ls %>% walk(print)
    dev.off()
    
    ## 2.4 CD8Tpex/CD8T, CAF/CAF ----------------
    cd8t_cells <- c(
      'CD8-Tem',
      'CD8-Tnaive',
      'CD8-Tpex',
      'CD8-Trm'
    )
    
    
    plot_dat4_freq_immune <-
      all_dat_used %>%
      dplyr::filter(Classification1 %in% cd8t_cells) %>%
      group_by(Slide_ID, Classification1) %>%
      summarise(count = n()) %>%
      mutate(count = count/sum(count)) %>%
      pivot_wider(
        names_from = "Classification1",
        values_from = "count",
        values_fill = 0
      ) %>%
      column_to_rownames('Slide_ID')
    
    plot_dat4_freq_caf <-
      all_dat_used %>%
      dplyr::filter(Classification1 %in% caf_cells) %>%
      group_by(Slide_ID, Classification1) %>%
      summarise(count = n()) %>%
      mutate(count = count/sum(count)) %>%
      pivot_wider(
        names_from = "Classification1",
        values_from = "count",
        values_fill = 0
      ) %>%
      column_to_rownames('Slide_ID')
    
    stopifnot(all(rownames(plot_dat4_freq_immune) == rownames(plot_dat4_freq_caf)))
    
    plot_dat4_freq <- cbind(
      plot_dat4_freq_caf,
      plot_dat4_freq_immune['CD8-Tpex']
    )
    
    p_cor_cd8tpex_caf_ls2 <-
      caf_cells %>%
      setdiff("CD8-Tpex") %>% # .[1] %>%
      setNames(., .) %>%
      map(function(cell){ # browser()
        dat <- plot_dat4_freq
        
        cor_test <- cor.test(dat[['CD8-Tpex']], dat[[cell]])
        cor_val <- cor_test$estimate %>% signif(2)
        cor_p <- cor_test$p.value %>% signif(2)
        cor_label <- str_glue("R = {cor_val} p = {cor_p}")
        
        p1 <- 
          dat %>%
          ggplot(aes(x = `CD8-Tpex`, y = .data[[cell]])) +
          geom_point() +
          geom_smooth(method = "lm", color = "black") +
          
          theme_classic() +
          annotate(
            'text',
            x = Inf,
            y = ifelse(cor_val > 0, 0, Inf),
            label = cor_label,
            vjust = ifelse(cor_val > 0, -2, 2),
            hjust = 1.2
          )
        
        p1
      })
    
    # export
    pdf(
      file.path(out_dir_cor, '4. correlation of CD8Tpex in CD8T cells and caf cells - cell freq.pdf'),
      height = 4,
      width = 5
    )
    p_cor_cd8tpex_caf_ls2 %>% walk(print)
    dev.off()
  })



# 3. correlation analysis - CN ------------------
all_data_meta_cn %>%
  group_by(CN) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  iwalk(function(dat, cn_name){
    if(! cn_name %in% c( "CN4")) return()
    dir.create(file.path(out_dir, 'correlation-CN', cn_name), recursive = TRUE)
    
    # run
    slide_ls %>%
      iwalk(function(slide, name){  
        out_dir_cor <- file.path(out_dir, 'correlation-CN', cn_name, name)
        dir.create(out_dir_cor, recursive = TRUE)
        
        if(name == "all_slide"){
          all_dat_used <- all_data_meta_pos_subset %>% dplyr::filter(cell_id %in% dat$cell_id)
        }else{
          all_dat_used <-
            all_data_meta_pos_subset %>%
            dplyr::filter(Slide_ID %in% slide) %>%
            dplyr::filter(cell_id %in% dat$cell_id)
        }
        
        # if(name == "all_slide") return()
        # if(name == "lym_meta") return()
        # if(name == "non_meta") browser()
        
        ## 3.1 CD8Tpex and others - counts ---------
        all_cells <- all_dat_used$Classification1 %>% unique()
        
        
        plot_dat1_count <-
          all_dat_used %>%
          group_by(Slide_ID, Classification1) %>%
          summarise(count = n()) %>%
          pivot_wider(
            names_from = "Classification1",
            values_from = "count",
            values_fill = 0
          ) %>%
          column_to_rownames('Slide_ID')
        
        if(nrow(plot_dat1_count) < 3) return()
        if(! "CD8-Tpex" %in% colnames(plot_dat1_count)) return()
        
        p_cor_cd8tpex_count_ls <-
          all_cells %>%
          setdiff("CD8-Tpex") %>% # .[1] %>%
          .[. %in% colnames(plot_dat1_count)] %>%
          setNames(., .) %>%
          map(function(cell){ # browser()
            dat <- plot_dat1_count
            # if( any(!c("CD8-Tpex", cell) %in% colnames(dat)) ) return()
            
            cor_test <- cor.test(dat[['CD8-Tpex']], dat[[cell]])
            cor_val <- cor_test$estimate %>% signif(2)
            cor_p <- cor_test$p.value %>% signif(2)
            cor_label <- str_glue("R = {cor_val} p = {cor_p}")
            
            p1 <- 
              dat %>%
              ggplot(aes(x = `CD8-Tpex`, y = .data[[cell]])) +
              geom_point() +
              geom_smooth(method = "lm", color = "black") +
              
              theme_classic() +
              annotate(
                'text',
                x = Inf,
                y = ifelse(cor_val > 0, 0, Inf),
                label = cor_label,
                vjust = ifelse(cor_val > 0, -2, 2),
                hjust = 1.2
              )
            p2 <- 
              p1 +
              ggtitle('log transfrom') +
              scale_x_continuous(trans =  "log1p") +
              scale_y_continuous(trans =  "log1p") 
            
            p1 + p2
          })
        
        # export
        pdf(
          file.path(out_dir_cor, '1. correlation of CD8Tpex and others cells - cell count.pdf'),
          height = 4,
          width = 9
        )
        p_cor_cd8tpex_count_ls %>% walk(print)
        dev.off()
        
        
        ## 3.2 CD8Tpex and others - freq ---------
        plot_dat2_freq <-
          all_dat_used %>%
          group_by(Slide_ID, Classification1) %>%
          summarise(count = n()) %>%
          mutate(count = count/sum(count)) %>%
          pivot_wider(
            names_from = "Classification1",
            values_from = "count",
            values_fill = 0
          ) %>%
          column_to_rownames('Slide_ID')
        
        
        p_cor_cd8tpex_freq_ls <-
          all_cells %>%
          setdiff("CD8-Tpex") %>% # .[1] %>%
          .[. %in% colnames(plot_dat2_freq)] %>%
          setNames(., .) %>%
          map(function(cell){ # browser()
            dat <- plot_dat2_freq
            
            cor_test <- cor.test(dat[['CD8-Tpex']], dat[[cell]])
            cor_val <- cor_test$estimate %>% signif(2)
            cor_p <- cor_test$p.value %>% signif(2)
            cor_label <- str_glue("R = {cor_val} p = {cor_p}")
            
            p1 <- 
              dat %>%
              ggplot(aes(x = `CD8-Tpex`, y = .data[[cell]])) +
              geom_point() +
              geom_smooth(method = "lm", color = "black") +
              
              theme_classic() +
              annotate(
                'text',
                x = Inf,
                y = ifelse(cor_val > 0, 0, Inf),
                label = cor_label,
                vjust = ifelse(cor_val > 0, -2, 2),
                hjust = 1.2
              )
            
            p1
          })
        
        # export
        pdf(
          file.path(out_dir_cor, '2. correlation of CD8Tpex and others cells - cell freq.pdf'),
          height = 4,
          width = 5
        )
        p_cor_cd8tpex_freq_ls %>% walk(print)
        dev.off()
        
        ## 3.3 CD8Tpex/immune cells, CAF/CAF ----------------
        immune_cells <- c(
          'Macrophage',
          'B-naive',
          'CD8-Tem',
          'CD4-Tmemory',
          'CD4-Tnaive',
          'CD8-Tnaive',
          'CD4T-Tcm',
          'Tfh',
          'Treg',
          'proliferating-T',
          'Mast',
          'pDCs',
          'NK',
          'cDC1',
          'B-memory',
          'CD8-Tpex',
          'B-GC',
          'Monocyte',
          'CD8-Trm'
        )
        
        caf_cells
        
        plot_dat3_freq_immune <-
          all_dat_used %>%
          dplyr::filter(Classification1 %in% immune_cells) %>%
          group_by(Slide_ID, Classification1) %>%
          summarise(count = n()) %>%
          mutate(count = count/sum(count)) %>%
          pivot_wider(
            names_from = "Classification1",
            values_from = "count",
            values_fill = 0
          ) %>%
          column_to_rownames('Slide_ID')
        
        plot_dat3_freq_caf <-
          all_dat_used %>%
          dplyr::filter(Classification1 %in% caf_cells) %>%
          group_by(Slide_ID, Classification1) %>%
          summarise(count = n()) %>%
          mutate(count = count/sum(count)) %>%
          pivot_wider(
            names_from = "Classification1",
            values_from = "count",
            values_fill = 0
          ) %>%
          column_to_rownames('Slide_ID')
        

        plot_dat3_freq <-
          merge(
            plot_dat3_freq_caf,
            plot_dat3_freq_immune['CD8-Tpex'],
            all = TRUE,
            by = "row.names"
          ) %>%
          column_to_rownames("Row.names")
        plot_dat3_freq[is.na(plot_dat3_freq)] <- 0
          
        
        p_cor_cd8tpex_caf_ls <-
          caf_cells %>%
          .[. %in% colnames(plot_dat3_freq)] %>%
          setdiff("CD8-Tpex") %>% # .[1] %>%
          setNames(., .) %>%
          map(function(cell){ # browser()
            dat <- plot_dat3_freq
            if(nrow(dat) < 3) return()
            
            cor_test <- cor.test(dat[['CD8-Tpex']], dat[[cell]])
            cor_val <- cor_test$estimate %>% signif(2)
            cor_p <- cor_test$p.value %>% signif(2)
            cor_label <- str_glue("R = {cor_val} p = {cor_p}")
            
            p1 <- 
              dat %>%
              ggplot(aes(x = `CD8-Tpex`, y = .data[[cell]])) +
              geom_point() +
              geom_smooth(method = "lm", color = "black") +
              
              theme_classic() +
              annotate(
                'text',
                x = Inf,
                y = ifelse(cor_val > 0, 0, Inf),
                label = cor_label,
                vjust = ifelse(cor_val > 0, -2, 2),
                hjust = 1.2
              )
            
            p1
          })
        
        # export
        pdf(
          file.path(out_dir_cor, '3. correlation of CD8Tpex in immunes cells and caf cells - cell freq.pdf'),
          height = 4,
          width = 5
        )
        p_cor_cd8tpex_caf_ls %>% walk(print)
        dev.off()
        
        ## 3.4 CD8Tpex/CD8T, CAF/CAF ----------------
        cd8t_cells <- c(
          'CD8-Tem',
          'CD8-Tnaive',
          'CD8-Tpex',
          'CD8-Trm'
        )
        
        
        plot_dat4_freq_immune <-
          all_dat_used %>%
          dplyr::filter(Classification1 %in% cd8t_cells) %>%
          group_by(Slide_ID, Classification1) %>%
          summarise(count = n()) %>%
          mutate(count = count/sum(count)) %>%
          pivot_wider(
            names_from = "Classification1",
            values_from = "count",
            values_fill = 0
          ) %>%
          column_to_rownames('Slide_ID')
        
        plot_dat4_freq_caf <-
          all_dat_used %>%
          dplyr::filter(Classification1 %in% caf_cells) %>%
          group_by(Slide_ID, Classification1) %>%
          summarise(count = n()) %>%
          mutate(count = count/sum(count)) %>%
          pivot_wider(
            names_from = "Classification1",
            values_from = "count",
            values_fill = 0
          ) %>%
          column_to_rownames('Slide_ID')
        
        # stopifnot(all(rownames(plot_dat4_freq_immune) == rownames(plot_dat4_freq_caf)))
        
        plot_dat4_freq <-
          merge(
            plot_dat4_freq_caf,
            plot_dat4_freq_immune['CD8-Tpex'],
            all = TRUE,
            by = "row.names"
          ) %>%
          column_to_rownames("Row.names")
        plot_dat4_freq[is.na(plot_dat4_freq)] <- 0
        
        p_cor_cd8tpex_caf_ls2 <-
          caf_cells %>% 
          .[. %in% colnames(plot_dat4_freq)] %>%
          setdiff("CD8-Tpex") %>% # .[1] %>%
          setNames(., .) %>%
          map(function(cell){ # browser()
            dat <- plot_dat4_freq
            if(nrow(dat) < 3) return()
            
            cor_test <- cor.test(dat[['CD8-Tpex']], dat[[cell]])
            cor_val <- cor_test$estimate %>% signif(2)
            cor_p <- cor_test$p.value %>% signif(2)
            cor_label <- str_glue("R = {cor_val} p = {cor_p}")
            
            p1 <- 
              dat %>%
              ggplot(aes(x = `CD8-Tpex`, y = .data[[cell]])) +
              geom_point() +
              geom_smooth(method = "lm", color = "black") +
              
              theme_classic() +
              annotate(
                'text',
                x = Inf,
                y = ifelse(cor_val > 0, 0, Inf),
                label = cor_label,
                vjust = ifelse(cor_val > 0, -2, 2),
                hjust = 1.2
              )
            
            p1
          })
        
        # export
        pdf(
          file.path(out_dir_cor, '4. correlation of CD8Tpex in CD8T cells and caf cells - cell freq.pdf'),
          height = 4,
          width = 5
        )
        p_cor_cd8tpex_caf_ls2 %>% walk(print)
        dev.off()
      })
    
  })
