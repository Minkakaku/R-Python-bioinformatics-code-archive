library(tidyverse)
library(Seurat)

CD8T <- readRDS("dat_analysis_cycle3/subtype/CD8T/CD8_subtype.rds")
CD8T$sub.cluster %>% is.na %>% table

out_dir_260417 <- "dat_analysis_cycle3/result_20260428"
dir.create(out_dir_260417)

# 1 CD8T subcluster adjustment --------------

# change cluster0 and cluster1 to CD8_Tpex
cell_cluster <- c(
  '0'  = "CD8_Tpex",
  '1'  = "CD8_Tpex",
  '2'  = "CD8_Tem",
  '3_0'  = "CD8_Tnaive",
  '3_1'  = "CD8_Tnaive",
  '3_2'  = "CD8_Trm",
  '3_3'  = "CD8_Tnaive",
  '3_4'  = "CD8_Tem",
  '3_5'  = "CD8_Tnaive",
  '4'  = "CD8_Tem",
  '5'  = "CD8_Tem",
  '6'  = "CD8_Tem",
  '7'  = "CD8_Tem",
  '8'  = "CD8_Tem",
  '9'  = "CD8_Tex",
  '10'  = "CD8_Tem",
  '11'  = "CD8_Tem",
  '12'  = "CD8_Tem"
)

Idents(CD8T) <- "sub.cluster"
CD8T <- RenameIdents(CD8T, cell_cluster)
CD8T$subtype2 <- Idents(CD8T)


CD8T <- 
  CD8T %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'harmonyumap
    dims = 1:20,  # <----
    
    refdata = list(subtype_full2 = "subtype2")
  )

# saveRDS(CD8T, "dat_analysis_cycle3/subtype/CD8T/CD8_subtype2.rds")

pdf(
  file.path(out_dir_260417, "1 umap plot of CD8T annotation.pdf"),
  height = 6,
  width = 8
)
# DimPlot(
#   CD8T,
#   group.by = "subtype_full2",
#   reduction = "full.harmonyumap"
# )
DimPlot(
  CD8T,
  group.by = "subtype_full2",
  reduction = "full.harmonyumap.1"
)
dev.off()

# 2 proximity analysis --------------
all_data <- readRDS(file.path("dat_analysis_cycle3", "subtype", "all_data_cycle3_250219.rds"))


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


# add new CD8T annotation
meta_anno_new <-
  left_join(
    rownames_to_column(all_data[['annotation5']], "cell_id"),
    rownames_to_column(CD8T[['subtype_full2']], "cell_id"),
    by = "cell_id"
  ) %>%
  mutate(annotation6 = ifelse(is.na(subtype_full2), annotation5, subtype_full2)) %>%
  column_to_rownames("cell_id") 

all_data <- all_data %>% AddMetaData(meta_anno_new['annotation6'])

# saveRDS(all_data, "dat_analysis_cycle3/subtype/all_data_cycle3_260417.rds")


# table(
#   meta_anno_new$annotation5[!is.na(meta_anno_new$subtype_full2)], 
#   meta_anno_new$subtype_full2[!is.na(meta_anno_new$subtype_full2)] 
# )
# 
# meta_anno_new$annotation5[is.na(meta_anno_new$subtype_full2)] %>% unique %>% str_subset("CD8")

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


## 2.1 each cluster of CD-T cell ------------------
# all cafs cells around CD8-T cell

### 1 input data ------------
caf_cells <- c(
  "apCAF",
  "iCAF",
  "mCAF",
  "rCAF",
  "vCAF"
)

cd8t_cells <- c(
  'CD8-Tem',
  'CD8-Tpex',
  'CD8-Tex',
  'CD8-Tnaive',
  'CD8-Trm'
) 

cell_annotaion_ls <- list(
  CAF = caf_cells,
  cd8t_cells = cd8t_cells
)


all_data_meta_pos_subset1 <- 
  all_data_meta_pos_subset %>%
  dplyr::filter(Slide_ID %in% c( group_primary, group_lym_meta, group_non_lym_meta)) %>% #  ###
  # mutate(slide_group = case_when(
  #   Slide_ID %in% group_primary ~ "Primary",
  #   Slide_ID %in% group_lym_meta ~ "Metastatic ",
  #   Slide_ID %in% group_non_lym_meta ~ "Metastatic ",
  #   TRUE ~ "Others"
  # )) %>%
  
  dplyr::filter(Classification1 %in% c(caf_cells, CD8_T_cells)) %>%
  mutate(Classification1 = fct_drop(Classification1)) 

# all_data_meta_pos_subset1_cd8tpex <- all_data_meta_pos_subset1[all_data_meta_pos_subset1$Classification1 == "CD8-Tpex", ]

### 2 proximity analysis --------------- 
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


# saveRDS(res_pa_200um, file = file.path(out_dir, str_glue("mcaf_res_pa_200um.rds")))


# cell_match_var_cd8t_caf <- structure(
#   as.character(all_data_meta_pos_subset1$Classification1),
#   names = all_data_meta_pos_subset1$cell_id
# )

cell_match_var_cd8t_caf_ls <- 
  all_data_meta_pos_subset1 %>% 
  with(split(Classification1, Slide_ID)) %>%
  map(as.character)

res_nn_df_ls_cd8t_caf <-
  res_pa_200um_1 %>%
  # .[1] %>%
  imap(function(nn, grp){ # browser()
    print(grp)
    print(Sys.time())
    res_ls <- list()
    
    cell_match_var_cd8t_caf <- cell_match_var_cd8t_caf_ls[[grp]]
    
    nn.idx  <- nn$nn.idx[, -1] + 1
    nn.dist <- nn$nn.dists[, -1]
    total.cell <- nn$total_cell_id
    
    cell_match_var <- factor( c("Others", cell_match_var_cd8t_caf) ) # <------
    
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
res_nn_count_df_ls_cd8t_center_ls <-
  cd8t_cells %>%
  setNames(., .) %>%
  map(function(cell){
    
    res_nn_count_df_ls_cd8t_caf <-
      res_nn_df_ls_cd8t_caf %>%
      imap(function(x, proj){ # browser()
        cell_match_var_cd8t_caf <- cell_match_var_cd8t_caf_ls[[proj]]
        center_cell_idx <- cell_match_var_cd8t_caf == cell
        if(!any(center_cell_idx)){
          print(str_glue("cencter {cell} in {proj} have no cells and is droped."))
          return(NULL)
        }
        
        xx <- x[center_cell_idx, caf_cells]
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
    
    res_nn_count_df_cd8t_caf <-
      res_nn_count_df_ls_cd8t_caf %>%
      enframe(name = "project") %>%
      unnest(value)  %>%
      # left_join(distinct(all_data_meta_cn_subset[c("project", "group")]), by = "project") %>%
      dplyr::filter(radius != "200-Inf")
    
    res_nn_count_df_cd8t_caf
  })
  


### 3 plot -----------
# plot  - group
p_proximity_count_cd8t_center_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(x, name){
    x %>%
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


p_proximity_freq_cd8t_center_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(x, name){
    x %>%
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
  file = file.path(out_dir_260417, str_glue("2 proximity analysis barplot of each cluster of CD8-T - count - all metastatic slide.pdf")),
  width  = 6,
  height = 8
)
p_proximity_count_cd8t_center_ls %>% walk(print)
dev.off()

pdf(
  file = file.path(out_dir_260417, str_glue("2 proximity analysis barplot of each cluster of CD8-T - freq - all metastatic slide.pdf")),
  width  = 6,
  height = 8
)
p_proximity_freq_cd8t_center_ls %>% walk(print)
dev.off()



# plot - project
p_proximity_proj_count_cd8t_center_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(x, name){
    x %>%
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

p_proximity_proj_freq_cd8t_center_ls <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  imap(function(x, name){
    x %>%
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
  file = file.path(out_dir_260417, str_glue("2 proximity analysis barplot of each cluster of CD8-T - project - count.pdf")),
  width  = 12,
  height = 24
)
p_proximity_proj_count_cd8t_center_ls %>% walk(print)
dev.off()

pdf(
  file = file.path(out_dir_260417, str_glue("2 proximity analysis barplot of each cluster of CD8-T - project - freq.pdf")),
  width  = 12,
  height = 24
)
p_proximity_proj_freq_cd8t_center_ls %>% walk(print)
dev.off()


# write out cell count of others cells of epithelial
res_pa_count_tidy_cd8t_center_all <-
  res_nn_count_df_ls_cd8t_center_ls %>%
  map(function(x){
    x %>%
      mutate(
        # project = str_c(group, project, sep = "__"),
        Others  = str_c(cell, radius, sep = "__")
      ) %>%
      arrange(project) %>%
      dplyr::select(project, Others, count) %>%
      pivot_wider(names_from = "project", values_from = "count")
    
  })
  

# export 
openxlsx::write.xlsx( 
  res_pa_count_tidy_cd8t_center_all, 
  file = file.path(out_dir_260417, str_glue("2 proximity analysis of each cluster of CD8-T - all cell counts table.xlsx")), 
  asTable = TRUE 
) 




### 4 density plot ----------------
res_nn_df_ls_cd8t_center_ls <-
  cd8t_cells %>%
  setNames(., .) %>%
  map(function(cell){
    res_nn_df_ls_cd8t_caf  %>% 
      imap(function(dat, proj){
        
        cell_match_var_cd8t_caf <- cell_match_var_cd8t_caf_ls[[proj]]
        center_cell_idx <- cell_match_var_cd8t_caf == cell
        if(!any(center_cell_idx)){
          print(str_glue("cencter {cell} in {proj} have no cells and is droped."))
          return(NULL)
        }
        
        dat[center_cell_idx, caf_cells]
      }) %>%
      purrr::discard(is.null)
  })


p_den_ls_cd8t_center_ls <-
  res_nn_df_ls_cd8t_center_ls %>%
  map(function(dat_ls){
    dat_ls  %>% 
      imap(function(dat, proj){
        
        dat2 <- dat %>%
          pivot_longer(
            everything(),
            names_to = 'CAF',
            values_to = "dist"
          )
        dat2$dist[is.na(dat2$dist)] <- 1000
        
        dat2 %>%
          ggplot(aes(x = dist, fill = CAF)) +
          # geom_density(position = "fill")
          geom_density(alpha = 0.5) +
          theme_classic() +
          ggtitle(proj)
      })
  })  




p_den_total_cd8t_center_ls <-
  res_nn_df_ls_cd8t_center_ls %>%
  imap(function(dat_ls, name){
    dat <- dat_ls  %>% 
      bind_rows() %>%
      pivot_longer(
        everything(),
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


p_den_total_cd8t_center_ls_slide_grp <-
  res_nn_df_ls_cd8t_center_ls %>%
  imap(function(dat_ls, name){
    slide_ls <-
      list(
        "Primary" = group_primary,
        "Lym_meta"= group_lym_meta,
        "Lym_non_meta" = group_non_lym_meta
      )
    
    slide_ls %>%
      imap(function(slide, name2){
        if(!any(slide %in% names(dat_ls))){
          print(str_glue("cencter {name} in {name2} have no cells and is droped."))
          return(NULL)
        }
        
        dat <- dat_ls[slide] %>% 
          bind_rows() %>%
          pivot_longer(
            everything(),
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
      }) %>%
      purrr::discard(is.null)
  })


# export
pdf(
  file = file.path(out_dir_260417, str_glue("2 proximity analysis denstiy plot - all.pdf")),
  width  = 8,
  height = 6
)
p_den_total_cd8t_center_ls %>% walk(print)
dev.off()

pdf(
  file = file.path(out_dir_260417, str_glue("2 proximity analysis denstiy plot - slide group.pdf")),
  width  = 8,
  height = 6
)
p_den_total_cd8t_center_ls_slide_grp %>% walk(print)
dev.off()

p_den_ls_cd8t_center_ls %>%
  iwalk(function(p_ls, name){
    pdf(
      file = file.path(out_dir_260417, str_glue("2 proximity analysis denstiy plot - project - {name}.pdf")),
      width  = 8,
      height = 6
    )
    p_ls %>% walk(print)
    dev.off()
    
  })


## 2.2 each cluster of caf cell ------------------
# all CD8T cells around caf cells


### 2 proximity analysis --------------- 
# cell count prop 
res_nn_count_df_ls_caf_center_ls <-
  caf_cells %>%
  setNames(., .) %>%
  map(function(cell){
    
    res_nn_count_df_ls_cd8t_caf <-
      res_nn_df_ls_cd8t_caf %>%
      imap(function(x, proj){ # browser()
        cell_match_var_cd8t_caf <- cell_match_var_cd8t_caf_ls[[proj]]
        center_cell_idx <- cell_match_var_cd8t_caf == cell
        others_cell_idx <- cd8t_cells[cd8t_cells %in% colnames(x)]
        if(!any(center_cell_idx)){
          print(str_glue("cencter {cell} in {proj} have no cells and is droped."))
          return(NULL)
        }
        
        xx <- x[center_cell_idx, others_cell_idx]
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
    
    res_nn_count_df_cd8t_caf <-
      res_nn_count_df_ls_cd8t_caf %>%
      enframe(name = "project") %>%
      unnest(value)  %>%
      # left_join(distinct(all_data_meta_cn_subset[c("project", "group")]), by = "project") %>%
      dplyr::filter(radius != "200-Inf")
    
    res_nn_count_df_cd8t_caf
  })



### 3 plot -----------
# plot  - group
p_proximity_count_caf_center_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(x, name){
    x %>%
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


p_proximity_freq_caf_center_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(x, name){
    x %>%
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
  file = file.path(out_dir_260417, str_glue("3 proximity analysis barplot of each cluster of caf - count - all metastatic slide.pdf")),
  width  = 6,
  height = 8
)
p_proximity_count_caf_center_ls %>% walk(print)
dev.off()

pdf(
  file = file.path(out_dir_260417, str_glue("3 proximity analysis barplot of each cluster of caf - freq - all metastatic slide.pdf")),
  width  = 6,
  height = 8
)
p_proximity_freq_caf_center_ls %>% walk(print)
dev.off()



# plot - project
p_proximity_proj_count_caf_center_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(x, name){
    x %>%
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

p_proximity_proj_freq_caf_center_ls <-
  res_nn_count_df_ls_caf_center_ls %>%
  imap(function(x, name){
    x %>%
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
  file = file.path(out_dir_260417, str_glue("3 proximity analysis barplot of each cluster of caf - project - count.pdf")),
  width  = 12,
  height = 24
)
p_proximity_proj_count_caf_center_ls %>% walk(print)
dev.off()

pdf(
  file = file.path(out_dir_260417, str_glue("3 proximity analysis barplot of each cluster of caf - project - freq.pdf")),
  width  = 12,
  height = 24
)
p_proximity_proj_freq_caf_center_ls %>% walk(print)
dev.off()


# write out cell count of others cells of epithelial
res_pa_count_tidy_caf_center_all <-
  res_nn_count_df_ls_caf_center_ls %>%
  map(function(x){
    x %>%
      mutate(
        # project = str_c(group, project, sep = "__"),
        Others  = str_c(cell, radius, sep = "__")
      ) %>%
      arrange(project) %>%
      dplyr::select(project, Others, count) %>%
      pivot_wider(names_from = "project", values_from = "count")
    
  })


# export 
openxlsx::write.xlsx( 
  res_pa_count_tidy_caf_center_all, 
  file = file.path(out_dir_260417, str_glue("3 proximity analysis of each cluster of caf - all cell counts table.xlsx")), 
  asTable = TRUE 
) 




### 4 density plot ----------------
res_nn_df_ls_caf_center_ls <-
  caf_cells %>%
  setNames(., .) %>%
  map(function(cell){
    res_nn_df_ls_cd8t_caf  %>% 
      imap(function(dat, proj){
        
        cell_match_var_cd8t_caf <- cell_match_var_cd8t_caf_ls[[proj]]
        center_cell_idx <- cell_match_var_cd8t_caf == cell
        others_cell_idx <- cd8t_cells[cd8t_cells %in% colnames(dat)]
        if(!any(center_cell_idx)){
          print(str_glue("cencter {cell} in {proj} have no cells and is droped."))
          return(NULL)
        }
        
        dat[center_cell_idx, others_cell_idx]
      }) %>%
      purrr::discard(is.null)
  })


p_den_ls_caf_center_ls <-
  res_nn_df_ls_caf_center_ls %>%
  map(function(dat_ls){
    dat_ls  %>% 
      imap(function(dat, proj){
        
        dat2 <- dat %>%
          pivot_longer(
            everything(),
            names_to = 'CD8T',
            values_to = "dist"
          )
        dat2$dist[is.na(dat2$dist)] <- 1000
        
        dat2 %>%
          ggplot(aes(x = dist, fill = CD8T)) +
          # geom_density(position = "fill")
          geom_density(alpha = 0.5) +
          theme_classic() +
          ggtitle(proj)
      })
  })  




p_den_total_caf_center_ls <-
  res_nn_df_ls_caf_center_ls %>%
  imap(function(dat_ls, name){
    dat <- dat_ls  %>% 
      bind_rows() %>%
      pivot_longer(
        everything(),
        names_to = 'CD8T',
        values_to = "dist"
      )  
    dat$dist[is.na(dat$dist)] <- 1000
    
    dat %>%
      ggplot(aes(x = dist, fill = CD8T)) +
      # geom_density(position = "fill")
      geom_density(alpha = 0.5) +
      theme_classic() +
      ggtitle(name)
  })


p_den_total_caf_center_ls_slide_grp <-
  res_nn_df_ls_caf_center_ls %>%
  imap(function(dat_ls, name){
    slide_ls <-
      list(
        "Primary" = group_primary,
        "Lym_meta"= group_lym_meta,
        "Lym_non_meta" = group_non_lym_meta
      )
    
    slide_ls %>%
      imap(function(slide, name2){
        if(!any(slide %in% names(dat_ls))){
          print(str_glue("cencter {name} in {name2} have no cells and is droped."))
          return(NULL)
        }
        
        dat <- dat_ls[slide] %>% 
          bind_rows() %>%
          pivot_longer(
            everything(),
            names_to = 'CD8T',
            values_to = "dist"
          )  
        dat$dist[is.na(dat$dist)] <- 1000
        
        dat %>%
          ggplot(aes(x = dist, fill = CD8T)) +
          # geom_density(position = "fill")
          geom_density(alpha = 0.5) +
          theme_classic() +
          ggtitle(str_glue("{name} - {name2}"))
      }) %>%
      purrr::discard(is.null)
  })


# export
pdf(
  file = file.path(out_dir_260417, str_glue("3 proximity analysis denstiy plot - all.pdf")),
  width  = 8,
  height = 6
)
p_den_total_caf_center_ls %>% walk(print)
dev.off()

pdf(
  file = file.path(out_dir_260417, str_glue("3 proximity analysis denstiy plot - slide group.pdf")),
  width  = 8,
  height = 6
)
p_den_total_caf_center_ls_slide_grp %>% walk(print)
dev.off()

p_den_ls_caf_center_ls %>%
  iwalk(function(p_ls, name){
    pdf(
      file = file.path(out_dir_260417, str_glue("3 proximity analysis denstiy plot - project - {name}.pdf")),
      width  = 8,
      height = 6
    )
    p_ls %>% walk(print)
    dev.off()
    
  })
