library(tidyverse)
library(Seurat)


all_data <- readRDS("R_OUTPUT/all_data_brian.rds")
all_data <- all_data %>% JoinLayers()
DefaultAssay(all_data) <- "RNA"

out_dir_v5 <- "spatial_domain_segmentation"
dir.create(out_dir_v5)

cells_boundary_ls <- readRDS('R_OUTPUT/data_for_analysis/cells_boundary_ls_brian.rds')

# path_ls <-
#   readxl::excel_sheets("Supplementary Table 2_41586_2024_7806_MOESM4_ESM.xlsx") %>%
#   setNames(., str_trim(str_extract(., "\\s[A-Za-z]*\\s"))) %>%
#   setNames(., str_replace(names(.), "Ischemic", "IZ")) %>%
#   map(function(x){
#     dat <- readxl::read_excel("Supplementary Table 2_41586_2024_7806_MOESM4_ESM.xlsx", sheet = x)
#     # if('avg_log2FC' %in% colnames(dat)){
#     #   dat <-
#     #     dat %>%
#     #     dplyr::filter(avg_log2FC > .25, p_val_adj < 0.05, pct.1 > pct.2)
#     # }
#     dat$gene
#   }) %>%
#   map(str_to_title)
# 
# path_ls <- path_ls[c('BZ', 'IZ')]
# 
# 
# # 0 findmarker of singlecell data and pathway score
# all_markers_sc <-
#   FindAllMarkers(
#     JoinLayers(all_data_sc2),
#     group.by = 'annotation',
#     return.thresh = 1.1
#   )
# 
# de_gene_cm <- 
#   all_markers_sc %>%
#   dplyr::filter(cluster == 'CM') %>%
#   dplyr::filter(
#     avg_log2FC > 2,
#     p_val_adj < 0.05,
#     pct.1 - pct.2 > 0.3
#   ) %>%
#   dplyr::pull('gene')
# 
# # signature of papers
# path_ls2 <-
#   readxl::read_excel("44161_2022_160_MOESM2_ESM.xlsx", sheet = 1) %>%
#   group_by(cluster) %>%
#   nest() %>%
#   deframe() %>%
#   map('gene') %>%
#   setNames(., fct_recode(names(.), RZ = 'Myh6', BZ1 = 'Ankrd1', BZ2 = 'Xirp2'))
# 
# # path_ls[['CM']] <- de_gene_cm
# path_ls <- c(path_ls, path_ls2)

PVN <- c('Ntrk1', 'Kcnj5', 'Gpx3', 'Thbs2', 'Slc41a3', 'Insrr', 'Fam65b', 'Fam102b', 'Cntnap3', 
         'Drc1', 'Cacng3', '4932418E24Rik', 'Prokr2', 'Satb2', 'Nenf', 'Zdbf2', 'Fam20a', 'Kctd6', 
         'Baalc', 'Fhod3', 'Tmem176a', 'Pcdh18', 'Gpr151', 'Rwdd2a', 'Sdc2', 'Maoa', 'Celf6', 
         'Tspan6', 'Tmem176a', 'Nts', 'Rftn1', 'Cda', 'Nek7', 'Pdzrn3', 'Rab3b', 'D8Ertd82e', 'Doc2a', 
         'Gdpd2', 'Baiap3', 'Galr1', 'Sema3f', 'Sulf1', 'LOC381076', 'Pdxk', 'Gsto1', 'Dlk1', 'Maob', 
         'N4bp2', 'Fnbp1l', 'Hpcal4')
SCN <- c('Lhx1', 'Vipr2', 'Rasl11b', 'Rgs16', 'Rorb', 'Prokr2', 'Avp', 'Rora', 'Nr1d1', 'Zim1', 'Flt3', 
         'Zic1', 'Adcy3', 'Btg1', 'Per2', 'Ntpcr', 'Spon1', 'Rftn1', 'Dlk1', 'Neurod2', 'Trp53i11', 
         'Igfbp5', 'Id4', 'Thsd7b', 'Plch2', 'Slc2a13', 'Adck4', 'Pdcd4', 'Msi2', 'Epha6', 'Nkd1', 
         'Myt1', 'Zfhx3', 'Vgf', 'Galnt6', 'Celf3', 'Slc5a3', 'Blcap', 'B630019K06Rik', 'Ppp1r17', 
         'Plk2', 'Rps6ka2', 'Cdh13', 'Iqsec3', 'Pck2')
path_ls <- list(PVN, SCN)
names(path_ls) <- c('PVN', 'SCN')


all_data_score_zone <- 
  all_data %>%
  AddModuleScore(path_ls)

colnames(all_data_score_zone@meta.data)[
  str_detect(colnames(all_data_score_zone@meta.data), 'Cluster')
] <- 
  names(path_ls)


# 0 score spatial plot ----------
all_data_score_zone_meta <- all_data_score_zone@meta.data %>% rownames_to_column('cell_id')
target_gene_ls <- c('PVN', 'SCN')
names(target_gene_ls) <- target_gene_ls

marker_p_ls <- 
  target_gene_ls %>%
  map(function(x) {
    cells_boundary_ls %>%
      imap(function(cells_boundary, sample) { # browser()
        dat <-
          cells_boundary %>%
          right_join(all_data_score_zone_meta, by = "cell_id")
        
        p <- dat %>%
          arrange(.data[[x]]) %>% 
          as.data.frame() %>%
          ggplot(aes(x = vertex_x, y = -vertex_y)) + 
          geom_polygon(aes(fill = .data[[x]], group = cell_id))  +
          # scattermore::geom_scattermore(aes(alpha = .data[[x]]), pointsize = 0.5, show.legend = TRUE) +
          guides(alpha = "none") +
          theme_void() +
          coord_equal() +
          ggtitle(sample)+ 
          scale_fill_gradientn(
            colours = c("blue", "white", "red"),
            # limits = c(0, ceiling(max(dat_expr[-1]))),  # 统一颜色范围
            limits = c(0, 0.1), 
            name = "Expression"
          ) + 
          guides(
            fill = guide_colorbar(
              title = str_glue(x, " Expression"),
              barwidth = 0.5,
              barheight = 10
            )
          )
        p
        
      })
  })


marker_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path(out_dir_v5,  
                  str_glue("spatial dim plot score - {grp}.pdf")), 
        height = 10, 
        width = 10)
    print(p)
    dev.off()
  }) 






# 1 nn analysis ----------
cell_coords_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_position_ls_brian.rds")
cell_coords_ls2 <- cell_coords_ls %>% map(~  .[.$cell_id %in% colnames(all_data), ])

## find nn point -----------
k_nn <- 25
res_nn_ls <- 
  cell_coords_ls2 %>%
  map(function(dat){ # browser()
    cn_nn <- 
      RANN::nn2(
        dat[c('x_centroid', 'y_centroid')], 
        k = k_nn
      )
    
    list(
      nn = cn_nn,
      cell_id = dat$cell_id
    ) 
  })




## combined nn expression -----------
dat_zone_sum_ls <- 
  res_nn_ls %>%
  map(function(res_nn){ # browser()
    # dat
    print('subset data')
    dat_used <- FetchData(
      all_data_score_zone,
      vars = names(path_ls),
      cells = res_nn$cell_id
    )
    dat_used <- dat_used[res_nn$cell_id, ]
    
    stopifnot(all(rownames(dat_used) == res_nn$cell_id))
    
    
    # idx
    nn_idx <- res_nn$nn$nn.idx %>% t %>% c()
    group_idx <- rep(res_nn$cell_id, each = k_nn)
    
    # sum
    print('sum of nn cells')
    dat_merged <- dat_used %>% apply(2, function(x) tapply(x[nn_idx], group_idx, sum) )
    
    print('done.')
    
    # return
    dat_merged
  })


# # Find Cluster
# all_data_path_ls <-
#   dat_zone_sum_ls %>%
#   map(function(dat){
#     sc <- 
#       CreateSeuratObject(
#         CreateAssayObject(
#           # counts = t(dat),
#           data = t(dat)
#         ),
#         meta.data = all_data@meta.data
#       )
#     sc[['RNA']]$scale.data <- sc[['RNA']]$data
#     sc[['RNA']] <- as(sc[['RNA']], "Assay5")
#     sc
#   })
#   
# ll <- all_data_path_ls$A02994A6 %>% RunPCA(features = rownames(.)) %>% RunUMAP(dims = 1:3)


## kmeans cluster -------------
res_kmeans_used_ls <- 
  dat_zone_sum_ls %>%
  map(function(dat){
    set.seed(123)
    kmeans(dat, centers = 3) 
  })


p_heat_ls <-
  res_kmeans_used_ls %>%
  map('centers') %>%
  map(function(x){
    pheatmap::pheatmap(
      t(x),
      scale = "row",
      cluster_cols = F,
      # cluster_rows = F,
      border_color = NA,
      silent = T
    )
  })


pdf(
  file = file.path(out_dir_v5, "1 kmeans cluster PVN and SCN - heatmap.pdf"),
  height = 5,
  width  = 5
)
p_heat_ls %>% map('gtable') %>% map(cowplot::plot_grid)
dev.off()

# df_kmeans_cluster <-
#   res_kmeans_used_ls %>%
#   map('cluster') %>%
#   map(enframe, 'cell_id', 'cluster') %>%
#   imap(function(dat, name){
#     # cluster annotation by heatmap
#     if(name %in% c('PSO_F2')){
#       dat <- dat %>%
#         mutate(
#           cluster = str_c('cluster', cluster),
#           cluster = fct_recode(
#             cluster, 
#             "PVN" = 'cluster1', 
#             "SCN" = 'cluster2',
#             "others" = 'cluster3'
#           )
#         )  
#     }else{
#       dat <- dat %>%
#         mutate(
#           cluster = str_c('cluster', cluster),
#           cluster = fct_recode(
#             cluster, 
#             "PVN" = 'cluster1', 
#             "SCN" = 'cluster2',
#             "others" = 'cluster3'
#           )
#         )  
#     }
#   }) %>%
#   bind_rows() %>%
#   column_to_rownames("cell_id")  

# 读取并预处理映射表
cluster_mapping <- read.csv("spatial_PVN_SCN.csv", stringsAsFactors = FALSE)

df_kmeans_cluster <-
  res_kmeans_used_ls %>%
  map('cluster') %>%
  map(enframe, 'cell_id', 'cluster') %>%
  imap(function(dat, sample_name){
    
    mapping_row <- cluster_mapping[cluster_mapping$sample == sample_name, ]
    
    if (nrow(mapping_row) == 0) {
      # 样本不在映射表中
      return(dat %>% mutate(cluster = paste0("cluster", cluster)))
    }
    
    # 检查是否需要映射
    if (mapping_row$cluster1 == "cluster1") {
      # 不需要映射，保持原样
      return(dat %>% mutate(cluster = paste0("cluster", cluster)))
    }
    
    # 需要映射的情况
    mapping_vec <- c(mapping_row$cluster1, mapping_row$cluster2, mapping_row$cluster3)
    names(mapping_vec) <- paste0("cluster", 1:3)
    
    dat %>%
      mutate(
        cluster = paste0("cluster", cluster),
        cluster = mapping_vec[cluster]
      )
  }) %>%
  bind_rows() %>%
  column_to_rownames("cell_id")

# # 准备映射表（长格式）
# mapping_long <- cluster_mapping %>%
#   pivot_longer(
#     cols = c(cluster1, cluster2, cluster3),
#     names_to = "cluster_num",
#     values_to = "new_cluster"
#   ) %>%
#   mutate(
#     cluster_num = str_remove(cluster_num, "cluster"),
#     cluster_id = paste0("cluster", cluster_num)
#   ) %>%
#   select(sample, cluster_id, new_cluster)
# 
# # 使用left_join合并
# df_kmeans_cluster <-
#   res_kmeans_used_ls %>%
#   map('cluster') %>%
#   map(enframe, 'cell_id', 'cluster') %>%
#   imap(function(dat, sample_name){
#     dat %>%
#       mutate(
#         sample = sample_name,
#         cluster_id = paste0("cluster", cluster)
#       )
#   }) %>%
#   bind_rows() %>%
#   # 合并映射表
#   left_join(mapping_long, by = c("sample", "cluster_id")) %>%
#   # 如果没有映射，设为"others"
#   mutate(
#     cluster = ifelse(is.na(new_cluster), "others", new_cluster),
#     cluster = factor(cluster, levels = c("PVN", "SCN", "others"))
#   ) %>%
#   select(cell_id, cluster) %>%
#   column_to_rownames("cell_id")


p_sp_iz_used_ls <-
  cell_coords_ls2 %>%
  imap(function(coord, name){
    coord %>%
      left_join(
        rownames_to_column(df_kmeans_cluster, "cell_id"),
        by = "cell_id"
      ) %>%
      
      slice_sample(n = 50000) %>%
      
      ggplot(aes(x = x_centroid, y = -y_centroid, color = as.character(cluster))) +
      geom_point(size = 0.3) +
      coord_equal() +
      theme_minimal() +
      # scale_color_manual(values = c(
      #   '1' = "red",
      #   '2' = 'blue'
      # )) +
      guides(color = guide_legend(override.aes = list(size = 2))) +
      labs(
        color = 'cluster',
        title = name
      )
  })


pdf(
  file = file.path(out_dir_v5, "1 kmeans cluster PVN and SCN - spatial plot.pdf"),
  height = 20,
  width  = 20
)
p_sp_iz_used_ls %>% walk(print)
dev.off()


## vlnplot of PVN and SCN ------------
all_data_score_zone <-
  all_data_score_zone %>%
  AddMetaData(df_kmeans_cluster)

p_vln_zone_ls <-
  path_ls %>% 
  names() %>%
  map(function(x){
    all_data_score_zone %>%
      VlnPlot(
        features = x,
        group.by = "project",
        split.by = "cluster",
        pt.size = 0
      )  +
      geom_boxplot(width = 0.15, outliers = F, position = position_dodge(width = 0.9)) +
      ggtitle(x)
    
  })



pdf(
  file = file.path(out_dir_v5, "1 kmeans cluster PVN and SCN - vln plot.pdf"),
  height = 6,
  width  = 12
)
p_vln_zone_ls %>% walk(print)
dev.off()


df_kmeans_cluster %>% dim()
saveRDS(df_kmeans_cluster, 'df_kmeans_cluster.rds')




# export
# saveRDS(all_data_score_zone, 'all_data_score_zone.rds')


# # test not used  --------------
# # combined data and cluster
# res_kmeans_used <- 
#   dat_bz_iz_sum_ls2$A02994A6 %>%
#   kmeans(centers = 2) 
# 
# 
# p_sp_iz_used <-
#   cell_coords_ls2$A02994A6 %>%
#   left_join(
#    enframe(res_kmeans_used$cluster, "cell_id", "cluster"),
#    by = "cell_id"
#   ) %>%
#   
#   slice_sample(n = 50000) %>%
#   
#   ggplot(aes(x = x, y = -y, color = as.character(cluster))) +
#   geom_point() +
#   theme_minimal() +
#   scale_color_manual(values = c(
#     '1' = "red",
#     '2' = 'blue'
#   ))
# 
# 
# ggsave(
#   p_sp_iz_used,
#   filename = file.path(out_dir_v5, "kmeans cluster IZ and BZ using marker gene - used.pdf"),
#   height = 20,
#   width  = 20
# )
#   
# 
# # not used combined data
# res_kmeans_not_used <- 
#   FetchData(
#     all_data,
#     vars = path_ls[1:2] %>% unlist,
#     cells = colnames(all_data)[all_data$project == 'A02994A6']
#   ) %>%
#   kmeans(centers = 2)
# 
# 
# 
# p_sp_iz_not_used <-
#   cell_coords_ls2$A02994A6 %>%
#   left_join(
#     enframe(res_kmeans_not_used$cluster, "cell_id", "cluster"),
#     by = "cell_id"
#   ) %>%
#   
#   slice_sample(n = 50000) %>%
#   
#   ggplot(aes(x = x, y = -y, color = as.character(cluster))) +
#   geom_point() +
#   # scale_color_continuous(trans = "reverse") +
#   theme_minimal() +
#   scale_color_manual(values = c(
#     '2' = "red",
#     '1' = 'blue'
#   ))
# 
# ggsave(
#   p_sp_iz_not_used,
#   filename = file.path(out_dir_v5, "kmeans cluster IZ and BZ using marker gene - not used combined.pdf"),
#   height = 20,
#   width  = 20
# )
