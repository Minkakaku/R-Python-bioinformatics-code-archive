

library(Seurat)
library(tidyverse)
source("spatial analysis.R")

CN_out_dir <- "spatial_analysis - 20251201/CN analysis"
# dir.create("spatial_analysis - 20251201", recursive = T)
dir.create(CN_out_dir, recursive = T)

all_data_filt <- readRDS("all_data_DC_annotation_1201.rds")
Idents(all_data_filt) <- "annotation4"

cells_position_ls <- readRDS("cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()


## CN analysis ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

all_data_submeta_for_cn_analysis <- 
  all_data_filt[[c("sample", "group", "annotation4")]] %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  dplyr::rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    annotation4 = factor(annotation4), # import!
    obj_id = cell_id # import!
  )

### CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "annotation4",
  k.nn           = 10,
  k.cluster      = 10
)


saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))

### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 7,
  width  = 10
)
cowplot::plot_grid(p_cn10_heat$gtable)
dev.off()


### CN k15 ---------
res_cn_15 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "sample",
  classification = "annotation4",
  k.nn           = 10,
  k.cluster      = 15
)

saveRDS(res_cn_15, file.path(CN_out_dir, "res_cn_15.rds"))

### heatmap ----------------------
# FoldChange cluster/allCell
cn15_cluster_mean <- res_cn_15$cluster_centers
cn15_baseMean     <- res_cn_15$base_mean

# log2(FC)
cn15_cluster_FC <- log2(1 + sweep(cn15_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn15_baseMean))

# heatmap
p_cn15_heat <-
  pheatmap::pheatmap(
    cn15_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn15_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_15_heatmap.pdf"),
  height = 9,
  width  = 10
)
cowplot::plot_grid(p_cn15_heat$gtable)
dev.off()



### CN annotaton --------------------
# # CN k10
# cn_annotation <- c(
#   "CN1"  = "Macrophage",
#   "CN2"  = "T",
#   "CN3"  = "B",
#   "CN4"  = "Fibro_Stroma",
#   "CN5"  = "Myo_Stroma",
#   "CN6"  = "Pan immune",
#   "CN7"  = "NK_Endothelial_Stroma",
#   "CN8"  = "Epithelial",
#   "CN9"  = "T_Epithelial",
#   "CN10" = "T_mast_Epithelial"
# )
# 
# p_cn_heat2 <-
#   cn_cluster_FC %>%
#   pheatmap::pheatmap(
#     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#     border_color = NA,
#     treeheight_row = 0,
#     labels_row = cn_annotation,
#     # clustering_method = 'ward.D2',
#     cluster_rows = F,
#     scale = "column",
#     silent = TRUE
#   )
# 
# pdf(
#   file.path(CN_out_dir, "CN_heatmap_annotation.pdf"),
#   height = 6,
#   width  = 11
# )
# cowplot::plot_grid(p_cn_heat2$gtable)
# dev.off()


### CN10 spatial plot -----------------
# cn10_sp_out_dir <- file.path(CN_out_dir, 'CN10 spatial plot') 
cn10_sp_out_dir <- file.path(CN_out_dir, 'CN10 spatial plot') 
dir.create(cn10_sp_out_dir)

cells_boundary_ls <- readRDS("cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn10 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn10, file.path(CN_out_dir, "all_data_meta_cn10.rds"))


All_type_meta_cn10_filt <- All_type_meta_cn10 %>% 
  select(-annotation, -annotation2, -annotation3) %>% 
  rename(annotation = annotation4)
colnames(All_type_meta_cn10_filt)


WriteXLS(
  x = "All_type_meta_cn10_filt", 
  ExcelFileName = "spatial_analysis - 20251201/CN analysis/All_type_meta_cn10.xlsx", 
  SheetNames = "CN10", 
  row.names = FALSE, 
  col.names = TRUE, 
  AdjWidth = TRUE, 
  BoldHeaderRow = TRUE
)





cn10_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    # image_wh_r <- image_wh_r_ls[[sample]]
    
    anno_data   <- All_type_meta_cn10[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(sample, group_info, sep = "_") )
    p
  })


# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn10_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

### cell props ------
library(ggpubr)
All_type_cn10_count_df <-
  All_type_meta_cn10 %>% 
  group_by(sample, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn10_props <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn10_props_group <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  stat_compare_means(label = "p.format",size=2,
                     aes(group = group))+
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn10_props %>% print
p_cn10_props_group %>% print
dev.off()



### CN15 spatial plot -----------------
cn15_sp_out_dir <- file.path(CN_out_dir, 'CN15 spatial plot') 
dir.create(cn15_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

All_type_meta_cn15 <- 
  all_data_filt[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_15$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn15, file.path(CN_out_dir, "all_data_meta_cn15.rds"))




All_type_meta_cn15_filt <- All_type_meta_cn15 %>% 
  select(-annotation, -annotation2, -annotation3) %>% 
  rename(annotation = annotation4)
colnames(All_type_meta_cn15_filt)


WriteXLS(
  x = "All_type_meta_cn15_filt", 
  ExcelFileName = "spatial_analysis - 20251201/CN analysis/All_type_meta_cn15.xlsx", 
  SheetNames = "CN10", 
  row.names = FALSE, 
  col.names = TRUE, 
  AdjWidth = TRUE, 
  BoldHeaderRow = TRUE
)




cn15_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, sample){ # browser()
    # image_wh_r <- image_wh_r_ls[[sample]]
    
    anno_data   <- All_type_meta_cn15[c("sample", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$sample == sample]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn15, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(sample, group_info, sep = "_") )
    p
  })


# export image
cn15_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn15_sp_out_dir, str_glue("CN15 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

### cell props ------
All_type_cn15_count_df <-
  All_type_meta_cn15 %>% 
  group_by(sample, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn15_props <-
  All_type_cn15_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn15_props_group <-
  All_type_cn15_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  stat_compare_means(label = "p.format",size=3,
                     aes(group = group))+
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN15_boxplot_and_barplot.pdf"),
  height = 6,
  width  = 11
)
p_cn15_props %>% print
p_cn15_props_group %>% print
dev.off()










# cellphone --------
library(Seurat)
library(tidyverse)
library(reticulate)
library(igraph)


# config python
use_condaenv(condaenv = "../szt/miniforge3/envs/rpy/bin/python")
py_config()
py_available()
py_module_available("pandas")

# os <- import('os')
# os$getcwd()
# os$listdir()



# test_var <- "variable created in r"
# 
# 
# repl_python() # start a python console
# test_var = "variable created in python"
# r.test_var
# 
# quit # return to R
# 
# py$test_var

# # test python code running ...
# source_python
# reticulate::py_run_file
# reticulate::py_run_string("  
# print(1)
# print(2)
# 
# for i in range(5):
#   print(i)
#     
# ")


### export cellphonedb files   ------------
### read data  -------
all_data <- readRDS('all_data_DC_annotation4_1202.rds')

dir.create('cellphonedb_result2')
all_data_meta <- all_data[[]]
all_data <- JoinLayers(all_data)
all_data_norm <- LayerData(all_data, layer ='counts')


human_mouse <- read_tsv('human_mouse_id.txt')
human_mouse <- human_mouse %>% select("Mouse gene name","Gene name")
human_mouse_clean <- drop_na(human_mouse)

new_rownames <- human_mouse_clean$`Gene name`[match(rownames(all_data_norm), human_mouse_clean$`Mouse gene name`)]
rownames(all_data_norm)<-new_rownames
all_data_norm2 <- all_data_norm[!is.na(rownames(all_data_norm)),]



all_grp_ls <-
  all_data_meta %>%
  dplyr::filter(nCount_RNA > 100) %>%
  rownames_to_column('temp') %>%
  group_by(group) %>%
  nest() %>%
  deframe()



all_grp_ls %>%
  iwalk(function(dat, name){ #browser()
    
    dat<-column_to_rownames(dat,'temp')
    cells <- rownames(all_data_meta) %in% rownames(dat)
    
    meta <- all_data_meta[cells, ] 
    dat  <- all_data_norm[, cells]
    
    out_dir <- file.path("cellphonedb_result2", 'in', name)
    dir.create(out_dir,recursive = TRUE)
    
    # write meta files
    meta['annotation2'] %>%
      data.table::fwrite(
        file = file.path(out_dir, "meta.csv"),
        row.names = TRUE
      )
    
    # write count files
    dat %>%
      as.matrix() %>% 
      as.data.frame() %>%
      data.table::fwrite(
        file = file.path(out_dir, "counts.csv"),
        row.names = TRUE
      )
    
  })








### download cellphone db files ------
reticulate::py_run_string("
import os
from cellphonedb.utils import db_utils

# Version of the databse
cpdb_version = 'v5.0.0'

# Path where the input files to generate the database are located
cpdb_target_dir = os.path.join('cellphonedb/database', cpdb_version)
db_utils.download_database(cpdb_target_dir, cpdb_version)

")



### running cellphondb ------
cellphonedb_path <- "cellphonedb/database/v5.0.0/cellphonedb.zip"

### ################### python code start ########################

repl_python() ## start python
import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
db_path = "cellphonedb/database/v5.0.0/cellphonedb.zip" 


in_dir = os.path.join('cellphonedb_result2', 'in')
os.makedirs(in_dir, exist_ok = True)
out_dir = os.path.join('cellphonedb_result2', 'out')
os.makedirs(out_dir, exist_ok = True)


all_dirs = os.listdir(in_dir)

for group in all_dirs:

  meta_path   = os.path.join(in_dir, group, 'meta.csv') # meta文件，两列：细胞id及细胞类型
counts_path = os.path.join(in_dir, group, 'counts.csv') # counts文件， norm文件，行为基因，列为细胞id

out_path = os.path.join(out_dir, group)
os.makedirs(out_path, exist_ok = True)

cpdb_results = cpdb_statistical_analysis_method.call(
  cpdb_file_path = db_path,
  meta_file_path = meta_path,
  counts_file_path = counts_path,
  counts_data = 'hgnc_symbol',
  output_path = out_path,
  output_suffix = group,
  score_interactions = True,
  threads = 10
)

quit

### ################### python code end. ########################


### plot cellphonedb ------------------
cpdb_plot_out_dir <- "cellphonedb_result"

all_files_ls <- 
  list.files(file.path(cpdb_plot_out_dir, "out")) %>%
  set_names(., basename(.)) %>%
  map(function(x){
    all_files <- list.files(str_glue('{cpdb_plot_out_dir}/out/{x}/'), full.names = TRUE)
    all_files
  })



### parse cellphonedb输出文件：means文件   -----
network_dat_ls <-
  all_files_ls %>%
  map(function(all_files){# browser()
    network_dat <- all_files  %>% 
      str_subset('_analysis_pvalues_') %>% 
      read_tsv() %>% 
      dplyr::select(interacting_pair, where(is.numeric)) %>%
      pivot_longer(
        -interacting_pair,
        names_to = 'cell_pair',
        values_to = 'p'
      ) %>% 
      dplyr::filter(p < 0.05) %>% 
      group_by(cell_pair) %>% 
      summarise(counts = n()) %>% 
      mutate(SOURCE = str_split(cell_pair, "\\|") %>% map_chr(1),
             TARGET = str_split(cell_pair, "\\|") %>% map_chr(2)) %>% 
      dplyr::select(SOURCE, TARGET, counts) %>% 
      dplyr::filter(counts != 0) 
    
    network_dat
  })


### full dotplot ----------
net_ls <-
  network_dat_ls %>%
  map(function(network_dat){
    all_nodes <- unique(c(network_dat$SOURCE, network_dat$TARGET))
    
    if(!exists('cell_color')){
      cell_color <- ggsci::pal_igv()(length(all_nodes))
      names(cell_color) <- all_nodes
    }
    
    
    net <- graph_from_data_frame(network_dat)
    
    E(net)$width <- E(net)$counts#  %>% {(. - min(.))/diff(range(.)) + 1}
    E(net)$label <- E(net)$counts
    
    color_match_table     <- cell_color[names(cell_color) %in% all_nodes] 
    E(net)$color <- attr(E(net), 'vnames' ) %>% str_remove("\\|.*$")  %>% { color_match_table[.]}
    V(net)$color <- color_match_table[names(V(net))] 
    
    # return
    net
  })


# plot
net_ls %>%
  iwalk(function(net, name){# browser()
    karate_groups <- cluster_optimal(net)
    # coords        <- layout_in_circle(net, order = order(membership(karate_groups)))  # 设置网络布局
    coords        <- layout_in_circle(net)  # 设置网络布局
    
    # 调整节点位置的线条角度
    edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
    angles <- atan(coords[igraph::V(net),2]/coords[igraph::V(net),1])
    
    loop.angle <- ifelse(coords[igraph::V(net), 1]>0, -angles, pi-angles)
    loop.edge.idx <- which(edge.start[,2] == edge.start[,1])
    
    igraph::E(net)$loop.angle <- 0
    igraph::E(net)$loop.angle[loop.edge.idx] <- loop.angle[edge.start[loop.edge.idx, 1]]
    
    
    # export pdf
    pdf(
      file.path(cpdb_plot_out_dir, str_glue("{name} network plot.pdf")),
      height = 10,
      width  = 10
    )
    plot(
      net,
      edge.curved=0.3,
      edge.arrow.size=0.3,
      layout = coords,
      margin = rep(0.2, 4),
      main = name
    )
    dev.off()
  })

### split dotplot ------
net_ls %>%
  iwalk(function(net, name){
    all_nodes <- attr(V(net), "names")
    
    # export pdf
    pdf(
      file.path(cpdb_plot_out_dir, str_glue("{name} network plot - single node.pdf")),
      height = 10,
      width  = 10
    )
    all_nodes %>%
      walk(function(x){ # browser()
        
        coords        <- layout_in_circle(net)  # 设置网络布局
        
        # 调整节点位置的线条角度
        edge.start <- igraph::ends(net, es=igraph::E(net), names=FALSE)
        angles <- atan(coords[igraph::V(net),2]/coords[igraph::V(net),1])
        
        loop.angle <- ifelse(coords[igraph::V(net), 1]>0, -angles, pi-angles)
        loop.edge.idx <- which(edge.start[,2] == edge.start[,1])
        
        igraph::E(net)$loop.angle <- 0
        igraph::E(net)$loop.angle[loop.edge.idx] <- loop.angle[edge.start[loop.edge.idx, 1]]
        
        # 只保留当前细胞发出的线
        E(net)$color[! str_starts(attr(E(net), 'vnames' ), fixed(x))] <- NA
        E(net)$label[! str_starts(attr(E(net), 'vnames' ), fixed(x))] <- NA
        
        plot(
          net,
          edge.curved=0.3,
          edge.arrow.size=0.3,
          # label.cex = 15,
          layout = coords,
          margin = rep(0.2, 4),
          main = name
        )
        
      })
    dev.off()
  })

### heatmap --------
p_cpdb_heatmap_ls2 <-
  network_dat_ls %>%
  imap(function(network_dat, name){
    #browser()
    dat <- network_dat %>% 
      pivot_wider(
        names_from  = "TARGET",
        values_from = "counts",
        values_fill = 0
      ) %>%
      column_to_rownames("SOURCE")
    
    # pheatmap::pheatmap(
    #   dat,
    #   scale = "none",
    #   color = colorRampPalette(c("#124E8B", "#FDD8B8", "#8C0C51"))(50),
    #   method = "ward.D2",
    #   border_color = NA,
    #   treeheight_row = 0,
    #   treeheight_col = 0,
    #   # main = name,
    #   silent = F,
    #   na_col = "white"
    # )
    cluster_row <- hclust(dist(dat), method = "ward.D2")
    cluster_col <- hclust(dist(t(dat)), method = "ward.D2")
    
    dat %>%
      rownames_to_column("SOURCE") %>%
      pivot_longer(
        -SOURCE,
        names_to = "TARGET",
        values_to = "counts"
      ) %>%
      mutate(
        SOURCE = factor(SOURCE, levels = rev(cluster_row$labels[cluster_row$order])),
        TARGET = factor(TARGET, levels = cluster_col$labels[cluster_col$order]),
      ) %>%
      ggplot(aes(x = TARGET, y = SOURCE, fill = counts)) +
      geom_tile() +
      scale_fill_gradientn(colours = c("#124E8B", "#FDD8B8", "#8C0C51")) +
      theme_minimal() +
      coord_fixed() +
      ggtitle(name) +
      theme(
        # axis.title = element_text(),
        axis.ticks.length.x = unit(0, 'mm'),
        axis.text.x = element_text(
          angle = -90,
          hjust = 0, 
          vjust = 0.5,
          size = 15
        ),
        axis.text.y = element_text(size = 15),
        axis.title  = element_text(size = 20)
      )
    
  })


p_cpdb_heatmap_combn_ls2 <-
  p_cpdb_heatmap_ls2 %>%
  enframe(name = "group", value = "p") %>%
  separate(group, into = c("Group", "group"), sep = "_") %>%
  dplyr::select(-group) %>%
  group_by(Group) %>%
  nest() %>%
  deframe() %>%
  map(1) %>%
  map(function(x){ # browser()
    rng_ls <- x %>% map('data') %>% map('counts') %>% map(range)
    min_val <- rng_ls %>% map_dbl(1) %>% min()
    max_val <- rng_ls %>% map_dbl(2) %>% max()
    
    x %>%
      map(~ .x + scale_fill_gradientn(colours = c("#124E8B", "#FDD8B8", "#8C0C51"), limits = c(min_val, max_val))) %>%
      patchwork::wrap_plots()
  })

# export heatmap plot
pdf(
  file.path(cpdb_plot_out_dir, str_glue("heatmap plot.pdf")),
  height = 12,
  width  = 12
)
p_cpdb_heatmap_combn_ls2 %>%
  walk(function(p){
    cowplot::plot_grid(p) %>% print()
  })
dev.off()


### dotplot ----------

# 解析 mean文件和pvalues文件
inter_df_mean_ls <- 
  all_files_ls %>%
  map(function(all_files){
    all_files  %>%
      str_subset('_analysis_means_') %>%
      read_tsv() %>%
      dplyr::select(interacting_pair, where(is.numeric)) %>%
      pivot_longer(
        -interacting_pair,
        names_to = 'cell_pair',
        values_to = 'mean'
      ) %>%
      dplyr::filter(!is.na(mean)) %>%
      dplyr::filter(!cell_pair == 'rank') %>% 
      mutate(help_col = str_c(interacting_pair, cell_pair)) %>% 
      dplyr::select(-c(1, 2)) 
  })



inter_df_p_ls <- 
  all_files_ls %>%
  map(function(all_files){ # browser()
    all_files  %>%
      str_subset('pvalues') %>%
      read_tsv() %>%
      dplyr::select(interacting_pair, where(is.numeric)) %>%
      pivot_longer(
        -interacting_pair,
        names_to = 'cell_pair',
        values_to = 'p'
      ) %>%
      # mutate(p = p.adjust(p)) %>%
      dplyr::filter(p < 0.05) %>% # <------ p值: 没有筛选的作用 
      mutate(help_col = str_c(interacting_pair, cell_pair))
  })

inter_df_inter_score_ls <- 
  all_files_ls %>%
  map(function(all_files){ # browser()
    all_files  %>%
      str_subset('_interaction_scores_') %>%
      read_tsv() %>%
      dplyr::select(interacting_pair, where(is.numeric)) %>%
      pivot_longer(
        -interacting_pair,
        names_to = 'cell_pair',
        values_to = 'score'
      ) %>%
      # mutate(p = p.adjust(p)) %>%
      # dplyr::filter(p < 0.001) %>% # <------ p值: 没有筛选的作用 
      mutate(help_col = str_c(interacting_pair, cell_pair))
  })


dot_dat_ls <-
  inter_df_p_ls %>%
  imap(function(inter_df_p, name){ # browser()
    inter_df_mean <- inter_df_mean_ls[[name]]
    inter_df_inter_score <- inter_df_inter_score_ls[[name]]
    
    inter_df_p %>%
      inner_join(inter_df_mean, by = 'help_col') %>%
      inner_join(inter_df_inter_score[c('score', 'help_col')], by = 'help_col') %>%
      mutate(p = str_replace(p,"^0$","1e-12")) %>%
      group_by(interacting_pair) %>%
      # dplyr::filter(p <= 0.05) %>%
      arrange(p) %>%
      mutate('-log10(pvalue)' = -log10(as.numeric(p))) %>%
      mutate(mean = log2(mean)) 
  })


# dotplot
p_dot_ls <-
  dot_dat_ls %>%
  imap(function(dot_dat, name){
    dot_dat  %>% 
      as.data.frame() %>% 
      mutate(
        cell_pair2 = 
          cell_pair %>%
          str_split("\\|") %>% 
          map(~ str_remove_all(.x, "\\s.*$")) %>% 
          map_chr(str_c, collapse = "|")
      ) %>% 
      
      # slice_sample(prop = 0.1) %>%
      
      ggplot(aes(
        x = cell_pair2,
        y = interacting_pair,
        color = `-log10(pvalue)`
      )) +
      geom_point(shape = 22, aes(fill = mean, size = score), stroke = 2) +
      theme_classic() +
      labs(
        x = NULL,
        y = NULL,
        title = name,
        fill = "Log2 mean(Molecule 1, Molecule 2)",
        # size = "-log10(p value)"
        size = "Interaction scores"
      ) +
      scale_fill_gradientn(colors = c("#fce5d7", "#f5ae8b",  "#d00f19" ,  "#6e020f")) +
      scale_size(range = c(0.1, 5)) +
      scale_color_gradient(low = "#386cb044", high = "#386cb0") + 
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ))
  })

# export dotplot
pdf(
  file.path(cpdb_plot_out_dir, str_glue("dotplot.pdf")),
  height = 30,
  width  = 50
)
p_dot_ls %>% walk(print)
dev.off()


### dotplot - split ------
# ligand
p_dot_ls2 <-
  dot_dat_ls %>%
  enframe(name = "group") %>%
  unnest(value) %>%
  mutate(group_ligand = str_remove(cell_pair, "\\|.*$")) %>%
  
  group_by(group_ligand) %>% 
  nest() %>%
  deframe() %>%
  
  imap(function(dat, name){ # browser()
    dat %>% 
      as.data.frame() %>% 
      # mutate(
      #   cell_pair2 = 
      #     cell_pair %>%
      #     str_split("\\|") %>% 
      #     map(~ str_remove_all(.x, "\\s.*$")) %>% 
      #     map_chr(str_c, collapse = "|")
      # ) %>% 
      
      # slice_sample(prop = 0.1) %>%
      
      ggplot(aes(
        x = cell_pair,
        y = interacting_pair,
        color = `-log10(pvalue)`
      )) +
      geom_point(shape = 22, aes(fill = mean, size = score), stroke = 2) +
      ggpubr::theme_pubclean() +
      labs(
        x = NULL,
        y = NULL,
        title = str_glue("{name} - ligand"),
        fill = "Log2 mean(Molecule 1, Molecule 2)",
        # size = "-log10(p value)"
        size = "Interaction scores"
      ) +
      facet_wrap(~ group, ncol = 3) +
      scale_fill_gradientn(colors = c("#fce5d7", "#f5ae8b",  "#d00f19" ,  "#6e020f")) +
      scale_size(range = c(0.1, 5)) +
      scale_color_gradient(low = "#386cb044", high = "#386cb0") + 
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ),
        legend.position = "right"
      )
  })

# receptor
p_dot_ls3 <-
  dot_dat_ls %>%
  enframe(name = "group") %>%
  unnest(value) %>%
  mutate(group_receptor = str_remove(cell_pair, "^.*\\|")) %>%
  
  group_by(group_receptor) %>% 
  nest() %>%
  deframe() %>%
  
  imap(function(dat, name){ # browser()
    dat %>% 
      as.data.frame() %>% 
      # mutate(
      #   cell_pair2 = 
      #     cell_pair %>%
      #     str_split("\\|") %>% 
      #     map(~ str_remove_all(.x, "\\s.*$")) %>% 
      #     map_chr(str_c, collapse = "|")
      # ) %>% 
      
      # slice_sample(prop = 0.1) %>%
      
      ggplot(aes(
        x = cell_pair,
        y = interacting_pair,
        color = `-log10(pvalue)`
      )) +
      geom_point(shape = 22, aes(fill = mean, size = score), stroke = 2) +
      ggpubr::theme_pubclean() +
      labs(
        x = NULL,
        y = NULL,
        title = str_glue("{name} - receptor"),
        fill = "Log2 mean(Molecule 1, Molecule 2)",
        # size = "-log10(p value)"
        size = "Interaction scores"
      ) +
      facet_wrap(~ group, ncol = 3) +
      scale_fill_gradientn(colors = c("#fce5d7", "#f5ae8b",  "#d00f19" ,  "#6e020f")) +
      scale_size(range = c(0.1, 5)) +
      scale_color_gradient(low = "#386cb044", high = "#386cb0") + 
      theme(
        axis.text.x = element_text(
          angle = -90,
          hjust = 0,
          vjust = 0.5
        ),
        legend.position = "right"
      )
  })


# export dotplot
pdf(
  file.path(cpdb_plot_out_dir, str_glue("dotplot2.pdf")),
  height = 20,
  width  = 15
)
p_dot_ls2 %>% walk(print)
dev.off()

## 
pdf(
  file.path(cpdb_plot_out_dir, str_glue("dotplot3.pdf")),
  height = 20,
  width  = 15
)
p_dot_ls3 %>% walk(print)
dev.off()



