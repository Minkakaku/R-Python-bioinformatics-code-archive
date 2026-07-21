library(Seurat)
library(tidyverse)

source("spatial analysis.R")

out_dir <- "dat_analysis_cycle3"
cn_out_dir <- "dat_analysis_cycle3/cn_analysis_250219" 
dir.create(cn_out_dir)

# all_data <- readRDS(file.path(out_dir, "all_data_cycle3.rds"))
# Idents(all_data) <- "annotation2"

all_data <- readRDS(file.path(out_dir, "subtype", "all_data_cycle3_250219.rds"))
Idents(all_data) <- "annotation5"

cells_position_ls <- readRDS("data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()



# implement CN ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id

all_data_meta_for_cn_analysis <- 
  all_data[[c("project", "group", "group2", "annotation5")]] %>% # 2025.2.19 used annotation5
  dplyr::rename(annotation2 = annotation5) %>%
  
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    annotation2 = factor(annotation2), # import!
    obj_id = cell_id # import!
  )

res_cn <- CN_compute_CN(
  all_data_meta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "annotation2",
  k.nn           = 50,
  k.cluster      = 10
)


# saveRDS(res_cn, file.path(cn_out_dir, "res_cn.rds"))

# heatmap ----------------------
# FoldChange cluster/allCell
cn_cluster_mean <- res_cn$cluster_centers
cn_baseMean     <- res_cn$base_mean

# log2(FC)
cn_cluster_FC <- log2(1 + sweep(cn_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn_baseMean))

# heatmap
p_cn_heat <-
  pheatmap::pheatmap(
    cn_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(cn_out_dir, "CN_heatmap_CN.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn_heat$gtable)
dev.off()

# CN annotaton --------------------
# CN k10
# cn_annotation <- c(
#   "CN1"  = "Plasma",
#   "CN2"  = "B boundary",
#   "CN3"  = "B core",
#   "CN4"  = "Tumor core",
#   "CN5"  = "Stroma_Fibroblast",
#   "CN6"  = "CD8T + pan-immune",
#   "CN7"  = "Epi_immune",
#   "CN8"  = "CD4T + pan-immune",
#   "CN9"  = "Tumor boundary",
#   "CN10" = "Stroma_Endothelial"
# )

# cn_annotation <- c(
#   "CN1"  = "Plasma",
#   "CN2"  = "B boundary",
#   "CN3"  = "B core",
#   "CN4"  = "Tumor core",
#   "CN5"  = "Endo + immnune",
#   "CN6"  = "CD8T + pan-immune",
#   "CN7"  = "Immune_epi0",
#   "CN8"  = "CD4T + pan-immune",
#   "CN9"  = "Tumor boundary",
#   "CN10" = "Stroma"
# )

cn_annotation <- c(
  "CN1"  = "Plasma",
  "CN2"  = "pan-immune1", 
  "CN3"  = "B core",
  "CN4"  = "Tumor core",
  "CN5"  = "B boundary",
  "CN6"  = "pan-immune2",
  "CN7"  = "CD4T_APC",
  "CN8"  = "CD4T + pan-immune",
  "CN9"  = "Tumor boundary",
  "CN10" = "Stroma"
)

p_cn_heat2 <-
  cn_cluster_FC %>%
  pheatmap::pheatmap(
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = cn_annotation,
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(cn_out_dir, "CN_heatmap_annotation.pdf"),
  height = 6,
  width  = 11
)
cowplot::plot_grid(p_cn_heat2$gtable)
dev.off()


# CN spatial plot -----------------
cn_sp_out_dir <- file.path(cn_out_dir, 'cn spatial plot') 
dir.create(cn_sp_out_dir)

cells_boundary_ls <- readRDS("data_for_analysis/cells_boundary_ls.rds")
image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

all_data_meta <- 
  all_data[[]] %>%
  rownames_to_column('cell_id') %>%
  left_join(rownames_to_column(res_cn$cn_cluster, "cell_id"), by = "cell_id") %>%
  mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))

saveRDS(all_data_meta, file.path(cn_out_dir, str_glue("meta_data_cn.rds")))

cn_sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN2, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() +
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


# export image
cn_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn_sp_out_dir, str_glue("CN spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

# cellphonedb -------------
## export cellphonedb files   -----
###  read data  -------
all_data_meta_cpdb <- all_data[[]] 
all_data_norm_cpdb <- all_data %>% JoinLayers( ) %>% LayerData(assay = "RNA", layer = "counts")

cpdb_dir <- file.path(cn_out_dir, 'cellphonedb_result')
dir.create(cpdb_dir)

c(
  'A','B','C'
) %>% 
  walk(function(grp){# browser()
    
    cells <- which(all_data_meta_cpdb$group == grp)
    # if(length(cells) < 1){
    #   next
    # }
    # set.seed(1234)
    # if(length(cells) > 30000){
    #   cells <- sample(cells, size = 30000)
    # }
    
    
    print("subset")
    meta <- all_data_meta_cpdb[cells, ]
    dat  <- all_data_norm_cpdb[, cells]
    
    out <- file.path(cpdb_dir, 'in', grp)
    dir.create(out, recursive = TRUE)
    
    # write meta files
    meta['annotation5'] %>%
      data.table::fwrite(
        file = file.path(out, "meta.csv"),
        row.names = TRUE
      )
    
    # write count files
    dat %>%
      as.matrix() %>%
      as.data.frame() %>%
      data.table::fwrite(
        file = file.path(out, "counts.csv"),
        row.names = TRUE
      )
    
    # print('parquet')
    # # write count parquet
    # dat %>%
    #   as.matrix() %>% 
    #   as.data.frame() %>%
    #   arrow::write_parquet(file.path(out, "counts.parquet"))
  })


## running cellphondb ------
cellphonedb_path <- "cellphonedb/database/v5.0.0/cellphonedb.zip"

### ################### python code start ########################
library(reticulate)
use_python("../szt/miniforge3/envs/rpy/bin/python")
py_module_available("cellphonedb")

repl_python() ## start python
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import os

db_path = "cellphonedb/database/v5.0.0/cellphonedb.zip" 


in_dir = os.path.join('dat_analysis_cycle3/result_20250220', 'cellphonedb_result', 'in')
os.makedirs(in_dir, exist_ok = True)
out_dir = os.path.join('dat_analysis_cycle3/result_20250220', 'cellphonedb_result', 'out')
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



# plot cellphonedb ------------------
library(igraph)

cpdb_dir <- file.path(cn_out_dir, 'cellphonedb_result')
dir.create(cpdb_dir)

all_files_ls <- 
  list.files(file.path(cpdb_dir, "out")) %>%
  set_names(., basename(.)) %>%
  map(function(x){
    all_files <- list.files(str_glue('{cpdb_dir}/out/{x}/'), full.names = TRUE)
    all_files
  })



## parse cellphonedb输出文件：means文件   -----
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


## full dotplot ----------
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
      file.path(cpdb_dir, str_glue("{name} network plot.pdf")),
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

## split dotplot ------
net_ls %>%
  iwalk(function(net, name){
    all_nodes <- attr(V(net), "names")
    
    # export pdf
    pdf(
      file.path(cpdb_dir, str_glue("{name} network plot - single node.pdf")),
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

## heatmap --------
network_dat_B <- network_dat_ls$B %>% 
  pivot_wider(
    names_from  = "TARGET",
    values_from = "counts",
    values_fill = 0
  ) %>%
  column_to_rownames("SOURCE")

p_cpdb_heatmap_ls2 <-
  network_dat_ls %>%
  imap(function(network_dat, name){
    # browser()
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
    # cluster_row <- hclust(dist(dat), method = "ward.D2")
    # cluster_col <- hclust(dist(t(dat)), method = "ward.D2")
    cluster_row <- hclust(dist(network_dat_B), method = "ward.D2")
    cluster_col <- hclust(dist(t(network_dat_B)), method = "ward.D2")
    
    levels_row <- rev(cluster_row$labels[cluster_row$order])
    levels_row <- c(levels_row, setdiff(rownames(dat), levels_row))
    
    levels_col <- cluster_col$labels[cluster_col$order]
    levels_col <- c(levels_col, setdiff(colnames(dat), levels_col))
    
    dat %>%
      rownames_to_column("SOURCE") %>%
      pivot_longer(
        -SOURCE,
        names_to = "TARGET",
        values_to = "counts"
      ) %>%
      mutate(
        SOURCE = factor(SOURCE, levels = levels_row),
        TARGET = factor(TARGET, levels = levels_col)
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

rng_ls <- p_cpdb_heatmap_ls2 %>% map('data') %>% map('counts') %>% unlist %>% range

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
    min_val <- rng_ls[1] # %>% map_dbl(1) %>% min()
    max_val <- rng_ls[2] # %>% map_dbl(2) %>% max()
    
    x %>%
      map(~ .x + scale_fill_gradientn(colours = c("#124E8B", "#FDD8B8", "#8C0C51"), limits = c(min_val, max_val))) %>%
      patchwork::wrap_plots()
  })


# export heatmap plot
pdf(
  file.path(cpdb_dir, str_glue("heatmap plot - 2.pdf")),
  height = 12,
  width  = 12
)
p_cpdb_heatmap_combn_ls2 %>%
  walk(function(p){
    cowplot::plot_grid(p) %>% print()
  })
dev.off()


## dotplot ----------

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
  file.path(cpdb_dir, str_glue("dotplot.pdf")),
  height = 6,
  width  = 50 * 2
)
p_dot_ls %>% walk(print)
dev.off()


# dotplot - split ------
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
  file.path(cpdb_dir, str_glue("dotplot2.pdf")),
  height = 8,
  width  = 25
)
p_dot_ls2 %>% walk(print)
dev.off()

## 
pdf(
  file.path(cpdb_dir, str_glue("dotplot3.pdf")),
  height = 8,
  width  = 25
)
p_dot_ls3 %>% walk(print)
dev.off()


# dotplot - user defined ----------------
# 对我们关注的以下细胞的受配体分析，类似如下输出
# 1. Epithelial：T、B、CAF、macrophage各亚群
# 2. CD8Tem、Treg、CD8Tex：DC、B、macrophage各亚群、NK
# 3. apCAF、iCAF、mCAF；其他所有细胞类型（图可能很长，我们后续根据结果再调整）
# 4. Macrophage各亚群：其他所有细胞类型（图可能很长，我们后续根据结果再调整）

df_cell_major_to_sub <- c(
  fibro = 'apCAF', 
  fibro = 'iCAF',
  fibro = 'rCAF', 
  fibro = 'vCAF',
  fibro = 'mCAF', 
  
  b =     'ASC', 
  b =     'B-GC', 
  b =     'B-memory', 
  b =     'B-naive',
  
  t =     'CD4-Tmemory', 
  t =     'CD4-Tnaive', 
  t =     'CD4T-Tcm', 
  t =     'CD8-Tem', 
  t =     'CD8-Tex', 
  t =     'CD8-Tnaive', 
  t =     'CD8-Trm', 
  t =     'Tfh', 
  t =     'Treg', 
  t =     'proliferating-T', 
  
  mye =   'cDC1', 
  mye =   'Macrophage', 
  mye =   'Mast', 
  mye =   'Monocyte', 
  mye =   'pDCs', 

  
  'NK', 
  'Fat cells', 
  'Epi0_Mast cell', 
  
  'Endothelial', 
  'Smooth muscle cells', 
  
  'Epithelial', 
  'Fibro_like Epithelial', 
  'Immune_like Epithelial'
)


## 1. Epithelial：T、B、CAF、macrophage各亚群 ----
p_dot_ls_userdefined_1 <-
  dot_dat_ls %>%
  enframe(name = "group") %>%
  unnest(value) %>%
  mutate(group_ligand = str_remove(cell_pair, "\\|.*$")) %>%
  mutate(group_receptor = str_remove(cell_pair, "^.*\\|")) %>% {
    dat <- .
    cell1 <- c('Epithelial')
    cell2 <- c(df_cell_major_to_sub[names(df_cell_major_to_sub) %in% c('t', 'b', 'fibro')], 'Macrophage')
    sub_dat1 <-
      dplyr::filter(
        dat,
        group_ligand %in% cell1,
        group_receptor %in% cell2
      ) %>%
      mutate(group_ligand = factor(group_ligand, levels = cell1)) %>%
      arrange(group_ligand) %>%
      as.data.frame() %>%
      mutate(cell_pair = factor(cell_pair,  levels = unique(cell_pair))) 
    sub_dat2 <-
      dplyr::filter(
        dat,
        group_ligand %in% cell2,
        group_receptor %in% cell1
      )%>%
      mutate(group_receptor = factor(group_receptor, levels = cell1)) %>%
      arrange(group_receptor) %>%
      as.data.frame() %>%
      mutate(cell_pair = factor(cell_pair,  levels = unique(cell_pair))) 
    
    # return
    list(
      ligand   = sub_dat1,
      receptor = sub_dat2
    )
  } %>% 
  imap(function(dat, name){ # browser()
    dat %>% 
      as.data.frame() %>% 
      
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
        title = str_glue("{name}"),
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
          hjust = ifelse(name == "ligand", 0, 1),
          vjust = 0.5
        ),
        legend.position = "right"
      )
  })

# export
pdf(
  file.path(cpdb_dir, str_glue("dotplot - user-defined1.pdf")),
  height = 8,
  width  = 20
)
p_dot_ls_userdefined_1 %>% walk(print)
dev.off()


## 2. CD8Tem、Treg、CD8Tex：DC、B、macrophage各亚群、NK ----
p_dot_ls_userdefined_2 <-
  dot_dat_ls %>%
  enframe(name = "group") %>%
  unnest(value) %>%
  mutate(group_ligand = str_remove(cell_pair, "\\|.*$")) %>%
  mutate(group_receptor = str_remove(cell_pair, "^.*\\|")) %>% {
    dat <- .
    cell1 <- c('CD8-Tem', 'Treg', 'CD8-Tex')
    cell2 <- c(df_cell_major_to_sub[names(df_cell_major_to_sub) %in% c('b')], 'cDC1', 'pDCs', 'macrophage', 'NK')
    sub_dat1 <-
      dplyr::filter(
        dat,
        group_ligand %in% cell1,
        group_receptor %in% cell2
      ) %>%
      mutate(group_ligand = factor(group_ligand, levels = cell1)) %>%
      arrange(group_ligand) %>%
      as.data.frame() %>%
      mutate(cell_pair = factor(cell_pair,  levels = unique(cell_pair))) 
    sub_dat2 <-
      dplyr::filter(
        dat,
        group_ligand %in% cell2,
        group_receptor %in% cell1
      )%>%
      mutate(group_receptor = factor(group_receptor, levels = cell1)) %>%
      arrange(group_receptor) %>%
      as.data.frame() %>%
      mutate(cell_pair = factor(cell_pair,  levels = unique(cell_pair))) 
    # return
    list(
      ligand   = sub_dat1,
      receptor = sub_dat2
    )
  } %>% 
  imap(function(dat, name){ # browser()
    dat %>% 
      as.data.frame() %>% 
      
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
        title = str_glue("{name}"),
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
          hjust = ifelse(name == "ligand", 0, 1),
          vjust = 0.5
        ),
        legend.position = "right"
      )
  })

# export
pdf(
  file.path(cpdb_dir, str_glue("dotplot - user-defined2.pdf")),
  height = 8,
  width  = 20
)
p_dot_ls_userdefined_2 %>% walk(print)
dev.off()


## 3. apCAF、iCAF、mCAF；其他所有细胞类型 ----
p_dot_ls_userdefined_3 <-
  dot_dat_ls %>%
  enframe(name = "group") %>%
  unnest(value) %>%
  mutate(group_ligand = str_remove(cell_pair, "\\|.*$")) %>%
  mutate(group_receptor = str_remove(cell_pair, "^.*\\|")) %>% {
    dat <- .
    cell1 <- c('apCAF', 'iCAF', 'mCAF')
    cell2 <- df_cell_major_to_sub
    sub_dat1 <-
      dplyr::filter(
        dat,
        group_ligand %in% cell1,
        group_receptor %in% cell2
      ) %>%
      mutate(group_ligand = factor(group_ligand, levels = cell1)) %>%
      arrange(group_ligand) %>%
      as.data.frame() %>%
      mutate(cell_pair = factor(cell_pair,  levels = unique(cell_pair))) 
    sub_dat2 <-
      dplyr::filter(
        dat,
        group_ligand %in% cell2,
        group_receptor %in% cell1
      )%>%
      mutate(group_receptor = factor(group_receptor, levels = cell1)) %>%
      arrange(group_receptor) %>%
      as.data.frame() %>%
      mutate(cell_pair = factor(cell_pair,  levels = unique(cell_pair))) 
    # return
    list(
      ligand   = sub_dat1,
      receptor = sub_dat2
    )
  } %>% 
  imap(function(dat, name){ # browser()
    dat %>% 
      as.data.frame() %>% 
      
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
        title = str_glue("{name}"),
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
          hjust = ifelse(name == "ligand", 0, 1),
          vjust = 0.5
        ),
        legend.position = "right"
      )
  })

# export
pdf(
  file.path(cpdb_dir, str_glue("dotplot - user-defined3.pdf")),
  height = 8,
  width  = 50
)
p_dot_ls_userdefined_3 %>% walk(print)
dev.off()

## 4. Macrophage各亚群：其他所有细胞类型 ----
p_dot_ls_userdefined_4 <-
  dot_dat_ls %>%
  enframe(name = "group") %>%
  unnest(value) %>%
  mutate(group_ligand = str_remove(cell_pair, "\\|.*$")) %>%
  mutate(group_receptor = str_remove(cell_pair, "^.*\\|")) %>% {
    dat <- .
    cell1 <- c('Macrophage')
    cell2 <- df_cell_major_to_sub
    sub_dat1 <-
      dplyr::filter(
        dat,
        group_ligand %in% cell1,
        group_receptor %in% cell2
      ) %>%
      mutate(group_ligand = factor(group_ligand, levels = cell1)) %>%
      arrange(group_ligand) %>%
      as.data.frame() %>%
      mutate(cell_pair = factor(cell_pair,  levels = unique(cell_pair))) 
    sub_dat2 <-
      dplyr::filter(
        dat,
        group_ligand %in% cell2,
        group_receptor %in% cell1
      )%>%
      mutate(group_receptor = factor(group_receptor, levels = cell1)) %>%
      arrange(group_receptor) %>%
      as.data.frame() %>%
      mutate(cell_pair = factor(cell_pair,  levels = unique(cell_pair))) 
    # return
    list(
      ligand   = sub_dat1,
      receptor = sub_dat2
    )
  } %>% 
  imap(function(dat, name){ # browser()
    dat %>% 
      as.data.frame() %>% 
      
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
        title = str_glue("{name}"),
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
          hjust = ifelse(name == "ligand", 0, 1),
          vjust = 0.5
        ),
        legend.position = "right"
      )
  })

# export
pdf(
  file.path(cpdb_dir, str_glue("dotplot - user-defined4.pdf")),
  height = 8,
  width  = 25
)
p_dot_ls_userdefined_4 %>% walk(print)
dev.off()

