
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


# export cellphonedb files   -----
## read data  -------
all_data <- readRDS('all_data_fb_rctd_1223.rds')
# cells <- !is.na(all_data@meta.data$annotation_fbsub)
# all_data_filt <- all_data[, cells]
all_data_filt <- all_data[, all_data@meta.data$annotation_fbsub != "Unknown"]

# dir.create('cellphonedb_result3')
all_data_meta <- all_data_filt[[]]
all_data_filt <- JoinLayers(all_data_filt)
all_data_norm <- LayerData(all_data_filt, layer ='counts')


human_mouse<-read_tsv('human_mouse_id.txt')
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
    
    out_dir <- file.path("cellphonedb_result_fb", 'in', name)
    dir.create(out_dir,recursive = TRUE)
    
    # write meta files
    meta['annotation_fbsub'] %>%
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








# # download cellphone db files ------
# reticulate::py_run_string("
# import os
# from cellphonedb.utils import db_utils
# 
# # Version of the databse
# cpdb_version = 'v5.0.0'
# 
# # Path where the input files to generate the database are located
# cpdb_target_dir = os.path.join('cellphonedb/database', cpdb_version)
# db_utils.download_database(cpdb_target_dir, cpdb_version)
# 
# ")



# running cellphondb ------
cellphonedb_path <- "cellphonedb/database/v5.0.0/cellphonedb.zip"

## ################### python code start ########################

repl_python() ## start python
import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
db_path = "cellphonedb/database/v5.0.0/cellphonedb.zip" 


in_dir = os.path.join('cellphonedb_result_fb', 'in')
os.makedirs(in_dir, exist_ok = True)
out_dir = os.path.join('cellphonedb_result_fb', 'out')
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
  threads = 32
)

quit

## ################### python code end. ########################


# plot cellphonedb ------------------
cpdb_plot_out_dir <- "cellphonedb_result_fb"

all_files_ls <- 
  list.files(file.path(cpdb_plot_out_dir, "out")) %>%
  set_names(., basename(.)) %>%
  map(function(x){
    all_files <- list.files(str_glue('{cpdb_plot_out_dir}/out/{x}/'), full.names = TRUE)
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

## split dotplot ------
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

## heatmap --------
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
  # separate(group, into = c("Group", "group"), sep = "_") %>%
  # dplyr::select(-group) %>%
  group_by(group) %>%
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
  file.path(cpdb_plot_out_dir, str_glue("dotplot.pdf")),
  height = 30,
  width  = 50
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






p_dot_ls4 <-
  dot_dat_ls %>%
  enframe(name = "group") %>%
  unnest(value) %>%
  mutate(
    ligand = str_extract(interacting_pair, "^[^_]+"), 
    group_ligand = str_remove(cell_pair, "\\|.*$"), 
    group_receptor = str_remove(cell_pair, "^.*\\|")) %>%
  # filter(group_ligand == "Fibroblasts", group_receptor %in% c("Macrophage", "DC"), ligand %in% c("VEGFA", "TGFB1")) %>% 
  filter(group_ligand == "Hypoxia_CAFs", group_receptor %in% c("Macrophage", "DC")) %>% 
  
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
          angle = 30,
          hjust = 1,
          vjust = 1
        ),
        legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 22)
      )
  })

pdf(
  file.path(cpdb_plot_out_dir, str_glue("dotplot_4.pdf")),
  height = 15,
  width  = 12
)
p_dot_ls4 %>% walk(print)
dev.off()





p_dot_ls5 <-
  dot_dat_ls %>%
  enframe(name = "group") %>%
  unnest(value) %>%
  mutate(
    ligand = str_extract(interacting_pair, "^[^_]+"), 
    group_ligand = str_remove(cell_pair, "\\|.*$"), 
    group_receptor = str_remove(cell_pair, "^.*\\|")) %>%
  filter(group_ligand == "Hypoxia_CAFs", group_receptor %in% c("Macrophage", "DC"), ligand %in% c("VEGFA", "TGFB1")) %>%
  # filter(group_ligand == "Hypoxia_CAFs", group_receptor %in% c("Macrophage", "DC")) %>% 
  
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
          angle = 30,
          hjust = 1,
          vjust = 1
        ),
        legend.position = "right", 
        plot.title = element_text(hjust = 0.5, size = 22)
      )
  })

pdf(
  file.path(cpdb_plot_out_dir, str_glue("dotplot_5.pdf")),
  height = 6,
  width  = 12
p_dot_ls5 %>% walk(print)










