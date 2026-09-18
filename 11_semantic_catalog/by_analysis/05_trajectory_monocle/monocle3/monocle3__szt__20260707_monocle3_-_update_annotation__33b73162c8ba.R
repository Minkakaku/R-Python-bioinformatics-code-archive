library(monocle3)
library(Seurat)
library(tidyverse)

out_dir_260707 <- "20260707 monocle3 - annov7"
dir.create(out_dir_260707)

all_data_full <- readRDS('all_data_full_v7.rds')
meta_roi_full <- all_data_full[[]]

# 1. gastric body ------------
cds3_body <- readRDS("monocle3 - v1 - gastric body/monocle3.rds")
stopifnot(all(colnames(cds3_body) %in% colnames(all_data_full)))

pData(cds3_body)$sub_celltype_v7  <- meta_roi_full$sub_celltype_v7[match(colnames(cds3_body), rownames(meta_roi_full))]
pData(cds3_body)$sub_celltype2_v7 <- meta_roi_full$sub_celltype2_v7[match(colnames(cds3_body), rownames(meta_roi_full))]

## sub_celltype -------
p_cds3_body <- plot_cells(
  cds3_body,
  reduction_method    = "UMAP",
  color_cells_by      = "sub_celltype_v7",
  group_cells_by      = "sub_celltype_v7",
  show_trajectory_graph = T,
  label_cell_groups     = T,
  label_leaves          = F,
  label_branch_points   = F,
  trajectory_graph_segment_size = 1,
  group_label_size = 5,  
  graph_label_size = 2.5, 
  cell_size = 0.3  
)

p_cds3_body_2 <-
  p_cds3_body +
  labs(
    title = "Monocle3",
    x = "harmonyumap_1",
    y = "harmonyumap_2"
  ) +
  ggsci::scale_color_igv() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )

sub_celltype_epi_levs <- 
  p_cds3_body_2$data$sub_celltype_v7 %>%
  unique() %>% 
  as.character() %>%
  sort()  
  # setdiff('Tumor_C1_Stem_like') %>%
  # c("Tumor_C1_Stem_like")

p_cds3_body_2$data <-
  p_cds3_body_2$data %>%
  dplyr::filter(!is.na(sub_celltype_v7)) %>%
  mutate(sub_celltype_v7 = factor(as.character(sub_celltype_v7), levels = sub_celltype_epi_levs)) %>%
  arrange(sub_celltype_v7)

ggsave(
  p_cds3_body_2,
  filename = file.path(out_dir_260707, "1. gastric body - monocle trajectory - annotation - v7 - sub_celltype.pdf"),
  width  = 8,
  height = 6
)

## sub_celltype2 -------
p_cds3_body2 <- plot_cells(
  cds3_body,
  reduction_method    = "UMAP",
  color_cells_by      = "sub_celltype2_v7",
  group_cells_by      = "sub_celltype2_v7",
  show_trajectory_graph = T,
  label_cell_groups     = T,
  label_leaves          = F,
  label_branch_points   = F,
  trajectory_graph_segment_size = 1,
  group_label_size = 5,  
  graph_label_size = 2.5, 
  cell_size = 0.3  
)

p_cds3_body2_2 <-
  p_cds3_body2 +
  labs(
    title = "Monocle3",
    x = "harmonyumap_1",
    y = "harmonyumap_2"
  ) +
  ggsci::scale_color_igv() +
  theme(
    plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    legend.position = "right"
  )


p_cds3_body2_2$data <-
  p_cds3_body2_2$data %>%
  dplyr::filter(!is.na(sub_celltype2_v7))

ggsave(
  p_cds3_body2_2,
  filename = file.path(out_dir_260707, "1. gastric body - monocle trajectory - annotation - v7 - sub_celltype2.pdf"),
  width  = 8,
  height = 6
)

## 胃体246283、255670 ----------
p_cds3_body_2sample_ls <-
  c('A246283', 'A255670', 'combined') %>%
  setNames(.,.) %>%
  map(function(x){
    p <- p_cds3_body_2
    if(x == "combined"){
      p$data <-
        p$data %>%
        dplyr::filter(project %in% c('A246283', 'A255670') )
    }else{
      p$data <-
        p$data %>%
        dplyr::filter(project == x)
    }
    
    p + ggtitle(x)
  })


pdf(
  file = file.path(out_dir_260707, "1. gastric body - monocle trajectory - annotation - v7 - sub_celltype - 2samples.pdf"),
  width  = 8,
  height = 6
)
p_cds3_body_2sample_ls %>% walk(print)
dev.off()

## sub_celltype2 
p_cds3_body2_2sample_ls <-
  c('A246283', 'A255670', 'combined') %>%
  setNames(.,.) %>%
  map(function(x){
    p <- p_cds3_body2_2
    if(x == "combined"){
      p$data <-
        p$data %>%
        dplyr::filter(project %in% c('A246283', 'A255670') )
    }else{
      p$data <-
        p$data %>%
        dplyr::filter(project == x)
    }
    
    p + ggtitle(x)
  })


pdf(
  file = file.path(out_dir_260707, "1. gastric body - monocle trajectory - annotation - v7 - sub_celltype2 - 2samples.pdf"),
  width  = 8,
  height = 6
)
p_cds3_body2_2sample_ls %>% walk(print)
dev.off()




# 2. gastric body ------------
cds3_antrum <- readRDS("monocle3 - v1 - gastric antrum/monocle3.rds")
stopifnot(all(colnames(cds3_antrum) %in% colnames(all_data_full)))

pData(cds3_antrum)$sub_celltype_v7  <- meta_roi_full$sub_celltype_v7[match(colnames(cds3_antrum), rownames(meta_roi_full))]
pData(cds3_antrum)$sub_celltype2_v7 <- meta_roi_full$sub_celltype2_v7[match(colnames(cds3_antrum), rownames(meta_roi_full))]


## sub_celltype -------
p_cds3_antrum <- plot_cells(
  cds3_antrum,
  reduction_method    = "UMAP",
  color_cells_by      = "sub_celltype_v7",
  group_cells_by      = "sub_celltype_v7",
  show_trajectory_graph = T,
  label_cell_groups     = T,
  label_leaves          = F,
  label_branch_points   = F,
  trajectory_graph_segment_size = 1,
  group_label_size = 5,  
  graph_label_size = 2.5, 
  cell_size = 0.3  
)

p_cds3_antrum_2 <-
  p_cds3_antrum +
  labs(
    title = "Monocle3",
    x = "harmonyumap_1",
    y = "harmonyumap_2"
  ) +
  ggsci::scale_color_igv() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )

sub_celltype_epi_levs2 <- 
  p_cds3_antrum_2$data$sub_celltype_v7 %>%
  unique() %>% 
  as.character() %>%
  sort()  


p_cds3_antrum_2$data <-
  p_cds3_antrum_2$data %>%
  dplyr::filter(!is.na(sub_celltype_v7)) %>%
  mutate(sub_celltype_v7 = factor(as.character(sub_celltype_v7), levels = sub_celltype_epi_levs2)) %>%
  arrange(sub_celltype_v7)

ggsave(
  p_cds3_antrum_2,
  filename = file.path(out_dir_260707, "2. gastric antrum - monocle trajectory - annotation - v7 - sub_celltype.pdf"),
  width  = 8,
  height = 6
)

## sub_celltype2 -------
p_cds3_antrum2 <- plot_cells(
  cds3_antrum,
  reduction_method    = "UMAP",
  color_cells_by      = "sub_celltype2_v7",
  group_cells_by      = "sub_celltype2_v7",
  show_trajectory_graph = T,
  label_cell_groups     = T,
  label_leaves          = F,
  label_branch_points   = F,
  trajectory_graph_segment_size = 1,
  group_label_size = 5,  
  graph_label_size = 2.5, 
  cell_size = 0.3  
)

p_cds3_antrum2_2 <-
  p_cds3_antrum2 +
  labs(
    title = "Monocle3",
    x = "harmonyumap_1",
    y = "harmonyumap_2"
  ) +
  ggsci::scale_color_igv() +
  theme(
    plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    legend.position = "right"
  )


p_cds3_antrum2_2$data <-
  p_cds3_antrum2_2$data %>%
  dplyr::filter(!is.na(sub_celltype2_v7))

ggsave(
  p_cds3_antrum2_2,
  filename = file.path(out_dir_260707, "2. gastric antrum - monocle trajectory - annotation - v7 - sub_celltype2.pdf"),
  width  = 8,
  height = 6
)

## 胃窦234244 ----------
p_cds3_antrum_1sample_ls <-
  c('A234244') %>%
  setNames(.,.) %>%
  map(function(x){
    p <- p_cds3_antrum_2

    p$data <-
      p$data %>%
      dplyr::filter(project == x)
    
    
    p + ggtitle(x)
  })


pdf(
  file = file.path(out_dir_260707, "2. gastric antrum - monocle trajectory - annotation - v7 - sub_celltype - 1samples.pdf"),
  width  = 8,
  height = 6
)
p_cds3_antrum_1sample_ls %>% walk(print)
dev.off()


## sub_celltype2 
p_cds3_antrum2_1sample_ls <-
  c('A234244') %>%
  setNames(.,.) %>%
  map(function(x){
    p <- p_cds3_antrum2_2

    p$data <-
      p$data %>%
      dplyr::filter(project == x)
    
    
    p + ggtitle(x)
  })


pdf(
  file = file.path(out_dir_260707, "2. gastric antrum - monocle trajectory - annotation - v7 - sub_celltype2 - 1samples.pdf"),
  width  = 8,
  height = 6
)
p_cds3_antrum2_1sample_ls %>% walk(print)
dev.off()



