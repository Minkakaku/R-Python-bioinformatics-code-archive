detach("package:monocle")
# detach("package:monocle3")


library(monocle3)
library(Seurat)
library(tidyverse)

outdir_monocle3 <- "monocle3 - v4"
dir.create(outdir_monocle3)


all_data_full <- readRDS('all_data_full_v6.rds')


# antrum ----------------
cds3_an <- readRDS('monocle3 - v1 - gastric antrum/monocle3.rds')

meta_cds3_an <- 
  data.frame(cell_id = colnames(cds3_an)) %>%
  left_join(
    rownames_to_column(all_data_full[[c('sub_celltype', 'sub_celltype3')]], "cell_id"),
    by = 'cell_id'
  )
table(is.na(meta_cds3_an$sub_celltype))
stopifnot(all(meta_cds3_an$cell_id == colnames(cds3_an)))

pData(cds3_an)$sub_celltype  <- meta_cds3_an$sub_celltype %>% fct_drop()
pData(cds3_an)$sub_celltype3 <- meta_cds3_an$sub_celltype3 %>% fct_drop()

remove_idx_an <- (!is.na(cds3_an$sub_celltype3)) & (cds3_an$sub_celltype3 != "endocrine_cells")

cds3_an_p1 <- plot_cells(
  cds3_an,
  reduction_method    = "UMAP",
  color_cells_by      = "sub_celltype",
  group_cells_by      = "sub_celltype",
  show_trajectory_graph = T,
  label_cell_groups     = T,
  label_leaves          = F,
  label_branch_points   = F,
  trajectory_graph_segment_size = 1,
  group_label_size = 5,  
  graph_label_size = 2.5, 
  cell_size = 0.3  
)

cds3_an_p1_2 <-
  cds3_an_p1 +
  labs(
    title = "Monocle3",
    x = "harmonyumap_1",
    y = "harmonyumap_2"
  ) +
  # scale_color_manual(values = fibro_cols_sc) +
  theme(
    plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )
cds3_an_p1_2$layers[[6]]$data <- cds3_an_p1_2$layers[[6]]$data %>% dplyr::filter(cell_group != "endocrine_cells")


ggsave(
  cds3_an_p1_2,
  filename = file.path(outdir_monocle3, "monocle trajectory - antrum - annotation.pdf"),
  width  = 8,
  height = 6
)


cds3_an_p2 <- plot_cells(
  cds3_an,
  reduction_method    = "UMAP",
  color_cells_by      = "sub_celltype3",
  group_cells_by      = "sub_celltype3",
  show_trajectory_graph = T,
  label_cell_groups     = T,
  label_leaves          = F,
  label_branch_points   = F,
  trajectory_graph_segment_size = 1,
  group_label_size = 5,  
  graph_label_size = 2.5, 
  cell_size = 0.3  
)

cds3_an_p2_2 <-
  cds3_an_p2 +
  labs(
    title = "Monocle3",
    x = "harmonyumap_1",
    y = "harmonyumap_2"
  ) +
  # scale_color_manual(values = fibro_cols_sc) +
  theme(
    plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )
cds3_an_p2_2$layers[[6]]$data <- cds3_an_p2_2$layers[[6]]$data %>% dplyr::filter(cell_group != "endocrine_cells")

ggsave(
  cds3_an_p2_2,
  filename = file.path(outdir_monocle3, "monocle trajectory - antrum - annotation - v2.pdf"),
  width  = 8,
  height = 6
)

# body ----------------
cds3_body <- readRDS('monocle3 - v1 - gastric body/monocle3.rds')

meta_cds3_body <- 
  data.frame(cell_id = colnames(cds3_body)) %>%
  left_join(
    rownames_to_column(all_data_full[[c('sub_celltype', 'sub_celltype3')]], "cell_id"),
    by = 'cell_id'
  )
table(is.na(meta_cds3_body$sub_celltype))
stopifnot(all(meta_cds3_body$cell_id == colnames(cds3_body)))

pData(cds3_body)$sub_celltype  <- meta_cds3_body$sub_celltype %>% fct_drop()
pData(cds3_body)$sub_celltype3 <- meta_cds3_body$sub_celltype3 %>% fct_drop()

# remove_idx_body1 <- (!is.na(cds3_body$sub_celltype)) & (cds3_body$sub_celltype != "endocrine_cells")
# remove_idx_body2 <- (!is.na(cds3_body$sub_celltype3)) & (cds3_body$sub_celltype3 != "endocrine_cells")


cds3_body_p1 <- plot_cells(
  cds3_body,
  reduction_method    = "UMAP",
  color_cells_by      = "sub_celltype",
  group_cells_by      = "sub_celltype",
  show_trajectory_graph = T,
  label_cell_groups     = T,
  label_leaves          = F,
  label_branch_points   = F,
  trajectory_graph_segment_size = 1,
  group_label_size = 5,  
  graph_label_size = 2.5, 
  cell_size = 0.3  
)

cds3_body_p1_2 <-
  cds3_body_p1 +
  labs(
    title = "Monocle3",
    x = "harmonyumap_1",
    y = "harmonyumap_2"
  ) +
  # scale_color_manual(values = fibro_cols_sc) +
  theme(
    plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )
# cds3_body_p1_2$data <- cds3_body_p1_2$data[remove_idx_body1, ]
cds3_body_p1_2$layers[[6]]$data <- cds3_body_p1_2$layers[[6]]$data %>% dplyr::filter(cell_group != "endocrine_cells")

ggsave(
  cds3_body_p1_2,
  filename = file.path(outdir_monocle3, "monocle trajectory - body - annotation.pdf"),
  width  = 8,
  height = 6
)


cds3_body_p2 <- plot_cells(
  cds3_body,
  reduction_method    = "UMAP",
  color_cells_by      = "sub_celltype3",
  group_cells_by      = "sub_celltype3",
  show_trajectory_graph = T,
  label_cell_groups     = T,
  label_leaves          = F,
  label_branch_points   = F,
  trajectory_graph_segment_size = 1,
  group_label_size = 5,  
  graph_label_size = 2.5, 
  cell_size = 0.3  
)

cds3_body_p2_2 <-
  cds3_body_p2 +
  labs(
    title = "Monocle3",
    x = "harmonyumap_1",
    y = "harmonyumap_2"
  ) +
  # scale_color_manual(values = fibro_cols_sc) +
  theme(
    plot.title = element_text(hjust = 0.5, size =  15, face = "bold"),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black")
  )
cds3_body_p2_2$layers[[6]]$data <- cds3_body_p2_2$layers[[6]]$data %>% dplyr::filter(cell_group != "endocrine_cells")

ggsave(
  cds3_body_p2_2,
  filename = file.path(outdir_monocle3, "monocle trajectory - body - annotation - v2.pdf"),
  width  = 8,
  height = 6
)
