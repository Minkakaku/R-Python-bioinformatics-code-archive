library(Seurat)
library(tidyverse)

source("spatial analysis.R")

proximity_out_dir <- "proximity_out_dir"
dir.create(proximity_out_dir)


# read data -------
meta_data <- data.table::fread("proximity_out_dir/adata_cellty_obs.csv", data.table = FALSE)
coords_df <- data.table::fread("proximity_out_dir/coords.csv", data.table = FALSE)

meta_data <- cbind(meta_data[-1], coords_df[-1])

# implement ---------
# options(warn = 0)
#  
# Classification1 must be a factor
# add a column: obj_id



all_data_meta_for_cn_analysis <-
  meta_data %>%
  dplyr::mutate(annotation_full = celltype) %>% 
  dplyr::mutate(cell_id = str_c(sample, cell_id, sep = "_")) %>%
  rename(X = X, Y = Y) %>%
  mutate(
    annotation_full = factor(annotation_full), # import!
    obj_id = cell_id # import!
  )


## NN: 10 cells ------------------
all_cell_counts_ls <-
  all_data_meta_for_cn_analysis %>%
  group_by(sample) %>%
  nest() %>%
  deframe() %>%
  map(function(dat){
    CN_compute_cellfreq(
      dat,
      x = 'X',
      y = 'Y',
      classification = "annotation_full",
      return.vars = "count",
      k.nn = 10
    )
  })

all_cell_counts_df <-
  all_cell_counts_ls %>%
  map(t) %>%
  map(as.data.frame) %>%
  bind_rows() %>%
  rownames_to_column("cell_id")
  

all_cell_counts_df <- 
  all_cell_counts_df %>%
  left_join(all_data_meta_for_cn_analysis[c("cell_id", "sample", "celltype")], by = "cell_id")

all_cell_counts_summary_df <-
  all_cell_counts_df %>%
  dplyr::select(-cell_id) %>%
  group_by(sample, celltype) %>%
  summarise(across(everything(), median)) %>%
  dplyr::rename(Target = celltype)


p_proximity_box <-
  all_cell_counts_summary_df %>%
  pivot_longer(
    -c(sample, Target),
    names_to  = "Others",
    values_to = "count"
  ) %>%

  ggplot(aes(x = Others, y = count)) +
  geom_boxplot(aes(fill = Others), outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target) +
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggsci::scale_fill_aaas() +
  labs(
    caption = "Median count of other cells surrounding the target cells per sample"
  )

# export
openxlsx::write.xlsx(
  all_cell_counts_summary_df,
  file = file.path("proximity_out_dir/Median cell count - nn10.xlsx"),
  asTable = TRUE
)
ggsave(
  p_proximity_box,
  filename = file.path(proximity_out_dir, "boxplot of others cells count surrounding the target cells - nn10.pdf"),
  width = 15,
  height = 12
)


## NN:30 cells ------------------
all_cell_counts_ls_nn30 <-
  all_data_meta_for_cn_analysis %>%
  group_by(sample) %>%
  nest() %>%
  deframe() %>%
  map(function(dat){
    CN_compute_cellfreq(
      dat,
      x = 'X',
      y = 'Y',
      classification = "annotation_full",
      return.vars = "count",
      k.nn = 30
    )
  })


all_cell_counts_df_nn30 <-
  all_cell_counts_ls_nn30 %>%
  map(t) %>%
  map(as.data.frame) %>%
  bind_rows() %>%
  rownames_to_column("cell_id")


all_cell_counts_df_nn30 <- 
  all_cell_counts_df_nn30 %>%
  left_join(all_data_meta_for_cn_analysis[c("cell_id", "sample", "celltype")], by = "cell_id")

all_cell_counts_summary_df_nn30 <-
  all_cell_counts_df_nn30 %>%
  dplyr::select(-cell_id) %>%
  group_by(sample, celltype) %>%
  summarise(across(everything(), median)) %>%
  dplyr::rename(Target = celltype)


p_proximity_box_nn30 <-
  all_cell_counts_summary_df_nn30 %>%
  pivot_longer(
    -c(sample, Target),
    names_to  = "Others",
    values_to = "count"
  ) %>%
  
  ggplot(aes(x = Others, y = count)) +
  geom_boxplot(aes(fill = Others), outlier.colour = NA) +
  geom_jitter(width = 0.25) +
  facet_wrap(~ Target) +
  ggpubr::theme_pubclean() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggsci::scale_fill_aaas() +
  labs(
    caption = "Median count of other cells surrounding the target cells per sample"
  )

# export
openxlsx::write.xlsx(
  all_cell_counts_summary_df_nn30,
  file = file.path("proximity_out_dir/Median cell count - nn30.xlsx"),
  asTable = TRUE
)
ggsave(
  p_proximity_box_nn30,
  filename = file.path(proximity_out_dir, "boxplot of others cells count surrounding the target cells - nn30.pdf"),
  width = 15,
  height = 12
)




