library(tidyverse)
library(Seurat)

all_data <- readRDS('all_data.epi.t.rds')

out_dir_export <- "export annotation for xenium explorer"
dir.create(out_dir_export)


all_data[[]] %>% 
  dplyr::select(
    fov, 
    cell_id = cell_id_sample, 
    group = annotation
  ) %>%
  mutate(cell_id = str_replace(cell_id, "-1-.*$", "-1")) %>%
  group_by(fov) %>%
  nest() %>%
  deframe() %>%
  iwalk(function(dat, name){
    write.csv(
      dat,
      file = file.path(out_dir_export, str_glue("{name}.csv")),
      row.names = FALSE
    )
  })
