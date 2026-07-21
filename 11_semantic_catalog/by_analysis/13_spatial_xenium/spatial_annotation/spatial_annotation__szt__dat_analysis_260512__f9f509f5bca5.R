# export cell annotation table

library(Seurat)
library(tidyverse)

all_data <- readRDS("dat_analysis_cycle3/subtype/all_data_cycle3_260417.rds")

all_data[[]] %>% head()

all_data[[c('group', 'group2', 'annotation6')]] %>%
  rownames_to_column("cell_id") %>%
  write.csv(
    file = file.path('dat_analysis_cycle3', '20260513 cell annotation - add CD8Tpex.csv'),
    row.names = FALSE
  )
