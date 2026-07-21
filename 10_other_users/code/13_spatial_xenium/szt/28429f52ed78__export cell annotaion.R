library(tidyverse)
library(Seurat)

all_data_xenium_pcf <- readRDS('all_data_xenium_pcf-v1-final.rds')


all_data_xenium_pcf[[c('project', 'annotation_full2')]] %>%
  dplyr::rename(annotation = annotation_full2) %>% 
  rownames_to_column('cell_id') %>%
  group_by(project) %>%
  nest() %>%
  deframe() %>%
  
  iwalk(function(dat, name){
    write.csv(
      dat, 
      file = str_glue("celltype - {name}.csv"),
      row.names = FALSE
    )
  })


# 20250914 export v2
all_data_xenium_pcf[[c('project', 'annotation_full2_2')]] %>%
  dplyr::rename(annotation = annotation_full2_2) %>% 
  rownames_to_column('cell_id') %>%
  group_by(project) %>%
  nest() %>%
  deframe() %>%
  
  iwalk(function(dat, name){
    write.csv(
      dat, 
      file = str_glue("celltype - {name} - 20250915.csv"),
      row.names = FALSE
    )
  })


DimPlot(
  all_data_xenium_pcf,
  reduction = "umap",
  group.by = "annotation_full2_2",
  label = TRUE
) +
  NoLegend()

# ll <- list.files(pattern = "ID[10]*.csv") %>%
#   map(read.csv) %>%
#   bind_rows()
# ll <- ll %>% column_to_rownames("cell_id")
# ll$annotation_full2 <- ll$annotation
# ll$annotation <- NULL
# head(ll)
# 
# all_data_xenium_pcf$annotation_full2_2 <- all_data_xenium_pcf$annotation_full2
# all_data_xenium_pcf <-
#   all_data_xenium_pcf %>%
#   AddMetaData(ll)



## only used xenium data
all_data <- readRDS('all_data_projected_subcluster.rds')

all_data[[c('project', 'annotation_full2')]] %>%
  dplyr::rename(annotation = annotation_full2) %>% 
  rownames_to_column('cell_id') %>%
  group_by(project) %>%
  nest() %>%
  deframe() %>%
  
  iwalk(function(dat, name){
    write.csv(
      dat, 
      file = str_glue("celltype - {name} - only xenium - v2.csv"),
      row.names = FALSE
    )
  })


## only used xenium data - 20250914
all_data <- readRDS('all_data_projected_subcluster_update250914.rds')

all_data[[c('project', 'annotation_full2')]] %>%
  dplyr::rename(annotation = annotation_full2) %>% 
  rownames_to_column('cell_id') %>%
  group_by(project) %>%
  nest() %>%
  deframe() %>%
  
  iwalk(function(dat, name){
    write.csv(
      dat, 
      file = str_glue("celltype - {name} - only xenium - 20250915.csv"),
      row.names = FALSE
    )
  })


## only used pcf data
all_data_pcf <- readRDS('all_data.pcf.rds')

all_data_pcf[[c('project', 'annotation')]] %>%
  # dplyr::rename(annotation = annotation_full2) %>% 
  rownames_to_column('cell_id') %>%
  group_by(project) %>%
  nest() %>%
  deframe() %>%
  
  iwalk(function(dat, name){
    write.csv(
      dat, 
      file = str_glue("celltype - {name} - only pcf.csv"),
      row.names = FALSE
    )
  })

