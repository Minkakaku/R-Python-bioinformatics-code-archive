library(Seurat)
library(tidyverse)

result_dir <- 'result_8um/inferCNVPy'  # dir.create(result_dir)
Luminal4 <- readRDS('result_8um/Luminal_reannotation3/Luminal4_final.rds')
adata_meta <- read_csv('result_8um/Luminal_reannotation3/adata_manual_obs_res1.csv')
adata_meta <- 
  adata_meta %>% 
  dplyr::rename(cell_id = ...1)

stopifnot(adata_meta['cell_id'] == Cells(Luminal4))

Luminal4_meta <- 
  Luminal4@meta.data %>%
  rownames_to_column('cell_id') %>% 
  left_join(adata_meta, by = 'cell_id') %>% 
  column_to_rownames('cell_id') %>% 
  select(cnv_leiden)

Luminal4 <- AddMetaData(Luminal4, Luminal4_meta['cnv_leiden'])




# spatial plot ------
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
names(allcolour) <- Luminal4$cnv_leiden %>% unique() %>% sort()

ggsave(
  SpatialDimPlot(Luminal4, group.by = 'cnv_leiden', pt.size.factor = 6, image.scale = 'hires', cols = allcolour, shape = 22, label = T, label.size = 3),
  filename = file.path(result_dir, "Saptial_luminal_type_by_cnv_leiden_res1.pdf"),
  height = 10,
  width  = 20
)




# Annotation info output -------------------------------------------------------
Luminal4$cnv_leiden %>% table()
cellType_anno_Luminal4 <- 
  Luminal4@meta.data[c('project', 'cnv_leiden')] %>% 
  rownames_to_column('cell_id') %>% 
  mutate(
    Barcode = str_replace(cell_id, '^.*?_', ''), 
    group = cnv_leiden
  ) %>% 
  dplyr::select(Barcode, group)

cellType_anno_Luminal4 %>% 
  write.csv(file.path(result_dir, 'luminal_type_by_cnv_leiden_res1.csv'), row.names = F)
