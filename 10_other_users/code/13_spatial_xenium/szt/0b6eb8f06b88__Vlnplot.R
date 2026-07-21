
# cell_id <- all_data_ls %>% 
#   map(function(all_data) {
#     # 提取列名
#     colnames(all_data)
#   }) %>% 
#   unlist() 
# length(cell_id)
# length(unique(cell_id))

library(Seurat)
library(tidyverse)
library(ggplot2)

saveRDS(all_data, "all_data_annotated.rds")
saveRDS(Duct_cells_1, "Duct_cell_1.rds")

VinPlot <- "VinPlot_0627"
dir.create(VinPlot)

# VlnPlot(Duct_cells_1,c('Csf1','Vegfa'),group.by = 'group')
# Idents(all_data)<- 'annotation'
# VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1')


all_data <- NormalizeData(all_data, scale.factor = median(all_data$nCount_RNA))
# VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1')
VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1',pt.size = 0,raster = FALSE)
# VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1',pt.size = 0,raster = FALSE,y.max = 2)

# Ldha,Ldhb,Eno1,Eno2,Pkm,Pfkp,Hk2,Slc16a1,Slc16a3(使用小提琴图可视化)
VlnPlot(all_data, features = "Ldha", group.by = 'group', idents = 'Duct_cells_1', pt.size = 0, raster = FALSE)
ggsave("VinPlot_0627/Vlnplot_Ldha.pdf", width=9, height=6)

VlnPlot(all_data, features = "Ldhb", group.by = 'group', idents = 'Duct_cells_1', pt.size = 0, raster = FALSE)
ggsave("VinPlot_0627/Vlnplot_Ldhb.pdf", width=9, height=6)

VlnPlot(all_data, features = "Eno1", group.by = 'group', idents = 'Duct_cells_1', pt.size = 0, raster = FALSE)
ggsave("VinPlot_0627/Vlnplot_Eno1.pdf", width=9, height=6)

VlnPlot(all_data, features = "Eno2", group.by = 'group', idents = 'Duct_cells_1', pt.size = 0, raster = FALSE)
ggsave("VinPlot_0627/Vlnplot_Eno2.pdf", width=9, height=6)

VlnPlot(all_data, features = "Pkm", group.by = 'group', idents = 'Duct_cells_1', pt.size = 0, raster = FALSE)
ggsave("VinPlot_0627/Vlnplot_Pkm.pdf", width=9, height=6)

VlnPlot(all_data, features = "Pfkp", group.by = 'group', idents = 'Duct_cells_1', pt.size = 0, raster = FALSE)
ggsave("VinPlot_0627/Vlnplot_Pfkp.pdf", width=9, height=6)

VlnPlot(all_data, features = "Hk2", group.by = 'group', idents = 'Duct_cells_1', pt.size = 0, raster = FALSE)
ggsave("VinPlot_0627/Vlnplot_Hk2.pdf", width=9, height=6)

VlnPlot(all_data, features = "Slc16a1", group.by = 'group', idents = 'Duct_cells_1', pt.size = 0, raster = FALSE)
ggsave("VinPlot_0627/Vlnplot_Slc16a1.pdf", width=9, height=6)

VlnPlot(all_data, features = "Slc16a3", group.by = 'group', idents = 'Duct_cells_1', pt.size = 0, raster = FALSE)
ggsave("VinPlot_0627/Vlnplot_Slc16a3.pdf", width=9, height=6)




# p<- VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1',pt.size = 0,raster = FALSE)
# plot_data <- VlnPlot(all_data,c('Csf1','Vegfa'),group.by = 'group',idents = 'Duct_cells_1',pt.size = 0,raster = FALSE,return.data = TRUE)
# expression_data <- all_data[which(meta_data$group == "Duct_cells_1"), c("Csf1", "Vegfa")]
# expression_data <- all_data[which(all_data@meta.data$group == "Duct_cells_1"), c("Csf1", "Vegfa")]
# expression_data <- all_data[which(all_data$annotation == "Duct_cells_1"), c("Csf1", "Vegfa")]


subset(all_data,annotation=='Duct_cells_1') %>% FetchData(vars = c("Csf1", "Vegfa",slot = 'data')) -> temp
subset(all_data,annotation=='Duct_cells_1') %>% FetchData(vars = c("Ldha","Ldhb","Eno1","Eno2","Pkm","Pfkp","Hk2","Slc16a1","Slc16a3",slot = 'data')) -> temp1

head(temp)
head(temp1)
Duct_cells_1 @meta.data %>% colnames()
head(Duct_cells_1)
Duct_cells_1


temp %>%
  rownames_to_column('cell_id') %>%
  left_join(Duct_cells_1[[]] %>%
              select(group) %>%
              column_to_rownames('cell_id'),by = 'cell_id')

temp1 %>% rownames_to_column('cell_id') %>%
  left_join(Duct_cells_1[[]] %>%
              select(group) %>%
              column_to_rownames('cell_id'),by = 'cell_id')


temp %>% rownames_to_column('cell_id')
Duct_cells_1[[]] %>% select(group)
Duct_cells_1[[]] %>% select(group) %>% 
  column_to_rownames('cell_id')

Duct_cells_1[[]] %>% select(group) %>% head()

temp1 %>% rownames_to_column('cell_id') %>% 
  left_join(Duct_cells_1[[]] %>% 
              select(group) %>% 
              rownames_to_column('cell_id'),by = 'cell_id')
getwd()
temp %>%
  rownames_to_column('cell_id') %>%
  left_join(Duct_cells_1[[]] %>%
              select(group) %>%
              rownames_to_column('cell_id'),by = 'cell_id') %>%
  column_to_rownames('cell_id') %>%
  write.csv('vlnplot_data.csv')

temp1 %>%
  rownames_to_column('cell_id') %>%
  left_join(Duct_cells_1[[]] %>%
              select(group) %>%
              rownames_to_column('cell_id'),by = 'cell_id') %>%
  column_to_rownames('cell_id') %>%
  write.csv('VinPlot_0627/vlnplot_data_20250627.csv')





