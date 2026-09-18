# data preparation -------

### read data ------
all_data <- readRDS("/mnt/NexData_file/home/liuran/projects/2_VisiumHD/20250630_FWSP20242387_hypopharyngeal/epi_subcluster/all_data/all_data_project.rds")

sc_epi <- readRDS("/mnt/NexData_file/home/liuran/projects/2_VisiumHD/20250630_FWSP20242387_hypopharyngeal/inferCNV_test/sc_epi.rds")

### 筛选细胞做inferCNV -------
all_data@meta.data$annotation_infercnv <- all_data@meta.data$annotation

all_data@meta.data$annotation_infercnv <- ifelse(
  all_data@meta.data$annotation_infercnv %in% c("Plasma", "B cell", "CD4+ T cell", "Macrophage", "CD8+ T cell"), # 将这些细胞合并成"Immune cells"作为参考细胞，避免单个细胞类型数目太少导致报错；
  "Immune cells", 
  all_data@meta.data$annotation_infercnv
  )
unique(all_data@meta.data$annotation_infercnv)

visim_ls4 <- all_data %>% SplitObject(split.by = "sample")

sample_names <- all_data@meta.data$sample %>% unique()
for (sn in sample_names) {
  saveRDS(visim_ls4[[sn]], 
          file = paste0("all_data_", sn, ".rds"))
}