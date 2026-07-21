library(Seurat)
library(infercnv)
library(tidyverse)
library(dplyr)


setwd("/mnt/NexData_file/home/liuran/projects/2_VisiumHD/20250630_FWSP20242387_hypopharyngeal/inferCNV_test")

all_data <- readRDS("/mnt/NexData_file/home/liuran/projects/2_VisiumHD/20250630_FWSP20242387_hypopharyngeal/epi_subcluster/all_data/all_data_project.rds")
all_data_G5 <- readRDS("/mnt/NexData_file/home/liuran/projects/2_VisiumHD/20250630_FWSP20242387_hypopharyngeal/inferCNV_test/all_data_G5.rds")
sc_epi_G5 <- readRDS("/mnt/NexData_file/home/liuran/projects/2_VisiumHD/20250630_FWSP20242387_hypopharyngeal/inferCNV_test/sc_epi_G5.rds")



# 筛选细胞做inferCNV -------
# all_data@meta.data$annotation_infercnv <- all_data@meta.data$annotation

# all_data@meta.data$annotation_infercnv <- ifelse(
  # all_data@meta.data$annotation_infercnv %in% c("Plasma", "B cell", "CD4+ T cell", "Macrophage", "CD8+ T cell"), 
  # "Immune cells", 
  # all_data@meta.data$annotation_infercnv
# )
# unique(all_data@meta.data$annotation_infercnv)



# G5 ---------
sc_epi_G5_anno_meta <- 
  sc_epi_G5[['Spatial.008um_snn_res.0.2']] %>% 
  rownames_to_column("cell_id") 

all_data_G5_anno_meta <- 
  all_data_G5[[]] %>% dplyr::filter(nCount_Spatial.008um > 100) %>% .['annotation_infercnv'] %>% 
  # all_data[['annotation_infercnv']] %>% 
  rownames_to_column('cell_id') %>% 
  dplyr::filter(annotation_infercnv %in% c("Epithelial", "Immune cells")) %>% 
  left_join(sc_epi_G5_anno_meta, by = "cell_id") %>% 
  map_dfc(as.character) %>% 
  mutate(annotation2 = ifelse(
    annotation_infercnv %in% "Epithelial", 
    Spatial.008um_snn_res.0.2, 
    annotation_infercnv
  )) %>% 
  dplyr::select(cell_id, annotation = annotation2)


# all_data <- all_data %>% AddMetaData(all_data_anno_meta)  
meta_data_infer_subset <- 
  all_data_G5_anno_meta %>% 
  group_by(annotation) %>% 
  slice_sample(prop = 0.5)



### inferCNV ---------------
# used to infer epethial cells
# immunes cell used as control

# 当前路径为pipeline文件夹
library(infercnv)
infercnv_out_dir = "/mnt/NexData_file/home/liuran/projects/2_VisiumHD/20250630_FWSP20242387_hypopharyngeal/inferCNV_test"
infercnv_G5_out_dir <- file.path(infercnv_out_dir, "infercnv_G5")
dir.create(infercnv_G5_out_dir) 


# export as file
meta_data_infer_subset %>% 
  column_to_rownames("cell_id") %>% 
  data.table::fwrite(
    file = file.path(infercnv_G5_out_dir, "meta_data_infer_subset.txt"), 
    row.names = TRUE, 
    col.names = FALSE, 
    sep = "\t"
  )


# undebug(CreateInfercnvObject)
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix = JoinLayers(all_data)@assays$Spatial.008um$counts[VariableFeatures(all_data)[1:500],  meta_data_infer_subset$cell_id], 
  annotations_file  = file.path(infercnv_G5_out_dir, "meta_data_infer_subset.txt"), 
  gene_order_file   = file.path('/mnt/NexData_file/home/liuran/projects/2_VisiumHD/20250630_FWSP20242387_hypopharyngeal/infercnv2', "hsa_gene_position.txt"), 
  
  min_max_counts_per_cell = c(30, Inf), 
  
  delim             = "\t", 
  ref_group_names   = "Immune cells"
)



# all_data_G5@assays$RNA@counts %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

infercnv_obj@count.data %>% rowMeans %>% { sum(. > 0.001) } # 394 
infercnv_obj@count.data %>% rowMeans %>% { sum(. < 0.005) } # 199


infercnv_obj2  = infercnv::run(
  infercnv_obj, 
  cutoff              =    0.001, 
  out_dir             =    infercnv_G5_out_dir, 
  cluster_by_groups   =    TRUE, 
  plot_steps          =    FALSE, 
  denoise             =    TRUE, 
  HMM                 =    FALSE, 
  no_prelim_plot      =    TRUE, 
  png_res             =    300, 
  # useRaster           =   FALSE, 
  output_format       =   'pdf', 
  no_plot = T ### figure margin too large ----
)

infercnv::plot_cnv(
  infercnv_obj2, 
  out_dir           =  infercnv_G5_out_dir, 
  output_filename   =  'infercnv', 
  output_format     =  "pdf", 
  png_res           =  300, 
  cluster_by_groups =  TRUE, 
  x.range           =  "auto", 
  x.center          =  1, 
  title             =  "infercnv", 
  color_safe_pal    =  FALSE
)

saveRDS(infercnv_obj2, file.path(infercnv_G5_out_dir, 'infercnv_obj2.rds'))
saveRDS(infercnv_obj, file.path(infercnv_G5_out_dir, 'infercnv_obj.rds'))


infercnv_obj2_default_median_filtering <- 
  infercnv_obj2 %>% 
  infercnv::apply_median_filtering()

infercnv_obj2_default_median_filtering %>% 
  infercnv::plot_cnv(
    out_dir           =  infercnv_G5_out_dir, 
    output_filename   =  'infercnv_median_filtering', 
    output_format     = "pdf", 
    cluster_by_groups =  TRUE, 
    x.range           =  "auto", 
    x.center          =  1, 
    title             =  "infercnv", 
    color_safe_pal    =  FALSE
  )
saveRDS(infercnv_obj2_default_median_filtering, file.path(infercnv_G5_out_dir, 'infercnv_obj2_default_median_filtering .rds'))


### cnv score ------------
library(scales)
cnvScore <- function(data){
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  rng <- quantile(data, c(0.01, 0.99))
  data <- pmax( pmin(data, rng[2]), rng[1] )
  
  data <- data %>% as.matrix() %>% 
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

# cnvScore <- function(data){
#   sd_num <- apply(data, 1, sd)
#   data <- data[sd_num > median(sd_num), ]
#   
#   data <- data %>% as.matrix() %>%
#     t() %>% 
#     scale() %>% 
#     t()
#   
#   cnv_score <- as.data.frame(colSums(data * data))
#   return(cnv_score)
#   }
cnv_score <- cnvScore(infercnv_obj2@expr.data)
colnames(cnv_score) <-'CNV_score'
temp <- 
  cnv_score %>% 
  rownames_to_column('cell_id') %>% 
  left_join(.,as.data.frame(meta_data_infer_subset),'cell_id')

p_infer_cnv_score <- 
  temp %>% 
  ggplot(aes(annotation,log(CNV_score)))+geom_violin(aes(fill= annotation),color="NA",drop = FALSE) + 
  ggsci::scale_fill_igv() + 
  # scale_fill_manual(values = color_v) + 
  theme_bw() -> p

ggsave(
  p_infer_cnv_score, 
  filename = file.path(infercnv_G5_out_dir, "cnv_vlnplot.pdf"), 
  height = 6, 
  width  = 12
)

cols_brewer = c("#f94144","#f25c27","#eb770a","#f89f12","#f9c74f","#90be6d",
                "#43aa8b","#4d908e","#577590","#277da1","#967aa1",
                "#F2698D","#F391AE","#F3B9CF","#AF8EB8","#7B84B7","#3F6FAF",
                "#035AA6","#0473B4","#8491B7","#737AAB","#5E5E97")

p_infer_cnv_score <- 
  temp %>% 
  ggplot(aes(annotation,log(CNV_score))) + 
  geom_boxplot(aes(fill= annotation),color="NA",outlier.colour = "black") + 
  # ggsci::scale_fill_igv() + 
  scale_fill_manual(values = cols_brewer) + 
  theme_bw() -> p

ggsave(
  p_infer_cnv_score, 
  filename = file.path(infercnv_G5_out_dir, "cnv_boxplot.pdf"), 
  height = 8, 
  width  = 10
)







