# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 150 * 1024 ^ 3) # 10Gb


library(infercnv)
library(Seurat)
library(dplyr)
library(tidyverse)
library(stringr)

all_data <- readRDS(file.path("data_for_analysis", "all_data_projected_rctd_20251230.rds"))

tumor_out_dir <- "tumor_subcluster"
Tumor <- readRDS(file.path(tumor_out_dir, "Tumor.rds"))
table(Tumor$Spatial.Polygons_snn_res.0.6)


### 筛选细胞做inferCNV -------
tumor_anno_meta <-
  Tumor[['annotation_sub']] %>% 
  rownames_to_column("cell_id") 

all_data$annotation <- all_data$first_type

all_data_anno_meta <- 
  all_data[[]] %>% dplyr::filter(nCount_Spatial.Polygons > 100) %>% .['annotation'] %>%
  # all_data[['annotation']] %>%
  rownames_to_column('cell_id') %>%
  dplyr::filter(annotation %in% c("Tumor", "Macrophages", "NK_T-cells", "Monocytes", "B-cells")) %>%
  left_join(tumor_anno_meta, by = "cell_id") %>% 
  map_dfc(as.character) %>%
  mutate(annotation2 = ifelse(
    annotation %in% "Tumor",
    annotation_sub,
    annotation
  )) %>% 
  dplyr::select(cell_id, annotation = annotation2)




# all_data <- all_data %>% AddMetaData(all_data_anno_meta)  
meta_data_infer_subset <- 
  all_data_anno_meta %>% 
  group_by(annotation) %>% 
  slice_sample(prop = 0.25)

### inferCNV ---------------
# used to infer epethial cells
# immunes cell used as control

# 当前路径为pipeline文件夹
library(infercnv)
infercnv_out_dir = "infercnv_tumor_sub_variable_all500"
dir.create(infercnv_out_dir) 


# meta_data_infer <-
#   all_data[[]] %>%
#   rownames_to_column("cell_id") %>%
#   dplyr::filter(
#     annotation %in% c("Club cells","TRB-SC","Basal","AT2",'Endothelial',"Plasma","Macrophage","T cells")
#   ) %>%
#   dplyr::select(cell_id, annotation, Spatial.008um_snn_res.0.8) %>%
#   map_dfc(as.character) 

# unique(meta_data_infer$annotation)
# "Acinar cells"     "Ductal cells"     "Endocrine cells"  "Epithelial cells" "NKT cells"


# cell is too much
# less than 50000
# 取1/20细胞： 1179497 --> 47179.88
# 1662182 --> 83000
# meta_data_infer_subset <-
#   meta_data_infer %>%
#   group_by(annotation) %>%
#   slice_sample(prop = 0.5)

# 合并annotation
# meta_data_infer_subset <-
#   meta_data_infer_subset %>%
#   mutate(Spatial.008um_snn_res.0.8 = str_c('Cluster', Spatial.008um_snn_res.0.8, sep = "_")) %>%
#   mutate(annotation2 = ifelse(
#     annotation %in% c("Plasma", "Macrophage","T cells","Endothelial"),
#     annotation,
#     Spatial.008um_snn_res.0.8
#   )) %>%
#   dplyr::select(
#     cell_id,
#     annotation = annotation2
#   )

# export as file
meta_data_infer_subset %>%
  column_to_rownames("cell_id") %>%
  data.table::fwrite(
    file = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
    row.names = TRUE,
    col.names = FALSE,
    sep = "\t"
  )

# undebug(CreateInfercnvObject)
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix = JoinLayers(all_data)@assays$Spatial.Polygons$counts[VariableFeatures(all_data)[1:500], meta_data_infer_subset$cell_id] ,
  annotations_file  = file.path(infercnv_out_dir, "meta_data_infer_subset.txt"),
  gene_order_file   = file.path('infercnv2', "hsa_gene_position.txt"),
  
  min_max_counts_per_cell = c(30, Inf),
  
  delim             = "\t",
  ref_group_names   = c("Macrophages", "NK_T-cells", "Monocytes", "B-cells")
)


# all_data@assays$RNA@counts %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

infercnv_obj@count.data %>% rowMeans %>% { sum(. > 0.001) } # 425 <--- 
infercnv_obj@count.data %>% rowMeans %>% { sum(. < 0.005) } # 137


infercnv_obj2  = infercnv::run(
  infercnv_obj,
  cutoff              =    0.001,
  out_dir             =    infercnv_out_dir,
  cluster_by_groups   =    TRUE,
  plot_steps          =    FALSE,
  denoise             =    TRUE,
  HMM                 =    FALSE,
  no_prelim_plot      =    TRUE,
  png_res             =    300,
  # useRaster         =    FALSE,
  output_format       =    'pdf',
  no_plot             =    T #  figure margin too large ----
)
infercnv::plot_cnv(
  infercnv_obj2,
  out_dir           =  infercnv_out_dir, 
  output_filename   =  'infercnv',
  output_format     =  "pdf",
  png_res           =  300, 
  cluster_by_groups =  TRUE,
  x.range           =  "auto",
  x.center          =  1,
  title             =  "infercnv",
  color_safe_pal    =  FALSE
)

saveRDS(infercnv_obj2, file.path(infercnv_out_dir, 'infercnv_obj2.rds'))
saveRDS(infercnv_obj, file.path(infercnv_out_dir, 'infercnv_obj.rds'))

# # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
# infercnv_obj_default = infercnv::run(
#   infercnv_obj,
#   cutoff              =    0.01, 
#   out_dir    =    infercnv_out_dir,
#   cluster_by_groups   =    TRUE, 
#   plot_steps          =    FALSE,
#   denoise             =    TRUE,
#   HMM                 =    FALSE,
#   no_prelim_plot      =    TRUE,
#   png_res             =    300,
#   # useRaster   = FALSE,
#   output_format = 'pdf',
#   no_plot = T #  figure margin too large ----
# )
# 
# # C stack error
# infercnv::plot_cnv(
#   infercnv_obj_default,
#   out_dir           =  infercnv_out_dir,
#   output_filename   =  'infercnv',
#   output_format     =  "pdf",
#   png_res = 60,
#   cluster_by_groups =  TRUE,
#   x.range           =  "auto",
#   x.center          =  1,
#   title             =  "infercnv",
#   color_safe_pal    =  FALSE
# )


# download inferCNV object from linux server
# plot inferCNV locally
# saveRDS(infercnv_obj_default, file.path(infercnv_out_dir, 'infercnv_obj_default.rds'))
# 
# infercnv_obj_default <- readRDS("inferCNV/infercnv_obj_default.rds")
# 
# infercnv::plot_cnv(
#   infercnv_obj_default,
#   out_dir           =  'inferCNV', 
#   output_filename   =  'infercnv',
#   output_format     =  "pdf",
#   png_res = 300, 
#   cluster_by_groups =  TRUE,
#   x.range           =  "auto",
#   x.center          =  1,
#   title             =  "infercnv",
#   color_safe_pal    =  FALSE
# )
# 
infercnv_obj2_default_median_filtering <-
  infercnv_obj2 %>%
  infercnv::apply_median_filtering()

infercnv_obj2_default_median_filtering %>%
  infercnv::plot_cnv(
    out_dir           =  infercnv_out_dir,
    output_filename   =  'infercnv_median_filtering',
    output_format     =  "pdf",
    cluster_by_groups =  TRUE,
    x.range           =  "auto",
    x.center          =  1,
    title             =  "infercnv",
    color_safe_pal    =  FALSE
  )
saveRDS(infercnv_obj2_default_median_filtering, file.path(infercnv_out_dir, 'infercnv_obj2_default_median_filtering.rds'))

# expFile='expFile.txt'
# groupFiles='groupinfo.txt'
# geneFile='gene_pos.txt'
# infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
#                                     annotations_file=groupFiles,
#                                     delim="\t",
#                                     gene_order_file= geneFile,
#                                     ref_group_names=c('B_cells','plasma_cells','T_cells','microglial')) 
# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics



# infercnv::plot_cnv(
#   infercnv_obj2,
#   out_dir           =  'infercnv_output_2',
#   output_filename   =  'infercnv',
#   output_format     =  "pdf",
#   png_res = 300,
#   cluster_by_groups =  TRUE,
#   x.range           =  "auto",
#   x.center          =  1,
#   title             =  "infercnv",
#   color_safe_pal    =  FALSE
# )
# infercnv_obj2 = infercnv::run(infercnv_obj,
#                               cutoff=0.01, 
#                               out_dir= "infercnv_output_2",  
#                               cluster_by_groups=T ,  
#                               output_format="pdf",
#                               hclust_method="ward.D2", plot_steps=F,no_plot = T)
# infercnv::plot_cnv(infercnv_obj2, 
#                    out_dir='infercnv_output_2',
#                    output_filename='infercnv', 
#                    x.range="auto",
#                    x.center=1,
#                    title = "infercnv", 
#                    output_format="pdf",
#                    color_safe_pal = FALSE)
# infercnv_obj2 = infercnv::run(infercnv_obj,
#                               cutoff=0.1, 
#                               out_dir= "infercnv_output_sd",  
#                               write_expr_matrix=T,
#                               cluster_by_groups=T ,  
#                               output_format="pdf",
#                               sd_amplifier=3,
#                               denoise=T,
#                               noise_logistic=TRUE,
#                               hclust_method="ward.D2", plot_steps=F)
# saveRDS(infercnv_obj2, file.path(infercnv_output, 'infercnv_obj2.rds'))
# saveRDS(infercnv_obj, file.path(infercnv_output, 'infercnv_obj.rds'))


# cnv score ------------
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
# }
cnv_score <- cnvScore(infercnv_obj2@expr.data)
colnames(cnv_score) <-'CNV_score'
temp <- 
  cnv_score %>% 
  rownames_to_column('cell_id') %>% 
  left_join(.,as.data.frame(meta_data_infer_subset),'cell_id')  

p_infer_cnv_score <-
  temp %>% 
  ggplot(aes(annotation,log(CNV_score))) + geom_violin(aes(fill = annotation), color = "NA", drop = FALSE) + 
  ggsci::scale_fill_igv() + 
  # scale_fill_manual(values = color_v)+
  theme_bw() -> p

ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_vlnplot.pdf"),
  height   = 5,
  width    = 15
)

cols_brewer = c("#f94144","#f25c27","#eb770a","#f89f12","#f9c74f","#90be6d",
                "#43aa8b","#4d908e","#577590","#277da1","#967aa1",
                "#F2698D","#F391AE","#F3B9CF","#AF8EB8","#7B84B7","#3F6FAF",
                "#035AA6","#0473B4","#8491B7","#737AAB","#5E5E97")

p_infer_cnv_score <-
  temp %>% 
  ggplot(aes(annotation,log(CNV_score)))+
  geom_boxplot(aes(fill= annotation),color="NA",outlier.colour = "black")+
  # ggsci::scale_fill_igv()+
  scale_fill_manual(values = cols_brewer)+
  theme_bw() -> p

ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_boxplot.pdf"),
  height = 8,
  width  = 12
)



# spatial plot --------------
# Idents(all_data) <- 'Spatial.008um'
# Idents(all_data) <- 'Spatial.qc.pass'

sc_epi_infer_anno <-
  sc_epi[['Spatial.008um_snn_res.0.4']] %>% 
  rownames_to_column("cell_id") %>%
  mutate(annotation_epi = ifelse(
    Spatial.008um_snn_res.0.4 %in% c('0', '1', '3' ,'4', '5', '6', '7', '8', '9'),
    "Tumor cells",
    "Epithelial cells"
    # Spatial.008um_snn_res.0.4
  )) 

all_data_infer_anno <-
  all_data[['annotation']] %>%
  rownames_to_column("cell_id") %>%
  left_join(sc_epi_infer_anno, by = "cell_id") %>%
  mutate(annotation2 = ifelse(
    is.na(annotation_epi),
    as.character(annotation),
    as.character(annotation_epi)
  )) %>%
  column_to_rownames("cell_id")


sp_p_ls <-
  Images(all_data) %>%
  set_names(., .) %>%
  map(function(img){
    all_data %>%
      AddMetaData(all_data_infer_anno['annotation2']) %>%
      SpatialDimPlot(
        images = img, 
        group.by = "annotation2",
        pt.size.factor	= 2,
        image.scale = "lowres",
        crop = FALSE
      ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img)
  })


# export
pdf(
  file.path(epi_out_dir, "spatial plot of infercnv analysis.pdf"),
  height = 8,
  width  = 8
)
sp_p_ls %>% walk(print)
dev.off()













