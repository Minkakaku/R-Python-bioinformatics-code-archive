# installed into linux server
# mamba install conda-forge::jags
# install.packages('rjags')
# BiocManager::install("infercnv")


# BiocManager::install("infercnv")
# https://sourceforge.net/projects/mcmc-jags/files
# If you are using R 4.1.3 or earlier then install JAGS-4.3.0.exe


library(Seurat)
library(tidyverse)
library(infercnv)


# inferCNV ---------------
# used to infer epethial cells
# immunes cell used as control

# 当前路径为pipeline文件夹

infercnv_in_dir  <- "infercnv/infercnv_files"
infercnv_out_dir <- "infercnv/infercnv_run"
dir.create(infercnv_in_dir) 
dir.create(infercnv_out_dir) 

# all_data_used <- merge(
#   sc_epi2, 
#   list(all_data[, all_data$annotation_full %in% c("Plasma", "Fibroblasts", "Myeloid cells")])
# )
# all_data_used$annotation_used <- ifelse(
#   all_data_used$annotation_full == "Epithelial",
#   all_data_used$annotation_epi,
#   all_data_used$annotation_full
# )

all_data_used <- readRDS("all_data_full2.rds")
all_data_used <- all_data_used[, !is.na(all_data_used$group_all)]

table( all_data_used$project, all_data_used$nCount_Spatial.008um > 50 )
#           FALSE   TRUE
# A210239   406399  98193
# A223925   255793 429956
# A224091   364135 291985
# A234244   161401  68013
# A240495   169799  86818
# A244628   130946  47125
# A246283   122615 191822
# A247950   100839  99406
# A251034   136856 205799
# A251692   119881 103563
# A255670   148606  99283
# A266027    87933 107024
# A271888    84703 126751
# A306965   211017 137740
# F550277   162622 162642
# LC_296844  85849 553822
# LC_297623 266593 289185
# LC_309016 131791 169817

# only used NORMAL and CANCER, ADVANCE_CANCER  
# QC: 50

DefaultAssay(all_data_used) <- "Spatial.008um"
all_data_used[['sketch']]   <- NULL


infercnv_used_cells <-
  c(
    # 'Chief_cell',
    'Endocrine',
    # 'Enterocytes',
    # 'Enteroendocrine',
    'GMC',
    # 'Goblet_cell',
    'PMC',
    'Tumor',
    'parietal_cells'
    # 'tuft_cells'
  )


all_cell_id_used_ls <-
  all_data_used$project %>%
  unique %>%
  setNames(., .) %>%
  map(function(x){
    idx  <- all_data_used$project == x & all_data_used$nCount_Spatial.008um > 50 & all_data_used$celltype %in% infercnv_used_cells
    meta <- all_data_used[[]] %>% .[idx, ] 
    
    norm_exists <- any( c('NORMAL') %in% meta$group_all ) 
    tumor_exists  <- any( c('CANCER', 'ADVANCE_CANCER') %in% meta$group_all )  
    
    if(!norm_exists || !norm_exists){
      return(NULL)
    }

    res <- colnames(all_data_used)[idx]
    res 
  })

all_cell_id_used_ls <- all_cell_id_used_ls %>% purrr::discard(is.null)


# infercnv obj
all_infercnv_obj_ls <-
  all_cell_id_used_ls %>%
  imap(function(cell_id, name){
    
    dir.create(file.path(infercnv_in_dir, name))
    
    all_data <- all_data_used[, cell_id] %>% JoinLayers()
    
    #  export meta
    all_data[['group_all']] %>% 
      data.table::fwrite(
        file = file.path(infercnv_in_dir, name, "meta_data_infer_subset.txt"),
        row.names = TRUE,
        col.names = FALSE,
        sep = "\t"
      )
    
    
    raw_counts_matrix_infercnv <- LayerData(all_data, "counts")
    
    infercnv_obj = infercnv::CreateInfercnvObject(
      raw_counts_matrix = raw_counts_matrix_infercnv,
      annotations_file  = file.path(infercnv_in_dir, name, "meta_data_infer_subset.txt"),
      gene_order_file   = file.path('infercnv', "GRCh38_gene_id_pos.txt"),
      
      min_max_counts_per_cell = c(50, Inf), # <--- qc prefilter ------
      delim             = "\t",
      ref_group_names   = c("NORMAL")
    )
    
    infercnv_obj
  })

saveRDS(all_infercnv_obj_ls, file.path('infercnv/', 'all_infercnv_obj_ls.rds') )



# raw_counts_matrix_infercnv %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

raw_counts_matrix_infercnv[, meta_data_infer_subset$cell_id] %>% rowMeans %>% { sum(. > 0.01) }   
raw_counts_matrix_infercnv[, meta_data_infer_subset$cell_id] %>% rowMeans %>% { sum(. > 0.05) }  

infercnv_obj@count.data %>% rowMeans %>% { sum(. > 0.01) }   
infercnv_obj@count.data %>% rowMeans %>% { sum(. > 0.05) }   

# options()$bitmapType
# [1] "cairo"
# options(bitmapType="Xlib")
# options(bitmapType="cairo")

options(scipen = 100)
# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
infercnv_obj_default = infercnv::run(
  infercnv_obj,
  cutoff              =    0.05, 
  out_dir             =    infercnv_out_dir,
  cluster_by_groups   =    TRUE, 
  plot_steps          =    FALSE,
  denoise             =    TRUE,
  HMM                 =    FALSE,
  no_prelim_plot      =    TRUE,
  write_expr_matrix   =    TRUE,
  png_res             =    300,
  no_plot = T #  figure margin too large  
)

# download inferCNV object from linux server
# plot inferCNV locally
saveRDS(infercnv_obj_default, file.path(infercnv_out_dir, 'infercnv_obj_default.rds'))

infercnv_obj_default <- readRDS("inferCNV/infercnv_obj_default.rds")

infercnv::plot_cnv(
  infercnv_obj_default,
  out_dir           =  infercnv_out_dir, 
  output_filename   =  'infercnv',
  output_format     =  "pdf",
  png_res = 300,
  cluster_by_groups =  TRUE,
  x.range           =  "auto",
  x.center          =  1,
  title             =  "infercnv",
  color_safe_pal    =  FALSE
)

infercnv_obj_default_median_filtering <-
  infercnv_obj_default %>%
  infercnv::apply_median_filtering()

infercnv_obj_default_median_filtering %>%
  infercnv::plot_cnv(
    out_dir           =  infercnv_out_dir, # 'inferCNV', 
    output_filename   =  'infercnv_median_filtering',
    output_format     =  "pdf",
    png_res = 300,
    cluster_by_groups =  TRUE,
    x.range           =  "auto",
    x.center          =  infercnv_obj_default_median_filtering@expr.data %>% median,
    title             =  "infercnv",
    color_safe_pal    =  FALSE
  )
saveRDS(infercnv_obj_default_median_filtering, file.path(infercnv_out_dir, 'infercnv_obj_default_median_filtering.rds'))


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

cnvScore <- function(data){
  sd_num <- apply(data, 1, sd)
  data <- data[sd_num > median(sd_num), ]
  
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}


cnvScore <- function(data){
  # sd_num <- apply(data, 1, sd)
  # data <- data[sd_num > median(sd_num), ]
  
  # rng <- quantile(data, c(0.01, 0.99))
  # data <- pmax( pmin(data, rng[2]), rng[1] )
  
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}

cnv_score <- cnvScore(infercnv_obj_default@expr.data)
# cnv_score <- cnvScore(infercnv_obj_default_median_filtering@expr.data)
colnames(cnv_score) <-'CNV_score'
cnv_score_tidy <-
  cnv_score %>% 
  rownames_to_column('cell_id') %>% 
  left_join(.,as.data.frame(meta_data_infer_subset),'cell_id')  

p_infer_cnv_score <-
  cnv_score_tidy %>% 
  ggplot(aes(annotation_full, log(CNV_score))) + geom_violin(aes(fill= annotation_full),color="NA",drop = FALSE)+
  ggsci::scale_fill_igv()+
  # scale_fill_manual(values = color_v)+
  theme_bw() +
  RotatedAxis()

ggsave(
  p_infer_cnv_score,
  filename = file.path(infercnv_out_dir, "cnv_vlnplot.pdf"),
  height = 6,
  width  = 12
)



# spatial plot --------------
# Idents(all_data) <- 'Spatial.008um'
# Idents(all_data) <- 'Spatial.qc.pass'

sc_epi_infer_anno <-
  sc_epi[['Spatial.008um_snn_res.0.4']] %>% 
  rownames_to_column("cell_id") %>%
  mutate(annotation_epi = ifelse(
    Spatial.008um_snn_res.0.4 %in% c('3', '6', '12' ,'13'),
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











