library(Seurat)
library(tidyverse)

out_dir <- 'luyichao'

all_data2 <- readRDS(file.path(out_dir, 'all_data2_TLS.rds'))

# inferCNV ---------------------------------------------------------------------
Epithelial <- all_data2[, all_data2$annotation == 'Epithelial']
Epithelial <- FindVariableFeatures(Epithelial)

## set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb

library(infercnv)
infercnv_out_dir = "luyichao/inferCNV_episub_allgenes" 
dir.create(infercnv_out_dir) 


## 映射-------------------------------
library(tibble)
library(data.table)

ref_cells <- c('B cell', 'Plasma', 'T cell') # , 'Endothelial', 'Fibroblast', 'Macrophage', 'Pericyte'
tumor_cells <- "Epithelial"

meta_data_infer_ls <- 
  all_data2[[c('annotation', 'sample')]] %>% 
  dplyr::filter(annotation %in% c(ref_cells, tumor_cells)) %>%
  rownames_to_column('cell_id') %>% 
  group_by(sample) %>%
  nest() %>% 
  deframe()


# cell is too much
# less than 50000
set.seed(1234)
meta_data_infer_subset_ls <-
  meta_data_infer_ls %>% 
  imap(function(meta_data_infer, sample){
    meta_data_infer2 <- 
      meta_data_infer %>% 
      group_by(annotation) %>%
      slice_sample(prop = 1)
  })


# export as file
meta_data_infer_subset_ls %>% 
  imap(function(meta_data_infer_subset, sample){
    if(!file.exists(file.path(infercnv_out_dir, sample))){
      dir.create(file.path(infercnv_out_dir, sample))
    }
    data.table::fwrite(
      meta_data_infer_subset, 
      file = file.path(infercnv_out_dir, sample, "meta_data_infer_subset.txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
  })



infercnv_obj_ls <- 
  meta_data_infer_subset_ls %>% 
  imap(function(meta_data_infer_subset, sample){ # browser()
    infercnv_obj <- CreateInfercnvObject(
      raw_counts_matrix = JoinLayers(all_data2)@assays$RNA$counts[, meta_data_infer_subset$cell_id],
      # raw_counts_matrix = JoinLayers(all_data2)@assays$RNA$counts[VariableFeatures(all_data2)[1:2000], meta_data_infer_subset$cell_id],
      annotations_file  = file.path(infercnv_out_dir, sample, "meta_data_infer_subset.txt"),
      gene_order_file   = "luyichao/gene_pos47.txt",
      delim             = "\t",
      min_max_counts_per_cell = c(50,Inf),
      ref_group_names = ref_cells
    )
    infercnv_obj
  })



infercnv_obj_ls[[2]]@count.data %>% rowMeans %>% { sum(. > 0.01) } # 3280 <--- 
infercnv_obj_ls[[2]]@count.data %>% rowMeans %>% { sum(. < 0.01) } # 1477

undebug(infercnv:::.single_tumor_leiden_subclustering)

options(scipen = 100)
infercnv_obj_ls2 <- 
  infercnv_obj_ls %>% 
  imap(function(infercnv_obj, sample){
    infercnv_obj2 <- infercnv::run(
      infercnv_obj,
      cutoff              =    0.01,
      out_dir    =   file.path(infercnv_out_dir, sample),
      cluster_by_groups   =    TRUE, 
      # k_nn                =    20,
      # leiden_method       =    "simple", # PCA 
      plot_steps          =    FALSE,
      denoise             =    TRUE,
      HMM                 =    FALSE,
      no_prelim_plot      =    TRUE,
      png_res             =    300,
      # useRaster   = FALSE,
      output_format = 'pdf',
      no_plot = T #  figure margin too large 
    )
    infercnv_obj2
  })

infercnv_obj_ls2 %>% 
  imap(function(infercnv_obj2, sample){
    infercnv::plot_cnv(
      infercnv_obj2,
      out_dir           =  file.path(infercnv_out_dir, sample),
      output_filename   =  'infercnv',
      output_format     =  "pdf",
      png_res = 300,
      cluster_by_groups =  TRUE,
      x.range           =  "auto",
      x.center          =  mean(infercnv_obj2@expr.data), # 1,
      title             =  "infercnv",
      color_safe_pal    =  FALSE
    )
  })



infercnv_obj2_default_median_filtering_ls <-
  infercnv_obj_ls2 %>%
  map(function(infercnv_obj2){
    infercnv::apply_median_filtering(infercnv_obj2)
  })


infercnv_obj2_default_median_filtering_ls %>% 
  imap(function(infercnv_obj2_default_median_filtering, sample){
    infercnv::plot_cnv(
      infercnv_obj2_default_median_filtering, 
      out_dir           =  file.path(infercnv_out_dir, sample),
      output_filename   =  'infercnv_median_filtering',
      output_format="pdf",
      cluster_by_groups =  TRUE,
      x.range           =  "auto",
      x.center          =  1,
      title             =  "infercnv",
      color_safe_pal    =  FALSE
    )
  })

saveRDS(infercnv_obj2_default_median_filtering_ls, file.path(infercnv_out_dir, 'infercnv_obj2_default_median_filtering_ls.rds'))

saveRDS(infercnv_obj_ls2, file.path(infercnv_out_dir, 'infercnv_obj_ls2.rds'))
saveRDS(infercnv_obj_ls, file.path(infercnv_out_dir, 'infercnv_obj_ls.rds'))


## plot -------------------------------------------------------------------------
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

cnv_score_ls <- 
  infercnv_obj_ls2 %>% 
  map(function(infercnv_obj2){
    cnv_score <- cnvScore(infercnv_obj2@expr.data)
    colnames(cnv_score)<-'CNV_score'
    cnv_score
  })


temp_ls <- 
  cnv_score_ls %>%
  imap(function(cnv_score, sample){
    cnv_score %>% 
      rownames_to_column('cell_id') %>% 
      left_join(.,as.data.frame(meta_data_infer_subset_ls[[sample]]),'cell_id') -> temp
    temp
  })


plot_ls <- 
  temp_ls %>% 
  map(function(temp){
    temp %>% 
      ggplot(aes(annotation, log(CNV_score)))+
      geom_violin(aes(fill= annotation), color="NA", drop = FALSE, trim = FALSE)+
      ggsci::scale_fill_igv()+
      # scale_fill_manual(values = color_v)+
      theme_bw()  +
      RotatedAxis()
  })

plot_ls %>% 
  iwalk(function(plot, sample){
    ggsave(
      plot,
      filename = file.path(infercnv_out_dir, sample, "CNV_score.pdf"),
      height = 8,
      width  = 15
    )
  })


plot_ls2 <-
  temp_ls %>%
  map(function(temp){
    temp %>% 
      group_by(annotation) %>%
      
      ggplot(aes(x = annotation, y = log(CNV_score), fill = annotation)) +
      geom_boxplot(outlier.colour = NA) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )  
  })


plot_ls2 %>% 
  iwalk(function(plot2, sample){
    ggsave(
      plot2,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_barplot.pdf"),
      height = 6,
      width  = 8
    )
  })


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

cnv_score_noQuantile_ls <- 
  infercnv_obj_ls2 %>% 
  map(function(infercnv_obj2){
    cnv_score <- cnvScore(infercnv_obj2@expr.data)
    colnames(cnv_score)<-'CNV_score'
    cnv_score
  })


temp_noQuantile_ls <- 
  cnv_score_noQuantile_ls %>%
  imap(function(cnv_score, sample){
    cnv_score %>% 
      rownames_to_column('cell_id') %>% 
      left_join(., as.data.frame(meta_data_infer_subset_ls[[sample]]), 'cell_id') -> temp
    temp
  })


plot_noQuantile_ls <- 
  temp_noQuantile_ls %>% 
  map(function(temp){
    temp %>% 
      ggplot(aes(annotation, log(CNV_score)))+
      geom_violin(aes(fill= annotation), color="NA", drop = FALSE, trim = FALSE)+
      ggsci::scale_fill_igv()+
      # scale_fill_manual(values = color_v)+
      theme_bw()  +
      RotatedAxis()
  })

plot_noQuantile_ls %>% 
  iwalk(function(plot, sample){
    ggsave(
      plot,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_noQuantile.pdf"),
      height = 8,
      width  = 15
    )
  })


plot_noQuantile_ls2 <-
  temp_noQuantile_ls %>%
  map(function(temp){
    temp %>% 
      group_by(annotation) %>%
      
      ggplot(aes(x = annotation, y = log(CNV_score), fill = annotation)) +
      geom_boxplot(outlier.colour = NA) +
      ggsci::scale_fill_igv() +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      )  
  })


plot_noQuantile_ls2 %>% 
  iwalk(function(plot2, sample){
    ggsave(
      plot2,
      filename = file.path(infercnv_out_dir, sample, "CNV_score_barplot_noQuantile.pdf"),
      height = 6,
      width  = 8
    )
  })
