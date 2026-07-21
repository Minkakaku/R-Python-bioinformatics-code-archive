library(Seurat)
library(tidyverse)

all_data <- readRDS(file.path(out_dir, '/mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/Epithelial_subtype/all_data_Epithelial_annotation.rds'))

# inferCNV ---------------------------------------------------------------------
Epithelial <- all_data %>% JoinLayers() %>% .[, .$annotation1 %in% c("Enterocyte","Goblet","Paneth cell","Cancer stem cell")]
Epithelial <- FindVariableFeatures(Epithelial)




## set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 100 * 1024 ^ 3) # 100Gb

library(infercnv)
infercnv_out_dir = "Epithelial_infercnv" 
dir.create(infercnv_out_dir) 


## 映射-------------------------------
library(tibble)
library(data.table)

ref_cells <- c('B cell', 'Plasma', 'T cell') # , 'Endothelial', 'Fibroblast', 'Macrophage', 'Pericyte'
tumor_cells <- c("Enterocyte","Goblet","Paneth cell","Cancer stem cell")

meta_data_infer_ls <- 
  all_data[[c('annotation1', 'sample')]] %>% 
  dplyr::filter(annotation1 %in% c(ref_cells, tumor_cells)) %>%
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
      group_by(annotation1) %>%
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
      #raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[, meta_data_infer_subset$cell_id],
      raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[VariableFeatures(Epithelial)[1:2000], meta_data_infer_subset$cell_id],
      annotations_file  = file.path(infercnv_out_dir, sample, "meta_data_infer_subset.txt"),
      gene_order_file   = "/mnt/NexData_file/home/lxy/case/XE2603162231-02/20260616_XE2603162231-02/Epithelial_infercnv/mouse_gene_pos_mm10.txt",
      delim             = "\t",
      min_max_counts_per_cell = c(50,Inf),
      ref_group_names = ref_cells
    )
    infercnv_obj
  })



infercnv_obj_ls[[2]]@count.data %>% rowMeans %>% { sum(. > 0.01) } # 3280 <--- 
infercnv_obj_ls[[2]]@count.data %>% rowMeans %>% { sum(. < 0.01) } # 1477

#undebug(infercnv:::.single_tumor_leiden_subclustering)

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


## plot -------------------------------------------------------------------------(按照样本)
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
      ggplot(aes(annotation1, log(CNV_score)))+
      geom_violin(aes(fill= annotation1), color="NA", drop = FALSE, trim = FALSE)+
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
      height = 6,
      width  = 8
    )
  })


plot_ls2 <-
  temp_ls %>%
  map(function(temp){
    temp %>% 
      group_by(annotation1) %>%
      
      ggplot(aes(x = annotation1, y = log(CNV_score), fill = annotation1)) +
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
      ggplot(aes(annotation1, log(CNV_score)))+
      geom_violin(aes(fill= annotation1), color="NA", drop = FALSE, trim = FALSE)+
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
      group_by(annotation1) %>%
      
      ggplot(aes(x = annotation1, y = log(CNV_score), fill = annotation1)) +
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


# 把打分映射回原亚群画图
Epithelial <- readRDS("Epithelial_subtype/Epithelial_clean.rds")

##########下面之前跑过了，可以先不跑
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

##########上面可以先不跑,从下面开始

# 合并所有样本得分，cell_id 为行名
all_cnv_df <- map_dfr(cnv_score_ls, ~rownames_to_column(.x, "cell_id"))

# 提取上皮已有cell_id
epi_cell_ids <- colnames(Epithelial)

# 筛选仅上皮细胞的得分
epi_cnv_df <- all_cnv_df %>%
  filter(cell_id %in% epi_cell_ids) %>%
  column_to_rownames("cell_id")

# 追加元数据
Epithelial <- AddMetaData(Epithelial, metadata = epi_cnv_df, col.name = "CNV_score")

saveRDS(Epithelial, file.path(infercnv_out_dir, "Epithelial_clean_CNVscore.rds"))

library(ggsci)
library(ggpubr)

# 1 提取元数据 + 去除CNV_score NA
meta_df <- Epithelial[[]] %>%
  filter(!is.na(CNV_score)) %>%
  mutate(log_CNV = log(CNV_score))

# 2 整理长数据，实现四个横坐标共用绘图结构
# 构造分组标签列，方便分面
plot_df <- pivot_longer(
  data = meta_df,
  cols = c(sample, group, RNA_snn_res.0.8, subtype),
  names_to = "X_var_type",
  values_to = "X_value"
) %>%
  mutate(
    X_var_type = factor(
      X_var_type,
      levels = c("sample", "group", "RNA_snn_res.0.8", "subtype"),
      labels = c("Sample", "Group", "Cluster", "Subtype")
    )
  )


p_violin <- ggplot(plot_df, aes(x = X_value, y = log_CNV, fill = X_value)) +
  geom_violin(color = NA, drop = FALSE, trim = FALSE) +
  scale_fill_igv() +
  facet_wrap(~ X_var_type, scales = "free_x", ncol = 2) +
  theme_bw() +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none",
    plot.margin = margin(10,10,10,10)
  ) +
  labs(y = "log(CNV_score)", x = NULL)

# 保存
ggsave(
  filename = file.path(infercnv_out_dir, "Epithelial_CNVscore_violin.pdf"),
  plot = p_violin,
  width = 12,
  height = 12
)


p_box <- ggplot(plot_df, aes(x = X_value, y = log_CNV, fill = X_value)) +
  geom_boxplot(outlier.colour = NA, width = 0.45) +
  scale_fill_igv() +
  facet_wrap(~ X_var_type, scales = "free_x", ncol = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none",
    plot.margin = margin(10,10,10,10)
  ) +
  labs(y = "log(CNV_score)", x = NULL)

# 保存
ggsave(
  filename = file.path(infercnv_out_dir, "Epithelial_CNVscore_boxplot.pdf"),
  plot = p_box,
  width = 12,
  height = 12
)

