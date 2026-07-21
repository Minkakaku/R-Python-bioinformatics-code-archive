infercnv_out_dir = "infercnv_group_cycling_result"
dir.create(infercnv_out_dir) 

all_data <- FindVariableFeatures(all_data)

Epi_anno_meta <-
  Epi[['RNA_snn_res.0.1']] %>% 
  rownames_to_column("cell_id") 


#B cells  pDCs  T cells Endothelial  Macrophage Fibroblasts        
#Plasma  MergeDC Epithelial Smooth muscle cells
all_data_anno_meta <- 
  all_data[[]] %>% 
  dplyr::filter(nCount_RNA > 100) %>%
  rownames_to_column('cell_id') %>%
  # 同时保留annotation和group.y列
  dplyr::select(cell_id, annotation, group.y) %>%
  dplyr::filter(annotation %in% c("Epithelial", "T cells")) %>%
  left_join(Epi_anno_meta, by = "cell_id") %>%
  map_dfc(as.character) %>%
  mutate(annotation2 = ifelse(
    annotation %in% "Epithelial",
    RNA_snn_res.0.1,
    annotation
  )) %>%
  # 选择最终需要的列：cell_id、重命名的annotation2和group.y
  dplyr::select(cell_id, annotation = annotation2, group.y)


# "N0_CA","N1_CA","Ni_CA","Nmi_CA"
# "N0_LN","Ni_LN","Nmi_LN","N1_LN"

# library(tidyverse)
# groups <- all_data@meta.data$group.y %>% unique()
# 
# for (group in groups) {
# # all_data <- all_data %>% AddMetaData(all_data_anno_meta)
#     meta_data_infer_subset <-
#       all_data_anno_meta %>%
#       filter(group.y == group) %>%
#       group_by(annotation) %>%
#       slice_sample(prop = 0.1)
# }

library(tidyverse)

groups <- all_data@meta.data$group.y %>% unique()

pb <- txtProgressBar(min = 0, max = length(groups), style = 3)  # [8](@ref)

for (i in seq_along(groups)) {
  group <- groups[i]
  
  #1. 动态创建当前group专属文件夹
  group_out_dir <- file.path("infercnv_group_cycling_result", group)
  if (!dir.exists(group_out_dir)) {
    dir.create(group_out_dir, recursive = TRUE, showWarnings = FALSE)  # [6,8](@ref)
  }
  
  #2. 筛选当前group的细胞子集
  meta_data_infer_subset <- all_data_anno_meta %>%
    filter(group.y == group) %>%
    group_by(annotation) %>%
    slice_sample(prop = 0.5)
  
  #3. 保存元数据到当前group文件夹
  meta_data_infer_subset %>%
    column_to_rownames("cell_id") %>%
    data.table::fwrite(
      file = file.path(group_out_dir, "meta_data_infer_subset.txt"),
      row.names = TRUE,
      col.names = FALSE,
      sep = "\t"
    )
  
  #4. 创建inferCNV对象（路径指向当前group文件夹）
  infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix = JoinLayers(all_data)@assays$RNA$counts[VariableFeatures(all_data)[1:500], meta_data_infer_subset$cell_id],
    annotations_file = file.path(group_out_dir, "meta_data_infer_subset.txt"),  # 关键修改
    gene_order_file = file.path('infercnv2', "gene_pos47.txt"),
    min_max_counts_per_cell = c(30, Inf),
    delim = "\t",
    ref_group_names = "T cells"
  )
  
  #5. 运行inferCNV分析
  infercnv_obj2 = infercnv::run(
    infercnv_obj,
    cutoff = 0.001,
    out_dir = group_out_dir,  # 输出到当前group文件夹
    cluster_by_groups = TRUE,
    plot_steps          =    FALSE,
    denoise             =    TRUE,
    HMM                 =    FALSE,
    no_prelim_plot      =    TRUE,
    png_res             =    300,
    # useRaster   = FALSE,
    output_format = 'pdf',
    no_plot = T #  figure margin too large ----
  )
  
  #6. 绘图并保存到当前group文件夹
  infercnv::plot_cnv(
    infercnv_obj2,
    out_dir = group_out_dir,  # 输出到当前group文件夹
    output_filename = 'infercnv',
    output_format = "pdf",
    png_res = 300,     
    cluster_by_groups = TRUE,
    x.range           =  "auto",
    x.center          =  1,
    title             =  "infercnv",
    color_safe_pal    =  FALSE
  )
  
  # 💾 7. 保存RDS对象到当前group文件夹
  saveRDS(infercnv_obj2, file.path(group_out_dir, 'infercnv_obj2.rds'))
  saveRDS(infercnv_obj, file.path(group_out_dir, 'infercnv_obj.rds'))
  
  
  # 1. 应用中值滤
  infercnv_obj2_default_median_filtering <- infercnv_obj2 %>%
    infercnv::apply_median_filtering()
  
  # 2. 绘制滤后CNV热图
  infercnv_obj2_default_median_filtering %>%
    infercnv::plot_cnv(
      out_dir = group_out_dir,
      output_filename = 'infercnv_median_filtering',
      output_format = "pdf",
      cluster_by_groups = TRUE,
      x.range           =  "auto",
      x.center          =  1,
      title             =  "infercnv",
      color_safe_pal    =  FALSE
    )
  
  # 3. 保存滤波后对象
  saveRDS(infercnv_obj2_default_median_filtering, 
          file.path(group_out_dir, 'infercnv_obj2_default_median_filtering.rds'))
  
  # 4. 计算CNV评分
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
  
  cnv_score <- cnvScore(infercnv_obj2@expr.data)
  colnames(cnv_score) <- 'CNV_score'
  
  # 5. 合并评分与元数据
  temp <- cnv_score %>% 
    rownames_to_column('cell_id') %>% 
    left_join(as.data.frame(meta_data_infer_subset), by = 'cell_id')
  
  # 6. 绘制CNV评分小提琴图
  p_violin <- ggplot(temp, aes(annotation, log(CNV_score))) +
    geom_violin(aes(fill = annotation), color = NA) +
    scale_fill_manual(values = cols_brewer) +
    theme_bw() +
    labs(title = paste("CNV Score Distribution -", group))
  
  ggsave(
    plot = p_violin,
    filename = file.path(group_out_dir, "cnv_vlnplot.pdf"),
    height = 6,
    width = 12
  )
  
  # 7. 绘制CNV评分箱线图
  p_boxplot <- ggplot(temp, aes(annotation, log(CNV_score))) +
    geom_boxplot(aes(fill = annotation), color = "black", outlier.colour = "red") +
    scale_fill_manual(values = cols_brewer) +
    theme_bw() +
    labs(title = paste("CNV Score Distribution -", group))
  
  ggsave(
    plot = p_boxplot,
    filename = file.path(group_out_dir, "cnv_boxplot.pdf"),
    height = 6,
    width = 12
  )
  
  # 更新进度条
  setTxtProgressBar(pb, i)  # [8](@ref)
}
close(pb)
