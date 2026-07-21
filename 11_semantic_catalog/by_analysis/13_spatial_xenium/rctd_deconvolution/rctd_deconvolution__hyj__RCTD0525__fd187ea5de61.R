#只对6-13区域所在样本RCTD-----

library(tidyverse)
library(Seurat)
getwd()
library(patchwork)
options(java.parameters = '-Xmx10g')  # options('java.parameters')
library(rJava)
library(RBioFormats)
library(ggpubr)
library(SeuratDisk)
library(ggrepel)
library(loupeR)

library(spacexr)

rctd_outdir <- "RCTD0526/"
dir.create(rctd_outdir)


#1.单细胞数据-----
scdata <- readRDS("/mnt/NexData_file/home/sh/Projects/AE2507033447-01/16vs16/01.Mission/Main_mapping.rds")

adata_sc <- sc$read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2507033447-01/13vs13/Main_mapping.h5ad")
adata_sc_meta <- adata_sc$obs

adata_sc_meta$sub_cell_type %>% table()
adata_sc_meta$cell_type %>% table()
adata_sc_meta$sub_cell_type <- gsub("T/NK","T_NK",adata_sc_meta$sub_cell_type)

scdata@meta.data <- adata_sc_meta

scdata1 <- scdata[,!scdata$sub_cell_type%in%c("Epithelial cells")]

anno_counts <- table(scdata1$sub_cell_type)
print(anno_counts)

Idents(scdata1) <- "sub_cell_type"
counts_matrix <- GetAssayData(scdata1, "RNA","counts")
cluster_sc <- as.factor(scdata1$sub_cell_type)
nUMI_sc    <- scdata1$total_counts

#sc_cell_filter_idx <-  cluster_sc != "Unknown"

# 2.create the RCTD reference object------
reference <- Reference(
  counts_matrix, 
  droplevels(cluster_sc), 
  nUMI_sc
)

saveRDS(reference,file = file.path(rctd_outdir, "reference.rds"))

#xenium数据
xenium_ROI <- readRDS("result/xenium_ROI_FOV.rds")
table(xenium_ROI$sample,xenium_ROI$sample_raw)

#XENIUM6-13在N16326-G1-1样本中--------

all_data <- readRDS("R_OUTPUT/all_data.rds")
all_data$sample %>% table()

#N16326_G1_1


# 3. create the RCTD query object-----
all_data_sp_rctd_ls <- all_data %>% 
  JoinLayers() %>% 
  SplitObject('sample')

cells_position_ls <- readRDS("~/R_OUTPUT/data_for_analysis/cells_position_ls.rds")

all_data_sp_rctd <- all_data_sp_rctd_ls["N16326_G1_1"]

rctd_query <-
  all_data_sp_rctd %>%
  imap(function(xenium, name){ # browser()
    counts_hd <- xenium %>% LayerData(layer = "counts")
    coords    <- 
      cells_position_ls[[name]] %>% 
      dplyr::select(
        cell_id,
        x = x_centroid, 
        y = y_centroid
      ) %>%
      column_to_rownames('cell_id')
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })

saveRDS(rctd_query,file = file.path(rctd_outdir, "rctd_query.rds"))

# 
# rctd_query_ls <-
#   all_data_sp_rctd_ls %>%
#   imap(function(xenium, name){ # browser()
#     counts_hd <- xenium %>% LayerData(layer = "counts")
#     coords    <- 
#       cells_position_ls[[name]] %>% 
#       dplyr::select(
#         cell_id,
#         x = x_centroid, 
#         y = y_centroid
#       ) %>%
#       column_to_rownames('cell_id')
#     
#     # create the RCTD query object
#     query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
#     query
#   })
# 
# saveRDS(rctd_query_ls,file = file.path(rctd_outdir, "rctd_query_ls.rds"))


#4. run RCTD-----
rctd_result <- 
  rctd_query %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 30, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

saveRDS(rctd_result, file = file.path(rctd_outdir, "rctd_result.rds"))

# rctd_result_ls <- 
#   rctd_query_ls %>%
#   imap(function(query, name){
#     # run RCTD
#     RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 30, UMI_min_sigma = 50)
#     RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
#     
#     RCTD
#   })
# 
# saveRDS(rctd_result_ls, file = file.path(rctd_outdir, "rctd_result_ls.rds"))

#rctd_result_ls <- readRDS("RCTD/rctd_result_ls.rds")

#5. addMetaData-----
all_data_sp_rctd_anno <-
  all_data_sp_rctd %>%
  imap(function(xenium, name){ # browser()
    rctd <- rctd_result[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)
    
    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })

saveRDS(all_data_sp_rctd_anno,file = file.path(rctd_outdir, "rctd.rds"))

# 
# all_data_sp_rctd_anno_ls <-
#   all_data_sp_rctd_ls %>%
#   imap(function(xenium, name){ # browser()
#     rctd <- rctd_result_ls[[name]]
#     xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)
#     
#     xenium$first_type <- as.character(xenium$first_type)
#     xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
#     xenium
#   })
# 
# saveRDS(all_data_sp_rctd_anno_ls,file = file.path(rctd_outdir, "rctd.rds"))

#all_data_sp_rctd_anno_ls <- readRDS("rctd.rds")


#6. dimplot 
# all_data_sp_anno <-
#   merge(
#     all_data_sp_rctd_anno_ls[[1]],
#     all_data_sp_rctd_anno_ls[-1]
#   )

meta_rctd_result <-
  rctd_result %>%
  map(~ .@results$results_df) %>%
  bind_rows()

all_data_rctd_N16326_G1_1 <-
  all_data[,all_data$sample %in% c("N16326_G1_1")] %>%
  AddMetaData(meta_rctd_result)

all_data_rctd_N16326_G1_1$first_type %>% table()

saveRDS(all_data_rctd_N16326_G1_1,file = file.path(rctd_outdir, "all_data_rctd_N16326_G1_1.rds"))


# meta_rctd_result <-
#   rctd_result_ls %>%
#   map(~ .@results$results_df) %>%
#   bind_rows()
# 
# all_data_rctd <-
#   all_data %>%
#   AddMetaData(meta_rctd_result)
# 
# saveRDS(all_data_rctd,file = file.path(rctd_outdir, "all_data_rctd.rds"))


#all_data_rctd_N16326_G1_1 <- readRDS(file = file.path(rctd_outdir, "all_data_rctd_N16326_G1_1.rds"))
all_data_rctd_N16326_G1_1 <- JoinLayers(all_data_rctd_N16326_G1_1)
all_data_rctd_N16326_G1_1$cell_id %>% head()
all_data_rctd_N16326_G1_1$cell_id_raw <- all_data_rctd_N16326_G1_1$cell_id
all_data_rctd_N16326_G1_1$cell_id <- paste0(all_data_rctd_N16326_G1_1$project,"_",all_data_rctd_N16326_G1_1$cell_id_raw)

##加入images的xenium数据----
xenium_ROI <- readRDS("result/xenium_ROI_FOV.rds")

xenium_ROI_RCTD_N16326_G1_1 <- xenium_ROI
xenium_ROI_RCTD_N16326_G1_1@meta.data %>% head()
xenium_ROI_RCTD_N16326_G1_1@meta.data$cell_id_raw %>% head()
xenium_ROI_RCTD_N16326_G1_1@meta.data$cell_id %>% head()
xenium_ROI_RCTD_N16326_G1_1@meta.data$cell_id <-paste0(xenium_ROI_RCTD_N16326_G1_1$project,"_",xenium_ROI_RCTD_N16326_G1_1$cell_id_raw)
xenium_ROI_RCTD_N16326_G1_1@meta.data$sample_cell_id <- rownames(xenium_ROI_RCTD_N16326_G1_1@meta.data)

##把all_data_rctd_N16326_G1_1的first_type加入xenium_ROI中-----
all_data_rctd_N16326_G1_1@meta.data %>% head()
source_meta <- all_data_rctd_N16326_G1_1@meta.data[, c("cell_id", "first_type")]
xenium_ROI_RCTD_N16326_G1_1@meta.data <- merge(
  xenium_ROI_RCTD_N16326_G1_1@meta.data,
  source_meta,
  by = "cell_id",
  all.x = TRUE  
)

rownames(xenium_ROI_RCTD_N16326_G1_1@meta.data) <- xenium_ROI_RCTD_N16326_G1_1@meta.data$sample_cell_id 

xenium_ROI_RCTD_N16326_G1_1@meta.data %>% head()

xenium_ROI_RCTD_N16326_G1_1$first_type <- as.character(xenium_ROI_RCTD_N16326_G1_1$first_type)
xenium_ROI_RCTD_N16326_G1_1$first_type[is.na(xenium_ROI_RCTD_N16326_G1_1$first_type)] <- "other"

xenium_ROI_RCTD_N16326_G1_1$first_type <- factor(xenium_ROI_RCTD_N16326_G1_1$first_type)

xenium_ROI_RCTD_N16326_G1_1$first_type %>% unique()

saveRDS(xenium_ROI_RCTD_N16326_G1_1,file = file.path(rctd_outdir, "xenium_ROI_RCTD_N16326_G1_1.rds"))


#RCTD空间映射-------
##Imagedimplot all----
xenium_ROI_RCTD_N16326_G1_1@meta.data$first_type %>% table()
# xenium_ROI_RCTD_N16326_G1_1@meta.data$annotation<- as.factor(xenium_ROI_RCTD_N16326_G1_1@meta.data$annotation)

all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928", # 原始21色
  "#FF6B6B", "#4ECDC4", "#FF9F1C", "#7209B7", "#4361EE", 
  "#F15BB5", "#00BBF9", "#00F5D4", "#9B5DE5", "#00F5A7",
  "#FEE440", "#9C89B8", "#06D6A0"  # 新增13色
)


all_colors10 <- c(
  "B_cell" = "#00F5A7", 
  "Endothelial" = "#98D277",
  "Epithelial" = "#9B5DE5",
  "Fibroblast" = "#FDBB69",
  "Myeloid" = "#feee7d", #"#9E8099",
  "Pericyte" = "#00BBF9",
  "Plasma" =  "#52AF43",
  "SMC" ="#A6CEE3",
  "T_cell" = "#D9A295",
  "Uncertain" = "gray90"
) 


all_colors_rctd <- c(
  # CAF 家族（暖橙/棕色系 —— 相近）
  "Antigen-presenting CAF"   = "#FDBB69",
  "Inflammatory CAF"         = "#FE982C",
  "hsp tCAF"                 = "#F78620",
  "Matrix CAF"               = "#DD9A88",
  "tCAF"                     = "#D9A295",
  "Interferon CAF"           = "#DBB466",
  "Vascular CAF"             = "#B15928",
  
  # 上皮细胞 Epi 家族（紫色/粉色系 —— 相近）
  "Epithelial cells"         = "#9B5DE5",
  "AAEs"                     = "#600473",
  "non-AAEs"                 = "#9E8099",
  # "FABP5 Epi"                = "#7D54A5",
  # "KRT5 Epi"                 = "#B294C7",
  # "MUC16 Epi"                = "#600473",
  # "SPINK1 Epi"               = "#7209B7",
  # "SPRR3 Epi"                = "#9E8099",
  # "WFDC2 Epi"                = "#F15BB5",
  
  # 免疫髓系 Myeloid（黄色系 —— 相近）
  "Macrophages"              = "#feee7d",
  "Neutrophils"              = "#F0EB99",
  
  # 免疫淋巴系（青色/蓝色系 —— 相近）
  "B cells"                  = "#00F5A7",
  "Plasma cells"             = "#52AF43",
  "T_NK"                     = "#4361EE",
  
  # 基质/血管（绿色系 —— 相近）
  "Endothelial cells"        = "#98D277",
  "Smooth muscle cells"       = "#A6CEE3",
  
  # 其他
  "other"                    = "gray90"
)


sample_label <- c("XENIUM6-13")#unique(xenium_ROI_RCTD_N16326_G1_1@meta.data$sample)

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium_ROI_RCTD_N16326_G1_1@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_ROI_RCTD_N16326_G1_1,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors_rctd,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "first_type"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0(rctd_outdir, "Imagedimplot/",i,"_rctd_spatial.png"),width = 15 ,height = 15,dpi = 300)
  ggsave(p,filename = paste0(rctd_outdir, "Imagedimplot/",i,"_rctd_spatial.pdf"),width = 15 ,height = 15)
}


##Imagedimplot 4types-----
all_colors_rctd_5types <- c(
  # CAF 家族（暖橙/棕色系 —— 相近）
  "Antigen-presenting CAF"   = "gray90",
  "Inflammatory CAF"         = "gray90",
  "hsp tCAF"                 = "#F78620",
  "Matrix CAF"               = "gray90",
  "tCAF"                     = "gray90",
  "Interferon CAF"           = "gray90",
  "Vascular CAF"             = "#B15928",
  
  # 上皮细胞 Epi 家族（紫色/粉色系 —— 相近）
  "Epithelial cells"         = "gray90",
  "AAEs"                     = "#600473",
  "non-AAEs"                 = "#9E8099",
  # "FABP5 Epi"                = "#7D54A5",
  # "KRT5 Epi"                 = "#B294C7",
  # "MUC16 Epi"                = "#600473",
  # "SPINK1 Epi"               = "#7209B7",
  # "SPRR3 Epi"                = "#9E8099",
  # "WFDC2 Epi"                = "#F15BB5",
  
  # 免疫髓系 Myeloid（黄色系 —— 相近）
  "Macrophages"              = "gray90",
  "Neutrophils"              = "gray90",
  
  # 免疫淋巴系（青色/蓝色系 —— 相近）
  "B cells"                  = "gray90",
  "Plasma cells"             = "gray90",
  "T_NK"                     = "gray90",
  
  # 基质/血管（绿色系 —— 相近）
  "Endothelial cells"        = "#98D277",
  "Smooth muscle cells"      = "gray90",
  
  # 其他
  "other"                    = "gray90"
)

Idents(xenium_ROI_RCTD_N16326_G1_1) <- "first_type"

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium_ROI_RCTD_N16326_G1_1@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_ROI_RCTD_N16326_G1_1,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors_rctd,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "first_type",
                    cells = WhichCells(xenium_ROI_RCTD_N16326_G1_1, idents = c("AAEs",     
                                                          "non-AAEs",
                                                          "hsp tCAF",      
                                                          "Vascular CAF",
                                                          "Interferon CAF",
                                                          "Endothelial cells"
                                                          ) )
                    
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  ggsave(p,filename = paste0(rctd_outdir,"Imagedimplot/",i,"_rctd_6types_spatial.png"),width = 15 ,height = 15,dpi = 300)
  ggsave(p,filename = paste0(rctd_outdir,"Imagedimplot/",i,"_rctd_6types_spatial.pdf"),width = 15 ,height = 15)
}


##细胞比例统计------
xenium_ROI_RCTD_N16326_G1_1$group %>% table()

library(tidyverse)
library(Seurat)
library(ggalluvial)
library(ggthemes)
library(ggpubr)

dat_props_immune <- xenium_ROI_RCTD_N16326_G1_1[[]] %>% 
  filter(!first_type %in% c("other")) %>% 
  filter(sample %in% c("XENIUM6-13")) %>% 
  group_by(sample, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

dat_props_immune %>% head
# # A tibble: 6 × 4
# # Groups:   group [1]
# group annotation  count   freq
# <fct> <fct>       <int>  <dbl>
#   1 Ctrl  B cell       1437 0.153 
# 2 Ctrl  Plasma cell  1515 0.162 
# 3 Ctrl  T cell       4222 0.451 
# 4 Ctrl  NK cell       152 0.0162
# 5 Ctrl  Macrophage   1192 0.127 
# 6 Ctrl  Neutrophil    253 0.0270


p_props_immune_flow <- 
  dat_props_immune %>% 
  # mutate(id = 1:n()) %>%
  
  ggplot(aes(x = sample, y = freq, stratum = first_type, alluvium = first_type, fill = first_type, label = first_type)) +
  geom_stratum(alpha = 1, color = NA) +
  geom_flow(alpha = 0.5) +
  #geom_text(stat = "stratum") +
  theme_classic2() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1#, 
                               #vjust = 0.5
    ) 
  ) +
  scale_fill_manual(values = all_colors_rctd) +
  #scale_x_discrete(expand = c(0.01, 0.01)) +
  #scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = "sample",
    y = "cell proportion"
  )


ggsave(file = file.path(rctd_outdir, "props_flow.pdf"), p_props_immune_flow, width = 6, height = 6)
ggsave(file = file.path(rctd_outdir, "props_flow.png"), p_props_immune_flow, width = 6, height = 6, dpi = 300)


#邻域分析N16326_G1_1内-----
library(Seurat)
library(tidyverse)
source("spatial analysis.R")

CN_out_dir <- file.path(rctd_outdir, "CN analysis/")
dir.create(CN_out_dir)

all_data_rctd_N16326_G1_1@meta.data %>% head()

cells_position_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_position_ls.rds")
cells_position_df <- cells_position_ls %>% bind_rows()

all_data_submeta_for_cn_analysis <- 
  all_data_rctd_N16326_G1_1[[c("project", "group", "first_type")]] %>% 
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    first_type = factor(first_type), # import!
    obj_id = cell_id # import!
  )

saveRDS(all_data_submeta_for_cn_analysis,file.path(CN_out_dir,"all_data_submeta_for_cn_analysis.rds"))

## CN k10 ---------
res_cn_10 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "first_type",
  k.nn           = 10,
  k.cluster      = 10
)

saveRDS(res_cn_10, file.path(CN_out_dir, "res_cn_10.rds"))

### heatmap ----------------------
# FoldChange cluster/allCell
cn10_cluster_mean <- res_cn_10$cluster_centers
cn10_baseMean     <- res_cn_10$base_mean

# log2(FC)
cn10_cluster_FC <- log2(1 + sweep(cn10_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn10_baseMean))

# heatmap
p_cn10_heat <-
  pheatmap::pheatmap(
    cn10_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn10_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_10_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn10_heat$gtable)
dev.off()

### CN10 spatial plot -----------------
cn10_sp_out_dir <- file.path(CN_out_dir, 'CN10 spatial plot/') 
dir.create(cn10_sp_out_dir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

all_data_rctd_N16326_G1_1@meta.data$cell_id <- NULL
rownames(all_data_rctd_N16326_G1_1@meta.data) %>% head()

All_type_meta_cn10 <- 
  all_data_rctd_N16326_G1_1[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_10$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn10, file.path(CN_out_dir, "all_data_meta_cn10.rds"))

cn10_sp_p_ls <-
  cells_boundary_ls["N16326_G1_1"] %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn10[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn10, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })


# export image
cn10_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn10_sp_out_dir, str_glue("CN10 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


### cell props ------
All_type_cn10_count_df <-
  All_type_meta_cn10 %>% 
  group_by(sample, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn10_props <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn10_props_group <-
  All_type_cn10_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN10_barplot.pdf"),#boxplot_and_
  height = 6,
  width  = 3
)
p_cn10_props %>% print
#p_cn10_props_group %>% print
dev.off()


###CN空间映射------
All_type_meta_cn10 %>% head()

source_meta <- All_type_meta_cn10[, c("cell_id", "CN")] %>% 
  rename(CN10 = CN)

xenium_ROI_RCTD_N16326_G1_1 <- readRDS(file.path(rctd_outdir, "xenium_ROI_RCTD_N16326_G1_1.rds"))

xenium_ROI_RCTD_N16326_G1_1@meta.data <- merge(
  xenium_ROI_RCTD_N16326_G1_1@meta.data,
  source_meta,
  by = "cell_id",
  all.x = TRUE  
)

rownames(xenium_ROI_RCTD_N16326_G1_1@meta.data) <- xenium_ROI_RCTD_N16326_G1_1@meta.data$sample_cell_id 

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium_ROI_RCTD_N16326_G1_1@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_ROI_RCTD_N16326_G1_1,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "CN10"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )+
    ggsci::scale_fill_igv() 
  
  ggsave(p,filename = paste0(cn10_sp_out_dir, i,"_CN10_spatial.png"),width = 15 ,height = 15,dpi = 300)
  ggsave(p,filename = paste0(cn10_sp_out_dir, i,"_CN10_spatial.pdf"),width = 15 ,height = 15)
}



## CN k50 ---------
res_cn_50 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "first_type",
  k.nn           = 50,
  k.cluster      = 10
)

saveRDS(res_cn_50, file.path(CN_out_dir, "res_cn_50.rds"))

### heatmap ----------------------
# FoldChange cluster/allCell
cn50_cluster_mean <- res_cn_50$cluster_centers
cn50_baseMean     <- res_cn_50$base_mean

# log2(FC)
cn50_cluster_FC <- log2(1 + sweep(cn50_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn50_baseMean))

# heatmap
p_cn50_heat <-
  pheatmap::pheatmap(
    cn50_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn50_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "CN_50_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn50_heat$gtable)
dev.off()

### CN50 spatial plot -----------------
cn50_sp_out_dir <- file.path(CN_out_dir, 'CN50 spatial plot/') 
dir.create(cn50_sp_out_dir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

all_data_rctd_N16326_G1_1@meta.data$cell_id <- NULL
rownames(all_data_rctd_N16326_G1_1@meta.data) %>% head()

All_type_meta_cn50 <- 
  all_data_rctd_N16326_G1_1[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn_50$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_50$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn50, file.path(CN_out_dir, "all_data_meta_cn50.rds"))

cn50_sp_p_ls <-
  cells_boundary_ls["N16326_G1_1"] %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn50[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn50, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })


# export image
cn50_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn50_sp_out_dir, str_glue("CN50 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


### cell props ------
All_type_cn50_count_df <-
  All_type_meta_cn50 %>% 
  group_by(sample, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn50_props <-
  All_type_cn50_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn50_props_group <-
  All_type_cn50_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "CN50_barplot.pdf"),#boxplot_and_
  height = 6,
  width  = 3
)
p_cn50_props %>% print
#p_cn50_props_group %>% print
dev.off()


###CN空间映射------
All_type_meta_cn50 %>% head()

source_meta <- All_type_meta_cn50[, c("cell_id", "CN")] %>% 
  rename(CN50 = CN)

xenium_ROI_RCTD_N16326_G1_1 <- readRDS(file.path(rctd_outdir, "xenium_ROI_RCTD_N16326_G1_1.rds"))

xenium_ROI_RCTD_N16326_G1_1@meta.data <- merge(
  xenium_ROI_RCTD_N16326_G1_1@meta.data,
  source_meta,
  by = "cell_id",
  all.x = TRUE  
)

rownames(xenium_ROI_RCTD_N16326_G1_1@meta.data) <- xenium_ROI_RCTD_N16326_G1_1@meta.data$sample_cell_id 

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium_ROI_RCTD_N16326_G1_1@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_ROI_RCTD_N16326_G1_1,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "CN50"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )+
    ggsci::scale_fill_igv() 
  
  ggsave(p,filename = paste0(cn50_sp_out_dir, i,"_CN50_spatial.png"),width = 15 ,height = 15,dpi = 300)
  ggsave(p,filename = paste0(cn50_sp_out_dir, i,"_CN50_spatial.pdf"),width = 15 ,height = 15)
}



## CN k100 ---------
res_cn_100 <- CN_compute_CN(
  all_data_submeta_for_cn_analysis,
  x = "X",
  y = "Y",
  group.by       = "project",
  classification = "first_type",
  k.nn           = 100,
  k.cluster      = 10
)

saveRDS(res_cn_100, file.path(CN_out_dir, "res_cn_100.rds"))

### heatmap ----------------------
# FoldChange cluster/allCell
cn100_cluster_mean <- res_cn_100$cluster_centers
cn100_baseMean     <- res_cn_100$base_mean

# log2(FC)
cn100_cluster_FC <- log2(1 + sweep(cn100_cluster_mean, MARGIN = 2, FUN = "/", STATS = cn100_baseMean))

# heatmap
p_cn100_heat <-
  pheatmap::pheatmap(
    cn100_cluster_FC,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color = NA,
    treeheight_row = 0,
    labels_row = str_c("CN", rownames(cn100_cluster_FC), sep = ""),
    # clustering_method = 'ward.D2',
    cluster_rows = F,
    scale = "column",
    silent = TRUE
  )

pdf(
  file.path(CN_out_dir, "cn_100_heatmap.pdf"),
  height = 6,
  width  = 10
)
cowplot::plot_grid(p_cn100_heat$gtable)
dev.off()

### cn100 spatial plot -----------------
cn100_sp_out_dir <- file.path(CN_out_dir, 'cn100 spatial plot/') 
dir.create(cn100_sp_out_dir)

cells_boundary_ls <- readRDS("R_OUTPUT/data_for_analysis/cells_boundary_ls.rds")
# image_wh_r_ls     <- readRDS('data_for_analysis/image_wh_r_ls.rds')

all_data_rctd_N16326_G1_1@meta.data$cell_id <- NULL
rownames(all_data_rctd_N16326_G1_1@meta.data) %>% head()

All_type_meta_cn100 <- 
  all_data_rctd_N16326_G1_1[[]] %>%
  rownames_to_column('cell_id') %>%
  # left_join(rownames_to_column(res_cn_100$cn_cluster, "cell_id"), by = "cell_id") %>%
  # mutate(CN2 = fct_recode(CN, !!!structure(names(cn_annotation), names = cn_annotation)))
  left_join(rownames_to_column(res_cn_100$cn_cluster, "cell_id"), by = "cell_id")

saveRDS(All_type_meta_cn100, file.path(CN_out_dir, "all_data_meta_cn100.rds"))

cn100_sp_p_ls <-
  cells_boundary_ls["N16326_G1_1"] %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- All_type_meta_cn100[c("project", "group")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(All_type_meta_cn100, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = CN, group = cell_id))  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal() + 
      ggsci::scale_fill_igv() +
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, sep = "_") )
    p
  })


# export image
cn100_sp_p_ls %>% 
  iwalk(function(p, name){
    pdf(
      file.path(cn100_sp_out_dir, str_glue("cn100 spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


### cell props ------
All_type_cn100_count_df <-
  All_type_meta_cn100 %>% 
  group_by(sample, group, CN) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count))

p_cn100_props <-
  All_type_cn100_count_df %>%
  
  ggplot(aes(x = sample, y = freq, fill = CN)) +
  geom_col() +
  # ggsci::scale_fill_aaas() +
  ggsci::scale_fill_igv() + 
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )


p_cn100_props_group <-
  All_type_cn100_count_df %>%
  
  ggplot(aes(x = CN, y = freq, fill = group)) +
  geom_boxplot() +
  ggsci::scale_fill_aaas() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL
  )

pdf(
  file.path(CN_out_dir, "cn100_barplot.pdf"),#boxplot_and_
  height = 6,
  width  = 3
)
p_cn100_props %>% print
#p_cn100_props_group %>% print
dev.off()


###CN空间映射------
All_type_meta_cn100 %>% head()

source_meta <- All_type_meta_cn100[, c("cell_id", "CN")] %>% 
  rename(cn100 = CN)

xenium_ROI_RCTD_N16326_G1_1 <- readRDS(file.path(rctd_outdir, "xenium_ROI_RCTD_N16326_G1_1.rds"))

xenium_ROI_RCTD_N16326_G1_1@meta.data <- merge(
  xenium_ROI_RCTD_N16326_G1_1@meta.data,
  source_meta,
  by = "cell_id",
  all.x = TRUE  
)

rownames(xenium_ROI_RCTD_N16326_G1_1@meta.data) <- xenium_ROI_RCTD_N16326_G1_1@meta.data$sample_cell_id 

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium_ROI_RCTD_N16326_G1_1@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_ROI_RCTD_N16326_G1_1,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "cn100"
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )+
    ggsci::scale_fill_igv() 
  
  ggsave(p,filename = paste0(cn100_sp_out_dir, i,"_cn100_spatial.png"),width = 15 ,height = 15,dpi = 300)
  ggsave(p,filename = paste0(cn100_sp_out_dir, i,"_cn100_spatial.pdf"),width = 15 ,height = 15)
}







