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

rctd_outdir <- "RCTD/"
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

# 2. create the RCTD reference object------
reference <- Reference(
  counts_matrix, 
  droplevels(cluster_sc), 
  nUMI_sc
)

saveRDS(reference,file = file.path(rctd_outdir, "reference.rds"))
reference <- readRDS

#xenium数据
all_data <- readRDS("result/all_data_Lgr5_M1M2.rds")
all_data$sample %>% table()

# 3. create the RCTD query object-----
all_data_sp_rctd_ls <- all_data %>% 
  JoinLayers() %>% 
  SplitObject('sample')

cells_position_ls <- readRDS("~/R_OUTPUT/data_for_analysis/cells_position_ls.rds")

rctd_query_ls <-
  all_data_sp_rctd_ls %>%
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

saveRDS(rctd_query_ls,file = file.path(rctd_outdir, "rctd_query_ls.rds"))

#4. run RCTD-----
rctd_result_ls <-
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 30, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
    
  })

saveRDS(rctd_result_ls, file = file.path(rctd_outdir, "rctd_result_ls.rds"))
rctd_result_ls <- readRDS(file = file.path(rctd_outdir,"rctd_result_ls.rds"))

#5. addMetaData-----
all_data_sp_rctd_anno_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)
    
    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })

saveRDS(all_data_sp_rctd_anno_ls,file = file.path(rctd_outdir, "rctd.rds"))

#all_data_sp_rctd_anno_ls <- readRDS("rctd.rds")


#6. dimplot 
# all_data_sp_anno <-
#   merge(
#     all_data_sp_rctd_anno_ls[[1]],
#     all_data_sp_rctd_anno_ls[-1]
#   )

meta_rctd_result <-
  rctd_result_ls %>%
  map(~ .@results$results_df) %>%
  bind_rows()

all_data_rctd <-
  all_data %>%
  AddMetaData(meta_rctd_result)

all_data_rctd$first_type %>% table()

all_data_rctd$first_type <- as.character(all_data_rctd$first_type)
all_data_rctd$first_type[is.na(all_data_rctd$first_type)] <- "other"
all_data_rctd$first_type <- factor(all_data_rctd$first_type)

saveRDS(all_data_rctd,file = file.path(rctd_outdir, "all_data_rctd.rds"))

all_data_rctd <- readRDS(file.path(rctd_outdir, "all_data_rctd.rds"))

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


#all_data_rctd <- readRDS(file = file.path(rctd_outdir, "all_data_rctd.rds"))
all_data_rctd <- JoinLayers(all_data_rctd)
all_data_rctd$cell_id %>% head()
all_data_rctd$cell_id_rownames <- rownames(all_data_rctd@meta.data)
all_data_rctd$cell_id_raw <- all_data_rctd$cell_id
all_data_rctd$cell_id <- paste0(all_data_rctd$project,"_",all_data_rctd$cell_id_raw)



##包含images的xenium数据-------
xenium <- readRDS("result/xenium_FOV.rds")

xenium_RCTD <- xenium
xenium_RCTD@meta.data %>% head()
xenium_RCTD@meta.data$cell_id_raw %>% head()
xenium_RCTD@meta.data$cell_id %>% head()
xenium_RCTD@meta.data$cell_id <-paste0(xenium_RCTD$project,"_",xenium_RCTD$cell_id_raw)
xenium_RCTD@meta.data$sample_cell_id <- rownames(xenium_RCTD@meta.data)

##rctd的first_type加入xenium中-----
all_data_rctd@meta.data %>% head()
source_meta <- all_data_rctd@meta.data[, c("cell_id", "first_type")]
xenium_RCTD@meta.data <- merge(
  xenium_RCTD@meta.data,
  source_meta,
  by = "cell_id",
  all.x = TRUE  
)

rownames(xenium_RCTD@meta.data) <- xenium_RCTD@meta.data$sample_cell_id 
xenium_RCTD@meta.data %>% head()
xenium_RCTD$first_type <- as.character(xenium_RCTD$first_type)
xenium_RCTD$first_type[is.na(xenium_RCTD$first_type)] <- "other"
xenium_RCTD$first_type <- factor(xenium_RCTD$first_type)
xenium_RCTD$first_type %>% unique()

saveRDS(xenium_RCTD,file = file.path(rctd_outdir, "xenium_RCTD.rds"))


##RCTD空间映射-------
####Imagedimplot all----
xenium_RCTD@meta.data$first_type %>% table()
# xenium_RCTD@meta.data$annotation<- as.factor(xenium_RCTD@meta.data$annotation)

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
  #"tCAF"                     = "#D9A295",
  "Interferon CAF"           = "#DBB466",
  "Vascular CAF"             = "#B15928",
  
  # 上皮细胞 Epi 家族（紫色/粉色系 —— 相近）
  #"Epithelial cells"         = "#9B5DE5",
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


sample_label <- unique(xenium_RCTD@meta.data$sample)

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium_RCTD@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_RCTD,
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


####Imagedimplot 4types-----
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

Idents(xenium_RCTD) <- "first_type"

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium_RCTD@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_RCTD,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = all_colors_rctd,
                    #shuffle.cols = F, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "first_type",
                    cells = WhichCells(xenium_RCTD, idents = c("AAEs",     
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
xenium_RCTD$group %>% table()

library(tidyverse)
library(Seurat)
library(ggalluvial)
library(ggthemes)
library(ggpubr)

dat_props_immune <- xenium_RCTD[[]] %>% 
  filter(!first_type %in% c("other")) %>% 
  #filter(sample %in% c("XENIUM6-13")) %>% 
  group_by(group, first_type) %>%
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
  
  ggplot(aes(x = group, y = freq, stratum = first_type, alluvium = first_type, fill = first_type, label = first_type)) +
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


