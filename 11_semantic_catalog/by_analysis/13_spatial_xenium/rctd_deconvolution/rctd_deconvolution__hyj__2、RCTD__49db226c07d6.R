# RCTD

library(tidyverse)
library(Seurat)
getwd()
library(patchwork)
options(java.parameters = "-Xmx10g") # options('java.parameters')
library(rJava)
library(RBioFormats)
library(ggpubr)
library(SeuratDisk)
library(ggrepel)
#library(loupeR)
library(spacexr)

rctd_outdir <- "RCTD"
dir.create(rctd_outdir)


# 1.单细胞数据-----
## 1.1、获取原始counts----
### h5ad文件转化成rds用conda env tools里的sceasy包convertFormat()------
# library('sceasy')
# python <- use_python("/mnt/NexData_file/home/liuran/rstudio/miniforge3/envs/rpy/bin/python")
#
# reticulate::py_config()
# sc <- import("scanpy")

adata <- sc$read("/mnt/NexData_file/home/zwl/Project/FWSC20252451/all_data.h5ad")

# sc_counts <- convertFormat('/mnt/NexData_file/scRNA/scVision_output/FWSC20253251/03-Cell_filter/adata_counts.h5ad', to = 'seurat')
# adata <- sc_counts
# saveRDS(adata, 'sc_counts.rds')

counts_x <- adata$X %>%
  t() %>%
  as.matrix()
colnames(counts_x) <- rownames(adata$obs)
rownames(counts_x) <- rownames(adata$var)

counts_x[-20:-1, -20:-1]

## 1.2、metadata-----
adata_meta <- adata$obs
adata_meta$celltype %>% table()
anno_counts <- table(adata_meta$cell_type)
print(anno_counts)

scRNA2 <- CreateSeuratObject(
  counts_x,
  meta.data = adata_meta
)

# scRNA2@meta.data <- all_data_sc@meta.data
saveRDS(scRNA2, "RCTD/scRNA2.rds")
scRNA2$cell_type %>% table()
scRNA2 <- subset(scRNA2, cell_type != "Neurons")


Idents(scRNA2) <- "cell_type"
counts_matrix <- GetAssayData(scRNA2, "RNA", "counts")
cluster_sc <- as.factor(scRNA2$cell_type)
nUMI_sc <- scRNA2$total_counts

# sc_cell_filter_idx <-  cluster_sc != "Unknown"

# 2.create the RCTD reference object------
reference <- Reference(
  counts_matrix,
  droplevels(cluster_sc),
  nUMI_sc
)

saveRDS(reference, "RCTD/reference.rds")
# xenium数据
all_data <- readRDS("R_OUTPUT/all_data.rds")
# 3. create the RCTD query object-----
all_data_sp_rctd_ls <- all_data %>%
  JoinLayers() %>%
  SplitObject("sample")

cells_position_ls <- readRDS("~/R_OUTPUT/data_for_analysis/cells_position_ls.rds")

rctd_query_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name) { # browser()
    counts_hd <- xenium %>% LayerData(layer = "counts")
    coords <-
      cells_position_ls[[name]] %>%
      dplyr::select(
        cell_id,
        x = x_centroid,
        y = y_centroid
      ) %>%
      column_to_rownames("cell_id")

    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })

saveRDS(rctd_query_ls, "RCTD/rctd_query_ls.rds")

rm(all_data_ls)
rm()

reference <- readRDS("/mnt/NexData_file/home/hyj/Rstudio/XE2505092086-01_yangrui/RCTD/reference.rds")
rctd_query_ls <- readRDS("/mnt/NexData_file/home/hyj/Rstudio/XE2505092086-01_yangrui/RCTD/rctd_query_ls.rds")

# 4. run RCTD-----
rctd_result_ls <-
  rctd_query_ls %>%
  imap(function(query, name) {
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 30, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

    RCTD
  })

saveRDS(rctd_result_ls, "RCTD/rctd_result_ls.rds")

rctd_result_ls <- readRDS("RCTD/rctd_result_ls.rds")

# 5. addMetaData-----
all_data_sp_rctd_anno_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name) { # browser()
    rctd <- rctd_result_ls[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)

    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })

saveRDS(all_data_sp_rctd_anno_ls, "RCTD/rctd.rds")

# all_data_sp_rctd_anno_ls <- readRDS("rctd.rds")


# 6. dimplot
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

saveRDS(all_data_rctd, "RCTD/all_data_rctd.rds")


# all_data_rctd <- readRDS("RCTD/all_data_rctd.rds")
all_data_rctd <- JoinLayers(all_data_rctd)
all_data_rctd$cell_id %>% head()
all_data_rctd$cell_id_raw <- all_data_rctd$cell_id
all_data_rctd$cell_id <- paste0(all_data_rctd$project, "_", all_data_rctd$cell_id_raw)

#加入images的xenium数据----
xenium <- readRDS("result/xenium.rds")

xenium_RCTD <- xenium
xenium_RCTD@meta.data %>% head()
xenium_RCTD@meta.data$cell_id_raw %>% head()
xenium_RCTD@meta.data$cell_id %>% head()
xenium_RCTD$cell_id_raw <- xenium_RCTD$cell_id
xenium_RCTD@meta.data$cell_id <- paste0(xenium_RCTD$project, "_", xenium_RCTD$cell_id_raw)
xenium_RCTD@meta.data$sample_cell_id <- rownames(xenium_RCTD@meta.data)

#把all_data_rctd的first_type加入xenium中-----
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

xenium_RCTD$first_type %>% table()

saveRDS(xenium_RCTD, "RCTD/xenium_RCTD.rds")


# RCTD空间映射-------
## Imagedimplot all----
xenium_RCTD@meta.data %>% head()
# xenium_RCTD@meta.data$annotation<- as.factor(xenium_RCTD@meta.data$annotation)

all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928", # 原始21色
  "#FF6B6B", "#4ECDC4", "#FF9F1C", "#7209B7", "#4361EE",
  "#F15BB5", "#00BBF9", "#00F5D4", "#9B5DE5", "#00F5A7",
  "#FEE440", "#9C89B8", "#06D6A0" # 新增13色
)


all_colors <- c(
  # 1. 神经元类（蓝色系）
  #"Neurons" = "#2D82AF", # 主神经元
  "Excitatory neurons" = "#5B9EC9", # 兴奋性神经元
  "Inhibitory neurons" = "#A6CEE3", # 抑制性神经元

  # 2. 胶质/支持细胞（绿色系）
  "Astrocyte" = "#7EBA98", # 星形胶质细胞
  "Oligodendrocyte" = "#98D277", # 少突胶质细胞
  "OPC" = "#52AF43", # 少突胶质前体细胞
  "Ependymal cells" = "#6F9E4C", # 室管膜细胞

  # 3. 小胶质/髓系（橙棕系）
  "M1-like Mic" = "#FDBB69", # M1型小胶质
  "M2-like Mic" = "#FE982C", # M2型小胶质
  "Nrg3+ Mic" = "#F78620", # Nrg3+小胶质
  "Pde10a+ Mic" = "#DD9A88", # Pde10a+小胶质
  "Pde4b+ Mic" = "#D9A295", # Pde4b+小胶质

  # 4. 基质/内皮/其他（紫/粉/杂色）
  "Endothelium" = "#600473", # 内皮细胞
  "Fibroblast" = "#B294C7" # 成纤维细胞
)

sample_label <- unique(xenium_RCTD@meta.data$sample)

for (i in sample_label) {
  # j = gsub("_","-",i)
  DefaultBoundary(xenium_RCTD@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium_RCTD,
    fov = i,
    nmols = 20000,
    dark.background = F,
    border.color = "black",
    border.size = 0.1,
    cols = all_colors,
    # shuffle.cols = F, #随机打乱颜色顺序
    # size = 1.2,
    group.by = "first_type"
    # molecules = c("INHBA"),
    # mols.size = 1,
    # mols.cols = "red",
    # mols.alpha = 0.5, #分子的透明度值
  )
  ggsave(p, filename = paste0("RCTD/Imagedimplot/", i, "_rctd_spatial.png"), width = 12, height = 12, dpi = 300)
  ggsave(p, filename = paste0("RCTD/Imagedimplot/", i, "_rctd_spatial.pdf"), width = 12, height = 12)
}


## Imagefeatureplot-----
SegImagefeature_dir <- 'RCTD/SegImagefeature'
dir.create(SegImagefeature_dir)

gene_list <- c('P2ry12','C1qb','Ntrk1','P2y12','Hork3','Sp1999')
valid_gene <- gene_list[gene_list %in% rownames(xenium_RCTD)]
cat("有效基因:", length(valid_gene), "/", length(gene_list))

missing_gene <- setdiff(gene_list, valid_gene)

gene_cus <- valid_gene

for (sample_id in unique(xenium_RCTD@meta.data$sample)) {
  # 智能匹配样本ID
  j <- names(xenium_RCTD@images)[grepl(sample_id, names(xenium_RCTD@images))]
  
  if (length(j) == 0) {
    warning("跳过未匹配样本: ", sample_id)
    next
  }
  
  for (gene in gene_cus) {
    tryCatch({
      # 空值检查
      if (!is.null(xenium_RCTD@images[[j]])) {
        DefaultBoundary(xenium_RCTD@images[[j]]) <- "segmentation"
      }
      
      p <- ImageFeaturePlot(
        object = xenium_RCTD,
        fov = j,
        nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
        border.color = "black",  # 细胞边界的颜色（深蓝色）
        border.size = 0.1,    # 细胞边界的线宽
        features = gene,
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('#003366', "#E7475E"),
        size = 0.8,
        dark.background = F
      )
      
      ggsave(
        filename = file.path(SegImagefeature_dir, paste0(sample_id, "_", gene, "_spatial.png")),
        plot = p,
        width = 12,
        height = 12,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(SegImagefeature_dir, paste0(sample_id, "_", gene, "_spatial.pdf")),
        plot = p,
        width = 12,
        height = 12
      )
      
      message("成功保存: ", sample_id, "_", gene)
      
    }, error = function(e) {
      warning("失败[", sample_id, ":", gene, "]: ", e$message)
    })
  }
}


source("ImageFeature-plot.R")

Imagefeature_dir <- 'RCTD/Boundary_Imagefeature' 
dir.create(Imagefeature_dir)


xenium_RCTD$first_type %>% table()
##Imagefeatureplot➕边细胞界--------
for (sample_id in unique(xenium@meta.data$sample)){
  p <- ImageFeaturePlot_wyt(xenium_RCTD, 
                            features = c("Ntrk1"), 
                            size = 0.5, 
                           # cols = c('#003366', '#FFD700'),
                            cols = c('#003366', "#E7475E"),
                            fov =sample_id,
                            nmols = 20000,
                            border.color = all_colors, 
                            border.size = 0.1,
                            border.by = "first_type",
                            #cells.use = c("Macrophage"),
                            axes = FALSE,
                            dark.background = F
  )
  
  ggsave(
    filename = file.path(Imagefeature_dir, paste0(sample_id, "_Ntrk1_spatial.png")),
    plot = p,
    width = 12,
    height = 12,
    dpi = 300
  )
  ggsave(
    filename = file.path(Imagefeature_dir, paste0(sample_id, "_Ntrk1_spatial.pdf")),
    plot = p,
    width = 12,
    height = 12
  )
  
  message("成功保存: ", sample_id)
}

###合并为大类--------
xenium_RCTD@meta.data <- xenium_RCTD@meta.data %>% 
  mutate(first_type_merge = case_when(
    first_type%in% c("Excitatory neurons","Inhibitory neurons") ~ "Neurons",
    first_type%in% c("M1-like Mic", "M2-like Mic", "Nrg3+ Mic", "Pde10a+ Mic", "Pde4b+ Mic" ) ~ "Microglias",
    TRUE ~ first_type
    ))
xenium_RCTD@meta.data$first_type_merge %>% table()
table(xenium_RCTD$first_type_merge, xenium_RCTD$first_type)

all_colors <- c(
  "#A6CEE3", "#5B9EC9", "#2D82AF", "#7EBA98", "#98D277",
  "#52AF43", "#6F9E4C", "#DD9A88", "#FDBB69", "#FE982C",
  "#F78620", "#D9A295", "#feee7d", "#02c9c9", "#600473",
  "#B294C7", "#7D54A5", "#9E8099", "#F0EB99", "#DBB466",
  "#B15928", # 原始21色
  "#FF6B6B", "#4ECDC4", "#FF9F1C", "#7209B7", "#4361EE",
  "#F15BB5", "#00BBF9", "#00F5D4", "#9B5DE5", "#00F5A7",
  "#FEE440", "#9C89B8", "#06D6A0" # 新增13色
)



all_colors <- c(
  # 1. 神经元类（蓝色系）
  "Neurons" = "#600473", # 主神经元
  # "Excitatory neurons" = "#5B9EC9", # 兴奋性神经元
  # "Inhibitory neurons" = "#A6CEE3", # 抑制性神经元
  
  # 2. 胶质/支持细胞（绿色系）
  "Astrocyte" = "#00F5A7", # 星形胶质细胞
  "Oligodendrocyte" = "#4ECDC9", # 少突胶质细胞
  "OPC" = "#52AF43", # 少突胶质前体细胞
  "Ependymal cells" = "#4361EE", # 室管膜细胞
  
  # 3. 小胶质/髓系（橙棕系）
  "Microglias" = "#F78620",
  # "M1-like Mic" = "#FDBB69", # M1型小胶质
  # "M2-like Mic" = "#FE982C", # M2型小胶质
  # "Nrg3+ Mic" = "#F78620", # Nrg3+小胶质
  # "Pde10a+ Mic" = "#DD9A88", # Pde10a+小胶质
  # "Pde4b+ Mic" = "#D9A295", # Pde4b+小胶质
  
  # 4. 基质/内皮/其他（紫/粉/杂色）
  "Endothelium" ="#D9A295",# 内皮细胞
  "Fibroblast" = "#B294C7" # 成纤维细胞
)

Imagefeature_dir <- 'RCTD/Boundary_Imagefeature_firsttype_merge' 
dir.create(Imagefeature_dir)

##Imagefeatureplot➕边细胞界--------
for (sample_id in unique(xenium@meta.data$sample)){
  p <- ImageFeaturePlot_wyt(xenium_RCTD, 
                            features = c("Ntrk1"), 
                            size = 0.5, 
                            #cols = c('#003366', '#FFD700'),
                            cols = c("#A6CEE3","#FF6B6B"),
                            fov =sample_id,
                            nmols = 20000,
                            border.color = all_colors, 
                            border.size = 0.1,
                            border.by = "first_type_merge",
                            #cells.use = c("Macrophage"),
                            axes = FALSE,
                            dark.background = F
  )
  
  ggsave(
    filename = file.path(Imagefeature_dir, paste0(sample_id, "_Ntrk1_spatial.png")),
    plot = p,
    width = 12,
    height = 12,
    dpi = 300
  )
  ggsave(
    filename = file.path(Imagefeature_dir, paste0(sample_id, "_Ntrk1_spatial.pdf")),
    plot = p,
    width = 12,
    height = 12
  )
  
  message("成功保存: ", sample_id)
}


## Ntrk1 高表达细胞-------
D <- DotPlot(xenium_RCTD,features = "Ntrk1", group.by = "first_type_merge")

V <- VlnPlot(xenium_RCTD,features = "Ntrk1", group.by = "first_type_merge", pt.size = 0.1)




#细胞比例统计------
xenium_RCTD$group %>% table()
group_level <- c("Ctrl_E","NGF_E","IR_E","IRD_E",
                 "Ctrl_L","NGF_L","IR_L","IRD_L")

xenium_RCTD$group <- factor(xenium_RCTD$group,level = group_level)


library(tidyverse)
library(Seurat)
library(ggalluvial)
library(ggthemes)
library(ggpubr)

##总的-------
dat_props_immune <- xenium_RCTD[[]] %>%
  #filter(group %in% c())
  group_by(group, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))

dat_props_immune %>% head()
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
  # geom_text(stat = "stratum") +
  theme_classic2() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45, hjust = 1 # ,
      # vjust = 0.5
    )
  ) +
  scale_fill_manual(values = all_colors) +
  # scale_x_discrete(expand = c(0.01, 0.01)) +
  # scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = "group",
    y = "cell proportion"
  )


ggsave(file = file.path("RCTD/cell_prop", "props_flow.pdf"), p_props_immune_flow, width = 8, height = 6)
ggsave(file = file.path("RCTD/cell_prop", "props_flow.png"), p_props_immune_flow, width = 8, height = 6, dpi = 300)

##早期-------
group_level <- c("Ctrl_E","NGF_E","IR_E","IRD_E",
                 "Ctrl_L","NGF_L","IR_L","IRD_L")


dat_props_immune <- xenium_RCTD[[]] %>%
  filter(group %in% c("Ctrl_E","NGF_E","IR_E","IRD_E")) %>% 
  group_by(group, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))

dat_props_immune %>% head()
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
  # geom_text(stat = "stratum") +
  theme_classic2() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45, hjust = 1 # ,
      # vjust = 0.5
    )
  ) +
  scale_fill_manual(values = all_colors) +
  # scale_x_discrete(expand = c(0.01, 0.01)) +
  # scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = "group",
    y = "cell proportion"
  )


ggsave(file = file.path("RCTD/cell_prop", "Early_props_flow.pdf"), p_props_immune_flow, width = 8, height = 6)
ggsave(file = file.path("RCTD/cell_prop", "Early_props_flow.png"), p_props_immune_flow, width = 8, height = 6, dpi = 300)



##晚期-----
group_level <- c("Ctrl_E","NGF_E","IR_E","IRD_E",
                 "Ctrl_L","NGF_L","IR_L","IRD_L")


dat_props_immune <- xenium_RCTD[[]] %>%
  filter(group %in% c("Ctrl_L","NGF_L","IR_L","IRD_L")) %>% 
  group_by(group, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count))

dat_props_immune %>% head()
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
  # geom_text(stat = "stratum") +
  theme_classic2() +
  theme(
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(
      angle = 45, hjust = 1 # ,
      # vjust = 0.5
    )
  ) +
  scale_fill_manual(values = all_colors) +
  # scale_x_discrete(expand = c(0.01, 0.01)) +
  # scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    x = "group",
    y = "cell proportion"
  )


ggsave(file = file.path("RCTD/cell_prop", "Later_props_flow.pdf"), p_props_immune_flow, width = 8, height = 6)
ggsave(file = file.path("RCTD/cell_prop", "Later_props_flow.png"), p_props_immune_flow, width = 8, height = 6, dpi = 300)

# 1、细胞通讯-------
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
# library(SeuratData)
library(patchwork)

source("result/cellchat_output/netVisual_circle_custom.R")

all_data_rctd$first_type <- as.character(all_data_rctd$first_type)
all_data_rctd$first_type[is.na(all_data_rctd$first_type)] <- "other"

all_data_rctd$first_type <- factor(all_data_rctd$first_type)

all_data_rctd$first_type %>% table()
all_data_rctd$first_type %>% unique()
all_data_rctd$first_type <- as.character(all_data_rctd$first_type)
all_data_rctd@meta.data <- all_data_rctd@meta.data %>%
  mutate(first_type_combine = case_when(
    grepl("neurons", first_type) ~ "Neuron",
    TRUE ~ first_type
  ))

all_data_rctd@meta.data <- all_data_rctd@meta.data %>%
  mutate(first_type_combine = case_when(
    grepl("Mic", first_type) ~ "Microglial",
    TRUE ~ first_type_combine
  ))

all_data_rctd$first_type_combine %>% table()
all_data_rctd$first_type_combine %>% unique()
table(all_data_rctd$first_type, all_data_rctd$first_type_combine)
all_data_rctd$first_type_combine <- as.character(all_data_rctd$first_type_combine)

all_data_rctd$Celltype <- all_data_rctd$first_type_combine

## 循环各分组(神经元、小胶质细胞和OPC)------
types3_cellchat_dir <- file.path("RCTD", "cellchat_output_types3/")
dir.create(types3_cellchat_dir)

cellchat_data <- all_data_rctd[, all_data_rctd$Celltype %in% c("Neuron", "Microglial", "OPC")]
cellchat_data <- JoinLayers(cellchat_data)
cellchat_data$Celltype %>% table()
cellchat_data$Celltype <- as.factor(cellchat_data$Celltype)
cellchat_data$group %>% table()
cellchat_data$group <- as.character(cellchat_data$group)

for (i in unique(cellchat_data$group)) {
  dat <- cellchat_data[, cellchat_data$group == i]
  dat$group %>% table()
  dat$Celltype %>% table()
  dat$samples <- dat$sample
  dat$samples <- factor(dat$samples)
  Idents(dat) <- "Celltype"
  
  
  all_data_meta <- dat[[]]
  all_data_norm <- LayerData(JoinLayers(dat), "data")
  stopifnot(all(rownames(all_data_meta) == colnames(all_data_norm)))
  
  
  cellchat <- createCellChat(
    object = all_data_norm,
    meta = all_data_meta,
    group.by = "Celltype"
  )
  
  
  cellchat <- setIdent(cellchat, ident.use = "Celltype")
  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents))
  
  CellChatDB.use <- CellChatDB.mouse %>%
    subsetDB(
      search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
      key = "annotation"
    )
  
  
  p <- showDatabaseCategory(CellChatDB.use)
  ggsave(file.path(types3_cellchat_dir, "CellChatDB.mmu.png"), width = 10, height = 5, dpi = 300)
  ggsave(file.path(types3_cellchat_dir, "CellChatDB.mmu.pdf"), width = 10, height = 5)
  
  # print(colnames(CellChatDB$interaction))
  # CellChatDB.use <-  # use     Secreted Signaling
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- smoothData(cellchat, adj = PPI.mouse)
  
  ### Part II: Inference of cell-cell communication network
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, population.size = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  ### 信息提取
  df.net <- subsetCommunication(cellchat)
  write.table(df.net,
              file.path(types3_cellchat_dir, paste0("group_", i, "_Communication_info.txt")),
              sep = "\t",
              quote = F,
              row.names = F
  )
  
  ## Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))
  head(cellchat@net$count)
  head(cellchat@net$weight)
  
  
  pdf(
    file = file.path(types3_cellchat_dir, paste0("group_", i, "_network_interactions_count.pdf")),
    height = 8,
    width = 8,
    family = "serif"
  )
  netVisual_circle_custom(cellchat@net$count,
                          vertex.weight = groupSize,
                          weight.scale = T,
                          label.edge = F,
                          title.name = "Number of interactions"
  )
  dev.off()
  
  
  pdf(
    file = file.path(types3_cellchat_dir, paste0("group_", i, "_network_interactions_weight.pdf")),
    height = 8,
    width = 8,
    family = "serif"
  )
  netVisual_circle_custom(cellchat@net$weight,
                          vertex.weight = groupSize,
                          weight.scale = T,
                          label.edge = F,
                          title.name = "Interaction weights/strength"
  )
  dev.off()
  
  saveRDS(cellchat, file = file.path(types3_cellchat_dir, paste0("group_", i, "_cellchat.rds")))
}




## 循环各分组(神经元、小胶质细胞)------
types2_cellchat_dir <- file.path("RCTD", "cellchat_output_types2/")
dir.create(types2_cellchat_dir)

cellchat_data <- all_data_rctd[, all_data_rctd$Celltype %in% c("Neuron", "Microglial")]
cellchat_data <- JoinLayers(cellchat_data)
cellchat_data$Celltype %>% table()
cellchat_data$Celltype <- as.factor(cellchat_data$Celltype)
cellchat_data$group %>% table()
cellchat_data$group <- as.character(cellchat_data$group)

for (i in unique(cellchat_data$group)) {
  dat <- cellchat_data[, cellchat_data$group == i]
  dat$group %>% table()
  dat$Celltype %>% table()
  dat$samples <- dat$sample
  dat$samples <- factor(dat$samples)
  Idents(dat) <- "Celltype"
  
  
  all_data_meta <- dat[[]]
  all_data_norm <- LayerData(JoinLayers(dat), "data")
  stopifnot(all(rownames(all_data_meta) == colnames(all_data_norm)))
  
  
  cellchat <- createCellChat(
    object = all_data_norm,
    meta = all_data_meta,
    group.by = "Celltype"
  )
  
  
  cellchat <- setIdent(cellchat, ident.use = "Celltype")
  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents))
  
  CellChatDB.use <- CellChatDB.mouse %>%
    subsetDB(
      search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
      key = "annotation"
    )
  
  
  p <- showDatabaseCategory(CellChatDB.use)
  ggsave(file.path(types2_cellchat_dir, "CellChatDB.mmu.png"), width = 10, height = 5, dpi = 300)
  ggsave(file.path(types2_cellchat_dir, "CellChatDB.mmu.pdf"), width = 10, height = 5)
  
  # print(colnames(CellChatDB$interaction))
  # CellChatDB.use <-  # use     Secreted Signaling
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- smoothData(cellchat, adj = PPI.mouse)
  
  ### Part II: Inference of cell-cell communication network
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, population.size = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  ### 信息提取
  df.net <- subsetCommunication(cellchat)
  write.table(df.net,
              file.path(types2_cellchat_dir, paste0("group_", i, "_Communication_info.txt")),
              sep = "\t",
              quote = F,
              row.names = F
  )
  
  ## Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))
  head(cellchat@net$count)
  head(cellchat@net$weight)
  
  
  pdf(
    file = file.path(types2_cellchat_dir, paste0("group_", i, "_network_interactions_count.pdf")),
    height = 8,
    width = 8,
    family = "serif"
  )
  netVisual_circle_custom(cellchat@net$count,
                          vertex.weight = groupSize,
                          weight.scale = T,
                          label.edge = F,
                          title.name = "Number of interactions"
  )
  dev.off()
  
  
  pdf(
    file = file.path(types2_cellchat_dir, paste0("group_", i, "_network_interactions_weight.pdf")),
    height = 8,
    width = 8,
    family = "serif"
  )
  netVisual_circle_custom(cellchat@net$weight,
                          vertex.weight = groupSize,
                          weight.scale = T,
                          label.edge = F,
                          title.name = "Interaction weights/strength"
  )
  dev.off()
  
  saveRDS(cellchat, file = file.path(types2_cellchat_dir, paste0("group_", i, "_cellchat.rds")))
}


## 循环各分组亚类(神经元、小胶质细胞和OPC)------
all_data_rctd$Celltype <- all_data_rctd$first_type

types3_sub_cellchat_dir <- file.path("RCTD", "cellchat_output_types3_sub/")
dir.create(types3_sub_cellchat_dir)

cellchat_data <- all_data_rctd[, all_data_rctd$first_type_combine %in% c("Neuron", "Microglial", "OPC")]
cellchat_data <- JoinLayers(cellchat_data)
cellchat_data$Celltype %>% table()
cellchat_data$Celltype <- as.factor(cellchat_data$Celltype)
cellchat_data$group %>% table()
cellchat_data$group <- as.character(cellchat_data$group)

for (i in unique(cellchat_data$group)) {
  dat <- cellchat_data[, cellchat_data$group == i]
  dat$group %>% table()
  dat$Celltype %>% table()
  dat$samples <- dat$sample
  dat$samples <- factor(dat$samples)
  Idents(dat) <- "Celltype"
  
  
  all_data_meta <- dat[[]]
  all_data_norm <- LayerData(JoinLayers(dat), "data")
  stopifnot(all(rownames(all_data_meta) == colnames(all_data_norm)))
  
  
  cellchat <- createCellChat(
    object = all_data_norm,
    meta = all_data_meta,
    group.by = "Celltype"
  )
  
  
  cellchat <- setIdent(cellchat, ident.use = "Celltype")
  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents))
  
  CellChatDB.use <- CellChatDB.mouse %>%
    subsetDB(
      search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
      key = "annotation"
    )
  
  
  p <- showDatabaseCategory(CellChatDB.use)
  ggsave(file.path(types3_sub_cellchat_dir, "CellChatDB.mmu.png"), width = 10, height = 5, dpi = 300)
  ggsave(file.path(types3_sub_cellchat_dir, "CellChatDB.mmu.pdf"), width = 10, height = 5)
  
  # print(colnames(CellChatDB$interaction))
  # CellChatDB.use <-  # use     Secreted Signaling
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- smoothData(cellchat, adj = PPI.mouse)
  
  ### Part II: Inference of cell-cell communication network
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, population.size = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  ### 信息提取
  df.net <- subsetCommunication(cellchat)
  write.table(df.net,
              file.path(types3_sub_cellchat_dir, paste0("group_", i, "_Communication_info.txt")),
              sep = "\t",
              quote = F,
              row.names = F
  )
  
  ## Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))
  head(cellchat@net$count)
  head(cellchat@net$weight)
  
  
  pdf(
    file = file.path(types3_sub_cellchat_dir, paste0("group_", i, "_network_interactions_count.pdf")),
    height = 8,
    width = 8,
    family = "serif"
  )
  netVisual_circle_custom(cellchat@net$count,
                          vertex.weight = groupSize,
                          weight.scale = T,
                          label.edge = F,
                          title.name = "Number of interactions"
  )
  dev.off()
  
  
  pdf(
    file = file.path(types3_sub_cellchat_dir, paste0("group_", i, "_network_interactions_weight.pdf")),
    height = 8,
    width = 8,
    family = "serif"
  )
  netVisual_circle_custom(cellchat@net$weight,
                          vertex.weight = groupSize,
                          weight.scale = T,
                          label.edge = F,
                          title.name = "Interaction weights/strength"
  )
  dev.off()
  
  saveRDS(cellchat, file = file.path(types3_sub_cellchat_dir, paste0("group_", i, "_cellchat.rds")))
}




## 循环各分组亚类(神经元、小胶质细胞)------
types2_sub_cellchat_dir <- file.path("RCTD", "cellchat_output_types2_sub/")
dir.create(types2_sub_cellchat_dir)

cellchat_data <- all_data_rctd[, all_data_rctd$first_type_combine %in% c("Neuron", "Microglial")]
cellchat_data <- JoinLayers(cellchat_data)
cellchat_data$Celltype %>% table()
cellchat_data$Celltype <- as.factor(cellchat_data$Celltype)
cellchat_data$group %>% table()
cellchat_data$group <- as.character(cellchat_data$group)

for (i in unique(cellchat_data$group)) {
  dat <- cellchat_data[, cellchat_data$group == i]
  dat$group %>% table()
  dat$Celltype %>% table()
  dat$samples <- dat$sample
  dat$samples <- factor(dat$samples)
  Idents(dat) <- "Celltype"
  
  
  all_data_meta <- dat[[]]
  all_data_norm <- LayerData(JoinLayers(dat), "data")
  stopifnot(all(rownames(all_data_meta) == colnames(all_data_norm)))
  
  
  cellchat <- createCellChat(
    object = all_data_norm,
    meta = all_data_meta,
    group.by = "Celltype"
  )
  
  
  cellchat <- setIdent(cellchat, ident.use = "Celltype")
  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents))
  
  CellChatDB.use <- CellChatDB.mouse %>%
    subsetDB(
      search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
      key = "annotation"
    )
  
  
  p <- showDatabaseCategory(CellChatDB.use)
  ggsave(file.path(types2_sub_cellchat_dir, "CellChatDB.mmu.png"), width = 10, height = 5, dpi = 300)
  ggsave(file.path(types2_sub_cellchat_dir, "CellChatDB.mmu.pdf"), width = 10, height = 5)
  
  # print(colnames(CellChatDB$interaction))
  # CellChatDB.use <-  # use     Secreted Signaling
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- smoothData(cellchat, adj = PPI.mouse)
  
  ### Part II: Inference of cell-cell communication network
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, population.size = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  ### 信息提取
  df.net <- subsetCommunication(cellchat)
  write.table(df.net,
              file.path(types2_sub_cellchat_dir, paste0("group_", i, "_Communication_info.txt")),
              sep = "\t",
              quote = F,
              row.names = F
  )
  
  ## Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  groupSize <- as.numeric(table(cellchat@idents))
  head(cellchat@net$count)
  head(cellchat@net$weight)
  
  
  pdf(
    file = file.path(types2_sub_cellchat_dir, paste0("group_", i, "_network_interactions_count.pdf")),
    height = 8,
    width = 8,
    family = "serif"
  )
  netVisual_circle_custom(cellchat@net$count,
                          vertex.weight = groupSize,
                          weight.scale = T,
                          label.edge = F,
                          title.name = "Number of interactions"
  )
  dev.off()
  
  
  pdf(
    file = file.path(types2_sub_cellchat_dir, paste0("group_", i, "_network_interactions_weight.pdf")),
    height = 8,
    width = 8,
    family = "serif"
  )
  netVisual_circle_custom(cellchat@net$weight,
                          vertex.weight = groupSize,
                          weight.scale = T,
                          label.edge = F,
                          title.name = "Interaction weights/strength"
  )
  dev.off()
  
  saveRDS(cellchat, file = file.path(types2_sub_cellchat_dir, paste0("group_", i, "_cellchat.rds")))
}




# ## 逐行跑-----
# library(CellChat)
# library(ggplot2)
# library(ggalluvial)
# library(svglite)
# library(Seurat)
# # library(SeuratData)
# library(patchwork)
# 
# source("result/cellchat_output/netVisual_circle_custom.R")
# 
# cellchat_dir <- file.path("RCTD", "cellchat_output/")
# dir.create(cellchat_dir)
# 
# all_data$annotation %>% table()
# all_data$annotation %>% unique()
# all_data$annotation <- as.character(all_data$annotation)
# all_data@meta.data <- all_data@meta.data %>%
#   mutate(annotation_combine = case_when(
#     grepl("Neuron", annotation) ~ "Neuron",
#     TRUE ~ annotation
#   ))
# 
# all_data$annotation_combine %>% table()
# all_data$annotation_combine %>% unique()
# table(all_data$annotation, all_data$annotation_combine)
# all_data$annotation_combine <- as.character(all_data$annotation_combine)
# 
# ## 循环各分组(神经元、小胶质细胞和OPC)------
# all_data$Celltype <- all_data$annotation_combine
# all_data$group %>% table()
# 
# cellchat_data <- all_data[, all_data$Celltype %in% c("Neuron", "Microglial", "OPC")]
# cellchat_data <- JoinLayers(cellchat_data)
# cellchat_data$Celltype %>% table()
# cellchat_data$Celltype <- as.factor(cellchat_data$Celltype)
# cellchat_data$group %>% table()
# cellchat_data$group <- as.character(cellchat_data$group)
# 
# 
# 
# 
# 
# 
# 
# i <- "Ctrl_E"
# 
# for (i in unique(cellchat_data$group)) {
#   dat <- cellchat_data[, cellchat_data$group == i]
#   dat$group %>% table()
#   dat$Celltype %>% table()
#   dat$samples <- dat$sample
#   dat$samples <- factor(dat$samples)
#   Idents(dat) <- "Celltype"
# 
#   # cellchat <- createCellChat(object = cellchat_data@assays$Spatial.008um@layers$data,meta=cellchat_data@meta.data    ,group.by='Celltype')
#   cellchat <- createCellChat(
#     object = dat,
#     group.by = "Celltype"
#   )
#   cellchat <- setIdent(cellchat, ident.use = "Celltype")
#   levels(cellchat@idents)
#   groupSize <- as.numeric(table(cellchat@idents))
# 
#   CellChatDB <- CellChatDB.mouse
# 
#   p <- showDatabaseCategory(CellChatDB)
#   ggsave(file.path(cellchat_dir, "CellChatDB.mmu.png"), width = 10, height = 5, dpi = 300)
#   ggsave(file.path(cellchat_dir, "CellChatDB.mmu.pdf"), width = 10, height = 5)
# 
#   print(colnames(CellChatDB$interaction))
#   CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")) # use     Secreted Signaling
#   cellchat@DB <- CellChatDB.use
#   cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#   # future::plan("multicore", workers = 1) # do parallel
#   cellchat <- identifyOverExpressedGenes(cellchat)
#   cellchat <- identifyOverExpressedInteractions(cellchat)
# 
#   # cellchat <- projectData(cellchat, PPI.mouse)
#   cellchat <- smoothData(cellchat, adj = PPI.human)
# 
#   ### Part II: Inference of cell-cell communication network
#   cellchat <- computeCommunProb(cellchat)
#   # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#   cellchat <- filterCommunication(cellchat, min.cells = 10)
#   ### 信息提取
#   df.net <- subsetCommunication(cellchat)
#   write.table(df.net,
#     file.path(cellchat_dir, paste0("group_", i, "_Communication_info.txt")),
#     sep = "\t",
#     quote = F,
#     row.names = F
#   )
# 
#   ## Infer the cell-cell communication at a signaling pathway level
#   cellchat <- computeCommunProbPathway(cellchat)
#   cellchat <- aggregateNet(cellchat)
# 
#   groupSize <- as.numeric(table(cellchat@idents))
#   head(cellchat@net$count)
#   head(cellchat@net$weight)
# 
# 
#   pdf(
#     file = file.path(cellchat_dir, paste0("group_", i, "_network_interactions_count.pdf")),
#     height = 9,
#     width = 7,
#     family = "serif"
#   )
#   netVisual_circle_custom(cellchat@net$count,
#     vertex.weight = groupSize,
#     weight.scale = T,
#     label.edge = F,
#     title.name = "Number of interactions"
#   )
#   dev.off()
# 
# 
#   pdf(
#     file = file.path(cellchat_dir, paste0("group_", i, "_network_interactions_weight.pdf")),
#     height = 9,
#     width = 7,
#     family = "serif"
#   )
#   netVisual_circle_custom(cellchat@net$weight,
#     vertex.weight = groupSize,
#     weight.scale = T,
#     label.edge = F,
#     title.name = "Interaction weights/strength"
#   )
# 
#   dev.off()
# 
#   # saveRDS(cellchat, file = file.path(cellchat_dir, paste0("group_", i, "_cellchat.rds")))
# }




all_data_original <- readRDS("R_OUTPUT/all_data.rds")

cellchat_data <- all_data_original[, all_data_original$Celltype %in% c("Neuron", "Microglial", "OPC")]
cellchat_data <- JoinLayers(cellchat_data)
cellchat_data$Celltype %>% table()
cellchat_data$Celltype <- as.factor(cellchat_data$Celltype)
cellchat_data$group %>% table()
cellchat_data$group <- as.character(cellchat_data$group)

for (i in unique(cellchat_data$group)) {
  all_data <- cellchat_data[, cellchat_data$group == i]
  all_data$group %>% table()
  all_data$Celltype %>% table()
  all_data$samples <- all_data$sample
  all_data$samples <- factor(all_data$samples)
  Idents(all_data) <- "Celltype"


  all_data_meta <- all_data[[]]
  all_data_norm <- LayerData(JoinLayers(all_data), "data")
  stopifnot(all(rownames(all_data_meta) == colnames(all_data_norm)))


  cellchat <- createCellChat(
    object = all_data_norm,
    meta = all_data_meta,
    group.by = "Celltype"
  )


  cellchat <- setIdent(cellchat, ident.use = "Celltype")
  levels(cellchat@idents)
  groupSize <- as.numeric(table(cellchat@idents))

  CellChatDB.use <- CellChatDB.mouse %>%
    subsetDB(
      search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"),
      key = "annotation"
    )


  p <- showDatabaseCategory(CellChatDB.use)
  ggsave(file.path(cellchat_dir, "CellChatDB.mmu.png"), width = 10, height = 5, dpi = 300)
  ggsave(file.path(cellchat_dir, "CellChatDB.mmu.pdf"), width = 10, height = 5)

  # print(colnames(CellChatDB$interaction))
  # CellChatDB.use <-  # use     Secreted Signaling
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- smoothData(cellchat, adj = PPI.mouse)

  ### Part II: Inference of cell-cell communication network
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = TRUE, population.size = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  ### 信息提取
  df.net <- subsetCommunication(cellchat)
  write.table(df.net,
    file.path(cellchat_dir, paste0("group_", i, "_Communication_info.txt")),
    sep = "\t",
    quote = F,
    row.names = F
  )

  ## Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  groupSize <- as.numeric(table(cellchat@idents))
  head(cellchat@net$count)
  head(cellchat@net$weight)


  pdf(
    file = file.path(cellchat_dir, paste0("group_", i, "_network_interactions_count.pdf")),
    height = 8,
    width = 8,
    family = "serif"
  )
  netVisual_circle_custom(cellchat@net$count,
    vertex.weight = groupSize,
    weight.scale = T,
    label.edge = F,
    title.name = "Number of interactions"
  )
  dev.off()


  pdf(
    file = file.path(cellchat_dir, paste0("group_", i, "_network_interactions_weight.pdf")),
    height = 8,
    width = 8,
    family = "serif"
  )
  netVisual_circle_custom(cellchat@net$weight,
    vertex.weight = groupSize,
    weight.scale = T,
    label.edge = F,
    title.name = "Interaction weights/strength"
  )
  dev.off()

  saveRDS(cellchat, file = file.path(cellchat_dir, paste0("group_", i, "_cellchat.rds")))
}





