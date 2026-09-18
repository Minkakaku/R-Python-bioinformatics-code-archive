#单细胞数据------
# FIB<- readRDS("customer_RDS/cd4_sp.rds")
# MYE <- readRDS("customer_RDS/cd8_sp.rds")
# T <- readRDS("customer_RDS/mTEC.rds")


sc_data <- readRDS("RCTD/scRNA_0227.rds")


#xenium数据-----
all_data <- readRDS("R_OUTPUT/all_data.rds")

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

# result_dir <- 'result' 
# dir.create(result_dir)


# RCTD ------------------
library(spacexr)
rctd_out_dir <- 'RCTD' 
dir.create(rctd_out_dir)

#1. sketch data  
# NOT RUN. ROI data, only 50k cells
# visium_roi_ls %>%
#   map(function(visium){
#     # sketch the cortical subset of the Visium HD dataset
#     DefaultAssay(visium) <- "Spatial.008um"
#     visium <- FindVariableFeatures(visium)
#     visium <- SketchData(
#       object = cortex,
#       ncells = 50000,
#       method = "LeverageScore",
#       sketched.assay = "sketch"
#     )
#     
#     DefaultAssay(visium) <- "sketch"
#     visium <- ScaleData(visium)
#     visium <- RunPCA(visium, assay = "sketch", reduction.name = "pca.sketch", verbose = T)
#     visium <- FindNeighbors(visium, reduction = "pca.sketch", dims = 1:50)
#     visium <- RunUMAP(visium, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50, verbose = T)
#     
#     visium
#   })

##2. singlecell data处理----
# all_data_rctd <- readRDS("all_data_sc_dim15_anno.rds")
# all_data_rctd <- all_data

sc_data
anno_counts <- table(sc_data@meta.data$Celltypes)
print(anno_counts)

# 保留细胞数 ≥25 的细胞类型
valid_anno <- names(anno_counts[anno_counts >= 25])
print(valid_anno)

# 筛选出属于 valid_anno 的细胞
valid_cells <- sc_data@meta.data$celltype %in% valid_anno
sc_data_filtered <- sc_data[, valid_cells]

# 验证过滤结果
print(table(sc_data_filtered@meta.data$celltype))
dim(sc_data_filtered)


# 从 SCT 提取原始计数矩阵
counts_matrix <- GetAssayData(sc_data_filtered, "RNA","counts")  
# 创建新的 RNA Assay
sc_data[["RNA"]] <- CreateAssayObject(counts = counts_matrix)  
# 设为默认 Assay
DefaultAssay(sc_data_filtered) <- "RNA"  

# all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")

sc_data_filtered$celltype <- sc_data_filtered$celltype 

Idents(sc_data_filtered) <- "celltype"
counts_sc  <- sc_data_filtered
cluster_sc <- as.factor(sc_data_filtered$celltype)
nUMI_sc    <- sc_data_filtered$nCount_RNA

sc_cell_filter_idx <-  cluster_sc != "Unknown"

# create the RCTD reference object
reference <- Reference(
  counts_matrix, 
  droplevels(cluster_sc), 
  nUMI_sc
)


cells_position_ls <- readRDS("~/R_OUTPUT/data_for_analysis/cells_position_ls.rds")

## 3. create the RCTD query object-----
#all_data$sample <- all_data$project
all_data_sp_rctd_ls <- all_data %>% JoinLayers() %>% SplitObject('sample')

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

##4. run RCTD-----
rctd_result_ls <- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 30, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

saveRDS(rctd_result_ls, 'RCTD/rctd_result_ls.rds')

#rctd_result_ls <- readRDS("RCTD/rctd_result_ls.rds")

##5. addMetaData-----
all_data_sp_rctd_anno_ls <-
  all_data_sp_rctd_ls %>%
  imap(function(xenium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    xenium <- AddMetaData(xenium, metadata = rctd@results$results_df)
    
    xenium$first_type <- as.character(xenium$first_type)
    xenium$first_type[is.na(xenium$first_type)] <- "Unknown"
    xenium
  })

saveRDS(all_data_sp_rctd_anno_ls,'RCTD/rctd.rds')

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

all_data_sp_rctd_anno <-
  all_data %>% 
  AddMetaData(meta_rctd_result)


all_data2 <- JoinLayers(all_data_sp_rctd_anno)


xenium <- readRDS("result_new/xenium.rds")
# all_data2_metadata <- all_data2[[]]
xenium_RCTD <- xenium
# xenium_RCTD$cell_id1 <- rownames(xenium_RCTD@meta.data)
# rownames(xenium_RCTD@meta.data) <- paste0(
#   xenium_RCTD@meta.data$orig.ident,"_",xenium_RCTD@meta.data$cell_id)
# 
# xenium_RCTD <- AddMetaData(all_data2_metadata)

all_data2@meta.data$cell_id <- gsub("TMA_","",rownames(all_data2@meta.data))

source_meta <- all_data2@meta.data[, c("cell_id", "first_type")]

xenium_RCTD@meta.data <- merge(
  xenium_RCTD@meta.data,
  source_meta,
  by = "cell_id",
  all.x = TRUE  # 保留所有目标对象细胞
)

##改样本名POST1改成POST2------
xenium_RCTD@meta.data$sample <- gsub("P2_POST1","P2_POST2",xenium_RCTD@meta.data$sample)
xenium_RCTD@meta.data$group1 <- gsub("POST1","POST2",xenium_RCTD@meta.data$group1)

xenium_RCTD@meta.data$cell_id1 <- paste0(xenium_RCTD@meta.data$sample,"_",xenium_RCTD@meta.data$cell_id)
rownames(xenium_RCTD@meta.data) <- xenium_RCTD@meta.data$cell_id1

head(xenium_RCTD@meta.data[, c("cell_id", "first_type")])
Idents(xenium_RCTD) <- "first_type"
xenium_RCTD@meta.data$first_type <- as.character(xenium_RCTD@meta.data$first_type)

##替换NA为other-----
xenium_RCTD@meta.data$first_type %>% unique()
is.na(xenium_RCTD@meta.data$first_type) %>% table()

xenium_RCTD$first_type <- 
  xenium_RCTD$first_type %>%
  as.character() %>% 
  str_replace_na("Others")
all_types <- unique(xenium_RCTD$first_type)

#确保 "Others" 在最后
if("Others" %in% all_types) {
  all_types <- c(setdiff(all_types, "Others"), "Others")
}

xenium_RCTD$first_type <- factor(
  xenium_RCTD$first_type,
  levels = all_types
)

table(xenium_RCTD$first_type, useNA = "ifany")

T1 <- data_frame(xenium_RCTD@meta.data$group,xenium_RCTD@meta.data$group4,xenium_RCTD@meta.data$group1,xenium_RCTD@meta.data$sample)

xenium_RCTD@meta.data$cell_id2 <- NULL

saveRDS(xenium_RCTD,"RCTD/xenium_RCTD.rds")


#RCTD结果绘图----
plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", 
          "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", 
          "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879", 
          "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", 
          "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", "#ffc300", 
          "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

##umap----
pdf(file.path(rctd_out_dir, 'RCTD_umap.pdf'), height = 12, width = 14)
DimPlot(xenium_RCTD, group.by = "first_type", label = F, cols = plan, reduction = 'umap')
dev.off()

pdf(file.path(rctd_out_dir, 'RCTD_group3_umap.pdf'), height = 10, width = 40)
DimPlot(xenium_RCTD, group.by = 'first_type', label = F, split.by = 'group3',cols = plan, reduction = 'umap')
dev.off()

##cell props barplot  --------
xenium_RCTD_meta <-
  xenium_RCTD[[]] #%>%
  #rownames_to_column("cell_id")

p_col_proj <-
  xenium_RCTD_meta %>%
  group_by(sample, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = sample, y = freq, fill = first_type)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj,
  filename = file.path(rctd_out_dir, "cell_freq_project.pdf"),
  height = 6,
  width  = 12
)

p_box_group  <-
  xenium_RCTD_meta %>%
  group_by(group3, sample, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = first_type, y = freq, fill = group3)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) +
  ggsci::scale_fill_igv() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave(
  p_box_group,
  filename = file.path(rctd_out_dir, "cell_freq_group3.pdf"),
  height = 6,
  width  = 20
)


##imagedimplot只画关注细胞-----
# xenium_RCTD@meta.data$first_type1 <- "other"
# xenium_RCTD@meta.data$first_type <- as.character(xenium_RCTD@meta.data$first_type)
# xenium_RCTD@meta.data$first_type1[which(xenium_RCTD@meta.data$first_type %in% c("ASPN_mCAF", "Cancer", "Myeloid"))] <- xenium_RCTD@meta.data$first_type[which(xenium_RCTD@meta.data$first_type %in% c("ASPN_mCAF", "Cancer", "Myeloid"))]
# 
# xenium_RCTD@meta.data$first_type1 %>% unique()
# xenium_RCTD@meta.data$first_type1 <- factor(xenium_RCTD@meta.data$first_type1, levels = c("ASPN_mCAF", "Cancer", "Myeloid","other"))
# 

# all_colors4 <- c("#FF6B6B","#7EBA98","#FDBB69","gray90")
# names(all_colors4) <- c("ASPN_mCAF", "Cancer", "Myeloid","other")
# sample_label <- unique(xenium@meta.data$sample)
# 
# for( i in sample_label){
#   #j = gsub("_","-",i)
#   DefaultBoundary(xenium_RCTD@images[[i]]) <- "segmentation"
#   p <- ImageDimPlot(xenium_RCTD,
#                     fov = i,
#                     nmols = 20000,
#                     dark.background = F,
#                     border.color = "black", 
#                     border.size = 0.1,
#                     cols = all_colors4,
#                     group.by = "first_type1"
#   )+
#     ggtitle(paste(i,"_annotation")) +
#     theme(plot.title = element_text(
#       size = 20, 
#       face = "bold", #加粗
#       hjust = 0.5, #居中
#       vjust = -3 #下移
#     ))#+
#   #coord_flip()+ # 交换X/Y轴
#   #scale_x_reverse()  # 反转x轴方向（匹配顺时针效果）
#   
#   ggsave(p,filename = paste0("result/RCTD/Imagedimplot/",i,"_Epithelial_annotation.png"),width = 12 ,height = 12,dpi = 300)
# }

##imagedimplot所有细胞-----
colors47 <- c(#"#CD5C5C",
  # "#FA8072",
  # "#FFA07A",
  "#DC143C",
  #"#B22222",
  "#FFC0CB",
  "#FF69B4",
  #"#FF1493",
  "#C71585",
  #"#DB7093",
  "#FFA07A",
  #"#FF7F50",
  "#FF6347",
  "#FF4500",
  "#FF8C00",
  "#FFD700",
  "#FFFF00",
  "#F0E68C",
  "#BDB76B",
  "#E6E6FA",
  "#EE82EE",
  "#FF00FF",
  "#9370DB",
  "#800080",
  "#ADFF2F",
  "#32CD32",
  "#00FA9A",
  "#3CB371",
  "#008000",
  "#808000",
  "#20B2AA",
  "#008080",
  "#00FFFF",
  "#AFEEEE",
  "#40E0D0",
  "#4682B4",
  "#ADD8E6",
  "#00BFFF",
  "#7B68EE",
  "#00008B",
  "#FFE4C4",
  "#DEB887",
  "#BC8F8F",
  "#D2691E",
  "#A52A2A",
  "#DCDCDC",
  "#808080",
  "#708090"
)

sample_label <- unique(xenium_RCTD@meta.data$sample)

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium_RCTD@images[[i]]) <- "segmentation"
  # i = "P2-POST1"
  # xenium_RCTD@images[["P2-POST1"]] <- xenium_RCTD@images[["P2_POST1"]] 
  p <- ImageDimPlot(xenium_RCTD,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "black", 
                    border.size = 0.1,
                    cols = colors47,
                    group.by = "first_type"
                    #shuffle.cols = F
  )+
    ggtitle(paste(i,"_annotation")) +
    theme(plot.title = element_text(
      size = 20, 
      face = "bold", #加粗
      hjust = 0.5, #居中
      vjust = -3 #下移
    ))#+
  #coord_flip()+ # 交换X/Y轴
  #scale_x_reverse()  # 反转x轴方向（匹配顺时针效果）
  
  ggsave(p,filename = paste0("RCTD/Imagedimplot/",i,"_annotation.pdf"),width = 12 ,height = 12,dpi = 300)
}



##SegImagefeature----
SegImagefeature_dir <- 'result/RCTD/SegImagefeature_ASPN'
dir.create(SegImagefeature_dir)

gene_list <- c('ASPN')
valid_gene <- gene_list[gene_list %in% rownames(xenium)]
cat("有效基因:", length(valid_gene), "/", length(gene_list))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(gene_list, valid_gene)

gene_cus <- valid_gene

for (sample_id in unique(xenium@meta.data$sample)) {
  # 智能匹配样本ID
  j <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  
  if (length(j) == 0) {
    warning("跳过未匹配样本: ", sample_id)
    next
  }
  
  for (gene in gene_cus) {
    tryCatch({
      # 空值检查
      if (!is.null(xenium@images[[j]])) {
        DefaultBoundary(xenium@images[[j]]) <- "segmentation"
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
      
      message("成功保存: ", sample_id, "_", gene)
      
    }, error = function(e) {
      warning("失败[", sample_id, ":", gene, "]: ", e$message)
    })
  }
}


###normalize ------
all_data2 <- NormalizeData(all_data2, scale.factor = median(all_data2$nCount_RNA))
all_data2 <- FindVariableFeatures(all_data2)

###dimension reduction - pca ---------
dr_dir <- "dimension reduction" 
dir.create(dr_dir)

DefaultAssay(all_data2) <- "RNA"

# all_data <- NormalizeData(all_data, scale.factor = median(all_data$nCount_RNA))

all_data2 <- FindVariableFeatures(all_data2, nfeatures = 1500)
VariableFeaturePlot(all_data2)

all_data2 <- ScaleData(all_data2, vars.to.regress = c('nCount_RNA'))
all_data2 <- RunPCA(all_data2)

pdf(file.path(dr_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data2, ndims = 30)
dev.off()


set.seed(1234)
suset_cells <- 
  all_data2[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data2) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(dr_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 1:5)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(dr_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 16:20)
dev.off()


#####dimension reduction --------
dims <- 16


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


all_data2 <- RunUMAP(all_data2, dims = 1:dims, return.model = TRUE, n.neighbors = 30)


pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 8, width = 8)
DimPlot(all_data2, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


all_data2 <- all_data2 %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data2 <- RunUMAP(
  all_data2,
  dims = 1:dims,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_", 
  return.model = TRUE
)

pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data2, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'


pdf(file.path(dr_dir, str_glue("Dimplot dim{dims} annotation {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data2, group.by = "first_type", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# pdf(file.path(dr_dir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
# DotPlot(all_data_rctd, group.by = "annotation", features = major_marker1)+ RotatedAxis()
# dev.off()

saveRDS(all_data2,'RCTD/rctd_umap_anno.rds')
all_data2 <- readRDS("~/RCTD/rctd_umap_anno.rds")

# cells_boundary_ls <- readRDS(file.path("~/data_for_analysis", "cells_boundary_ls.rds"))
# cells_position_ls <- readRDS(file.path("~/data_for_analysis", "cells_position_ls.rds"))
# image_wh_r_ls <- readRDS(file.path("~/data_for_analysis", "image_wh_r_ls.rds"))

all_data2@meta.data %>% head()
#all_data2@meta.data$cell_id <- gsub("1_","",rownames(all_data2@meta.data))
# 使用正则表达式移除所有行名前缀
all_data2@meta.data$cell_id <- gsub("^[0-9]+_", "", rownames(all_data2@meta.data))
# 提取所需的列并重命名
meta_subset <- all_data2@meta.data[, c("cell_id", "first_type")]
colnames(meta_subset) <- c("cell_id", "group")
# 获取样本列表（根据orig.ident列）
samples <- unique(all_data2@meta.data$orig.ident)
# 为每个样本创建CSV文件
for(sample_id in samples) {
  # 筛选当前样本的数据
  sample_data <- meta_subset[all_data2@meta.data$orig.ident == sample_id, ]
  
  # 保存为CSV文件
  write.csv(sample_data, 
            file = paste0(sample_id, ".csv"), 
            row.names = FALSE)
  
  # 输出确认信息
  message(paste("已创建", sample_id, ".csv (包含", nrow(sample_data), "个细胞)"))
}


plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", 
          "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", 
          "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879", 
          "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", 
          "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", "#ffc300", 
          "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")

#RCTD结果绘图----
###umap----
pdf(file.path(rctd_out_dir, 'RCTD_umap.pdf'), height = 12, width = 14)
DimPlot(xenium_RCTD, group.by = 'first_type', label = F, cols = plan, reduction = 'umap') 
dev.off()

pdf(file.path(rctd_out_dir, 'RCTD_group_umap.pdf'), height = 10, width = 18)
DimPlot(xenium_RCTD, group.by = 'first_type', label = F, split.by = 'group',cols = plan, reduction = 'umap')
dev.off()

###cell props barplot  --------
all_data_meta <-
  all_data2[[]] %>%
  rownames_to_column("cell_id")


p_col_proj <-
  all_data_meta %>%
  group_by(project, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = first_type)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


p_box_group  <-
  all_data_meta %>%
  group_by(group, project, first_type) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = first_type, y = freq, fill = group)) +
  geom_boxplot(outlier.colour = NA) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) +
  ggsci::scale_fill_igv() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


# export
ggsave(
  p_col_proj,
  filename = file.path(rctd_out_dir, "cell_freq_project.pdf"),
  height = 6,
  width  = 6
)

ggsave(
  p_box_group,
  filename = file.path(rctd_out_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 10
)

## spatial plot(annotation) -----
###将对象转为样本级列表
all_data2_list <- Seurat::SplitObject(all_data2, split.by = "project")
sp_p_ls_subset <-
  cells_boundary_ls %>%
  .[names(all_data2_list)] %>%
  imap(function(cells_boundary, project){  #browser()
    image_wh_r <- image_wh_r_ls[[project]]
    all_data_meta <- all_data2_list[[project]][[]] %>% rownames_to_column('cell_id')
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = first_type, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      # scale_fill_gradient(low = "gray90", high = "red") +
      theme_void() +
      labs(
        title = project
      )
    p
  })

# export
sp_p_ls_subset %>% 
  iwalk(function(p, grp){
    pdf(file.path(result_dir, str_glue("spatial plot annotation {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })



## spatial plot(probes) -----
gene_probes <- list(
  "mTEC" = c("Tbata", "Gpbar1", "Shp", "Nr1h4", "Ahr"),
  "SP4-0" = c("Ccr9", "Itm2a", "Ldhb"),
  "SP4-1" = c("Klf2", "S1pr1", "Ccr7"),
  "Sp8-3" = c("Klf2", "Ccr7"),
  "Sp8-6" = c("Bst2", "Isg15")
)

genes <- c(
  "Tbata", "Gpbar1", "Shp", "Nr1h4", "Ahr",
  "Ccr9", "Itm2a", "Ldhb",
  "Klf2", "S1pr1", "Ccr7",
  "Klf2", "Ccr7",
  "Bst2", "Isg15"
)
valid_genes <- genes[genes %in% rownames(xenium)]
cat("有效基因:", length(valid_genes), "/", length(genes))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_genes <- setdiff(genes, valid_genes)

# 定义基因探针列表（基于图片中的表格）
gene_probes <- list(
  mTEC = c("Gpbar1", "Nr1h4", "Ahr"),
  `SP4-0` = c("Ccr9", "Ldhb"),
  `SP4-1` = c("S1pr1", "Ccr7"),
  `Sp8-3` = c("Ccr7"),
  `Sp8-6` = c("Isg15")
)
# 提取所有唯一基因名
all_genes <- unique(unlist(gene_probes))

project_dir <- c("output-XETG00305__0041785__K01__20241226__052857",
                 "output-XETG00305__0041785__K02__20241226__052857",
                 "output-XETG00305__0041785__WT1__20241226__052858",
                 "output-XETG00305__0041785__WT3__20241226__052858"
                 )

# 修改绘图函数
sp_p_ls_subset <-
  cells_boundary_ls %>%
  .[names(all_data2_list)] %>%
  imap(function(cells_boundary, project){  #browser()
    image_wh_r <- image_wh_r_ls[[project]]
    all_data_meta <- all_data2_list[[project]][[]] %>% rownames_to_column('cell_id')
    
    # 创建基础数据框
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    # 获取当前样本的Seurat对象
    sample_obj <- all_data2_list[[project]]
    
    # 创建一个空列表存储所有图
    all_plots <- list()
    
    # # 1. 添加细胞类型注释图
    # p_cell <- dat %>%
    #   ggplot(aes(x = vertex_x, y = vertex_y)) +
    #   geom_polygon(aes(fill = first_type, group = cell_id), alpha = 0.8)  +
    #   scale_y_reverse() +
    #   coord_fixed(ratio = 1/image_wh_r) +
    #   ggsci::scale_fill_igv() +
    #   theme_void() +
    #   labs(title = str_glue("{project} - 细胞类型"))
    # all_plots[["Cell Types"]] <- p_cell
    # 
    # # 2. 添加探针基因组平均表达热图
    # for (probe_name in names(gene_probes)) {
    #   genes <- gene_probes[[probe_name]]
    # 
    #   # 计算探针基因组的平均表达量
    #   probe_expr <- FetchData(sample_obj, vars = genes) %>%
    #     rowMeans(na.rm = TRUE)
    # 
    #   # 创建表达热图
    #   p_probe <- dat %>%
    #     mutate(probe_expr = probe_expr[cell_id]) %>%
    #     ggplot(aes(x = vertex_x, y = vertex_y)) +
    #     geom_polygon(aes(fill = probe_expr, group = cell_id)) +
    #     scale_y_reverse() +
    #     coord_fixed(ratio = 1/image_wh_r) +
    #     scale_fill_gradientn(
    #       colors = c("navy", "blue", "cyan", "yellow", "red"),
    #       name = "平均表达量"
    #     ) +
    #     theme_void() +
    #     labs(title = str_glue("{project} - {probe_name}表达热图"))
    # 
    #   all_plots[[probe_name]] <- p_probe
    # }

    # # 3. 添加基因分子点位置图
    # for (gene in all_genes) {
    #   # 尝试寻找分子数据文件
    #   mol_file_found <- FALSE
    #   
    #   for(i in project_dir){
    #     mol_file <- file.path(data_dir, i, "transcripts.parquet")
    #     
    #     if (file.exists(mol_file)) {
    #       mol_file_found <- TRUE
    #       
    #       # 使用arrow包读取Parquet文件
    #       require(arrow)
    #       tryCatch({
    #         molecules <- read_parquet(mol_file) %>%
    #           dplyr::filter(gene == !!gene) %>%
    #           dplyr::sample_n(min(10000, n()))  # 限制点数量
    #         
    #         # 创建分子点图
    #         p_mol <- dat %>%
    #           ggplot(aes(x = vertex_x, y = vertex_y)) +
    #           geom_polygon(aes(group = cell_id), fill = "gray90", alpha = 0.3) +
    #           geom_point(
    #             data = molecules,
    #             aes(x = x_location, y = y_location),
    #             color = "red", size = 0.3, alpha = 0.5
    #           ) +
    #           scale_y_reverse() +
    #           coord_fixed(ratio = 1/image_wh_r) +
    #           theme_void() +
    #           labs(title = str_glue("{project} - {gene}分子位置"))
    #         
    #         all_plots[[gene]] <- p_mol
    #       }, error = function(e) {
    #         warning("读取文件失败: ", mol_file, " - ", e$message)
    #       })
    #       break  # 找到文件后退出循环
    #     }
    #   }
    #   
    #   if (!mol_file_found) {
    #     warning("未找到分子文件: 检查了 ", paste(project_dir, collapse = ", "))
    #   }
    # }
    
    
    # 3. 添加基因分子点位置图（包含细胞分割边界）
    for (gene in all_genes) {
      # 尝试寻找分子数据文件
      mol_file_found <- FALSE
      
      for(i in project_dir){
        mol_file <- file.path(data_dir, i, "transcripts.parquet")
        
        if (file.exists(mol_file)) {
          mol_file_found <- TRUE
          
          # 使用arrow包读取Parquet文件
          require(arrow)
          tryCatch({
            molecules <- read_parquet(mol_file) %>%
              dplyr::filter(gene == !!gene) %>%
              dplyr::sample_n(min(10000, n()))  # 限制点数量
            
            # 创建分子点图 - 添加细胞边界线
            p_mol <- ggplot() +
              # 绘制细胞分割边界（黑色轮廓）
              geom_path(
                data = cells_boundary,
                aes(x = vertex_x, y = vertex_y, group = cell_id),
                color = "black", size = 0.1, alpha = 0.5
              ) +
              # 绘制细胞内部半透明填充
              geom_polygon(
                data = dat,
                aes(x = vertex_x, y = vertex_y, group = cell_id),
                fill = "white", alpha = 0.2
              ) +
              # 绘制基因分子点
              geom_point(
                data = molecules,
                aes(x = x_location, y = y_location),
                color = "red", size = 0.01, alpha = 0.7
              ) +
              scale_y_reverse() +
              coord_fixed(ratio = 1/image_wh_r) +
              theme_void() +
              labs(title = str_glue("{project} - {gene}分子位置 (含细胞边界)"))
            
            all_plots[[gene]] <- p_mol
          }, error = function(e) {
            warning("读取文件失败: ", mol_file, " - ", e$message)
          })
          break  # 找到文件后退出循环
        }
      }
      
      if (!mol_file_found) {
        warning("未找到分子文件: 检查了 ", paste(project_dir, collapse = ", "))
      }
    }
    
    # 返回所有图的列表
    return(all_plots)
  })

# 创建嵌套输出目录
probe_dir <- file.path(result_dir, "gene_probes")
dir.create(probe_dir, showWarnings = FALSE, recursive = TRUE)

# 输出所有图像
sp_p_ls_subset %>%
  iwalk(function(plot_list, project) {
    # 1. 单独输出每个图
    for (plot_name in names(plot_list)) {
      clean_name <- gsub("/", "_", plot_name)  # 处理特殊字符
      pdf_file <- file.path(probe_dir, str_glue("spatial_{project}_{clean_name}.pdf"))
      pdf(pdf_file, width = 12, height = 10)
      print(plot_list[[plot_name]])
      dev.off()
      
      # png_file <- file.path(probe_dir, str_glue("spatial_{project}_{clean_name}.png"))
      # ggsave(png_file, plot_list[[plot_name]], width = 12, height = 10, dpi = 300)
    }
    
    # # 2. 输出组合对比图（每行3个）
    # pdf(file.path(probe_dir, str_glue("spatial_{project}_combined.pdf")), 
    #     width = 18, height = 6 * ceiling(length(plot_list)/3))
    # print(wrap_plots(plot_list, ncol = 3))
    # dev.off()
  })







#加入FOV----
library(Seurat)
library(SeuratDisk)
library(sceasy)
library(reticulate)
library(dplyr)
library(ggplot2)
library(argparser)
library(patchwork)
library(RColorBrewer)
library(arrow)

source("~/h5ad2rds(1).R")

all_data2@assays$RNA@layers %>% class


all_data3 <- all_data2
###改行名
# 使用正则表达式移除所有行名前缀
all_data3@meta.data$cell_id <- gsub("^[0-9]+_", "", rownames(all_data3@meta.data))
all_data3@meta.data$sample <- all_data3@meta.data$project
rownames(all_data3@meta.data) <- paste0(all_data3@meta.data$project,"_",all_data3@meta.data$cell_id)

#all_data2 <- JoinLayers(all_data2)
rownames(all_data3@assays$RNA@layers$counts) <- rownames(all_data3)
colnames(all_data3@assays$RNA@layers$counts) <- colnames(all_data3)
xenium <- CreateSeuratObject(
  all_data3@assays$RNA@layers$counts,
  assay = "Xenium")

xenium@meta.data <- all_data3@meta.data

xenium@meta.data$cell_id1 <- rownames(xenium@meta.data)

umap <- all_data3@reductions$umap@cell.embeddings
rownames(umap) <- rownames(xenium@meta.data)
metadata <- xenium@meta.data 
xenium <- Add_meta_umapData_Xenium(object = xenium,umapfile = umap,metedata = metadata)
colnames(xenium@assays$Xenium) <-  colnames(xenium)
row.names(xenium@reductions$umap@cell.embeddings) <- colnames(xenium)

xenium <- RenameCells(xenium,new.names = colnames(xenium))
label <- unique(xenium@meta.data$project)

colnames(xenium) #%>% head()
rownames(xenium) #%>% head()
rownames(xenium@assays$Xenium@layers$counts) #%>% head()
colnames(xenium@assays$Xenium@layers$counts) #%>% head()
rownames(xenium@assays$Xenium@layers$counts) <- rownames(xenium)
colnames(xenium@assays$Xenium@layers$counts) <- colnames(xenium)

xenium <- NormalizeData(xenium)

###加入1个FOV测试----
library(arrow)
coord1 <- Read_fov_spatial(data.dir = "/data/xenium/20241226__052526__FWSP20242241_1491/output-XETG00305__0041785__K01__20241226__052857",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "K01",object = xenium,cov = NULL)
xenium@images[["K01"]] = coord1

options(future.globals.maxSize = 10 * 1024^3)
#xenium <- NormalizeData(xenium)

p <- ImageFeaturePlot(xenium, fov = "K01",features = c("Cd3d"),max.cutoff = "q90",axes = FALSE,
                      cols = c('#003366', '#FFD700'))
ggsave(p,filename = paste0("_feature_spatial.png"),width = 12 ,height = 12,dpi = 300)


sample_label <- unique(xenium@meta.data$sample)
for( i in sample_label){
  j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[j]]) <- "segmentation"
  p <- ImageDimPlot(xenium,fov = j,nmols = 20000,border.color = "#003366", border.size = 0.1,
                    cols = c('#003366', "#E7475E"), group.by = "first_type")
  ggsave(p,filename = paste0(i,"_spatial.png"),width = 12 ,height = 12,dpi = 300)
}

#saveRDS(xenium,paste0('xenium_coor',"_final2.rds")) 

# data_path <-"/data/xenium/20241226__052526__FWSP20242241_1491/output-XETG00305__0041785__K01__20241226__052857"
# dir.exists(data_path) 
# required_files <- c("cells.csv.gz", "transcripts.csv.gz", "cell_polygons.csv.gz")
# all(file.exists(file.path(data_path, required_files)))

###加入4个样本FOV(正式)----

coord1 <- Read_fov_spatial(data.dir = "/data/xenium/20241226__052526__FWSP20242241_1491/output-XETG00305__0041785__K01__20241226__052857",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "K01",object = xenium,cov = NULL)
xenium@images[["K01"]] = coord1

coord2 <- Read_fov_spatial(data.dir = "/data/xenium/20241226__052526__FWSP20242241_1491/output-XETG00305__0041785__K02__20241226__052857",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "K02",object = xenium,cov = NULL)
xenium@images[["K02"]] = coord2

coord3 <- Read_fov_spatial(data.dir = "/data/xenium/20241226__052526__FWSP20242241_1491/output-XETG00305__0041785__WT1__20241226__052858",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "WT1",object = xenium,cov = NULL)
xenium@images[["WT1"]] = coord3

coord4 <- Read_fov_spatial(data.dir = "/data/xenium/20241226__052526__FWSP20242241_1491/output-XETG00305__0041785__WT3__20241226__052858",assay = 'Xenium',
                           mols.qv.threshold = 20,label = "WT3",object = xenium,cov = NULL)
xenium@images[["WT3"]] = coord4

options(future.globals.maxSize = 10 * 1024^3)
#xenium <- NormalizeData(xenium)

#saveRDS(xenium,paste0('xenium_fov',"_anno.rds")) 
xenium <- readRDS("~/RCTD/xenium_fov_anno.rds")

immune_genes <- c(
  "Cd3d", "Cd3e", "Cd4", "Cd8a"
  # "CD3E", "CD4", "CD68", "CD69", "CD74", "CD8A", 
  # "CD163", "CD274", "CIITA", "CTLA4", "CXCL9", "CXCL10", "CXCL13", "FOXP3", 
  # "GZMA", "GZMB", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "ICOS", 
  # "IDO1", "IFNG", "IL1B", "IL6", "IL10", "IL12B", "LAG3", "NOS2", "PDCD1", 
  # "PRF1", "TGFB1", "TIGIT", "TNF", "TNFRSF9", "TOX"
)
valid_immune_genes <- immune_genes[immune_genes %in% rownames(xenium)]
cat("有效基因:", length(valid_immune_genes), "/", length(immune_genes))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_gene <- setdiff(immune_genes, valid_immune_genes)

gene_cus <- valid_immune_genes

###ImageFeaturePlot----
SegImagefeature_dir <- 'result/SegImagefeature' 
dir.create(SegImagefeature_dir)

for (sample_id in unique(xenium@meta.data$sample)) {
  # 智能匹配样本ID
  j <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  
  if (length(j) == 0) {
    warning("跳过未匹配样本: ", sample_id)
    next
  }
  
  for (gene in gene_cus) {
    tryCatch({
      # 空值检查
      if (!is.null(xenium@images[[j]])) {
        DefaultBoundary(xenium@images[[j]]) <- "segmentation"
      }
      
      p <- ImageFeaturePlot(
        object = xenium,
        fov = j,
        nmols = 20000,       # 每个图像中显示的最大分子数（随机子抽样）
        border.color = "white",  # 细胞边界的颜色（深蓝色）
        border.size = 0,    # 细胞边界的线宽
        features = gene,
        #group.by = "first_type",
        max.cutoff = "q90",
        axes = FALSE,
        cols = c('#003366', "#E7475E"),
        #size = 0.8,
        dark.background = F
      )
      
      # ggsave(
      #   filename = file.path(SegImagefeature_dir, paste0(sample_id, "_", gene, "_spatial.png")),
      #   plot = p,
      #   width = 12,
      #   height = 12,
      #   dpi = 300
      # )
      
      ggsave(
        filename = file.path(SegImagefeature_dir, paste0(sample_id, "_", gene, "_spatial.pdf")),
        plot = p,
        width = 12,
        height = 12
        #dpi = 300
      )
      
      message("成功保存: ", sample_id, "_", gene)
      
    }, error = function(e) {
      warning("失败[", sample_id, ":", gene, "]: ", e$message)
    })
  }
}


###ImageDimPlot----
ImagedimPlot_dir <- 'result/ImagedimPlot' 
dir.create(ImagedimPlot_dir)

sample_label <- unique(xenium@meta.data$sample)

for( i in sample_label){
  #j = gsub("_","-",i)
  DefaultBoundary(xenium@images[[i]]) <- "segmentation"
  p <- ImageDimPlot(xenium,
                    fov = i,
                    nmols = 20000,
                    dark.background = F,
                    border.color = "white", 
                    border.size = 0.05,
                    cols = plan,
                    #shuffle.cols = T, #随机打乱颜色顺序
                    #size = 1.2,
                    group.by = "first_type",
                    #molecules = c("INHBA"),
                    #mols.size = 1,
                    #mols.cols = "red",
                    #mols.alpha = 0.5, #分子的透明度值
                    
  )
  #ggsave(p,filename = paste0("result/Imagedimplot/",i,"_spatial.png"),width = 12 ,height = 12,dpi = 300)
  ggsave(p,filename = paste0("result/ImagedimPlot/",i,"_spatial.pdf"),width = 12 ,height = 12)
}


###各亚群单独空间映射图----
# 创建按亚群分组的输出目录
by_group_dir <- file.path(ImagedimPlot_dir, "spatial_by_first_type")
dir.create(by_group_dir)

# 获取所有亚群类型
#all_types <- na.omit(unique(xenium@meta.data$first_type))
all_types <- c("mTEC","CD4_C0","CD4_C1","CD8_C3","CD8_C6")

# 循环处理每个样本
for(sample_id in sample_label) {
  # 匹配样本对应的图像名称
  img_name <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  if(length(img_name) == 0) next
  
  # 设置当前样本的边界类型
  DefaultBoundary(xenium@images[[img_name]]) <- "segmentation"
  
  # 循环绘制每个亚群
  for(cell in all_types) {
    tryCatch({
      # 创建临时元数据列标记目标细胞
      xenium$celltype <- ifelse(
        xenium$first_type == cell & xenium$sample == sample_id,
        cell,
        "Other"
      )
      
      # 动态颜色映射
      color_map <- c("#1F77B4", "gray")
      names(color_map) <- c(cell, "Other")  # 键名=实际分组值
      
      # 确保分组为因子
      xenium$celltype <- factor(
        xenium$celltype, 
        levels = names(color_map)  # 与颜色向量顺序一致
      )
      
      # 核心绘图函数（使用group.by替代cells.highlight）
      p <- ImageDimPlot(
        xenium,
        fov = img_name,
        group.by = "celltype",  # 使用新创建的元数据列
        cols = color_map, 
        dark.background = F,
        border.color = "white",
        border.size = 0.05,
        size = 0.8
      ) +
        ggtitle(paste(sample_id, "-", cell)) +
        theme(
          plot.title = element_text(size=14, face="bold", hjust = 0.5),
          legend.position = "right"
        )
      
      # 文件名安全处理
      safe_name <- gsub("[/\\]", "_", cell) %>% gsub(" ", "_", .)
      
      # 保存高分辨率图像
      ggsave(
        filename = file.path(by_group_dir, paste0(sample_id, "_", safe_name, ".pdf")),
        plot = p,
        width = 10,
        height = 10
        #,
        #dpi = 300
      )
      message("成功保存: ", sample_id, " | ", cell)
      
      # 移除临时列
      xenium$celltype <- NULL
    }, error = function(e) {
      warning("失败[", sample_id, ":", cell, "]: ", e$message)
    })
  }
}

xenium$celltype <- NULL


###各亚群和mTEC共现空间映射图-----
# 创建按亚群分组的输出目录
by_group_dir <- file.path(ImagedimPlot_dir, "spatial_by_first_type_with_mTEC")
dir.create(by_group_dir, showWarnings = FALSE, recursive = TRUE)

# 获取所有亚群类型
#all_types <- na.omit(unique(xenium@meta.data$first_type))
all_types <- c("mTEC","CD4_C0","CD4_C1","CD8_C3","CD8_C6")

# 循环处理每个样本
for(sample_id in sample_label) {
  # 匹配样本对应的图像名称
  img_name <- names(xenium@images)[grepl(sample_id, names(xenium@images))]
  if(length(img_name) == 0) next
  
  # 设置当前样本的边界类型
  DefaultBoundary(xenium@images[[img_name]]) <- "segmentation"
  
  # 循环绘制每个亚群
  for(cell in all_types) { 
    tryCatch({
      # 创建临时元数据列标记目标细胞
      xenium$celltype <- "Other"  # 初始化所有细胞为Other
      
      # 标记当前亚群细胞 (蓝色)
      idx_target <- which(
        xenium$first_type == cell & 
          xenium$sample == sample_id
      )
      xenium$celltype[idx_target] <- cell
      
      # 标记mTEC细胞 (红色)
      idx_mtec <- which(
        xenium$first_type == "mTEC" & 
          xenium$sample == sample_id
      )
      xenium$celltype[idx_mtec] <- "mTEC"
      
      # 动态颜色映射
      color_map <- c("#1F77B4", "#FF7F0E", "gray")
      names(color_map) <- c(cell, "mTEC", "Other")  # 键名=实际分组值
      
      # 确保分组为因子
      xenium$celltype <- factor(
        xenium$celltype, 
        levels = names(color_map)  # 与颜色向量顺序一致
      )
      
      # 核心绘图函数
      p <- ImageDimPlot(
        xenium,
        fov = img_name,
        group.by = "celltype",
        cols = color_map,
        dark.background = FALSE,
        border.color = "white",
        border.size = 0.05,
        size = 0.8
      ) +
        ggtitle(paste(sample_id, "-", cell, "with mTEC")) +
        theme(
          plot.title = element_text(
            size = 14, 
            face = "bold", 
            hjust = 0.5
            ),
          legend.position = "right"
        )
        
      
      # 文件名安全处理
      safe_name <- gsub("[/\\]", "_", cell) %>% gsub(" ", "_", .)
      
      # 保存高分辨率PDF
      ggsave(
        filename = file.path(by_group_dir, paste0(sample_id, "_", safe_name, "_with_mTEC.pdf")),
        plot = p,
        width = 10,
        height = 10
      )
      message("成功保存: ", sample_id, " | ", cell, " with mTEC")
      
      # 移除临时列
      xenium$celltype <- NULL
    }, error = function(e) {
      warning("失败[", sample_id, ":", cell, "]: ", e$message)
    })
  }
}




###放大一个区域(xenium explorer获取坐标)------
options(future.globals.maxSize = 10 * 1024^3)
crop1 <- Crop(xenium[["K01"]], x = c(1366.4283131414966,
                                     1632.2729147238167
                                     ), 
                               y = c(1457.1216307457019,
                                     1756.848387431651
                                     ))

xenium@images[["zoom1"]] <- crop1
DefaultBoundary(xenium@images[["zoom1"]]) <- "segmentation"

p <- ImageDimPlot(xenium, fov = "zoom1", 
                  cols = plan, 
                  nmols = 20000,
                  dark.background = F,
                  border.color = "white", 
                  border.size = 0.05,
                  group.by = "first_type",
                  coord.fixed = FALSE,
                  flip_xy = FALSE
                  )
                  #+coord_flip()  # 交换X/Y轴
                  #+scale_x_reverse()  # 反转x轴方向（匹配顺时针效果）

p
ggsave(p,filename = "result/ImagedimPlot/zoom1_spatial.pdf",width = 12 ,height = 12)



genes <- c(
  "Tbata", "Gpbar1", "Shp", "Nr1h4", "Ahr",
  "Ccr9", "Itm2a", "Ldhb",
  "Klf2", "S1pr1", "Ccr7",
  "Klf2", "Ccr7",
  "Bst2", "Isg15"
)
valid_genes <- genes[genes %in% rownames(xenium)]
cat("有效基因:", length(valid_genes), "/", length(genes))

# 提取差集：virus_genes中存在但valid_Vriusgenes中缺失的基因
missing_genes <- setdiff(genes, valid_genes)

# note the clouds of TPSAB1 molecules denoting mast cells
p <- ImageDimPlot(xenium, 
                  fov = "zoom1", 
                  cols = plan, 
                  group.by = "first_type",
                  alpha = 0.6, 
                  molecules = valid_genes, 
                  mols.size = 0.6, 
                  nmols = 20000, 
                  border.color = "white", 
                  dark.background = F,
                  coord.fixed = FALSE
                )+
  ggtitle(paste("k01", "-", "zoom1", "gene")) +
  theme(plot.title = element_text(
                     size = 30, 
                     face = "bold", #加粗
                     hjust = 0.5, #居中
                     vjust = -3 #下移
                     ))+
  theme(legend.position = "bottom" #,  # 图例置底
        #legend.title = element_blank() # 隐藏图例标题
        )+
  guides(color = guide_legend(ncol = 3))  # 图例分2列显示

              

ggsave(p,filename = "result/ImagedimPlot/k01_zoom1_mole_spatial.pdf",width = 12 ,height = 12)



