library(tidyverse)
library(Seurat)
library(patchwork)

# setting dir -------
img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)
mk_output <- "4. marker gene" # dir.create(mk_output)
er_output <- "5. enriment analysis" # dir.create(er_output)

dir.create(img_output)
dir.create(qc_outdir)
dir.create(dr_outdir)
dir.create(sp_output)
dir.create(mk_output)
dir.create(er_output)


list.files("/data/xenium/XE2601270675-01")

# library --------------
# copy jar file.
bf_cache_dir <- tools::R_user_dir("RBioFormats", which = "cache")
dir.create(bf_cache_dir)
bf_jar_path <- list.files(system.file("java", package = "RBioFormats"), pattern = "jar", full.names = TRUE)
bf_ver <- read.dcf(file.path(system.file(package = "RBioFormats"), "DESCRIPTION"), "BioFormats")
bf_cache_name <- sprintf("bioformats_package_%s.jar", bf_ver)
file.copy(bf_jar_path, file.path(bf_cache_dir, bf_cache_name), overwrite = TRUE)

options(timeout = 10 ^ 4)
options(java.parameters = '-Xmx10g') # options('java.parameters')
library(rJava)
library(RBioFormats)
source("xenium.R")
options(timeout = 60)


# multi workers  -----
# library(future)
# plan("multisession", workers = 10, gc = TRUE)
# options(future.globals.maxSize = 10 * 1024 ^ 3)
# 
# nbrOfFreeWorkers()

# options(timeout = 10000)

# setting variable   --------------------
data_dir <- "/data/xenium/XE2601270675-01"
anno_file <- "data/Annotation.csv"
# ROI_anno_file <- "../data/ROI_Annotation.xlsx"
# panel_file <- "SAM000626 xenium-human-gene-expression-panel-gene-list.xlsx"


# reading data -------------
## data path and project name  -----
parse_data_path <-
  dir(data_dir, full.names = TRUE) %>% 
  parseProjName_xenium(pattern = "output[-]?") %>%
  mutate(.counts_path = file.path(.path, "cell_feature_matrix")) %>%
  identity()

## add id for parse_data_path --------
parse_data_path <- 
  parse_data_path %>%
  mutate(.id = seq_along(.proj))



## image w/h ratio ------
# used to fix width and heigth ratio in In Situ Spatial plot.
# 这段代码的主要目的是计算每个图像的长宽比 (width-to-height ratio)，并将其存储在列表 image_wh_r_ls 中
# 现在可以用 coord_equal 调整宽高比，这段代码可以不跑
# image_wh_r_ls <-
#   with(parse_data_path, split(.path, .proj)) %>%
#   map(function(x){
#     imag_wh <-
#       RBioFormats::read.metadata(
#         file.path(x, "morphology.ome.tif")
#       ) %>%
#       .[[1]] %>%
#       .$coreMetadata %>%
#       .[c("sizeX", "sizeY")]
# 
#     image_wh_r <- imag_wh$sizeX/imag_wh$sizeY
#     image_wh_r
#   })


## read into annotation -------
if(file.exists(anno_file)){
  anno_data <- readr::read_csv(anno_file)
  User_columns <- colnames(anno_data)
}else{
  anno_data <- NULL
  User_columns <- "project"
}


## seurat object ------
all_data_ls <-
  parse_data_path %>% 
  dplyr::select(.counts_path, .proj, .id) %>%
  pmap(function(.counts_path, .proj, .id){
    # merge anno_data
    if(!is.null(anno_data)){
      row_idx <- which(anno_data$project == .proj)
      col_idx <- anno_data %>% colnames() 
      
      # 没有找到当前样本的注释数据
      #<<<<<< 所以Annotation.csv也有筛选文件的作用 >>>>>#
      if(length(row_idx) < 1){
        message(str_glue("Not found sample {.proj} in column: project of Annotation table.") )
        return(NULL)
      }
      # 没有任何数据，报错
      if(length(col_idx) < 1) stop("Not found any data in Annotation table.")
      
      meta_dat <- anno_data %>% .[row_idx, col_idx] %>% as.list()
    }else{
      meta_dat <- list(project = .proj)
    }
    
    # 构建Seurat
    sc <-
      Read10X(.counts_path) %>% 
      .[["Gene Expression"]] %>%
      CreateSeuratObject(
        project = .proj,
        min.cells    = 5,
        min.features = 5
      ) %>%
      AddMetaData(meta_dat) %>%
      RenameCells(add.cell.id = .id)
    sc
  }) 

all_data_ls <- all_data_ls %>% purrr::discard(is.null)
names(all_data_ls) <- map_chr(all_data_ls, "project.name")


# parse ROI----------------
# sample_df <-  
#   data.frame(
#     sample = c("kxz1","kxz6","kxz7","kcz2","kcz3","kcz8","gxz10","gxz11","gxz12","gcz13","gcz14","gcz15"),
#     group  = c("AirMale","AirMale","AirMale","AirFemale","AirFemale","AirFemale","HyperoxiaMale","HyperoxiaMale","HyperoxiaMale","HyperoxiaFemale","HyperoxiaFemale","HyperoxiaFemale")
#   )

### sample_df中的内容：sample为样本名，group为对应FOV名
sample_df <- 
  data.frame(
    sample = c("GPHO1", "LXYE1", "DAYI1", "LJME1", "ZXLU1", "GPHO2", "LXYE2", "DAYI2", "ZXLU2", "GNJI2", "WJXI2"),
    group  = c("TMA", "TMA", "TMA", "TMA", "TMA", "TMA", "TMA", "TMA", "TMA", "TMA", "TMA")
  )


roi_ls <- list.files(
  'ROI',
  full.names = TRUE
) 

names(roi_ls) <- 
  roi_ls %>%
  basename() %>%
  str_remove("_cells_stats.csv")


roi_coords_ls <-
  roi_ls %>%
  map(function(path){
    # 读取文件内容，跳过前两行
    dat <- read.csv(path,skip = 2)
    # 将列名转换为小写
    colnames(dat) <- dat %>% colnames() %>% str_to_lower()
    # 添加一个新的列 'cell'，其值与 'cell.id' 列相同
    dat<-dat %>% mutate(cell=cell.id)%>%
      # 将 'cell.id' 列转换为行名
      column_to_rownames('cell.id')
    dat
  }) 

### 根据每个样本所在芯片的顺序编号(names(roi_coords_ls)的顺序)，来对应all_data_ls中的cell id
.id <- c("1_", "1_", "1_", "1_", "1_", "1_", "1_", "1_", "1_", "1_", "1_")
# .id<- str_c(names(roi_coords_ls), "_")
# .id
names(.id) <- names(roi_coords_ls)
roi_coords_ls <- roi_coords_ls %>%
  map2(.id, ~ {
    .x$cell <- str_c(.y, .x$cell)
    return(.x)
  })


xenium_roi_ls <-
  roi_coords_ls %>%
  imap(function(roi, name){ # browser()
    # print(name)
    # return()
    # get Visium Object
    xenium_name <- sample_df$group[sample_df$sample == name]
    xenium <- all_data_ls[[xenium_name]]
    
    # subset visium
    cells <- roi$cell
    xenium_roi <- subset(xenium, cells = cells)
    xenium_roi
  })


# 计算 roi_coords_ls 中每个样本的细胞数
cell_counts <- sapply(roi_coords_ls, function(x) length(x$cell))
cell_counts


# xenium_roi_ls$`F21-6989A`<-all_data_ls$`F21-6989A`
# xenium_roi_ls$`F22-2468A`<-all_data_ls$`F22-2468A`

# add sample_id or sample_name -------------------
xenium_roi_ls <-
  xenium_roi_ls %>%
  imap(function(xenium, name){ #browser()
    xenium$project <- name
    Idents(xenium) <- 'project'
    # RenameCells(visium, add.cell.id = name)
    xenium
  })

# visium_roi override visium
all_data_ls <- xenium_roi_ls





# project number, used as condition judges
proj_count <- length(all_data_ls)
proj_name  <- names(all_data_ls)



sample_df_fix <- 
  data.frame(
    sample = c("GPHO1", "LXYE1", "DAYI1", "LJME1", "ZXLU1", "GPHO2", "LXYE2", "DAYI2", "ZXLU2", "GNJI2", "WJXI2"),
    group  = c("P1", "P1", "P1", "P1", "NP1", "P2", "P2", "P2", "NP2", "NP2", "NP2")
  )

all_data_ls_fix <- all_data_ls %>%
  imap(function(data, name) {
    # 获取正确的组名
    correct_group <- sample_df_fix$group[sample_df_fix$sample == name]
    correct_project <- sample_df_fix$sample[sample_df_fix$sample == name]
    # 检查是否找到正确的组名
    if (length(correct_group) == 0) {
      warning(paste("No group found for sample:", name))
      return(data)
    }
    
    # 添加或更新组名到元数据
    data@meta.data$group <- correct_group
    data@meta.data$sample <- correct_project
    # 返回更新后的数据
    data
  })

# 验证组名是否更新正确
print(unique(all_data_ls_fix$GPHO1@meta.data$group))
print(unique(all_data_ls_fix$GPHO1@meta.data$project))


all_data_ls <- all_data_ls_fix


# merged object
all_data <- merge(all_data_ls[[1]], all_data_ls[-1])





## read coords -------------------
# cell coordinate and cell bundary
cells_position_ls <- 
  # with(parse_data_path, split(.path, .proj)) %>% 
  parse_data_path %>% 
  group_by(.proj) %>% 
  nest %>% 
  deframe() %>% 
  map(function(x){
    xx <-
      readr::read_csv(
        file.path(x$.path, "cells.csv.gz"), 
        col_types = list(cell_id = "c")
      ) %>% 
      mutate(cell_id = str_c(x$.id, cell_id, sep = "_"))
    xx
  })

roi_position_ls <- 
  roi_coords_ls %>% 
  imap(function(roi, name){  #browser()
    # get Visium Object
    xenium_name <- sample_df$group[sample_df$sample == name]
    xenium <- cells_position_ls[[xenium_name]]
    
    # subset visium
    cells <- roi$cell
    position_roi <-  xenium[which(xenium$cell_id %in% cells),]
    position_roi
  })

cells_position_ls <- cells_position_ls[proj_name]
cells_position_ls <- cells_position_ls %>% purrr::discard(is.null)
cells_position_ls <- c(roi_position_ls,cells_position_ls)


cells_boundary_ls <- 
  # with(parse_data_path, split(.path, .proj)) %>%
  parse_data_path %>%
  group_by(.proj) %>%
  nest %>%
  deframe() %>%
  map(function(x){
    xx <-
      readr::read_csv(
        file.path(x$.path, "cell_boundaries.csv.gz"),
        col_types = list(cell_id = "c")
      ) %>%
      mutate(cell_id = str_c(x$.id, cell_id, sep = "_"))
    xx
  })

roi_boundary_ls <-
  roi_coords_ls %>%
  imap(function(roi, name){  #browser()
    # get Visium Object
    xenium_name <- sample_df$group[sample_df$sample == name]
    xenium <- cells_boundary_ls[[xenium_name]]
    
    # subset visium
    cells <- roi$cell
    boundary_roi <- xenium[which(xenium$cell_id %in% cells),]
    boundary_roi
  })

cells_boundary_ls <- cells_boundary_ls[proj_name]
cells_boundary_ls <- cells_boundary_ls %>% purrr::discard(is.null)
cells_boundary_ls <- c(roi_boundary_ls,cells_boundary_ls)



# export data
dir.create("data_for_analysis")
saveRDS(cells_boundary_ls, file.path("data_for_analysis", "cells_boundary_ls.rds"))
saveRDS(cells_position_ls, file.path("data_for_analysis", "cells_position_ls.rds"))
saveRDS(all_data_ls, file.path("data_for_analysis", "all_data_ls.rds"))
# saveRDS(image_wh_r_ls, file.path("data_for_analysis", "image_wh_r_ls.rds"))






# FORMAT CONVERSION ---------
# config python
library(reticulate)
# use_condaenv(condaenv = "/home/szt/miniforge3/envs/rpy/bin/python")
use_condaenv(condaenv = "/home/liuran/software/miniconda3/envs/rpy/bin/python")
py_config()
py_available()
py_module_available("pandas")


sc <- import("scanpy")
adata <- sc$read_h5ad("adata_celltype.h5ad")

adata_meta <- adata$obs
adata_umap <- adata$obsm['X_umap']

# # sample match
# if("sample_raw" %in% colnames(adata_meta)){
#   adata_meta_new <-
#     adata_meta %>% 
#     # left_join(sample_match_df2, by = c(sample = "to"))
#     mutate(
#       from = sample_raw
#     )
# }else{
#   # find original sample_id
#   sample_match_df <- 
#     data.frame(
#       from = c( # project name
#         "NC1", "NC2", "NC3", "NC4", "CTE1", "CTE2", "CTE3", "CTE4"
#       ),
#       to   = c( # sample name
#         "NC1", "NC2", "NC3", "NC4", "CTE1", "CTE2", "CTE3", "CTE4"
#       )
#     ) 
#   
#   if(nrow(sample_match_df) > 0){
#     sample_match_df2 <-
#       sample_match_df %>%
#       right_join(parse_data_path['.proj'], by = c(from = '.proj')) %>%
#       mutate(to = ifelse(is.na(to), from, to))
#   }else{
#     sample_match_df2 <- 
#       parse_data_path %>%
#       dplyr::select(
#         from = .proj, 
#         to = .proj,
#       )
#   }
#   
#   
#   adata_meta_new <-
#     adata_meta %>% 
#     left_join(sample_match_df2, by = c(sample = "to"))
# }




# stopifnot(all(!is.na(adata_meta$sample)))
# 
# if("cell_id_full" %in% colnames(adata_meta_new)){
#   adata_meta_new$cell_id2 <- adata_meta_new$cell_id_full
# }else{
#   adata_meta_new$cell_id2 <- str_c(adata_meta_new$from, adata_meta_new$cell_id, sep = "_")
# }

# # some cells duplicated 
# dup_cell_ids <- adata_meta_new$cell_id2 %>% .[duplicated(.)]
# dup_cell_ids <- adata_meta$cell_id_fov  %>% .[duplicated(.)]
# 
# dup_cell_idx <- adata_meta$cell_id_fov %in% dup_cell_ids
# 
# if(sum(dup_cell_idx) > 0){
#   warning(str_glue("{sum(dup_cell_idx)} cell_ids are duplicated and have been removed"))
#   adata_meta_new <- adata_meta_new[!dup_cell_idx, ]
#   adata_umap <- adata_umap[!dup_cell_idx, ]
# }





library(dplyr)
adata_meta <- adata_meta %>%
  mutate(cell_id2 = paste0("1_", cell_id))
head(adata_meta)




# add METADATA
all_data2 <- all_data %>% subset(cells = adata_meta$cell_id2)
all_data2@meta.data <- 
  all_data2@meta.data %>%
  rownames_to_column("cell_id2") %>%
  left_join(adata_meta, by = "cell_id2") %>%
  column_to_rownames("cell_id2")


# add UMAP
adata_umap2 <- adata_umap
colnames(adata_umap2) <- str_c("UMAP", 1:ncol(adata_umap2), sep = "_")
rownames(adata_umap2) <- adata_meta$cell_id2
adata_umap2 <- adata_umap2[colnames(all_data2), ]

all_data2[['umap']] <- 
  CreateDimReducObject(
    adata_umap2,
    key    = "UMAP_",
    assay  = DefaultAssay(all_data2),
    global = TRUE
  )

# add useful columns
cluster_res_cols <- all_data2[[]] %>% colnames() %>% str_subset("^res_")
for(col_name in cluster_res_cols){
  col_name_new <- col_name %>% str_replace("res_", str_glue("{DefaultAssay(all_data2)}_snn_res.")) 
  all_data2[[col_name_new]] <- all_data2[[col_name]][[1]]
}
all_data2$annotation <- all_data2$celltype


# plot
pdf(file.path(dr_outdir, str_glue("Dimplot all_data2.pdf")), height = 6, width = 8)
all_data2 %>% DimPlot(reduction = "umap", group.by = "celltype", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "annotation", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "res_0.6", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "RNA_snn_res.0.6", label = TRUE) + NoLegend()
dev.off()



# split ROI ------------
colnames(all_data2@meta.data)[colnames(all_data2@meta.data) == "group.x"] <- "group"
colnames(all_data2@meta.data)[colnames(all_data2@meta.data) == "sample.x"] <- "sample"

#1 split layers
all_data2 <- all_data2 %>% JoinLayers() %>% split(all_data2$sample)


# #2 cells_position_ls
# #3 cells_boundary_ls
# roi_match_ls <- 
#   all_data2[[c('from', 'sample')]] %>%
#   distinct() %>%
#   group_by(from) %>%
#   nest() %>%
#   deframe() %>%
#   map(unlist, use.names = FALSE)
# 
# roi_cell_id_ls <-
#   all_data2[['sample']] %>% 
#   rownames_to_column("cell_id") %>%
#   group_by(sample) %>%
#   nest() %>%
#   deframe() %>%
#   map(unlist, use.names = FALSE)
# 
# 
# cells_position_ls2  <- list()
# cells_boundary_ls2  <- list()
# 
# for(idx in seq_along(roi_match_ls)){
#   from_name <- names(roi_match_ls)[idx]
#   to_name   <- roi_match_ls[[idx]]
#   
#   print(str_glue("###### {from_name} ->"))
#   for(sub in to_name){
#     print(str_glue("  -> {sub} ######"))
#     cell_id_sub <- roi_cell_id_ls[[sub]]
#     
#     cells_position_new <- cells_position_ls[[from_name]] %>% dplyr:::filter(cell_id %in% cell_id_sub)
#     cells_boundary_new <- cells_boundary_ls[[from_name]] %>% dplyr:::filter(cell_id %in% cell_id_sub)
#     
#     cells_position_ls2[[sub]] <- cells_position_new
#     cells_boundary_ls2[[sub]] <- cells_boundary_new
#   }
# }
# 
# # export data
# saveRDS(cells_boundary_ls2, file.path("R_OUTPUT", "data_for_analysis", "cells_boundary_ls.rds"))
# saveRDS(cells_position_ls2, file.path("R_OUTPUT", "data_for_analysis", "cells_position_ls.rds"))


p1 <- cells_boundary_ls[['GPHO1']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p2 <- cells_boundary_ls[['LXYE1']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p3 <- cells_boundary_ls[['DAYI1']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p4 <- cells_boundary_ls[['LJME1']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p <- cells_boundary_ls[['GPHO1']] %>%
  slice_sample(prop = 0.25) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

library(patchwork)
((p2 | p1)/ p3)  / p
p1 / p2 / p3  / p4





# normalize ----------------------------------------------------------------------
all_data2 <- NormalizeData(all_data2, scale.factor = median(all_data2$nCount_RNA))
all_data2 <- FindVariableFeatures(all_data2)

# PCA ----------------------------------------------------------------------------
## PCA --------
#DefaultAssay(all_data2) <- "sketch"
#all_data2 <- FindVariableFeatures(all_data2)
VariableFeaturePlot(all_data2)

all_data2 <- ScaleData(all_data2, vars.to.regress = 'nCount_RNA')
all_data2 <- RunPCA(all_data2)

# dr_outdir <- '2. dimension reduction and cluster/'
### pca inspector ------
pdf(file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data2, ndims = 30)
dev.off()

all_data2@meta.data %>% head()

all_data2$cell_id <- NULL

set.seed(1234)
suset_cells <-
  all_data2[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  slice_sample(prop = 50000/ncol(all_data2) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(dr_outdir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 1:5)
dev.off()

pdf(file.path(dr_outdir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data2[, suset_cells], dims = 16:20)
dev.off()


pdf(file.path(dr_outdir, str_glue("Dimplot all_data2 - sample.pdf")), height = 6, width = 8)
DimPlot(all_data2, group.by = "sample", label = TRUE, reduction = "umap") + NoLegend()
dev.off()



# ## dimension reduction -----------------------------------------------------------------
# dims <- 20
# 
# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
# 
# 
# #all_data2 <- RunUMAP(all_data2, dims = 1:dims, return.model = TRUE, n.neighbors = 30)
# 
# all_data2 <- all_data2 %>%
#   harmony::RunHarmony(
#     dims.use = 1:dims,
#     group.by.vars = "sample",
#     reduction.save = "harmony"
#   )
# 
# all_data2 <- RunUMAP(
#   all_data2,
#   dims = 1:dims,
#   reduction = "harmony",
#   reduction.name = "harmonyumap",
#   reduction.key  = "harmonyUMAP_", 
#   return.model = TRUE
# )
# 
# pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} sample - harmony.pdf")), height = 6, width = 8)
# DimPlot(all_data2, group.by = "sample", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data2, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
# dev.off()
# 
# pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 8)
# DimPlot(all_data2, group.by = "group", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data2, group.by = "group", label = TRUE, reduction = "harmonyumap") + NoLegend()
# dev.off()

# export 
saveRDS(all_data2, file = file.path("data_for_analysis", "all_data.rds"))













