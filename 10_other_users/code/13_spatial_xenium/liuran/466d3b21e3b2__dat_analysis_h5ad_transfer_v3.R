library(tidyverse)
library(Seurat)
library(patchwork)


# setting variable   --------------------
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)

dir.create(dr_outdir)
dir.create(sp_output)

data_dir  <- "/mnt/NexData_file/home/liuran/projects/3_Xenium/XE2512247768-01/XE2512247768-01"
list.files(data_dir)
anno_file <- "data/Annotation.csv" # 非必须

# dir.create("R_OUTPUT")
dir.create("data_for_analysis")


# library --------------
options(java.parameters = '-Xmx10g') # options('java.parameters')
library(rJava)
library(RBioFormats)
source("xenium.R")

# multi workers  -----
# library(future)
# plan("multisession", workers = 10, gc = TRUE)
# 
# options(future.globals.maxSize = 10 * 1024 ^ 3)
# 
# nbrOfFreeWorkers()



# reading data -------------
## data path and project name  -----
parse_data_path <-
  dir(data_dir, full.names = TRUE) %>% 
  parseProjName_xenium(pattern = "output[-]?") %>%
  mutate(.counts_path = file.path(.path, "cell_feature_matrix")) %>%
  identity()

## image w/h ratio ------
# used to fix width and heigth ratio in In Situ Spatial plot.
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
  dplyr::select(.counts_path, .proj) %>%
  pmap(function(.counts_path, .proj){
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
      # RenameCells(add.cell.id = .id)
      RenameCells(add.cell.id = .proj)
    sc
  }) 

all_data_ls <- all_data_ls %>% purrr::discard(is.null)
names(all_data_ls) <- map_chr(all_data_ls, "project.name")

# project number, used as condition judges
proj_name  <- names(all_data_ls)


# merged object
all_data <- merge(all_data_ls[[1]], all_data_ls[-1])


## read coords -------------------
# cell coordinate and cell boundary
cells_position_ls <- 
  # with(parse_data_path, split(.path, .proj)) %>%
  parse_data_path %>%
  group_by(.proj) %>%
  nest %>%
  deframe() %>%
  imap(function(x, proj){
    xx <-
      readr::read_csv(
        file.path(x$.path, "cells.csv.gz"),
        col_types = list(cell_id = "c")
      ) %>%
      mutate(cell_id = str_c(proj, cell_id, sep = "_"))
    xx
  })
cells_position_ls <- cells_position_ls[proj_name]

cells_boundary_ls <- 
  # with(parse_data_path, split(.path, .proj)) %>%
  parse_data_path %>%
  group_by(.proj) %>%
  nest %>%
  deframe() %>%
  imap(function(x, proj){
    xx <-
      readr::read_csv(
        file.path(x$.path, "cell_boundaries.csv.gz"),
        col_types = list(cell_id = "c")
      ) %>%
      mutate(cell_id = str_c(proj, cell_id, sep = "_"))
    xx
  })
cells_boundary_ls <- cells_boundary_ls[proj_name]

# export data
saveRDS(cells_boundary_ls, file.path("data_for_analysis", "cells_boundary_ls_raw.rds"))
saveRDS(cells_position_ls, file.path("data_for_analysis", "cells_position_ls_raw.rds"))




# FORMAT CONVERSION ---------
# config python
library(reticulate)
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
#         'F2410395-A3_F2410396-A3-1_F2410396-A3-2', 'F2410395-A3_F2410396-A3-1_F2410396-A3-2', 'F2410395-A3_F2410396-A3-1_F2410396-A3-2'
#       ),
#       to   = c( # sample name
#         "PC", "PA_1", "PA_2"
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
# 
# 
# 
# 
# stopifnot(all(!is.na(adata_meta_new$from)))

# if("cell_id_full" %in% colnames(adata_meta_new)){
#   adata_meta_new$cell_id2 <- adata_meta_new$cell_id_full
# }else{
#   adata_meta_new$cell_id2 <- str_c(adata_meta_new$from, adata_meta_new$cell_id, sep = "_")
# }

adata_meta$cell_id2 <- adata_meta$cell_id

id_df <- data.frame(
  id = c("TMA1_", "TMA1_", "TMA1_", "TMA1_", "TMA1_", "TMA1_", "TMA1_", "TMA1_", "TMA1_", "TMA1_", "TMA1_", "TMA1_", "TMA2_", "TMA2_", "TMA2_", "TMA2_", "TMA2_", "TMA2_", "TMA2_", "TMA2_", "TMA2_", "TMA2_", "TMA2_", "TMA2_", "TMA2_", "TMA2_", "TMA3_", "TMA3_", "TMA3_", "TMA3_", "TMA3_", "TMA3_", "TMA3_", "TMA3_", "TMA3_", "TMA3_", "TMA3_", "TMA3_", "TMA3_", "TMA4_", "TMA4_", "TMA4_", "TMA4_", "TMA4_", "TMA4_", "TMA4_", "TMA4_", "TMA4_", "TMA4_", "TMA4_", "TMA4_", "TMA4_", "TMA4_"),
  sample = c("Ps_1", "Ps_2", "Ps_3", "Ps_4", "Ps_5", "Ps_6", "Ps_7", "Ps_8", "Ps_9", "Ps_10", "Ps_11", "Ps_12", "Ps_13", "Ps_14", "Ps_15", "Ps_16", "Ps_17", "Ps_18", "Ps_19", "Ps_20", "Ps_21", "Ps_22", "Ps_23", "Ps_24", "Ps_25", "Ps_26", "Ps_27", "Ps_28", "Ps_29", "Ps_30", "Ps_31", "Ps_32", "Ps_33", "Ps_34", "Ps_35", "Ps_36", "Ps_37", "Ps_38", "Ps_39", "NC_1", "NC_2", "NC_3", "NC_4", "NC_5", "NC_6", "NC_7", "NC_8", "NC_9", "NC_10", "NC_11", "NC_12", "NC_13", "NC_14"))

adata_meta_new <- adata_meta %>% 
  left_join(id_df, by = "sample") %>% 
  mutate(cell_id2 = ifelse(!is.na(id), paste0(id, cell_id2), cell_id2)) %>% 
  select(-id)
head(adata_meta_new)


# # some cells duplicated 
# dup_cell_ids <- adata_meta_new$cell_id2 %>% .[duplicated(.)]
# dup_cell_idx <- adata_meta_new$cell_id2 %in% dup_cell_ids
# 
# if(sum(dup_cell_idx) > 0){
#   warning(str_glue("{sum(dup_cell_idx)} cell_ids are duplicated and have been removed"))
#   adata_meta_new <- adata_meta_new[!dup_cell_idx, ]
#   adata_umap <- adata_umap[!dup_cell_idx, ]
# }



# add METADATA
# all_data@meta.data$sample <- gsub("-", "_", all_data@meta.data$project)

all_data2 <- all_data %>% subset(cells = adata_meta_new$cell_id2)
all_data2@meta.data$group <- NULL

dim(all_data)
dim(all_data2)
dim(adata_meta)

all_data2@meta.data <- 
  all_data2@meta.data %>%
  rownames_to_column("cell_id2") %>%
  left_join(adata_meta_new, by = "cell_id2") %>%
  column_to_rownames("cell_id2")


# add UMAP
adata_umap2 <- adata_umap
colnames(adata_umap2) <- str_c("UMAP", 1:ncol(adata_umap2), sep = "_")
rownames(adata_umap2) <- adata_meta_new$cell_id2
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
all_data2 %>% DimPlot(reduction = "umap", group.by = "celltype", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "annotation", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "res_0.8", label = TRUE) + NoLegend()
all_data2 %>% DimPlot(reduction = "umap", group.by = "RNA_snn_res.0.8", label = TRUE) + NoLegend()




# split ROI ------------
# 20250421 add

#1 split layers
all_data2 <- all_data2 %>% JoinLayers() %>% split(all_data2$sample)


#2 cells_position_ls
#3 cells_boundary_ls
roi_match_ls <- 
  all_data2[[c('orig.ident', 'sample')]] %>%
  distinct() %>%
  group_by(orig.ident) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.names = FALSE)

roi_cell_id_ls <-
  all_data2[['sample']] %>% 
  rownames_to_column("cell_id") %>%
  group_by(sample) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.names = FALSE)


cells_position_ls2  <- list()
cells_boundary_ls2  <- list()

for(idx in seq_along(roi_match_ls)){
  from_name <- names(roi_match_ls)[idx]
  to_name   <- roi_match_ls[[idx]]
  
  print(str_glue("###### {from_name} ->"))
  for(sub in to_name){
    print(str_glue("  -> {sub} ######"))
    cell_id_sub <- roi_cell_id_ls[[sub]]
    
    cells_position_new <- cells_position_ls[[from_name]] %>% dplyr:::filter(cell_id %in% cell_id_sub)
    cells_boundary_new <- cells_boundary_ls[[from_name]] %>% dplyr:::filter(cell_id %in% cell_id_sub)
    
    cells_position_ls2[[sub]] <- cells_position_new
    cells_boundary_ls2[[sub]] <- cells_boundary_new
  }
}

# export data
saveRDS(cells_boundary_ls2, file.path("data_for_analysis", "cells_boundary_ls.rds"))
saveRDS(cells_position_ls2, file.path("data_for_analysis", "cells_position_ls.rds"))


p1 <- cells_boundary_ls2[['Ps2344413']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p2 <- cells_boundary_ls2[['Ps20230296']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p3 <- cells_boundary_ls2[['Ps20230169']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()


p <- cells_boundary_ls[['TMA']] %>%
  slice_sample(prop = 0.25) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

library(patchwork)
((p2 | p1)/ p3)  / p



p1 <- cells_boundary_ls[['TMA1']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p2 <- cells_boundary_ls[['TMA2']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p3 <- cells_boundary_ls[['TMA3']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()

p4 <- cells_boundary_ls[['TMA4']] %>%
  slice_sample(prop = 0.5) %>%
  
  ggplot(aes(x = vertex_x, y = vertex_y)) +
  geom_polygon(aes(group = cell_id)) +
  scale_y_reverse() +
  coord_equal() +
  theme_void()
p1 | p2 | p3 | p4


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


pdf(file.path(dr_outdir, "Dimplot sample.pdf"), height = 6, width = 7)
DimPlot(all_data2, group.by = "sample", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, "Dimplot annotation.pdf"), height = 6, width = 7)
DimPlot(all_data2, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
dev.off()

pdf(file.path(dr_outdir, "Dimplot annotation group.pdf"), height = 6, width = 12)
DimPlot(all_data2, group.by = "annotation", label = TRUE, reduction = "umap", split.by = "group")
dev.off()

pdf(file.path(dr_outdir, "Dimplot group.pdf"), height = 6, width = 7)
DimPlot(all_data2, group.by = "group", label = TRUE, reduction = "umap") + NoLegend()
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
saveRDS(all_data2, file = file.path("data_for_analysis", "all_data_h5ad.rds"))





# spatial plot --------------
all_data_meta <- 
  all_data2[[]] %>% 
  rownames_to_column("cell_id")

plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#600473",
          "#7F95D1", "#dee2d1", "#ef5285", "#a5dff9", "#e8a0b8", "#ffeead", 
          "#737495", "#bccf3d", "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", 
          "#fed95c", "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", 
          "#FFE869", "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", 
          "#F0D879", "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", 
          "#f06966", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#60c5ba", "#feee7d", "#ffc300", "#02c9c9", "#fbf579", "#4695d6", 
          "#005995", "#fa625f", "#f17d80")

sp_p_ls <- 
  cells_boundary_ls2 %>% 
  # .[1] %>% 
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    dat <- 
      cells_boundary %>% 
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>% 
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation, group = cell_id)) + 
      scale_y_reverse() + 
      coord_equal() + 
      scale_fill_manual(values = plan) + 
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) + 
      theme_void()
    p
  })

# export image
sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path(sp_output, str_glue("spatial plot {grp}.pdf")), 
      height = 12, 
      width  = 12
    )
    print(p)
    dev.off()
  })


sp_p_TMA1_3 <- 
  cells_boundary_ls %>% 
  .[c("TMA1", "TMA2", "TMA3")] %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    dat <- 
      cells_boundary %>% 
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>% 
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation, group = cell_id)) + 
      scale_y_reverse() + 
      coord_equal() + 
      scale_fill_manual(values = plan) + 
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) + 
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_text(vjust = 0, hjust = 0.5, size = 25)) +
      labs(title = project)
    p
  })

sp_p_TMA4 <- 
  cells_boundary_ls %>% 
  .["TMA4"] %>%
  imap(function(cells_boundary, project){ # browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    dat <- 
      cells_boundary %>% 
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>% 
      ggplot(aes(x = vertex_x, y = vertex_y)) + 
      geom_polygon(aes(fill = annotation, group = cell_id)) + 
      scale_y_reverse() + 
      coord_equal() + 
      scale_fill_manual(values = plan) + 
      guides(fill = guide_legend(override.aes = list(size = 5, color = "black"))) + 
      theme_void() +
      theme(legend.position = "right",
            plot.title = element_text(vjust = 0, hjust = 0.5, size = 25)) +
      labs(title = project)
    p
  })

p_combined <- wrap_plots(
  sp_p_TMA1_3$TMA1, sp_p_TMA1_3$TMA2, sp_p_TMA1_3$TMA3, sp_p_TMA4$TMA4,
  nrow = 1, ncol = 4
)
ggsave(file.path(sp_output, "spatial plot combined.pdf"), p_combined, width = 40, height = 20)



