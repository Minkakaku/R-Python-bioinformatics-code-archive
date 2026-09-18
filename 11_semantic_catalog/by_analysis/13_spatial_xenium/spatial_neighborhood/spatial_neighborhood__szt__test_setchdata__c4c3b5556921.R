library(Seurat)
library(tidyverse)
library(BPCells)


out_dir <- 'skech_data_analysis'
dir.create(out_dir)



# 官网教程中，使用130万细胞，索引了5万个细胞
# 其数据是多样本数据，已经是cellranger aggr合并的数据了
# 注意，SketchData是每个data layer进行一次操作，因此最终的细胞数应该是ncells * length(data layer)

# 多样本数据
# 策略1 直接JointLayers，索引细胞即可, 数据太大，可能报错
# Error in object %*% (R.inv %*% JL) : 
#   CHOLMOD error 'problem too large' at file '../Core/cholmod_sparse.c', line 89

# 策略2 每个样本单独索引细胞，根据获得细胞ID后，单独索引细胞并构建数据 】
# 参考SketchData的后半部分代码



all_data <- readRDS("all_data_merge.rds") 

# method1 ------
# 放弃，改为单个执行
all_data <- all_data %>% JoinLayers()
all_data <- all_data %>% NormalizeData()
all_data <- all_data %>% FindVariableFeatures()


all_data <- all_data %>% 
  SketchData(
    ncells = 50000, # 50000 per layer
    method = 'LeverageScore', # 'LeverageScore',  
    sketched.assay = "sketch"
  )


# method2 ------

# all_data_diet <- all_data 
# 
# all_data_diet[['RNA']]$data        <- NULL # new(Class = "dgCMatrix")
# all_data_diet[['RNA']]$scale.data  <- NULL # matrix(nrow = 0, ncol = 0)

# format(object.size(all_data), 'MB')
# format(object.size(all_data_diet), 'MB')

# Layers(all_data_diet)


# all_data_diet <- all_data_diet %>% split(all_data_diet$project)
# all_data_diet <- all_data_diet %>% NormalizeData()

# LayerData(
#   all_data_diet,
#   layer = 'counts.00037G2', 
#   features = VariableFeatures(all_data_diet, layer = 'counts.00037G2')
# ) %>%
# LeverageScore 

all_data_layer <- Layers(all_data, search = 'data') %>% setNames(., .)
all_lev_score <- 
  all_data_layer %>%
  map(function(lyr){
    dat <-
      LayerData(
        all_data,
        layer = lyr, 
        features = VariableFeatures(all_data, layer = lyr)
      )
    lev_score <- LeverageScore(dat)
    lev_score
  })

number_of_skech <- 50000
number_of_skech_scale_factor <-  number_of_skech / ncol(all_data)


all_sketch_cells <- 
  all_lev_score %>%
  imap(function(score, lyr){ # browser()
    # if(lyr == 'data.09925E1') browser()
     
    dat <-
      LayerData(
        all_data,
        layer = lyr, 
        features = VariableFeatures(all_data, layer = lyr)
      )
    all_cells <- colnames(dat)
    
    stopifnot(length(score) == ncol(dat))
    if(!is.null(names(score))){
      stopifnot(all(names(score) == colnames(dat)))
    }
    
    set.seed(1234)
    subset_cells <- sample(all_cells, size = round(ncol(dat) * number_of_skech_scale_factor), prob = score)
    
    subset_cells
  })

# subset assays
sketched <- subset(
  all_data[['RNA']], 
  cells  = unlist(all_sketch_cells), 
  layers = Layers(all_data[['RNA']], search = c("counts", "data"))
)

# set VariableFeature
for (lyr in all_data_layer) {
  try(
    VariableFeatures(sketched, method = "sketch", layer = lyr) <- 
      VariableFeatures(all_data[['RNA']], layer = lyr), 
    silent = TRUE
    )
}

if (inherits(sketched, "Assay5")) {
  sketched <- CastAssay(sketched, to = 'dgCMatrix')
}

# others slot setting
Key(sketched) <- Key('sketch', quiet = TRUE)
all_data[['sketch']] <- sketched
DefaultAssay(all_data) <- 'sketch'


format(object.size(all_data), 'GB')

# saveRDS(all_data, 'all_data_sketch.rds')

# dimension reduction -----------
all_data <- ScaleData(all_data, vars.to.regress = 'nCount_RNA')
all_data <- RunPCA(all_data)


## pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()

pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data, dims = 16:20)
dev.off()


 
## multi core -------
library(future)
plan(multicore, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


##  remove batch effect ------
dims <- 13

all_data <- RunUMAP(all_data, dims = 1:dims)

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE) + NoLegend()
dev.off()


all_data <- all_data %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

all_data <- RunUMAP(
  all_data, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


## add metadata ------
anno_data <- read.table("Annotation.txt", header = TRUE)

all_data_meta <-
  all_data@meta.data %>%
  rownames_to_column("cell_id") %>%
  left_join(anno_data, by = c(project = 'sample')) %>% 
  column_to_rownames('cell_id') 

stopifnot(all(rownames(all_data@meta.data) == rownames(all_data_meta)))
all_data@meta.data <- all_data_meta


## 可视化

reduction <- "harmony" # pca harmony
reduction_show <- ifelse(reduction == 'pca', "uamp", "harmonyumap")



pdf(file.path(out_dir, str_glue("Dimplot dim{dims} group - harmony.pdf")), height = 6, width = 15)
DimPlot(all_data, group.by = "group", label = TRUE, reduction = reduction_show) + NoLegend()
dev.off()



## cluster ------
### cluster
all_data <-
  FindNeighbors(
    all_data,
    reduction = reduction,
    dims = 1:dims,
    k.param = 20
  )


all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = reduction_show) + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = reduction_show) + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = reduction_show) + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = reduction_show) + NoLegend()
dev.off()



all_data[[]] %>% 
  .[!is.na(.$sketch_snn_res.0.4), ] %>%
  group_by(sketch_snn_res.0.4, group) %>% 
  summarise(count =n()) %>%
  ggplot(aes(x = group, y = count, fill = sketch_snn_res.0.4)) +
  geom_col() +
  ggsci::scale_fill_igv()



## dotplot -------
major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  
  "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  # "Fibroblasts" = "ACTA2",
  
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "MYH11", 
  "Smooth muscle cells" = "ACTA2",
  
  "Breast glandular cells" = "SERPINA3", 
  "Breast glandular cells" = "DSC2",
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  # "Plasma" = "JCHAIN",
  
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "CD4",
  "T cells" = "CD8A",
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  "Myeloid" = "ITGAM",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Dendritic cell" = "CD83",
  'Dendritic cell' = "LILRA4",  
  
  'Dendritic cell' = "CLEC9A",
  'Dendritic cell' = "IDO1",
  
  'Dendritic cell' = "CD1C",
  'Dendritic cell' = "FCER1A"
  
)


# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()


## annotation ----

cell_cluster <- c(
  '0'  = "CD4T",
  '1'  = "B",
  '2'  = "Macrophage",
  '3'  = "Epithelial",
  '4'  = "Fibroblast",
  '5'  = "CD8T",
  '6'  = "Endothelial",
  '7'  = "Plasma",
  '8'  = "Endothelial - macrophage",
  '9'  = "B",
  '10' = "Fibroblast",
  '11' = "DC"
)


Idents(all_data) <- "sketch_snn_res.0.4"

all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = reduction_show) + NoLegend()
dev.off()


# project data -----------
plan(multicore, workers = 1)


all_data <- 
  all_data %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "pca.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'
    
    dims = 1:13, # 50 # <---- harmony只有13个dims
    refdata = list(annotation_full = "annotation")
  )

# test
DefaultAssay(all_data) <- 'RNA'

pdf(file.path(out_dir, str_glue("Dimplot dim{dims} annotation - harmony -projectdata.pdf")), height = 6, width = 8)
DimPlot(all_data, reduction = "full.harmonyumap", group.by = "annotation_full", label = TRUE) + NoLegend()
dev.off()


# 官方教程中 
# 将integrated的ProjectData分为两步骤
# ProjectIntegration返回的DimReducObject的Loading是来源于reference（sketch）的pca/harmony的loading
# 这个loading会被ProjectData用于创建所有细胞的Embedding
# so, 
# ProjectIntegration + ProjectData == ProjectData ??


all_data2 <- all_data
DefaultAssay(all_data2) <- 'sketch'

# add reduction: harmony_full
all_data2 <-
  all_data2 %>% 
  ProjectIntegration(
    sketched.assay = "sketch",
    assay = "RNA",
    reduction = "harmony"
  )


undebug(ProjectIntegration)

all_data2 <- 
  all_data2 %>%
  ProjectData(
    sketched.assay = "sketch",
    sketched.reduction = "harmony.full", 
    assay = "RNA",
    full.reduction = "harmony.full.projected",
    umap.model = "harmonyumap",
    dims = 1:13,
    refdata = list(annotation_full = "annotation")
  )


all_data2 <- 
  all_data2 %>%
  RunUMAP(
    reduction = "harmony.full.projected",
    dims = 1:13,
    reduction.name = "umap.harmony.full",
    reduction.key  = "UMAP_full_"
  )


DimPlot(all_data2, reduction = "full.harmonyumap", group.by = "annotation_full")
DimPlot(all_data2[, all_data2$annotation_full.score > 0.5], reduction = "umap.harmony.full", group.by = "annotation_full", label = TRUE) + NoLegend()


# test SCTranform-----
# SCT 效果不佳，直接放弃 ---

all_data3 <- all_data
DefaultAssay(all_data3) <- "sketch"
all_data3 <- SCTransform(all_data3, assay = "sketch", new.assay.name = 'SCT.sketch')

DefaultAssay(all_data3) <- 'SCT.sketch'
all_data3 <- RunPCA(all_data3)


## dimension reduction ------
dims_sct <- 10

all_data3 <- RunUMAP(all_data3, dims = 1:dims_sct, return.model = TRUE)


DimPlot(all_data3, group.by = "project", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(all_data3, group.by = "group",   reduction = "umap", label = TRUE) + NoLegend()
DimPlot(all_data3, reduction = "umap", label = TRUE) + NoLegend()

all_data3 <- all_data3 %>% 
  harmony::RunHarmony(
    dims.use = 1:dims_sct,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

all_data3 <- RunUMAP(
  all_data3, 
  dims = 1:dims_sct, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims_sct} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data3, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data3, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


## cluster ------
all_data3 <-
  FindNeighbors(
    all_data3,
    reduction = reduction,
    dims = 1:dims_sct,
    k.param = 20
  )


all_data3 <-
  FindClusters(
    all_data3,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims_sct} cluster - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data3, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data3, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data3, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data3, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(all_data3, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data3, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data3, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(all_data3, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()


## dotplot -------
pdf(file.path(out_dir, str_glue("dotplot dims{dims_sct} major_marker.pdf")), width = 25, height = 6)
DotPlot(all_data3, group.by = "sketch_snn_res.0.2", features = major_marker_ls)
DotPlot(all_data3, group.by = "sketch_snn_res.0.4", features = major_marker_ls)
DotPlot(all_data3, group.by = "sketch_snn_res.0.6", features = major_marker_ls)
DotPlot(all_data3, group.by = "sketch_snn_res.0.8", features = major_marker_ls)
dev.off()


# adjust UMAP spread ------------
all_data4 <- all_data
DefaultAssay(all_data4) <- "sketch"

all_data4 <- 
  all_data4 %>%
  RunUMAP(
    dims = 1:dims_sct, 
    reduction = "harmony", 
    reduction.name = "harmonyumap", 
    reduction.key  = "harmonyUMAP_",
    return.model = TRUE,
    
    min.dist = 0.3,
    spread   = 1.2,
    local.connectivity = 100
  )

DimPlot(all_data4, group.by = "annotation", label = TRUE, reduction = "harmonyumap") + NoLegend()


#  test cluster
### cluster
DefaultAssay(all_data4) <- "sketch"

all_data4 <-
  FindNeighbors(
    all_data4,
    reduction = 'harmony',
    dims = 1:dims,
    annoy.metric = "cosine",
    k.param = 30,
    prune.SNN = 1/10
  )


all_data4 <-
  FindClusters(
    all_data4,
    resolution = 0.2, #  c(0.2, 0.4, 0.6, 0.8)
  )
  

DimPlot(all_data4, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data4, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data4, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data4, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data4, group.by = "sketch_snn_res.1",   label = TRUE, reduction = "harmonyumap") + NoLegend()



major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  Epithelial  = "EPCAM", 
  Epithelial  = "KRT5",
  Epithelial  = "CEACAM6",
  
  "Fibroblasts" = "SFRP4", 
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  # "Fibroblasts" = "ACTA2",
  
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "MYH11", 
  "Smooth muscle cells" = "ACTA2",
  
  "Breast glandular cells" = "SERPINA3", 
  "Breast glandular cells" = "DSC2",
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "TNFRSF17",
  
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  "Myeloid" = "ITGAM",
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  
  "Dendritic cell" = "CD83",
  "Dendritic cell" = "FCER1A",
  'Dendritic cell' = "LILRA4",  
  'Dendritic cell' = "CLEC9A"
)


DotPlot(all_data4, features = major_marker_ls, group.by = "sketch_snn_res.0.2")
