library(tidyverse)
library(Seurat)

# conda activate R
# srun --pty /bin/bash
# conda activate R
# R


# setwd('pipeline')
getwd() # "pipeline"
all_data_ls_pass <- readRDS("../all_data_ls_pass.rds")


# 多进程
library(future)
plan("multisession", workers = 6, gc =  TRUE)

options(future.globals.maxSize = 10 * 1024 ^ 3)

nbrOfFreeWorkers()



# 常规merge
all_data <- merge(all_data_ls_pass[[1]], all_data_ls_pass[-1])

all_data <- 
  NormalizeData(
    object = all_data, 
    normalization.method = "LogNormalize", 
    scale.factor = 10000
  )

all_data <- 
  FindVariableFeatures(
    object = all_data,
    nfeatures = round(0.9 * nrow(all_data))
  )

all_data <- 
  all_data %>%
  ScaleData(
    vars.to.regress = c(
      "nCount_RNA"
    )
  )



# pca -----------
all_data <- RunPCA(all_data)


pdf('elbow.pdf', height = 5, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()


set.seed(1234)
dim_heatmap_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 0.07)

pdf("dimheatmap 1-5.pdf", height = 5, width = 10)
DimHeatmap(all_data[, dim_heatmap_cells$cell_id], dims = 1:5)
dev.off()

pdf("dimheatmap 6-10.pdf", height = 5, width = 10)
DimHeatmap(all_data[, dim_heatmap_cells$cell_id], dims = 6:10)
dev.off()

pdf("dimheatmap 11-15.pdf", height = 5, width = 10)
DimHeatmap(all_data[, dim_heatmap_cells$cell_id], dims = 11:15)
dev.off()



# dims: 6 10
# 测试harmony --------------------

# umap
all_data <- all_data %>% RunUMAP(dims = 1:6)
all_data <- all_data %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:6)
all_data <- all_data %>% RunUMAP(dims = 1:6, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")




pdf('dimplot.pdf')
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "project") + NoLegend()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "project") + NoLegend()
dev.off()


all_data <- FindNeighbors(object = all_data, reduction = "harmony", dims = 1:6)
all_data <- FindClusters(object = all_data, resolution = c(0.2, 0.4, 0.6, 0.8))


pdf('dimplot - cluster.pdf')
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()

DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()

dev.off()

saveRDS(all_data, "all_data_dim6.rds")


# 测试harmony dim10, used --------------------

# umap
all_data <- all_data %>% RunUMAP(dims = 1:10)
all_data <- all_data %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:10)
all_data <- all_data %>% RunUMAP(dims = 1:10, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")




pdf('dimplot - dim10.pdf')
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "project") + NoLegend()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "project") + NoLegend()
dev.off()


all_data <- FindNeighbors(object = all_data, reduction = "harmony", dims = 1:10)
all_data <- FindClusters(object = all_data, resolution = c(0.2, 0.4, 0.6, 0.8))


pdf('dimplot - cluster - dim10.pdf')
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()

DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()

dev.off()

saveRDS(all_data, "all_data_dim10.rds")


# annotation --------------
cell_marker <- 
  c(Epithelial = "EPCAM", Epithelial = "KRT18", Epithelial = "KRT20", 
    Fibroblasts = "ACTA2", Fibroblasts = "COL5A2", Fibroblasts = "FGF7", 
    Fibroblasts = "MME", Fibroblasts = "COL1A1", Fibroblasts = "COL1A2", 
    `Endothelial cells` = "PECAM1", `Endothelial cells` = "VWF", `Endothelial cells` = "MCAM",
    Macrophage = "CD68", Macrophage = "CD163", Macrophage = "MARCO", Macrophage = "MRC1",
    Monocytes = "CD14",  Monocytes = "FCGR3A", Monocytes = "FCN1", Monocytes = "LYZ", 
    Monocytes = "VCAN", `Dendritic cells1` = "CLEC9A", `Dendritic cells1` = "XCR1", 
    `Dendritic cells2` = "FCER1A", `Dendritic cells2` = "CD1C", `Dendritic cells2` = "CD1E", 
    `pDendritic cells` = "IRF4", `pDendritic cells` = "IRF7", `pDendritic cells` = "IRF8", 
    `pDendritic cells` = "LAMP5", `pDendritic cells` = "SCT", `pDendritic cells` = "SLC32A1", 
    `pDendritic cells` = "SHD", `pDendritic cells` = "LRRC26", `pDendritic cells` = "PACSIN1", 
    `pDendritic cells` = "LILRA4", 
    `T cells` = "CD2", `T cells` = "CD3D", 
    `T cells` = "CD3E", `T cells` = "CD4", `T cells` = "CD8A", `T cells` = "PRF1", 
    `T cells` = "GZMB", `T cells` = "CD69", `T cells` = "FOXP3", `T cells` = "TOP2A",
    `NK cells` = "NKG7", `NK cells` = "GNLY", 
    `NK cells` = "KLRB1", `NK cells` = "KLRC1", `NK cells` = "KLRD1", 
    `NK cells` = "NCR1",  `NK cells` = "FGFBP2", `NK cells` = "CX3CR1",
    `B cells` = "CD79A", `B cells` = "MS4A1", `B cells` = "CD19", `B cells` = "CXCR4", `B cells` = "CD27", 
    `B cells` = "IGHG1", `B cells` = "IGHA1", `B cells` = "IGKC", `B cells` = "MZB1",
    Granulocytes = "CPA3", Granulocytes = "CTSG", Granulocytes = "HPGDS", 
    Granulocytes = "MS4A2", Granulocytes = "VWA5A", Granulocytes = "CCR3", 
    Granulocytes = "CD33", Granulocytes = "IL5RA",
    `mast cell` = "TPSAB1", `mast cell` = "TPSB2"
  )



pdf(
  "dotplot.pdf",
  height = 10,
  width = 25
)

dot_theme <-  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) 
dot_scale <- scale_colour_gradient2(
  low = "white",
  mid = "gray90",
  high = "blue",
  midpoint = 0
)
DotPlot(all_data, features = cell_marker, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
DotPlot(all_data, features = cell_marker, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
DotPlot(all_data, features = cell_marker, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
DotPlot(all_data, features = cell_marker, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()

# anno cluster
pdf('dimplot - cluster - dim10 - split.pdf', height = 5, width = 25)
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "RNA_snn_res.0.4", split.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4", split.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
dev.off()


Idents(all_data) <- "RNA_snn_res.0.4"

major_cluster <- c(
  '0' = "Epithelial",
  '1' = "Epithelial",
  '2' = "Epithelial",
  '3' = "Endothelial",
  '4' = "Myeloid",
  '5' = "Epithelial",
  '6' = "Fibroblast",
  '7' = "Lymphocyte"
)

all_data <- RenameIdents(all_data, major_cluster)
all_data$major_annotation <- Idents(all_data)



# anno cluster
pdf('dimplot - cluster - dim10 - annotation.pdf')
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "major_annotation") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "major_annotation") + NoLegend() + ggsci::scale_color_igv()
dev.off()

