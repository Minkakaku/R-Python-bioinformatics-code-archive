library(tidyverse)
library(Seurat)
library(harmony)
source("xenium.R")


setwd('pipeline')

options(future.globals.maxSize = 100 * 1024 ^ 3)

# 普通合并 ------------------
all_data_ls_pass <- readRDS("all_data_ls_pass.rds")

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
    nfeatures = round(477 * 0.5),
    # mean.cutoff = c(0.1, 8),
    # dispersion.cutoff = c(1, Inf),
    # mean.function = ExpMean, 
    # dispersion.function = LogVMR
  )




pdf("VariableFeaturePlot.pdf")
VariableFeaturePlot(all_data)
dev.off()


all_data <- ScaleData(
  object = all_data, 
  vars.to.regress = "nCount_RNA",
  features = rownames(all_data)
)


# pca umap tsne
all_data <- RunPCA(object = all_data, features = VariableFeatures(object = all_data))

(all_data@reductions$pca@stdev) %>% {cumsum(.)/sum(.) } %>% {which(. > 0.9)}



pdf("merge data - ElbowPlot.pdf")
ElbowPlot(all_data, ndims = 50)
dev.off()



# all_data <- FindNeighbors(object = all_data, dims = 1:PCNum)
# all_data <- FindClusters(object = all_data, resolution = 1.2)


all_data <- RunUMAP(object = all_data, dims = 1:30)


pdf("merge data - umap - dim30.pdf")
DimPlot(object = all_data, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()



# all_data <- RunUMAP(object = all_data, dims = 1:50)
# 
# 
# pdf("merge data - umap - dim50.pdf")
# DimPlot(object = all_data, reduction = "umap", label = TRUE, pt.size = 0.5)
# dev.off()

# 测试harmony --------------------
all_data <- all_data %>% RunHarmony(group.by.vars = 'project')
all_data <- RunUMAP(object = all_data, dims = 1:30, reduction = "harmony")



pdf("merge data - umap - harmony - dim30.pdf")
DimPlot(object = all_data, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

saveRDS(all_data, "all_data_harmony.rds")



# 测试integrated data ----------------
all_data <- readRDS("all_data.rds")

# pdf("integrated data - umap - dim30.pdf")
# DimPlot(object = all_data, reduction = "umap", label = TRUE, pt.size = 0.5)
# dev.off()


DefaultAssay(all_data) <- "integrated"

# pca umap tsne
all_data <- ScaleData(all_data, vars.to.regress = "nCount_RNA")
all_data <- RunPCA(object = all_data, features = VariableFeatures(object = all_data))

(all_data@reductions$pca@stdev) %>% {cumsum(.)/sum(.) } %>% {which(. > 0.9)}



pdf("integrated data - ElbowPlot.pdf")
ElbowPlot(all_data, ndims = 50)
dev.off()


all_data <- RunUMAP(object = all_data, dims = 1:30)

all_data$Group <- ifelse(all_data$Tissue == "T", all_data$TP53, all_data$Tissue)

pdf("integrated data - umap - dim30.pdf")
DimPlot(object = all_data, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(object = all_data, group.by = 'Tissue', reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(object = all_data, group.by = 'TP53', reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(object = all_data, group.by = 'Group', reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

pdf("integrated data - umap - split - dim30.pdf", height = 4, width = 35)
DimPlot(object = all_data, reduction = "umap", split.by = "project", label = TRUE, pt.size = 0.5)
dev.off()


saveRDS(all_data, "all_data_dr.rds")

# 重新整合 ------------------------------------------------------------------
# 用于合并的anchor gene没有变化，但是用于合并数据集中的基因是全部基因

IntegrateReferId <-
  findReference(
    all_data_ls_pass,
    name =  c(
      "20231128__L01-001", # WT
      "20231128__R01-001", # MUT
      "20240112__L02-001"  # P
    )
  )


features <- SelectIntegrationFeatures(
  object.list = all_data_ls_pass, 
  assay = rep("RNA", times = length(all_data_ls_pass)),
  nfeatures = round(nrow(all_data_ls_pass[[1]]) * 0.5 )
)


for(i in seq_along(all_data_ls_pass) ){
  all_data_ls_pass[[i]] <- 
    all_data_ls_pass[[i]] %>%
    FindVariableFeatures(nfeatures = round(nrow(all_data_ls_pass[[1]]) * 0.5 )) %>%
    NormalizeData() %>%
    ScaleData(vars.to.regress = "nCount_RNA")
  
  all_data_ls_pass[[i]] <- all_data_ls_pass[[i]] %>% RunPCA()
}


AnchorSet <- 
  FindIntegrationAnchors(
    object.list = all_data_ls_pass,
    reference   = IntegrateReferId,
    normalization.method = 'LogNormalize', # !important
    anchor.features      = features,
    reduction            = 'rpca' #CCA往往无法运行成功
  )

all_data <- 
  IntegrateData(
    AnchorSet, 
    normalization.method = "LogNormalize", # !important
    features.to.integrate = rownames(all_data_ls_pass[[1]])# NULL
  )



# add annotation -------------------
anno_data <- read.csv('../data/Annotation.csv')
meta_anno_data <-
  all_data[['project']] %>%
  rownames_to_column("obj_id") %>%
  left_join(anno_data, by = "project") %>%
  column_to_rownames("obj_id") 
meta_anno_data$Group <- ifelse(meta_anno_data$Tissue == "T", meta_anno_data$TP53, all_data$Tissue)

all_data <- all_data %>% AddMetaData(metadata = meta_anno_data)


saveRDS(all_data, 'all_data_integrated.rds')


# all_data <- 
#   FindVariableFeatures(
#     object = all_data,
#     nfeatures = round(477 * 0.5),
#     # mean.cutoff = c(0.1, 8),
#     # dispersion.cutoff = c(1, Inf),
#     # mean.function = ExpMean, 
#     # dispersion.function = LogVMR
#   )

all_data <- ScaleData(
  object = all_data,
  vars.to.regress = "nCount_RNA"
)

# pca umap tsne
all_data <- RunPCA(object = all_data, features = VariableFeatures(object = all_data))

(all_data@reductions$pca@stdev) %>% {cumsum(.)/sum(.) } %>% {which(. > 0.9)}



pdf("integrated - 16 - data - ElbowPlot.pdf")
ElbowPlot(all_data, ndims = 50)
dev.off()


all_data <- RunUMAP(object = all_data, dims = 1:30)


pdf("integrated - 16 - data - umap - dim30.pdf")
DimPlot(object = all_data, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()


saveRDS(all_data, 'all_data_integrated_dr.rds')


all_data <- FindNeighbors(object = all_data, dims = 1:30)

library(future)
plan("multisession", workers = 6, gc = TRUE)
options(future.globals.maxSize = 100 * 1024 ^ 3)
nbrOfFreeWorkers()

all_data <- FindClusters(object = all_data, resolution = 1)

plan("sequential")
future:::ClusterRegistry("stop")

saveRDS(all_data, 'all_data_integrated_cluster_res1.rds')

all_data <- FindClusters(object = all_data, resolution = 0.2)

saveRDS(all_data, "all_data_integrated_cluster_res1_res0.2.rds")

# integrated -> harmony --------------------
all_data <- all_data %>% RunHarmony(group.by.vars = 'project')
all_data <- RunUMAP(object = all_data, dims = 1:30, reduction = "harmony")

# export
pdf("integrated - harmony - 16 - data - umap - dim30.pdf")
DimPlot(object = all_data, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

saveRDS(all_data, "all_data_integrated_harmony.rds")


# 放弃integrated -> harmony



# 细胞注释 ------------------
all_data <- readRDS("all_data_integrated_cluster_res1_res0.2.rds")


# cell_marker <- 
#   readxl::read_excel(
#     "marker gene list.xlsx",
#     sheet = "Sheet2"
#   ) %>%
#   deframe()


cell_marker <- 
  c(Epithelial = "EPCAM", Epithelial = "KRT18", Epithelial = "KRT20", 
    Fibroblasts = "ACTA2", Fibroblasts = "COL5A2", Fibroblasts = "FGF7", 
    Fibroblasts = "MME", fibroblast = "COL1A1", fibroblast = "COL1A2", 
    `Endothelial cells` = "PECAM1", `Endothelial cells` = "VWF", 
    `Endothelial cells` = "MCAM", Macrophage = "CD68", Macrophage = "CD163", 
    Macrophage = "MARCO", Macrophage = "MRC1", Monocytes = "CD14", 
    Monocytes = "FCGR3A", Monocytes = "FCN1", Monocytes = "LYZ", 
    Monocytes = "VCAN", `Dendritic cells1` = "CLEC9A", `Dendritic cells1` = "XCR1", 
    `Dendritic cells2` = "FCER1A", `Dendritic cells2` = "CD1C", `Dendritic cells2` = "CD1E", 
    `pDendritic cells` = "IRF4", `pDendritic cells` = "IRF7", `pDendritic cells` = "IRF8", 
    `pDendritic cells` = "LAMP5", `pDendritic cells` = "SCT", `pDendritic cells` = "SLC32A1", 
    `pDendritic cells` = "SHD", `pDendritic cells` = "LRRC26", `pDendritic cells` = "PACSIN1", 
    `pDendritic cells` = "LILRA4", `T cells` = "CD2", `T cells` = "CD3D", 
    `T cells` = "CD3E", `T cells` = "CD4", `T cells` = "CD8A", `T cells` = "PRF1", 
    `T cells` = "GZMB", `T cells` = "CD69", `T cells` = "FOXP3", 
    `T cells` = "TOP2A", `NK cells` = "NKG7", `NK cells` = "GNLY", 
    `NK cells` = "KLRB1", `NK cells` = "KLRC1", `NK cells` = "KLRD1", 
    `NK cells` = "NCR1",  `NK cells` = "FGFBP2", 
    `NK cells` = "CX3CR1", `B cells` = "CD79A", `B cells` = "MS4A1", 
    `B cells` = "CD19", `B cells` = "CXCR4", `B cells` = "CD27", 
    Granulocytes = "CPA3", Granulocytes = "CTSG", Granulocytes = "HPGDS", 
    Granulocytes = "MS4A2", Granulocytes = "VWA5A", Granulocytes = "CCR3", 
    Granulocytes = "CD33", Granulocytes = "IL5RA", `mast cell` = "TPSAB1", 
    `mast cell` = "TPSB2")

DefaultAssay(all_data) <- "integrated"
p_dot <- 
  DotPlot(all_data, features = cell_marker) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )


# export 
pdf(
  "integrated - 16 - data - dotplot.pdf",
  height = 15,
  width = 25
)
p_dot %>% print 
dev.off()


# "Epithelial"        "Fibroblasts"       "fibroblast"        "Endothelial cells" "Macrophage"        "Monocytes"        
# "Dendritic cells1"  "Dendritic cells2"  "pDendritic cells"  "T cells"           "NK cells"          "B cells"          
# "Granulocytes"      "mast cell"  

cell_cluster <- c(
  '0'  = 'Epithelial', 
  '1'  = 'Fibroblasts', 
  '2'  = 'Epithelial', 
  '3'  = 'T cells', 
  '4'  = 'Endothelial', 
  '5'  = 'NK cells', 
  '6'  = 'Monocytes', 
  '7'  = 'Epithelial', 
  '8'  = 'Dendritic cells', 
  '9'  = 'Granulocytes', 
  '10' = 'B cells', 
  '11' = 'NA', # "Epithelial"
  '12' = 'Fibroblasts', 
  '13' = 'B cells', 
  '14' = 'NA', # "Epithelial"
  '15' = 'Fibroblasts', 
  '16' = 'Fibroblasts', 
  '17' = 'T cells', 
  '18' = 'NA', # "Epithelial"
  '19' = 'Fibroblasts', 
  '20' = 'Fibroblasts', 
  '21' = 'Monocytes',  # NA
  '22' = 'Granulocytes', # "Epithelial"
  '23' = 'pDendritic cells', # NA
  '24' = 'Epithelial', 
  
  '25' = 'pDendritic cells', 
  '26' = 'Epithelial', 
  '27' = 'Epithelial', 
  '28' = 'NA', # "Fibroblasts"
  '29' = 'NA', # "Epithelial"
  '30' = 'T cells', 
  '31' = 'Macrophage', 
  '32' = 'NA', # 'T cells'
  '33' = 'NA', 
  '34' = 'pDendritic cells', # 'T cells' 
  '35' = 'NA', 
  '36' = 'NK', 
  '37' = 'NK', # NA
  '38' = 'Macrophage', # NA
  '39' = 'Fibroblasts', # 'T cells' 
  '40' = 'Macrophage', # NA
  '41' = 'Monocytes', # 'T cells' 
  '42' = 'Dendritic cells', 
  '43' = 'Epithelial', 
  '44' = 'Epithelial', 
  '45' = 'Monocytes', 
  '46' = 'NA', 
  '47' = 'Monocytes'
)



all_data <- eval(expr(RenameIdents(all_data, !!!cell_cluster)))

# export 
pdf(
  "integrated - 16 - cell - annotation.pdf",
  height = 6,
  width = 8
)
DimPlot(all_data) %>% print
dev.off()

pdf(
  "integrated - 16 - cell - annotation - splitby - project.pdf",
  height = 25,
  width = 12
)
DimPlot(all_data, split.by = "project", ncol = 3) +
  ggsci::scale_color_igv()
dev.off()


# merge 
DefaultAssay(all_data) <- "RNA"
all_data <- all_data %>% NormalizeData()

all_data$integrated_snn_res.0.2 <-
  factor(
    all_data$integrated_snn_res.0.2,
    levels = as.character(sort(as.numeric(unique(all_data$integrated_snn_res.0.2))))
  )
p_dot2 <- 
  DotPlot(all_data, features = cell_marker, group.by = "integrated_snn_res.0.2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - data - dotplot - using - RNA-assay.pdf",
  height = 15,
  width = 25
)
p_dot2 %>% print 
dev.off()


# 细胞注释 - 胰腺细胞 -------------------------

# cell_marker2 <-
#   readxl::read_excel(
#     "marker gene list - pancreas.xlsx"
#   ) %>%
#   deframe()

cell_marker2 <-
  c(`Ductal cells` = "CFTR", `Ductal cells` = "KRT18", `Ductal cells` = "ADAMTS1", 
    `Ductal cells` = "ALDH1A3", `Ductal cells` = "C5orf46", `Ductal cells` = "CCL5", 
    `Ductal cells` = "CFB", `Ductal cells` = "CNN1", `Ductal cells` = "CXCL10", 
    `Ductal cells` = "CXCL2", `Ductal cells` = "CXCL6", `Ductal cells` = "EPCAM", 
    `Ductal cells` = "FXYD2", `Ductal cells` = "GLIPR1", `Ductal cells` = "GPX2", 
    `Ductal cells` = "HES4", `Ductal cells` = "LIF", `Ductal cells` = "MYC", 
    `Ductal cells` = "NPDC1", `Ductal cells` = "TCIM", `Ductal cells` = "TNC", 
    `Endocrine cells` = "CHGA", `Endocrine cells` = "GKN2", `Endocrine cells` = "GPRC5A", 
    `Endocrine cells` = "IL1R2", `Endocrine cells` = "SFTA2", `Endocrine cells` = "TFF2", 
    `Exocrine glandular cells` = "AMY2A", `Exocrine glandular cells` = "DPEP1", 
    `Exocrine glandular cells` = "FSTL3", `Exocrine glandular cells` = "HAMP", 
    `Exocrine glandular cells` = "KLK11", `Alpha cells` = "GCG", 
    `Beta cells` = "INS", `Delta cells` = "SST", `Epsilon cells` = "GHRL", 
    `Gamma cells` = "PPY")

# 
DefaultAssay(all_data) <- "integrated"

all_data$integrated_snn_res.0.2 <-
  factor(
    all_data$integrated_snn_res.0.2,
    levels = as.character(sort(as.numeric(unique(as.character(all_data$integrated_snn_res.0.2)))))
  )

p_dot3 <- 
  DotPlot(all_data, features = cell_marker2, group.by = "integrated_snn_res.0.2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - data - dotplot - pdac - marker.pdf",
  height = 15,
  width = 25
)
p_dot3 %>% print 
dev.off()



# cell_cluster <- c(
#   '0'  = 'Epithelial', 
#   '2'  = 'Epithelial', 
#   '7'  = 'Epithelial', 
#   '11' = 'NA', # "Epithelial"
#   '14' = 'NA', # "Epithelial"
#   '18' = 'NA', # "Epithelial"
#   '24' = 'Epithelial', 
#   '26' = 'Epithelial', 
#   '27' = 'Epithelial', 
#   '29' = 'NA', # "Epithelial"
#   '43' = 'Epithelial', 
#   '44' = 'Epithelial'
# )


sc_epic <- all_data[, WhichCells(all_data, idents = 'Epithelial')]

DefaultAssay(sc_epic) <- "integrated"

sc_epic <-
  sc_epic %>%
  ScaleData() %>%
  RunPCA() 

sc_epic <- FindNeighbors(object = sc_epic, dims = 1:30)
sc_epic <- FindClusters(object = sc_epic, resolution = 0.2)


p_dot_epic <- 
  DotPlot(sc_epic, features = cell_marker2, group.by = "integrated_snn_res.0.2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - epic data - dotplot - pdac - marker.pdf",
  height = 12,
  width = 25
)
p_dot_epic %>% print 
dev.off()



cell_marker3 <-
  c(
  `Epithelial cell of pancreas` = "EPCAM",
  `Epithelial cell of pancreas` = "KRT19",
  `Epithelial cell of pancreas` = "KRT20",
  
  `Ductal cells` = "CFTR", 
  `Ductal cells` = "KRT18", 
  `Ductal cells` = "SPP1", 
  `Ductal cells` = "AQP1", 


  # `Endocrine cells` = "CHGA", `Endocrine cells` = "GKN2", `Endocrine cells` = "GPRC5A", 

  `Exocrine glandular cells` = "AMY2A", 
  # `Exocrine glandular cells` = "PRSS1", 
  # `Exocrine glandular cells` = "CELA3A", 
  # `Exocrine glandular cells` = "KLK1", 
  # `Exocrine glandular cells` = "PNLIP", 
  # `Exocrine glandular cells` = "CTRC", 
  # `Exocrine glandular cells` = "AKR1C3", 
  # `Exocrine glandular cells` = "REG3A", 
  # `Exocrine glandular cells` = "CTRB1", 
  
  # `Exocrine glandular cells` = "CPA1", 
  # `Exocrine glandular cells` = "CPA2", 
  # `Exocrine glandular cells` = "REG1A", 
  # `Exocrine glandular cells` = "REG1B", 
  # `Exocrine glandular cells` = "DPEP1", 
  # `Exocrine glandular cells` = "FSTL3", # fibro ductal-cells
  # `Exocrine glandular cells` = "HAMP", 
  # `Exocrine glandular cells` = "KLK11",  # ductal-cells
  
  `Alpha cells` = "GCG", 
  `Beta cells` = "INS", 
  `Delta cells` = "SST", 
  `Epsilon cells` = "GHRL", 
  `Gamma cells` = "PPY")


p_dot_epic2 <- 
  DotPlot(sc_epic, features = cell_marker3, group.by = "integrated_snn_res.0.2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - epic data - dotplot - pdac - marker2.pdf",
  height = 8,
  width = 25
)
p_dot_epic2 %>% print 
dev.off()

epic_anno <- c(
  '0'  = 'Ductal cells',
  '1'  = 'Ductal cells',
  '2'  = 'Ductal cells',
  '3'  = 'Ductal cells',
  '4'  = 'Ductal cells',
  '5'  = 'Epsilon cells',
  '6'  = 'Ductal cells',
  '7'  = 'Ductal cells',
  '8'  = 'Epsilon cells',
  '9'  = 'Ductal cells',
  '10' = 'Ductal cells',
  '11' = 'Beta cells',
  '12' = 'Ductal cells',
  '13' = 'Ductal cells',
  '14' = 'Ductal cells',
  '15' = 'Epithelial',
  '16' = 'Ductal cells',
  '17' = 'Epithelial',
  '18' = 'Alpha cells',
  '19' = 'Ductal cells',
  '20' = 'Acinar cells',
  '21' = 'Ductal cells',
  '22' = 'Epithelial',
  '23' = 'Delta cells'
)


sc_epic <- eval(expr(RenameIdents(sc_epic, !!!epic_anno)))

p_dot_epic_cell <- 
  DotPlot(sc_epic, features = cell_marker3) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - epic data - cellcluster - dotplot - pdac - marker2.pdf",
  height = 8,
  width = 25
)
p_dot_epic_cell %>% print 
dev.off()


saveRDS(sc_epic, "sc_epic.rds")



# 细胞注释 - T细胞 --------------------------------------------------------------

sc_t <- all_data[, WhichCells(all_data, idents = 'T cells')]

DefaultAssay(sc_t) <- "integrated"

sc_t <-
  sc_t %>%
  ScaleData() %>%
  RunPCA() 

sc_t <- FindNeighbors(object = sc_t, dims = 1:30)
sc_t <- FindClusters(object = sc_t, resolution = 0.2)





cell_marker_t <-
  c(
    "CD2",
    "CD3D",
    "CD3E",
    "CD4",
    "CD8A",
    
    # "PRF1",
    "GZMB",
    # "CD69",
    "FOXP3",
    # "TOP2A",
    "NKG7",
    # "GNLY", 
    # "KLRB1",
    # "KLRC1",
    # "KLRD1",
    # "NCR1",
    "FGFBP2"
  )


p_dot_t <- 
  DotPlot(sc_t, features = cell_marker_t, group.by = "integrated_snn_res.0.2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - t data - dotplot.pdf",
  height = 8,
  width = 25
)
p_dot_t %>% print 
dev.off()




t_anno <- c(
  '0'  = "CD8T",
  '1'  = "Ohters T",
  '2'  = "CD4T",
  '3'  = "CD4T",
  '4'  = "Ohters T",
  '5'  = "Tregs",
  '6'  = "CD4T",
  '7'  = "Ohters T",
  '8'  = "CD4T",
  '9'  = "Tregs",
  '10' = "NK/NKT",
  '11' = "Ohters T",
  '12' = "Tregs",
  '13' = "Ohters T",
  '14' = "Ohters T",
  '15' = "NK/NKT"
)


sc_t <- eval(expr(RenameIdents(sc_t, !!!t_anno)))

p_dot_t_cell2 <- 
  DotPlot(sc_t, features = cell_marker_t) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - t data - cellcluster - dotplot.pdf",
  height = 6,
  width = 8
)
p_dot_t_cell2 %>% print 
dev.off()


saveRDS(sc_t, "sc_t.rds")


# final annotation ###########
epic_anno_df <- sc_epic %>%  Idents %>% enframe("cell_id", "annotation_epic")
t_anno_df <- sc_t %>%  Idents %>% enframe("cell_id", "annotation_t")

# all_data$annotation <- Idents(all_data)

all_anno_df <- all_data %>% Idents %>% enframe("cell_id", "annotation_all")

annotation_df <-
  all_anno_df %>%
  left_join(epic_anno_df, by = "cell_id") %>%
  left_join(t_anno_df, by = "cell_id") %>%
  map_dfc(as.character) %>%
  mutate(
    annotation = ifelse(annotation_all == "Epithelial", annotation_epic, annotation_all),
    annotation = ifelse(annotation == "T cells", annotation_t, annotation)
  )

all_data <-
  all_data %>%
  AddMetaData(column_to_rownames(annotation_df[c('cell_id', "annotation")], "cell_id"))

all_data$annotation[all_data$annotation == "NA"] <- "Undefined"
all_data$annotation[all_data$annotation == "NK"] <- "NK/NKT"
all_data$annotation[all_data$annotation == "NK cells"] <- "NK/NKT"


Idents(all_data) <- "annotation"


saveRDS(all_data, "all_data_final.rds")



# export
pdf("integrated - 16 - data - umap - dim30 - res0.2.pdf")
DimPlot(object = all_data, reduction = "umap", group.by = "integrated_snn_res.0.2", label = TRUE, pt.size = 0.5)
dev.off()

pdf(
  "integrated - 16 - data - umap - split - dim30 - res0.2.pdf",
  height = 25,
  width = 12
)
DimPlot(all_data, reduction = "umap", group.by = "integrated_snn_res.0.2", label = TRUE, pt.size = 0.5, split.by = "project", ncol = 3) 
dev.off()

# 2 降维聚类图 -----------------------------------------------------------------
all_data$Group <- ifelse(all_data$Tissue == "T", all_data$TP53, all_data$Tissue)


pdf("integrated - 16 - data - umap - dim30.pdf")
DimPlot(
  object = all_data,
  reduction = "umap", 
  label = TRUE, 
  pt.size = 0.5,
  repel = TRUE
)+
  ggsci::scale_color_igv() +
  coord_fixed(ratio = 0.9)
DimPlot(object = all_data, group.by = 'Tissue', reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'TP53', reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'Group', reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
dev.off()

pdf("integrated - 16 - data - umap - split - dim30.pdf", height = 4, width = 35)
  DimPlot(
    object = all_data, 
    reduction = "umap",
    split.by = "project", 
    label = TRUE, 
    pt.size = 0.5,
    repel = TRUE
  ) +
    ggsci::scale_color_igv()
  DimPlot(object = all_data, group.by = 'Tissue', split.by = "project", reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
  DimPlot(object = all_data, group.by = 'TP53', split.by = "project", reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
  DimPlot(object = all_data, group.by = 'Group', split.by = "project", reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
dev.off()



# export 
pdf(
  "integrated - 16 - data - umap - dim30 - annotation.pdf",
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "annotation") +
  ggsci::scale_color_igv()
dev.off()

pdf(
  "integrated - 16 - data - umap - dim30 - project.pdf",
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "project") +
  ggsci::scale_color_igv()
dev.off()

pdf(
  "integrated - 16 - data - umap - split - dim30 - annotation.pdf",
  height = 25,
  width = 12
)
DimPlot(all_data, group.by = "annotation", split.by = "project", ncol = 3) +
  ggsci::scale_color_igv()
dev.off()

pdf(
  "integrated - 16 - data - umap - split - dim30 - project.pdf",
  height = 25,
  width = 12
)
DimPlot(all_data, group.by = "project", split.by = "project", ncol = 3) +
  ggsci::scale_color_igv()
dev.off()

# 3 spatial plot ----------------------------------------------------------

data_dir <- "../data"
anno_file <- "../data/Annotation.csv"
# ROI_anno_file <- "../data/ROI_Annotation.xlsx"
# panel_file <- "Xenium_hMulti_v1_metadata.xlsx"

# data path and project name
parse_data_path <-
  dir(data_dir, full.names = TRUE) %>% 
  map(dir, full.names = TRUE) %>%
  unlist %>% 
  
  parseProjName_xenium(pattern = "output[-]?", extract_postion = c(4, 3)) %>%
  mutate(.counts_path = file.path(.path, "cell_feature_matrix")) %>%
  identity()

# add id for parse_data_path
parse_data_path <- 
  parse_data_path %>%
  mutate(.id = seq_along(.proj))


# image w/h ratio
# used to fix width and heigth ratio in In Situ Spatial plot.
image_wh_r_ls <-
  with(parse_data_path, split(.path, .proj)) %>%
  map(function(x){
    imag_wh <-
      RBioFormats::read.metadata(
        file.path(x, "morphology_mip.ome.tif")
      ) %>%
      .[[1]] %>%
      .$coreMetadata %>%
      .[c("sizeX", "sizeY")]
    
    image_wh_r <- imag_wh$sizeX/imag_wh$sizeY
    
    # try another method
    if(length(image_wh_r) <1){
      imag_wh <-
        RBioFormats::read.metadata(
          file.path(x, "morphology_mip.ome.tif")
        ) %>%
        .$coreMetadata %>%
        .[c("sizeX", "sizeY")]
      
      image_wh_r <- imag_wh$sizeX/imag_wh$sizeY
    }
    
    # length of iamge_wh_r should equal 1
    stopifnot(length(image_wh_r) == 1)
    
    image_wh_r
  }) 

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
# cells_position_ls <- cells_position_ls[proj_name]

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
# cells_boundary_ls <- cells_boundary_ls[proj_name]


# plot
all_data_meta <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })


# export image
pdf(
  "integrated - 16 - data - spatial plot.pdf",
  height = 15,
  width  = 15
)
sp_p_ls %>% iwalk(~ {print(.x + ggtitle(.y))})
dev.off()




# 2024.3.22 细胞注释 校对 --------------------------------------------------------
##### Undefined 
# 11   14   18   28   29   32   33   35   46 

Idents(all_data) <- 'integrated_snn_res.0.2'
na_markers <-
  c(
    11,14,18
    # 28,29,32,
    # 33,35,46 
  ) %>%
  map(function(x){
    FindMarkers(
      all_data,
      ident.1 = as.character(x),
      ident.2 = NULL,
      assay = "integrated"
    )
  })


na_markers %>%
  set_names(c(11,14,18)) %>%
  map(function(x){
    x %>%
      arrange(p_val_adj, -avg_log2FC) %>%
      head(10) 
  }) %>%
  map(rownames) %>%
  map(str_c, collapse = ",") %>%
  enframe() %>% 
  unnest(value) %>%
  apply(1, function(x) str_glue("cluster_{x[1]}:{x[2]}")) %>%
  walk(cat, "\n")

# cluster_11:GCG,CHGA,INS,SST,PCSK2,SCGN,HEPACAM2,DPP6,PRG4,TM4SF4
# cluster_14:PLIN4,PEBP4,LPL,CFTR,KCNK3,AQP8,CA4,GDF15,DIRAS3,CFB
# cluster_18:UMOD,AQP8,KCNK3,GDF15,DIRAS3,CFB,DPEP1,COCH,ARFGEF3,HAMP

# 三个群均设置为Epethelial
# 另第24-47群的细胞每群均只有2个，不再做注释，命名为Undefined，后面绘图时丢弃

cell_cluster_cycle2 <- c(
  '0'  = 'Epithelial', 
  '1'  = 'Fibroblasts', 
  '2'  = 'Epithelial', 
  '3'  = 'T cells', 
  '4'  = 'Endothelial', 
  '5'  = 'NK cells', 
  '6'  = 'Monocytes', 
  '7'  = 'Epithelial', 
  '8'  = 'Dendritic cells', 
  '9'  = 'Granulocytes', 
  '10' = 'B cells', 
  '11' = 'Epithelial', # "Epithelial"
  '12' = 'Fibroblasts', 
  '13' = 'B cells', 
  '14' = 'Epithelial', # "Epithelial"
  '15' = 'Fibroblasts', 
  '16' = 'Fibroblasts', 
  '17' = 'T cells', 
  '18' = 'Epithelial', # "Epithelial"
  '19' = 'Fibroblasts', 
  '20' = 'Fibroblasts', 
  '21' = 'Monocytes',  # NA
  '22' = 'Granulocytes', # "Epithelial"
  '23' = 'pDendritic cells', # NA
  '24' = 'Epithelial', 
  
  '25' = 'NA', 
  '26' = 'NA', 
  '27' = 'NA', 
  '28' = 'NA', # "Fibroblasts"
  '29' = 'NA', # "Epithelial"
  '30' = 'NA', 
  '31' = 'NA', 
  '32' = 'NA', # 'T cells'
  '33' = 'NA', 
  '34' = 'NA', # 'T cells' 
  '35' = 'NA', 
  '36' = 'NA', 
  '37' = 'NA', # NA
  '38' = 'NA', # NA
  '39' = 'NA', # 'T cells' 
  '40' = 'NA', # NA
  '41' = 'NA', # 'T cells' 
  '42' = 'NA', 
  '43' = 'NA', 
  '44' = 'NA', 
  '45' = 'NA', 
  '46' = 'NA', 
  '47' = 'NA'
)


all_data <- eval(expr(RenameIdents(all_data, !!!cell_cluster_cycle2)))
all_data$annotation2 <- Idents(all_data)


#### ductal cells ----------------
cell_marker2_cycle2 <-
  c(`Ductal cells` = "CFTR", 
    `Ductal cells` = "GPX2",  # Acinar cells/Ductal cells
    `Ductal cells` = "HES4",  # Pancreatic stellate cells/Ductal cells
    `Ductal cells` = "FSTL3", 
    
    `Exocrine glandular cells` = "AMY2A", 
    `Exocrine glandular cells` = "DPEP1",
    
    `Endocrine cells` = "CHGA", 
    
    `Alpha cells` = "GCG", 
    `Beta cells` = "INS", 
    `Delta cells` = "SST", 
    `Epsilon cells` = "GHRL", 
    `Gamma cells` = "PPY"
  )

sc_epic_cycle2 <- all_data[, WhichCells(all_data, idents = 'Epithelial')]

DefaultAssay(sc_epic_cycle2) <- "integrated"

sc_epic_cycle2 <-
  sc_epic_cycle2 %>%
  ScaleData() %>%
  RunPCA() 

sc_epic_cycle2 <- FindNeighbors(object = sc_epic_cycle2, dims = 1:30) #, k.param = 50)
sc_epic_cycle2 <- FindClusters(object = sc_epic_cycle2, resolution = 0.2)

table(sc_epic_cycle2$integrated_snn_res.0.2)

p_dot_epic_cycle2 <- 
  DotPlot(sc_epic_cycle2, features = cell_marker2_cycle2, group.by = "integrated_snn_res.0.2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - epic data - dotplot - pdac - marker - refined - k50.pdf",
  height = 8,
  width = 16
)
p_dot_epic_cycle2 %>% print 
dev.off()

# epic umap
sc_epic_cycle2 <- 
  sc_epic_cycle2 %>%
  RunUMAP(
    dims = 1:30,
    min.dist = 0.001,
    spread = 5,
    local.connectivity = 15,
    n.neighbors = 100 
  )


# epic annotation
epic_anno_cycle2 <- c(
  '0'  = 'Ductal cells',
  '1'  = 'Ductal cells',
  '2'  = 'Ductal cells',
  '3'  = 'Endocrine cells',
  '4'  = 'Ductal cells',
  '5'  = 'Ductal cells',
  '6'  = 'Ductal cells',
  '7'  = 'Endocrine cells',
  '8'  = 'Ductal cells',
  '9'  = 'Ductal cells',
  '10' = 'Acinar cells',
  '11' = 'Ductal cells',
  '12' = 'Ductal cells',
  '13' = 'Acinar cells',
  '14' = 'Ductal cells',
  '15' = 'Acinar cells',
  '16' = 'Ductal cells',
  '17' = 'Ductal cells',
  '18' = 'Ductal cells',
  '19' = 'Endocrine cells'
)

sc_epic_cycle2 <- eval(expr(RenameIdents(sc_epic_cycle2, !!!epic_anno_cycle2)))
sc_epic_cycle2$annotation_epic <- Idents(sc_epic_cycle2)



pdf(
  "integrated - 16 - epic data - umap - pdac.pdf",
  height =6,
  width = 7
)
DimPlot(
  sc_epic_cycle2,
  group.by = "integrated_snn_res.0.2",
  reduction = "umap"
)

DimPlot(
  sc_epic_cycle2,
  group.by = "ident",
  reduction = "umap"
)
dev.off()


#### Endocrine cells ----------------
sc_endocrine <- sc_epic_cycle2[, WhichCells(sc_epic_cycle2, idents = 'Endocrine cells')]

DefaultAssay(sc_endocrine) <- "integrated"

sc_endocrine <-
  sc_endocrine %>%
  ScaleData() %>%
  RunPCA() 

sc_endocrine <- FindNeighbors(object = sc_endocrine, dims = 1:20, k.param = 50)
sc_endocrine <- FindClusters(object = sc_endocrine, resolution = 0.2)

sc_endocrine <- 
  sc_endocrine %>%
  RunUMAP(
    dims = 1:30,
    min.dist = 0.01,
    spread = 3,
    local.connectivity = 7,
    n.neighbors = 50 
  )



# plot
table(sc_endocrine$integrated_snn_res.0.2)

p_dot_endocrine_cycle2 <-
  DotPlot(sc_endocrine, features = c(    `Alpha cells` = "GCG",
                                         `Beta cells` = "INS",
                                         `Delta cells` = "SST",
                                         `Epsilon cells` = "GHRL",
                                         `Gamma cells` = "PPY"
  ), group.by = "integrated_snn_res.0.2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) 
  # scale_colour_gradient2(
  #   low = "white",
  #   mid = "gray90",
  #   high = "blue",
  #   midpoint = 0
  # )

# export
pdf(
  "integrated - 16 - epic data - dotplot - pdac - endocrine - marker - refined.pdf",
  height = 8,
  width = 12
)
p_dot_endocrine_cycle2 %>% print
DimPlot(
  sc_endocrine,
  group.by = "ident",
  label = TRUE
)
dev.off()


saveRDS(sc_endocrine, 'sc_endocrine.rds')

## T细胞 - cycle2 -------------------------
# 2024.4.7
# t_cells_id <-
#   all_data[[c('integrated_snn_res.0.2', 'annotation_cycle2')]] %>%
#   rownames_to_column("cell_id") %>%
#   mutate(annotation_cycle2 = ifelse(integrated_snn_res.0.2 == '23', "T cells", annotation_cycle2)) %>%
#   dplyr::filter(annotation_cycle2 %in% c('T cells', "NK cells")) %>%
#   dplyr::pull(cell_id)
  



  
sc_t_cycle2 <- all_data[, WhichCells(all_data, idents = c('T cells', "NK cells"))]
# sc_t_cycle2 <- all_data[, t_cells_id]

DefaultAssay(sc_t_cycle2) <- "integrated"

sc_t_cycle2 <-
  sc_t_cycle2 %>%
  ScaleData() %>%
  RunPCA() 

sc_t_cycle2 <- FindNeighbors(object = sc_t_cycle2, dims = 1:30)
sc_t_cycle2 <- FindClusters(object = sc_t_cycle2, resolution = 0.2)

table(sc_t_cycle2$integrated_snn_res.0.2)

# remove small cluster
sc_t_cycle2 <- sc_t_cycle2[, !sc_t_cycle2$integrated_snn_res.0.2 %in% c(17, 18)]
sc_t_cycle2 <-
  sc_t_cycle2 %>%
  ScaleData() %>%
  RunPCA() 

sc_t_cycle2 <- FindNeighbors(object = sc_t_cycle2, dims = 1:30)
sc_t_cycle2 <- FindClusters(object = sc_t_cycle2, resolution = 0.2)

table(sc_t_cycle2$integrated_snn_res.0.2)



cell_marker_t <-
  c(
    "CD2",
    "CD3D",
    "CD3E",
    "CD4",
    "CD8A",
    
    # "PRF1",
    "GZMB",
    # "CD69",
    "FOXP3",
    # "TOP2A",
    "NKG7",
    # "GNLY", 
    # "KLRB1",
    # "KLRC1",
    # "KLRD1",
    # "NCR1",
    "FGFBP2"
  )


p_dot_t_cycle2 <- 
    DotPlot(sc_t_cycle2, features = cell_marker_t, group.by = "integrated_snn_res.0.2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - t data - dotplot - cycle2.pdf",
  height = 8,
  width = 18
)
p_dot_t_cycle2 %>% print 
dev.off()




t_anno_cycle2 <- c(
  '0'  = "NK/NKT",
  '1'  = "NA",
  '2'  = "Tregs",
  '3'  = "NA",
  '4'  = "CD8Ts",
  '5'  = "CD4Ts",
  '6'  = "CD8Ts",
  '7'  = "Tregs",
  '8'  = "CD8Ts",
  '9'  = "Tregs",
  '10' = "NK/NKT",
  '11' = "CD8Ts",
  '12' = "NA",
  '13' = "Tregs",
  '14' = "NA",
  '15' = "NA",
  '16' = "CD8Ts",
  '17' = "CD4Ts",
  '18' = "NA",
  '19' = "NA",
  '20' = "NA",
  '21' = "NA",
  '22' = "NA",
  '23' = "NA",
  '24' = "CD8Ts",
  '25' = "NA"
)


Idents(sc_t_cycle2) <- "integrated_snn_res.0.2"
sc_t_cycle2 <- eval(expr(RenameIdents(sc_t_cycle2, !!!t_anno_cycle2)))

p_dot_t_cell_cycle2 <- 
  DotPlot(sc_t_cycle2, features = cell_marker_t) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - t data - cellcluster - dotplot - cycle2.pdf",
  height = 6,
  width = 8
)
p_dot_t_cell_cycle2 %>% print 
dev.off()


# saveRDS(sc_t_cycle2, "sc_t_cycle2.rds")


## monocyte亚分群 -------------------------------------------------------------
# 2024.4.7 尝试分群

sc_mono <- all_data[, WhichCells(all_data, idents = c("Monocytes", "Dendritic cells", "pDendritic cells"))]

DefaultAssay(sc_mono) <- "integrated"

sc_mono <-
  sc_mono %>%
  ScaleData() %>%
  RunPCA() 

sc_mono <- FindNeighbors(object = sc_mono, dims = 1:30)
sc_mono <- FindClusters(object = sc_mono, resolution = 0.2)

table(sc_mono$integrated_snn_res.0.2)

# # remove small cluster
# sc_mono <- sc_mono[, !sc_mono$integrated_snn_res.0.2 %in% c(17, 18)]
# sc_mono <-
#   sc_mono %>%
#   ScaleData() %>%
#   RunPCA() 
# 
# sc_mono <- FindNeighbors(object = sc_mono, dims = 1:30)
# sc_mono <- FindClusters(object = sc_mono, resolution = 0.2)
# 
# table(sc_mono$integrated_snn_res.0.2)


cell_marker_mye <-
  c(
    # 'T' = c('CD4', 'CD8A', 'FOXP3', 'NKG7'),
    # 'B' = c('CD79A'),
    'cDC' =  c('FCER1A','CST3', 'CD1C', "CD1E", "CCL19", "LAMP3", "DNASE1L3", "HLA-DQB2"),
    'pDC' =  c('IL3RA','GZMB','SERPINF1','ITM2C'),
    'Monocyte'  =  c('CD14','LYZ','S100A8','S100A9','LST1', 'FCN1', 'FCGR3A'),
    'Macrophage' =  c('CD86','CD68', 'CD163', 'MARCO', 'MRC1'),
    'Proliferative signal' =  c('MKI67','TOP2A','STMN1')
  )


# umap
sc_mono <- 
  sc_mono %>%
  RunUMAP(
    dims = 1:30,
    min.dist = 0.1,
    spread = 2,
    # local.connectivity = 7,
    n.neighbors = 20
  )

ggsave(
  DimPlot(sc_mono, group.by = "integrated_snn_res.0.2", label = TRUE),
  file = "integrated - 16 - monocyte data - dimplot - cycle2.pdf"
)

p_dot_mye_cycle2 <- 
  DotPlot(sc_mono, features = cell_marker_mye, group.by = "integrated_snn_res.0.2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - monocyte data - dotplot - cycle2.pdf",
  height = 8,
  width = 25
)
p_dot_mye_cycle2 %>% print 
dev.off()

mono_anno_cycle2 <- c(
  '0'  = "cDCs", # cDC2
  '1'  = "Macrophage",
  '2'  = "NA", # NKT
  '3'  = "cDCs",
  '4'  = "NA",
  '5'  = "NA",
  '6'  = "Monocyte",
  '7'  = "NA", # B
  '8'  = "cDCs",
  '9'  = "NA"
)


Idents(sc_mono) <- "integrated_snn_res.0.2"
sc_mono <- eval(expr(RenameIdents(sc_mono, !!!mono_anno_cycle2)))

p_dot_mye_anno_cycle2 <- 
  DotPlot(sc_mono, features = cell_marker_mye) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - monocyte data - cellcluster - dotplot - cycle2.pdf",
  height = 6,
  width = 25
)
p_dot_mye_anno_cycle2 %>% print 
dev.off()


# saveRDS(sc_mono, "sc_mono_cycle2.rds")


## macropage ----------------------
# 2024.4.8 尝试分群

sc_mac <- sc_mono[, WhichCells(sc_mono, idents = c("Macrophage"))]
DefaultAssay(sc_mac) <- "integrated"

sc_mac <-
  sc_mac %>%
  ScaleData() %>%
  RunPCA() 

sc_mac <- FindNeighbors(object = sc_mac, dims = 1:20)
sc_mac <- FindClusters(object = sc_mac, resolution = 0.2)

table(sc_mac$integrated_snn_res.0.2)

# # only use macphage
# sc_mac <- sc_mac[, sc_mac$integrated_snn_res.0.2 %in% c("0")]
# sc_mac <-
#   sc_mac %>%
#   ScaleData() %>%
#   RunPCA() 
# 
# sc_mac <- FindNeighbors(object = sc_mac, dims = 1:10)
# sc_mac <- FindClusters(object = sc_mac, resolution = 1)


cell_marker_mac <-
  c(
    # 'T' = c('CD4', 'CD8A', 'FOXP3', 'NKG7'),
    # 'B' = c('CD79A'),
    # 'cDC' =  c('FCER1A','CST3', 'CD1C', "CD1E", "CCL19", "LAMP3", "DNASE1L3", "HLA-DQB2"),
    # 'pDC' =  c('IL3RA','GZMB','SERPINF1','ITM2C'),
    # 'Monocyte'  =  c('CD14','LYZ','S100A8','S100A9','LST1', 'FCN1', 'FCGR3A'),
    'Macrophage' =  c('CD86','CD68', 'CD163', 'MARCO', 'MRC1')
    # 'Proliferative signal' =  c('MKI67','TOP2A','STMN1')
  )


# # umap
# sc_mac <- 
#   sc_mac %>%
#   RunUMAP(
#     dims = 1:30,
#     min.dist = 0.1,
#     spread = 2,
#     # local.connectivity = 7,
#     n.neighbors = 20
#   )
# 
# ggsave(
#   DimPlot(sc_mac, group.by = "integrated_snn_res.0.2", label = TRUE),
#   file = "integrated - 16 - macropahge data - dimplot - cycle2.pdf"
# )

p_dot_mac_cycle2 <- 
  DotPlot(sc_mac, features = cell_marker_mac, group.by = "integrated_snn_res.0.2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

# export 
pdf(
  "integrated - 16 - macrophage data - dotplot - cycle2.pdf",
  height = 8,
  width = 25
)
p_dot_mac_cycle2 %>% print 
dev.off()




## final annotation - 2024.3.24 ------------
epic_anno_df <- sc_epic_cycle2 %>%  Idents %>% enframe("cell_id", "annotation_epic")
# t_anno_df <- sc_t %>%  Idents %>% enframe("cell_id", "annotation_t")

# all_data$annotation <- Idents(all_data)

all_anno_df <- all_data %>% Idents %>% enframe("cell_id", "annotation_all")

annotation_df <-
  all_anno_df %>%
  left_join(epic_anno_df, by = "cell_id") %>%
  # left_join(t_anno_df, by = "cell_id") %>%
  map_dfc(as.character) %>%
  mutate(
    annotation_cycle2 = ifelse(annotation_all == "Epithelial", annotation_epic, annotation_all),
    # annotation = ifelse(annotation == "T cells", annotation_t, annotation)
  )

all_data <-
  all_data %>%
  AddMetaData(column_to_rownames(annotation_df[c('cell_id', "annotation_cycle2")], "cell_id"))

all_data$annotation_cycle2[all_data$annotation_cycle2 == "NK"] <- "NA"

Idents(all_data) <- "annotation_cycle2"

# saveRDS(all_data, "all_data_final.rds")

## final annotation - 2024.4.8 ------------
# epic_anno_df <- sc_epic_cycle2 %>%  Idents %>% enframe("cell_id", "annotation_epic")
epic_anno_df <- sc_epic_cycle2[["annotation_epic"]] %>% rownames_to_column("cell_id")
t_anno_df <- sc_t_cycle2 %>%  Idents %>% enframe("cell_id", "annotation_t")
mono_anno_df <- sc_mono %>%  Idents %>% enframe("cell_id", "annotation_mono")

# all_data$annotation <- Idents(all_data)

# all_anno_df <- all_data %>% Idents %>% enframe("cell_id", "annotation_all")
all_data$annotation_all <- all_data$annotation2
all_data$annotation_all[all_data$annotation_all == "NK"] <- "NK cells"
all_anno_df <- all_data[["annotation_all"]] %>% rownames_to_column("cell_id") 

# t在进行呀分群时中有少数细胞删除了，
# 因此annotation_cycle2_v2有几个细胞是NA
# 无法注释时的注释名称是"NA"
annotation_df <-
  all_anno_df %>%
  left_join(epic_anno_df, by = "cell_id") %>%
  left_join(t_anno_df, by = "cell_id") %>%
  left_join(mono_anno_df, by = "cell_id") %>%
  map_dfc(as.character) %>%
  mutate(
    annotation_cycle2_v2 = ifelse(annotation_all == "Epithelial", annotation_epic, annotation_all),
    annotation_cycle2_v2 = ifelse(annotation_all %in% c('T cells', "NK cells"), annotation_t, annotation_cycle2_v2),
    annotation_cycle2_v2 = ifelse(annotation_all %in% c("Monocytes", "Dendritic cells", "pDendritic cells"), annotation_mono, annotation_cycle2_v2)
  )

all_data <-
  all_data %>%
  AddMetaData(column_to_rownames(annotation_df[c('cell_id', "annotation_cycle2_v2")], "cell_id"))
Idents(all_data) <- "annotation_cycle2_v2"

all_data$annotation_final <- Idents(all_data)

# saveRDS(all_data, "all_data_final2.rds")



## drop NA -----------------------------------------------------------------
# 2024.4.9
not_na_cells <- all_data$annotation_final != 'NA'
  
all_data2 <- 
  CreateSeuratObject(
    counts = all_data@assays$integrated@data[, not_na_cells],
    meta.data = all_data[[]] %>% .[not_na_cells, ],
    assay = "integrated"
  )

DefaultAssay(all_data2) <- "integrated"
Idents(all_data2) <- "annotation_final"


all_data2 <- 
  FindVariableFeatures(
    object = all_data2,
    nfeatures = round(477 * 0.5),
    # mean.cutoff = c(0.1, 8),
    # dispersion.cutoff = c(1, Inf),
    # mean.function = ExpMean, 
    # dispersion.function = LogVMR
  )


all_data2 <- ScaleData(
  object = all_data2, 
  vars.to.regress = "nCount_RNA"
  # features = rownames(all_data2)
)


# pca umap tsne
all_data2 <- RunPCA(object = all_data2, features = VariableFeatures(object = all_data2))

(all_data2@reductions$pca@stdev) %>% {cumsum(.)/sum(.) } %>% {which(. > 0.9)}

all_data2 <- RunUMAP(
  object = all_data2, 
  dims = 1:50,
  min.dist = 0.01,
  spread = 2,
  local.connectivity = 5,
  n.neighbors = 50 
)


all_data2$annotation_final <- factor(as.character(all_data2$annotation_final))


# pdf('test.pdf')
# DimPlot(all_data2)
# dev.off()

saveRDS(all_data2, "../material_for_analysis/all_data2.rds")

all_data2 %>%
  Seurat::DietSeurat(
    .,
    scale.data = FALSE,
    assays = names(.@assays),
    dimreducs = Reductions(.)
  ) %>%
  saveRDS("../material_for_analysis/all_data2_diet.rds")


## export R objects ---------------
all_data_meta <- all_data2[[]]
saveRDS(all_data_meta, "../material_for_analysis/all_data_meta.rds")

saveRDS(parse_data_path, "../material_for_analysis/parse_data_path.rds")

avg_expr <- 
  AverageExpression(all_data2, slot = "data")
saveRDS(avg_expr, "../material_for_analysis/avg_expr.rds")

norm_data <- all_data2@assays$integrated@data
saveRDS(norm_data, "../material_for_analysis/norm_data.rds")


saveRDS(image_wh_r_ls, "../material_for_analysis/image_wh_r_ls.rds")
saveRDS(cells_boundary_ls, "../material_for_analysis/cells_boundary_ls.rds")
saveRDS(cells_position_ls, "../material_for_analysis/cells_position_ls.rds")


# all cells
cell_colors <- c(
  'Acinar cells'     =  "#FB117E",
  'Ductal cells'     =  "#D1207C",
  'Endocrine cells'  =  "#BB1117",
   
  'B cells'          =  "#4180BD",
  'CD4Ts'            =  "#21B746",
  'CD8Ts'            =  "#006B14",
  'Tregs'            =  "#5F8D66",
  'NK/NKT'           =  "#9329AE", 
  
  'cDCs'             =  "#954200",

  'Granulocytes'     =  "#F59F0C",
  'Macrophage'       =  "#EF7111",
  'Monocyte'         =  "#F9CC01",
  
  'Endothelial'      =  "#D56683",
  'Fibroblasts'      =  "#008E8C"
)


## export dimplot cycle2 - 2024.4.9 ------------------

dir_20240409 <- '2024.4.9 cycle2'
dir.create(dir_20240409)

pdf(file.path(dir_20240409, "integrated - 16 - data - umap - dim30 - cycle2.pdf"))
DimPlot(
  object = all_data2,
  reduction = "umap", 
  label = TRUE, 
  pt.size = 0.5,
  repel = TRUE
)+
  ggsci::scale_color_igv() +
  coord_fixed(ratio = 0.9)
DimPlot(object = all_data2, group.by = 'integrated_snn_res.0.2', reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data2, group.by = 'Tissue', reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data2, group.by = 'TP53', reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data2, group.by = 'Group', reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
dev.off()

pdf(file.path(dir_20240409, "integrated - 16 - data - umap - split - dim30.pdf"), height = 25, width = 12)
DimPlot(
  object = all_data2, 
  reduction = "umap",
  split.by = "project", 
  label = TRUE, 
  ncol = 3,
  pt.size = 0.5,
  repel = TRUE
) +
  ggsci::scale_color_igv()
DimPlot(object = all_data2, group.by = 'integrated_snn_res.0.2', split.by = "project", ncol = 3, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data2, group.by = 'Tissue', split.by = "project", reduction = "umap", ncol = 3, label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data2, group.by = 'TP53', split.by = "project", reduction = "umap", ncol = 3, label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data2, group.by = 'Group', split.by = "project", reduction = "umap", ncol = 3, label = TRUE, pt.size = 0.5, repel = TRUE)
dev.off()



## featureplot - marker gene ---------------
# 2024.4.10

p_feature_ls <-
  c(
    "EPCAM", "ACTA2", "PECAM1", "FCN1",
    "CD3D",  "CD3E",  "NKG7", "CD79A", "MS4A1",
    "CPA3",
    "CFTR", "AMY2A", "CHGA"
  ) %>% 
  map(function(x){
    FeaturePlot(all_data2, features = x) +
      scale_colour_gradient2(
        low = "white",
        mid = "gray90",
        high = "blue",
        midpoint = 0
      )
  })

pdf(
  file.path(dir_20240409, "integrated - 16 - data - featureplot - dim30 - cycle2.pdf")
)
p_feature_ls %>% walk(print)
dev.off()



## DoHeatmap ----------------------------------------------------
data_do_heatmap <-
  all_data2[c(
    "EPCAM", "ACTA2", "PECAM1", "FCN1",
    "CD3D",  "CD3E",  "NKG7", "CD79A", "MS4A1",
    "CPA3",
    "CFTR", "AMY2A", "CHGA"
  ), ]

data_do_heatmap <- 
  data_do_heatmap %>%
  ScaleData(features = rownames(data_do_heatmap))


ggsave(
  filename = "test doheatmap.pdf",
  height = 8,
  width  = 30,
  plot = DoHeatmap(
    data_do_heatmap,
    features = c(
      "ACTA2",
      "CD79A",
      "MS4A1",
      "CHGA",
      "PECAM1",
      "CD3D", 
      "CD3E",
      "CFTR", 
      "AMY2A", 
      "FCN1",
      "CPA3",
      "NKG7",
      "EPCAM"
    )
  )
)


scaled_expr_dat <-
  data_do_heatmap@assays$integrated@scale.data %>%
  apply(1, function(x) tapply(x, Idents(data_do_heatmap), mean))


marker_order <- c(
  "ACTA2",
  "PECAM1",
  "EPCAM",
  "CFTR", 
  "CHGA",
  "AMY2A", 
  
  "CD79A",
  "MS4A1",
  
  "CD3D", 
  "CD3E",
  "NKG7",
  
  "FCN1",
  
  "CPA3"
)
cell_order <- c(
  'Fibroblasts',
  'Endothelial', 
  
  'Ductal cells',
  'Endocrine cells',
  'Acinar cells',
  
  'B cells',
  
  'CD8Ts',
  'CD4Ts',
  'Tregs',
  'NK/NKT',
  
  'Macrophage',
  'Monocyte',
  'cDCs',
  
  'Granulocytes'
)

pdf(
  file.path(dir_20240409, "integrated - 16 - data - doHeatmap - dim30 - cycle2.pdf"),
  width = 10,
  height = 5
)
pheatmap::pheatmap(
  scaled_expr_dat[cell_order, marker_order],
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  border_color = NA
)
dev.off()


### test export data ----- data.table读写很强大 ----------
system.time(
  sc_norm[1:5, 1:1000] %>%
  as.data.frame() %>%
  rownames_to_column('Gene') %>%
  write.table(file = str_glue("test1.txt"), sep = "\t", quote=F, row.names = F)
)

# # sc_norm[1:5, 1:5] %>%
# #   as.data.frame() %>%
# #   write.table(file = str_glue("test2.txt"), sep = "\t", quote=F, row.names = T)
# 

system.time(
  sc_norm[1:5, 1:1000] %>%
  as.data.frame() %>%
  rownames_to_column('Gene') %>%
  write_delim(file = str_glue("test3.txt"), delim = "\t",  quote = "none", num_threads = 2)
)

system.time(
  sc_norm[1:5, 1:1000] %>%
    as.data.frame() %>%
    rownames_to_column('Gene') %>%
    data.table::fwrite(file = str_glue("test4.txt"), sep = "\t",  quote = FALSE, nThread = 2)
)



sc_meta <- readRDS("../material_for_analysis/all_data_meta.rds")
sc_norm <- readRDS('../material_for_analysis/norm_data.rds')
sc_meta$Group %>% 
  unique() %>% 
  walk(function(group){#browser()
    new_sc_meta <- sc_meta[sc_meta$Group == group, ]
    select_cell <- new_sc_meta %>% rownames()
    
    # str(sc_norm)
    # new_sc_norm <- sc_norm[, select_cell] %>%
    #   as.data.frame() %>%
    #   rownames_to_column('Gene')
    # str(new_sc_norm)
    
    print(group)
    dir.create(str_glue('Cellphone/V1/{group}'), recursive = T)
    sc_norm[, select_cell] %>%
      as.data.frame() %>%
      rownames_to_column('Gene') %>% 
      data.table::fwrite(file = str_glue("Cellphone/V1/{group}/counts.txt"), sep = "\t",  quote = FALSE, nThread = 5)
    
    new_sc_meta %>% 
      as.data.frame() %>%
      rownames_to_column('Cell') %>%
      dplyr::select(Cell, annotation_final) %>%
      data.table::fwrite(file = str_glue("Cellphone/V1/{group}/meta.txt"), sep = "\t",  quote = FALSE, nThread = 2)
  })






## export dimplot cycle2 - 不完全 ------------------
dir_20240324 <- '2024.3.24 cycle2'
dir.create(dir_20240324)

pdf(file.path(dir_20240324, "integrated - 16 - data - umap - dim30 - cycle2.pdf"))
DimPlot(
  object = all_data,
  reduction = "umap", 
  label = TRUE, 
  pt.size = 0.5,
  repel = TRUE
)+
  ggsci::scale_color_igv() +
  coord_fixed(ratio = 0.9)
DimPlot(object = all_data, group.by = 'integrated_snn_res.0.2', reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'Tissue', reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'TP53', reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'Group', reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
dev.off()

pdf(file.path(dir_20240324, "integrated - 16 - data - umap - split - dim30.pdf"), height = 25, width = 12)
DimPlot(
  object = all_data, 
  reduction = "umap",
  split.by = "project", 
  label = TRUE, 
  ncol = 3,
  pt.size = 0.5,
  repel = TRUE
) +
  ggsci::scale_color_igv()
DimPlot(object = all_data, group.by = 'integrated_snn_res.0.2', split.by = "project", ncol = 3, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'Tissue', split.by = "project", reduction = "umap", ncol = 3, label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'TP53', split.by = "project", reduction = "umap", ncol = 3, label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'Group', split.by = "project", reduction = "umap", ncol = 3, label = TRUE, pt.size = 0.5, repel = TRUE)
dev.off()



# 2024.5.13 细胞注释 - cycle3 -------------------------------------------------

all_data2 <- readRDS("../material_for_analysis/all_data2.rds")

# # pca umap tsne
# all_data2 <- RunPCA(object = all_data2, features = VariableFeatures(object = all_data2))
pdf('cycle3 - elbow.pdf', height = 5, width = 8)
ElbowPlot(all_data2, ndims = 30)
dev.off()


set.seed(1234)
dim_heatmap_cells <- 
  all_data2[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(annotation_final) %>%
  slice_sample(prop = 0.02)

pdf("cycle3 -  dimheatmap 1-5.pdf", height = 5, width = 10)
DimHeatmap(all_data2[, dim_heatmap_cells$cell_id], dims = 1:5)
dev.off()

pdf("cycle3 -  dimheatmap 6-10.pdf", height = 5, width = 10)
DimHeatmap(all_data2[, dim_heatmap_cells$cell_id], dims = 6:10)
dev.off()

pdf("cycle3 -  dimheatmap 11-15.pdf", height = 5, width = 10)
DimHeatmap(all_data2[, dim_heatmap_cells$cell_id], dims = 11:15)
dev.off()


(all_data2@reductions$pca@stdev) %>% {cumsum(.)/sum(.) } %>% {which(. > 0.9)}

all_data2 <- RunUMAP(
  object = all_data2, 
  dims = 1:11,
  min.dist = 0.01,
  spread = 2,
  local.connectivity = 5,
  n.neighbors = 50
)


pdf('test.pdf')
DimPlot(all_data2, label = TRUE) + NoLegend()
dev.off()


cell_marker <- 
  c(Epithelial = "EPCAM", Epithelial = "KRT18", Epithelial = "KRT20", 
    Fibroblasts = "ACTA2", Fibroblasts = "COL5A2", Fibroblasts = "FGF7", 
    Fibroblasts = "MME", fibroblast = "COL1A1", fibroblast = "COL1A2", 
    `Endothelial cells` = "PECAM1", `Endothelial cells` = "VWF", 
    `Endothelial cells` = "MCAM", Macrophage = "CD68", Macrophage = "CD163", 
    Macrophage = "MARCO", Macrophage = "MRC1", Monocytes = "CD14", 
    Monocytes = "FCGR3A", Monocytes = "FCN1", Monocytes = "LYZ", 
    Monocytes = "VCAN", `Dendritic cells1` = "CLEC9A", `Dendritic cells1` = "XCR1", 
    `Dendritic cells2` = "FCER1A", `Dendritic cells2` = "CD1C", `Dendritic cells2` = "CD1E", 
    `pDendritic cells` = "IRF4", `pDendritic cells` = "IRF7", `pDendritic cells` = "IRF8", 
    `pDendritic cells` = "LAMP5", `pDendritic cells` = "SCT", `pDendritic cells` = "SLC32A1", 
    `pDendritic cells` = "SHD", `pDendritic cells` = "LRRC26", `pDendritic cells` = "PACSIN1", 
    `pDendritic cells` = "LILRA4", `T cells` = "CD2", `T cells` = "CD3D", 
    `T cells` = "CD3E", `T cells` = "CD4", `T cells` = "CD8A", `T cells` = "PRF1", 
    `T cells` = "GZMB", `T cells` = "CD69", `T cells` = "FOXP3", 
    `T cells` = "TOP2A", `NK cells` = "NKG7", `NK cells` = "GNLY", 
    `NK cells` = "KLRB1", `NK cells` = "KLRC1", `NK cells` = "KLRD1", 
    `NK cells` = "NCR1",  `NK cells` = "FGFBP2", 
    `NK cells` = "CX3CR1", `B cells` = "CD79A", `B cells` = "MS4A1", 
    `B cells` = "CD19", `B cells` = "CXCR4", `B cells` = "CD27", 
    Granulocytes = "CPA3", Granulocytes = "CTSG", Granulocytes = "HPGDS", 
    Granulocytes = "MS4A2", Granulocytes = "VWA5A", Granulocytes = "CCR3", 
    Granulocytes = "CD33", Granulocytes = "IL5RA", `mast cell` = "TPSAB1", 
    `mast cell` = "TPSB2")

DefaultAssay(all_data2) <- "integrated"
p_dot <- 
  DotPlot(all_data2, features = cell_marker, group.by = "integrated_snn_res.1") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )


# export 
pdf(
  "cycle3 - integrated - 16 - data - dotplot.pdf",
  height = 15,
  width = 25
)
p_dot %>% print 
dev.off()

all_data2 <- FindNeighbors(object = all_data2, dims = 1:11)
all_data2 <- FindClusters(object = all_data2, resolution = c(0.2, 0.4))
all_data2 <- FindClusters(object = all_data2, resolution = c(0.6))
all_data2 <- FindClusters(object = all_data2, resolution = c(0.8))

table(all_data2$integrated_snn_res.0.2)
table(all_data2$integrated_snn_res.0.4)
table(all_data2$integrated_snn_res.0.6)
table(all_data2$integrated_snn_res.0.8)


pdf(
  "cycle3 - integrated - 16 - data - dimplot.pdf",
  height = 6,
  width = 7
)
DimPlot(all_data2, label = TRUE, group.by = "integrated_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data2, label = TRUE, group.by = "integrated_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data2, label = TRUE, group.by = "integrated_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data2, label = TRUE, group.by = "integrated_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()
dev.off()



pdf(
  "cycle3 - integrated - 16 - data - dotplot.pdf",
  height = 15,
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
DotPlot(all_data2, features = cell_marker, group.by = "integrated_snn_res.0.2") + dot_theme + dot_scale
DotPlot(all_data2, features = cell_marker, group.by = "integrated_snn_res.0.4") + dot_theme + dot_scale
DotPlot(all_data2, features = cell_marker, group.by = "integrated_snn_res.0.6") + dot_theme + dot_scale
DotPlot(all_data2, features = cell_marker, group.by = "integrated_snn_res.0.8") + dot_theme + dot_scale
  
dev.off()

### 使用RNA assay进行dotplot绘制
all_data <- readRDS('all_data_integrated.rds')
all_data <- all_data[, colnames(all_data2)] # 少量细胞被去除了
stopifnot(all(colnames(all_data) == colnames(all_data2)))

all_data2@assays[['RNA']] <- all_data@assays[['RNA']]

pdf(
  "cycle3 - integrated - 16 - data - dotplot - RNA - slot.pdf",
  height = 15,
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
DotPlot(all_data2, assay = 'RNA', features = cell_marker, group.by = "integrated_snn_res.0.2") + dot_theme + dot_scale
DotPlot(all_data2, assay = 'RNA', features = cell_marker, group.by = "integrated_snn_res.0.4") + dot_theme + dot_scale
DotPlot(all_data2, assay = 'RNA', features = cell_marker, group.by = "integrated_snn_res.0.6") + dot_theme + dot_scale
DotPlot(all_data2, assay = 'RNA', features = cell_marker, group.by = "integrated_snn_res.0.8") + dot_theme + dot_scale

dev.off()


saveRDS(all_data2, '../all_data2_cycle3.rds')

# library(future)
# plan("multisession", workers = 6, gc = TRUE)
# options(future.globals.maxSize = 100 * 1024 ^ 3)
# nbrOfFreeWorkers()


# plan("sequential")
# future:::ClusterRegistry("stop")



## 2024.5.14 harmony  ------------------------------------------------------

all_data <- readRDS("all_data_harmony.rds")

# add annotation 
# anno_data <- read.csv('../data/Annotation.csv')
# meta_anno_data <-
#   all_data[['project']] %>%
#   rownames_to_column("obj_id") %>%
#   left_join(anno_data, by = "project") %>%
#   column_to_rownames("obj_id") 
all_data$Group <- ifelse(all_data$Tissue == "T", all_data$TP53, all_data$Tissue)


# normalize
all_data <- 
  all_data %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = 'vst', nfeatures = nrow(.) * 0.9) %>%
  ScaleData(
    vars.to.regress = c("nCount_RNA")
  )


# PCA
all_data <- all_data %>% RunPCA()

pdf("cycle3 - harmony - qc pca.pdf", height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()


set.seed(1234)
dim_heatmap_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 0.02)

pdf("cycle3 - harmony -  dimheatmap 1-5.pdf", height = 5, width = 10)
DimHeatmap(all_data[, dim_heatmap_cells$cell_id], dims = 1:5)
dev.off()

pdf("cycle3 - harmony -  dimheatmap 6-10.pdf", height = 5, width = 10)
DimHeatmap(all_data[, dim_heatmap_cells$cell_id], dims = 6:10)
dev.off()

pdf("cycle3 - harmony -  dimheatmap 11-15.pdf", height = 5, width = 10)
DimHeatmap(all_data[, dim_heatmap_cells$cell_id], dims = 11:15)
dev.off()


# dims: 8 or 13

# umap
all_data <- all_data %>% RunUMAP(dims = 1:8)
all_data <- all_data %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:8)
all_data <- all_data %>% RunUMAP(dims = 1:8, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")



pdf('test.pdf')
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "Group") + NoLegend()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP") + NoLegend()
dev.off()


all_data <- FindNeighbors(object = all_data, reduction = "harmony", dims = 1:8)
all_data <- FindClusters(object = all_data, resolution = c(0.4))

saveRDS(all_data, "all_data_harmony_dim8")

pdf(
  "cycle3 - harmony - integrated - 16 - data - dimplot.pdf",
  height = 6,
  width = 7
)
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
dev.off()


pdf(
  "cycle3 - harmony - integrated - 16 - data - dimplot - split.pdf",
  height = 20,
  width = 15
)
DimPlot(all_data, label = TRUE, reduction = "umap", split.by = "RNA_snn_res.0.4", ncol = 3) + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", split.by = "RNA_snn_res.0.4", ncol = 3) + NoLegend() + ggsci::scale_color_igv()
dev.off()



pdf(
  "cycle3 - harmony - integrated - 16 - data - dotplot - RNA - slot.pdf",
  height = 15,
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
DotPlot(all_data, assay = 'RNA', features = cell_marker, group.by = "RNA_snn_res.0.4") + dot_theme + dot_scale
dev.off()


## 2024.5.14 important: harmony dim 13  ------------
## major annotation: 使用的dim13 ---------------
# umap
all_data <- all_data %>% RunUMAP(dims = 1:13)
all_data <- all_data %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:13)
all_data <- all_data %>% RunUMAP(dims = 1:13, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")
all_data <- all_data %>%
  RunUMAP(
    dims = 1:13, 
    min.dist = 0.01,
    spread = 2,
    reduction = "harmony", 
    reduction.name = "harmony_uamp", 
    reduction.key = "harmonyumap_"
  )



pdf('cycle3 - harmony - integrated - 16 - data - dimplot -  dim 13 - project.pdf')
DimPlot(all_data, label = TRUE, reduction = "umap", group.by = "project") + NoLegend()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "project") + NoLegend()
dev.off()


all_data <- FindNeighbors(object = all_data, reduction = "harmony", dims = 1:13)
all_data <- FindClusters(object = all_data, resolution = c(0.2, 0.4, 0.6, 0.8))


# dimplot
pdf(
  "cycle3 - harmony - integrated - 16 - data - dimplot -  dim 13.pdf",
  height = 6,
  width = 7
)
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()
dev.off()

pdf(
  "cycle3 - harmony - integrated - 16 - data - dimplot -  dim 13 - spread2 - min_dist0.1.pdf",
  height = 6,
  width = 7
)
DimPlot(all_data, label = TRUE, reduction = "harmony_uamp", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmony_uamp", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmony_uamp", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmony_uamp", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()
dev.off()

# dotplot
pdf(
  "cycle3 - harmony - integrated - 16 - data - dotplot - RNA - slot - dim13.pdf",
  height = 15,
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
DotPlot(all_data, assay = 'RNA', features = cell_marker, group.by = "RNA_snn_res.0.2") + dot_theme + dot_scale
DotPlot(all_data, assay = 'RNA', features = cell_marker, group.by = "RNA_snn_res.0.4") + dot_theme + dot_scale
DotPlot(all_data, assay = 'RNA', features = cell_marker, group.by = "RNA_snn_res.0.6") + dot_theme + dot_scale
DotPlot(all_data, assay = 'RNA', features = cell_marker, group.by = "RNA_snn_res.0.8") + dot_theme + dot_scale

dev.off()

# cell_annotation
cell_annotation_res.0.4 <- c(
  '0'  = 'Epithelial cells', 
  '1'  = 'Myeloid cells', 
  '2'  = 'Fibroblasts', 
  '3'  = 'T/NK cells', 
  '4'  = 'Epithelial cells', 
  '5'  = 'Epithelial cells', 
  '6'  = 'Endothelial cells', 
  '7'  = 'Granulocytes', 
  '8'  = 'Fibroblasts', 
  '9'  = 'Epithelial cells', 
  '10' = 'B cells', 
  '11' = 'Epithelial cells'
)

Idents(all_data) <- "RNA_snn_res.0.4"
all_data <- RenameIdents(all_data, cell_annotation_res.0.4)
all_data$major_annotation <- Idents(all_data) 


# annotation

pdf(
  "cycle3 - harmony - integrated - 16 - data - dimplot -  dim 13 - annotation.pdf",
  height = 6,
  width = 7
)
DimPlot(all_data, label = TRUE, reduction = "harmonyUMAP", group.by = "major_annotation") + NoLegend() + ggsci::scale_color_igv()
DimPlot(all_data, label = TRUE, reduction = "harmony_uamp", group.by = "major_annotation") + NoLegend() + ggsci::scale_color_igv()
dev.off()

# dotplot
pdf(
  "cycle3 - harmony - integrated - 16 - data - dotplot - RNA - slot - dim13 - annotation.pdf",
  height = 8,
  width = 18
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
DotPlot(all_data, assay = 'RNA', features = cell_marker, group.by = "major_annotation") + dot_theme + dot_scale
dev.off()

# save rds
saveRDS(all_data, '../all_data.harmony.re-analysis.rds')

## epic cell cluster  --------------------
sc_epi <- 
  CreateSeuratObject(
    counts = all_data@assays$RNA@counts[, all_data$major_annotation == "Epithelial cells"],
    meta.data = all_data[[]] %>% .[all_data$major_annotation == "Epithelial cells", ]
  )


# normalize
sc_epi <- 
  sc_epi %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = 'vst', nfeatures = nrow(.) * 0.9) %>%
  ScaleData(
    vars.to.regress = c("nCount_RNA")
  )


# PCA
sc_epi <- sc_epi %>% RunPCA()

# find optimal PCA dimensions
pdf("cycle3 - harmony - epi - qc pca.pdf", height = 6, width = 8)
ElbowPlot(sc_epi, ndims = 30)
dev.off()


set.seed(1234)
dim_heatmap_cells_epi <- 
  sc_epi[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 0.05)

pdf("cycle3 - harmony - epi -  dimheatmap 1-5.pdf", height = 5, width = 10)
DimHeatmap(sc_epi[, dim_heatmap_cells_epi$cell_id], dims = 1:5)
dev.off()

pdf("cycle3 - harmony - epi -  dimheatmap 6-10.pdf", height = 5, width = 10)
DimHeatmap(sc_epi[, dim_heatmap_cells_epi$cell_id], dims = 6:10)
dev.off()

pdf("cycle3 - harmony - epi -  dimheatmap 11-15.pdf", height = 5, width = 10)
DimHeatmap(sc_epi[, dim_heatmap_cells_epi$cell_id], dims = 11:15)
dev.off()



### dims6 - 使用这个 ---------------

# umap
sc_epi <- sc_epi %>% RunUMAP(dims = 1:6)
sc_epi <- sc_epi %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:6)
sc_epi <- sc_epi %>% RunUMAP(dims = 1:6, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")


pdf(
  "cycle3 - harmony - epi - integrated - 16 - data - dimplot - project.pdf",
  height = 6,
  width = 9
)
DimPlot(sc_epi, reduction = "umap", group.by = "project") + ggsci::scale_color_igv()
DimPlot(sc_epi, reduction = "harmonyUMAP", group.by = "project") + ggsci::scale_color_igv()
dev.off()



sc_epi <- FindNeighbors(object = sc_epi, reduction = "harmony", dims = 1:6)
sc_epi <- FindClusters(object = sc_epi, resolution = c(0.2, 0.4, 0.6, 0.8))


# dimplot
pdf(
  "cycle3 - harmony - epi - integrated - 16 - data - dimplot.pdf",
  height = 6,
  width = 7
)
DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()


DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP", group.by = "Tissue") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP", group.by = "TP53") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP", group.by = "Group") + NoLegend() + ggsci::scale_color_igv()

dev.off()


epi_marker2 <-
  c(`Ductal cells` = "CFTR", 
    # `Ductal cells` = "GPX2",  # Acinar cells/Ductal cells
    `Ductal cells` = "HES4",  # Pancreatic stellate cells/Ductal cells
    `Ductal cells` = "FSTL3", 
    
    `Exocrine glandular cells` = "AMY2A", 
    `Exocrine glandular cells` = "DPEP1",
    
    `Endocrine cells` = "CHGA", 
    
    `Alpha cells` = "GCG", 
    `Beta cells` = "INS", 
    `Delta cells` = "SST", 
    `Epsilon cells` = "GHRL", 
    `Gamma cells` = "PPY",
    
    Epithelial = "EPCAM", 
    Epithelial = "KRT18", 
    Epithelial = "KRT20"
  )

pdf(
  "cycle3 - harmony - epi - integrated - 16 - data - dotplot - RNA - slot.pdf",
  height = 10,
  width = 20
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
DotPlot(sc_epi, assay = 'RNA', features = epi_marker2, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
DotPlot(sc_epi, assay = 'RNA', features = epi_marker2, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
DotPlot(sc_epi, assay = 'RNA', features = epi_marker2, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
DotPlot(sc_epi, assay = 'RNA', features = epi_marker2, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale


dev.off()

# saveRDS(sc_epi, "sc_epi.dim6.rds")


# annotation
epi_anno_cycle3 <- c(
  '0'  = "Epithelial cells",
  '1'  = "Acinar cells", # "Exocrine glandular cells",
  '2'  = "Epithelial cells",
  '3'  = "Ductal cells",
  '4'  = "Acinar cells", # "Exocrine glandular cells",
  '5'  = "Endocrine cells"
)


Idents(sc_epi) <- "RNA_snn_res.0.2"
sc_epi <- RenameIdents(sc_epi, epi_anno_cycle3)
sc_epi$epi_annotation <- Idents(sc_epi)


# export 
pdf(
  "cycle3 - harmony - epi - integrated - 16 - data - dimplot - dim6 - annotation.pdf",
  height = 6,
  width = 7
)
DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP") + NoLegend() + ggsci::scale_color_igv()
dev.off()

saveRDS(sc_epi, "sc_epi.dim6.annotation.rds")


### dims10 ---------------

# umap
sc_epi <- sc_epi %>% RunUMAP(dims = 1:10)
sc_epi <- sc_epi %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:10)
sc_epi <- sc_epi %>% RunUMAP(dims = 1:10, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")


pdf(
  "cycle3 - harmony - epi - integrated - 16 - data - dimplot - project - dim10.pdf",
  height = 6,
  width = 9
)
DimPlot(sc_epi, reduction = "umap", group.by = "project") + ggsci::scale_color_igv()
DimPlot(sc_epi, reduction = "harmonyUMAP", group.by = "project") + ggsci::scale_color_igv()
dev.off()



sc_epi <- FindNeighbors(object = sc_epi, reduction = "harmony", dims = 1:10)
sc_epi <- FindClusters(object = sc_epi, resolution = c(0.2, 0.4, 0.6, 0.8))


# dimplot
pdf(
  "cycle3 - harmony - epi - integrated - 16 - data - dimplot - dim10.pdf",
  height = 6,
  width = 7
)
DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_epi, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()

dev.off()


pdf(
  "cycle3 - harmony - epi - integrated - 16 - data - dotplot - RNA - slot - dim10.pdf",
  height = 10,
  width = 20
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
DotPlot(sc_epi, assay = 'RNA', features = epi_marker2, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
DotPlot(sc_epi, assay = 'RNA', features = epi_marker2, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
DotPlot(sc_epi, assay = 'RNA', features = epi_marker2, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
DotPlot(sc_epi, assay = 'RNA', features = epi_marker2, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()

# saveRDS(sc_epi, "sc_epi.dim10.rds")


### inferCNV ---------------
# used to infer epethial cells
# immunes cell used as control

# 当前路径为pipeline文件夹
dir.create('../inferCNV') 
out_dir = "../inferCNV"

meta_data_infer <-
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  dplyr::filter(major_annotation %in% c("Epithelial cells", "T/NK cells", "B cells", "Fibroblasts")) %>%
  dplyr::select(cell_id, major_annotation, final_annotation) %>%
  map_dfc(as.character) %>%
  dplyr::mutate(annotation = ifelse(major_annotation == "Epithelial cells", final_annotation, major_annotation)) %>%
  dplyr::select(cell_id, annotation) %>%
  
  mutate(annotation = str_replace(annotation, "T\\/NK", "NKT"))

# unique(meta_data_infer$annotation)
# "Acinar cells"     "Ductal cells"     "Endocrine cells"  "Epithelial cells" "NKT cells"


# cell is too much
# less than 50000
# 取1/20细胞： 1179497 --> 47179.88
# 1662182 --> 83000
meta_data_infer_subset <-
  meta_data_infer %>%
  group_by(annotation) %>%
  slice_sample(prop = 0.05)

# export as file
meta_data_infer_subset %>%
  column_to_rownames("cell_id") %>%
  data.table::fwrite(
  file = file.path(out_dir, "meta_data_infer_subset.txt"),
  row.names = TRUE,
  col.names = FALSE,
  sep = "\t"
)

infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix = all_data@assays$RNA@counts[, meta_data_infer_subset$cell_id],
  annotations_file  = file.path(out_dir, "meta_data_infer_subset.txt"),
  gene_order_file   = file.path(out_dir, "gtf_infercnv_dropdup"),
  
  delim             = "\t",
  ref_group_names   = c("NKT cells", "Fibroblasts", "Acinar cells")
)


# all_data@assays$RNA@counts %>% rowMeans %>% quantile
# 0%          25%          50%          75%         100%
# 0.0008099304 0.0203149125 0.0807838652 0.2460118495 7.3282037744

all_data@assays$RNA@counts[, meta_data_infer_subset$cell_id] %>% rowMeans %>% { sum(. < 0.01) } # 87 <--- 
all_data@assays$RNA@counts[, meta_data_infer_subset$cell_id] %>% rowMeans %>% { sum(. < 0.005) } # 34



# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
infercnv_obj_default = infercnv::run(
  infercnv_obj,
  cutoff              =    0.01, 
  out_dir             =    out_dir,
  cluster_by_groups   =    TRUE, 
  plot_steps          =    FALSE,
  denoise             =    TRUE,
  HMM                 =    FALSE,
  no_prelim_plot      =    TRUE,
  png_res             =    300,
  no_plot = T # <--- figure margin too large ----
)

# C stack error
infercnv::plot_cnv(
  infercnv_obj_default,
  out_dir           =  out_dir,
  output_filename   =  'infercnv',
  # output_format     =  "pdf",
  png_res = 60,
  cluster_by_groups =  TRUE,
  x.range           =  "auto",
  x.center          =  1,
  title             =  "infercnv",
  color_safe_pal    =  FALSE
)


# download inferCNV object from linux server
# plot inferCNV locally
saveRDS(infercnv_obj_default, file.path(out_dir, 'infercnv_obj_default.rds'))

infercnv_obj_default <- readRDS("inferCNV/infercnv_obj_default.rds")

infercnv::plot_cnv(
  infercnv_obj_default,
  out_dir           =  'inferCNV', 
  output_filename   =  'infercnv',
  # output_format     =  "pdf",
  png_res = 300,
  cluster_by_groups =  TRUE,
  x.range           =  "auto",
  x.center          =  1,
  title             =  "infercnv",
  color_safe_pal    =  FALSE
)

infercnv_obj_default_median_filtering <-
  infercnv_obj_default %>%
  infercnv::apply_median_filtering()

infercnv_obj_default_median_filtering %>%
  infercnv::plot_cnv(
    out_dir           =  'inferCNV', 
    output_filename   =  'infercnv_median_filtering',
    # output_format     =  "pdf",
    png_res = 300,
    cluster_by_groups =  TRUE,
    x.range           =  "auto",
    x.center          =  1,
    title             =  "infercnv",
    color_safe_pal    =  FALSE
 )

# options()$bitmapType
# [1] "cairo"
# options(bitmapType="Xlib")
# options(bitmapType="cairo")


# tutorials
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix="D:/Program Files/R/R-4.1.1/library/infercnv/extdata/oligodendroglioma_expression_downsampled.counts.matrix.gz",
  annotations_file="D:/Program Files/R/R-4.1.1/library/infercnv/extdata/oligodendroglioma_annotations_downsampled.txt",
  delim="\t",
  gene_order_file="D:/Program Files/R/R-4.1.1/library/infercnv/extdata/gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt",
  ref_group_names= c("Microglia/Macrophage"), #c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")
)

# dir.create("test infercnv")
# out_dir = "test infercnv"

dir.create("test infercnv2")
out_dir = "test infercnv2"

infercnv_obj_default = infercnv::run(
  infercnv_obj,
  cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=out_dir,
  cluster_by_groups=TRUE, 
  plot_steps=FALSE,
  denoise=TRUE,
  HMM=FALSE,
  no_prelim_plot=TRUE,
  png_res=300
)



## T cell cluster  --------------------
sc_t <- 
  CreateSeuratObject(
    counts = all_data@assays$RNA@counts[, all_data$major_annotation == "T/NK cells"],
    meta.data = all_data[[]] %>% .[all_data$major_annotation == "T/NK cells", ]
  )


# normalize
sc_t <- 
  sc_t %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = 'vst', nfeatures = nrow(.) * 0.9) %>%
  ScaleData(
    vars.to.regress = c("nCount_RNA")
  )


# PCA
sc_t <- sc_t %>% RunPCA()

# find optimal PCA dimensions
pdf("cycle3 - harmony - t - qc pca.pdf", height = 6, width = 8)
ElbowPlot(sc_t, ndims = 30)
dev.off()


set.seed(1234)
dim_heatmap_cells_t <- 
  sc_t[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 0.2)

pdf("cycle3 - harmony - t -  dimheatmap 1-5.pdf", height = 5, width = 10)
DimHeatmap(sc_t[, dim_heatmap_cells_t$cell_id], dims = 1:5)
dev.off()

pdf("cycle3 - harmony - t -  dimheatmap 6-10.pdf", height = 5, width = 10)
DimHeatmap(sc_t[, dim_heatmap_cells_t$cell_id], dims = 6:10)
dev.off()

pdf("cycle3 - harmony - t -  dimheatmap 11-15.pdf", height = 5, width = 10)
DimHeatmap(sc_t[, dim_heatmap_cells_t$cell_id], dims = 11:15)
dev.off()



### dims6 ---------------

# umap
sc_t <- sc_t %>% RunUMAP(dims = 1:6)
sc_t <- sc_t %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:6)
sc_t <- sc_t %>% RunUMAP(dims = 1:6, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")


pdf(
  "cycle3 - harmony - t - integrated - 16 - data - dimplot - project.pdf",
  height = 6,
  width = 9
)
DimPlot(sc_t, reduction = "umap", group.by = "project") + ggsci::scale_color_igv()
DimPlot(sc_t, reduction = "harmonyUMAP", group.by = "project") + ggsci::scale_color_igv()
dev.off()



sc_t <- FindNeighbors(object = sc_t, reduction = "harmony", dims = 1:6)
sc_t <- FindClusters(object = sc_t, resolution = c(0.2, 0.4, 0.6, 0.8))


# dimplot
pdf(
  "cycle3 - harmony - t - integrated - 16 - data - dimplot.pdf",
  height = 6,
  width = 7
)
DimPlot(sc_t, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_t, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_t, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_t, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()
dev.off()




cell_marker_t <-
  c(
    "CD2",
    "CD3D",
    "CD3E",
    "CD4",
    "CD8A",
    
    # "PRF1",
    "GZMB",
    # "CD69",
    "FOXP3",
    # "TOP2A",
    "NKG7",
    # "GNLY", 
    # "KLRB1",
    # "KLRC1",
    # "KLRD1",
    # "NCR1",
    "FGFBP2"
  )

pdf(
  "cycle3 - harmony - t - integrated - 16 - data - dotplot - RNA - slot.pdf",
  height = 15,
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
DotPlot(sc_t, assay = 'RNA', features = cell_marker_t, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
DotPlot(sc_t, assay = 'RNA', features = cell_marker_t, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
DotPlot(sc_t, assay = 'RNA', features = cell_marker_t, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
DotPlot(sc_t, assay = 'RNA', features = cell_marker_t, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()

saveRDS(sc_t, "sc_t.dim6.rds")

### dims8 - 使用这个 ---------------

# umap
sc_t <- sc_t %>% RunUMAP(dims = 1:8)
sc_t <- sc_t %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:8)
sc_t <- sc_t %>% RunUMAP(dims = 1:8, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")


pdf(
  "cycle3 - harmony - t - integrated - 16 - data - dimplot - project - dim8.pdf",
  height = 6,
  width = 9
)
DimPlot(sc_t, reduction = "umap", group.by = "project") + ggsci::scale_color_igv()
DimPlot(sc_t, reduction = "harmonyUMAP", group.by = "project") + ggsci::scale_color_igv()
dev.off()



sc_t <- FindNeighbors(object = sc_t, reduction = "harmony", dims = 1:8)
sc_t <- FindClusters(object = sc_t, resolution = c(0.2, 0.4, 0.6, 0.8))


# dimplot
pdf(
  "cycle3 - harmony - t - integrated - 16 - data - dimplot - dim8.pdf",
  height = 6,
  width = 7
)
DimPlot(sc_t, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_t, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_t, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_t, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()
dev.off()




cell_marker_t <-
  c(
    "CD2",
    "CD3D",
    "CD3E",
    "CD4",
    "CD8A",
    
    # "PRF1",
    "GZMB",
    # "CD69",
    "FOXP3",
    # "TOP2A",
    "NKG7",
    "GNLY",
    "KLRB1",
    "KLRC1",
    "KLRD1",
    # "NCR1",
    "FGFBP2"
  )

pdf(
  "cycle3 - harmony - t - integrated - 16 - data - dotplot - RNA - slot - dim8.pdf",
  height = 15,
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
DotPlot(sc_t, assay = 'RNA', features = cell_marker_t, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
DotPlot(sc_t, assay = 'RNA', features = cell_marker_t, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
DotPlot(sc_t, assay = 'RNA', features = cell_marker_t, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
DotPlot(sc_t, assay = 'RNA', features = cell_marker_t, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()


t_marker_2 <- c(
  CD8 = "CD8A",
  NK = "NKG7",
  NK = "KLRB1",
  CD4 = "CD4",
  Cycling = "MKI67", 
  Cycling = "TOP2A", 
  B = "MS4A1", 
  B = "CD79A",
  ILC = "IL7R"
)

pdf(
  "cycle3 - harmony - t - integrated - 16 - data - dotplot - 2 - RNA - slot - dim8.pdf",
  height = 15,
  width = 25
)

dot_theme <-  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) 

DotPlot(sc_t, assay = 'RNA', features = t_marker_2, group.by = "RNA_snn_res.0.2") + dot_theme 
DotPlot(sc_t, assay = 'RNA', features = t_marker_2, group.by = "RNA_snn_res.0.4") + dot_theme 
DotPlot(sc_t, assay = 'RNA', features = t_marker_2, group.by = "RNA_snn_res.0.6") + dot_theme 
DotPlot(sc_t, assay = 'RNA', features = t_marker_2, group.by = "RNA_snn_res.0.8") + dot_theme 

dev.off()

# annotation
# t_anno_cycle3 <- c(
#   '0'  = "CD8Ts",
#   '1'  = "CD4Ts",
#   '2'  = "T/NK cells",
#   '3'  = "CD8Ts",
#   '4'  = "Tregs",
#   '5'  = "T/NK cells",
#   '6'  = "T/NK cells",
#   '7'  = "CD8Ts",
#   '8'  = "T/NK cells",
#   '9'  = "CD8Ts"
# )

t_anno_cycle3 <- c(
  '0'  = "CD8Ts",
  '1'  = "CD4Ts",
  '2'  = "NKT cells",
  '3'  = "CD8Ts",
  '4'  = "Tregs",
  '5'  = "NKT cells",
  '6'  = "NKT cells",
  '7'  = "CD8Ts",
  '8'  = "NKT cells",
  '9'  = "CD8Ts"
)

Idents(sc_t) <- "RNA_snn_res.0.6"
sc_t <- RenameIdents(sc_t, t_anno_cycle3)
sc_t$t_annotation <- Idents(sc_t)


# export 
pdf(
  "cycle3 - harmony - t - integrated - 16 - data - dimplot - dim8 - annotation.pdf",
  height = 6,
  width = 7
)
DimPlot(sc_t, label = TRUE, reduction = "harmonyUMAP") + NoLegend() + ggsci::scale_color_igv()
dev.off()

saveRDS(sc_t, "sc_t.dim8.annotation.rds")


## myeloid cluster ----------------------
sc_mye <- 
  CreateSeuratObject(
    counts = all_data@assays$RNA@counts[, all_data$major_annotation == "Myeloid cells"],
    meta.data = all_data[[]] %>% .[all_data$major_annotation == "Myeloid cells", ]
  )


# normalize
sc_mye <- 
  sc_mye %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = 'vst', nfeatures = nrow(.) * 0.9) %>%
  ScaleData(
    vars.to.regress = c("nCount_RNA")
  )


# PCA
sc_mye <- sc_mye %>% RunPCA()

# find optimal PCA dimensions
pdf("cycle3 - harmony - myeloid - qc pca.pdf", height = 6, width = 8)
ElbowPlot(sc_mye, ndims = 30)
dev.off()


set.seed(1234)
dim_heatmap_cells_myeloid <- 
  sc_mye[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 0.1)

pdf("cycle3 - harmony - myeloid -  dimheatmap 1-5.pdf", height = 5, width = 10)
DimHeatmap(sc_mye[, dim_heatmap_cells_myeloid$cell_id], dims = 1:5)
dev.off()

pdf("cycle3 - harmony - myeloid -  dimheatmap 6-10.pdf", height = 5, width = 10)
DimHeatmap(sc_mye[, dim_heatmap_cells_myeloid$cell_id], dims = 6:10)
dev.off()

pdf("cycle3 - harmony - myeloid -  dimheatmap 11-15.pdf", height = 5, width = 10)
DimHeatmap(sc_mye[, dim_heatmap_cells_myeloid$cell_id], dims = 11:15)
dev.off()



### dims5 ---------------

# umap
sc_mye <- sc_mye %>% RunUMAP(dims = 1:5)
sc_mye <- sc_mye %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:5)
sc_mye <- sc_mye %>% RunUMAP(dims = 1:5, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")
# sc_mye <- sc_mye %>%
#   RunUMAP(
#     dims = 1:13, 
#     min.dist = 0.01,
#     spread = 2,
#     reduction = "harmony", 
#     reduction.name = "harmony_uamp", 
#     reduction.key = "harmonyumap_"
#   )

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dimplot - project.pdf",
  height = 6,
  width = 9
)
DimPlot(sc_mye, reduction = "umap", group.by = "project") + ggsci::scale_color_igv()
DimPlot(sc_mye, reduction = "harmonyUMAP", group.by = "project") + ggsci::scale_color_igv()
dev.off()



sc_mye <- FindNeighbors(object = sc_mye, reduction = "harmony", dims = 1:5)
sc_mye <- FindClusters(object = sc_mye, resolution = c(0.2, 0.4, 0.6, 0.8))


# dimplot
pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dimplot.pdf",
  height = 6,
  width = 7
)
DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()
dev.off()


# dotplot
mye_marker <- c(
  pDC = "CLEC4C",
  pDC = "IRF7",
  pDC = "TCF4",
  pDC = "GZMB",
  
  cDC1 = "XCR1",
  cDC1 = "CLNK",
  cDC1 = "CLEC9A",
  
  cDC2 = "FCER1A",
  cDC2 = "HLA-DPB1",
  cDC2 = "HLA-DQB1",
  cDC2 = "CD1E",
  cDC2 = "CD1C",
  cDC2 = "CLEC10A",
  cDC2 = "HLA-DQA2",
  
  cDC3 = "CCL19",
  cDC3 = "LAMP3",
  cDC3 = "IDO1",
  cDC3 = "IDO2",
  cDC3 = "LAD1",
  cDC3 = "FSCN1",
  cDC3 = "CCR7",
  cDC3 = "LY75",
  cDC3 = "CCL22",
  cDC3 = "CD40",
  cDC3 = "BIRC3",
  cDC3 = "NFKB2",
  
  macrophage = "APOC1",
  macrophage = "APOE",
  macrophage = "HLA-DRB5",
  macrophage = "C1QA",
  macrophage = "C1QB",
  macrophage = "CD14",
  
  monocyte = "LYZ",
  monocyte = "HLA-DRB1",
  monocyte = "TIMP1",
  monocyte = "S100A11",
  monocyte = "CXCL8",
  monocyte = "IL1B",
  monocyte = "PTGS2",
  monocyte = "S100A9",
  monocyte = "S100A8",
  monocyte = "MMP19",
  
  # 'CD14+ monocytes' =	'CD14',
  # 'CD14+ monocytes'	= 'LYZ',
  
  'FCGR3A+ monocytes'	= 'FCGR3A',
  'FCGR3A+ monocytes'	= 'MS4A7',
  
  # Macrophages  = "MARCO",
  Macrophages  = "ITGAM", 
  Macrophages  = "ADGRE1",
  Macrophages  = "CD68",
  Macrophages  = "CD163",
  
  # M1 = "S100A8",
  M1 = "CXCL10",
  M1 = "SELENOP",
  
  M2 = "MARCO",
  M2 = "COL1A1",
  M2 = "MMP9"
)

mye_marker <- c(
  'CD14+ monocytes' =	'CD14',
  'CD14+ monocytes'	= 'LYZ',
  
  'FCGR3A+ monocytes'	= 'FCGR3A',
  'FCGR3A+ monocytes'	= 'MS4A7',
  
  'Conventional dendritic cells'	= 'FCER1A',
  'Conventional dendritic cells'	=  'CST3',
  
  # Macrophages  = "MARCO",
  Macrophages  = "ITGAM", 
  Macrophages  = "ADGRE1",
  Macrophages  = "CD68",
  Macrophages  = "CD163",
  
  M1 = "S100A8",
  M1 = "CXCL10",
  M1 = "SELENOP",
  
  M2 = "MARCO",
  M2 = "COL1A1",
  M2 = "MMP9"
  
)


# M1.S100A8,M2.CXCL10, M2.SELENOP, M2.MARCO, M2.COL1A1, M2.MMP9
# M1 CD80
# M2 CD163
# M3 ARG1

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dotplot - RNA - slot.pdf",
  height = 15,
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
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - 2 - dotplot - RNA - slot.pdf",
  height = 15,
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
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()

# saveRDS(sc_mye, "sc_mye.dim5.rds")

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - 2 - dotplot - RNA - slot.pdf",
  height = 15,
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
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()


### dims8 - 确定使用这个 -------

# umap
sc_mye <- sc_mye %>% RunUMAP(dims = 1:8)
sc_mye <- sc_mye %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:8)
sc_mye <- sc_mye %>% RunUMAP(dims = 1:8, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")
# sc_mye <- sc_mye %>%
#   RunUMAP(
#     dims = 1:13, 
#     min.dist = 0.01,
#     spread = 2,
#     reduction = "harmony", 
#     reduction.name = "harmony_uamp", 
#     reduction.key = "harmonyumap_"
#   )

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dimplot - project - dim8.pdf",
  height = 6,
  width = 9
)
DimPlot(sc_mye, reduction = "umap", group.by = "project") + ggsci::scale_color_igv()
DimPlot(sc_mye, reduction = "harmonyUMAP", group.by = "project") + ggsci::scale_color_igv()
dev.off()



sc_mye <- FindNeighbors(object = sc_mye, reduction = "harmony", dims = 1:8)
sc_mye <- FindClusters(object = sc_mye, resolution = c(0.2, 0.4, 0.6, 0.8))

sc_mye <- FindClusters(object = sc_mye, resolution = c(0.1))

# dimplot
pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dimplot - dim8.pdf",
  height = 6,
  width = 7
)
DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.1") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()
dev.off()


pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dotplot - RNA - slot - dim8.pdf",
  height = 10,
  width = 20
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
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()

# saveRDS(sc_mye, "sc_mye.dim8.rds")

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - 2 - dotplot - RNA - slot - dim8.pdf",
  height = 10,
  width = 20
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
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()



# findmarkers
Idents(sc_mye) <- "RNA_snn_res.0.2"
all_markers_dim8 <-
  FindAllMarkers(
    sc_mye,
    assay = "RNA",
    only.pos = TRUE
  )

all_markers_ls_dim8 <-
  sc_mye$RNA_snn_res.0.2 %>%
  unique %>%
  as.character() %>%
  set_names(., .) %>%
  map(function(x){
    FindConservedMarkers(
      sc_mye,
      ident.1 = x,
      grouping.var = "project",
      assay = "RNA",
      only.pos = TRUE
    )
  })
# FindConservedMarkers(
#   sc_mye,
#   ident.1 = "2",
#   grouping.var = "project",
#   assay = "RNA",
#   logfc.threshold = 0.1,
#   only.pos = TRUE
# )

# used for ACT 
all_markers_dim8 %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01) %>% 
  # dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)
  dplyr::slice_max(avg_log2FC, n = 10)  %>%
  group_by(cluster) %>%
  nest() %>% 
  deframe() %>%
  map(~ paste(.x$gene, collapse = ",")) %>%
  enframe() %>%
  unnest(value) %>%
  apply(1, function(x) str_glue('cluster_{x[1]}:{x[2]}')) %>%
  walk(cat, "\n")

# doheatmap
mye_gene_dim8 <-
  all_markers_dim8 %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01) %>% 
  # dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)
  dplyr::slice_max(avg_log2FC, n = 5) %>%
  dplyr::pull(gene)

mye_gene2_dim8 <- 
  all_markers_ls_dim8 %>%
  map(function(x){
    x %>%
      rownames_to_column("gene") %>%
      dplyr::filter(minimump_p_val < 0.05) %>% 
      rowwise() %>%
      mutate(avg_log2FC = mean(c_across(ends_with("avg_log2FC")))  )  %>%
      ungroup() %>%
      dplyr::slice_max(avg_log2FC, n = 5) %>%
      dplyr::pull('gene')
  }) 
mye_gene2_dim8 <- mye_gene2_dim8[as.character(seq(0, length(mye_gene2_dim8) - 1))]

set.seed(1234)
do_heatmap_cells_myeloid_dim8 <- 
  sc_mye[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(RNA_snn_res.0.2) %>%
  slice_sample(prop = 0.1)

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - doheatmap - RNA - slot - dim8.pdf",
  height = 10,
  width = 20
)
DoHeatmap(sc_mye[, do_heatmap_cells_myeloid_dim8$cell_id], features = mye_gene_dim8, group.by = "RNA_snn_res.0.2")
dev.off()

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - doheatmap - 2 - RNA - slot - dim8.pdf",
  height = 10,
  width = 20
)
DoHeatmap(sc_mye[, do_heatmap_cells_myeloid_dim8$cell_id], features = unlist(mye_gene2_dim8), group.by = "RNA_snn_res.0.2")
dev.off()

# 
# saveRDS(all_markers_dim8, 'all_markers_dim8.rds')
# saveRDS(all_markers_ls_dim8, 'all_markers_ls_dim8.rds')

### dims4 -------

# umap
# sc_mye <- sc_mye %>% RunUMAP(dims = 1:4)
sc_mye <- sc_mye %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:4)
sc_mye <- sc_mye %>% RunUMAP(dims = 1:4, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")
# sc_mye <- sc_mye %>%
#   RunUMAP(
#     dims = 1:13, 
#     min.dist = 0.01,
#     spread = 2,
#     reduction = "harmony", 
#     reduction.name = "harmony_uamp", 
#     reduction.key = "harmonyumap_"
#   )

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dimplot - project - dim4.pdf",
  height = 6,
  width = 9
)
# DimPlot(sc_mye, reduction = "umap", group.by = "project") + ggsci::scale_color_igv()
DimPlot(sc_mye, reduction = "harmonyUMAP", group.by = "project") + ggsci::scale_color_igv()
dev.off()



sc_mye <- FindNeighbors(object = sc_mye, reduction = "harmony", dims = 1:4)
sc_mye <- FindClusters(object = sc_mye, resolution = c(0.2))


# dimplot
pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dimplot - dim4.pdf",
  height = 6,
  width = 7
)
DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
# DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
# DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
# DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()
dev.off()


pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dotplot - RNA - slot - dim4.pdf",
  height = 10,
  width = 20
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
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
# DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
# DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
# DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()

# saveRDS(sc_mye, "sc_mye.dim8.rds")


mye_marker <- c(
  cDC1 = "CLEC9A",
  cDC1 = "CPNE3",
  cDC1 = "CPVL",
  cDC1 = "IDO1",
  cDC1 = "LGALS2",
  
  cDC2 = "AREG",
  cDC2 = "CD1C",
  cDC2 = "CLEC10A",
  cDC2 = "FCER1A",
  cDC2 = "S100B",

  mDC = "BIRC3",
  mDC = "CCL19",
  mDC = "CCR7",
  mDC = "MARCKSL1",
  mDC = "TXN",
  
  pDC = "BCL11A",
  pDC = "CLIC3",
  pDC = "GZMB",
  pDC = "JCHAIN",
  pDC = "PTGDS",
  
  M1 = "FCN1",
  M1 = "S100A12",
  M1 = "S100A8",
  M1 = "S100A9",
  M1 = "VCAN",
  
  M2 = "CCL8",
  M2 = "CXXL10",
  M2 = "CXCL11",
  M2 = "IFIT1",
  M2 = "IFIT2",
  
  M2 = "CCL3",
  M2 = "CCL3L1",
  M2 = "CCL4L2",
  M2 = "SELENOP",
  M2 = "SLC40A1",

  M2 = "MARCO",
  M2 = "MT1G",
  M2 = "MT1H",
  M2 = "MT1X",
  M2 = "PLIN2",
  
  M2 = "CALD1",
  M2 = "COL1A1",
  M2 = "COL1A2",
  M2 = "COL3A1",
  M2 = "SPARC",
  
  M2 = "APC5",
  M2 = "CKB",
  M2 = "CTSK",
  M2 = "MMP9",
  M2 = "SLC9B2"
)

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dotplot - 2 - RNA - slot - dim4.pdf",
  height = 10,
  width = 20
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
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
# DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
# DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
# DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()

saveRDS(sc_mye, "sc_mye.dim4.rds")



# findmarkers
Idents(sc_mye) <- "RNA_snn_res.0.2"
all_markers <-
  FindAllMarkers(
    sc_mye,
    assay = "RNA",
    only.pos = TRUE
  )

all_markers_ls <-
  sc_mye$RNA_snn_res.0.2 %>%
  unique %>%
  as.character() %>%
  set_names(., .) %>%
  map(function(x){
    FindConservedMarkers(
      sc_mye,
      ident.1 = x,
      grouping.var = "project",
      assay = "RNA",
      only.pos = TRUE
    )
  })



# doheatmap
mye_gene <-
  all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01) %>% 
  # dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)
  dplyr::slice_max(avg_log2FC, n = 5) %>%
  dplyr::pull(gene)

mye_gene2 <- 
  all_markers_ls %>%
  map(function(x){
    x %>%
      rownames_to_column("gene") %>%
      dplyr::filter(minimump_p_val < 0.05) %>% 
      rowwise() %>%
      mutate(avg_log2FC = sum(c_across(ends_with("c_across")))) %>%
      ungroup() %>%
      dplyr::slice_max(avg_log2FC, n = 5) %>%
      dplyr::pull('gene')
  }) 


set.seed(1234)
do_heatmap_cells_myeloid <- 
  sc_mye[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(RNA_snn_res.0.2) %>%
  slice_sample(prop = 0.1)

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - doheatmap - RNA - slot - dim4.pdf",
  height = 10,
  width = 20
)
DoHeatmap(sc_mye[, do_heatmap_cells_myeloid$cell_id], features = mye_gene, group.by = "RNA_snn_res.0.2")
dev.off()

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - doheatmap - 2 - RNA - slot - dim4.pdf",
  height = 10,
  width = 20
)
DoHeatmap(sc_mye[, do_heatmap_cells_myeloid$cell_id], features = unlist(mye_gene2), group.by = "RNA_snn_res.0.2")
dev.off()

# used for ACT 
all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.01) %>% 
  # dplyr::arrange(desc(avg_log2FC), .by_group = TRUE)
  dplyr::slice_max(avg_log2FC, n = 10)  %>%
  group_by(cluster) %>%
  nest() %>% 
  deframe() %>%
  map(~ paste(.x$gene, collapse = ",")) %>%
  enframe() %>%
  unnest(value) %>%
  apply(1, function(x) str_glue('cluster_{x[1]}:{x[2]}')) %>%
  walk(cat, "\n")


### dims14 -------

# umap
# sc_mye <- sc_mye %>% RunUMAP(dims = 1:14)
sc_mye <- sc_mye %>% harmony::RunHarmony(group.by.vars = "project", dims.use = 1:14)
sc_mye <- sc_mye %>% RunUMAP(dims = 1:14, reduction = "harmony", reduction.name = "harmonyUMAP", reduction.key = "harmonyUMAP_")
# sc_mye <- sc_mye %>%
#   RunUMAP(
#     dims = 1:13, 
#     min.dist = 0.01,
#     spread = 2,
#     reduction = "harmony", 
#     reduction.name = "harmony_uamp", 
#     reduction.key = "harmonyumap_"
#   )

pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dimplot - project - dim14.pdf",
  height = 6,
  width = 9
)
# DimPlot(sc_mye, reduction = "umap", group.by = "project") + ggsci::scale_color_igv()
DimPlot(sc_mye, reduction = "harmonyUMAP", group.by = "project") + ggsci::scale_color_igv()
dev.off()



sc_mye <- FindNeighbors(object = sc_mye, reduction = "harmony", dims = 1:14)
sc_mye <- FindClusters(object = sc_mye, resolution = c(0.2))


# dimplot
pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dimplot - dim14.pdf",
  height = 6,
  width = 7
)
DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.2") + NoLegend() + ggsci::scale_color_igv()
# DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.4") + NoLegend() + ggsci::scale_color_igv()
# DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.6") + NoLegend() + ggsci::scale_color_igv()
# DimPlot(sc_mye, label = TRUE, reduction = "harmonyUMAP", group.by = "RNA_snn_res.0.8") + NoLegend() + ggsci::scale_color_igv()
dev.off()


mye_marker <- c(
  DC = 'LAMP3',
  DC = 'IDO1',
  DC = 'IDO2',
  DC = 'CD1E',
  DC = 'CD1C',
  Macphage = "CD68",
  Macphage = "CD14",
  Macphage = "CD163",
  Macphage = "MARCO",
  Macphage = "MMP9"
)


pdf(
  "cycle3 - harmony - myeloid - integrated - 16 - data - dotplot - RNA - slot - dim14.pdf",
  height = 10,
  width = 20
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
DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.2") + dot_theme # + dot_scale
# DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.4") + dot_theme # + dot_scale
# DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.6") + dot_theme # + dot_scale
# DotPlot(sc_mye, assay = 'RNA', features = mye_marker, group.by = "RNA_snn_res.0.8") + dot_theme # + dot_scale

dev.off()

# saveRDS(sc_mye, "sc_mye.dim14.rds")


## cycle3 final annotation ---------------------

sc_t_anno <- sc_t[['t_annotation']] %>% rownames_to_column("cell_id") %>% mutate(t_annotation = as.character(t_annotation))
sc_epi_anno <- sc_epi[['epi_annotation']] %>% rownames_to_column("cell_id") %>% mutate(epi_annotation = as.character(epi_annotation))

all_data_anno <-
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  left_join(sc_t_anno, by = "cell_id") %>%
  left_join(sc_epi_anno, by = "cell_id") %>%
  mutate(
    major_annotation2 = major_annotation %>% as.character(),
    final_annotation = ifelse(major_annotation2 == "T/NK cells", t_annotation, major_annotation2),
    final_annotation = ifelse(major_annotation2 == "Epithelial cells", epi_annotation, final_annotation)
  ) %>%
  column_to_rownames("cell_id") 


all_data <- all_data %>% AddMetaData(all_data_anno['final_annotation'] )
Idents(all_data) <- "final_annotation"

## cycle3 export R objects ---------------
file.rename('../material_for_analysis', '../material_for_analysis_V1')
dir.create('../material_for_analysis')


all_data_meta <- all_data[[]]
saveRDS(all_data_meta, "../material_for_analysis/all_data_meta.rds")

avg_expr <- 
  AverageExpression(all_data, slot = "data")
saveRDS(avg_expr, "../material_for_analysis/avg_expr.rds")

norm_data <- all_data@assays$RNA@data
saveRDS(norm_data, "../material_for_analysis/norm_data_rna.rds")


# saveRDS(parse_data_path, "../material_for_analysis/parse_data_path.rds")
# saveRDS(image_wh_r_ls, "../material_for_analysis/image_wh_r_ls.rds")
# saveRDS(cells_boundary_ls, "../material_for_analysis/cells_boundary_ls.rds")
# saveRDS(cells_position_ls, "../material_for_analysis/cells_position_ls.rds")

file.copy("../material_for_analysis_V1/parse_data_path.rds", "../material_for_analysis/parse_data_path.rds")
file.copy("../material_for_analysis_V1/image_wh_r_ls.rds", "../material_for_analysis/image_wh_r_ls.rds")
file.copy("../material_for_analysis_V1/cells_boundary_ls.rds", "../material_for_analysis/cells_boundary_ls.rds")
file.copy("../material_for_analysis_V1/cells_position_ls.rds", "../material_for_analysis/cells_position_ls.rds")
file.copy("../material_for_analysis_V1/norm_data.rds", "../material_for_analysis/norm_data_integrated.rds")


saveRDS(all_data, "../material_for_analysis/all_data.rds")

all_data %>%
  Seurat::DietSeurat(
    .,
    scale.data = FALSE,
    assays = names(.@assays),
    dimreducs = Reductions(.)
  ) %>%
  saveRDS("../material_for_analysis/all_data_diet.rds")


# all cells
# cell_colors <- c(
#   'Acinar cells'     =  "#FB117E",
#   'Ductal cells'     =  "#D1207C",
#   'Endocrine cells'  =  "#BB1117",
# 
#   'B cells'          =  "#4180BD",
#   'CD4Ts'            =  "#21B746",
#   'CD8Ts'            =  "#006B14",
#   'Tregs'            =  "#5F8D66",
#   'NK/NKT'           =  "#9329AE",
# 
#   'cDCs'             =  "#954200",
# 
#   'Granulocytes'     =  "#F59F0C",
#   'Macrophage'       =  "#EF7111",
#   'Monocyte'         =  "#F9CC01",
# 
#   'Endothelial'      =  "#D56683",
#   'Fibroblasts'      =  "#008E8C"
# )

cell_colors <- c(
  # 'Acinar cells'     =  "#FB117E",
  'Exocrine glandular cells'     =  "#FB117E",
  'Ductal cells'     =  "#D1207C",
  'Endocrine cells'  =  "#BB1117",
  'Epithelial cells' =  "#F9CC01",
  
  
  'B cells'          =  "#4180BD",
  'CD4Ts'            =  "#21B746",
  'CD8Ts'            =  "#006B14",
  'Tregs'            =  "#5F8D66",
  'NKT cells'        =  "#9329AE", 
  
  'Granulocytes'     =  "#F59F0C",
  'Myeloid cells'    =  "#EF7111",

  'Endothelial'      =  "#D56683",
  'Fibroblasts'      =  "#008E8C"
)

## cycle3  export dimplot cycle3 - 2024.5.22 ------------------

dir_20240514 <- '2024.5.14 cycle3'
dir.create(dir_20240514)

pdf(file.path(dir_20240514, "integrated - 16 - data - umap - cycle3.pdf"))
DimPlot(
  object = all_data,
  reduction = "harmonyUMAP", 
  label = TRUE, 
  pt.size = 0.5,
  repel = TRUE
)+
  ggsci::scale_color_igv() +
  coord_fixed(ratio = 0.9)

DimPlot(
  object = all_data,
  reduction = "harmonyUMAP", 
  group.by = 'major_annotation',
  label = TRUE, 
  pt.size = 0.5,
  repel = TRUE
)+
  ggsci::scale_color_igv() +
  coord_fixed(ratio = 0.9)

DimPlot(object = all_data, group.by = 'Tissue', reduction = "harmonyUMAP", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'TP53', reduction = "harmonyUMAP", label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'Group', reduction = "harmonyUMAP", label = TRUE, pt.size = 0.5, repel = TRUE)
dev.off()

pdf(file.path(dir_20240514, "integrated - 16 - data - umap - split - dim30.pdf"), height = 25, width = 12)
DimPlot(
  object = all_data, 
  reduction = "harmonyUMAP",
  split.by = "project", 
  label = TRUE, 
  ncol = 3,
  pt.size = 0.5,
  repel = TRUE
) +
  ggsci::scale_color_igv()
DimPlot(
  object = all_data, 
  reduction = "harmonyUMAP",
  split.by = "final_annotation", 
  label = TRUE, 
  ncol = 3,
  pt.size = 0.5,
  repel = TRUE
) +
  ggsci::scale_color_igv()
DimPlot(object = all_data, group.by = 'Tissue', split.by = "project", reduction = "harmonyUMAP", ncol = 3, label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'TP53', split.by = "project", reduction = "harmonyUMAP", ncol = 3, label = TRUE, pt.size = 0.5, repel = TRUE)
DimPlot(object = all_data, group.by = 'Group', split.by = "project", reduction = "harmonyUMAP", ncol = 3, label = TRUE, pt.size = 0.5, repel = TRUE)
dev.off()



## cycle3 featureplot - marker gene ---------------
# 2024.5.22

p_feature_ls <-
  # c(
  #   "EPCAM", "ACTA2", "PECAM1", "FCN1",
  #   "CD3D",  "CD3E",  "NKG7", "CD79A", "MS4A1",
  #   "CPA3",
  #   "CFTR", "AMY2A", "CHGA"
  # ) %>% 
  c(
    "ACTA2", 
    "PECAM1", 
    "FCN1",
    "CD3D",  "CD3E", "CD8A", "CD4", "FOXP3",  "NKG7",
    "CD79A", "MS4A1",
    "CPA3",
    "EPCAM", "CFTR", "AMY2A", "CHGA"
  ) %>%
  map(function(x){
    FeaturePlot(all_data, features = x, reduction = "harmonyUMAP") + # 《---- reuducitoon ----------
      scale_colour_gradient2(
        low = "white",
        mid = "gray90",
        high = "blue",
        midpoint = 0
      )
  })

pdf(
  file.path(dir_20240514, "integrated - 16 - data - featureplot - dim30 - cycle3.pdf")
)
p_feature_ls %>% walk(print)
dev.off()


## cycle3 heatmapPlot - marker gene ---------------
# 2024.6.3

avg_expr <- all_data %>% AverageExpression %>% .[[1]]
cell_markers_for_heatmap <- c(
  "ACTA2", 
  "PECAM1", 
  "FCN1",
  "CD3D",  "CD3E", "CD8A", "CD4", "FOXP3",  "NKG7",
  "CD79A", "MS4A1",
  "CPA3",
  "EPCAM", "CFTR", "AMY2A", "CHGA"
)

cell_for_heatmap <- c(
  "Fibroblasts",
  "Endothelial cells",
  "Myeloid cells",
  "CD8Ts",
  "CD4Ts",
  "Tregs",
  "NKT cells",
  "B cells",
  "Granulocytes",
  "Epithelial cells",
  "Ductal cells",
  "Acinar cells",
  "Endocrine cells"
  
)



pdf(
  file.path(dir_20240514, "cycle3 - harmony - integrated - 16 - data - heatplot .pdf"),
  height = 5,
  width = 8
)
pheatmap::pheatmap(
  avg_expr[cell_markers_for_heatmap, cell_for_heatmap] %>% t,
  scale = "column",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA
)
dev.off()






# Spatial plot ---------------------------
##spatial plot  2024.3.24 - 不完全注释 ---------------------------
d
data_dir <- "../data"
anno_file <- "../data/Annotation.csv"
# ROI_anno_file <- "../data/ROI_Annotation.xlsx"
# panel_file <- "Xenium_hMulti_v1_metadata.xlsx"

# data path and project name
parse_data_path <-
  dir(data_dir, full.names = TRUE) %>%
  map(dir, full.names = TRUE) %>%
  unlist %>%

  parseProjName_xenium(pattern = "output[-]?", extract_postion = c(4, 3)) %>%
  mutate(.counts_path = file.path(.path, "cell_feature_matrix")) %>%
  identity()

# add id for parse_data_path
parse_data_path <-
  parse_data_path %>%
  mutate(.id = seq_along(.proj))


# image w/h ratio
# used to fix width and heigth ratio in In Situ Spatial plot.
image_wh_r_ls <-
  with(parse_data_path, split(.path, .proj)) %>%
  map(function(x){
    imag_wh <-
      RBioFormats::read.metadata(
        file.path(x, "morphology_mip.ome.tif")
      ) %>%
      .[[1]] %>%
      .$coreMetadata %>%
      .[c("sizeX", "sizeY")]

    image_wh_r <- imag_wh$sizeX/imag_wh$sizeY

    # try another method
    if(length(image_wh_r) <1){
      imag_wh <-
        RBioFormats::read.metadata(
          file.path(x, "morphology_mip.ome.tif")
        ) %>%
        .$coreMetadata %>%
        .[c("sizeX", "sizeY")]

      image_wh_r <- imag_wh$sizeX/imag_wh$sizeY
    }

    # length of iamge_wh_r should equal 1
    stopifnot(length(image_wh_r) == 1)

    image_wh_r
  })

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
# cells_position_ls <- cells_position_ls[proj_name]

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
# cells_boundary_ls <- cells_boundary_ls[proj_name]


# plot
all_data_meta_cycle2 <-
  all_data[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls_cycle2 <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]

    dat <-
      cells_boundary %>%
      right_join(all_data_meta_cycle2, by = "cell_id")

    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })


# export image
pdf(
  file.path(dir_20240324, "integrated - 16 - data - spatial plot.pdf"),
  height = 15,
  width  = 15
)
sp_p_ls_cycle2 %>% iwalk(~ {print(.x + ggtitle(.y))})
dev.off()


##spatial plot  2024.4.18 ---------------------------

cell_colors <- c(
  'Acinar cells'     =  "#FB117E",
  'Ductal cells'     =  "#D1207C",
  'Endocrine cells'  =  "#BB1117",
  
  'B cells'          =  "#4180BD",
  'CD4Ts'            =  "#21B746",
  'CD8Ts'            =  "#006B14",
  'Tregs'            =  "#5F8D66",
  'NK/NKT'           =  "#9329AE", 
  
  'cDCs'             =  "#954200",
  
  'Granulocytes'     =  "#F59F0C",
  'Macrophage'       =  "#EF7111",
  'Monocyte'         =  "#F9CC01",
  
  'Endothelial'      =  "#D56683",
  'Fibroblasts'      =  "#008E8C"
)

data_dir <- "../data"
anno_file <- "../data/Annotation.csv"
# ROI_anno_file <- "../data/ROI_Annotation.xlsx"
# panel_file <- "Xenium_hMulti_v1_metadata.xlsx"

# data path and project name
parse_data_path <- readRDS("../material_for_analysis/parse_data_path.rds")

# image w/h ratio
image_wh_r_ls <- readRDS("../material_for_analysis/image_wh_r_ls.rds")
  

# cell coordinate and cell bundary
cells_position_ls <- readRDS("../material_for_analysis/cells_position_ls.rds")
cells_boundary_ls <- readRDS("../material_for_analysis/cells_boundary_ls.rds")

# plot
all_data_meta_cycle2_v2 <-
  all_data2[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls_cycle2_v2 <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta_cycle2_v2, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = annotation_final, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      scale_fill_manual(values = cell_colors) +
      theme_void()
    p
  })


# export image
pdf(
  file.path(dir_20240409, "integrated - 16 - data - spatial plot.pdf"),
  height = 15,
  width  = 15
)
sp_p_ls_cycle2_v2 %>% iwalk(~ {print(.x + ggtitle(.y))})
dev.off()


# adipotal cells - marker不存在 ----------------

# ASC-1: SLC7A10
# PAT2: SLC36A2
# P2RX5: P2RX5

  DotPlot(
    all_data,
    features = c('SLC7A10', 'SLC36A2', 'P2RX5')
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.7, vjust = 0.5)
  ) +
  scale_colour_gradient2(
    low = "white",
    mid = "gray90",
    high = "blue",
    midpoint = 0
  )

  

# 测试umap -------------
# 先保存当前的数据
# umap 30, res0.2, annotation已是最新

all_data <- RunUMAP(
  object = all_data, 
  dims = 1:30,
  min.dist = 0.001,
  spread = 5,
  local.connectivity = 15,
  n.neighbors = 100 
)

# export
pdf(
  "integrated - 16 - data - umap - dim30 - adjust umap.pdf",
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "annotation") +
  ggsci::scale_color_igv()
DimPlot(all_data, group.by = "integrated_snn_res.0.2", label = TRUE) +
  ggsci::scale_color_igv()
dev.off()


pdf("integrated - 16 - data - umap - split - dim30  - adjust umap.pdf", height = 20, width = 20)
DimPlot(all_data, split.by = "project", group.by = "integrated_snn_res.0.2", label = TRUE, ncol = 4) +
  ggsci::scale_color_igv()
dev.off()



