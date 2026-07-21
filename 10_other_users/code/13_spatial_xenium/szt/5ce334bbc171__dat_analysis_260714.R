library(tidyverse)
library(Seurat)


out_dir_260714 <- "dat_analysis_cycle3/result_20260714"
dir.create(out_dir_260714)



# 1. fibroblast subcluster -----------
Fibroblast <- readRDS("dat_analysis_cycle3/subtype/Fibroblast/Fibroblast_subtype.rds")
DefaultAssay(Fibroblast) <- "sketch"

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
  'Dendritic cell' = "FCER1A",
  
  
  'Fat cell' = 'ADIPOQ',
  'Fat cell' = 'UCP1'
)


fibro_marker_ls <- c(
  Epithelial  = "EPCAM", 
  # Epithelial  = "KRT5",
  Epithelial  = "KRT8",
  Epithelial  = "S100A8",
  Epithelial  = "KRT18",
  
  "Mast cell" = "CPA3",
  "Mast cell" = "CTSG",
  "Mast cell" = "TPSAB1",
  
  'mCAF'='MMP11',
  'mCAF'='COL1A1',
  'mCAF'='POSTN',
  'mCAF'='COMP',
  'mCAF'='COL10A1',
  'mCAF'='COL11A1',
  'mCAF'='DCN',
  'mCAF'='CDH11',
  'mCAF'='LUM',
  'iCAF'='PLA2G2A',
  'iCAF'='CFD',
  'iCAF'='CD34',
  'iCAF'='CXCL14',
  'iCAF'='PI16',
  'iCAF'='PRG4',
  'iCAF'='APOD',
  'iCAF'='DPT',
  'iCAF'='CXCL12',
  'iCAF'='CXCL1',
  'iCAF'='IL6',
  'iCAF'='CCL2',
  'apCAF'='HLA-DRA',
  'apCAF'='HLA-DRB1',
  'apCAF'='CD74',
  'apCAF'='IGLC2',
  'apCAF'='IGKC',
  'apCAF'='CXCR4',
  'vCAF'='ACTA2',
  'vCAF'='MCAM',
  'vCAF'='RGS5',
  'vCAF'='MYH11',
  'vCAF'='ADIRF',
  'tCAF'='ENO1',
  'tCAF'='GAPDH',
  'tCAF'='NDRG1',
  # 'tCAF'='HSPH1',
  # 'tCAF'='VEGFA',
  'tCAF'='MME',
  'dCAF'='MKI67',
  'dCAF'='STMN1',
  'dCAF'='PTTG1',
  'dCAF'='H2AFZ',
  'dCAF'='CENPF',
  'infCAF'='IDO1',
  'infCAF'='CXCL9',
  'infCAF'='CXCL10',
  'infCAF'='CXCL11',
  'rCAF'='CCL19',
  'rCAF'='CCL21',
  'pericyte'='COL4A1',
  'pericyte'='COL4A2',
  'hsptcaf'='HSPH1',
  'hsptcaf'='VEGFA'
  
)


fibro_dims <- 10

## dotplot -------
pdf(file.path(out_dir_260714, str_glue("1 fibro - dotplot dims{fibro_dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.2", features = major_marker_ls) + RotatedAxis() + ggtitle('sketch_snn_res.0.2')
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.4", features = major_marker_ls) + RotatedAxis() + ggtitle('sketch_snn_res.0.4')
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.6", features = major_marker_ls) + RotatedAxis() + ggtitle('sketch_snn_res.0.6')
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.8", features = major_marker_ls) + RotatedAxis() + ggtitle('sketch_snn_res.0.8')
dev.off()

pdf(file.path(out_dir_260714, str_glue("1 fibro - dotplot dims{fibro_dims} fibro_marker.pdf")), width = 25, height = 6)
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.2", features = fibro_marker_ls) + RotatedAxis() + ggtitle('sketch_snn_res.0.2')
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.4", features = fibro_marker_ls) + RotatedAxis() + ggtitle('sketch_snn_res.0.4')
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.6", features = fibro_marker_ls) + RotatedAxis() + ggtitle('sketch_snn_res.0.6')
DotPlot(Fibroblast, group.by = "sketch_snn_res.0.8", features = fibro_marker_ls) + RotatedAxis() + ggtitle('sketch_snn_res.0.8')
dev.off()

## dimplot -------------------
pdf(file.path(out_dir_260714, str_glue("1 fibro - Dimplot dim{fibro_dims} cluster.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(Fibroblast, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.5", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Fibroblast, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

## find markers --------------
Idents(Fibroblast)<-'sketch_snn_res.0.6'
Fibroblast_markers <-
  FindAllMarkers(
    JoinLayers(Fibroblast),
    test.use = "wilcox",
    logfc.threshold = 1,
    min.pct = 0.1,
    only.pos = TRUE
  )

# export marker gene
Fibroblast_markers %>%
  dplyr::select(gene, everything()) %>%
  write.csv(
    file = file.path(out_dir_260714, "1 fibro - marker.csv"), 
    row.names = FALSE
  )

# vlnplot 
pdf(
  file.path(out_dir_260714, str_glue("1 fibro - vlnplot of ncount_rna - cluster.pdf")), 
  height = 6,
  width = 8
)
VlnPlot(Fibroblast, features = "nCount_RNA", group.by = "sketch_snn_res.0.2", pt.size = 0) + ggtitle('sketch_snn_res.0.2')
VlnPlot(Fibroblast, features = "nCount_RNA", group.by = "sketch_snn_res.0.4", pt.size = 0) + ggtitle('sketch_snn_res.0.4')
VlnPlot(Fibroblast, features = "nCount_RNA", group.by = "sketch_snn_res.0.6", pt.size = 0) + ggtitle('sketch_snn_res.0.6')
VlnPlot(Fibroblast, features = "nCount_RNA", group.by = "sketch_snn_res.0.8", pt.size = 0) + ggtitle('sketch_snn_res.0.8')
dev.off()

## annotation - v1 -------------
cell_cluster <- c(
  '0'  = "T-like fibrobalst",
  '1'  = "iCAF",
  '2'  = "mCAF",
  '3'  = "T-like fibrobalst",
  '4'  = "vCAF",
  '5'  = "low-quality cells",
  '6'  = "Epi-like fibroblast",
  '7'  = "T-like fibrobalst",
  '8'  = "low-quality cells",
  '9'  = "apCAF",
  '10' = "Mast-like fibroblast",
  '11' = "Epi-like fibroblast"
)

Idents(Fibroblast) <- "sketch_snn_res.0.6"

Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype2 <- Idents(Fibroblast)

# dimplot
pdf(file.path(out_dir_260714, str_glue("1 fibro - Dimplot dim{fibro_dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype2", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

# dotplot 
pdf(
  file.path(out_dir_260714, str_glue("1 fibro - dotplot dims{fibro_dims} major_major - annotation.pdf")), 
  width = 18,
  height = 4
)
  DotPlot(
    Fibroblast, 
    group.by = "subtype2",
    features = c(
      'T-like fibrobalst' = "CD3E",
      'T-like fibrobalst' = "CD3G",
      'Epi-like fibroblast' = 'EPCAM',
      'Mast-like fibroblast' = 'CPA3',
        
      'iCAF' = 'CD34',
      'iCAF' = 'DPT',
      'iCAF' = 'CXCL12',
      'iCAF' = 'IL6',
      
      'mCAF' = 'COL1A1',
      'mCAF' = 'POSTN',
      'mCAF' = 'LUM',
      
      'vCAF' = 'ACTA2',
      'vCAF' = 'MYH11',
      
      'apCAF' = "HLA-DRA",
      'apCAF' = "HLA-DRB1",
      'apCAF' = "CXCR4"
    )
  ) +
  RotatedAxis() 
dev.off()

pdf(
  file.path(out_dir_260714, str_glue("1 fibro - dotplot dims{fibro_dims} major_major - annotation - clean.pdf")), 
  width = 15,
  height = 4
)
DotPlot(
  Fibroblast[, Fibroblast$subtype2 %in% c('iCAF', 'mCAF', 'vCAF', 'apCAF')], 
  group.by = "subtype2",
  features = c(
    # 'T-like fibrobalst' = "CD3E",
    # 'T-like fibrobalst' = "CD3G",
    # 'Epi-like fibroblast' = 'EPCAM',
    # 'Mast-like fibroblast' = 'CPA3',
    
    'iCAF' = 'CD34',
    'iCAF' = 'DPT',
    'iCAF' = 'CXCL12',
    'iCAF' = 'IL6',
    
    'mCAF' = 'COL1A1',
    'mCAF' = 'POSTN',
    'mCAF' = 'LUM',
    
    'vCAF' = 'ACTA2',
    'vCAF' = 'MYH11',
    
    'apCAF' = "HLA-DRA",
    'apCAF' = "HLA-DRB1",
    'apCAF' = "CXCR4"
  )
) +
  RotatedAxis() 
dev.off()

## project data - v1 ------
Fibroblast <- 
  Fibroblast %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",#pca.full
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap' harmonyumap
    dims = 1:fibro_dims, 
    
    refdata = list(subtype_full2 = "subtype2")
  )

DefaultAssay(Fibroblast) <- 'RNA'

pdf(
  file.path(out_dir_260714, str_glue("Dimplot dim{fibro_dims} annotation - projected.pdf")), 
  height = 6, 
  width = 8
)
  DimPlot(
    Fibroblast, 
    group.by = "subtype_full2", 
    label = TRUE, 
    reduction = 'full.harmonyumap.2'
  ) + 
  NoLegend()
dev.off()

DefaultAssay(Fibroblast) <- 'sketch'
pdf(
  file.path(out_dir_260714, "dotplot celltype main_marker.pdf"),
  width = 15, 
  height = 6
)
  DotPlot(
    Fibroblast, 
    group.by = "subtype_full2",
    features = c(
      'T-like fibrobalst' = "CD3E",
      'T-like fibrobalst' = "CD3G",
      'Epi-like fibroblast' = 'EPCAM',
      'Mast-like fibroblast' = 'CPA3',
      
      'iCAF' = 'CD34',
      'iCAF' = 'DPT',
      'iCAF' = 'CXCL12',
      'iCAF' = 'IL6',
      
      'mCAF' = 'COL1A1',
      'mCAF' = 'POSTN',
      'mCAF' = 'LUM',
      
      'vCAF' = 'ACTA2',
      'vCAF' = 'MYH11',
      
      'apCAF' = "HLA-DRA",
      'apCAF' = "HLA-DRB1",
      'apCAF' = "CXCR4"
    ),
    dot.scale = 10
  ) + 
  RotatedAxis()
dev.off()


# saveRDS(Fibroblast, "dat_analysis_cycle3/subtype/Fibroblast/Fibroblast_subtype_260714.rds")


## annotation - v2 -------------
cell_cluster <- c(
  '0'  = "rCAF",
  '1'  = "iCAF",
  '2'  = "mCAF",
  '3'  = "rCAF",
  '4'  = "vCAF",
  '5'  = "low-quality cells",
  '6'  = "Epi-like fibroblast",
  '7'  = "T-like fibrobalst",
  '8'  = "low-quality cells",
  '9'  = "apCAF",
  '10' = "Mast-like fibroblast",
  '11' = "Epi-like fibroblast"
)

Idents(Fibroblast) <- "sketch_snn_res.0.6"

Fibroblast <- RenameIdents(Fibroblast, cell_cluster)
Fibroblast$subtype2 <- Idents(Fibroblast)

# dimplot
pdf(file.path(out_dir_260714, str_glue("1 fibro - Dimplot dim{fibro_dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Fibroblast, group.by = "subtype2", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
dev.off()

# dotplot 
pdf(
  file.path(out_dir_260714, str_glue("1 fibro - dotplot dims{fibro_dims} major_major - annotation.pdf")), 
  width = 18,
  height = 4
)
DotPlot(
  Fibroblast, 
  group.by = "subtype2",
  features = c(
    'T-like fibrobalst' = "CD3E",
    'T-like fibrobalst' = "CD3G",
    'Epi-like fibroblast' = 'EPCAM',
    'Mast-like fibroblast' = 'CPA3',
    
    'iCAF' = 'CD34',
    'iCAF' = 'DPT',
    'iCAF' = 'CXCL12',
    'iCAF' = 'IL6',
    
    'mCAF' = 'COL1A1',
    'mCAF' = 'POSTN',
    'mCAF' = 'LUM',
    
    'vCAF' = 'ACTA2',
    'vCAF' = 'MYH11',
    
    'apCAF' = "HLA-DRA",
    'apCAF' = "HLA-DRB1",
    'apCAF' = "CXCR4",
    
    'rCAF'='CCL19',
    'rCAF'='CCL21'
  )
) +
  RotatedAxis() 
dev.off()

pdf(
  file.path(out_dir_260714, str_glue("1 fibro - dotplot dims{fibro_dims} major_major - annotation - clean.pdf")), 
  width = 15,
  height = 4
)
DotPlot(
  Fibroblast[, Fibroblast$subtype2 %in% c('iCAF', 'mCAF', 'vCAF', 'apCAF', 'rCAF')], 
  group.by = "subtype2",
  features = c(
    # 'T-like fibrobalst' = "CD3E",
    # 'T-like fibrobalst' = "CD3G",
    # 'Epi-like fibroblast' = 'EPCAM',
    # 'Mast-like fibroblast' = 'CPA3',
    
    'iCAF' = 'CD34',
    'iCAF' = 'DPT',
    'iCAF' = 'CXCL12',
    'iCAF' = 'IL6',
    
    'mCAF' = 'COL1A1',
    'mCAF' = 'POSTN',
    'mCAF' = 'LUM',
    
    'vCAF' = 'ACTA2',
    'vCAF' = 'MYH11',
    
    'apCAF' = "HLA-DRA",
    'apCAF' = "HLA-DRB1",
    'apCAF' = "CXCR4",
    
    'rCAF'='CCL19',
    'rCAF'='CCL21'
  )
) +
  RotatedAxis() 
dev.off()

## project data - v2 ------
Fibroblast <- 
  Fibroblast %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",#pca.full
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap' harmonyumap
    dims = 1:fibro_dims, 
    
    refdata = list(subtype_full2 = "subtype2")
  )

DefaultAssay(Fibroblast) <- 'RNA'

pdf(
  file.path(out_dir_260714, str_glue("Dimplot dim{fibro_dims} annotation - projected.pdf")), 
  height = 6, 
  width = 8
)
DimPlot(
  Fibroblast, 
  group.by = "subtype_full2", 
  label = TRUE, 
  reduction = 'full.harmonyumap.2'
) + 
  NoLegend()
dev.off()

DefaultAssay(Fibroblast) <- 'sketch'
pdf(
  file.path(out_dir_260714, "dotplot celltype main_marker.pdf"),
  width = 15, 
  height = 6
)
DotPlot(
  Fibroblast, 
  group.by = "subtype_full2",
  features = c(
    'T-like fibrobalst' = "CD3E",
    'T-like fibrobalst' = "CD3G",
    'Epi-like fibroblast' = 'EPCAM',
    'Mast-like fibroblast' = 'CPA3',
    
    'iCAF' = 'CD34',
    'iCAF' = 'DPT',
    'iCAF' = 'CXCL12',
    'iCAF' = 'IL6',
    
    'mCAF' = 'COL1A1',
    'mCAF' = 'POSTN',
    'mCAF' = 'LUM',
    
    'vCAF' = 'ACTA2',
    'vCAF' = 'MYH11',
    
    'apCAF' = "HLA-DRA",
    'apCAF' = "HLA-DRB1",
    'apCAF' = "CXCR4",
    
    'rCAF'='CCL19',
    'rCAF'='CCL21'
  ),
  dot.scale = 10
) + 
  RotatedAxis()
dev.off()


# saveRDS(Fibroblast, "dat_analysis_cycle3/subtype/Fibroblast/Fibroblast_subtype_260717.rds")

## clean celltype, re-run umap ------------
fibro_others_cells <- Fibroblast$subtype2 %>% unique %>%  str_subset("-like|low-qua")

Fibroblast2 <- Fibroblast[,
  (!Fibroblast$subtype2 %in% fibro_others_cells) &
  (!Fibroblast$subtype_full2 %in% fibro_others_cells)
]

DefaultAssay(Fibroblast2) <- 'sketch'
Fibroblast2 <- RunUMAP(
  Fibroblast2, 
  dims = 1:fibro_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)

pdf(
  file.path(out_dir_260714, str_glue("Dimplot dim{fibro_dims} annotation - clean - rerun umap.pdf")), 
  height = 6, 
  width = 8
)
DimPlot(
  Fibroblast2, 
  group.by = "subtype2", 
  label = TRUE, 
  reduction = 'harmonyumap'
) + 
  NoLegend()
dev.off()


# project 
Fibroblast2 <- 
  Fibroblast2 %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",#pca.full
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap' harmonyumap
    dims = 1:fibro_dims, 
    
    refdata = list(subtype_full2 = "subtype2")
  )

DefaultAssay(Fibroblast) <- 'RNA'

pdf(
  file.path(out_dir_260714, str_glue("Dimplot dim{fibro_dims} annotation - projected - clean - rerun umap.pdf")), 
  height = 6, 
  width = 8
)
DimPlot(
  Fibroblast2, 
  group.by = "subtype_full2", 
  label = TRUE, 
  reduction = 'full.harmonyumap.2'
) + 
  NoLegend()
dev.off()

# saveRDS(Fibroblast2, "dat_analysis_cycle3/subtype/Fibroblast/Fibroblast_subtype_cleacn_260717.rds")
