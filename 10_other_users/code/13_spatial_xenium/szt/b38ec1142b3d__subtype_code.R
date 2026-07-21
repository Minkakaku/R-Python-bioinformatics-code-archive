library(Seurat)
library(tidyverse)
library(BPCells)
#All_T subtype-----------------------
all_data<-readRDS('/home/rstudio/R_OUTPUT/all_data.rds')
All_T<-subset(all_data,annotation=='T cells')
All_T<-JoinLayers(All_T)
### PCA --------
DefaultAssay(All_T) <- "RNA"
All_T <- NormalizeData(All_T, scale.factor = median(All_T$nCount_RNA))
All_T <- FindVariableFeatures(All_T)

All_T <- ScaleData(All_T, vars.to.regress = 'nCount_RNA')

All_T <- RunPCA(All_T)

out_dir<-'T_subtype'
dir.create(out_dir)
### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(All_T, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  All_T[[]] %>%
  rownames_to_column("cell_id2") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(All_T)) %>%
  dplyr::pull("cell_id2")


pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(All_T[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(All_T[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(All_T[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(All_T[, suset_cells], dims = 16:20)
dev.off()
### dimension reduction --------
dims <- 22


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


All_T <- RunUMAP(All_T, dims = 1:dims,return.model = TRUE)



pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


All_T <- All_T %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

All_T <- RunUMAP(
  All_T, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
All_T  <-
  FindNeighbors(
    All_T ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


All_T <-
  FindClusters(
    All_T,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
# DimPlot(All_T, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(All_T, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(All_T, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(All_T, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(All_T, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(All_T, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

### dotplot -------
major_marker_ls <- c(
  'CD8'='CD8A',
  'CD8'='CD8B',
  'CD8'='GZMB',
  'CD8'='GZMA',
  'CD8'='GZMK',
  'CD8'='GZMH',
  'CD8'='CCL5',
  'Naive'='CCR6',
  'Naive'='LEF1',
  'Naive'='TCF7',
  'Naive'='SELL',
  'Naive'='CCR7',
  'Naive'='YBX3',
  'Naive'='TGFB1',
  'Naive'='CCl5',
  'Naive'='CXCR3',
  'Memory'='GPR183',
  'Memory'='IL7R',
  'Memory'='FOS',
  'Memory'='JUN',
  'Memory'='KLRB1',
  'TRM'='CD89',
  'TRM'='CXCR8',
  # 'TEM'='IFNG',
  'TEM'='OCL4',
  'TEM'='OCL5',
  'TEM'='JAK3',
  # 'TCM'='IDO1',
  # 'TCM'='SKIV2L',
  # 'TCM'='TRAF2',
  # 'TCM'='SLC2A1',
  # 'TCM'='CD40LG',
  'TCM'='DUSP2',
  'TCM'='LYAR',
  'TCM'='ZNF683',
  'TCM'='IL32',
  'Effector'='FGFBP2',
  'Effector'='TYROBP',
  'Effector'='KLRC2',
  'Effector'='CX3CR1',
  'Effector'='FCGR3A',
  'Effector'='GNLY',
  'cytotoxic'='XCL1',
  'cytotoxic'='IFNG',
  'cytotoxic'='CCL4',
  # 'cytotoxic'='CCL5',
  'Exhausted'='CXCL13',
  'Exhausted'='LAYN',
  'Exhausted'='LAG3',
  'Exhausted'='HAVCR2',
  'Exhausted'='TNFRSF9',
  'Exhausted'='CTLA4',
  'Exhausted'='TIGIT',
  'Exhausted'='TOX',
  'Exhausted'='PDCD1',
  'CD4'='FOXP3',
  'CD4'='CD4',
  'CD4'='TRDC',
  'CD4'='IL2RA',
  'CD4'='CD25',
  'CD4'='STAT1',
  'CD4'='IL2',
  'Tfh'='CD200',
  'Tfh'='CD83',
  'Tfh'='CXCL3',
  'Tfh'='ICOS',
  'Tfh'='BATF',
  'Tfh'='PDCD1',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",


  'proliferation'='MKI67',
  'proliferation'='H4C3',
  'proliferation'='TUBB',
  'proliferation'='STMN1',
  'proliferation'='TOP1',
  'proliferation'='ANAPC11',
  # 'proliferation'='MCM7',
  # 'proliferation'='DNMT1',
  'proliferation'='TOP2A'
)

major_marker_ls_CD4 <- c(
  # "B cells" = "CD19",
  # "B cells" = "CD79A",
  # "B cells" = "MS4A1",
  'T' ='CD3E',
  'CD4'='CD4',
  # 'CD4'='TRDC',
  # 'CD4'='IL2RA',
  # 'CD4'='CD25',
  # 'CD4'='STAT1',
  # 'CD4'='IL2',
  'CD4_Treg'='TNFRSF4',
  'CD4_Treg'='TNFRSF9',
  'CD4_Treg'='FOXP3',
  'CD4_Treg'='IL2RA',
  'CD4_Treg'='LAYN',
  'CD4_Treg'='CTLA4',
  'CD4_Treg'='TIGIT',
  'CD4_naive'='LEF1',
  'CD4_Tem'='GZMK',
  'CD4_Tem'='GZMA',
  'CD4_Tem'='NKG7',
  'CD4_Tem'='CCL4',
  'CD4_Tem'='CCL5',
  'CD4_Tem'='CAPG',
  'CD4_Th1'='ANXA1',
  'CD4_Th1'='CXCL13',
  'CD4_Th1'='BHLHE40',
  'CD4_Th1'='IFNG',
  'CD4_Tfh'='TOX2',
  'CD4_Tfh'='PDCD1',
  'CD4_Tfh'='CXCR5',
  'CD4_Th17'='IL17A',
  'CD4_Th17'='IL17F',
  'CD4_Th17'='IL23R',
  'CD4_Th17'='RORC',
  'CD4_Th17'='CCR6',
  'CD4_Tcm'='TCF7',
  'CD4_Tcm'='CCR7',
  'CD4_Tcm'='S1PR1',
  'CD4_Tcm'='TXK',
  'CD4_Tcm'='ICAM2',
  'CD4_Tcm'='CDC25B',
  'CD4_Tcm'='PTGER2'
)


major_marker_ls_CD8 <- c(  
    'CD8'='CD8A',
    'CD8'='CD8B',
    'CD8_Tex'='HAVCR2',
    'CD8_Tex'='CTLA4',
    'CD8_Tex'='LAYN',
    'CD8_Tex'='LAG3',
    'CD8_Tex'='TIGIT',
    'CD8_Tex'='ENPD1',
    'CD8_Tem'='GZMK',
    'CD8_Tem'='GZMA',
    'CD8_Tem'='IFNG',
    'CD8_Tem'='KLRG1',
    'CD8_Tem'='GZMB',
    'CD8_Tem'='PRF1',
    'CD8_Tem'='NKG7',
    'CD8_Tnaive'='SELL',
    'CD8_Tnaive'='LEF1',
    'CD8_Tnaive'='TCF7',
    'CD8_Tnaive'='CCR7',
    'CD8_MAIT'='NCR3',
    'CD8_MAIT'='SLC4A10',
    'CD8_MAIT'='KLRB1',
    'CD8_Trm'='MYADM',
    'CD8_Trm'='CD69',
    'CD8_Trm'='ITGAE',
    'CD8_Trm'='CXCR6',
    'CD8_Trm'='ITGA1',
    'CD8_Trm'='ZNF683',
    'CD8_Tcm'='IL7R',
    'CD8_Tcm'='CD28',
    'CD8_Tcm'='GPR183'
    )

major_marker_ls_CD4_given <- c(
  "T_cells" = "CD3D",
  "T_cells" = "CD3E",
  "T_cells" = "CD3G",
  "CD4_T" = "CD4",
  "CD8_T" = "CD8A",
  "CD8_T" = "CD8B",
  "CD4_Naive" = "CCR7",
  "CD4_Naive" = "SELL",
  "CD4_Naive" = "TCF7",
  "CD4_Naive" = "LEF1",
  "CD4_Naive" = "IL7R",
  "CD4_Naive" = "CD28",
  "CD4_Naive" = "CD27",
  "CD4_Naive" = "TNFRSF25",
  "Th1" = "TBX21",
  "Th1" = "IFNG",
  "Th1" = "IL12RB2",
  "Th1" = "CXCR3",
  "Th1" = "STAT1",
  "Th1" = "STAT4",
  "Th1" = "TNF",
  "Th1" = "IL2",
  "Th2" = "GATA3",
  "Th2" = "IL4",
  "Th2" = "IL5",
  "Th2" = "IL13",
  "Th2" = "CCR4",
  "Th2" = "STAT6",
  "Th2" = "PTGDR2",
  "Th2" = "IL1RL1",
  "Th17" = "RORC",
  "Th17" = "IL17A",
  "Th17" = "IL17F",
  "Th17" = "IL23R",
  "Th17" = "CCR6",
  "Th17" = "STAT3",
  "Th17" = "IL22",
  "Th17" = "IL21",
  "Tfh" = "CXCR5",
  "Tfh" = "BCL6",
  "Tfh" = "PDCD1",
  "Tfh" = "ICOS",
  "Tfh" = "SH2D1A",
  # "Tfh" = "IL21",
  "Tfh" = "CD40LG",
  "Tfh" = "BTLA",
  "Treg" = "FOXP3",
  "Treg" = "IL2RA",
  "Treg" = "CTLA4",
  "Treg" = "TNFRSF18",
  "Treg" = "IKZF2",
  "Treg" = "ENTPD1",
  "Treg" = "TIGIT",
  "Treg" = "IL10",
  "proliferating T" = 'TOP2A',
  "proliferating T" = 'MKI67'
)



major_marker_ls_CD8_given<-c(
  "T_CD8_naive" = "CCR7",
  "T_CD8_naive" = "SELL",
  "T_CD8_naive" = "TCF7",
  "T_CD8_naive" = "LEF1",
  "T_CD8_naive" = "IL7R",
  "T_CD8_naive" = "CD28",
  "T_CD8_naive" = "CD27",
  "T_CD8_central_memory" = "CD44",
  "T_CD8_central_memory" = "CXCR3",
  "T_CD8_central_memory" = "BCL2",
  "T_CD8_effector_memory" = "GZMK",
  "T_CD8_effector_memory" = "GZMA",
  "T_CD8_effector_memory" = "KLRG1",
  "T_CD8_effector_memory" = "CX3CR1",
  "T_CD8_effector_memory" = "EOMES",
  "T_CD8_effector" = "GZMB",
  "T_CD8_effector" = "GNLY",
  "T_CD8_effector" = "PRF1",
  "T_CD8_effector" = "NKG7",
  "T_CD8_effector" = "TBX21",
  "T_CD8_exhausted" = "PDCD1",
  "T_CD8_exhausted" = "HAVCR2",
  "T_CD8_exhausted" = "LAG3",
  "T_CD8_exhausted" = "CTLA4",
  "T_CD8_exhausted" = "TOX",
  "T_CD8_exhausted" = "TIGIT",
  "T_CD8_exhausted" = "ENTPD1"
)

celltype_marker<-c(
  "B cell"= "CD19", 
  "B cell"= "CD79A", 
  "B cell"= "CD79B", 
  "B cell"= "MS4A1",
  "Plasma" = "TNFRSF17", 
  "Plasma" = "MZB1",
  # "Plasma" = "XBP1",
  "T cell" = "CD3E", 
  "T cell" = "CD3G", 
  "T cell" = "CD4",
  "T cell" = "FOXP3", 
  "T cell" = "CD8A", 
  "T cell" = "GZMB",
  "Endothelial" = "PECAM1",
  "Endothelial" = "CLDN5", 
  "Epithelial" = "EPCAM",
  "Epithelial" = "KRT5",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "PLN",
  "Smooth muscle cells" = "MYLK",
  "Monocytes" = "ITGAM",
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  "Monocytes" = "CD68",
  "Macrophage" = "CSF1R",
  "Macrophage" = "MARCO", 
  "DC"= "FCER1A", 
  "DC"= "CD1C", 
  "DC"= "CD83", 
  "DC"= "LILRA4", 
  "DC"= "CLEC9A", 
  "DC"= "IDO1",
  # 'pDC' = 'CLEC4C',
  # 'pDC' = "LILRA4",
  # 'pDC' = 'IL3RA',
  # 'pDC' = 'IRF8',
  # 'mregDC' = 'CCR7',
  # 'mregDC' ='CCL19',
  # 'mregDC' ='LAMP3',
  # 'mregDC' ='CD83',
  "NK" = "NCAM1",
  "Mast cell" = "CTSG", 
  "Mast cell" = "MS4A2", 
  "Mast cell" = "KIT", 
  "Mast cell" = 'HDC', 
  'Neutrophil'= "CSF3R", 
  'Neutrophil'= "FCGR3B", # "S100A8", "S100A9" 
  "Melanocy" = "MITF",
  "Melanocy" = "MLANA"
)







# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
# DotPlot(All_T, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
# DotPlot(All_T, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.8", features = celltype_marker)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.8", features = major_marker_ls_CD4)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.8", features = major_marker_ls_CD8)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.8", features = major_marker_ls_CD4_given)+RotatedAxis()
DotPlot(All_T, group.by = "RNA_snn_res.0.8", features = major_marker_ls_CD8_given)+RotatedAxis()
dev.off()
# find markers -------
Idents(All_T) <- "RNA_snn_res.0.6"
all_markers <-
  FindAllMarkers(
    JoinLayers(All_T),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  # arrange(desc(pct.1)) %>%
  slice_head(n = 30) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

### annotation ----

cell_cluster <- c(
  '0'  = "CD4_Tfh",
  '1'  = "Naive_T",
  '2'  = "Naive_T",
  '3'  = "CD8_Tem",
  '4'  = "Treg",
  '5'  = "CD8_naive",
  '6'  = "Naive_T",
  '7'  = "Tfh",
  '8'  = "Naive_T",
  '9'  = "Tfh",
  '10'  = "proliferating T",
  '11'  = "CD8_MAIT",
  '12'  = "CD4_naive",
  '13'  = "proliferating T",
  '14'  = "Naive_T",
  '15'  = "CD4_naive",
  '16'  = "Th2",
  '17'  = "proliferating T",
  '18'  = "proliferating T"
)


Idents(All_T) <- "RNA_snn_res.0.8"

All_T <- RenameIdents(All_T, cell_cluster)
All_T$subtype <- Idents(All_T)


pdf(file.path(out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf(file.path(out_dir,"dotplot celltype main_marker.pdf"), width = 20, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(All_T, group.by = "subtype", features = major_marker_ls_CD4_given, dot.scale = 10) + theme_90angle
DotPlot(All_T, group.by = "subtype", features = major_marker_ls_CD8, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()


### project data ------
All_T <- 
  All_T %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'harmonyumap
    dims = 1:dims, 
    
    refdata = list(subtype_full = "subtype")
  )


DefaultAssay(All_T) <- 'RNA'

pdf(file.path("/home/rstudio/dat_analysis_cycle3/subtype/All_T",str_glue("Dimplot dim{dims} annotation - projected.pdf")), height = 6, width = 8)
DimPlot(All_T, group.by = "subtype_full", label = TRUE, reduction = 'full.harmonyumap.1') + NoLegend()
dev.off()

# saveRDS(All_T,"T_subtype/All_T_subtype.rds")
# all_data <- readRDS("CD8_subtype.rds")

#spatial plot------------
All_T_meta <-
  All_T[[]] %>%
  rownames_to_column("cell_id")
all_data_meta <-
  all_data2[[]] %>%
  rownames_to_column("cell_id")


sp_p_ls <-
  cells_boundary_ls2 %>% 
  imap(function(cells_boundary, project){  #browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat2 <-
      cells_boundary %>%
      right_join(All_T_meta, by = "cell_id")
    
    p <- dat[which(dat$annotation != 'T cells'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray80")  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal()+
      theme_void() +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat2, aes(x = vertex_x, y = vertex_y, fill = subtype, group = cell_id)) +
      scale_fill_manual(values = plan)+
      theme_void()
    p
  })


plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", 
          "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", 
          "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879", 
          "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", 
          "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", "#ffc300", 
          "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")


sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("T_subtype/spatial_plot", str_glue("T subtype spatial plot {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(All_T_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })

all_cells <- unique(All_T_meta$subtype)
cell_cols <- structure(
  ggsci::pal_igv()(length(all_cells)),
  names = all_cells
)

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){  #browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray90")  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      theme_void()
    
    p_ls <-
      all_cells %>%
      set_names(., .) %>%
      map(function(cell){ #browser()
        dat_sub <- dat %>% dplyr::filter(annotation3 == .env[['cell']])
        p +
          geom_polygon(
            data = dat_sub,
            aes(group = cell_id), 
            fill = cell_cols[cell]
          )  +
          ggtitle(cell)
      })
    p_ls
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path("T_subtype/spatial_plot", str_glue("spatial plot {grp} celltype split.pdf")),
      height = 12,
      width  = 8
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })
# cell props barplot  --------

All_T_meta <-
  All_T[[]] %>%
  rownames_to_column("cell_id")

All_T_meta$project <- factor(All_T_meta$project,levels=c("GX_1","GYF_1","HXF_1","LS_1","MXH_1","GX_2","GYF_2","HXF_2","LS_2","MXH_2","GX_3","GYF_3","HXF_3","LS_3","MXH_3"))
p_col_proj <-
  All_T_meta %>%
  group_by(project, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


p_box_group  <-
  all_data_meta %>%
  group_by(group, project, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation, y = freq, fill = group)) +
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
  filename = file.path("T_subtype", "cell_freq_group.pdf"),
  height = 6,
  width  = 8
)

ggsave(
  p_box_group,
  filename = file.path(dr_dir, "cell_freq_group.pdf"),
  height = 6,
  width  = 10
)

##映射------
library(tibble)
library(dplyr)
meta.data.t <-
  All_T[['subtype']] %>%
  rownames_to_column('cells')%>%
  right_join(
    all_data[['annotation']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation)=='T cells'),as.character(subtype), as.character(annotation)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation2'])
pdf(file.path("mapping_annotation T umap.pdf"), height = 6, width = 12)
DimPlot(all_data, group.by = "annotation3", label = TRUE, reduction = 'harmonyumap',pt.size = 2,repel = TRUE) 
#DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()





#spatial plot------------
Fibroblast_meta <-
  Fibroblast[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){  #browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat2 <-
      cells_boundary %>%
      right_join(Fibroblast_meta, by = "cell_id")
    
    p <- dat[which(dat$annotation2 != 'Fibroblast'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray60")  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      theme_void() +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat2, aes(x = vertex_x, y = vertex_y, fill = subtype, group = cell_id)) +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("Fibroblast_subtype/spatial_plot", str_glue("Fibroblast subtype spatial plot {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })
# 
# sp_p_ls <-
#   cells_boundary_ls %>%
#   imap(function(cells_boundary, project){ # browser()
#     image_wh_r <- image_wh_r_ls[[project]]
#     
#     dat <-
#       cells_boundary %>%
#       right_join(Macrophage_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = vertex_x, y = vertex_y)) +
#       geom_polygon(aes(fill = subtype, group = cell_id))  +
#       scale_y_reverse() +
#       coord_fixed(ratio = 1/image_wh_r) +
#       ggsci::scale_fill_igv() +
#       theme_void()
#     p
#   })

all_cells <- unique(Fibroblast_meta$subtype)
cell_cols <- structure(
  ggsci::pal_igv()(length(all_cells)),
  names = all_cells
)

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){  #browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray60")  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      theme_void()
    
    p_ls <-
      all_cells %>%
      set_names(., .) %>%
      map(function(cell){ #browser()
        dat_sub <- dat %>% dplyr::filter(annotation3 == .env[['cell']])
        p +
          geom_polygon(
            data = dat_sub,
            aes(group = cell_id), 
            fill = cell_cols[cell]
          )  +
          ggtitle(cell)
      })
    p_ls
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path("Fibroblast_subtype/spatial_plot", str_glue("spatial plot {grp} celltype split.pdf")),
      height = 12,
      width  = 8
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })

##映射------
library(tibble)
library(data.table)
meta.data.t <-
  Fibroblasts[['subtype']] %>%
  rownames_to_column('cells')%>%
  right_join(
    all_data[['annotation4']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation4 = ifelse((as.character(annotation4)=='Fibroblast'),as.character(subtype), as.character(annotation4)))%>%
  dplyr::select(cells, annotation4) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation4'])
# meta.data <-
#   rbindlist(list(CD4T[['subtype_full']]%>%rownames_to_column('cells'),
#                  CD8T[['subtype_full']]%>%rownames_to_column('cells'),
#                  B[['subtype_full']]%>%rownames_to_column('cells'), 
#                  Fibroblast[['subtype_full']]%>%rownames_to_column('cells')))%>%as.data.frame()%>%
#   right_join(
#     all_data_cycle3[['annotation3']] %>% rownames_to_column('cells'),
#     by = "cells"
#   ) %>% 
#   mutate(
#     annotation4 = ifelse((as.character(annotation3)=='CD4-T'|as.character(annotation3)=='CD8-T'|as.character(annotation3)=='B'|as.character(annotation3)=='Epi0_B'|as.character(annotation3)=='CD83-dim B'|as.character(annotation3)=='CD83 B'|as.character(annotation3)=='Fibroblast'|as.character(annotation3)=='Epi0_Fibroblast'),as.character(subtype_full), as.character(annotation3)))%>%
#   dplyr::select(cells, annotation4) %>%
#   column_to_rownames('cells')
# 
# 
# all_data_cycle3 <- AddMetaData(all_data_cycle3, meta.data['annotation4'])

# cell props barplot  --------
Fibroblast_meta <-
  Fibroblast[[]] %>%
  rownames_to_column("cell_id")

Fibroblast_meta$project <- factor(Fibroblast_meta$project,levels=c("GX_1","GYF_1","HXF_1","LS_1","MXH_1","GX_2","GYF_2","HXF_2","LS_2","MXH_2","GX_3","GYF_3","HXF_3","LS_3","MXH_3"))
p_col_proj <-
  Fibroblast_meta %>%
  group_by(project, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = subtype)) + 
  geom_col() +
  ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

# export
ggsave(
  p_col_proj,
  filename = file.path("Fibroblast_subtype", "cell_freq_group.pdf"),
  height = 6,
  width  = 8
)

#macrophage subtype--------------------
Macrophage_out_dir <- "Macrophage_subtype"
dir.create(Macrophage_out_dir)
Macrophage <- subset(all_data, annotation=='Macrophage')
### PCA --------
Macrophage <- JoinLayers(Macrophage)
Macrophage <- NormalizeData(Macrophage, scale.factor = median(Macrophage$nCount_RNA))
Macrophage <- FindVariableFeatures(Macrophage)

Macrophage <- ScaleData(Macrophage, vars.to.regress = 'nCount_RNA')

Macrophage <- RunPCA(Macrophage)


### pca inspector ------
pdf(file.path(Macrophage_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(Macrophage, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  Macrophage[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(Macrophage)) %>%
  dplyr::pull("cell_id")


pdf(file.path(Macrophage_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(Macrophage_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(Macrophage[, suset_cells], dims = 16:20)
dev.off()
### dimension reduction --------
dims <- 16 # 18


library(future)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb
plan(multisession, workers = 10)


Macrophage <- RunUMAP(Macrophage, dims = 1:dims,return.model = TRUE)



pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


Macrophage <- Macrophage %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

Macrophage <- RunUMAP(
  Macrophage, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
Macrophage  <-
  FindNeighbors(
    Macrophage ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


Macrophage <-
  FindClusters(
    Macrophage,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
# DimPlot(Macrophage, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(Macrophage, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(Macrophage, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(Macrophage, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# Idents(B) <- "sketch_snn_res.0.6"
# B <- FindSubCluster(B, cluster = 7, graph.name = "sketch_snn", algorithm = 1,resolution = 0.3)
# Idents(B) <- "sub.cluster"
# DimPlot(B, group.by = "sub.cluster", label = TRUE, reduction = 'harmonyumap') + NoLegend()

### dotplot -------
major_marker <- c(
  'SPP1',
  M1 = 'NOS2',
  M1 = 'IL12',
  M1 = 'FCGR1A',
  M1 = 'FCGR1B',
  M1 = 'FCGR1C',
  M1 = 'CD80',
  M1 = 'CXCR10',
  M1 = 'IL23',
  M1 = 'CXCL9',
  M1 = 'CXCL10',
  M1 = 'CXCL11',
  M1 = 'CD86',
  M1 = 'IL1A',
  M1 = 'IL1B',
  M1 =  'IL6',
  M1 =  'TNF',
  M1 = 'MHCII',
  M1 = 'CCL5',
  M1 = 'IRF5',
  M1 = 'IRF1',
  M1 = 'CD40',
  M1 = 'KYNU',
  M1 = 'CCR7',
  M1 =  'VCAN',
  M1 =  'CD300E',
  M1 =  'INHBA',
  M1 =  'TNIP3',
  M2 = 'ARG1',
  M2 = 'IL10',
  M2 = 'CD32',
  M2 = 'CD163',
  M2 = 'FCER2',
  M2 = 'CD200R1',
  M2 = 'PDCD1LG2',
  M2 = 'CD274',
  M2 = 'MARCO',
  M2 = 'CSF1R',
  M2 = 'MRC1',
  M2 = 'IL1RN',
  M2 = 'Il1R2',
  M2 = 'IL4R',
  M2 = 'CCL4',
  M2 = 'CCL13',
  M2 = 'CCL20',
  M2 = 'CCL17',
  M2 = 'CCL18',
  M2 = 'CCL22',
  M2 = 'CCL24',
  M2 = 'LYVE1',
  M2 = 'VEGFA',
  M2 = 'VEGFB',
  M2 = 'VEGFC',
  M2 = 'VEGFD',
  M2 = 'EGF',
  M2 = 'CTSA',
  M2 = 'CTSB',
  M2 = 'CSTC',
  M2 = 'CTSD',
  M2 = 'TGFB1',
  M2 = 'CD209',
  M2 = 'CD163L1',
  M2 = 'FOLR2',
  cycling = "MKI67",
  cycling = 'TOP2A'
)




# dotplot
pdf(file.path(Macrophage_out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
# DotPlot(B, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
# DotPlot(B, group.by = "RNA_snn_res.0.4", features = major_marker_ls)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = celltype_marker)+RotatedAxis()
DotPlot(Macrophage, group.by = "RNA_snn_res.0.4", features = major_marker)+RotatedAxis()
# DotPlot(B, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
dev.off()



Idents(Macrophage)<-'RNA_snn_res.0.4'
markers <-
  FindAllMarkers(
    JoinLayers(Macrophage),
    test.use = "wilcox",
    logfc.threshold = 1,
    min.pct = 0.1,
    only.pos = TRUE
  )


markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.1) %>%
  dplyr::filter(abs(avg_log2FC) > 0.75) %>%
  # arrange(desc(avg_log2FC)) %>%
  arrange(desc(pct.1),.by_group = TRUE) %>%
  slice_head(n =20) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")


cell_cluster <- c(
  '0'  = "T_like_Macrophage",
  '1'  = "M2",
  '2'  = "M1",
  '3'  = "B_like_Macrophage",
  '4'  = "M1",
  '5'  = "Fibroblast_like_Macrophage",
  '6'  = "cycling macrophage",
  '7'  = "DC_like_Macrophage",
  '8'  = "Epithelial_like_Macrophage",
  '9'  = "Mast_like_Macrophage"
)

Idents(Macrophage) <- "RNA_snn_res.0.4"

Macrophage <- RenameIdents(Macrophage, cell_cluster)
Macrophage$subtype <- Idents(Macrophage)


pdf(file.path(Macrophage_out_dir, str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(Macrophage, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf(file.path(Macrophage_out_dir, "dotplot celltype main_marker.pdf"), width = 10, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(Macrophage, group.by = "subtype", features = major_marker, dot.scale = 10) + theme_90angle
DotPlot(Macrophage, group.by = "subtype", features = celltype_marker, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()

#spatial plot------------
Macrophage_meta <-
  Macrophage[[]] %>%
  rownames_to_column("cell_id")

sp_p_ls <-
  cells_boundary_ls2 %>% 
  imap(function(cells_boundary, project){  #browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat2 <-
      cells_boundary %>%
      right_join(Macrophage_meta, by = "cell_id")
    
    p <- dat[which(dat$annotation != 'Macrophage'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray80")  +
      scale_y_reverse() +
      coord_equal()+
      theme_void() +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat2, aes(x = vertex_x, y = vertex_y, fill = subtype, group = cell_id)) +
      scale_fill_manual(values = plan)+
      theme_void()
    p
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("Macrophage_subtype/spatial_plot", str_glue("Macrophage subtype spatial plot {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })
# 
# sp_p_ls <-
#   cells_boundary_ls %>%
#   imap(function(cells_boundary, project){ # browser()
#     image_wh_r <- image_wh_r_ls[[project]]
#     
#     dat <-
#       cells_boundary %>%
#       right_join(Macrophage_meta, by = "cell_id")
#     
#     p <- dat %>%
#       ggplot(aes(x = vertex_x, y = vertex_y)) +
#       geom_polygon(aes(fill = subtype, group = cell_id))  +
#       scale_y_reverse() +
#       coord_fixed(ratio = 1/image_wh_r) +
#       ggsci::scale_fill_igv() +
#       theme_void()
#     p
#   })

all_cells <- unique(Macrophage_meta$subtype)
cell_cols <- structure(
  ggsci::pal_igv()(length(all_cells)),
  names = all_cells
)

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){  #browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray60")  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      theme_void()
    
    p_ls <-
      all_cells %>%
      set_names(., .) %>%
      map(function(cell){ #browser()
        dat_sub <- dat %>% dplyr::filter(annotation3 == .env[['cell']])
        p +
          geom_polygon(
            data = dat_sub,
            aes(group = cell_id), 
            fill = cell_cols[cell]
          )  +
          ggtitle(cell)
      })
    p_ls
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path("Macrophage_subtype/spatial_plot", str_glue("spatial plot {grp} celltype split.pdf")),
      height = 12,
      width  = 8
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })
# cell props barplot  --------

Macrophage_meta <-
  Macrophage[[]] %>%
  rownames_to_column("cell_id")

Macrophage_meta$project <- factor(Macrophage_meta$project,levels=c("GX_1","GX_2","GX_3","MXH_1","MXH_2","MXH_3","HXF_1","HXF_2","HXF_3","LS_1","LS_2","LS_3","GYF_1","GYF_2","GYF_3"))

Macrophage_meta$project <- factor(Macrophage_meta$project,levels=c("GX_1","GYF_1","HXF_1","LS_1","MXH_1","GX_2","GYF_2","HXF_2","LS_2","MXH_2","GX_3","GYF_3","HXF_3","LS_3","MXH_3"))
p_col_proj <-
  Macrophage_meta %>%
  group_by(project, subtype) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = project, y = freq, fill = subtype)) + 
  geom_col() +
  scale_fill_manual(values = plan)+
  # ggsci::scale_fill_igv() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


p_box_group  <-
  all_data_meta %>%
  group_by(group, project, annotation_full) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation, y = freq, fill = group)) +
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
  filename = file.path("Macrophage_subtype", "cell_freq_group.pdf"),
  height = 6,
  width  = 8
)


##映射------
library(tibble)
library(dplyr)
meta.data.t <-
  Macrophage[['subtype']] %>%
  rownames_to_column('cells')%>%
  right_join(
    all_data[['annotation2']] %>% rownames_to_column('cells'),
    by = "cells"
  ) %>% 
  mutate(
    annotation2 = ifelse((as.character(annotation2)=='Macrophage'),as.character(subtype), as.character(annotation2)))%>%
  dplyr::select(cells, annotation2) %>%
  column_to_rownames('cells')


all_data<-AddMetaData(all_data,meta.data.t['annotation2'])

#EMT subtype-----------------------
EMT<-subset(all_data,annotation=='Epithelial')
EMT<-JoinLayers(EMT)
### PCA --------
DefaultAssay(EMT) <- "RNA"
EMT <- NormalizeData(EMT, scale.factor = median(EMT$nCount_RNA))
EMT <- FindVariableFeatures(EMT)

EMT <- ScaleData(EMT, vars.to.regress = 'nCount_RNA')

EMT <- RunPCA(EMT)

out_dir<-'EMT_subtype'
dir.create(out_dir)
### pca inspector ------
pdf(file.path(out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(EMT, ndims = 50)
dev.off()

set.seed(1234)
suset_cells <- 
  EMT[[]] %>%
  rownames_to_column("cell_id2") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(EMT)) %>%
  dplyr::pull("cell_id2")


pdf(file.path(out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(EMT[, suset_cells], dims = 1:5)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(EMT[, suset_cells], dims = 6:10)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(EMT[, suset_cells], dims = 11:15)
dev.off()
pdf(file.path(out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(EMT[, suset_cells], dims = 16:20)
dev.off()
### dimension reduction --------
dims <- 20


library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


EMT <- RunUMAP(EMT, dims = 1:dims,return.model = TRUE)



pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(EMT, group.by = "sample", label = TRUE) + NoLegend()
dev.off()


EMT <- EMT %>% 
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "sample",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

EMT <- RunUMAP(
  EMT, 
  dims = 1:dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(EMT, group.by = "sample", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



reduction.use <- 'harmony' # 'harmony' 'pca'


### cluster -------
EMT  <-
  FindNeighbors(
    EMT ,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )


EMT <-
  FindClusters(
    EMT,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )



# dimplot
pdf(file.path(out_dir, str_glue("Dimplot dim{dims} cluster.pdf")), height = 6, width = 8)
# DimPlot(All_T, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(All_T, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(All_T, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(All_T, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(EMT, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(EMT, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(EMT, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = 'harmonyumap') + NoLegend()
DimPlot(EMT, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

### dotplot ------
major_marker_ls<-c(
  "epithelial state" = "CDH1",
  "epithelial state" = "EPCAM",
  "epithelial state" = "CLDN3",
  "epithelial state" = "OCLN",
  "mesenchymal state" = "VIM",
  "mesenchymal state" = "CDH2",
  "mesenchymal state" = "FN1",
  "EMT" = "SNAI1",
  "EMT" = "SNAI2",
  "EMT" = "TWIST1",
  "EMT" = "TWIST2",
  "EMT" = "ZEB1",
  "EMT" = "ZEB2",
  "EMT" = "TGFB1",
  "EMT" = "MMPS"
)

celltype_marker<-c(
  "B cell"= "CD19", 
  "B cell"= "CD79A", 
  "B cell"= "CD79B", 
  "B cell"= "MS4A1",
  "Plasma" = "TNFRSF17", 
  "Plasma" = "MZB1",
  # "Plasma" = "XBP1",
  "T cell" = "CD3E", 
  "T cell" = "CD3G", 
  "T cell" = "CD4",
  "T cell" = "FOXP3", 
  "T cell" = "CD8A", 
  "T cell" = "GZMB",
  "Endothelial" = "PECAM1",
  "Endothelial" = "CLDN5", 
  "Epithelial" = "EPCAM",
  "Epithelial" = "KRT18",
  "Epithelial" = "CDH1",
  "Fibroblasts" = "SFRP4",
  "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "PLN",
  "Smooth muscle cells" = "MYLK",
  "Monocytes" = "ITGAM",
  "Monocytes" = "CD14", 
  "Monocytes" = "FCGR3A",
  "Monocytes" = "CD68",
  "Macrophage" = "CSF1R",
  "Macrophage" = "MARCO", 
  "DC"= "FCER1A", 
  "DC"= "CD1C", 
  "DC"= "CD83", 
  "DC"= "LILRA4", 
  "DC"= "CLEC9A", 
  "DC"= "IDO1",
  # 'pDC' = 'CLEC4C',
  # 'pDC' = "LILRA4",
  # 'pDC' = 'IL3RA',
  # 'pDC' = 'IRF8',
  # 'mregDC' = 'CCR7',
  # 'mregDC' ='CCL19',
  # 'mregDC' ='LAMP3',
  # 'mregDC' ='CD83',
  "NK" = "NCAM1",
  "Mast cell" = "CTSG", 
  "Mast cell" = "MS4A2", 
  "Mast cell" = "KIT", 
  "Mast cell" = 'HDC', 
  'Neutrophil'= "CSF3R", 
  'Neutrophil'= "FCGR3B", # "S100A8", "S100A9" 
  "Melanocy" = "MITF",
  "Melanocy" = "MLANA"
)







# dotplot
pdf(file.path(out_dir, str_glue("dotplot dims{dims} major_marker.pdf")), width = 25, height = 6)
DotPlot(EMT, group.by = "RNA_snn_res.0.2", features = celltype_marker)+RotatedAxis()
DotPlot(EMT, group.by = "RNA_snn_res.0.2", features = major_marker_ls)+RotatedAxis()
DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls)+RotatedAxis()
dev.off()
# find markers -------
Idents(All_T) <- "RNA_snn_res.0.6"
all_markers <-
  FindAllMarkers(
    JoinLayers(All_T),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  # arrange(desc(pct.1)) %>%
  slice_head(n = 30) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

### annotation ----

cell_cluster <- c(
  '0'  = "CD4_Tfh",
  '1'  = "Naive_T",
  '2'  = "Naive_T",
  '3'  = "CD8_Tem",
  '4'  = "Treg",
  '5'  = "CD8_naive",
  '6'  = "Naive_T",
  '7'  = "Tfh",
  '8'  = "Naive_T",
  '9'  = "Tfh",
  '10'  = "proliferating T",
  '11'  = "CD8_MAIT",
  '12'  = "CD4_naive",
  '13'  = "proliferating T",
  '14'  = "Naive_T",
  '15'  = "CD4_naive",
  '16'  = "Th2",
  '17'  = "proliferating T",
  '18'  = "proliferating T"
)


Idents(EMT) <- "RNA_snn_res.0.8"

EMT <- RenameIdents(EMT, cell_cluster)
EMT$subtype <- Idents(EMT)


pdf(file.path(out_dir,str_glue("Dimplot dim{dims} annotation.pdf")), height = 6, width = 8)
DimPlot(EMT, group.by = "subtype", label = TRUE, reduction = 'harmonyumap',repel = TRUE) + NoLegend()
# DimPlot(CD8T, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


pdf(file.path(out_dir,"dotplot celltype main_marker.pdf"), width = 20, height = 6)
# DotPlot(all_data, group.by = "RNA_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(EMT, group.by = "subtype", features = major_marker_ls_CD4_given, dot.scale = 10) + theme_90angle
DotPlot(EMT, group.by = "subtype", features = major_marker_ls_CD8, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "RNA_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle
dev.off()


### project data ------
EMT <- 
  EMT %>% 
  ProjectData(
    assay = "RNA",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony", # 'pca'
    
    umap.model = "harmonyumap", # 'umap'harmonyumap
    dims = 1:dims, 
    
    refdata = list(subtype_full = "subtype")
  )


DefaultAssay(EMT) <- 'RNA'

pdf(file.path("/home/rstudio/dat_analysis_cycle3/subtype/All_T",str_glue("Dimplot dim{dims} annotation - projected.pdf")), height = 6, width = 8)
DimPlot(EMT, group.by = "subtype_full", label = TRUE, reduction = 'full.harmonyumap.1') + NoLegend()
dev.off()

# saveRDS(All_T,"T_subtype/All_T_subtype.rds")
# all_data <- readRDS("CD8_subtype.rds")

#spatial plot------------
EMT_meta <-
  EMT[[]] %>%
  rownames_to_column("cell_id")


sp_p_ls <-
  cells_boundary_ls2 %>% 
  imap(function(cells_boundary, project){  #browser()
    # image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    dat2 <-
      cells_boundary %>%
      right_join(EMT_meta, by = "cell_id")
    
    p <- dat[which(dat$annotation != 'EMT'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray80")  +
      scale_y_reverse() +
      # coord_fixed(ratio = 1/image_wh_r) +
      coord_equal()+
      theme_void() +
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat2, aes(x = vertex_x, y = vertex_y, fill = subtype, group = cell_id)) +
      scale_fill_manual(values = plan)+
      theme_void()
    p
  })


plan <- c("#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
          "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", 
          "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", 
          "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879", 
          "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", 
          "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
          "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", "#ffc300", 
          "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
          "#f17d80", "#737495")


sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(file.path("EMT/spatial_plot", str_glue("EMT subtype spatial plot {grp}.pdf")), height = 18, width = 18)
    print(p)
    dev.off()
  })

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(All_T_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id))  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      ggsci::scale_fill_igv() +
      theme_void()
    p
  })

all_cells <- unique(All_T_meta$subtype)
cell_cols <- structure(
  ggsci::pal_igv()(length(all_cells)),
  names = all_cells
)

sp_p_ls <-
  cells_boundary_ls %>%
  imap(function(cells_boundary, project){  #browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(group = cell_id), fill = "gray90")  +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      theme_void()
    
    p_ls <-
      all_cells %>%
      set_names(., .) %>%
      map(function(cell){ #browser()
        dat_sub <- dat %>% dplyr::filter(annotation3 == .env[['cell']])
        p +
          geom_polygon(
            data = dat_sub,
            aes(group = cell_id), 
            fill = cell_cols[cell]
          )  +
          ggtitle(cell)
      })
    p_ls
  })

sp_p_ls %>% 
  iwalk(function(p, grp){
    pdf(
      file.path("T_subtype/spatial_plot", str_glue("spatial plot {grp} celltype split.pdf")),
      height = 12,
      width  = 8
    )
    # cat(paste0("#### ", grp, " {.unnumbered .tabshow} \n\n"))
    print(p)
    # cat("\n\n")
    dev.off()
  })




###score--------------------
M1<-c( 'NOS2',
       'IL12',
       'FCGR1A',
       'FCGR1B',
       'FCGR1C',
       'CD80',
       'CXCR10',
       'IL23',
       'CXCL9',
       'CXCL10',
       'CXCL11',
       'CD86',
       'IL1A',
       'IL1B',
        'IL6',
        'TNF',
       'MHCII',
       'CCL5',
       'IRF5',
       'IRF1',
       'CD40',
       'KYNU',
       'CCR7',
        'VCAN',
        'CD300E',
        'INHBA',
        'TNIP3')
M2<-c( 'ARG1',
       'IL10',
       'CD32',
       'CD163',
       'FCER2',
       'CD200R1',
       'PDCD1LG2',
       'CD274',
       'MARCO',
       'CSF1R',
       'MRC1',
       'IL1RN',
       'Il1R2',
       'IL4R',
       'CCL4',
       'CCL13',
       'CCL20',
       'CCL17',
       'CCL18',
       'CCL22',
       'CCL24',
       'LYVE1',
       'VEGFA',
       'VEGFB',
       'VEGFC',
       'VEGFD',
       'EGF',
       'CTSA',
       'CTSB',
       'CSTC',
       'CTSD',
       'TGFB1',
       'CD209',
       'CD163L1',
       'FOLR2')

pathway<-list(M1,M2)
names(pathway)<-c('M1','M2')

# M1_M2<-read.table("/home/lr/ATAC_thirteen/Immune/MG/M1_M2.txt",header=TRUE,fill=TRUE)
n_na_idx <- pathway %>% map_lgl(function(x) all(is.na(x))) %>% {! .}%>%  which
pathway <- pathway[n_na_idx] %>% map(function(x) x[!is.na(x)])

for ( i in seq_along(names(pathway))){
  Macrophage <- AddModuleScore(JoinLayers(Macrophage), pathway[i], name =  names(pathway)[i],nbin=5,ctrl=50)
}

Macrophage <- Macrophage[[]] %>% rownames_to_column('cell_id')
p<-VlnPlot(Macrophage,features = c('M1','M2'),group.by='RNA_snn_res.0.4',pt.size = 0)

ggsave(
  p,
  filename = 'Macrophage_subtype/M1_M2 score.pdf',
  height = 6,
  width  = 8
)
all_data_T <- all_data_T %>%
  mutate(subtype4 = case_when(
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))
sp_p_ls_tumor <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_T[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_T, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype4, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#CE3D32FF",'grey90'))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


# export image
sp_p_ls_tumor %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })
all_data_meta <- all_data_cycle3[[]] %>% rownames_to_column('cell_id')
#HLA--------
# sp_out_dir <- "dat_analysis_cycle3/subtype/newplot"

all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "CD4-Tmemory" | 
      as.character(annotation4) == "CD4-Tnaive" | 
      as.character(annotation4) == "CD8-Tem" | 
      as.character(annotation4) == "CD8-Tex" | 
      as.character(annotation4) == "CD8-Tnaive" | 
      as.character(annotation4) == "CD8-Trm" | 
      as.character(annotation4) == "Tfh" | 
      as.character(annotation4) == "Treg" | 
      as.character(annotation4) == "proliferating-T" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_HLA <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 'grey90',"#802268FF","#6BD76BFF","#D595A7FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = HLA, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_HLA %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot T {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "B-memory" | 
      as.character(annotation4) == "B-naive" | 
      as.character(annotation4) == "B-GC" | 
      as.character(annotation4) == "ASC" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_HLA <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF",'grey90',"#802268FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = HLA, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_HLA %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot B_plasma {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
      as.character(annotation4) == "apCAF" | 
      as.character(annotation4) == "iCAF" | 
      as.character(annotation4) == "mCAF"| 
      as.character(annotation4) == "vCAF"| 
      as.character(annotation4) == "rCAF"| 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_HLA <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF",'grey90',"#802268FF","#466983FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = HLA, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_HLA %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot fibroblast {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

#tumor_cell_LN_metastasis----------
all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "CD4-Tmemory" | 
      as.character(annotation4) == "CD4-Tnaive" | 
      as.character(annotation4) == "CD8-Tem" | 
      as.character(annotation4) == "CD8-Tex" | 
      as.character(annotation4) == "CD8-Tnaive" | 
      as.character(annotation4) == "CD8-Trm" | 
      as.character(annotation4) == "Tfh" | 
      as.character(annotation4) == "Treg" | 
      as.character(annotation4) == "proliferating-T" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_tumor_cell_LN_metastasis <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 'grey90',"#802268FF","#6BD76BFF","#D595A7FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = tumor_cell_LN_metastasis, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-1,0,1,2))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_tumor_cell_LN_metastasis %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot T {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "B-memory" | 
      as.character(annotation4) == "B-naive" | 
      as.character(annotation4) == "B-GC" | 
      as.character(annotation4) == "ASC" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_tumor_cell_LN_metastasis <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF",'grey90',"#802268FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = tumor_cell_LN_metastasis, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_tumor_cell_LN_metastasis %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot B_plasma {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "apCAF" | 
      as.character(annotation4) == "iCAF" | 
      as.character(annotation4) == "mCAF"| 
      as.character(annotation4) == "vCAF"| 
      as.character(annotation4) == "rCAF"| 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_tumor_cell_LN_metastasis <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF",'grey90',"#802268FF","#466983FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = tumor_cell_LN_metastasis, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_tumor_cell_LN_metastasis %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot fibroblast {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })






all_data_meta2$Inflammation[which(all_data_meta2$subtype4=='Others')]<-NA
sp_p_ls_Inflammation <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta2[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta2, by = "cell_id")
    
    p <- dat %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = Inflammation, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      # scale_fill_manual(values=c("#CE3D32FF",'grey90'))+
      scale_fill_gradient(low = "#FFB800", high ="#CC00FF",na.value = 'grey90')+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })

sp_p_ls_Inflammation %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


#IFN-------------
all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "CD4-Tmemory" | 
      as.character(annotation4) == "CD4-Tnaive" | 
      as.character(annotation4) == "CD8-Tem" | 
      as.character(annotation4) == "CD8-Tex" | 
      as.character(annotation4) == "CD8-Tnaive" | 
      as.character(annotation4) == "CD8-Trm" | 
      as.character(annotation4) == "Tfh" | 
      as.character(annotation4) == "Treg" | 
      as.character(annotation4) == "proliferating-T" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_IFN <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 'grey90',"#802268FF","#6BD76BFF","#D595A7FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = IFN, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-1,0,1,2))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_IFN %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot T {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "B-memory" | 
      as.character(annotation4) == "B-naive" | 
      as.character(annotation4) == "B-GC" | 
      as.character(annotation4) == "ASC" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_IFN <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF",'grey90',"#802268FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = IFN, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_IFN %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot B_plasma {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "apCAF" | 
      as.character(annotation4) == "iCAF" | 
      as.character(annotation4) == "mCAF"| 
      as.character(annotation4) == "vCAF"| 
      as.character(annotation4) == "rCAF"| 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_IFN <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF",'grey90',"#802268FF","#466983FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = IFN, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_IFN %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot fibroblast {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })

#Inflammation----------
all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "CD4-Tmemory" | 
      as.character(annotation4) == "CD4-Tnaive" | 
      as.character(annotation4) == "CD8-Tem" | 
      as.character(annotation4) == "CD8-Tex" | 
      as.character(annotation4) == "CD8-Tnaive" | 
      as.character(annotation4) == "CD8-Trm" | 
      as.character(annotation4) == "Tfh" | 
      as.character(annotation4) == "Treg" | 
      as.character(annotation4) == "proliferating-T" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_Inflammation <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 'grey90',"#802268FF","#6BD76BFF","#D595A7FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = Inflammation, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-1,0,1,2))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_Inflammation %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot T {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "B-memory" | 
      as.character(annotation4) == "B-naive" | 
      as.character(annotation4) == "B-GC" | 
      as.character(annotation4) == "ASC" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_Inflammation <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF",'grey90',"#802268FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = Inflammation, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_Inflammation %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot B_plasma {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "apCAF" | 
      as.character(annotation4) == "iCAF" | 
      as.character(annotation4) == "mCAF"| 
      as.character(annotation4) == "vCAF"| 
      as.character(annotation4) == "rCAF"| 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_Inflammation <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF",'grey90',"#802268FF","#466983FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = Inflammation, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_Inflammation %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot fibroblast {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })
#cancer_stem_cell_and_EMT-----------
all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "CD4-Tmemory" | 
      as.character(annotation4) == "CD4-Tnaive" | 
      as.character(annotation4) == "CD8-Tem" | 
      as.character(annotation4) == "CD8-Tex" | 
      as.character(annotation4) == "CD8-Tnaive" | 
      as.character(annotation4) == "CD8-Trm" | 
      as.character(annotation4) == "Tfh" | 
      as.character(annotation4) == "Treg" | 
      as.character(annotation4) == "proliferating-T" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_cancer_stem_cell_and_EMT <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF",
                                 'grey90',"#802268FF","#6BD76BFF","#D595A7FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = cancer_stem_cell_and_EMT, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4,6,8))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_cancer_stem_cell_and_EMT %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot T {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "B-memory" | 
      as.character(annotation4) == "B-naive" | 
      as.character(annotation4) == "B-GC" | 
      as.character(annotation4) == "ASC" | 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_cancer_stem_cell_and_EMT <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#749B58FF","#CE3D32FF",'grey90',"#802268FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = cancer_stem_cell_and_EMT, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4,6,8))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_cancer_stem_cell_and_EMT %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot B_plasma {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })


all_data_meta <- all_data_meta %>%
  mutate(subtype = case_when(
    as.character(annotation4) == "apCAF" | 
      as.character(annotation4) == "iCAF" | 
      as.character(annotation4) == "mCAF"| 
      as.character(annotation4) == "vCAF"| 
      as.character(annotation4) == "rCAF"| 
      as.character(annotation4) == "Epithelial" ~ as.character(annotation4),
    TRUE ~ "Others"
  ))

sp_p_ls_cancer_stem_cell_and_EMT <-
  cells_boundary_ls %>% 
  imap(function(cells_boundary, project){ # browser()
    image_wh_r <- image_wh_r_ls[[project]]
    
    anno_data   <- all_data_meta[c("project", "group", "group2")] %>% distinct()
    group_info  <- anno_data$group[anno_data$project == project]
    group_info2 <- anno_data$group2[anno_data$project == project]
    
    dat <-
      cells_boundary %>%
      right_join(all_data_meta, by = "cell_id")
    
    p <- dat[which(dat$subtype!='Epithelial'),] %>%
      ggplot(aes(x = vertex_x, y = vertex_y)) +
      geom_polygon(aes(fill = subtype, group = cell_id)) +
      scale_y_reverse() +
      coord_fixed(ratio = 1/image_wh_r) +
      # ggsci::scale_fill_igv() +
      scale_fill_manual(values=c("#5050FFFF","#CE3D32FF","#749B58FF",'grey90',"#802268FF","#466983FF"))+
      ggnewscale::new_scale_fill() +
      geom_polygon(data = dat[which(dat$subtype=='Epithelial'),],aes(fill = cancer_stem_cell_and_EMT, group = cell_id),color='black',linewidth=0.1) +
      scale_fill_gradient2(low = "cyan", high ="#FF00DB",midpoint = 0,mid = "white",breaks = c(-2,0,2,4,6,8))+
      theme_void() +
      labs(fill = NULL, title = paste(project, group_info, group_info2, sep = "_") )
    p
  })


sp_p_ls_cancer_stem_cell_and_EMT %>% 
  iwalk(function(p, name){
    pdf(
      file.path(sp_out_dir, str_glue("spatial plot fibroblast {name}.pdf")),
      width  = 20,
      height = 20
    )
    print(p)
    dev.off()
  })