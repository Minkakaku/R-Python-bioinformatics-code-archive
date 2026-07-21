library(Seurat)
library(tidyverse)

# setting dir -------
img_output <- "materials_for_rmd_cg" # dir.create(img_output)
qc_outdir <- "1. QC_cg" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster_cg"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot_cg" # dir.create(sp_output)
mk_output <- "4. marker gene_cg" # dir.create(mk_output)
er_output <- "5. enriment analysis_cg" # dir.create(er_output)

# data import -----------
all_cell_seg_ls <- 
  list.files("cell_segment", full.names = TRUE) %>% 
  list.files(full.names = TRUE)

all_data_ls <-  
  all_cell_seg_ls %>% 
  str_subset("rds$") %>%
  set_names(., basename(dirname(.))) %>%
  map(readRDS)


all_gdf_ls <- 
  all_cell_seg_ls %>% 
  str_subset("gdf") %>%
  set_names(., basename(dirname(.))) %>%
  map(read.csv)


all_gdf_ls <-
  all_gdf_ls %>%
  imap(function(gdf, name){ # browser()
    sc <- all_data_ls[[name]]
    
    gdf[gdf$id %in% colnames(sc), ]
  })


# cell area filter ---------------
all_gdf_ls %>%
  map('area') %>%
  enframe(name = "sample") %>%
  unnest(value) %>%
  ggplot(aes(x = value, fill = sample)) +
  geom_density(alpha = 0.5) +
  theme_bw()
  
# "838375-1-181" 2000
# "838375-2-91"  2000

cell_area_cutoff_df <- 
  data.frame(
    sample = c("838375-1-181", "838375-2-91"),
    cutoff = c(2000, 2000)
  )


all_gdf_filter_ls <- 
  all_gdf_ls %>%
  imap(function(gdf, name){
    cutoff <- cell_area_cutoff_df$cutoff[cell_area_cutoff_df$sample == name]
    gdf[gdf$area < cutoff, ]
  })
  

all_data_filter_ls <-
  all_data_ls %>%
  imap(function(sc, name){
    gdf <- all_gdf_ls[[name]]
    sc[, colnames(sc) %in% gdf$id]
  })


# parse cell_coords  ---------------
parse_boundary <- function(x){ # browser()
  stopifnot(length(x) == 1)
  stopifnot(str_detect(x, "^POLYGON"))
  
  xx <- x %>% str_remove("POLYGON \\(\\(") %>% str_remove("\\)\\)")
  
  res <- xx %>% 
    str_split_1(",") %>%
    str_trim %>% 
    str_split(pattern = "\\s+") %>% 
    map(as.numeric) %>% 
    unlist %>%
    matrix(byrow = TRUE, ncol = 2)
  

  colnames(res) <- c("x", "y")
  as.data.frame(res)
}

parse_boundary_multi <- function(x){
  res_ls <- list()
  for( i in seq_along(x)){
    res <- parse_boundary(x[i])
  }
  res_ls[[i]] <- res
  res_ls
}

all_gdf_ls$`838375-1-181`$geometry %>% head(1)  %>% parse_boundary()

coords_1_181 <- all_gdf_ls$`838375-1-181`
parsered_coords_1_181 <- coords_1_181$geometry %>% parse_boundary_multi()

coords_1_181$boundary <- coords_1_181$geometry %>% parse_boundary_multi()


visium <- Load10X_Spatial("/home/ztt/spatial/visium_HD/FWSP20240131", bin.size = 8, slice = name)


# QC -----------

all_data_filter_ls <- 
  all_data_filter_ls %>%
  map(function(x){# browser()
    x[["nFeature_RNA"]] <- colSums(x[['RNA']]@counts != 0)
    x[["nCount_RNA"]]   <- colSums(x[['RNA']]@counts)
    
    
    x[["percent.mt"]]    <- PercentageFeatureSet(x, pattern = "^MT-")
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = "^RP[SL]")
    x$log10GenesPerUMI   <- log10(x$nFeature_RNA + 1) / log10(x$nCount_RNA + 1)
    x
  })


p_qc_ls <-
  all_data_filter_ls %>%
  imap(function(x, name){ # browser()
    # only subset 10,000 cells
    x_sub <- subset(x, cells = Cells(x[['RNA']]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI"), 
      ncol     = 4, 
      pt.size  = 0,
      log      = FALSE, 
      raster   = TRUE
    ) +
      patchwork::plot_annotation(title = name)
    
    # 精细划分y轴标度
    p & scale_y_continuous(
      breaks = function(x){
        xx <- round(x)
        
        # 5代表分为5个breaks
        # 10代表以10的倍数进行划分breaks
        steps <- if(xx[2] < 100){
          scaled_times <- ceiling(100/xx[2])
          
          floor((scaled_times * xx[2] / 5 ) / 10 ) * 10 / scaled_times
        }else{
          floor((xx[2] / 5 ) / 10 ) * 10
        }
        
        # 最小值都是0
        bk <- seq(0, xx[2], by = steps)
        
        # 第一间隔细分
        bk_minor <- seq(bk[1], bk[2], length.out = 5)
        bk <- c(bk_minor, bk[-(1:2)])
        
        # 最后间隔细分
        bk <- rev(bk)
        bk_minor <- seq(bk[1], bk[2], length.out = 5)
        bk <- c(bk_minor, bk[-(1:2)])
        bk <- rev(bk)
        
        bk
      }
    )
  })

p_qc_log_ls <-
  all_data_filter_ls %>%
  imap(function(x, name){
    # only subset 10,000 cells
    x_sub <- subset(x, cells = Cells(x[['RNA']]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI"), 
      ncol     = 4, 
      pt.size  = 0,
      log      = TRUE,
      raster   = TRUE
    ) +
      patchwork::plot_annotation(title = name)
    
  })


# export
# saveRDS(p_qc_ls, file = file.path(rds_output, "1. p_qc_ls.rds"))

p_qc_ls %>%
  iwalk(function(p, name){
    png(file.path(img_output, str_glue("1 qc plot - {name}.png")), width = 480 * 3, height = 480 )
    print(p)
    dev.off()
  })


pdf(
  file.path(qc_outdir, "qc plot.pdf"),
  height = 6, 
  width = 20
)
p_qc_ls %>% walk(print)
dev.off()

pdf(
  file.path(qc_outdir, "qc plot - log.pdf"),
  height = 6, 
  width = 20
)
p_qc_log_ls %>% walk(print)
dev.off()

# qc filter ----------
qc_metric <- read.csv("qc_metric_cg.csv") %>% map_dfc(str_trim)

qc_metric_variable <- c(
  "nFeature_RNA", "nCount_RNA",
  "percent.mt", "log10GenesPerUMI"
)
qc_metric_ls <-
  qc_metric %>% 
  apply(1, name = names(.), simplify = FALSE, function(x, name = name){
    qc_ls <- map2(x, name, function(num, variable){
      if(!variable %in% qc_metric_variable) return()
      
      if(str_detect(num, ":")){
        num <- str_split(num, pattern = ":", n = 2) %>% .[[1]]
        qc <- str_glue("{variable} {num[1]} & {variable} {num[2]}")
      }else{
        qc <- str_glue("{variable} {num}")
      }
      qc
    })
    
    res <- qc_ls %>% discard(function(x) length(x) < 1) %>% str_c(collapse = " & ")
    res
  })

names(qc_metric_ls) <- qc_metric$sample



# filter
all_data_filter_ls_pass <-
  names(all_data_filter_ls) %>%
  set_names(.,.) %>%
  map(function(sample){ #browser()
    all_data <- all_data_filter_ls[[sample]]
    qc <- qc_metric_ls[[sample]]
    
    e <- str_glue("subset(
      all_data,
      {qc}
    )")
    
    eval(parse(text = e))
  })


# qc plot
qc_plot_filter_ls <-
  names(all_data_filter_ls) %>%
  setNames(., .) %>%
  map(function(sample){
    all_data <- all_data_filter_ls[[sample]] %>% subset(downsample = 100000, seed = 1)
    all_data_pass <- all_data_filter_ls_pass[[sample]] %>% subset(downsample = 100000, seed  = 1)
    
    # all_data
    plot1 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = "nCount_RNA", feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + patchwork::plot_layout(ncol = 3)
    
    # all_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_RNA", feature2 = "log10GenesPerUMI")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    p2 <- plot1_pass + plot2_pass + plot3_pass  + patchwork::plot_layout(ncol = 3)
    
    p1 / p2
  })



# export for rmd
qc_plot_filter_ls %>%
  iwalk(function(p, name){
    png(file.path(img_output, str_glue("1 qc plot - pass - {name}.png")), width = 360 * 3, height = 360 * 2 )
    print(p)
    dev.off()
  })


# export 
pdf(
  file.path(qc_outdir, "qc plot - pass.pdf"),
  height = 12, 
  width  = 18
)
qc_plot_filter_ls %>% walk(print)
dev.off()

qc_metric_ls %>%
  enframe("sample", "qc") %>%
  unnest(qc) %>% 
  openxlsx::write.xlsx(
    file = file.path(qc_outdir, "qc_metric.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )

# cells filter
cell_counts_df <-
  names(all_data_filter_ls) %>%
  set_names(., .) %>%
  map(function(sample){
    all_data <- all_data_filter_ls[[sample]]
    all_data_pass <- all_data_filter_ls_pass[[sample]]
    
    structure(
      c(ncol(all_data), ncol(all_data_pass)),
      names = c("Before filter", "After filter")
    )
  }) %>%
  as.data.frame(check.names = FALSE)


# export
cell_counts_df %>% 
  rownames_to_column(" ") %>%
  openxlsx::write.xlsx(
    file = file.path(qc_outdir, "cell count after filter.xlsx"),
    asTable = TRUE,
    withFilter = FALSE
  )


# dimension reduction and cluster and annotation ------
all_data_filter_ls_pass <- 
  all_data_filter_ls_pass %>%
  imap(function(sc, name){
    sc$project <- name
    sc
  })

# all_data <- merge(all_data_filter_ls_pass[[1]], all_data_filter_ls_pass[-1])
all_data <- all_data_filter_ls_pass$`838375-2-91`

## normalize for sketch data ------
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(all_data)


## sketch data --------
all_data <-
  SketchData(
    all_data,
    ncells         = 50000L,
    sketched.assay = "sketch",
    method         = "LeverageScore",
    var.name       = "leverage.score",
    seed = 123L
  )


### dimension reduction
DefaultAssay(all_data) <- "sketch"
all_data <- FindVariableFeatures(all_data, nfeatures = 3000)
all_data <- ScaleData(
  all_data,
  vars.to.regress = c(
    'nCount_RNA'
  )
)

# all_data <- SCTransform(all_data, assay = "sketch")

all_data <- RunPCA(all_data, features = VariableFeatures(all_data))

pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()


set.seed(1234)
suset_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(dr_outdir, 'pca dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 1:5)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(dr_outdir, 'pca dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 16:20)
dev.off()


## dimension reduction --------
dims <- 6

all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE, assay  = 'sketch')


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()
dev.off()


all_data <- all_data %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data <- RunUMAP(
  all_data,
  dims = 1:dims,
  return.model = TRUE,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'pca' # 'harmony' 'pca'
## cluster ------
all_data <-
  FindNeighbors(
    all_data,
    reduction = reduction.use,
    dims = 1:dims,
    k.param = 20
  )

all_data <-
  FindClusters(
    all_data,
    resolution = c(0.2, 0.4, 0.6, 0.8)
  )


# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()

DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()

# DimPlot(all_data, group.by = "SCT_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()



## dotplot -------
major_marker_ls2 <- c(
  Hep = 'CYP3A7',
  Hep = 'CYP2A7',
  Hep = 'CYP2A6',
  
  Hep = 'SCD',
  Hep = 'HMGCS1',
  Hep = 'ACSS2',
  Hep = 'TM7SF2',
  
  Hep = 'SEC16B',
  Hep = 'SLBP',
  Hep = 'RND3',
  Hep = 'PCK1',
  
  Hep = 'BCHE',
  Hep = 'G6PC',
  Hep = 'GHR',
  Hep = 'ALDH6A1',
  
  Hep = 'APP25L',
  Hep = 'HSD11B1',
  Hep = 'HAMP',
  # Hep = 'GHR',
  
  Hep = 'HPR',
  Hep = 'GSTA2',
  Hep = 'AKR1C1',
  Hep = 'MASP2',
  
  Zone1 = 'MGP',
  Zone1 = 'SPARCL1',
  Zone1 = 'TM4SF1',
  Zone1 = 'CLEC14A',
  
  Zone2_3 = 'CCL14',
  Zone2_3 = 'CLEC1B',
  Zone2_3 = 'FCN2',
  Zone2_3 = 'S100A13',
  
  Portal = 'RAMP3',
  Portal = 'INMT',
  Portal = 'DNASE1L3',
  Portal = 'LIFR',
  
  Chol = 'KRT7',
  Chol = 'KET19',
  Chol = 'SOX9',
  Chol = 'EPCAM',
  
  Stellate = 'ACTA2',
  Stellate = 'COL1A1',
  Stellate = 'ABP1',
  
  Inf_macs = 'S100A8',
  Inf_macs = 'LYZ',
  Inf_macs = 'S100A9',
  Inf_macs = 'HLA-DPB1',
  
  Non_inf_macs = 'CD5L',
  Non_inf_macs = 'MARCO',
  Non_inf_macs = 'VSIG4',
  
  CD3T = 'CD2',
  CD3T = 'CD3D',
  CD3T = 'TRAC',
  CD3T = 'GZMK',
  
  gamma_delta_T = 'GNLY',
  gamma_delta_T = 'PTGDS',
  gamma_delta_T = 'GZMB',
  gamma_delta_T = 'TRDC',
  
  gamma_delta_T = 'STMN1',
  gamma_delta_T = 'HMGB2',
  gamma_delta_T = 'TYMS',
  
  B = 'MS4A1',
  B = 'LTB',
  B = 'CD37',
  B = 'CD79B',
  
  Plasma = 'IGLC2',
  Plasma = 'IGHG1',
  Plasma = 'IGKC',
  
  NK = 'CD7',
  NK = 'KLRB1',
  NK = 'NKG7',
  
  Erythroid = 'HBB',
  Erythroid = 'CA1',
  Erythroid = 'ALAS2'
)


major_marker_ls <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  Epithelial  = "EPCAM", 
  Epithelial  = "CEACAM5", 
  Epithelial  = "CEACAM6", 
  Cholangiocyte  = "KRT7",
  Cholangiocyte  = "KRT19",
  Hepatocyte = 'ALB',
  Hepatocyte = 'APOC3',
  
  "Fibroblasts" = "SFRP4", 
  # "Fibroblasts" = "PDGFRA", 
  "Fibroblasts" = "CCDC80",
  "Fibroblasts" = "COL1A1",
  CAF = "ACTA2",
  CAF = "COL1A2",
  iCAF = "CXCL12",
  iCAF = 'PDGFRA',
  mCAF = 'RGS5',
  mCAF = 'TAGLN',
  
  # "Pericytes" = "PDGFRB",
  
  
  "Smooth muscle cells" = "PDGFRB", 
  "Smooth muscle cells" = "MYH11", 
  # "Smooth muscle cells" = "ACTA2",
  
  # "Breast glandular cells" = "SERPINA3", 
  # "Breast glandular cells" = "DSC2",
  
  Immune = 'PTPRC',
  
  "B cells" = "CD19",
  "B cells" = "CD79A",
  "B cells" = "MS4A1",
  "Plasma" = "IGHG1",
  "Plasma" = "IGHA1",
  "Plasma" = "IGHG3",
  "Plasma" = "IGKC",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  
  Plasmablast = "MZB1",
  Plasmablast = "XBP1",
  
  "T cells" = "CD3E",
  "T cells" = "CD3G",
  "T cells" = "CD4",
  "T cells" = "CD8A",
  "T cells" = "GZMK",
  
  "NK" = "GNLY",
  "NK" = "NKG7",
  "NK" = "NCAM1",
  
  "Proliferating" = "MKI67",
  "Proliferating" = "TOP2A",
  
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
  
  'Fat cell' = 'ADIPOQ'
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} main_marker.pdf")), width = 25, height = 6)
# DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker_ls, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle

DotPlot(all_data, group.by = "RNA_snn_res.0.2", dot.scale = 10, features = major_marker_ls) + theme_90angle
DotPlot(all_data, group.by = "RNA_snn_res.0.4", dot.scale = 10, features = major_marker_ls) + theme_90angle
DotPlot(all_data, group.by = "RNA_snn_res.0.6", dot.scale = 10, features = major_marker_ls) + theme_90angle
DotPlot(all_data, group.by = "RNA_snn_res.0.8", dot.scale = 10, features = major_marker_ls) + theme_90angle

# DotPlot(all_data, group.by = "SCT_snn_res.0.2", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.4", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.6", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.8", dot.scale = 10, features = major_marker_ls) + theme_90angle
dev.off()




pdf(
  file.path(dr_outdir, str_glue("featureplot dims{dims} {reduction.use} main_marker.pdf")), 
  width = 25, 
  height = 50
)
FeaturePlot(
  all_data,
  features = major_marker_ls,
  ncol = 5,
  raster = TRUE
)  
dev.off()


## annotation ----
Idents(all_data) <- "sketch_snn_res.0.4"


# cell_cluster <- c(
#   '0'  = "Hepatocyte",
#   '1'  = "B",
#   '2'  = "Hepatocyte",
#   '3'  = "Fibroblasts",
#   '4'  = "T",
#   '5'  = "Plasma",
#   '6'  = "Hepatocyte",
#   '7'  = "Hepatocyte"
# )

cell_cluster <- c(
  '0'   = "Hepatocyte",
  '1'   = "B",
  '2'   = "Fibroblasts",
  '3'   = "Hepatocyte",
  '4'   = "T",
  '5'   = "Fibroblasts",
  '6'   = "Plasma",
  '7'   = "Fibroblasts",
  '8'   = "Cholangiocyte",
  '9'   = "Hepatocyte",
  '10'  = "Hepatocyte"
)


all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)


pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation {reduction.use} .pdf")), height = 6, width = 8)
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# pdf(str_glue(str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
# DotPlot(all_data, group.by = "annotation", features = major_cell_marker)
# dev.off()

