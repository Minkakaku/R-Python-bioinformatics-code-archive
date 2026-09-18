library(Seurat)
library(tidyverse)

# setting dir -------
img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)
mk_output <- "4. marker gene" # dir.create(mk_output)
er_output <- "5. enriment analysis" # dir.create(er_output)

# 4 11 20 hvg5000
# 12 hvg2000 regression out H cell


#### flatDoubleList
flatDoubleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  res_ls <- list()
  
  # layer 1
  name_1 <- names(x)
  id_1   <- seq_along(x)
  if(is.null(name_1)) name_1 <- id_1
  
  for(i in id_1){
    y <- x[[i]]
    stopifnot(is.vector(y))
    
    # layer 2
    name_2 <- names(y)
    id_2   <- seq_along(y)
    if(is.null(name_2)) name_2 <- id_2
    
    for(j in id_2){
      res  <- y[[j]]
      name <- paste0(name_1[i], sep, name_2[j])
      res_ls[[name]] <- res
    }
  }
  res_ls
}

##### nestSingleList
nestSingleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  name <- names(x)
  id_ <- seq_along(name)
  stopifnot(all(str_detect(name, sep)))
  name_1 <- str_split(name, pattern = sep) %>% map_chr(~ .x[1])
  name_2 <- str_split(name, pattern = sep) %>% map(~ .x[2:length(.x)]) %>% map_chr(str_c, collapse = sep)
  
  res_ls <- list()
  for(i in id_){    
    n1 <- name_1[i]
    n2 <- name_2[i]
    res <- x[[i]]
    res_ls[[n1]][[n2]] <- res
  }
  res_ls
}

# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


# set default -------------
bin.size     <- '008um' # '016um' 
assay.name   <- str_glue("Spatial.{bin.size}")
image.name   <- "tissue_hires_image.png" # "tissue_lowres_image.png"

stopifnot(bin.size %in% c('008um', '016um' ))

# read data -------
sample_df <- read.csv("data/Annoation.csv")
  # data.frame(
  #   sample = c('838375-1', '838375-2'),
  #   group  = c('Tumor', 'Handover')
  # )
  
file_list <- 
  list.files(
    'data/FWSP20240200', # "/home/ztt/spatial/visium_HD/FWSP20240200",
    pattern = paste(sample_df$sample, collapse = "|"), 
    full.names = TRUE
  ) %>%
  str_subset("_[no|old|remove]", negate = TRUE)

names(file_list) <- basename(file_list) # sample_df$group # basename(file_list)


visim_ls <-
  file_list %>%
  imap(function(path, name){
    path <- file.path(path, "outs")
    spatial_path <- file.path(path, "binned_outputs", str_glue("square_{bin.size}"), 'spatial')
    
    # read hi-res/low-res
    image <- Read10X_Image(
      spatial_path,
      assay      = assay.name,
      slice      = name,
      image.name = image.name,
    )
    
    # visium
    Load10X_Spatial(
      path,
      bin.size = 8, 
      image = image,
      slice = name
    )
  })

# parse ROI ---------------------
#1. read coordinate of  all cells  in all ROIs
roi_ls <- list.files(
  'data/ROI_coordinates/',
  pattern    = "Projection", 
  full.names = TRUE
) 

names(roi_ls) <- 
  roi_ls %>%
  basename() %>%
  str_remove("-.*$")

roi_coords_ls <-
  roi_ls %>%
  map(function(path){
    dat <- read.csv(path)
    dat <- 
      dat %>%
      dplyr::rename(
        X = X.Coordinate,
        Y = Y.Coordinate
      )
    dat
  }) 

# compute coords
compute_coords <- FALSE

if(compute_coords){
  #2. compute ROI boundary from ROI coordinate
  roi_boundary_ls <-
    roi_coords_ls %>%
    map(function(x){
      idx <- chull(x[c("X", "Y")])
      dat <- x[idx, c("X", "Y")]
      rownames(dat) <- NULL
      dat
    })
  
  
  #3. get all coordinate of all ROIs
  # roi_ls have all coordinates, why select all points' coordinate manually?
  # 3.1 only 008um have coordinate, not 016um
  # 3.2 coordinates manualy computed were not all identical to roi_ls which download from Loupe
  
  roi_coords_recomputed_ls <-
    roi_boundary_ls %>%
    imap(function(bound_coords, name){ # browser()
      visium_name <- sample_df$sample[sample_df$group == name]
      visium <- visim_ls[[visium_name]]
      
      # all points
      coords <- GetTissueCoordinates(visium)
      
      # compute all_points
      # sp_idx: 
      # 0: point exterior to pol 
      # 1: point interior to pol
      # 2: point lies on a line
      # 3: point lies on a vertex
      sp_idx <-
        sp::point.in.polygon(
          point.x = coords$x,
          point.y = coords$y,
          pol.x = bound_coords$X, 
          pol.y = bound_coords$Y
        )
      coords[as.logical(sp_idx), ]
    })
}else{
  roi_coords_recomputed_ls <- 
    roi_coords_ls %>%
    map(function(x){
      colnames(x) <- x %>% colnames() %>% str_to_lower()
      x$cell <- x$barcode # rownames(x)
      x
    })
}


#4. subset Visium Seurat Object
visium_roi_ls <-
  roi_coords_recomputed_ls %>%
  imap(function(roi, name){ # browser()
    # get Visium Object
    visium_name <- sample_df$sample[sample_df$group == name]
    visium <- visim_ls[[visium_name]]
    
    # subset visium
    cells <- roi$cell
    visium_roi <- subset(visium, cells = cells)
    
    # subset image
    image_name <- make.names(str_glue("{visium_name}.{bin.size}"))
    image <- visium_roi@images[[image_name]]
    stopifnot(!is.null(image))
    Key(image) <- suppressWarnings(Key(name))
    
    # change image names
    visium_roi@images[[image_name]] <- NULL
    visium_roi@images[[name]] <- image
    
    visium_roi
  })

# # test spatial plot
# SpatialFeaturePlot(
#   visium_roi_ls[[1]], 
#   image.scale = 'hires', 
#   pt.size.factor = 5,
#   shape = 22,
#   alpha = c(0.5, 0.8),
#   features = "nFeature_Spatial.008um"
# )

# normalize
visium_roi_ls <- 
  visium_roi_ls %>%
  map(~ {DefaultAssay(.x) <- 'Spatial.008um'; .x}) %>%
  map(NormalizeData) 

# add sample_id or sample_name -------------------
visium_roi_ls <-
  visium_roi_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })

# visium_roi override visium ------------
visim_ls <- visium_roi_ls

# add sample_id or sample_name -------------------
# visim_ls <-
#   visim_ls %>%
#   imap(function(visium, name){
#     visium$project <- name
#     Idents(visium) <- 'project'
#     RenameCells(visium, add.cell.id = name)
#   })


visim_ls <- 
  visim_ls %>%
  map(~ {DefaultAssay(.x) <- 'Spatial.008um'; .x}) %>%
  map(NormalizeData) 


# qc -----
visim_ls <- 
  visim_ls %>%
  map(function(visium){
    subset(visium, nCount_Spatial.008um > 0)
  })

visim_ls <- 
  visim_ls %>%
  map(function(x){# browser()
    x[["percent.mt"]]    <- PercentageFeatureSet(x, pattern = "^MT-")
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = "^RP[SL]")
    x$log10GenesPerUMI   <- log10(x$nFeature_Spatial.008um) / log10(x$nCount_Spatial.008um)
    x
  })


p_qc_ls <-
  visim_ls %>%
  imap(function(x, name){ # browser()
    # only subset 10,000 cells
    x_sub <- subset(x, cells = Cells(x[['Spatial.008um']]), downsample = 100000, seed = 1)
    p <- VlnPlot(
        x_sub, 
        features = c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mt", "log10GenesPerUMI"), 
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
  visim_ls %>%
  imap(function(x, name){
    # only subset 10,000 cells
    x_sub <- subset(x, cells = Cells(x[['Spatial.008um']]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c("nFeature_Spatial.008um", "nCount_Spatial.008um", "percent.mt", "log10GenesPerUMI"), 
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
qc_metric <- read.csv("qc_metric.csv") %>% map_dfc(str_trim)

qc_metric_variable <- c(
  "nFeature_Spatial.008um", "nCount_Spatial.008um",
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
visim_ls_pass <-
  names(visim_ls) %>%
  set_names(.,.) %>%
  map(function(sample){ #browser()
    all_data <- visim_ls[[sample]]
    qc <- qc_metric_ls[[sample]]
    
    e <- str_glue("subset(
      all_data,
      {qc}
    )")
    
    eval(parse(text = e))
  })


# qc plot
qc_plot_filter_ls <-
  names(visim_ls) %>%
  setNames(., .) %>%
  map(function(sample){
    all_data <- visim_ls[[sample]] %>% subset(downsample = 100000, seed = 1)
    all_data_pass <- visim_ls_pass[[sample]] %>% subset(downsample = 100000, seed  = 1)
    
    # all_data
    plot1 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.008um", feature2 = "nFeature_Spatial.008um") +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.008um", feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.008um", feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + patchwork::plot_layout(ncol = 3)
    
    # all_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.008um", feature2 = "nFeature_Spatial.008um")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.008um", feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.008um", feature2 = "log10GenesPerUMI")+
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
  names(visim_ls) %>%
  set_names(., .) %>%
  map(function(sample){
    all_data <- visim_ls[[sample]]
    all_data_pass <- visim_ls_pass[[sample]]
    
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
all_data <- merge(visim_ls_pass[[1]], visim_ls_pass[-1])

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
    'nCount_Spatial.008um'
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
dims <- 30

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


reduction.use <- 'harmony' # 'harmony' 'pca'
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
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()

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
major_marker_ls <- c(
  Myoid = "ACTA2",
  Myoid = "MYH11",
  Myoid = "PTCH1",
  Myoid = "PTCH2",
  
  Sertoli = "WT1",
  Sertoli = "SOX9",
  Sertoli = "BEX2",
  Sertoli = "GATA4",
  
  Macrophage = "CD68",
  Macrophage = "MS4A6A",
  Macrophage = "LYZ",
  Macrophage = "IL3RA",
  
  Leydig = "CYP11A1",
  Leydig = "STAR",
  Leydig = "IGF1",
  Leydig = "IGFBP5",
  
  Endothelial = "VWF",
  Endothelial = "PECAM1",
  Endothelial = "CD34",
  Endothelial = "TGFBR2",
  
  'T cell' = "CD3D",
  'T cell' = "CD3E",
  'T cell' = "CD3G",
  'T cell' = "CD2",
  
  Spermatogonia = "FGFR3",
  Spermatogonia = "STRA8",
  Spermatogonia = "DMRT1",
  Spermatogonia = "UCHL1",
  Spermatogonia = "TCF3",
  
  Spermatocytes = "SYCP1",
  Spermatocytes = "SMC1B",
  Spermatocytes = "SYCP3",
  Spermatocytes = "OVOL2",
  
  Spermatids = "TNP1",
  Spermatids = "PRM1",
  Spermatids = "TXNDC2",
  Spermatids = "PRM3"
)

theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} main_marker.pdf")), width = 25, height = 6)
DotPlot(all_data, group.by = "sketch_snn_res.0.2", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.4", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker_ls, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls, dot.scale = 10) + theme_90angle


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


## project data ------
all_data <- 
  all_data %>% 
  ProjectData(
    assay = "Spatial.008um",
    full.reduction = "pca.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony",
    
    umap.model = "harmonyumap", 
    dims = 1:dims, 
    
    refdata = list(
      annotation_full = "annotation",
      sketch_snn_res.0.4_full = "sketch_snn_res.0.4"
    )
  )

Idents(all_data) <- 'annotation_full'
DefaultAssay(all_data) <- "Spatial.008um"

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "sketch_snn_res.0.4_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

# saveRDS(all_data, "all_data_projected.rds")
# all_data <- readRDS("all_data_projected.rds")

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2 projected.pdf")), width = 25, height = 6)
DotPlot(all_data, assay = "sketch", group.by = "annotation_full", features = major_marker_ls) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()

# export for rmd
png(file.path(img_output, str_glue("2 dimplot - cluster.png")), width = 600, height = 480)
  DimPlot(
    all_data, 
    group.by = "sketch_snn_res.0.4_full", 
    label = TRUE, 
    reduction = 'full.harmonyumap'
  ) +
  NoLegend()
dev.off()

## find markers -------
Idents(all_data) <- "sketch_snn_res.0.4"
DefaultAssay(all_data) <- "sketch"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

# export marker gene
all_markers %>%
  dplyr::select(gene, everything()) %>%
  write.csv(
    file = file.path(mk_output, "marker.csv"), 
    row.names = FALSE
  )

all_de_genes <-
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  pull("gene")

all_de_genes_for_enrichmenet <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% 
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map("gene")


all_top5_gene <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  with(structure(
    gene,
    names = as.character(cluster)
  ))
  

# top5 dotplot
p_dot_top5 <-
  DotPlot(
    all_data, 
    features = all_top5_gene
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave(
  p_dot_top5,
  filename = file.path(mk_output, "top 5 gene dotplot.pdf"),
  height = 6,
  width  = 12
)

# export for rmd
png(file.path(img_output, str_glue("4 top5 gene dotplot.png")), width = 960, height = 480)
print(p_dot_top5)
dev.off()


# top1 gene
all_top1_gene <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>% ##
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 1)%>%
  with(structure(
    gene,
    names = as.character(cluster)
  ))

# top1 gene umap
pdf(
  file.path(mk_output, str_glue("top 1 gene featureplot.pdf")), 
  width = 25, 
  height = 5 * ceiling(length(all_top1_gene) / 5)
)
FeaturePlot(
  all_data,
  features = all_top1_gene,
  reduction = "harmonyumap",
  ncol = 5,
  raster = TRUE
)  
dev.off()

# export for rmd
png(file.path(img_output, str_glue("4 top1 gene umap.png")), width = 400 * 5, height = 320 * ceiling(length(all_top1_gene) / 5))
FeaturePlot(
  all_data,
  features = all_top1_gene,
  reduction = "harmonyumap",
  ncol = 5,
  raster = TRUE
) 
dev.off()


# spatial feature plot
DefaultAssay(all_data) <- "Spatial.008um"
Images(all_data) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  iwalk(function(img, img_name){
    
    png(file.path(img_output, str_glue("4 top1 gene spatial - {img_name}.png")), width = 400 * 5, height = 320 * ceiling(length(all_top1_gene) / 5))
    print(SpatialFeaturePlot(
      all_data,
      features = all_top1_gene,
      images = img,
      ncol = 5
    ))
    dev.off()
  })




# DoHeatmap(sc, features = all_de_genes, group.by = "RNA_snn_res.0.8")
if(!exists('avg_expr')) avg_expr <- AverageExpression(JoinLayers(all_data))[['sketch']]


all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")

# avg_expr[all_de_genes, ] %>% 
#   as.matrix() %>% 
#   as.data.frame() %>% 
#   setNames(., paste('cluster', (1:ncol(.) - 1) )) %>%  t %>%
#   pheatmap::pheatmap(
#     scale = "column",
#     cluster_row = F, # avg_expr[VariableFeatures(all_data), ] %>% t %>% dist %>% hclust,
#     cluster_col = FALSE,
#     color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#     border_col = NA
#   )

# BANKSY ----------------


# RCTD ------------------
library(spacexr)
rctd_out_dir <- 'RCTD' # dir.create(rctd_out_dir)

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

#2. single data
all_data_sc <- readRDS("all_data_sc_dim15_anno.rds")
# all_data_sc <- all_data


Idents(all_data_sc) <- "annotation"
counts_sc  <- all_data_sc %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation)
nUMI_sc    <- all_data_sc$nCount_RNA

sc_cell_filter_idx <- ! cluster_sc == "Unknown"

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
rctd_query_ls <-
  visium_roi_ls %>%
  map(function(visium){browser()
    counts_hd <- visium %>% LayerData(layer = "counts")
    coords    <- visium %>% GetTissueCoordinates() %>% .[, 1:2]
    
    # create the RCTD query object
    query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
    query
  })

#4. run RCTD
rctd_result_ls<- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

# saveRDS(rctd_result_ls, 'rctd_result_ls.rds')

#5. addMetaData
visium_roi_rctd_ls <-
  visium_roi_ls %>%
  imap(function(visium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    visium <- AddMetaData(visium, metadata = rctd@results$results_df)
    
    visium$first_type <- as.character(visium$first_type)
    visium$first_type[is.na(visium$first_type)] <- "Unknown"
    visium
  })

#6. spatial plot 
# SpatialFeaturePlot(
#   visium_roi_ls[[1]], 
#   image.scale = 'hires', 
#   pt.size.factor = 5,
#   shape = 22,
#   alpha = c(0.5, 0.8),
#   features = "nFeature_Spatial.008um"
# )

image_scale_ls <-
  visium_roi_rctd_ls %>%
  # visium_roi_ls %>% 
  map(GetTissueCoordinates) %>%
  map('x') %>%
  map(range) %>%
  map(diff) %>%
  map(log) %>%
  {map(., function(x, max_val)x/max_val, max_val = max(unlist(.)))}

pdf(
  file.path(rctd_out_dir, 'spatial plot of RCTD.pdf'), 
  height = 15, 
  width = 15
)
visium_roi_rctd_ls %>%
  imap(function(visium, name){
    SpatialDimPlot(
      visium,
      group.by       = 'first_type',
      image.scale    = 'hires',
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      alpha          = 0.5
    ) +
    ggtitle(name) +
    ggsci::scale_fill_igv() +
    guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()

# pahway score ------
geneset_ls <- list(
  Secretory_protein = c('SHBG', 'INHBB', 'FGF2', 'BMP4', 'KITLG'),
  Antioxidation_genes = c('SOD2', 'GPX4', 'JUND'),
  BTB_genes = c('CTNNB1', 'CDH1', 'TJP1', 'OCLN', 'CLDN1', 'LAMB2', 'GJA1', 'NCAM1')
)

# pathway score annotation -------
visium_pathway_ls <-
  visium_roi_ls %>%
  map(function(visium){ # browser()
    res_ls <-
      major_marker_ls %>% 
      enframe() %>% 
      with(split(value, name)) %>%
      map(function(pathway){ # browser()
        visium <- visium %>% AddModuleScore(features = list(pathway))
        visium$Cluster1
      })
    res_df <- as.data.frame(res_ls)
    rownames(res_df) <- colnames(visium)
    # visium <- visium %>% AddMetaData(res_df)
    # 
    # visium
    res_df
  })

# manual annotation ----
visium_annotation_ls <-
  visium_pathway_ls %>%
  map(function(x){ # browser()
    min_max <- function(x){
      if(all(near(x, x[1]))) return(x)
      (x-min(x))/(max(x)-min(x))
    }
    xx <- apply(x, 1, min_max)
    all_cells <- rownames(xx)
    sort2 <- function(x, ...) sort(x, decreasing = TRUE, ...)
    
    cell <- all_cells[apply(xx, 2, which.max)]
    score.1 <- apply(xx, 2, function(x) sort2(x)[1]/sort2(x)[2])
    score.2 <- apply(xx, 2, function(x) sum(sort2(x)[1:2])/sum(x))
    
    data.frame(
      annotation = cell,
      score.1 = score.1,
      score.2 = score.2,
      row.names = rownames(x)
    )
  })


# spatial plot 
pdf(
  file.path(rctd_out_dir, 'spatial plot of pathway score annotation.pdf'), 
  height = 15, 
  width = 15
)
visium_roi_rctd_ls %>%
  imap(function(visium, name){ # browser()
    annotation <- visium_annotation_ls[[name]]
    visium <- visium %>% AddMetaData(annotation)
    SpatialDimPlot(
      visium,
      group.by       = 'annotation',
      image.scale    = 'hires',
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      alpha          = 0.5
    ) +
      ggtitle(name) +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()


# spatial featureplot 

visium_roi_rctd_ls %>%
  imap(function(visium, name){ # browser()
    pathway_score <- visium_pathway_ls[[name]]
    visium <- visium %>% AddMetaData(pathway_score)
    
    pdf(
      file.path(rctd_out_dir, str_glue('spatial plot of pathway score - {name}.pdf')), 
      height = 15, 
      width = 15
    )
    
    pathway_score %>% 
      colnames() %>%
      map(function(x){
        SpatialFeaturePlot(
          visium,
          features = x,
          image.scale    = 'hires',
          pt.size.factor = 5 * image_scale_ls[[name]],
          shape          = 22,
          alpha          = c(0.1, 0.5)
        ) +
        ggtitle(name) 
      }) %>%
      walk(print)
    
    dev.off()
  })


# LSGI --------
library(RcppML)  # https://github.com/zdebruine/RcppML
library(LSGI)
library(Matrix)

lsgi_out_dir <- "LSGI"
dir.create(lsgi_out_dir)

scan.nmf.mse <- function(obj, ranks = seq(1, 30, 2), tol = 1e-04) {
    # users can customize the scan by changing 'ranks'
    # dat <- obj@assays$RNA@data
    dat <- LayerData(obj, layer = "data")
    errors <- c()
    ranks <- seq(1, 30, 2)
    for (i in ranks) {
        # cat('rank: ', i, '\n')
        mod <- RcppML::nmf(dat, i, tol = 1e-04, verbose = F)
        mse_i <- mse(dat, mod$w, mod$d, mod$h)
        errors <- c(errors, mse_i)
    }
    results <- data.frame(rank = ranks, MSE = errors)
    return(results)
}

sr.nmf <- function(obj, k = 10, tol = 1e-06, assay = "RNA") {
    # dat <- obj@assays$RNA@data
    dat <- LayerData(obj, layer = "data")
    nmf_model <- RcppML::nmf(dat, k = k, tol = tol, verbose = F)
    embeddings <- t(nmf_model$h)
    rownames(embeddings) <- colnames(obj)
    colnames(embeddings) <- paste0("nmf_", 1:k)
    loadings <- nmf_model$w
    rownames(loadings) <- rownames(obj)
    obj@reductions$nmf <- CreateDimReducObject(embeddings = embeddings,
        loadings = loadings, key = "nmf_", assay = assay)
    return(obj)
}

p_scan_nmf_ls <- visium_roi_ls %>%
  imap(function(visium, name){
    scan.nmf.res <- scan.nmf.mse(obj = visium)
    ggplot(scan.nmf.res, aes(x = rank, y = MSE)) + geom_point(size = 0.7) +
      geom_smooth(method = "loess", span = 0.2, color = "black",
                  linewidth = 1, se = F) + labs(x = "NMF rank", y = "MSE") +
      theme_classic() + scale_y_continuous(expand = c(0.01, 0)) +
      theme(aspect.ratio = 1) +
      ggtitle(name)
  })

# run nmf

# k_ls <- c(
#   AZFb   = 10,
#   AZFc1  = 10,
#   AZFc2  = 10,
#   N      = 10
# )

# visium_roi_lsgi_ls <-
#   visium_roi_ls[1] %>%
#   imap(function(visium, name){
#     k <- k_ls[name]
#     visium <- sr.nmf(obj = visium, k = k, tol = 1e-05)
#     visium
#   })

# run nmf on 4samples
# scan nmf
visium_roi <- 
  merge(visium_roi_ls[[1]], visium_roi_ls[-1]) %>%
  JoinLayers()

res.scan.nmf <- scan.nmf.mse(visium_roi)
saveRDS(res.scan.nmf, 'res.scan.nmf.4samples.rds')
plot(res.scan.nmf)


# run nmf
visium_roi_lsgi_ls <-
  visium_roi %>%
  sr.nmf(k = 10, tol = 1e-05) %>% # joint data used for nmf
  SplitObject(split.by = 'project')

# filter out some cells (the smallest nCount cells)
# make sure count of cells canbe divided by a number between 5-10
visium_roi_lsgi_ls <-
  visium_roi_lsgi_ls %>%
  map(function(visium){ # browser()
    spatial_coords <- GetTissueCoordinates(visium) %>% .[1:2]
    colnames(spatial_coords) <- c("X", "Y")
    
    embeddings <- visium@reductions$nmf@cell.embeddings
    
    # filter out some cells
    filter_out_number <- min(nrow(spatial_coords) %% 5:10)
    print(filter_out_number)
    
    if(filter_out_number == 0) return(visium)
    
    ncount_cols <- colnames(visium[[]]) %>% str_subset('nCount') %>% .[1]
    filter_out_order <- visium[[ncount_cols]][[1]] %>% order() %>% head(filter_out_number)
    
    filter_cells <- colnames(visium) %>% .[-filter_out_order]
    visium <- subset(visium, cells = filter_cells)
    
    visium
  })



# run lsgi
lsgi_res_ls <-
  visium_roi_lsgi_ls %>%
  imap(function(visium, name){ # browser()
    print(name)
    spatial_coords <- GetTissueCoordinates(visium) %>% .[1:2]
    colnames(spatial_coords) <- c("X", "Y")
    
    embeddings <- visium@reductions$nmf@cell.embeddings
    
    divided_number_ls <- DescTools::Divisors(nrow(spatial_coords)) %>% .[[1]] %>% .[-1]
    divided_number    <- divided_number_ls[which.min(abs(divided_number_ls - 8))]
    
    print(divided_number)
    # return()
    
    lsgi.res <-
      local.traj.preprocessing(
        spatial_coords    = spatial_coords,
        n.grids.scale     = divided_number ,
        embeddings        = embeddings, 
        n.cells.per.meta  = 25
    )
    
    lsgi.res
  })

saveRDS(lsgi_res_ls, 'lsgi_res_ls.4samples.rds')

# plot
cell_colors <- structure(
  ggsci::pal_igv()(10),
  names = str_c("nmf_", 1:10)
)

p_lsgi_res_ls <-
  lsgi_res_ls %>%
  imap(function(lsgi.res, name){
    # plot multiple factors
    plt.factors.gradient.ind(
      info = lsgi.res, 
      r_squared_thresh = 0.6,
      arrow.length.scale = 0.5, 
      minimum.fctr = 10) +  # plot gradient (had to appear in at least 10 grids)
    scale_color_manual(values = cell_colors, breaks = str_c("nmf_", 1:10))
  })

p_lsgi_res_he_ls <-
  lsgi_res_ls %>%
  imap(function(lsgi.res, name){
    visium <- visium_roi_lsgi_ls[[name]]
    
    p <- plot.overlay.factor.self(
      visium = visium,
      info   = lsgi.res, 
      img    = name,
      scale  =  "hires"
    ) 
    p$patches$plots[[1]] <-
      p$patches$plots[[1]] +
      scale_color_manual(values = cell_colors, breaks = str_c("nmf_", 1:10))
    p
  })


# export
pdf(
  file.path(lsgi_out_dir, "LSGI multi factor.pdf"),
  width  = 15,
  height = 15
)
p_lsgi_res_ls %>% imap(~ .x + ggtitle(.y)) %>% walk(print)
dev.off()

pdf(
  file.path(lsgi_out_dir, "LSGI multi factor with HE.pdf"),
  width  = 15,
  height = 15
)
p_lsgi_res_he_ls %>% imap(~ .x + ggtitle(.y))  %>% walk(print)
dev.off()



# gene functionnal enrichment 
get.nmf.info <- function(obj, top.n = 50) {
  feature.loadings <- as.data.frame(obj@reductions$nmf@feature.loadings)
  
  top.gene.list <- list()
  for (i in 1:ncol(feature.loadings)) {
    o <- order(feature.loadings[, i], decreasing = T)[1:top.n]
    features <- rownames(feature.loadings)[o]
    top.gene.list[[colnames(feature.loadings)[i]]] <- features
  }
  nmf.info <- list(feature.loadings = feature.loadings, top.genes = top.gene.list)
  return(nmf.info)
}


# nmf_info <- get.nmf.info(visium_roi_lsgi_ls$N)
nmf_info_ls <- visium_roi_lsgi_ls %>% map(get.nmf.info)


# export gene list
nmf_info_ls %>%
  map('top.genes') %>%
  map(as.data.frame) %>%
  .[[1]] %>% # 4个样本是一样的NMF结果
  openxlsx::write.xlsx(
    file = file.path(lsgi_out_dir, "NMF top 50 gene sorted by feature loading.xlsx"),
    asTable = TRUE
  )

nmf_info_ls$AZFb$feature.loadings %>% 
  rownames_to_column("gene") %>%
  openxlsx::write.xlsx(
    file = file.path(lsgi_out_dir, "NMF gene feature loading.xlsx"),
    asTable = TRUE
  )

# obtain gene sets
library(msigdbr)
library(hypeR)

mdb_h <- msigdbr(species = "Homo sapiens", category = "H")

gene.set.list <-
  mdb_h %>%
  dplyr::select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  nest() %>%
  deframe() %>%
  map(unlist, use.names = FALSE)

# run hypeR test
enrich_ls <-
  nmf_info_ls %>%
  map(function(nmf_info){
    mhyp <- hypeR(
      signature = nmf_info$top.genes,
      genesets  = gene.set.list,
      test = "hypergeometric", 
      background = rownames(nmf_info[["feature.loadings"]])
    )
    hyper.res.list <-
      mhyp$data %>%
      map(function(data){
        as.data.frame(data$data)
      })
    
    hyper.res.list
  })




# Visualize annotation results
p_enrich_ls_ls <-
  enrich_ls %>%
  map(function(enrich){
    p_enrich_ls <-
      enrich %>%
      imap(function(dat, name){
        dat[1:5, ] %>%
          ggplot(aes(x = reorder(label, -log10(fdr)), y = overlap/signature, fill = -log10(fdr))) + 
          geom_bar(stat = "identity", show.legend = T) + 
          xlab("Gene Set") + 
          ylab("Gene Ratio") +
          ggtitle(name) +
          viridis::scale_fill_viridis() + theme_classic() +
          coord_flip()  +
          theme(
            axis.text.x = element_text(color = "black", size = 12,angle = 45, hjust = 1), 
            axis.text.y = element_text(color = "black", size = 8, angle = 0),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
          )
      })
    
    p_enrich_ls
  })

# export
p_enrich_ls_ls %>%
  iwalk(function(p_ls, name){
    pdf(
      file.path(lsgi_out_dir, str_glue("LSGI enrichment analysis {name}.pdf")),
      width  = 8,
      height = 6
    )
    p_ls %>% walk(print)
    dev.off()
  })




# spatial plot --------------
p_sp_plot_anno <-
  Images(all_data) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    p <- SpatialDimPlot(
        all_data,
        images = img, 
        group.by = "annotation_full",
        pt.size.factor	= 4
      ) +
      ggsci::scale_fill_aaas() +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

p_sp_plot_cluster <-
  Images(all_data) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    p <- SpatialDimPlot(
      all_data,
      images = img, 
      group.by = "sketch_snn_res.0.4_full",
      pt.size.factor	= 4
    ) +
      # ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

# export for rmd
p_sp_plot_anno %>%
  iwalk(function(p, name){
    png(file.path(img_output, str_glue("3 spatial plot - annotation - {name}.png")), width = 720, height = 720)
    print(p)
    dev.off()
  })

p_sp_plot_cluster %>%
  iwalk(function(p, name){
    png(file.path(img_output, str_glue("3 spatial plot - cluster - {name}.png")), width = 720, height = 720)
    print(p)
    dev.off()
  })

# export
pdf(
  file.path(sp_output, "spatial plot annotation.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_anno %>% walk(print)
dev.off()

pdf(
  file.path(sp_output, "spatial plot cluster.pdf"),
  height = 8, 
  width  = 8
)
p_sp_plot_cluster %>% walk(print)
dev.off()

# enrichemt analysis -----
library(clusterProfiler)

all_de_genes_id_for_enrichmenet <-
  all_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  discard(function(x) length(x) < 1)


enrich_ls <-
  all_de_genes_id_for_enrichmenet %>%
  map(function(x){
    enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data = FALSE
    )
  })

# saveRDS(enrich_ls, 'enrich_ls.rds')

# expot enrichment table
enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    cluster_name <- paste0("cluster_", name)
    dir.create(file.path(er_output, cluster_name))
    
    write.csv(
      x@result,
      file = file.path(er_output, cluster_name, "enriment table.csv"),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  map(function(x){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    )
  })

# export dotplot
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    cluster_name <- paste0("cluster_", name)
    dir.create(file.path(er_output, cluster_name))
    
    ggsave(
      x,
      filename = file.path(er_output, cluster_name, "dotplot.pdf"),
      height = 6,
      width = 8
    )
  })

# export for rmd
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    cluster_name <- paste0("cluster_", name)
    
    png(file.path(img_output, str_glue("5 enrichment dotplot - {cluster_name}.png")), width = 720, height = 720)
    print(x)
    dev.off()
  })



# patch code -----------------------
dr_outdir2 <- str_glue("{dr_outdir} - 2") # dir.create(dr_outdir2)
all_files_dr_output <- list.files(dr_outdir, full.names = TRUE)

# dimplot
file.copy(
  all_files_dr_output %>% .[str_detect(basename(.), str_glue("Dimplot dim{dims}"))],
  dr_outdir2
)
# dotplot
file.copy(
  all_files_dr_output %>% .[str_detect(basename(.), str_glue("dotplot dims{dims}"))],
  dr_outdir2
)
# featureplot
file.copy(
  all_files_dr_output %>% .[str_detect(basename(.), str_glue("featureplot dims{dims}"))],
  dr_outdir2
)
# pca
file.copy(
  all_files_dr_output %>% .[str_detect(basename(.), "^pca")],
  dr_outdir2
)

