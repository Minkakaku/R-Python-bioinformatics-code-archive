# conda activate R_spaceranger4
# /mnt/NexData_file/home/liuran/projects/2_VisiumHD/FWSP20240360_20251106/analysis_spaceranger4
library(Seurat)
library(tidyverse)
library(patchwork)
library(sf)
library(ggplot2)
library(harmony)

# setting dir -------
img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)
mk_output <- "4. marker gene" # dir.create(mk_output)
er_output <- "5. enriment analysis" # dir.create(er_output)
data_output <- "data_for_analysis" # dir.create(data_output)

dir.create(img_output)
dir.create(qc_outdir)
dir.create(dr_outdir)
dir.create(sp_output)
dir.create(mk_output)
dir.create(er_output)
dir.create(data_output)


# set future -------------
library(future)
plan(multisession, workers = 4)
options(future.globals.maxSize = 200 * 1024 ^ 3) # 100Gb


# sample_df <- 
#   # read.csv("data/Annoation.csv")
#   data.frame(
#     sample = c("Z00", "Z22", "Z24", "Z88"),
#     group  = c("Z00", "Z22", "Z24", "Z88")
#   )

sample_df <-
  # read.csv("data/Annoation.csv")
  data.frame(
    sample = c("Z0_1","Z0_2","Z0_3","Z72_1",
               'Z2_1','Z2_2','Z2_3','Z72_2',
               'Z8_1','Z8_2','Z8_3','Z72_3',
               'Z24_1','Z24_2','Z24_3'),
    group  = c("Z00","Z00","Z00","Z00",
               "Z22","Z22","Z22","Z22",
               "Z88","Z88","Z88","Z88",
               "Z24","Z24","Z24")
  )


file_list <- 
  list.files(
    # "/home/szt/20240924_visiumHD_FWSP20240360/Spaceranger_out",
    "/mnt/NexData_file/home/hyj/Spatial_project/visium_hd/FWSP20253489-03_wcs_idney_20251230/analysis_spaceranger4/data", 
    pattern = paste(sample_df$group, collapse = "|"), 
    full.names = TRUE
  ) 
# %>%
#   str_subset("\\..*$", negate = TRUE) %>%
#   file.path(., basename(.))

stopifnot(nrow(file_list) == nrow(sample_df))

# TODO: fix bugs
# Samples order in file_list may not identical to sample_df
# names(file_list) <- sample_df$group # basename(file_list)
names(file_list) <- sample_df$group[match(basename(file_list), sample_df$group)]

# Load10X_Spatial

visim_ls <-
  file_list %>% 
  imap(function(path, name){
    path <- file.path(path, "outs")
    # print(name)
    # return()
    Load10X_Spatial(
      path, 
      bin.size = "polygons", 
      slice = name, 
      filter.matrix = TRUE,
      image.name = "tissue_hires_image.png"
    )
  }) %>% 
  setNames(names(file_list))

# visim_ls <-
  # visim_ls %>%
  # imap(function(visium, name){
    # visium$project <- name
    # Idents(visium) <- 'project'
    # RenameCells(visium, add.cell.id = name)
  # })

# parse ROI----------------
roi_ls <- list.files(
  'ROI',
  full.names = TRUE
) 

names(roi_ls) <- 
  roi_ls %>%
  basename() %>%
  str_remove("\\..*")

roi_coords_ls <-
  roi_ls %>%
  map(function(path){
    dat <- read.csv(path)
    colnames(dat) <- dat %>% colnames() %>% str_to_lower()
    dat<-dat %>% mutate(cell=barcode)%>%
      column_to_rownames('barcode')
    dat
  }) 

visium_roi_ls <-
  roi_coords_ls %>%
  imap(function(roi, name){ # browser()
    # get Visium Object
    visium_name <- sample_df$group[sample_df$sample == name]
    visium <- visim_ls[[visium_name]]
    
    # subset visium
    cells <- roi$cell
    visium_roi <- subset(visium, cells = cells)
    
    visium_roi
  })

# add sample_id or sample_name
visium_roi_ls <-
  visium_roi_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })


# visium_roi override visium
visim_ls <- visium_roi_ls

# visim_ls <- 
#   visim_ls %>%
#   map(~ {DefaultAssay(.x) <- 'Spatial.008um'; .x}) %>%
#   map(NormalizeData) 


# visim_ls <-
#   visim_ls %>%
#   imap(function(visium, name){
#     cat("Processing:", name, "\n")  # 添加进度提示
#     flush.console()
    
#     visium$project <- name
#     Idents(visium) <- 'project'
    
#     # 添加错误处理
#     tryCatch({
#       visium <- RenameCells(visium, add.cell.id = name)
#       cat("✓ Completed:", name, "\n")
#     }, error = function(e) {
#       cat("✗ Error in", name, ":", e$message, "\n")
#     })
    
#     return(visium)
#   })

visim_ls <- 
  visim_ls %>%
  map(~ {DefaultAssay(.x) <- 'Spatial.Polygons'; .x}) %>%
  map(NormalizeData) 


saveRDS(visim_ls, file = file.path(data_output, "visim_ls.rds"))

####qc -----
visim_ls <- 
  visim_ls %>%
  map(function(visium){
    subset(visium, nCount_Spatial.Polygons > 0)
  })

visim_ls$Z0_1@meta.data %>% head()

visim_ls <- 
  visim_ls %>%
  map(function(x){# browser()
    x[["percent.mt"]]    <- PercentageFeatureSet(x, pattern = "^mt-")
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = "^Rp[SL]")
    x$log10GenesPerUMI   <- log10(x$nFeature_Spatial.Polygons + 1) / log10(x$nCount_Spatial.Polygons + 1)
    x
  })

visim_ls$merge <-  merge(visim_ls[[1]], visim_ls[-1])

saveRDS(visim_ls, "./data_for_analysis/visim_ls.rds")
visim_ls <- readRDS("./data_for_analysis/visim_ls.rds")

p_qc_ls <-
  visim_ls %>%
  imap(function(x, name){ # browser()
    # only subset 10,000 cells
    x_sub <- subset(x, cells = Cells(x[['Spatial.Polygons']]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c("nFeature_Spatial.Polygons", "nCount_Spatial.Polygons", "percent.mt", "log10GenesPerUMI"), 
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
        # print(x)
        print(xx)
        
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
    x_sub <- subset(x, cells = Cells(x[['Spatial.Polygons']]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c("nFeature_Spatial.Polygons", "nCount_Spatial.Polygons", "percent.mt", "log10GenesPerUMI"), 
      ncol     = 4, 
      pt.size  = 0,
      log      = TRUE,
      raster   = TRUE
    ) +
      patchwork::plot_annotation(title = name)
    
  })

# export
saveRDS(p_qc_ls, file = file.path(data_output, "1. p_qc_ls.rds"))
p_qc_ls <- readRDS("./data_for_analysis/1. p_qc_ls.rds")

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
  "nFeature_Spatial.Polygons", "nCount_Spatial.Polygons",
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
    
    res <- qc_ls %>% purrr::discard(function(x) length(x) < 1) %>% str_c(collapse = " & ")
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

visim_ls_pass$merge<-NULL
visim_ls$merge<-NULL

# qc plot
qc_plot_filter_ls <-
  names(visim_ls) %>%
  setNames(., .) %>%
  map(function(sample){
    all_data <- visim_ls[[sample]] %>% subset(downsample = 100000, seed = 1)
    all_data_pass <- visim_ls_pass[[sample]] %>% subset(downsample = 100000, seed  = 1)
    
    # all_data
    plot1 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.Polygons", feature2 = "nFeature_Spatial.Polygons") +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.Polygons", feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = "nCount_Spatial.Polygons", feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + patchwork::plot_layout(ncol = 3)
    
    # all_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.Polygons", feature2 = "nFeature_Spatial.Polygons")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.Polygons", feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = "nCount_Spatial.Polygons", feature2 = "log10GenesPerUMI")+
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

all_data <- merge(visim_ls_pass[[1]], visim_ls_pass[-1])

## normalize for sketch data ------
all_data <- NormalizeData(all_data, scale.factor = floor(mean(all_data$nCount_Spatial.Polygons)))
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

sketch_cells_pass_qc <- all_data[['sketch']] %>% colnames() %>% .[. %in%  unlist(map(visim_ls_pass, colnames))]
all_data[['sketch']]  <- all_data[['sketch']] %>% subset(cells = sketch_cells_pass_qc)

### dimension reduction
DefaultAssay(all_data) <- "sketch"
all_data <- FindVariableFeatures(all_data, nfeatures = 2000)
all_data <- ScaleData(
  all_data,
  vars.to.regress = c(
    'nCount_Spatial.Polygons'
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
dims <- 15

all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE, assay = 'sketch')

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

# validate QC metric
pdf(file.path(dr_outdir, str_glue("Featureplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
FeaturePlot(
  all_data,
  features  = c('nCount_Spatial.Polygons', 'nFeature_Spatial.Polygons', 
                'percent.mt', 'log10GenesPerUMI'),
  reduction = "harmonyumap"
)
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
    algorithm = 1, 
    resolution = c(0.6, 0.8, 1.0)
  )

# dimplot
pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "sketch_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
#DimPlot(all_data, group.by = "sketch_snn_res.0.3", label = TRUE, reduction = "harmonyumap") + NoLegend()
#DimPlot(all_data, group.by = "sketch_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "sketch_snn_res.1", label = TRUE, reduction = "harmonyumap") + NoLegend()

# DimPlot(all_data, group.by = "SCT_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "SCT_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(all_data, "./data_for_analysis/all_data_clustered_dim15.rds")

## dotplot -------
major_marker_ls <- c(
  neron = 'Rbfox3', # 'NEUN',
  neron = 'Tubb3',
  neron = 'Map2',
  neron = 'Gap43',
  neron = 'Pcp4',

  'Glutamatergic neurons' = 'Slc17a7',
  'GABAergic neurons' = 'Gad2',
  'GABAergic neurons' = 'Gad1',
  'Dopaminergic neurons' = 'Th',
  'Dopaminergic neurons' = 'Aldh1a1',
  'Serotonergic neurons' = 'Tph1',
  'Serotonergic neurons' = 'Tph2',
  'Cholinergic neurons' = 'Chat',

  Astro = 'Gfap',
  Astro = 'Aldh1l1',
  Astro = 'S100b',

  Oligo = 'Mbp',
  Oligo = 'Plp1',
  Oligo = 'Mobp',
  Oligo = 'Cnp',

  'Microglial cells' = 'Dock8',
  'Microglial cells' = 'Sfmbt2',
  'Microglial cells' = 'Cx3cr1',
  'Microglial cells' = 'C1qc',

  'Endothelial cells' = 'Flt1',
  'Endothelial cells' = 'Vwf',
  'Endothelial cells' = 'Pecam1',
  'Endothelial cells' = 'Cdh5',
  'T cells' = 'Cd96',
  'T cells' = 'Skap1',
  'T cells' = 'Cd3g',
  'T cells' = 'Cd3e',
  'T cells' = 'Cd3d',
  Macrophages = 'Plxdc2',
  Macrophages = 'Cybb',
  Macrophages = 'Apbb1ip'
)


major_marker_ls1 <- c(
  MES = "Fhl2",
  MES = "Myl9",
  EC = "Pecam1",
  EC = "Vwf",
  EC = "Cldn5",
  EC = "Ackr1",
  POD = "Nphs2",
  POD = "Podxl",
  POD = "Ptpro",
  PT = "Cubn",
  PT = "Slc13a1",
  PT = "Lrp2",
  PT = "Aldob",
  DT = "Slc12a3",
  DT = "Calb1",
  DT = "Slc8a1",
  LOH = "Umod",
  LOH = "Slc12a1",
  LOH = "Cldn16",
  PC = "Aqp2",
  PC = "Aqp3",
  IC = "Slc26a7",
  IC = "Slc4a1",
  IC = "Atp6v1g3",
  MC = "Cd68",
  MC = "Mrc1",
  MC = "C1qa",
  MC = "C1qb",
  MON = "Vcan",
  MON = "Cd14",
  DC = "Fcer1a",
  DC = "Clec10a",
  FIB = "Lum",
  FIB = "Dcn",
  SMC = "Acta2",
  SMC = "Tagln",
  SMC = "Myh11",
  SMC = "Mylk",
  CC = "Top2a",
  CC = "Mki67"
  )

major_marker_ls2 <- c(
  EC = "Pecam1",
  EC = "Kdr",
  EC = "Nrp1",
  Podocyte = "Nphs2",
  Podocyte = "Cdkn1c",
  Podocyte = "Bcam",
  PT = "Cubn",
  PT = "Slc13a1",
  PT = "Lrp2",
  PT = "Aldob",
  DT = "Slc12a3",
  DT = "Calb1",
  DT = "Slc8a1",
  LOH = "Slc12a1",
  LOH = "Umod",
  LOH = "Cldn16",
  PC = "Aqp2",
  PC = "Aqp3",
  IC = "Slc26a7",
  IC = "Slc4a1",
  IC = "Atp6v1g3",
  FIB = "Lum",
  FIB = "Col1a1",
  FIB = "Dcn",
  FIB = "Plac8",
  FIB = "S100a4",
  FIB = "Pdgfrb",
  FIB = "Thy1",
  MC = "Lyz1",
  MC = "Cd68",
  MC = "Mrc1",
  MC = "C1qa",
  MC = "C1qb",
  Neutro = "S100a8",
  Neutro = "S100a9",
  BC = "Cd79a",
  BC = "Cd79b",
  TC = "Cd3d",
  TC = "Cd3e",
  TC = "Cd3g",
  NK = "Gzma",
  NK = "Nkg7",
  MES = "Fn1",
  MES = "Fhl2",
  MES = "Myl9"
)

major_marker_ls3 <- c(
  Endothelial = "Nrp1",
  Endothelial = "Kdr",
  Podocyte = "Nphs2",
  Proximal_tubule = "Lrp2",
  Proximal_tubule = "Slc27a2",
  ascendingloopof_Henle = "Umod",
  ascendingloopof_Henle = "Slc12a1",
  Distal_convoluted_tubule = "Slc12a3",
  Distal_convoluted_tubule = "Pvalb",
  collecting_duct_principal_cell = "Aqp2",
  collecting_duct_principal_cell = "Hsd11b2",
  CD_intercalated_cell = "Atp6v1g3",
  CD_intercalated_cell = "Atp6v0d2",
  CD_transitional_cell = "Insrr",
  CD_transitional_cell = "Rhbg",
  Fibroblast = "Plac8",
  Fibroblast = "S100a4",
  Fibroblast = "Mki67",
  Macrophage = "C1qa",
  Macrophage = "C1qb",
  Neutrophil = "S100a8",
  Neutrophil = "S100a9",
  B_cell = "Cd79a",
  B_cell = "Cd79b",
  T_cell = "Ltb",
  NK = "Gzma",
  NK = "Nkg7",
  Epi = "Epcam",
  Epi = "Krt19",
  LIP = "Fabp4",
  LIP = "Plin4",
  LIP = "Scd1"
)


theme_90angle <-
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
)

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} main_marker.pdf")), width = 25, height = 10)
#DotPlot(all_data, group.by = "sketch_snn_res.0.6", features = major_marker_ls1, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls1, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls2, dot.scale = 10) + theme_90angle
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = major_marker_ls3, dot.scale = 10) + theme_90angle

# DotPlot(all_data, group.by = "SCT_snn_res.0.3", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.4", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.6", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.8", dot.scale = 10, features = major_marker_ls) + theme_90angle
dev.off()

pdf(
  file.path(dr_outdir, str_glue("featureplot dims{dims} {reduction.use} main_marker.pdf")), 
  width = 50, 
  height = 50
)
FeaturePlot(
  all_data,
  features = major_marker_ls,
  ncol = 5,
  raster = TRUE
  )  
dev.off()


#find markers-------
Idents(all_data) <- "sketch_snn_res.0.8"
all_markers <-
  FindAllMarkers(
    JoinLayers(all_data),
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    only.pos = TRUE
  )

top10_markers <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>% 
  map(pull, 'gene') %>%
  imap(function(genes, name){
    re <- genes %>%
      map(function(gene){
        # browser()
        names(gene) = str_glue("cluster_{name}")
        gene
      })
    re
  }) %>%
  unlist() %>% 
  .[!duplicated(.)]

names(top10_markers) <- sub('.*\\.','',names(top10_markers))

top10_markers1 <- 
  all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 30) %>%
  arrange(desc(pct.1)) %>%
  slice_head(n = 10)

write.csv(top10_markers1,"top10_markers1.csv")

# top10_markers1 <- read.csv("top10_markers1.csv")

# top10_markers <-
#   top10_markers1%>%
#   nest() %>%
#   deframe()

# top10_markers <-
#   lapply(top10_markers, function(df) {
#     df %>% select(-X)  # -X 表示删除X列
#   })

# top10_markers <- top10_markers %>%
#   bind_rows() %>%  #合并列表为数据框
#   split(.$cluster) %>%  # 按cluster列拆分列表
#   map(pull, 'gene') %>%
#   imap(function(genes, name){
#     re <- genes %>%
#       map(function(gene){
#         # browser()
#         names(gene) = str_glue("cluster_{name}")
#         gene
#       })
#     re
#   }) %>%
#   unlist() %>% 
#   .[!duplicated(.)]

# names(top10_markers) <- sub('.*\\.','',names(top10_markers))

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} top10_marker.pdf")), width = 32, height = 10)
DotPlot(all_data, group.by = "sketch_snn_res.0.8", features = top10_markers, dot.scale = 10) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.3", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.4", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.6", dot.scale = 10, features = major_marker_ls) + theme_90angle
# DotPlot(all_data, group.by = "SCT_snn_res.0.8", dot.scale = 10, features = major_marker_ls) + theme_90angle
dev.off()

## annotation ----
cell_cluster <- c(
  '0'   = "LOH",
  '1'   = "PC",
  '2'   = "PT",
  '3'   = "PC",
  '4'   = "PT",
  '5'   = "DT",
  '6'   = "PT",
  '7'   = "LOH",
  '8'   = "FIB"  ,
  '9'   = "ENDO",
  '10'  = "POD",
  '11'  = "LOH",
  '12'  = "SMC",
  '13'  = "PT",
  '14'  = "PC",
  '15'  = "EPI",
  '16'  = "PT",
  '17'  = "NEU",
  '18'  = "MoMA",
  '19'  = "FIB",
  '20'  = "PT",
  '21'  = "delete",
  '22'  = "delete"
)

Idents(all_data) <- "sketch_snn_res.0.8"

all_data <- RenameIdents(all_data, cell_cluster)
all_data$annotation <- Idents(all_data)

pdf(file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation {reduction.use} .pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'umap') + NoLegend()
DimPlot(all_data, group.by = "annotation", label = TRUE, reduction = 'harmonyumap') + NoLegend()
dev.off()

# pdf(str_glue(str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2.pdf")), width = 25, height = 6)
# DotPlot(all_data, group.by = "annotation", features = major_cell_marker)
# dev.off()

## project data ------
all_data <- 
  all_data %>% 
  ProjectData(
    assay = "Spatial.Polygons",
    full.reduction = "harmony.full",
    
    sketched.assay = "sketch",
    sketched.reduction = "harmony",  # "harmony",
    
    umap.model = "harmonyumap", 
    dims = 1:dims, 
    
    refdata = list(
      annotation_full = "annotation",
      sketch_snn_res.0.8_full = "sketch_snn_res.0.8"
    )
  )

colnames(all_data@meta.data)
colnames(all_data@meta.data) <- gsub(
  pattern = "sketch_snn_res.0.6_full",  # 原列名（要替换的内容）
  replacement = "sketch_snn_res.0.8_full",  # 新列名（替换后的内容）
  x = colnames(all_data@meta.data)  # 待替换的列名向量
)

all_data$annotation_full %>% table()

##剔除注释为delete的细胞
all_data1 <- subset(all_data, annotation_full != c("delete"))
all_data1$annotation_full %>% table()
all_data <- all_data1
rm(all_data1)

Idents(all_data) <- 'annotation_full'
DefaultAssay(all_data) <- "Spatial.Polygons"
dims <- 15
reduction.use <- 'harmony'

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
  )
DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected group.pdf")),
  height = 6,
  width = 20
  )
DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap', split.by = "group") + NoLegend()
dev.off()

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
  )
DimPlot(all_data, group.by = "sketch_snn_res.0.8_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker projected.pdf")), width = 25, height = 6)
DotPlot(all_data, assay = "Spatial.Polygons", group.by = "annotation_full", features = top10_markers) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
DotPlot(all_data, assay = "Spatial.Polygons", group.by = "annotation_full", features = major_marker_ls1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
DotPlot(all_data, assay = "Spatial.Polygons", group.by = "annotation_full", features = major_marker_ls2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
DotPlot(all_data, assay = "Spatial.Polygons", group.by = "annotation_full", features = major_marker_ls3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()


saveRDS(all_data, "data_for_analysis/all_data_projected.rds")
#all_data <- readRDS("data_for_analysis/all_data_projected.rds")


# # ROI分析 ---------------------
# roi_df <-
#   list.files('ROI_2', full.names = TRUE, pattern = "csv") %>%
#   # setNames(., str_remove_all(basename(.), "^.*?-|-.*?$")) %>%
#   setNames(., str_extract(., "(?<=ROI_2/)[^.]+")) %>%
#   map(~ read.csv(.x, col.names = c('Barcode', 'ROI'))) %>%
#   enframe("group") %>%
#   unnest(value) %>%
#   mutate(Barcode = str_c(group, Barcode, sep = "_")) %>%
#   # mutate(roi_group = str_c(group, ROI, sep = "_")) %>%
#   dplyr::select(Barcode, ROI) %>%
#   column_to_rownames("Barcode")


# all_data <- 
#   all_data %>% 
#   AddMetaData(roi_df)

# all_data$ROI[is.na(all_data$ROI)] <- "Other"


# # ROI_3分析 ---------------------
# roi_df <- 
#   list.files('ROI_3', full.names = TRUE, pattern = "csv") %>%
#   # setNames(., str_remove_all(basename(.), "^.*?-|-.*?$")) %>%
#   setNames(., str_extract(., "(?<=ROI_3/)[^.]+")) %>%
#   map(~ read.csv(.x, col.names = c('Barcode', 'ROI'))) %>% 
#   enframe("group") %>%
#   unnest(value) %>%
#   mutate(Barcode = str_c(group, Barcode, sep = "_")) %>%
#   # mutate(roi_group = str_c(group, ROI, sep = "_")) %>%
#   dplyr::select(Barcode, ROI) %>%
#   column_to_rownames("Barcode")


# all_data <- 
#   all_data %>% 
#   AddMetaData(roi_df)

# all_data$ROI[is.na(all_data$ROI)] <- "Other"
# table(all_data$ROI)

# #saveRDS(all_data, "data_for_analysis/all_data_roi_1230.rds")



## annotation 空间映射图-----
#images改名
name_mapping <- c(
  "Z00.polygons"       = "Z0_1.polygons",
  "Z00.polygons.2"     = "Z0_2.polygons",
  "Z00.polygons.3"     = "Z0_3.polygons",
  "Z22.polygons"       = "Z2_1.polygons",
  "Z22.polygons.5"     = "Z2_2.polygons",
  "Z22.polygons.6"     = "Z2_3.polygons",
  "Z24.polygons"       = "Z24_1.polygons",
  "Z24.polygons.8"     = "Z24_2.polygons",
  "Z24.polygons.9"     = "Z24_3.polygons",
  "Z00.polygons.10"    = "Z72_1.polygons",
  "Z22.polygons.11"    = "Z72_2.polygons",
  "Z88.polygons"       = "Z72_3.polygons",
  "Z88.polygons.13"    = "Z8_1.polygons",
  "Z88.polygons.14"    = "Z8_2.polygons",
  "Z88.polygons.15"    = "Z8_3.polygons"
)

names(all_data@images) <- name_mapping[names(all_data@images)]

for (i in Images(all_data)){
  p1 <- SpatialDimPlot(all_data, 
  group.by = "annotation_full", images = i, image.scale = "hires", plot_segmentations = TRUE, image.alpha = 0.2, stroke = 0, pt.size.factor = 0.5)
  p2 <- p1 + scale_fill_manual(values = plan) + theme_bw() + xlab("")+ ylab("") + theme(axis.text = element_blank(), axis.ticks = element_blank())
  
  ggsave(p2, filename = paste0("3. spatial plot/", str_remove(i, "\\.polygons$"), "_annotation_spatial.pdf"), width = 20, height = 20)
  ggsave(p2, filename = paste0("3. spatial plot/", str_remove(i, "\\.polygons$"), "_annotation_spatial.png"), width = 20, height = 20, dpi = 300)
}
###加入group-------
all_data@meta.data$sample <- all_data@meta.data$project
group_map <- c(
  "Z0_1" = "Z0",
  "Z0_2" = "Z0",
  "Z0_3" = "Z0",
  "Z72_1" = "Z72",
  "Z2_1" = "Z2",
  "Z2_2" = "Z2",
  "Z2_3" = "Z2",
  "Z72_2" = "Z72",
  "Z8_1" = "Z8",
  "Z8_2" = "Z8",
  "Z8_3" = "Z8",
  "Z72_3" = "Z72",
  "Z24_1" = "Z24",
  "Z24_2" = "Z24",
  "Z24_3" = "Z24"
)
sample_col <- "sample" 
samples <- as.character(all_data@meta.data[[sample_col]])
group_vector <- group_map[samples]
names(group_vector) <- colnames(all_data)  # 确保分组向量与细胞barcode一一对应

all_data <- AddMetaData(
  object = all_data,
  metadata = group_vector,
  col.name = "group"
)

##细胞比例统计图-----
all_meta <- all_data@meta.data
all_meta$annotation_full %>% table()

celltype_counts <- all_meta %>%
  mutate(
    annotation_full = factor(annotation_full)
  ) %>%
  group_by(annotation_full, sample) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(sample) %>%
  mutate(prop = n / sum(n))

p <- ggplot(celltype_counts, aes(x = .data[['sample']], y = prop, fill = .data[['annotation_full']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = plan) +
  # 坐标轴与图例优化
  labs(x = "sample", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

ggsave(
  filename = paste0("2. dimension reduction and cluster/","sample_Celltype_proportions.pdf"),  # 保存文件名
  plot = p,                                
  width = 12,                              
  height = 5,                             
  dpi = 300,                               
  bg = "white"                             
)

celltype_counts <- all_meta %>%
  mutate(
    annotation_full = factor(annotation_full)
  ) %>%
  group_by(annotation_full, group) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(group) %>%
  mutate(prop = n / sum(n))

p <- ggplot(celltype_counts, aes(x = .data[['group']], y = prop, fill = .data[['annotation_full']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = plan) +
  # 坐标轴与图例优化
  labs(x = "group", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

ggsave(
  filename = paste0("2. dimension reduction and cluster/","group_Celltype_proportions.pdf"),  # 保存文件名
  plot = p,                                
  width = 5,                              
  height = 5,                             
  dpi = 300,                               
  bg = "white"                             
)

# 输出loupe browsers文件------
slide_mapping <- c(
  # Z00对应的sample
  "Z0_1" = "Z00",
  "Z0_2" = "Z00",
  "Z0_3" = "Z00",
  "Z72_1" = "Z00",
  # Z2对应的sample
  "Z2_1" = "Z22",
  "Z2_2" = "Z22",
  "Z2_3" = "Z22",
  "Z72_2" = "Z22",
  # Z8对应的sample
  "Z8_1" = "Z88",
  "Z8_2" = "Z88",
  "Z8_3" = "Z88",
  "Z72_3" = "Z88",
  # Z24对应的sample
  "Z24_1" = "Z24",
  "Z24_2" = "Z24",
  "Z24_3" = "Z24"
)

# 在meta.data中添加slide列
all_data@meta.data$slide <- slide_mapping[all_data@meta.data$sample]
table(all_data@meta.data$slide, all_data@meta.data$sample)

meta.data <- all_data@meta.data %>% 
  rownames_to_column("Barcode")

prefixes_to_remove <- c("Z0_1_", "Z0_2_", "Z0_3_", "Z72_1_", 
                        "Z2_1_", "Z2_2_", "Z2_3_", "Z72_2_",
                        "Z8_1_", "Z8_2_", "Z8_3_", "Z72_3_",
                        "Z24_1_", "Z24_2_", "Z24_3_")

# 拼接为正则表达式（|表示“或”，\\Q\\E避免特殊字符（如_）被解析为正则）
prefix_pattern <- paste0("\\Q", paste(prefixes_to_remove, collapse = "\\E|\\Q"), "\\E")

# 批量删除所有前缀
meta.data$Barcode <- gsub(
  pattern = prefix_pattern, 
  replacement = "", 
  x = meta.data$Barcode,
  fixed = FALSE  # 启用正则匹配，实现多前缀批量替换
)

head(meta.data$Barcode, 10)

for (s in unique(meta.data$slide)){
  anno_data <- meta.data%>% 
  filter(slide == s) %>% 
  #rownames_to_column("Barcode") %>% 
  #mutate(Barcode = str_remove(Barcode, paste0(s, "_"))) %>% 
  dplyr::select(Barcode, annotation_full) %>% 
  dplyr::rename(cell_type = annotation_full)

    # 打印去重前后行数
  cat("Slide", s, "去重前行数：", nrow(anno_data), "\n")
  
  anno_data <- anno_data %>% 
    dplyr::distinct(Barcode, .keep_all = TRUE)
  
  cat("Slide", s, "去重后行数：", nrow(anno_data), "\n\n")
  
  write.csv(anno_data, file = paste0("annocsv/", s, "_annotation.csv"), row.names = FALSE)
}





plan <- c("#EA7070","#FDC4B6","#E59572","#2694AB","#96CEB4","#ED5485",
        "#D9534F","#FFAD60","#57D1C9","#FFEEAD","#7F95D1","#FFE869",
        "#DE4307","#F29C2B","#F6D04D","#8BC24C","#4695D6","#FED95C",
        "#FA6E57","#F69E53","#248888","#999999","#E7475E","#F0D879",
        "#29A2C6","#FFCB18","#73B66B","#FF6D31","#F1AC9D","#F06966",
        "#DEE2D1","#6ABE83","#FCF4D9","#FFB85F","#FF7A5A","#8ED2C9",
        "#A5DFF9","#EF5285","#60C5BA","#FEEE7D","#E8A0B8","#FFC300",
        "#BCCF3D","#02C9C9","#FBF579","#005995","#FA625F","#600473",
        "#F17D80","#737495")

#把RCTD的first_type加入all_data@meta.data中

all_data_rctd <- readRDS("data_for_analysis/all_data_rctd.rds")

# 确保 all_data 的 meta.data 包含 cell_id 列（你已添加）
all_data@meta.data$cell_id <- rownames(all_data@meta.data)  # 更稳健的做法：meta.data 行名 = cell_id
#all_data$first_type <- "Other"

all_data_rctd@meta.data$cell_id <- rownames(all_data_rctd@meta.data)
first_type_mapping <- unique(all_data_rctd@meta.data[, c("cell_id", "first_type")])

# 在 all_data@meta.data 中按 cell_id 匹配并添加 first_type
all_data@meta.data$first_type <- first_type_mapping$first_type[
  match(rownames(all_data@meta.data), first_type_mapping$cell_id)
]

# 检查结果
table(all_data$first_type)

saveRDS(all_data, "data_for_analysis/all_data_merge_rctd.rds")


# RCTD结果绘制空间图 --------

all_data$first_type[all_data$first_type == "Other"] <- NA

for (i in Images(all_data)){
  p1 <- SpatialDimPlot(all_data, 
  group.by = "first_type", images = i, image.scale = "hires", plot_segmentations = TRUE, image.alpha = 0.2, stroke = 0, pt.size.factor = 0.5)
  p2 <- p1 + scale_fill_manual(values = plan, na.value = NA, labels = function(x) str_replace_na(x, "")) + 
  guides(fill = guide_legend(override.aes = list(size = 5, color = NA))) + 
  theme_bw() + xlab("")+ ylab("") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank())
  
  ggsave(p2, filename = paste0("3. spatial plot/", str_remove(i, "\\.polygons$"), "_RCTD_spatial.pdf"), width = 20, height = 20)
  ggsave(p2, filename = paste0("3. spatial plot/", str_remove(i, "\\.polygons$"), "_RCTD_spatial.png"), width = 20, height = 20, dpi = 300)
}


# RCTD结果输出loupe browsers文件------
meta.data <- all_data@meta.data %>% 
  rownames_to_column("Barcode")

prefixes_to_remove <- c("Z0_1_", "Z0_2_", "Z0_3_", "Z72_1_", 
                        "Z2_1_", "Z2_2_", "Z2_3_", "Z72_2_",
                        "Z8_1_", "Z8_2_", "Z8_3_", "Z72_3_",
                        "Z24_1_", "Z24_2_", "Z24_3_")

# 拼接为正则表达式（|表示“或”，\\Q\\E避免特殊字符（如_）被解析为正则）
prefix_pattern <- paste0("\\Q", paste(prefixes_to_remove, collapse = "\\E|\\Q"), "\\E")

# 批量删除所有前缀
meta.data$Barcode <- gsub(
  pattern = prefix_pattern, 
  replacement = "", 
  x = meta.data$Barcode,
  fixed = FALSE  # 启用正则匹配，实现多前缀批量替换
)

head(meta.data$Barcode, 10)

library(dplyr)

for (s in unique(meta.data$slide)){
  anno_data <- meta.data%>% 
  filter(slide == s) %>% 
  #rownames_to_column("Barcode") %>% 
  #mutate(Barcode = str_remove(Barcode, paste0(s, "_"))) %>% 
  dplyr::select(Barcode, first_type) %>% 
  dplyr::rename(cell_type = first_type)
  
  # 打印去重前后行数
  cat("Slide", s, "去重前行数：", nrow(anno_data), "\n")
  
  anno_data <- anno_data %>% 
    dplyr::distinct(Barcode, .keep_all = TRUE)
  
  cat("Slide", s, "去重后行数：", nrow(anno_data), "\n\n")

  write.csv(anno_data, file = paste0("annocsv/", s, "_rctd_annotation.csv"), row.names = FALSE)
}


##RCTD细胞比例统计图-----
all_meta <- all_data@meta.data
#all_meta$first_type[is.na(all_meta$first_type)] <- "Other"
all_meta$first_type %>% table()

celltype_counts <- all_meta %>%
  mutate(
    first_type = factor(first_type)
  ) %>%
  group_by(first_type, sample) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(sample) %>%
  mutate(prop = n / sum(n))

p <- ggplot(celltype_counts, aes(x = .data[['sample']], y = prop, fill = .data[['first_type']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = plan) +
  # 坐标轴与图例优化
  labs(x = "sample", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

ggsave(
  filename = paste0("2. dimension reduction and cluster/","sample_RCTD_Celltype_proportions.pdf"),  # 保存文件名
  plot = p,                                
  width = 12,                              
  height = 5,                             
  dpi = 300,                               
  bg = "white"                             
)

celltype_counts <- all_meta %>%
  mutate(
    first_type = factor(first_type)
  ) %>%
  group_by(first_type, group) %>%
  tally() %>%
  # 计算比例（按样本分组）
  group_by(group) %>%
  mutate(prop = n / sum(n))

p <- ggplot(celltype_counts, aes(x = .data[['group']], y = prop, fill = .data[['first_type']])) +
  geom_bar(stat = "identity", position = "stack") +
  # 应用自定义颜色方案
  scale_fill_manual(values = plan) +
  # 坐标轴与图例优化
  labs(x = "group", y = "Cell Type Proportion", fill = "Cell Type") +
  theme_minimal() +
  # 若样本名多，旋转 x 轴标签避免重叠
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

ggsave(
  filename = paste0("2. dimension reduction and cluster/","group_RCTD_Celltype_proportions.pdf"),  # 保存文件名
  plot = p,                                
  width = 5,                              
  height = 5,                             
  dpi = 300,                               
  bg = "white"                             
)




# Slc6a11 Astro表达量
astro_slc6a11 <- FetchData(
  object = all_data,
  vars = c("Slc6a11", "first_type"),
  cells = WhichCells(all_data, expression = first_type == "Astro")
)

astro_slc6a11$Cell_ID <- rownames(astro_slc6a11)
astro_slc6a11 <- astro_slc6a11[, c("Cell_ID", "Slc6a11", "first_type")]

write.csv(astro_slc6a11, "Astro_cells_Slc6a11_expression.csv", row.names = FALSE)


# 根据Astro中Slc6a11表达量分为两群；
Astro_anno <- read.table(file = "Astro_cells_Slc6a11_annotation.csv", header = TRUE, sep = ",")
head(Astro_anno)
nrow(Astro_anno)
colnames(Astro_anno) <- c("cell_id", "annotation_sub")

all_data@meta.data$Cells <- all_data@meta.data$cell_id
head(all_data@meta.data)

all_data@meta.data <- all_data@meta.data %>%
  left_join(Astro_anno, by = "cell_id") %>%
  mutate(
    annotation_updated = case_when(
      # 有 Astro_anno 注释的细胞
      !is.na(annotation_sub) ~ annotation_sub,
      
      # 默认情况：未分类的 Astro 细胞
      TRUE ~ first_type
    )
  ) 

unique(all_data$annotation_updated)
all_data@meta.data$annotation_sub <- all_data@meta.data$annotation_updated
head(all_data@meta.data)

all_data@meta.data <- 
  all_data@meta.data %>% 
  select(-c(annotation_updated, Cells))

Idents(all_data) <- "annotation_sub"
rownames(all_data@meta.data) <- all_data@meta.data$cell_id
table(all_data@meta.data$annotation_sub)

saveRDS(all_data, file = "data_for_analysis/all_data_sub.rds")
