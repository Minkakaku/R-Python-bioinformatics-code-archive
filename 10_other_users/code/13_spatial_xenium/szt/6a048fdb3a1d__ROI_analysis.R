library(Seurat)
library(tidyverse)

# setting dir -------
img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot" # dir.create(sp_output)
mk_output <- "4. marker gene" # dir.create(mk_output)
er_output <- "5. Calb1 gene" # dir.create(er_output)


# set default ------------
orginism <- "mm" # "hsa" # "mm"
orginism_kegg <- "mmu" # "hsa" # "mmu"

mt_pattern   <- if(orginism == 'hsa') "^MT-" else "^mt-"
ribo_pattern <- if(orginism == 'hsa') "^RP[SL]" else "^Rp[sl]"

org.db <- if(orginism == 'hsa') "org.Hs.eg.db" else "org.Mm.eg.db"

bin_size <- 8
assay_name <- case_when(
  bin_size == 8  ~ 'Spatial.008um',
  bin_size == 16 ~ 'Spatial.016um',
  TRUE ~ 'wrong bin size'
)

nfeature_name <- str_glue("nFeature_{assay_name}")
ncount_name   <- str_glue("nCount_{assay_name}")

# 4 11 20 hvg5000
# 12 hvg2000 regression out H cell


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
sample_df <-  #read.csv("data/Annoation.csv")
  data.frame(
    sample = c('DZZ_Saline_34_0um','DZZ_Saline_34_100um','SYZ_NH4_44_0um','SYZ_NH4_44_100um'),
    group  = c('DZZ_Saline','DZZ_Saline','SYZ_NH4','SYZ_NH4')
  )

file_list <- 
  list.files(
    # "/home/ztt/scRNA/spaceranger_out/FWSP20240508/spaceranger_output",
    "/data/outputdata/outputdata-sp/FWSP20240508/spaceranger_output",
    pattern = paste(sample_df$sample, collapse = "|"), 
    full.names = TRUE
  ) 

names(file_list) <- sample_df$group # basename(file_list)


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
# parse ROI----------------
roi_ls <- list.files(
  'data/ROI_coordinates',
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


# add sample_id or sample_name -------------------
visium_roi_ls <-
  visium_roi_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })

# visium_roi override visium
visim_ls <- visium_roi_ls


visim_ls <- 
  visim_ls %>%
  map(~ {DefaultAssay(.x) <- 'Spatial.008um'; .x}) %>%
  map(NormalizeData) 

# qc -----
visim_ls <- 
  visim_ls %>%
  map(function(visium){
    if(assay_name == "Spatial.008um"){
      subset(visium, nCount_Spatial.008um > 0)
    }else{
      subset(visium, nCount_Spatial.016um > 0)
    }
    
  })

visim_ls <- 
  visim_ls %>%
  map(function(x){# browser()
    x[["percent.mt"]]    <- PercentageFeatureSet(x, pattern = mt_pattern)
    x[["percent.ribo"]]  <- PercentageFeatureSet(x, pattern = ribo_pattern)
    x$log10GenesPerUMI   <- log10(x[[nfeature_name]][[1]] + 1) / log10(x[[ncount_name]][[1]] + 1)
    x
  })


p_qc_ls <-
  visim_ls %>%
  imap(function(x, name){ # browser()
    # only subset 10,000 cells
    x_sub <- subset(x, cells = Cells(x[[assay_name]]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c(nfeature_name, ncount_name, "percent.mt", "log10GenesPerUMI"), 
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
    x_sub <- subset(x, cells = Cells(x[[assay_name]]), downsample = 100000, seed = 1)
    p <- VlnPlot(
      x_sub, 
      features = c(nfeature_name, ncount_name, "percent.mt", "log10GenesPerUMI"), 
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
qc_metric <- read.csv("qc_metric_008um.csv") %>% map_dfc(str_trim)

qc_metric_variable <- c(
  nfeature_name, ncount_name,
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
    plot1 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = nfeature_name) +
      patchwork::plot_annotation(title = sample)
    plot2 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = "percent.mt") +
      patchwork::plot_annotation(title = sample)
    plot3 <- FeatureScatter(all_data, feature1 = ncount_name, feature2 = "log10GenesPerUMI") +
      patchwork::plot_annotation(title = sample)
    p1 <- plot1 + plot2 + plot3 + patchwork::plot_layout(ncol = 3)
    
    # all_data_pass
    plot1_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = nfeature_name)+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot2_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = "percent.mt")+
      patchwork::plot_annotation(title = str_glue('{sample} + log'))
    plot3_pass <- FeatureScatter(all_data_pass, feature1 = ncount_name, feature2 = "log10GenesPerUMI")+
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
all_data <- NormalizeData(all_data,scale.factor=median(all_data$nCount_Spatial.008um),verbose=FALSE)
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
    ncount_name
  )
)

# all_data <- SCTransform(all_data, assay = "sketch")

all_data <- RunPCA(all_data, features = VariableFeatures(all_data))

pdf(file = file.path(dr_outdir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 50)
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
dims <- 10

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
major_marker_ls%>%str_to_title()->major_marker_ls
major_marker_ls <- c(
  # neron = 'Rbfox3',
  # neron = 'Tubb3',
  # neron = 'Map2',
  # neron = 'Gap43',
  'inhibitory'='Sst',
  'inhibitory'='Pvalb',
  'inhibitory'='Lamp5',
  'inhibitory'='Vip',
  'excitatory'='SLC17A6',
  'Glutamatergic neurons' = 'Slc17a7',
  'GABAergic neurons' = 'Gad2',
  'GABAergic neurons' = 'Gad1',
  'MSN'='Syt6',
  'MSN'='Otof',
  'MSN'='Drd1',
  'MSN'='Adora2a',
  
  'Dopaminergic neurons' = 'Th',
  'Dopaminergic neurons' = 'Aldh1a1',
  'Serotonergic neurons' = 'Tph1',
  'Serotonergic neurons' = 'Tph2',
  'Cholinergic neurons' = 'Chat',


  Astro = 'Aqp4',
  Astro = 'Rorb',
  Astro = 'Slc7a10',
  Asrto = 'Gja1',
  
  OPC = 'Pdgfra',
  OEG = 'Npy',
  mNEUR = 'Syt1',
  Nendc = 'Baiap3',
  # HypEPC = 'Ccdc153',
  # HypEPC = 'Sspo',
  # HypEPC = 'Rax',
  # HypEPC = 'Ttr',
  
  
  
  Oligo = 'Mog',
  Oligo = 'Sox10',
  Oligo = 'Olig1',
  Oligo = 'Cldn11',
  Oligo = 'Mbp1',
  
  'Microglial cells' = 'Dock8',
  'Microglial cells' = 'Sfmbt2',
  'Microglial cells' = 'Cx3cr1',
  'Microglial cells' = 'C1qc',
  'Microglial cells' = 'P2ry12',
  'Microglial cells' = 'Aif1',
  'Microglial cells' = 'Tmem119',
  'Microglial cells' = 'Ctss',
  
  'Epen'='Foxj1',
  'Endothelial cells' =	'Flt1',
  'Endothelial cells' =	'Vwf',
  'Endothelial cells' =	'Pecam1',
  'Endothelial cells' =	'Cldn5',
  'Endothelial cells' =	'Cdh5',
  
  "Pericytes" = "Pdgfrb",
  "Pericytes" = "Rgs5",
  "Pericytes" = "Mcam",
  "Pericytes" = "Kcnj8",
  
  "Smooth muscle cells" = "Pln",
  "Smooth muscle cells" = "Mylk",
  "Smooth muscle cells" = "Myh11",
  "Smooth muscle cells" = "Acta2",
  
  "Hb_VC" = "Alas2",
  "VLMC" = "Slc6a13",
  "ABC" = "Slc47a1",
  'MNC' = "Plac8",
  'Mac' = "Pf4",
  'DC' = 'Cd209a',
  
  'Neutrophil' = "S100a8",
  'Neutrophil' = "S100a9",
  'Neutrophil' = "Csf3r",
  'Neutrophil' = "Fcgr3b",
  
  'T cells' =	'Cd3e',
  Macrophages	= 'Plxdc2',
  Macrophages	= 'Cybb',
  Macrophages	= 'Apbb1ip'
  
)


major_marker_ls <- c('Gapgh','Pdgfra','Cldn11','Npy','Thbs4','Cd44','Gja1','Cdk1','Sox11',
                     'Syt1','Baiap3','Ccdc153','Sspo','Rax','Ttr','Clan5','Kcnj8','Acta2',
                     'Alas2','Slc6a13','Slc47a1','Tmem119','Plac8','Pf4','Cd209a','S100a9')

major_marker_ls <- c('Cux2','Rspo1','Scube1','Fezf2','Ndst4','Nxph4','Hs3st4','Tshz2','Chat','Ptpru',
                     'Sst','Pyalb','Syt6','Cpne7','Lamp5','Lhx6','Vip','Adarb2','Calb2','Otof','Drd1',
                     'Adora2a','Pdgfra','Olig1','Rorb','Aqp4','Foxj1','Cspg4','Vtn','Cldn5','F13a1',
                     'Cd3e','Ctss')

major_marker_ls <- c(
  'excitatory'='Slc17a7',
  'excitatory'='SLC17A6',
  'inhibitory'='Gad2',
  'inhibitory'='Gad1',
  'inhibitory'='Lhfpl3',
  
  Astro = 'Aqp4',
  Astro = 'Sparcl1',
  
  Oligo = 'Mog',
  Oligo = 'Mag',
  Oligo = 'Apod',
  Oligo = 'Olig1',
  Oligo = 'Mbp1',
  # Oligo = 'Mobp',
  # Oligo = 'Cnp',
  
  'Microglial cells' = 'Ctss',
  'Microglial cells' = 'C1qb',
  
  'Vascular cells' =	'Ltih5',
  
  'Oligodendrocyte Precursors' = 'Pdgfra'
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
Idents(all_data) <- "sketch_snn_res.0.6"

cell_cluster <- c(
  '0'   = "InN",
  '1'   = "ExN",
  '2'   = "Oligo",
  '3'   = "Astro",
  '4'   = "mNEUR",
  '5'   = "Endo",
  '6'   = "Oligo",
  '7'   = "Th+ NendC",
  '8'   = "Micro",
  '9'   = "ExN_InN",
  '10'  = "Astro",
  '11'  = "Cholinergic neurons",
  '12'  = "SMCs",
  '13'  = "NendC",
  '14'  = "Epen"
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
full_umap <- ifelse(reduction.use == "pca", "umap", "harmonyumap")
all_data <- 
  all_data %>% 
  ProjectData(
    assay = "Spatial.008um",
    full.reduction = "pca.full",
    
    sketched.assay = "sketch",
    sketched.reduction = reduction.use, # "harmony",
    
    umap.model = full_umap, 
    dims = 1:dims, 
    
    refdata = list(
      annotation_full = "annotation",
      sketch_snn_res.0.6_full = "sketch_snn_res.0.6"
    )
  )

Idents(all_data) <- 'annotation_full'
DefaultAssay(all_data) <- "Spatial.008um"

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} annotation - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "annotation_full", label = TRUE, reduction = str_glue('full.{full_umap}')) + NoLegend()
dev.off()

pdf(
  file.path(dr_outdir, str_glue("Dimplot dim{dims} cluster - {reduction.use} - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data, group.by = "sketch_snn_res.0.6_full", label = TRUE, reduction = str_glue('full.{full_umap}')) + NoLegend()
dev.off()

# saveRDS(all_data, "all_data_projected_1028.rds")
# all_data <- readRDS("all_data_projected.rds")

pdf(file.path(dr_outdir, str_glue("dotplot dims{dims} {reduction.use} annotation major_marker2 projected.pdf")), width = 25, height = 6)
DotPlot(all_data, assay = "sketch", group.by = "annotation_full", features = major_marker_ls) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()

# export for rmd
png(file.path(img_output, str_glue("2 dimplot - cluster.png")), width = 600, height = 480)
DimPlot(
  all_data, 
  group.by = "sketch_snn_res.0.6_full", 
  label = TRUE, 
  reduction = full_umap
) +
  NoLegend()
dev.off()

#prop barplot-------------
coords_enrichment <-
  all_data %>%
  {left_join(
    .[[c('annotation_full', "project")]] %>% rownames_to_column("cell"),
    GetTissueCoordinates(.),
    by = 'cell'
  )}
props_df <-
  coords_enrichment %>%
  group_by(project, annotation_full) %>%
  summarise(count = n()) %>%
  group_by(project) %>%
  mutate(freq = count/sum(count) * 100) 
p_prop_celltype <-
  props_df %>%
  
  ggplot(aes(x = project, y = freq, fill = annotation_full)) +
  geom_col() +
  ggpubr::theme_pubclean() +
  ggsci::scale_fill_igv() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = "Cell Proportion(%)",
    fill = "Cell Type"
  )
# export
ggsave(
  p_prop_celltype,
  filename = file.path("Celltype Props distribution.pdf"),
  height = 6,
  width  = 5
)


## find celltype markers -------
Idents(all_data) <- "sketch_snn_res.0.6"
DefaultAssay(all_data) <- "sketch"
all_data$celltype<-paste0(all_data$project,'_',all_data$annotation_full)
index<-all_data%>%
  .[[c('annotation_full','celltype')]]%>% 
  {rownames(.)<-NULL;.}%>%distinct()%>%
  group_by(annotation_full)%>%
  nest()%>%
  deframe()%>%map(1)

# index%>%imap(data = all_data, function(sample, name, all_data){
#   temp<-JoinLayers(all_data)%>%
#     subset(annotation_full == name)%>%
#     {Idents(.)<-'celltype';.}
#   sample %>% 
#     map(function(cell){
#       marker_df<-FindMarkers(temp,ident.1 = cell, ident.2 = setdiff(sample, cell))
#       write.csv(marker_df,file.path(mk_output,str_glue('{cell}_diff_marker.csv')))
#     })
# })

index%>%iwalk(data = JoinLayers(all_data), function(sample,name,data){
  Idents(data) <- 'celltype'
  sample %>% 
    walk(function(cell){
      marker_df<-FindMarkers(data,ident.1 = cell, ident.2 = setdiff(sample, cell),
                             test.use = "wilcox",
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             only.pos = TRUE)
      write.csv(marker_df,file.path(mk_output,str_glue('{cell}_diff_marker.csv')))
    })
})

## find markers -------
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


all_markers %>% 
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.25) %>%
  dplyr::filter(abs(avg_log2FC) > 1) %>%
  # arrange(desc(avg_log2FC)) %>%
  arrange(desc(pct.1),.by_group = TRUE) %>%
  slice_head(n = 10) %>%
  nest() %>%
  deframe() %>%
  map(pull, 'gene') %>%
  imap(~{
    xx <- str_c(.x, collapse = ",")
    str_glue("cluster_{.y}:{xx}")
  }) %>%
  walk(cat, "\n")


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
    features = all_top5_gene %>% unique()
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
  # arrange(desc(avg_log2FC)) %>%
  arrange(desc(pct.1),.by_group = TRUE) %>%
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


# spatial plot --------------
##ROI spatial-------------------
coords_enrichment <-
  all_data %>%
  {left_join(
    .[[c('annotation_full', "project")]] %>% rownames_to_column("cell"),
    GetTissueCoordinates(.),
    by = 'cell'
  )}

p_spatial_enrichment_ls <-
  coords_enrichment %>%
  group_by(project) %>%
  nest() %>%
  deframe() %>%
  imap(function(sub_coords, name){ # browser()
    ratio <- sub_coords[c('x', 'y')] %>% map(range) %>% map(diff) %>% {.[[1]] / .[[2]]}
    
    sub_coords %>%
      ggplot(aes(x = y, y = -x)) +
      geom_point(aes(color = annotation_full), shape = 15, size = 0.2) +
      # facet_wrap(~ Group, scales = "free") +
      coord_fixed(ratio = 1) +
      
      ggsci::scale_color_igv() +
      theme_void() +
      ggtitle(name)
  })


# all_data_roi %>%
#   .[, all_data_roi$Group == 'Tcell_riched_area1'] %>%
#   SpatialDimPlot(group.by = "annotation_full", alpha = 1, pt.size.factor = 50) +
#   ggsci::scale_fill_igv() 


# export plot
pdf(
  file.path(output2, "3. ROI spatial plot.pdf"),
  height = 8,
  width  = 8
)
p_spatial_enrichment_ls %>% walk(print)
dev.off()

##spatial plot-----------------------
p_sp_plot_anno <-
  Images(all_data) %>%
  setNames(., c('DZZ_Saline_0um','DZZ_Saline_100um','SYZ_NH4_0um','SYZ_NH4_100um')) %>%
  imap(function(img, img_name){
    p <- SpatialDimPlot(
      all_data,
      images = img, 
      image.scale = "hires", 
      group.by = "annotation_full",
      pt.size.factor	= 10
    ) +
      #ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })

ll<-all_data[,!is.na(all_data$first_type)]
ll2<-all_data[,!is.na(all_data$roi_group)]
p_sp_plot_anno <-
  Images(ll2)[5:6] %>%
  setNames(., c('DZZ_Saline','SYZ_NH4')) %>%
  imap(function(img, img_name){
    ll2@images[[img_name]] <- 
      ll2@images[[img_name]] %>%
      subset(cells = colnames(ll2))
    p <- SpatialDimPlot(
      ll2,
      images = img, 
      image.scale = "hires", 
      group.by = "first_type",
      pt.size.factor	= 3,
      crop  = FALSE
    ) +
      #ggsci::scale_fill_aaas() +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values=col_vector)+
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })
pdf(
  file.path(rctd_out_dir, 'ROI spatial plot of RCTD.pdf'),
  height = 10, width = 13
)
p_sp_plot_anno %>% walk(print)
dev.off()
#### Test by SZT 2024.11.18

# p$coordinates$ratio <- p$data[c("x", 'y')] %>% map(range) %>% map(diff) %>% {.[[2]]/.[[1]]}
# p$coordinates <- coord_cartesian()

 



p <- p_sp_plot_anno$DZZ_Saline_100um
g <- ggplotGrob(p)

xy_range <- p$data[c('x', 'y')] %>% map(range) %>% map(diff)

# x: height
# y: width

w <- xy_range[[2]] / max(unlist(xy_range)) 
h <- xy_range[[1]] / max(unlist(xy_range))  



g2 <- g %>% gtable::gtable_filter("guide-box") %>% .$grobs %>% discard(function(x) 'zeroGrob' %in% class(x)) %>% .[[1]] %>% .$grobs %>% .[[1]]
g1 <- g %>% gtable::gtable_filter("panel") %>% .$grobs %>% .[[1]] 


g1$children$geom_spatial.gTree.547$children$GRID.rastergrob.544

getChildOfGrob <- function(x, pattern){
  all_grobs <- names(x$children)
  idx <- str_detect(all_grobs, pattern = pattern)

  if(sum(idx) > 1){
    warning(str_glue("Found more than 1 grob which match the pattern {pattern}"), immediate. = TRUE)
    idx <- which(idx)[1]
  }
  if(sum(idx) < 1) stop(str_glue("No one grob match the pattern {pattern}") )

  grob_<- all_grobs[idx]

  x$children[[grob_]]
}

vp_g1 <-
  g1 %>%
  getChildOfGrob("geom_spatial") %>%
  getChildOfGrob("rastergrob") %>%
  .$vp



x0 <- as.numeric(vp_g1$x)
y0 <- as.numeric(vp_g1$y)

w0 <- as.numeric(vp_g1$width) - x0
h0 <- as.numeric(vp_g1$height) - y0



p_g1 <- cowplot::ggdraw() +
  cowplot::draw_grob(
    g1,
    x = 0,
    y = 0,
    width  = w, 
    height = h,
    hjust = x0/w0,
    vjust = y0/h0,
    scale  = min(c(1/w0, 1/h0)) * 2 - 0.1
  )  +
  coord_fixed()

p_g2 <- cowplot::ggdraw() +
  cowplot::draw_grob(
    g2,
    x = 0,
    y = 0,
    width = w,
    height = h 
  )

cowplot::plot_grid(p_g1, p_g2, rel_widths = c(2, 1))

#### ####


p_sp_plot_cluster <-
  Images(all_data) %>%
  setNames(., c('DZZ_Saline_0um','DZZ_Saline_100um','SYZ_NH4_0um','SYZ_NH4_100um')) %>%
  imap(function(img, img_name){
    p <- SpatialDimPlot(
      all_data,
      images = img, 
      image.scale = "hires", 
      group.by = "sketch_snn_res.0.6_full",
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
  height = 5, 
  width  = 5
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

##Calb1 gene-----------------

Images(all_data) %>%
  setNames(.,c('DZZ_Saline_0um','DZZ_Saline_100um','SYZ_NH4_0um','SYZ_NH4_100um')) %>%
  iwalk(function(img, img_name){
    pdf(file.path(er_output, str_glue("Calb1 gene spatial - {img_name}_2.pdf")),width = 8, height = 8)
    print(SpatialFeaturePlot(
      all_data,
      features = 'Calb1',
      images = img,
      image.scale = "hires",
      pt.size.factor	= 4,
      alpha = c(0,1),
      image.alpha=0.5
    ))
    dev.off()
  })
pdf(
  file.path(er_output, "Calb1 gene umap.pdf"),
  height = 8, 
  width  = 8
)
FeaturePlot(
  all_data,
  features = 'Calb1',
  reduction = "harmonyumap",
  raster = TRUE,
  max.cutoff = 2
)  
dev.off()

##celltype spatial-------------
Idents(all_data) <- "annotation"
cells <- CellsByIdentities(all_data)
Images(all_data) %>%
  setNames(.,c('DZZ_Saline_0um','DZZ_Saline_100um','SYZ_NH4_0um','SYZ_NH4_100um')) %>%
  iwalk(function(img, img_name){
    pdf(file.path(er_output, str_glue("celltype spatial - {img_name}2.pdf")),width = 8, height = 8)
    print(SpatialDimPlot(all_data,
                         cells.highlight = cells[setdiff(names(cells), "NA")],
                         cols.highlight = c("#FFFF00", "grey50"),
                         facet.highlight = T, combine = T) + NoLegend())
    dev.off()
  })


table(all_data$annotation)%>%names()%>%
 walk(function(name){
    pdf(
      file.path(er_output, str_glue("celltype umap - {name}.pdf")),
      height = 8, 
      width  = 8
    )
    DimPlot(all_data, group.by = "annotation",
            cells.highlight = WhichCells(all_data, 
                                         idents = name), 
            cols.highlight = c("#DC0000B2"), cols = "grey",label=F,reduction = 'harmonyumap')+NoLegend() + labs(title=name) 
    dev.off()
  })

DimPlot(all_data, group.by = "annotation",  
        cells.highlight = WhichCells(all_data, idents = names(table(all_data$annotation))[1]),  
        cols.highlight = c("#DC0000B2"), cols = "grey", label = FALSE, reduction = 'harmonyumap') +  
  NoLegend() +  
  labs(title = names(table(all_data$annotation))[1]) 


# 遍历每个注释类别并生成 DimPlot  
plots<-map(names(table(all_data$annotation)), function(name) {  
  p <- DimPlot(all_data, group.by = "annotation",  
               cells.highlight = WhichCells(all_data, idents = name),  
               cols.highlight = c("#DC0000B2"), cols = "grey", label = FALSE, reduction = 'harmonyumap') +  
    NoLegend() +  
    labs(title = name)
  return(p)
})  
combined_plot <- patchwork::wrap_plots(plots[1:13])
ggsave(file.path(er_output, "celltype_split_umap.pdf"), combined_plot, width = 16, height = 16)

##第二次分析
# ROI分析 ---------------------
roi_out_dir <- "./output_sample/ROI_analysis" # dir.create(roi_out_dir)

all_data <- readRDS("all_data_projected.rds")

roi_df <- list.files(
  '/home/rstudio/data/ROI_new',
  full.names = TRUE
) 

names(roi_df) <- 
  roi_df %>%
  basename() %>%
  str_remove("\\..*")

roi_df <-
  list.files('/home/rstudio/data/ROI_new', full.names = TRUE, pattern = "csv") %>%
  setNames(., c("DZZ_Saline_34_0um","DZZ_Saline_34_100um","SYZ_NH4_44_0um","SYZ_NH4_44_100um")) %>%
  map(read.csv) %>%
  enframe("sample") %>%
  unnest(value) %>%
  mutate(Barcode = str_c(sample, Barcode, sep = "_")) %>%
  # mutate(roi_group = str_c(group, ROI, sep = "_")) %>%
  dplyr::select(group, Barcode)  %>%
  column_to_rownames("Barcode")  %>%
  rename(roi_group=group)

roi_ls <- 
  split(
    rownames(roi_df),
    roi_df$roi_group
  )

all_data <- all_data %>% AddMetaData(roi_df)

roi_assay <- all_data[['Spatial.008um']] %>% subset(cells = rownames(roi_df))
all_data[['ROI']] <- roi_assay


DefaultAssay(all_data) <- "ROI"



# Images(all_data)
# 
# SpatialDimPlot(all_data, images = "slice1.008um.2", cells.highlight = roi_ls$KA_ROI_2)
# 
# 
# 
# GetTissueCoordinates(all_data[['slice1.008um.2']]) %>%
#   dplyr::filter(cell %in% roi_ls$KA_ROI_2) %>%
#   ggplot(aes(x = x, y = y)) +
#   geom_point()
# 

## norm annd pca ------
all_data <- NormalizeData(all_data)
all_data <- FindVariableFeatures(all_data)
all_data <- ScaleData(all_data, vars.to.regress = c('nCount_Spatial.008um','percent.mt'))


all_data <- RunPCA(all_data)

pdf(file.path(roi_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data, ndims = 30)
dev.off()


set.seed(1234)
suset_cells <- 
  all_data[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data) ) %>%
  dplyr::pull("cell_id")

pdf(file.path(roi_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 1:5)
dev.off()

pdf(file.path(roi_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 6:10)
dev.off()

pdf(file.path(roi_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 11:15)
dev.off()

pdf(file.path(roi_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data[, suset_cells], dims = 16:20)
dev.off()


## dr ------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

dims <- 10
all_data <- RunUMAP(all_data, dims = 1:dims, return.model = TRUE)


pdf(file.path(roi_out_dir, str_glue("Dimplot dim{dims} project.pdf")), height = 6, width = 8)
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
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_"
)

pdf(file.path(roi_out_dir, str_glue("Dimplot dim{dims} project - harmony.pdf")), height = 6, width = 8)
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
pdf(file.path(roi_out_dir, str_glue("Dimplot dim{dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
# DimPlot(all_data, group.by = "ROI_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "ROI_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "ROI_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "ROI_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data, group.by = "ROI_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data, group.by = "ROI_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "ROI_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data, group.by = "ROI_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


#prop barplot-------------
coords_enrichment <-
  all_data %>%
  {left_join(
    .[[c('first_type', "roi_group")]] %>% rownames_to_column("cell"),
    GetTissueCoordinates(.),
    by = 'cell'
  )}
coords_enrichment<-
  coords_enrichment%>% 
  filter(!is.na((.$roi_group)))
props_df <-
  coords_enrichment %>%
  group_by(roi_group, first_type) %>%
  summarise(count = n()) %>%
  group_by(roi_group) %>%
  mutate(freq = count/sum(count) * 100) 
p_prop_celltype <-
  props_df %>%
  
  ggplot(aes(x = roi_group, y = freq, fill = first_type)) +
  geom_col() +
  ggpubr::theme_pubclean() +
  ggsci::scale_fill_igv() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = "Cell Proportion(%)",
    fill = "Cell Type"
  )
# export
ggsave(
  p_prop_celltype,
  filename = file.path("./RCTD/ROI Celltype Props distribution.pdf"),
  height = 6,
  width  = 8
)



# patch code -----------------------
dr_outdir2 <- str_glue("{dr_outdir} - 2") 
dir.create(dr_outdir2)

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

