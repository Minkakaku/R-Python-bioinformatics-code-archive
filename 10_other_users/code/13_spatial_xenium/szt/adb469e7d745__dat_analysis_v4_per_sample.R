library(Seurat)
library(tidyverse)

# setting dir -------
img_output <- "materials_for_rmd" # dir.create(img_output)
qc_outdir <- "1. QC_per_sample" # dir.create(qc_outdir)
dr_outdir <- "2. dimension reduction and cluster_per_sample"  # dir.create(dr_outdir)
sp_output <- "3. spatial plot_per_sample" # dir.create(sp_output)
mk_output <- "4. marker gene_per_sample" # dir.create(mk_output)
er_output <- "5. enriment analysis_per_sample" # dir.create(er_output)

# 4 11 20 hvg5000
# 12 hvg2000 regression out H cell


# set future -------------
library(future)
plan(multisession, workers = 10)
options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


# read data -------
sample_df <- 
  # read.csv("data/Annoation.csv")
  data.frame(
    sample = c('2310163', '2305735', '2304480', '2304434'),
    group  = c('HC', 'BP1', 'BP2', 'BP3')
  )

file_list <- 
  list.files(
    "/home/ztt/spatial/visium_HD/Human_Skin_FFPE/spaceranger_output",
    pattern = paste(sample_df$sample, collapse = "|"), 
    full.names = TRUE
  ) %>%
  str_subset("\\..*$", negate = TRUE)

stopifnot(nrow(file_list) == nrow(sample_df))

# TODO: fix bugs
# Samples order in file_list may not identical to sample_df
# names(file_list) <- sample_df$group # basename(file_list)
names(file_list) <- sample_df$group[match(basename(file_list), sample_df$sample)]




visim_ls <-
  file_list %>% 
  imap(function(path, name){
    path <- file.path(path, "outs")
    
    Read10X_Image <-
      function (image.dir, image.name = "tissue_lowres_image.png", 
                assay = "Spatial", slice = "slice1", filter.matrix = TRUE) {
        image <- png::readPNG(source = file.path(image.dir, image.name))
        scale.factors <- Read10X_ScaleFactors(filename = file.path(image.dir, 
                                                                   "scalefactors_json.json"))
        coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(image.dir, 
                                                                         "*tissue_positions*parquet")), filter.matrix)
        fov <- CreateFOV(coordinates[, c("imagerow", "imagecol")], 
                         type = "centroids", radius = scale.factors[["spot"]], 
                         assay = assay, key = Key(slice, quiet = TRUE))
        visium.fov <- new(Class = "VisiumV2", boundaries = fov@boundaries, 
                          molecules = fov@molecules, assay = fov@assay, key = fov@key, 
                          image = image, scale.factors = scale.factors)
        return(visium.fov)
      }
    
    Load10X_Spatial(
      path, 
      bin.size = 8, 
      slice = name, 
      image = Read10X_Image(file.path(path, "binned_outputs/square_008um/spatial"), assay = 'Spatial.008um')
    )
  })

# add sample_id or sample_name
visim_ls <-
  visim_ls %>%
  imap(function(visium, name){
    visium$project <- name
    Idents(visium) <- 'project'
    RenameCells(visium, add.cell.id = name)
  })


visim_ls <- 
  visim_ls %>%
  map(~ {DefaultAssay(.x) <- 'Spatial.008um'; .x}) %>%
  map(NormalizeData) 


# read high resolution HE ---------
visim_image_he_ls <-
  file_list %>% 
  imap(function(path, name){
    path <- file.path(path, "outs")
    
    Read10X_Image <-
      function (image.dir, image.name = "tissue_lowres_image.png", 
                assay = "Spatial", slice = "slice1", filter.matrix = TRUE) {
        image <- png::readPNG(source = file.path(image.dir, image.name))
        scale.factors <- Read10X_ScaleFactors(filename = file.path(image.dir, 
                                                                   "scalefactors_json.json"))
        coordinates <- Read10X_Coordinates(filename = Sys.glob(file.path(image.dir, 
                                                                         "*tissue_positions*parquet")), filter.matrix)
        fov <- CreateFOV(coordinates[, c("imagerow", "imagecol")], 
                         type = "centroids", radius = scale.factors[["spot"]], 
                         assay = assay, key = Key(slice, quiet = TRUE))
        visium.fov <- new(Class = "VisiumV2", boundaries = fov@boundaries, 
                          molecules = fov@molecules, assay = fov@assay, key = fov@key, 
                          image = image, scale.factors = scale.factors)
        return(visium.fov)
      }
    
    Read10X_Image(
      file.path(path, "binned_outputs/square_008um/spatial"),
      image.name = "tissue_hires_image.png"
    )
    
  })


# read coordinate ------------

visim_coords_ls <-
  file_list %>% 
  imap(function(path, name){ 
    proj <- basename(path)
    coord_file_name <- paste0("coords_", proj, ".csv")
    
    df <- read.csv(file.path("data", coord_file_name), header = TRUE)
    df$Barcode
  })



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
# all_data <- merge(visim_ls_pass[[1]], visim_ls_pass[-1])
all_data <- merge(visim_ls[[1]], visim_ls[-1])

## normalize for sketch data ------
all_data <- NormalizeData(all_data, scale.factor = floor(mean(all_data$nCount_Spatial.008um)))
all_data <- FindVariableFeatures(all_data)
