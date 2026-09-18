library(tidyverse)
library(Seurat)


out_dir <- "output" 
dir.create(out_dir)

# all_data <- readRDS("all_data_projected_v2.rds")
all_data <- readRDS("all_data_projected_v2_hi_res.rds")


# roi coords
coords_roi_ls <-
  list.files("data", full.names = TRUE, pattern = "csv") %>%
  {
    x <- .
    sample <- x  %>% basename() %>% str_remove_all("^coords_|\\.csv$")
    # print(sample)
    names(x) <- sample_df$group[match(sample, sample_df$sample)]
    x
  } %>%
  map(read.csv) %>%
  map(1) %>%
  imap(function(x, name){
    x <- str_c(name, x, sep = "_")
    x
  })


all_data_tissue_area <-
  all_data %>%
  subset(cells = unlist(coords_roi_ls))

DefaultAssay(all_data_tissue_area) <- "Spatial.008um"


#1 QC  plot -------------------

## 1.1 umapplot---------------

pdf(
  file.path(out_dir, str_glue("1. Dimplot annotation - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data_tissue_area, group.by = "annotation_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
DimPlot(all_data_tissue_area[, !is.na(all_data_tissue_area$annotation_full2)], group.by = "annotation_full2", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()



## 1.2 spatial plot --------------
sp_anno_plot_ls <-
  Images(all_data_tissue_area) %>%
  setNames(., str_remove(., "\\..*$")) %>%
  imap(function(image_name, name){ # browser()

    if(name == "NML") name <- "HC"
    
    flip_x  <- FALSE
    flip_y  <- FALSE
    flip_xy <- FALSE
    rotate_r90 <- FALSE
    rotate_90  <- FALSE
    
    if(name == "BP3"){
      rotate_r90 <- TRUE
    }
    
    # plot
    p <-
      all_data_tissue_area %>%
      SpatialPlot_manual(
        image       = image_name,
        image.scale = "hires",
        # features    = 'CD4',
        group.by    = "annotation_full",
        
        flip_x  = flip_x,
        flip_y  = flip_y,
        flip_xy = flip_xy,
        rotate_r90 = rotate_r90,
        rotate_90  = rotate_90,
        plot_he = TRUE,
        alpha = c(1, 1),
        pt.size.factor = 0.5
      ) +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 3))) +
      ggtitle(name)
    
    p
  })

# sp_anno_plot_ls

pdf(
  file.path(out_dir, '1. spatial plot of cell annotation.pdf'), 
  height = 10, 
  width = 10
)
sp_anno_plot_ls %>% walk(print)
dev.off()


### remove Other cells ----------
sp_anno_plot_remove_other_ls <-
  Images(all_data_tissue_area) %>%
  setNames(., str_remove(., "\\..*$")) %>%
  imap(function(image_name, name){ # browser()
    
    if(name == "NML") name <- "HC"
    
    flip_x  <- FALSE
    flip_y  <- FALSE
    flip_xy <- FALSE
    rotate_r90 <- FALSE
    rotate_90  <- FALSE
    
    if(name == "BP3"){
      rotate_r90 <- TRUE
    }
    
    remove_other_cells <- Cells(all_data_tissue_area)[all_data_tissue_area$annotation_full != "Other"]
    
    # plot
    p <-
      all_data_tissue_area %>%
      subset(cells = remove_other_cells) %>%
      SpatialPlot_manual(
        image       = image_name,
        image.scale = "hires",
        # features    = 'CD4',
        group.by    = "annotation_full",
        
        flip_x  = flip_x,
        flip_y  = flip_y,
        flip_xy = flip_xy,
        rotate_r90 = rotate_r90,
        rotate_90  = rotate_90,
        plot_he = TRUE,
        alpha = c(1, 1),
        pt.size.factor = 0.5
      ) +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 3))) +
      ggtitle(name)
    
    p
  })

# sp_anno_plot_ls

pdf(
  file.path(out_dir, '1. spatial plot of cell annotation - remove other cells.pdf'), 
  height = 10, 
  width = 10
)
sp_anno_plot_remove_other_ls %>% walk(print)
dev.off()


### remove Other cells - ggplot2 ----------
remove_other_cells <- Cells(all_data_tissue_area)[all_data_tissue_area$annotation_full != "Other"]

sp_anno_ls2 <-
  all_data_tissue_area@images %>%
  imap(v = subset(all_data_tissue_area, cells = remove_other_cells), 
       function(img, img_name, v = v){ # browser()
    coords  <- GetTissueCoordinates(img)
    anno_df <- v[['annotation_full']]
    
    inter_cells <- intersect(rownames(coords), rownames(anno_df))
    
    plot_data <-
      cbind(
        coords[inter_cells,  ,drop=FALSE],
        anno_df[inter_cells, ,drop=FALSE]
      ) %>%
      as.data.frame()
    
    # plot
    plot_data %>%
      # slice_sample(prop = 0.1) %>%
      mutate(annotation_full = fct_recode(annotation_full, Endothelial = "FABP4 Cells")) %>%
      ggplot(aes(x = y, y = -x)) +
      geom_point(aes(color = annotation_full), shape = 15, size = 1) +
      # geom_bin2d(aes(fill = annotation_full), binwidth = 100, size = 1) +
      theme_void() +
      ggsci::scale_color_igv()
      # ggsci::scale_fill_igv()
  })



pdf(
  file.path(out_dir, '1. spatial plot of cell annotation - only points - remove other cells.pdf'), 
  height = 10, 
  width = 10
)
sp_anno_ls2 %>% walk(print)
dev.off()



## 1.3 Dotplot -----------------
skin_marker_ordered <- c(
  'Eccrine gland' = "SCGB1D2",
  'Eccrine gland' = "PIP",
  'Eccrine gland' = "SCGB1B2P",
  'Eccrine gland' = "SCGB2A2",
  'Eccrine gland' = "MUCL1",
  'Eccrine gland' = "DCD",
  
  
  Endothelial = "CLDN5",
  # Endothelial = "FABP4",
  Endothelial = "PECAM1",
  Endothelial = "CDH5",
  Endothelial = "TM4SF1",
  Endothelial = "CCL21",
  
  'FABP4 Cells' = "FABP4",
  
  Fibroblast = "DCN",
  Fibroblast = "COL1A1",
  Fibroblast = "CFD",
  Fibroblast = "COL1A2",
  Fibroblast = "COL3A1",
  Fibroblast = "APOD",

  
  Keratinocyte = 'KRT14',
  Keratinocyte = 'KRT1',
  Keratinocyte = 'DMKN',
  Keratinocyte = 'KRT10',
  Keratinocyte = 'KRT5',
  Keratinocyte = 'KRTDAP',
  

  'Myeloid cell' = 'HLA-DRA',
  'Myeloid cell' = 'HLA-DPB1',
  'Myeloid cell' = 'CD74',
  'Myeloid cell' = 'LYZ',
  'Myeloid cell' = 'HLA-DPA1',
  'Myeloid cell' = 'HLA-DRB1',

  
  "Plasma" = "IGHG1",
  "Plasma" = "IGHA1",
  "Plasma" = "IGHG3",
  "Plasma" = "IGKC",
  "Plasma" = "TNFRSF17",
  "Plasma" = "JCHAIN",
  
  
  'Smooth muscle' = "TAGLN",
  'Smooth muscle' = "ACTA2",
  'Smooth muscle' = "MYL9",
  'Smooth muscle' = "RGS55",
  'Smooth muscle' = "TPM2",
  'Smooth muscle' = "CALD1",
  
  'Nerve cell' = 'MPZ',
  'Nerve cell' = 'PLP1',
  'Nerve cell' = 'S100B',
  'Nerve cell' = 'PMP22',
  'Nerve cell' = 'CRYAB',
  'Nerve cell' = 'MBP',
  
  'T cell' = 'IL32',
  'T cell' = 'CD52',
  'T cell' = 'CXCR4',
  'T cell' = 'CD3E',
  'T cell' = 'CD3D',
  'T cell' = 'TRAC',
  
  
  Melanocyte = 'DCT',
  Melanocyte = 'TYRP1',
  Melanocyte = 'PMEL',
  Melanocyte = 'MLANA',
  Melanocyte = 'QPCT',
  Melanocyte = 'MIFT',
  
  'Mast cell' = "TPSAB1",
  'Mast cell' = "CPA3",
  'Mast cell' = "CTSG",
  'Mast cell' = "HPGD",
  'Mast cell' = "HPGDS",
  
  Eosinophil = "ITGA4",
  Eosinophil = "IL5RA", 
  Eosinophil = 'CEACAM8',
  Eosinophil = 'FUT4',
  Eosinophil = 'CCR3',
  Eosinophil = 'CD125',
  Eosinophil = "PTGDR2",
  Eosinophil = "SIGLEC8",
  Eosinophil = "ITGAM"
)


pdf(
  file.path(out_dir, str_glue("1. dotplot annotation - projected.pdf")),
  width = 25,
  height = 4
)
  DotPlot(
    all_data_tissue_area, 
    assay = "Spatial.008um", # 'Spatial.008um', # "sketch", 
    group.by = "annotation_full", 
    features = skin_marker_ordered
  ) +
  coord_cartesian(clip = "off") +
  scale_color_gradient(low = "gray90", high = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()


# dotplot
pdf(
  file.path(out_dir, str_glue("1. dotplot annotation - projected - 2.pdf")),
  width = 15,
  height = 4
)
  DotPlot(
    all_data_tissue_area, 
    assay = 'Spatial.008um', # "sketch", 
    group.by = "annotation_full", 
    features = skin_marker_ordered %>%
      .[-seq(which(names(.) == "Nerve cell")[1], length(.))] %>%
      unname()
  ) +
  scale_color_gradient(low = "gray90", high = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
dev.off()




#3  RCTD -----------------------
library(spacexr)
rctd_out_dir <- 'RCTD' # dir.create(rctd_out_dir)

# all_data <- readRDS("all_data.full.annotation.rds")

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
# all_data_sc <- readRDS('BP_AD_Pso-v2.rds')
# all_data_sc <- readRDS('BP_AD_Pso_20240930.rds')
all_data_sc <- readRDS('BP_AD_Pso_20241028.rds')

all_data_sc[['RNA']] <- all_data_sc[['RNA']] %>% as("Assay5")
all_data_sc$annotation <- all_data_sc$celltype 

Idents(all_data_sc) <- "annotation"
counts_sc  <- all_data_sc %>% JoinLayers() %>% LayerData(layer = "counts")
cluster_sc <- as.factor(all_data_sc$annotation)
nUMI_sc    <- all_data_sc$nCount_RNA

sc_cell_filter_idx <-  cluster_sc != "Unknown"

# create the RCTD reference object
reference <- Reference(
  counts_sc[, sc_cell_filter_idx], 
  droplevels(cluster_sc[sc_cell_filter_idx]), 
  nUMI_sc[sc_cell_filter_idx]
)

# 3. create the RCTD query object
rctd_query_ls <-
  visim_ls %>%
  map(function(visium){# browser()
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
    RCTD <- create.RCTD(
      query, reference, max_cores = 10,
      UMI_min = 50
      # UMI_max
    )
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

# saveRDS(rctd_result_ls, 'rctd_result_ls.rds')

#5. addMetaData
visium_rctd_ls <-
  visim_ls %>%
  imap(function(visium, name){ # browser()
    rctd <- rctd_result_ls[[name]]
    visium <- AddMetaData(visium, metadata = rctd@results$results_df)
    
    visium$first_type <- as.character(visium$first_type)
    visium$first_type[is.na(visium$first_type)] <- "Other" # "Unknown"
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
  visium_rctd_ls %>%
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
visium_rctd_ls %>%
  imap(function(visium, name){
    if(name == "HC") name_fix <- "NML" else name_fix <- name
    # visium$first_type[visium$first_type == "Unknown"] <- NA
    visium$first_type[visium$first_type == "Other"] <- NA
    
    SpatialDimPlot(
      visium,
      group.by       = 'first_type',
      image.scale    = 'lowres',
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()

## RCTD remove non-tissue area -------
pdf(
  file.path(out_dir, '3. spatial plot of RCTD - only tissue area.pdf'), 
  height = 15, 
  width = 15
)
visium_rctd_ls %>%
  imap(function(visium, name){ # browser()
    # if(name == "HC") name_fix <- "NML" else name_fix <- name
    name_fix <- name
    
    visium$first_type[visium$first_type == "Other"] <- NA
    
    
    image  <- visim_image_he_ls[[name]]
    image  <- RenameCells(image, new.names = paste0(name, "_", Cells(image)))

    # coords <- visim_coords_ls[[name]]
    # coords_fixed <- paste(name, coords, sep = "_")

    # image_sub  <- subset(image,  cells = coords_fixed)
    # visium_sub <- subset(visium, cells = coords_fixed)
    # visium_sub@images[[1]] <- image_sub
    
    visium@images[[1]] <- image
    
    # all_celltype <- visium_sub$first_type %>% unique
    # cell_colors  <- structure(
    #   names = all_celltype,
    #   ggsci::pal_igv()(length(all_celltype))
    # )
    # cell_colors <- c(cell_colors, Other = NA)
    
    # remove others 
    # visium_sub <- visium_sub[, !is.na(visium_sub$first_type)]
    visium_sub <- visium[, !is.na(visium$first_type)]
    
    SpatialDimPlot(
      visium_sub,
      group.by       = 'first_type',
      image.scale    = 'hires',
      pt.size.factor = 5 * image_scale_ls[[name]],
      shape          = 22,
      # alpha          = 0.5
    ) +
      ggtitle(name_fix) +
      ggsci::scale_fill_igv() +
      # scale_fill_manual(values = cell_colors) +
      guides(fill = guide_legend(override.aes = list(size = 5)))
  })
dev.off()


# 4 gene feature plot ----------------
# roi coords
coords_roi_ls2 <-
  list.files("ROI - adjust name", full.names = TRUE) %>%
  setNames(., str_remove(basename(.), "\\.csv")) %>%
  map(read.csv) %>%
  map(1) %>%
  imap(function(x, name){
    if(name == 'NML') name <- "HC"
    x <- str_c(name, x, sep = "_")
    x
  })

# if(! 'HC' %in% names(coords_roi_ls2)){
#   coords_roi_ls2[['HC']] <- 
#     all_data %>%
#     .[[]] %>%
#     {rownames(.)[.$project == "HC"]}
# }

all_data_roi <- 
  all_data %>% 
  subset(cells = unlist(coords_roi_ls2))

DefaultAssay(all_data_roi) <- "Spatial.008um"


all_data %>% 
  GetTissueCoordinates(image = "HC.008um") %>%
  slice_sample(prop = 0.01) %>%
  ggplot(aes(x = y, y = -x)) +
  geom_point() +
  geom_point(
    data = all_data %>% 
      GetTissueCoordinates(image = "HC.008um") %>%
      dplyr::filter(cell %in% coords_roi_ls2$HC) %>%
      dplyr::slice_sample(prop = 0.01),
    aes(x = y, y = -x),
    color = "red"
  )
  


# jpg image

gene_to_plot <- c(
  'CD3D',
  'CD69',
  # 'ITGAE',
  'CXCR6',
  
  'CD4',
  'CD8A',
  'IFNG',
  "IL13",
  "IL17A",
  "IL9",
  "IL22"
)

first_genes <- c(
  'CD3D',
  'CD69',
  # 'ITGAE',
  'CXCR6'
)

other_genes <- c(
  'CD4',
  'CD8A',
  'IFNG',
  "IL13",
  "IL17A",
  "IL9",
  "IL22"
)


# SpatialPlot_manual(
#   image.scale = "hires",
#   features = c("CD3D", "CD69", "CXCR6"),
#   
#   flip_x = FALSE,
#   flip_y = FALSE,
#   flip_xy = FALSE,
#   plot_he = TRUE,
#   alpha = c(0, 1),
#   pt.size.factor = 10,
# ) +
#   scale_fill_gradient(low = "#FFFFFF11", high = "red", trans = 'log1p') +
#   guides(shape = guide_legend(override.aes = list(size = 5)))

all_genes_for_ploting <-
  all_data_roi %>%
  FetchData(vars = gene_to_plot, slot = 'data')
max_val_all_plot <- max(all_genes_for_ploting) %>% ceiling()

data_first_plot <- 
  all_data_roi %>%
  FetchData(vars = first_genes, slot = 'data')

max_val_first_plot <- 
  data_first_plot %>%
  dplyr::summarise(across(everything(), range)) %>%
  max() %>%
  ceiling()

first_plot_ls <-
  Images(all_data_roi) %>%
  sort() %>%
  setNames(., str_remove(., "\\..*$")) %>%
  
  imap(function(image_name, name){ # browser()
    # NML 左旋90度
    # BP3 垂直翻转
    
    if(name == "NML") name <- "HC"
    
    flip_x  <- FALSE
    flip_y  <- FALSE
    flip_xy <- FALSE
    rotate_r90 <- FALSE
    rotate_90  <- FALSE
    
    if(name == "BP3"){
      rotate_r90 <- TRUE
    }
    
    # plot
    p <-
      all_data_roi %>%
      SpatialPlot_manual(
        image       = image_name,
        image.scale = "hires",
        features    = first_genes,
        
        flip_x  = flip_x,
        flip_y  = flip_y,
        flip_xy = flip_xy,
        rotate_r90 = rotate_r90,
        rotate_90  = rotate_90,
        plot_he = TRUE,
        alpha = c(0, 1),
        pt.size.factor = 8
      ) +
      scale_fill_gradient(
        low = "#FFFFFF11", high = "red", trans = 'log1p', 
        # limits = c(NA, max_val_first_plot),
        limits = c(NA, max_val_all_plot),
        guide = guide_colorbar(theme = theme(legend.key.width = unit(10, "lines")))
      ) +
      guides(shape = guide_legend(override.aes = list(size = 3))) +
      ggtitle(name)
    
    p
  })


# Others genes plot
other_plot_ls <-
  other_genes %>% 
  set_names(., .) %>%
  map(function(gene){ # browser()
    data_plot <- 
      all_data_roi %>%
      FetchData(vars = gene, slot = 'data')
    
    max_val_df <- 
      data_plot %>%
      cbind(all_data_roi[['project']]) %>%
      group_by(project) %>%
      dplyr::summarise(across(everything(), max)) 
    
    max_val <-
      max_val_df[[gene]] %>%
      max() %>%
      ceiling()
    
    
    cat('max_val: ', max_val, "\n")
    
    Images(all_data_roi) %>%
      sort() %>%
      setNames(., str_remove(., "\\..*$")) %>%
      imap(function(image_name, name){  # browser()
        # NML 左旋90度
        # BP3 垂直翻转
        
        max_val_grp <- max_val_df[[gene]] %>% .[max_val_df$project == name]
        if(max_val_grp == 0){
          alpha_val <- c(0, 0)
        }else{
          alpha_val <- c(0, 1)
        }
        
        
        if(name == "NML") name <- "HC"
        
        flip_x  <- FALSE
        flip_y  <- FALSE
        flip_xy <- FALSE
        rotate_r90 <- FALSE
        rotate_90  <- FALSE
        
        
        if(name == "BP3"){
          rotate_r90 <- TRUE
        }
        
        # plot
        p <-
          all_data_roi %>%
          SpatialPlot_manual(
            image       = image_name,
            image.scale = "hires",
            features    = gene,
            
            flip_x  = flip_x,
            flip_y  = flip_y,
            flip_xy = flip_xy,
            rotate_r90 = rotate_r90,
            rotate_90  = rotate_90,
            plot_he = TRUE,
            alpha   = alpha_val,
            pt.size.factor = 5
          ) +
          scale_fill_gradient(
            low = "#FFFFFF11", high = "red", trans = 'log1p', 
            # limits = c(NA, max_val),
            limits = c(NA, max_val_all_plot),
            # guide = guide_colorbar(theme = theme(legend.key.width = unit(10, "lines")))
          ) +
          # guides(shape = guide_legend(override.aes = list(size = 5))) +
          # ggtitle(name) +
          theme()
        
        p
      })
  })


all_p_ls <- c(
  first_plot_ls, 
  purrr::list_c(other_plot_ls)
)

# # IL22 BP2, 
# all_p_ls[[idx_row2col2[24]]] 
# iris %>% 
#   ggplot(aes(x = Sepal.Length, y = Sepal.Width, fill = Petal.Length)) +
#   geom_point() + 
#   guides(fill = guide_colorbar(theme = theme(legend.key.width  = unit(10, "lines"))) )


math_image_layout <- matrix(seq_along(all_p_ls), byrow = FALSE, nrow = 4)
idx_row2col2 <-
  math_image_layout %>% 
  t %>% 
  c()


# 第一列只保留shape legend，取消颜色
all_p_ls[[1]] <-  
  all_p_ls[[1]] +
  guides(
    fill = FALSE,
    shape = guide_legend(override.aes = list(size = 1.5))
  ) +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 5),
    legend.key.spacing = unit(2, "mm")
  )

all_p_ls[2:4] <- 
  all_p_ls[2:4] %>%
  map(~ .x + guides(fill = FALSE, shape = FALSE))


# 其他列只有最右上角保留legend，位置在右侧，其他都取消legend，添加title
for(i in 2:ncol(math_image_layout)){
  p_idx <- math_image_layout[, i]
  subtitle <- other_genes[i - 1]
  
  # first row
  if(i == ncol(math_image_layout)){
    # print(i)
    all_p_ls[[ p_idx[1] ]]  <-
      all_p_ls[[ p_idx[1] ]] +
      theme(
        legend.position = "right",
        legend.title = element_blank(),
      ) +
      guides(fill = guide_colorbar(theme = theme(legend.key.width  = unit(1, "lines"))) ) +
      labs(subtitle = subtitle)
    
  }else{
    
    all_p_ls[[ p_idx[1] ]] <-
      all_p_ls[[ p_idx[1] ]] +
      theme(legend.position = "none") +
      labs(subtitle = subtitle)
    
  }
  
  # others row
  all_p_ls[ p_idx[-1] ] <-
    all_p_ls[ p_idx[-1] ] %>%
    map(~ .x + theme(legend.position = "none"))
  
}



p_spatial_combined_all2 <- 
  patchwork::wrap_plots(
    all_p_ls[idx_row2col2],
    nrow = 4,
    # widths = c(2.5, rep(1, 7))
    widths = c(1, rep(1, 6), 1)
    # guides = "collect"
  ) 


# pdf(
#   file.path(sp_output, 'spatial plot of genes_all_combined_ROI_and_rotation.pdf'), 
#   height = 15, 
#   width = 25
# )
# print(p_spatial_combined_all2)
# dev.off()



tiff(file.path(out_dir, '4. spatial plot of genes_all_combined_ROI_and_rotation.tiff'),
     width  = 720 * 2 * 3.2, 
     height = 480 * 2 * 3, 
     units = "px", 
     pointsize = 12,
     compression = "lzw", # c("none", "rle", "lzw", "jpeg", "zip", "lzw+p", "zip+p"),
     bg   = "transparent",
     res = 300,  
     type = "cairo" # c("cairo", "Xlib", "quartz")
)
print(p_spatial_combined_all2)
# print(p)
dev.off()



all_data %>% 
  SpatialFeaturePlot(
    features = c("KRT5", "KRT10", "KRT14"),
    images = "HC.008um", # "BP2.008um",
    image.scale = "hires"
  ) &
  scale_fill_gradient(
    low = "#FFFFFF11", high = "red", trans = 'log1p'
  )
  



# p <- 
#   iris %>%
#   ggplot(aes(x = Sepal.Length, y = Sepal.Width)) +
#   geom_point() +
#   labs(title = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
# p
# 
# (p + p )&
#   coord_cartesian(clip = "off")

