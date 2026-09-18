# 1. 圈选ROI
# 2. 改为出图jpg
# 3. NML左旋90
# 4. BP3垂直翻转
# 5. 细胞注释不需要太细，只需要主群
# 6. dotplot也是只需要大类即可


library(tidyverse)
library(Seurat)


out_dir <- "output" 
dir.create(out_dir)


all_data.full.annotation <- readRDS("all_data.full.annotation.rds")

# visim_image_he_ls是hires HE
# 来自于dat_analysis.R文件的read high resolution HE
for(i in seq_along(visim_image_he_ls )){
  hi_he_img <- visim_image_he_ls[[i]]
  hi_he_img_name <- names(visim_image_he_ls)[i]
  
  if(hi_he_img_name == "NML"){
    hi_he_img_name <- "HC"
  }
  
  hi_he_img <- hi_he_img %>% RenameCells(str_c(hi_he_img_name, "_", Cells(.)))
  Key(hi_he_img) <- suppressWarnings(Key(hi_he_img_name))
  
  if(hi_he_img_name == "HC" || hi_he_img_name == "NML"){
    img_name <- "HC.008um"
  }else{
    img_name <- paste0(hi_he_img_name, ".008um")
  }
  
  all_data.full.annotation@images[[img_name]] <- hi_he_img
}


# roi coords
coords_roi_ls <-
  list.files("ROI", full.names = TRUE) %>%
  setNames(., str_remove(basename(.), "\\.csv")) %>%
  map(read.csv) %>%
  map(1) %>%
  imap(function(x, name){
    if(name == 'NML') name <- "HC"
    x <- str_c(name, x, sep = "_")
    x
  })

if(! 'BP3' %in% names(coords_roi_ls)){
  coords_roi_ls[['BP3']] <- 
    all_data.full.annotation %>%
    .[[]] %>%
    {rownames(.)[.$project == "BP3"]}
}

all_data_roi <- 
  all_data.full.annotation %>% 
  subset(cells = unlist(coords_roi_ls))

DefaultAssay(all_data_roi) <- "Spatial.008um"

# 0. qc metric and gene ncount -------------
all_data.full.annotation$project %>% table
# BP1    BP2    BP3     HC 
# 554582 487245 526413 252852 

all_data.full.annotation[[]] %>%
  group_by(project) %>%
  summarise(
    nCount = median(nCount_Spatial.008um),
    nFeature = median(nFeature_Spatial.008um)
  )

# project nCount nFeature
# <chr>    <dbl>    <dbl>
#   1 BP1       35.3     28.6
# 2 BP2       29.9     26.0
# 3 BP3       20.1     16.7
# 4 HC        40.6     33.1

# 1. dotplot及annotation只需要主群 -----------------

all_data_roi$annotation3_full <-
  all_data_roi$annotation2_full %>%
  fct_recode(
    Keratinocytes = 'Keratinocytes diff',
    Keratinocytes = 'Keratinocytes undiff',
    Other = 'TMSB4X supporting cells',
    Other = 'Undefined'
  )

all_data_roi2 <- all_data_roi[, all_data_roi$annotation3_full != "Other"]


major_marker_tidy <- c(
  
  Endothelial = "PECAM1",
  Endothelial = "VWF",
  Endothelial = "AQP1",
  
  
  'Fat cell' = 'ADIPOQ',
  'Fat cell' = 'FABP4',
  'Fat cell' = 'G0S2',
  
  
  "Fibroblasts" = "DCN",
  "Fibroblasts" = "VIM",
  "Fibroblasts" = "COL1A1",
  "Fibroblasts" = "COL1A2",
  
  
  'Follicular epithelial cell' =  'KRT17',
  'Follicular epithelial cell' =  'S100A2',
  
  
  # undiff
  'Keratinocytes' = 'KRT5',
  'Keratinocytes' = 'KRT14',
  # 'Keratinocytes' = 'TP63',
  # 'Keratinocytes' = 'ITGB1',
  # 'Keratinocytes' = 'ITGA6',
  
  # diff
  'Keratinocytes' = 'KRT1',
  'Keratinocytes' = 'KRT10',
  'Keratinocytes' = 'SBSN',
  'Keratinocytes' = 'KRTDAP',

  # Melanocytes = 'PMEL',
  # Melanocytes = 'MLANA',
  # Melanocytes = 'TYRP1',
  # Melanocytes = 'DCT',
  
  
  "Macrophage" = "CD68",
  "Macrophage" = "CD163",
  "Macrophage" = "C1QA",
  "Macrophage" = "C1QC",
  
  
  Plasma = 'IGHG1',
  Plasma = 'IGKC',
  
  
  'Sebaceous gland cell' = 'DCD',
  # 'Sebaceous gland cell' = 'FST',
  'Sebaceous gland cell' = 'MUCL1',
  
  'PTPRC',
  'T' = 'CD3D',
  'T' = 'CD3E',
  'T' = 'TRAC'
)




# dimplot
pdf(
  file.path(out_dir, str_glue("Dimplot dim{dims} annotation3 - harmony - projected.pdf")),
  height = 6,
  width = 8
)
DimPlot(all_data_roi2, group.by = "annotation3_full", label = TRUE, reduction = 'full.harmonyumap') + NoLegend()
dev.off()


# dotplot
pdf(
  file.path(out_dir, str_glue("dotplot dims{dims} - harmony - annotation3 unname-major_marker2.pdf")), 
  width = 16, 
  height = 4
)
all_data_roi2 %>%
  {DefaultAssay(.) <- 'sketch'; .} %>%
  DotPlot(group.by = "annotation3_full",
        features = unname(major_marker_tidy) %>% c("CD19", "MS4A1", "CD79A"),
        dot.scale = 6) +
  theme_90angle
dev.off()


# 1.1. 参考文献中的marker绘制dotplot ----------------
skin_marker_ls <- c(
  Keratinocyte = 'KRT14',
  Keratinocyte = 'KRT1',
  Keratinocyte = 'DMKN',
  Keratinocyte = 'KRT10',
  Keratinocyte = 'KRT5',
  Keratinocyte = 'KRTDAP',
  
  Melanocyte = 'DCT',
  Melanocyte = 'TYRP1',
  Melanocyte = 'PMEL',
  Melanocyte = 'MLANA',
  Melanocyte = 'QPCT',
  Melanocyte = 'MIFT',
  
  'Eccrine gland' = "SCGB1D2",
  'Eccrine gland' = "PIP",
  'Eccrine gland' = "SCGB1B2P",
  'Eccrine gland' = "SCGB2A2",
  'Eccrine gland' = "MUCL1",
  'Eccrine gland' = "DCD",
  
  Endothelial = "CLDN5",
  Endothelial = "FABP4",
  Endothelial = "PECAM1",
  Endothelial = "CDH5",
  Endothelial = "TM4SF1",
  Endothelial = "CCL21",
  
  Fibroblast = "DCN",
  Fibroblast = "COL1A1",
  Fibroblast = "CFD",
  Fibroblast = "COL1A2",
  Fibroblast = "COL3A1",
  Fibroblast = "APOD",
  
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
  
  'Myeloid cell' = 'HLA-DRA',
  'Myeloid cell' = 'HLA-DPB1',
  'Myeloid cell' = 'CD74',
  'Myeloid cell' = 'LYZ',
  'Myeloid cell' = 'HLA-DPA1',
  'Myeloid cell' = 'HLA-DRB1',
  
  'Mast cell' = "TPSAB1",
  'Mast cell' = "CPA3",
  'Mast cell' = "CTSG",
  'Mast cell' = "HPGD",
  'Mast cell' = "HPGDS"
)


all_data_roi2 %>%
  {DefaultAssay(.) <- 'sketch'; .} %>%
  DotPlot(group.by = "annotation3_full", # "immune.sub.cluster_full",
          features = skin_marker_ls, # %>% c("CD19", "MS4A1", "CD79A"),
          dot.scale = 6) +
  theme_90angle



# 1.2 Fat cell spatial plot -----------
p_sp_plot_fat_cell <-
  Images(all_data) %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    if(img_name == "HC") img_name <- "NML"
    
    all_data$annotation2_full[all_data$annotation2_full != "Fat cell"] <- NA
    
    p <- SpatialDimPlot(
      all_data,
      images = img, 
      group.by = "annotation2_full",
      image.scale = "hires",
      pt.size.factor	= 4
    ) +
      ggsci::scale_fill_aaas() +
      # ggsci::scale_fill_igv() +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(img_name)
    p
  })


# export 
pdf(
  'spatial plot of fat cell.pdf',
  height = 15,
  width  = 15
)
p_sp_plot_fat_cell %>% walk(print)
dev.off()



# 2. spatia plot of some genes ---------------
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
      guides(shape = guide_legend(override.aes = list(size = 5))) +
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
            guide = guide_colorbar(theme = theme(legend.key.width = unit(10, "lines")))
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



idx_row2col2 <-
  matrix(seq_along(all_p_ls), byrow = FALSE, nrow = 4) %>% 
  t %>% 
  c()


p_spatial_combined_all2 <- 
  patchwork::wrap_plots(
    all_p_ls[idx_row2col2],
    nrow = 4,
    widths = c(2.5, rep(1, 7))
    # guides = "collect"
  ) 


pdf(
  file.path(sp_output, 'spatial plot of genes_all_combined_ROI_and_rotation.pdf'), 
  height = 15, 
  width = 25
)
print(p_spatial_combined_all2)
dev.off()


# jpeg(
#   file.path(sp_output, 'spatial plot of genes_all_combined_ROI_and_rotation.jpg'),
#   width = 480, 
#   height = 480, units = "px", pointsize = 12,
#   quality = 75,
#   bg = "transparent", res = 300, 
#   type = c("cairo", "Xlib", "quartz")
# )

tiff(file.path(sp_output, 'spatial plot of genes_all_combined_ROI_and_rotation.tiff'),
     width  = 720 * 2 * 8, 
     height = 480 * 2 * 4, 
     units = "px", pointsize = 12,
     compression = "lzw", # c("none", "rle", "lzw", "jpeg", "zip", "lzw+p", "zip+p"),
     bg   = "transparent", res = 300,  
     type = "cairo" # c("cairo", "Xlib", "quartz")
)
print(p_spatial_combined_all2)
# print(p)
dev.off()




# all_data_roi %>%
#   SpatialPlot_manual(
#     image       = 'HC.008um',
#     image.scale = "hires",
#     features    = "CD4",
#     
#     flip_x  = F,
#     flip_y  = F,
#     flip_xy = F,
#     rotate_r90 = T,
#     rotate_90  = F,
#     plot_he = TRUE,
#     alpha   = c(0, 1),
#     pt.size.factor = 5
#   ) 


# undebug(SpatialFeaturePlot)

# pt.size.factor = 100
# Images(all_data_roi)
# # "HC.008um"  "BP1.008um" "BP2.008um" "BP3.008um"
# 
# all_data_roi %>% SpatialFeaturePlot(features = "CD4", images = "HC.008um", pt.size.factor = 100, alpha = c(0, 1) )
# all_data_roi %>% SpatialFeaturePlot(features = "CD4", images = "BP1.008um", pt.size.factor = 100, alpha = c(0, 1) )
# all_data_roi %>% SpatialFeaturePlot(features = "CD4", images = "BP2.008um", pt.size.factor = 100, alpha = c(0, 1) )
# all_data_roi %>% SpatialFeaturePlot(features = "CD4", images = "BP3.008um", pt.size.factor = 100, alpha = c(0, 1) )




# 2. spatia plot of some genes ---------------
# jpg image

gene_to_plot <- c(
  'CD3D',
  'CD69',
  'ITGAE',
  'CXCR6',
  
  'CD4',
  'CD8A',
  "IL13",
  "IL17A",
  "IL9",
  "IL22"
)

# table(gene_to_plot %in% rownames(all_data_roi2))
# 
# pdf(
#   file.path(sp_output, 'spatial plot of genes v2.pdf'), 
#   height = 15, 
#   width = 15
# )



# change Image in all_data_roi2 to hires.
for(img_idx in seq_along(visim_image_he_ls)){
  
  image_new <- visim_image_he_ls[[img_idx]]
  name <- names(visim_image_he_ls)[[img_idx]]
  
  image_name <- paste0(name, '.008um')
  all_image_name <- Images(all_data_roi2)
  
  if(!image_name %in% all_image_name && name == "NML"){
    name <- "HC"
    image_name <- paste0(name, '.008um')
  }
  
  image <- all_data_roi2@images[[image_name]]
  
  image_new  <- RenameCells(image_new, new.names = paste0(name, "_", Cells(image_new)))
  coords <- image_new %>% GetTissueCoordinates()
  cells <- coords$cell %>% .[. %in% colnames(all_data_roi2)]
  
  image_new  <- image_new %>% subset(cells = cells)
  Key(image_new) <- Key(image)
  
  all_data_roi2@images[[image_name]] <- image_new
}


# value range
norm_data <-
  all_data_roi2[['Spatial.008um']] %>%
  JoinLayers() %>%
  LayerData(layer = "data") 

norm_data <- norm_data %>% .[rownames(.) %in% gene_to_plot, ]
range_val <- norm_data %>% range %>% ceiling()



p_all_spatial_ls <-
  Images(all_data_roi2) %>%
  setNames(., str_remove(., "\\..*$")) %>%
  imap(function(image_name, name){ # browser()
    if(name == "HC") name <- "NML"
    
    flip_x <- FALSE
    flip_y <- FALSE
    
    if(name == 'HC' | name == "NML" | name == "BP3"){
      flip_y <- TRUE
    }
    
    p <- SpatialPlot2(
      all_data_roi2,
      images         = image_name,
      features       = gene_to_plot,
      image.scale    = 'hires',
      pt.size.factor = 5 * 2, # * image_scale_ls[[name]],
      shape          = 22,
      flip_x         = flip_x,
      flip_y         = flip_y,
      keep.scale     = 'all',
      combine        = FALSE,
      alpha          = c(0, 1) # 0.5
    )
    
    # NML: rotate 90 degree
    # flip_y在SpatialPlot2中设置
    if(name == "HC" | name == "NML"){
      p <- p %>% map( ~{ .x + coord_flip() } )
    }
    
    # BP3: flip vertical 
    # flip_y在SpatialPlot2中设置

    
    # 统一数据大小
    p <- 
      p %>%
      map(~{
        .x + scale_fill_gradient(low = "#FFFFFF11", high = "red", trans = 'log1p', limits = range_val) 
      })
    
    p
    
  }) %>%
  purrr::list_flatten()



g <- p[[1]] %>% ggplotGrob()
gg <- g %>% gtable::gtable_filter(pattern = "panel") %>% .$grobs %>% .[[1]]

grid.newpage()
gg$children[[1]] %>% grid.draw()

gg$children[[3]] %>% grid.draw()

gg$children[[3]]$children$GRID.rastergrob.23788 %>% grid.draw()

g2 <- gg$children[[3]]$children$GRID.rastergrob.23788
g2$vp$width <- unit(1, "npc")
g2$vp$height <- unit(1, "npc")

grid.draw(g2)


idx_row2col <- matrix(seq_along(p_all_spatial_ls), byrow = TRUE, ncol = 4) 

# images <-
#   Images(all_data_roi2) %>%
#   str_remove("\\..*$")

# # NML: rotate 90 degree
# nml_idx <- which(images == "HC" | images == "NML")
# 
# p_all_spatial_ls[idx_row2col[, nml_idx]] <- 
#   p_all_spatial_ls[idx_row2col[, nml_idx]] %>%
#   map(function(x){
#     x + 
#       coord_flip() +
#       scale_y_continuous(transform = "reverse")
#   })
#   
#   
# # BP3: flip vertical
# bp3_idx <- which(images == "BP3")
# 
# p_all_spatial_ls[idx_row2col[, bp3_idx]] <- 
#   p_all_spatial_ls[idx_row2col[, bp3_idx]] %>%
#   map(function(x){
#     x + 
#       scale_y_continuous(transform = "reverse")
#   })

p_spatial_combined_all <- 
  patchwork::wrap_plots(
    p_all_spatial_ls[idx_row2col],
    nrow = 4,
    # guides = "collect"
  ) &
  scale_fill_gradient(low = "#FFFFFF11", high = "red", trans = 'log1p') 


pdf(
  file.path(sp_output, 'spatial plot of genes_all_combined - v2.pdf'), 
  height = 15, 
  width = 25
)
print(p_spatial_combined_all)
dev.off()









