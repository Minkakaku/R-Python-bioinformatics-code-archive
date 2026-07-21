library(tidyverse)
library(Seurat)

all_data_full <- readRDS("all_data_full_v6.rds")

# output_dir <- "20260227 spatial plot"
output_dir <- "20260228 spatial plot"
dir.create(output_dir)


acp5_macro_gene <- c(
  'ACP5',
  'MMP12',
  'ATP6V0D2',
  'CCL24',
  'DNASE1L3',
  'ACE',
  'RARRES1',
  'OTOA',
  'CD276',
  'RAB7B',
  'SIGLEC12',
  'CD209',
  'SEL1L3',
  'BOLA3',
  'LSS',
  'SLC7A8',
  'ARPIN',
  'SRD5A3'
  # 'IGF1',
  # 'LTC4S',
  # 'CD163L1'
)

selenop_macro_gene <- c(
  'SELENOP',
  'IGF1',
  'LTC4S',
  'CD163L1',
  'GPR155',
  'MVB12B',
  'MS4A4E',
  'IQCN',
  'BNC2',
  'NFATC2',
  'GAL3ST4',
  'WLS',
  'EMP2',
  'DNASE1L3',
  'CXCL12'
)


# ACP5 ------------
p_sp_plot_acp5_macro_gene <-
  Images(all_data_full) %>% #  .[1] %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    if(img_name == "HC") img_name <- "NML"
    
    p <- SpatialFeaturePlot(
      all_data_full,
      images = img, 
      features = acp5_macro_gene, # "MMP12",
      image.alpha = 0.5,
      alpha = c(0, 1),
      pt.size.factor	= 10,
      combine = FALSE
    ) 
      # ggsci::scale_fill_aaas() +
      # ggsci::scale_fill_igv() +
      # guides(fill = guide_legend(override.aes = list(size = 5))) +
      # ggtitle(img_name)
    names(p) <- acp5_macro_gene
    p
  })

# p_sp_plot_macro_gene



# export for rmd
p_sp_plot_acp5_macro_gene %>%
  iwalk(function(p_ls, name){
    dir.create(file.path(output_dir, name))
    p_ls %>% iwalk(function(p, gene){
      png(file.path(output_dir, name, str_glue("spatial plot - marker gene of ACP5 macrophage - {gene}.png")), width = 720, height = 720)
      print(p)
      dev.off()
    })
   
  })



# VlnPlot(all_data_full, group.by = "sub_celltype2", features = "LRP1")

# SELENOP  ------------
p_sp_plot_selenop_macro_gene <-
  Images(all_data_full) %>% #  .[1] %>%
  setNames(., str_remove(., "\\.[0-9]{3}um$")) %>%
  imap(function(img, img_name){
    if(img_name == "HC") img_name <- "NML"
    
    p <- SpatialFeaturePlot(
      all_data_full,
      images = img, 
      features = selenop_macro_gene, # "MMP12",
      image.alpha = 0.5,
      alpha = c(0, 1),
      pt.size.factor	= 10,
      combine = FALSE
    ) 
    # ggsci::scale_fill_aaas() +
    # ggsci::scale_fill_igv() +
    # guides(fill = guide_legend(override.aes = list(size = 5))) +
    # ggtitle(img_name)
    names(p) <- selenop_macro_gene
    p
  })




# export for rmd
p_sp_plot_selenop_macro_gene %>%
  iwalk(function(p_ls, name){
    dir.create(file.path(output_dir, name))
    p_ls %>% iwalk(function(p, gene){
      png(file.path(output_dir, name, str_glue("spatial plot - marker gene of SELENOP macrophage - {gene}.png")), width = 720, height = 720)
      print(p)
      dev.off()
    })
    
  })

