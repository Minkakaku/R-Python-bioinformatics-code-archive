library(tidyverse)
library(Seurat)
getwd()
library(patchwork)
options(java.parameters = '-Xmx10g')  # options('java.parameters')
library(rJava)
library(RBioFormats)
library(ggpubr)
library(SeuratDisk)
library(ggrepel)
library(spacexr)

reference <- readRDS("/mnt/NexData_file/home/szt/hyj/20250729_FWSP20242506-02_mouse_colon/result/RCTD20260604/reference.rds")
rctd_query_ls <- readRDS("/mnt/NexData_file/home/szt/hyj/20250729_FWSP20242506-02_mouse_colon/result/RCTD20260604/rctd_query_ls.rds")

#4. run RCTD-----
rctd_result_ls <- 
  rctd_query_ls %>%
  imap(function(query, name){
    # run RCTD
    RCTD <- create.RCTD(query, reference, max_cores = 10, UMI_min = 30, UMI_min_sigma = 50)
    RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
    
    RCTD
  })

saveRDS(rctd_result_ls, '/mnt/NexData_file/home/szt/hyj/20250729_FWSP20242506-02_mouse_colon/result/RCTD20260604/rctd_result_ls.rds')