# # terminal
# conda activate tools
# cd /mnt/NexData_file/home/szt/20260105_XE2511186721-01_liver/output/data
# 
# R



# R env
library(reticulate)
use_python("/home/szt/miniforge3/envs/tools/bin/python")

args <- commandArgs(trailingOnly = T)
print(args)


main_layer = 'counts'
if(length(args) > 0){
  infile = args[1]
  stopifnot(grepl('h5ad$', infile))
  stopifnot(file.exists(infile))
  out_file = gsub('h5ad$', 'rds', infile)
  
  if (length(args) > 1){
    main_layer = args[2]
  }
}


adata <-
  sceasy::convertFormat(
    infile, 
    from = "anndata", 
    to = "seurat",
    main_layer = main_layer # <----- 
  )

adata$project <- adata$sample
adata$nCount_RNA <- adata$total_counts
adata$nFeature_RNA <- adata$nFeaturess_RNA_by_counts

saveRDS(adata, out_file)
