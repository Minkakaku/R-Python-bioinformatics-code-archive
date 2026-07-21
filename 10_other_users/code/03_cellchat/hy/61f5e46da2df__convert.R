library(reticulate)
use_python("/home/szt/miniforge3/envs/tools/bin/python")

adata_celltype_anno <-
  sceasy::convertFormat(
    "all_pCR_NpCR.h5ad", 
    from = "anndata", 
    to = "seurat",
    main_layer = "counts" # <----- 
  )

saveRDS(adata_celltype_anno, 'all_pCR_NpCR.rds')

