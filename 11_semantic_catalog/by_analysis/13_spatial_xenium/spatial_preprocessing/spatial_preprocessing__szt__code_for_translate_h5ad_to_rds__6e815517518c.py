# cd /mnt/NexData_file/home/szt/others/20260131_hd_spatial_temp
# conda activate tools
# python

import scanpy as sc
import os

# norm
os.system('/home/szt/miniforge3/envs/tools/bin/Rscript --no-save --no-restore code_for_translate_h5ad_to_rds.R adata_celltype.h5ad data')

# raw data
os.system('/home/szt/miniforge3/envs/tools/bin/Rscript --no-save --no-restore code_for_translate_h5ad_to_rds.R adata_filter.h5ad counts')


