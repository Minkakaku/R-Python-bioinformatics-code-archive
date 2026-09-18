# /home/ztt/projects/FWSP20241029_xenium5k_liver_cancer_report

import os
import scanpy as sc
import pandas as pd

file_dir = "/home/ztt/projects/FWSP20241029_xenium5k_liver_cancer_report/data"
os.listdir(file_dir)

adata = sc.read(os.path.join(file_dir, 'adata_celltype.h5ad'))

adata.obs.to_csv("proximity_out_dir/adata_cellty_obs.csv")

coords_df = pd.DataFrame(adata.obsm['spatial'], columns = ['X', 'Y'])
coords_df.to_csv("proximity_out_dir/coords.csv")


