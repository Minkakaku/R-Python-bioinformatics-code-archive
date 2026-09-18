import anndata as ad
import diopy
import os
import scanpy as sc

adata = ad.read_h5ad('/home/yrc/trans/result.h5ad')
adata = adata.raw.to_adata()
adata.obs_names_make_unique()
diopy.output.write_h5(adata,file = '/home/yrc/result.h5')