import scanpy as sc

adata = sc.read("/mnt/NexData_file/home/lst/project/AE2511186665/售后1/adata_death_sig_new.h5ad")

adata

adata.X.max()
adata.X.min()


