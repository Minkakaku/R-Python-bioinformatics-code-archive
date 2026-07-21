import scanpy as sc
import os

adata = sc.read_h5ad("/mnt/NexData_file/home/wmx/Project/SC2601080130/shouhou6/anno/adata_anno.h5ad")

print(adata.obs['cell_type'].value_counts())
print(adata.obs['group'].value_counts())

adata_tmp = adata.raw.to_adata()

adata_count = sc.read_h5ad('/mnt/NexData_file/scRNA/scVision_output/SC2601080130/03-Cell_filter/adata_counts.h5ad')

missing_cells = [c for c in adata_tmp.obs_names if c not in adata_count.obs_names]
missing_genes = [g for g in adata_tmp.var_names if g not in adata_count.var_names]

print("missing cells:", len(missing_cells))
print("missing genes:", len(missing_genes))


adata_count_sub = adata_count[adata_tmp.obs_names, adata_tmp.var_names].copy()

adata_tmp.X = adata_count_sub.X.copy()

adata_tmp.write_h5ad("/mnt/NexData_file/home/wmx/Project/SC2601080130/shouhou6/anno/adata_anno_count.h5ad")
outfile="/mnt/NexData_file/home/wmx/Project/SC2601080130/shouhou6/anno/adata_anno_count.h5ad"
rdsfile="/mnt/NexData_file/home/wmx/Project/SC2601080130/shouhou6/anno/adata_anno_count.rds"
os.system(f'/mnt/NexData_file/share/scRNA-tools/envs/pip_SCENICR/bin/Rscript /mnt/NexData_file/share/zwl/adata2seurat.R  --h5ad {outfile} --rds {rdsfile}')
