import os, sys
import loompy as lp
import numpy as np
import anndata as ad
import scanpy as sc
output = "/mnt/NexData_file/home/zwl/Project/personal/demo/SCENIC/output"
ad_all = ad.read_h5ad('/mnt/NexData_file/scRNA/scVision_output/lung_cancer/anno.h5ad')
# 设置取样比例
sample_fraction = 0.3

# 计算要选取的细胞数量
n_samples = int(ad_all.n_obs * sample_fraction)

# 随机生成细胞索引
random_indices = np.random.choice(ad_all.n_obs, n_samples, replace=False)

# 根据随机索引选取细胞
adata = ad_all[random_indices, :]

sc.pp.filter_genes(adata,min_cells = 10)
adata.write('sub_anno.h5ad', compression="gzip")

row_attrs = {"Gene": np.array(adata.var_names),}
col_attrs = {"CellID": np.array(adata.obs_names)}
lp.create(os.path.join(output,"sce.loom"),adata.X.transpose(),row_attrs,col_attrs)

os.system(f'sh /mnt/NexData_file/home/zwl/Project/personal/demo/SCENIC/pyscenic.sh {os.path.join(output,"sce.loom")} {os.path.join(output,"grn_result.csv")} {os.path.join(output,"ctx_result.csv")} {os.path.join(output,"aue_result.loom")}')
