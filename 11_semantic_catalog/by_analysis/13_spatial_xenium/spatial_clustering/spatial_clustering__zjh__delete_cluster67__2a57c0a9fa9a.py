import scanpy as sc

adata1=sc.read('adata_cluster_clustered.h5ad')

res_cluster_to_remove = ['6','7']

if 'res_0.3' not in adata1.obs.columns:
    raise KeyError("错误: adata.obs 中不存在指定列。请检查列名是否正确。")

# **核心修改：使用 ~ 取反操作符，去除指定值的细胞，保留其余细胞**
adata = adata1[~adata1.obs['res_0.3'].isin(res_cluster_to_remove), :]

adata.write('adata_cluster_rm_cluster67.h5ad')
