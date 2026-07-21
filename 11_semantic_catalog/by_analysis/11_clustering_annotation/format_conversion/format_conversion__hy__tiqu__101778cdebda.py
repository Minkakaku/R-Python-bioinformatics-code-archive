import anndata as ad

# 读取原始数据
adata = ad.read_h5ad('All_filter_state_samp_new_with_CNs10-15.h5ad')

# 过滤：保留 celltype 为 'WNT+ ISC' 的细胞
adata_subset = adata[adata.obs['celltype'] == 'WNT+ ISC'].copy()

# 保存为新的 h5ad 文件
adata_subset.write('wnt_isc.h5ad')

print(f"原始细胞数: {adata.n_obs}")
print(f"过滤后细胞数: {adata_subset.n_obs}")
print("已保存为 wnt_isc.h5ad")
