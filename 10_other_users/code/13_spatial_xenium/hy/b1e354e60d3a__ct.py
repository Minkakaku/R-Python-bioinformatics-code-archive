import anndata

# 读取h5ad文件
adata = anndata.read_h5ad('TMA.h5ad')

# 提取obs中的celltype列
celltype_series = adata.obs['celltype']

# 保存为csv，保留行索引（barcode）作为第一列
celltype_series.to_csv('ct.csv', index=True, header=['celltype'])
