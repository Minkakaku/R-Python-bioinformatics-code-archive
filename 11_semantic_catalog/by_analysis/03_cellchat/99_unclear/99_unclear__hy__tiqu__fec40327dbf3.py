import scanpy as sc

# 1. 读取h5ad文件
adata = sc.read_h5ad("adata_with_band_by_sample.h5ad")

# 2. 查看数据的维度
print(f"原始数据形状: {adata.shape}")
print(f"site_band列的唯一值: {adata.obs['site_band'].unique()}")

# 3. 提取out_a和out_d的数据
# 方法A：使用布尔索引

# 方法B：同时提取两者到一个新adata
adata_combined = adata[adata.obs['site_band'].isin(['out_a', 'in_d']), :]

# 4. 查看提取结果
print(f"合并数据形状: {adata_combined.shape}")

# 5. 保存提取的数据
adata_combined.write("adata_out_a_and_in_d.h5ad")
