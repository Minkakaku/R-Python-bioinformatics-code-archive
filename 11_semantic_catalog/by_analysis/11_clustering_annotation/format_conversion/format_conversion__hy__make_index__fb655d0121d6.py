import scanpy as sc

# 读取数据
adata = sc.read_h5ad('Epi_new.h5ad')
print(f"原始数据: {adata.shape[0]} 个细胞")

# 1. 保存原始索引
if 'original_index' not in adata.obs.columns:
    adata.obs['original_index'] = adata.obs.index
    print("已保存原始索引到'original_index'列")

# 2. 将cell_id列设置为索引
if 'cell_id' in adata.obs.columns:
    adata.obs.index = adata.obs['cell_id']
    print("已将cell_id列设置为索引")
else:
    print("警告: 未找到'cell_id'列")
    print(f"可用列: {list(adata.obs.columns)}")

# 3. 检查索引是否唯一
if not adata.obs.index.is_unique:
    print("发现重复索引，开始去重...")
    adata_nodup = adata[~adata.obs.index.duplicated(keep='first'), :].copy()
    print(f"去重后: {adata_nodup.shape[0]} 个细胞")
    print(f"移除了 {adata.shape[0] - adata_nodup.shape[0]} 个重复细胞")
else:
    print("没有重复索引")
    adata_nodup = adata

# 4. 保存结果
adata_nodup.write('Epi_new_nodup_a.h5ad')
print("已保存为 Epi_new_nodup_a.h5ad")
