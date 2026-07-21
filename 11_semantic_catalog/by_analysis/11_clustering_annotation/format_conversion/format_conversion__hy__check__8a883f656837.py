import scanpy as sc

# 读取数据
adata = sc.read_h5ad('Epi_new.h5ad')
print(f"原始数据: {adata.shape[0]} 个细胞")

# 检查索引是否唯一
if not adata.obs.index.is_unique:
    print("发现重复索引，开始去重...")
    # 去除重复索引，保留第一个
    adata_nodup = adata[~adata.obs.index.duplicated(keep='first'), :].copy()
    print(f"去重后: {adata_nodup.shape[0]} 个细胞")
    print(f"移除了 {adata.shape[0] - adata_nodup.shape[0]} 个重复细胞")
else:
    print("没有重复索引")
    adata_nodup = adata

# 保存结果
adata_nodup.write('Epi_new_nodup.h5ad')
print("已保存为 Epi_new_nodup.h5ad")
