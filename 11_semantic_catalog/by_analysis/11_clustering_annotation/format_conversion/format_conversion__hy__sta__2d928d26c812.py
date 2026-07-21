import scanpy as sc

# 1. 读取数据
adata = sc.read_h5ad("adata_celltype.h5ad")

# 2. 统计所有细胞类型的细胞数量
print("=== 所有细胞类型统计 ===")
celltype_counts = adata.obs['celltype'].value_counts()
print(celltype_counts)
print()

# 3. 统计特定六种细胞类型的总细胞数
target_types = ["CD4+ T", "γδ T", "CD8+ T", "DC", "Treg", "MKI67+ T"]
total_target_cells = 0

print("=== 六种特定细胞类型统计 ===")
for cell_type in target_types:
    if cell_type in celltype_counts.index:
        count = celltype_counts[cell_type]
        print(f"{cell_type}: {count} 个细胞")
        total_target_cells += count
    else:
        print(f"{cell_type}: 0 个细胞 (在数据中未找到)")

print(f"\n六种特定细胞类型总计: {total_target_cells} 个细胞")
print(f"数据中总细胞数: {adata.n_obs} 个细胞")
print(f"特定细胞占比: {total_target_cells/adata.n_obs*100:.2f}%")
