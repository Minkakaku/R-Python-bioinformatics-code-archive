import scanpy as sc

# 1. 读取原始数据文件
# 读取 combined_data.h5ad 以获取目标细胞ID列表
combined_adata = sc.read_h5ad('combined_data.h5ad')
# 读取 adata_counts.h5ad 以进行细胞筛选
adata_counts = sc.read_h5ad('adata_counts.h5ad')

# 2. 提取细胞ID
# 假设细胞ID存储在 combined_adata 的观察索引（obs.index）中
# 如果细胞ID在obs的某一列中（如 combined_adata.obs['cell_id']），请替换为对应的列名
target_cell_ids = combined_adata.obs.index.tolist()  # 获取所有细胞ID列表

# 3. 检查细胞ID在 adata_counts 中的存在情况
# 计算两个数据集的交集细胞ID
existing_cells = [cell_id for cell_id in target_cell_ids if cell_id in adata_counts.obs.index]
print(f"目标细胞数量: {len(target_cell_ids)}")
print(f"在 adata_counts 中匹配的细胞数量: {len(existing_cells)}")

# 4. 根据细胞ID筛选数据
# 使用布尔索引选择存在的细胞
adata_filtered = adata_counts[existing_cells, :]  # 筛选细胞（保留所有基因）
print(f"筛选后数据形状: {adata_filtered.n_obs} 个细胞 × {adata_filtered.n_vars} 个基因")

# 5. 保存筛选结果
adata_filtered.write('adata_filter.h5ad')  # 输出为H5AD文件
print("筛选数据已保存为 adata_filter.h5ad")
