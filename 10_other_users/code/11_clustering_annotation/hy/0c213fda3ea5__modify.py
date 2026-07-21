import scanpy as sc

# 读取h5ad文件
adata = sc.read_h5ad('merged_single_cell_data.h5ad')  # 请替换为您的文件路径

# 检查 celltype 列是否存在
if 'celltype' not in adata.obs.columns:
    raise KeyError("在 adata.obs 中未找到 'celltype' 列。可用的列有: {}".format(list(adata.obs.columns)))

# 显示修改前的唯一细胞类型
print("修改前的唯一细胞类型:")
print(adata.obs['celltype'].value_counts())

# 使用 replace 方法进行替换
# 此方法会确保只替换完全匹配 'Plasma' 的条目，不会影响包含 'Plasma' 的其他字符串（如 'Plasma Cell'）
adata.obs['celltype'] = adata.obs['celltype'].replace('Plasma', 'Plasmas')

# 显示修改后的唯一细胞类型，确认更改
print("\n修改后的唯一细胞类型:")
print(adata.obs['celltype'].value_counts())

# 可选：保存修改后的数据到新文件，避免覆盖原文件
adata.write('merged_single_cell_data_modified.h5ad')
