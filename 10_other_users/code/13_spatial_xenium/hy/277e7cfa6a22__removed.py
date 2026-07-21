import scanpy as sc

# 读取h5ad文件（请将路径替换为您的实际文件路径）
adata = sc.read_h5ad('adata_celltype_RCTD.h5ad')  # 替换为您的h5ad文件路径

# 检查annotation列是否存在
if 'annotation' not in adata.obs.columns:
    raise KeyError("在adata.obs中未找到'annotation'列。可用的列有: " + ", ".join(adata.obs.columns))

# 定义要删除的细胞类型列表
remove_celltypes = [
    "Endo-CAF",
    "B-T", 
    "Unknown",
    "Plasma-CAF",
    "B-CAF",
    "T-Mast",
    "Endo-Immu", 
    "B-Epi",
    "T-Epi",
    "T-CAF"
]

# 显示删除前的细胞统计
print("删除前的细胞类型分布:")
print(adata.obs['annotation'].value_counts())
print(f"删除前总细胞数: {adata.n_obs}")

# 创建布尔掩码：选择annotation列不在remove_celltypes中的细胞
keep_mask = ~adata.obs['annotation'].isin(remove_celltypes)

# 应用过滤
adata_filtered = adata[keep_mask, :]

# 显示删除后的细胞统计
print("\n删除后的细胞类型分布:")
print(adata_filtered.obs['annotation'].value_counts())
print(f"删除后总细胞数: {adata_filtered.n_obs}")
print(f"删除了 {adata.n_obs - adata_filtered.n_obs} 个细胞")

# 保存过滤后的数据
output_filename = "adata_celltype_RCTD_removed.h5ad"
adata_filtered.write(output_filename)
print(f"\n过滤后的数据已保存为: {output_filename}")
