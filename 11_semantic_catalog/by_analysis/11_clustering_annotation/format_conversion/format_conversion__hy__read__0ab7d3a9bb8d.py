import anndata as ad

# 读取 h5ad 文件
adata = ad.read_h5ad("All_filter_state_samp_new_with_CNs10-15.h5ad")

# 获取 celltype 列的唯一值（去除 NaN）
unique_celltypes = adata.obs['celltype'].dropna().unique()

# 按字母顺序排序（可选，便于查看）
unique_celltypes = sorted(unique_celltypes)

# 保存到 CSV 文件（无行号、无表头）
with open("ct.csv", "w") as f:
    for ct in unique_celltypes:
        f.write(f"{ct}\n")

print(f"已保存 {len(unique_celltypes)} 个唯一细胞类型到 ct.csv")
