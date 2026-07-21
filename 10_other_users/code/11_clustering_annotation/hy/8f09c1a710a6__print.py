import anndata

# 1. 读取文件
adata = anndata.read_h5ad("merged_single_cell_data_with_group_new.h5ad")

# 2. 提取 sample 和 group_new 列
df = adata.obs[['sample', 'group_new']]

# 3. 打印到终端
print("=== sample 和 group_new 列 ===")
print(df.to_string())  # 打印全部行，若数据量大可改用 .head(20)

# 4. 保存为 CSV 文件
df.to_csv("sample_group_new.csv", index=False)
print("\n已保存至 sample_group_new.csv")
