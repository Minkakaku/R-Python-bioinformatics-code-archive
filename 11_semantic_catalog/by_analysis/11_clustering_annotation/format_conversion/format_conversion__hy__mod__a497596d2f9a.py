import scanpy as sc

# 1. 读取h5ad文件，请将‘your_data.h5ad’替换为您的实际文件路径
adata = sc.read_h5ad('adata_cluster.h5ad')

# 2. 检查obs中是否包含所需的列
print("adata.obs 中的所有列：")
print(adata.obs.columns.tolist())

# 3. 输出sample和group列
# 确保列名与您的数据中完全一致，注意大小写
if 'sample' in adata.obs and 'group' in adata.obs:
    print("\n=== sample 和 group 列的内容 ===")
    # 直接打印这两列
    selected_columns = adata.obs[['sample', 'group']]
    print(selected_columns)
    
    # 如果需要获取唯一值，例如查看有哪些样本和组别
    print("\n=== 唯一的 sample ===")
    print(adata.obs['sample'].unique())
    print("\n=== 唯一的 group ===")
    print(adata.obs['group'].unique())
else:
    print("\n警告：在adata.obs中未找到名为‘sample’或‘group’的列。")
    print("请检查上一步输出的列名列表，并使用正确的列名。")

adata.obs[['sample', 'group']].to_csv('sample_group_info.csv')
print("信息已保存到 sample_group_info.csv")
