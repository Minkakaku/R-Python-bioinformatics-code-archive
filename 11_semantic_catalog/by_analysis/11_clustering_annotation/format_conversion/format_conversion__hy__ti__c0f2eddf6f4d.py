import scanpy as sc

# 读取h5ad文件
adata = sc.read_h5ad('merged_celltype_stage.h5ad')  # 请替换为您的h5ad文件路径

# 查看adata.obs中的列名，确认包含'stage'列
print("adata.obs的列名:", list(adata.obs.columns))

# 查看stage列中有哪些类别
if 'stage' in adata.obs.columns:
    print("stage列中的类别:", adata.obs['stage'].unique())
else:
    print("警告: 在adata.obs中未找到'stage'列。")
    # 打印所有列名以便用户检查
    print("可用的列名:", list(adata.obs.columns))

# 提取stage为'Adenoma'的数据
# 确保类别名称完全匹配（包括大小写）
adata_adenoma = adata[adata.obs['stage'] == 'Tumor', :]

# 输出原始数据和提取后数据的信息
print("\n原始数据形状 (细胞数, 基因数):", adata.shape)
print("Adenoma数据形状 (细胞数, 基因数):", adata_adenoma.shape)

# 可选：保存提取后的数据到新的h5ad文件
adata_adenoma.write('tumor_data.h5ad')

print("提取完成！")
