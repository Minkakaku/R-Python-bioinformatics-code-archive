import scanpy as sc
import pandas as pd

# 1. 读取h5ad文件
input_file = "adata_cluster.h5ad"  # 替换为您的文件路径
adata = sc.read_h5ad(input_file)

# 2. 检查必要的列是否存在
required_columns = ['cell_id', 'res_0.5']  # 根据您的实际列名调整

for col in required_columns:
    if col not in adata.obs.columns:
        raise ValueError(f"列 '{col}' 在adata.obs中不存在")

# 3. 提取所需数据
cell_ids = adata.obs['cell_id']
cluster_results = adata.obs['res_0.5']

# 4. 创建DataFrame
result_df = pd.DataFrame({
    'cell_id': cell_ids,
    'res_0.5': cluster_results
})

# 5. 保存结果
output_file = "cell_clusters.csv"  # 输出文件路径
result_df.to_csv(output_file, index=False)

print(f"成功提取并保存了 {len(result_df)} 个细胞的结果到 {output_file}")
print("文件格式:")
print(result_df.head())  # 显示前几行示例
