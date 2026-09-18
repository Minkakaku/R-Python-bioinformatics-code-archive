import scanpy as sc
import numpy as np
import pandas as pd

# 设置随机种子确保结果可复现
np.random.seed(42)

# 1. 读取h5ad文件
input_file = "celltype_TNK.h5ad"  # 替换为您的文件路径
adata = sc.read_h5ad(input_file)

# 2. 检查细胞总数
total_cells = adata.n_obs
print(f"原始数据细胞总数: {total_cells}")

# 3. 计算需要抽取的细胞数
sample_size = min(10000, total_cells)
print(f"将抽取 {sample_size} 个细胞")

# 4. 随机抽样并保留细胞ID
if total_cells > sample_size:
    # 创建随机索引
    random_indices = np.random.choice(total_cells, size=sample_size, replace=False)
    
    # 获取抽样细胞的原始ID
    sampled_cell_ids = adata.obs.index[random_indices].tolist()
    
    # 创建子集
    adata_sampled = adata[random_indices].copy()
else:
    print("警告：细胞总数不足10000，将使用全部细胞")
    adata_sampled = adata.copy()
    sampled_cell_ids = adata.obs.index.tolist()

# 5. 保存细胞ID到CSV文件
cell_id_df = pd.DataFrame({'cell_id': sampled_cell_ids})
cell_id_df.to_csv('sampled_cell_ids.csv', index=False)
print(f"抽样细胞ID已保存至: sampled_cell_ids.csv")

# 6. 验证结果
print(f"抽样后细胞数: {adata_sampled.n_obs}")

# 7. 保存抽样结果
output_file = "sampled_data.h5ad"  # 输出文件路径
adata_sampled.write_h5ad(output_file)
print(f"抽样数据已保存至: {output_file}")

# 8. 检查抽样是否保持原始分布
if 'res_0.5' in adata.obs.columns:
    print("\n原始分类比例:")
    print(adata.obs['res_0.5'].value_counts(normalize=True))
    print("\n抽样后分类比例:")
    print(adata_sampled.obs['res_0.5'].value_counts(normalize=True))
