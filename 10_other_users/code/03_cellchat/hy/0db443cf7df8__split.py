import scanpy as sc
import anndata as ad

# 输入文件路径
input_file = "adata_with_band_by_sample.h5ad"
output_file = "adata_with_band_by_sample_core_NR.h5ad"

# 读取数据
adata = sc.read_h5ad(input_file)
print(f"原始数据：{adata.n_obs} 个细胞，{adata.n_vars} 个基因")

# 检查所需列是否存在
required_cols = ['region', 'group_final']
for col in required_cols:
    if col not in adata.obs.columns:
        raise ValueError(f"obs 中缺少列 '{col}'，请检查数据")

# 创建布尔掩码：同时满足两个条件
mask = (adata.obs['region'] == 'tumor') & (adata.obs['group_final'].isin(['POST-NR', 'PRE-NR']))

# 提取子集（使用 .copy() 避免后续修改影响原对象）
adata_subset = adata[mask].copy()
print(f"筛选后：{adata_subset.n_obs} 个细胞，{adata_subset.n_vars} 个基因")

# 保存结果
adata_subset.write_h5ad(output_file)
print(f"已保存至：{output_file}")
