import pandas as pd
import scanpy as sc

# 1. 读取 h5ad 文件
adata = sc.read_h5ad("All_filter_with_region.h5ad")

# 2. 读取制表符分隔的 CSV 文件（加 low_memory=False 避免 DtypeWarning）
df_cns = pd.read_csv("CNs_10_all_cell_features.csv", sep='\t', low_memory=False)

# 3. 提取需要的列，并根据 (cell_id, sample) 去重
df_cns_sub = df_cns[['cell_id', 'sample', 'CNs_10']].drop_duplicates(subset=['cell_id', 'sample'])

# 4. 为 adata.obs 添加临时索引
adata.obs['_temp_index'] = range(adata.obs.shape[0])

# 5. 左连接：将 CNs_10 信息合并到 obs
obs_merged = adata.obs.reset_index().merge(
    df_cns_sub,
    on=['cell_id', 'sample'],
    how='left'
)

# 6. 按临时索引排序，确保与 adata 顺序一致
obs_merged = obs_merged.sort_values('_temp_index')

# 7. 将 CNs_10 列直接赋给 adata.obs（不再使用 set_index）
adata.obs['CNs_10'] = obs_merged['CNs_10'].values

# 8. 删除临时列
adata.obs.drop(columns=['_temp_index'], inplace=True, errors='ignore')

# 检查匹配情况
matched_count = adata.obs['CNs_10'].notna().sum()
print(f"成功匹配 {matched_count} / {adata.n_obs} 个细胞")

# 可选：保存为新的 h5ad 文件
adata.write_h5ad("All_filter_with_region_CNs10.h5ad")
