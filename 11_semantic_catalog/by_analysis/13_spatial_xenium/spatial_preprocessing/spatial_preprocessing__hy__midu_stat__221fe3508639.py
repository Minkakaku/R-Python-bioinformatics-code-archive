import scanpy as sc
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull

# ==================== 读取数据 ====================
adata = sc.read_h5ad('All_filter_with_region_CNs10.h5ad')

# ==================== 指定要保留的样本列表 ====================
samples_to_keep = [
    'PRE_10-1', 'PRE_10-2', 'PRE_15-1', 'PRE_15-3', 'PRE_15-4', 'PRE_16',
    'PRE_17-1', 'PRE_17-3', 'PRE_18-1', 'PRE_18-2', 'PRE_19', 'PRE_20-1',
    'PRE_20-2', 'PRE_20-3', 'PRE_20-4', 'PRE_21', 'PRE_22-1', 'PRE_22-2',
    'PRE_23-1', 'PRE_23-2', 'PRE_25-1', 'PRE_25-2', 'PRE_25-3', 'PRE_25-4',
    'PRE_26-1', 'PRE_26-3', 'PRE_27-1', 'PRE_27-2', 'PRE_3-1', 'PRE_3-2',
    'PRE_5-1', 'PRE_5-2', 'PRE_6-1', 'PRE_6-2', 'PRE_6-3', 'PRE_8'
]

adata = adata[adata.obs['sample'].isin(samples_to_keep)].copy()
print(f"筛选后剩余细胞数: {adata.shape[0]}")

# ==================== 提取必要信息并转换类型 ====================
obs = adata.obs[['sample', 'celltype', 'region']].copy()
obs['x'] = adata.obsm['spatial'][:, 0]
obs['y'] = adata.obsm['spatial'][:, 1]
# 转为字符串，避免后续拼接错误
obs['celltype'] = obs['celltype'].astype(str)
obs['region'] = obs['region'].astype(str)

# ==================== 按样本 + region 分组，计算每个组的凸包面积 ====================
area_data = []
for (sample, region), group in obs.groupby(['sample', 'region'], observed=True):
    coords = group[['x', 'y']].values
    if len(coords) >= 3:
        hull = ConvexHull(coords)
        area_um2 = hull.volume   # 凸包面积（平方微米）
    else:
        area_um2 = np.nan
    area_data.append([sample, region, area_um2])
area_df = pd.DataFrame(area_data, columns=['sample', 'region', 'area_um2'])
area_df['area_mm2'] = area_df['area_um2'] / 1e6   # 转换为平方毫米

# ==================== 统计每个样本+region+celltype 的计数 ====================
counts = obs.groupby(['sample', 'region', 'celltype'], observed=True).size().reset_index(name='count')

# 计算每个样本+region 的总细胞数
region_total = obs.groupby(['sample', 'region'], observed=True).size().reset_index(name='region_total_cells')

# 合并
counts = counts.merge(region_total, on=['sample', 'region'])
counts = counts.merge(area_df, on=['sample', 'region'], how='left')

# 比例密度 = 计数 / 该 region 总细胞数
counts['proportion_density'] = counts['count'] / counts['region_total_cells']
# 空间密度 = 计数 / 面积_mm2
counts['spatial_density_per_mm2'] = counts['count'] / counts['area_mm2']

# 创建组合列名： celltype_region
counts['celltype_region'] = counts['celltype'] + '_' + counts['region']

# ==================== 转换为宽格式并保存 ====================
# 比例密度宽表
prop_wide = counts.pivot(index='sample', columns='celltype_region', values='proportion_density')
prop_wide = prop_wide.fillna(0)
prop_wide.to_csv('region_proportion_density_wide.csv')
print("已保存比例密度矩阵（每个 region 内细胞类型占比）：region_proportion_density_wide.csv")

# 空间密度宽表
spatial_wide = counts.pivot(index='sample', columns='celltype_region', values='spatial_density_per_mm2')
spatial_wide = spatial_wide.fillna(0)
spatial_wide.to_csv('region_spatial_density_wide.csv')
print("已保存空间密度矩阵（每个 region 内细胞类型每平方毫米个数）：region_spatial_density_wide.csv")
