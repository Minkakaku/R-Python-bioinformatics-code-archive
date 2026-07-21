import scanpy as sc
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull

# 读取数据
adata = sc.read_h5ad('All_filter_with_region_CNs10.h5ad')

# 提取 obs 中的 sample 和 CNs_10 列
obs = adata.obs[['sample', 'CNs_10']].copy()
obs['x'] = adata.obsm['spatial'][:, 0]
obs['y'] = adata.obsm['spatial'][:, 1]

# ==================== 1. 比例密度（细胞数占比） ====================
counts = pd.crosstab(obs['sample'], obs['CNs_10'])
sample_total = counts.sum(axis=1)
prop_density = counts.div(sample_total, axis=0)
prop_density.to_csv('CN_proportion_density_wide.csv')
print("已保存比例密度矩阵：CN_proportion_density_wide.csv")

# ==================== 2. 空间密度（每平方毫米细胞数） ====================
areas = {}
for sample, group in obs.groupby('sample'):
    coords = group[['x', 'y']].values
    if len(coords) >= 3:
        hull = ConvexHull(coords)
        area_um2 = hull.volume
    else:
        area_um2 = np.nan
    areas[sample] = area_um2

area_df = pd.DataFrame(list(areas.items()), columns=['sample', 'area_um2'])
total_df = sample_total.reset_index()
total_df.columns = ['sample', 'total_cells']
info = total_df.merge(area_df, on='sample', how='left')
info['area_mm2'] = info['area_um2'] / 1e6

counts_long = counts.reset_index().melt(id_vars='sample', var_name='CNs_10', value_name='count')
counts_long = counts_long.merge(info[['sample', 'area_mm2']], on='sample', how='left')
counts_long['spatial_density_per_mm2'] = counts_long['count'] / counts_long['area_mm2']
spatial_density_wide = counts_long.pivot(index='sample', columns='CNs_10', values='spatial_density_per_mm2')
spatial_density_wide.to_csv('CN_spatial_density_wide.csv')
print("已保存空间密度矩阵：CN_spatial_density_wide.csv")
