import scanpy as sc
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull

# 读取数据
adata = sc.read_h5ad('All_filter_with_region_CNs10.h5ad')

# 提取 obs 和 spatial 坐标
obs = adata.obs[['sample', 'celltype']].copy()
obs['x'] = adata.obsm['spatial'][:, 0]
obs['y'] = adata.obsm['spatial'][:, 1]

# ==================== 1. 比例密度（细胞数占比） ====================
# 交叉表：样本 × 细胞类型 的计数
counts = pd.crosstab(obs['sample'], obs['celltype'])
# 每个样本的总细胞数
sample_total = counts.sum(axis=1)
# 比例密度 = 计数 / 总细胞数
prop_density = counts.div(sample_total, axis=0)
# 保存
prop_density.to_csv('ct_proportion_density_wide.csv')
print("已保存比例密度矩阵：ct_proportion_density_wide.csv")

# ==================== 2. 空间密度（每平方毫米细胞数） ====================
# 计算每个样本的凸包面积（平方微米）
areas = {}
for sample, group in obs.groupby('sample'):
    coords = group[['x', 'y']].values
    if len(coords) >= 3:
        hull = ConvexHull(coords)
        area_um2 = hull.volume  # 凸包面积（二维面积）
    else:
        area_um2 = np.nan
    areas[sample] = area_um2

# 转换为 DataFrame
area_df = pd.DataFrame(list(areas.items()), columns=['sample', 'area_um2'])

# 每个样本的总细胞数（从 counts 获取）
total_df = sample_total.reset_index()
total_df.columns = ['sample', 'total_cells']

# 合并面积信息
info = total_df.merge(area_df, on='sample', how='left')
# 计算面积（平方毫米）
info['area_mm2'] = info['area_um2'] / 1e6

# 计算每种细胞类型的空间密度（细胞数 / 平方毫米）
# 先将 counts 转换为长格式
counts_long = counts.reset_index().melt(id_vars='sample', var_name='celltype', value_name='count')
# 合并面积
counts_long = counts_long.merge(info[['sample', 'area_mm2']], on='sample', how='left')
# 计算密度
counts_long['spatial_density_per_mm2'] = counts_long['count'] / counts_long['area_mm2']
# 转换回宽格式
spatial_density_wide = counts_long.pivot(index='sample', columns='celltype', values='spatial_density_per_mm2')
# 保存
spatial_density_wide.to_csv('ct_spatial_density_wide.csv')
print("已保存空间密度矩阵：ct_spatial_density_wide.csv")
