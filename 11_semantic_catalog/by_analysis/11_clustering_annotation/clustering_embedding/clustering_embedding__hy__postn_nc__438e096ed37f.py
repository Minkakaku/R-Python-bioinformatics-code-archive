import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# 1. 加载数据并确保obs唯一
adata = sc.read_h5ad('integrated_data.h5ad')
adata.obs_names_make_unique()  # 修复细胞名重复警告

# 2. 提取POSTN+ FB细胞
postn_cells = adata[adata.obs['celltype_final'] == 'POSTN+ FB'].copy()

# 3. 计算最近邻
nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree').fit(adata.obsm['spatial'])
_, indices = nbrs.kneighbors(postn_cells.obsm['spatial'])
postn_cells.obs['nearest_type'] = adata.obs.iloc[indices[:, 1]]['celltype_final'].values

# 4. 修正的统计函数
def calculate_proportions(adata_subset):
    """计算各细胞类型占比，返回带索引的Series"""
    counts = adata_subset.obs['nearest_type'].value_counts(normalize=True) * 100
    # 确保包含所有细胞类型
    all_types = adata.obs['celltype_final'].unique()
    return counts.reindex(all_types, fill_value=0)

# 5. 分组统计
disease_subset = postn_cells[postn_cells.obs['group'] == 'disease'].copy()
normal_subset = postn_cells[postn_cells.obs['group'] == 'normal'].copy()

prop_disease = calculate_proportions(disease_subset)
prop_normal = calculate_proportions(normal_subset)  # 修正变量名

# 6. 准备绘图数据
cell_types = prop_disease.index.tolist()
colors = plt.cm.tab20(np.linspace(0, 1, len(cell_types)))

# 7. 创建堆叠柱状图
fig, ax = plt.subplots(figsize=(12, 8))
ax.set_facecolor('#f5f5f5')
plt.grid(axis='y', linestyle='--', alpha=0.7)

bar_width = 0.4
x_pos = {'disease': 0, 'normal': 1}
cum_disease, cum_normal = 0, 0
patches = []

for i, cell_type in enumerate(cell_types):
    # 获取当前类型高度
    h_disease = prop_disease.loc[cell_type]
    h_normal = prop_normal.loc[cell_type]  # 使用loc按索引取值
    
    # 绘制分组柱状图
    ax.bar(x_pos['disease'], h_disease, bottom=cum_disease, 
           width=bar_width, color=colors[i], edgecolor='white', 
           label=cell_type if i == 0 else "")
    ax.bar(x_pos['normal'], h_normal, bottom=cum_normal, 
           width=bar_width, color=colors[i], edgecolor='white')
    
    # 创建连接线多边形
    polygon = Polygon([
        [x_pos['disease'] + bar_width/2, cum_disease],
        [x_pos['disease'] + bar_width/2, cum_disease + h_disease],
        [x_pos['normal'] - bar_width/2, cum_normal + h_normal],
        [x_pos['normal'] - bar_width/2, cum_normal]
    ], closed=True, alpha=0.3, color=colors[i])
    patches.append(polygon)
    
    # 更新累积高度
    cum_disease += h_disease
    cum_normal += h_normal

# 8. 添加图形元素
ax.add_collection(PatchCollection(patches, match_original=True))
ax.set_xticks([x_pos['disease'], x_pos['normal']])
ax.set_xticklabels(['Disease Group', 'Normal Group'], fontsize=12)
ax.set_ylabel('Nearest Neighbor Proportion (%)', fontsize=12)
ax.set_title('POSTN+ FB Cell Nearest Neighbor Distribution', fontsize=14)
ax.legend(loc='upper right', bbox_to_anchor=(1.15, 1), title='Cell Types')

# 9. 添加统计标注
mean_diff = (prop_disease - prop_normal).mean()
plt.figtext(0.5, 0.01, 
            f"Total POSTN+ FB cells: {postn_cells.n_obs} | Mean Difference: {mean_diff:.1f}%", 
            ha='center', fontsize=10, style='italic')

plt.tight_layout()
plt.savefig('postn_fb_neighbor_proportions.pdf', dpi=300)
plt.show()
