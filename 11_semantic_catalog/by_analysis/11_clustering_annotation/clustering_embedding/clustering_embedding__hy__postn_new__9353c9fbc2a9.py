import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import os

# 加载数据
adata = sc.read_h5ad("modified_data_with_ct.h5ad")  # 替换为你的h5ad文件路径
print(adata.obs)
#celltypes = [
#    "CD4+ CD247+T", "CD4+ Tem", "CD4+ Treg", "CD8+ Tex", "Disease Stem",
#    "FAMS", "Kera Basal", "Kera Granular", "Kera Inflammatory", "Kera Spinous",
#    "LAMP3+ DC", "M2", "MARCO+ Mø", "Normal Stem", "SBC Proliferation",
#    "SBC Spinous", "SIGLEC1+ Mø", "SOCS3+ CTMCS", "SOCS3+ FB", "SOCS3+ M1",
#    "THY1+ FB", "cDC1", "cDC2", "iDC", "CTMCS"
#]
#
## 创建ct列并初始化为"others"
#adata.obs['ct'] = "others"
#
## 遍历每个细胞类型，将目标细胞类型的ct设为原名称
#for ct in target_celltypes:
#    # 检查当前细胞类型是否存在于数据中
#    if ct in adata.obs['ct'].unique():
#        # 创建布尔掩码标识该类型的所有细胞
#        mask = adata.obs['ct'] == ct
#        # 只修改存在此类型的细胞
#        adata.obs.loc[mask, 'ct'] = ct
#    else:
#        print(f"警告: 细胞类型 '{ct}' 在数据中未找到，跳过处理")
#ID到obs
#adata.obs['cell_id'] = adata.obs.index

# 提取POSTN+ FB细胞
postn_fb_cells = adata[adata.obs['ct'] == 'POSTN+ FB'].copy()
other_cells = adata[adata.obs['ct'] != 'POSTN+ FB'].copy()

# 检查空间坐标
if 'spatial' not in adata.obsm:
    raise ValueError("Spatial coordinates not found in .obsm['spatial']")
    
# 计算最近邻居
def find_nearest_neighbor(postn_cell, other_cells):
    postn_coord = postn_cell.obsm['spatial']
    other_coords = other_cells.obsm['spatial']
    
    # 计算所有距离
    dists = distance.cdist(postn_coord, other_coords, 'euclidean')[0]
    
    # 找到最近的非POSTN+FB细胞
    nearest_idx = np.argmin(dists)
    nearest_cell_id = other_cells.obs.index[nearest_idx]
    nearest_cell_type = other_cells.obs['ct'].iloc[nearest_idx]
    
    return nearest_cell_id, nearest_cell_type

# 应用函数到每个POSTN+ FB细胞
results = []
for i, cell_id in enumerate(postn_fb_cells.obs.index):
    if i % 100 == 0:
        print(f"Processing cell {i+1}/{len(postn_fb_cells)}")
    
    cell_data = postn_fb_cells[postn_fb_cells.obs.index == cell_id]
    group = cell_data.obs['group'].values[0]
    
    nearest_id, nearest_type = find_nearest_neighbor(cell_data, other_cells)
    
    results.append({
        'postn_cell_id': cell_id,
        'nearest_cell_id': nearest_id,
        'nearest_cell_type': nearest_type,
        'group': group
    })

# 创建结果DataFrame
results_df = pd.DataFrame(results)

# 文件1: 最近邻细胞列表
nearest_neighbor_df = results_df[['nearest_cell_id', 'nearest_cell_type']]
nearest_neighbor_df.to_csv("nearest_neighbors.csv", index=False)

# 文件2: 统计结果
# 按组和细胞类型计数
count_df = results_df.groupby(['group', 'nearest_cell_type']).size().reset_index(name='count')

# 计算比例
total_counts = results_df.groupby('group').size()
count_df['proportion'] = count_df.apply(
    lambda row: row['count'] / total_counts[row['group']], axis=1
)
count_df.to_csv("nearest_neighbor_stats.csv", index=False)

# 准备绘图数据
pivot_df = count_df.pivot(index='nearest_cell_type', columns='group', values='proportion').fillna(0)
pivot_df = pivot_df.sort_index(axis=1)  # 确保顺序一致

# 创建绘图
plt.figure(figsize=(12, 8))
ax = plt.gca()

# 颜色设置
unique_types = pivot_df.index
cmap = plt.cm.get_cmap('tab20', len(unique_types))
type_colors = {t: cmap(i) for i, t in enumerate(unique_types)}
light_colors = {t: (*cmap(i)[:3], 0.4) for i, t in enumerate(unique_types)}  # 稍淡颜色

# 分组位置
groups = sorted(results_df['group'].unique())
bar_width = 0.35
x_pos = np.arange(len(groups))

# 绘制堆积柱状图
bottom = np.zeros(len(groups))
bars = []

for i, cell_type in enumerate(unique_types):
    proportions = pivot_df.loc[cell_type].values
    bar = ax.bar(
        x_pos, 
        proportions, 
        bottom=bottom, 
        width=bar_width,
        color=type_colors[cell_type],
        label=cell_type,
        edgecolor='white'
    )
    bars.append(bar)
    bottom += proportions

# 添加梯形连接
trapezoids = []
for i, cell_type in enumerate(unique_types):
    # 获取当前类型在两个组中的比例和累积高度
    disease_prop = pivot_df.loc[cell_type, 'disease']
    normal_prop = pivot_df.loc[cell_type, 'normal']
    
    disease_index = groups.index('disease')
    normal_index = groups.index('normal')
    
    # 计算累积高度
    disease_height = 0
    normal_height = 0
    for t in unique_types:
        if t == cell_type:
            break
        disease_height += pivot_df.loc[t, 'disease']
        normal_height += pivot_df.loc[t, 'normal']
    
    # 计算梯形顶点坐标 (按组的位置)
    x_disease = x_pos[disease_index] + bar_width/2
    x_normal = x_pos[normal_index] - bar_width/2
    
    y_disease_bottom = disease_height
    y_disease_top = y_disease_bottom + disease_prop
    y_normal_bottom = normal_height
    y_normal_top = y_normal_bottom + normal_prop
    
    # 创建梯形顶点
    vertices = np.array([
        [x_disease, y_disease_bottom],  # 疾病组右下
        [x_disease, y_disease_top],    # 疾病组右上
        [x_normal, y_normal_top],      # 正常组左上
        [x_normal, y_normal_bottom]    # 正常组左下
    ])
    
    trapezoid = patches.Polygon(vertices, closed=True)
    trapezoids.append(trapezoid)

# 添加梯形集合
pc = PatchCollection(
    trapezoids, 
    facecolor=[light_colors[t] for t in unique_types],
    edgecolor='none',
    alpha=0.6
)
ax.add_collection(pc)

# 完善图形
plt.xticks(x_pos, groups, fontsize=12)
plt.ylabel('Proportion', fontsize=14)
plt.title('Nearest Neighbor Cell Type Composition by Group', fontsize=16)
plt.grid(axis='y', linestyle='--', alpha=0.7)

# 图例处理
handles, labels = ax.get_legend_handles_labels()
plt.legend(
    handles[:len(unique_types)], 
    labels[:len(unique_types)], 
    title='Cell Types',
    bbox_to_anchor=(1.05, 1), 
    loc='upper left'
)

plt.tight_layout()
plt.savefig('neighbor_composition_trend.png', dpi=300, bbox_inches='tight')
plt.show()
