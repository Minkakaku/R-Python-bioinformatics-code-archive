import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import distance

# 1. 加载数据并确保名称唯一性
adata = sc.read_h5ad('integrated_data.h5ad')
adata.obs_names_make_unique()

# 2. 筛选目标细胞类型
target_celltypes = [
    "CD4+ TEM", "CD4+ Treg", "CD8+ Tex", "M2", "iDC",
    "POSTN+ FB", "THY1+ FB", "Kera Spinous", "Kera Inflammatory", "CTMCs"
]
adata = adata[adata.obs['celltype_final'].isin(target_celltypes)].copy()

# 3. 配置参数
tumor_cell_type = 'POSTN+ FB'
cell_type_col = 'celltype_final'
group_col = 'group'
distances = [10, 50, 100, 150, 200]  # 距离环定义
group_colors = {'disease': '#FF0000', 'normal': '#0000FF'}  # 红/蓝分组

# 4. 定位中心细胞
tumor_idx = adata.obs[adata.obs[cell_type_col] == tumor_cell_type].index[0]
center_coords = adata.obsm['spatial'][adata.obs.index == tumor_idx][0]

# 5. 计算距离并分区
distances_to_center = distance.cdist([center_coords], adata.obsm['spatial'])[0]
adata.obs['distance_to_center'] = distances_to_center

# 距离分区函数
def assign_zone(d):
    if d < 10: return "0-10"
    elif d < 50: return "10-50"
    elif d < 100: return "50-100"
    elif d < 150: return "100-150"
    elif d < 200: return "150-200"
    else: return "≥200"
adata.obs['distance_zone'] = adata.obs['distance_to_center'].apply(assign_zone)

# 6. 创建极坐标系
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(111, polar=True)
ax.set_theta_offset(np.pi/2)
ax.set_theta_direction(-1)
ax.set_ylim(0, 1.2)

# 7. 定义环形半径位置 [3,5](@ref)
ring_radii = {
    "0-10": 0.15, "10-50": 0.35, 
    "50-100": 0.55, "100-150": 0.75, 
    "150-200": 0.95, "≥200": 1.15
}

# 8. 绘制距离环标注线 [3](@ref)
for i, dist in enumerate(distances):
    r_pos = ring_radii[f"{dist}-{distances[i+1]}" if i<4 else "≥200"]
    ax.plot(np.linspace(0, 2*np.pi, 100), [r_pos]*100, 
            color='gray', linestyle='--', linewidth=0.7, alpha=0.7)
    ax.text(0, r_pos, f"{dist}μm", ha='left', va='center', 
            fontsize=10, bbox=dict(facecolor='white', alpha=0.8))

# 9. 设置细胞类型扇区 [7](@ref)
cell_types = adata.obs[cell_type_col].unique()
n_types = len(cell_types)
theta_step = 2 * np.pi / n_types

# 10. 细胞类型标签（最外围固定位置）[8](@ref)
label_radius = 1.25
for i, ct in enumerate(cell_types):
    angle = i * theta_step + theta_step/2
    ax.text(angle, label_radius, ct, ha='center', va='center', 
            fontsize=10, rotation=0,  # 取消旋转
            bbox=dict(boxstyle="round,pad=0.2", fc='white', ec='none', alpha=0.8))

# 11. 扇区分割线
for i in range(n_types + 1):
    ax.plot([i * theta_step, i * theta_step], [0, 1.15], 
            color='lightgray', linewidth=0.8)

# 12. 聚合数据（直接使用细胞数量）
agg_df = adata.obs.groupby(
    ['celltype_final', 'distance_zone', group_col]
).size().reset_index(name='cell_count')

# 13. 绘制分组气泡 [3,4](@ref)
max_count = agg_df['cell_count'].max()
for _, row in agg_df.iterrows():
    ct, zone, group = row['celltype_final'], row['distance_zone'], row[group_col]
    
    # 计算气泡位置
    angle_idx = list(cell_types).index(ct)
    angle = angle_idx * theta_step + theta_step/2
    r = ring_radii[zone]
    
    # 动态气泡大小（基于细胞数量）
    bubble_size = 50 + (row['cell_count'] / max_count) * 1500
    
    # 绘制分组气泡
    ax.scatter(
        angle, r, 
        s=bubble_size, 
        color=group_colors[group], 
        alpha=0.7,
        edgecolors='w',
        zorder=5
    )
    
    # 添加数量标签
    ax.text(angle, r, str(row['cell_count']), 
            ha='center', va='center', fontsize=8)

# 14. 中心标注
ax.text(0, 0, tumor_cell_type, ha='center', va='center', fontsize=14,
        bbox=dict(boxstyle="circle,pad=0.3", fc='white', ec='black', alpha=0.8))

# 15. 图例系统 [4](@ref)
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', label='Disease',
               markerfacecolor='#FF0000', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='Normal',
               markerfacecolor='#0000FF', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='Cell Count',
               markerfacecolor='grey', markersize=5)
]
ax.legend(handles=legend_elements, loc='lower right', 
          bbox_to_anchor=(1.15, 0.15), title="Group Legend")

# 16. 美化输出
ax.set_xticks([])
ax.set_yticks([])
ax.spines['polar'].set_visible(False)
plt.title('Spatial Distribution of Cell Types', fontsize=16, pad=20)
plt.tight_layout()
plt.savefig('cell_spatial_bubble.png', dpi=300, bbox_inches='tight')
plt.show()
