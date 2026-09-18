import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, Wedge
import pandas as pd
from scipy.spatial import distance

# 加载数据
adata = sc.read_h5ad('integrated_data.h5ad')
target_celltypes = {"CD4+ TEM", "CD4+ Treg", "CD8+ Tex", "M2", "CTMCs", 
                    "iDC", "POSTN+ FB", "THY1+ FB", "Kera Spinous", "Kera Inflammatory"}
cell_mask = adata.obs['celltype_final'].isin(target_celltypes)
adata = adata[cell_mask].copy()

# 配置参数
tumor_cell_type = 'POSTN+ FB'
cell_type_col = 'celltype_final'
group_col = 'group'
distances = [10, 50, 100, 150, 200]  # 距离环
group_colors = {'disease': '#FF0000', 'normal': '#0000FF'}  # 红蓝分组色

# 1. 定位中心细胞
tumor_idx = adata.obs[adata.obs[cell_type_col] == tumor_cell_type].index[0]
center_coords = adata.obsm['spatial'][adata.obs.index == tumor_idx][0]

# 2. 计算所有细胞到中心的距离
distances_to_center = distance.cdist([center_coords], adata.obsm['spatial'])[0]
adata.obs['distance_to_center'] = distances_to_center

# 3. 定义距离区域
distance_zones = [f"<{distances[0]}"] + [f"{distances[i-1]}-{distances[i]}" for i in range(1, len(distances))] + [f">={distances[-1]}"]
adata.obs['distance_zone'] = adata.obs['distance_to_center'].apply(
    lambda d: next((zone for i, zone in enumerate(distance_zones) 
                  if (i==0 and d<distances[0]) or (i==len(distances) and d>=distances[-1]) or (distances[i-1] <= d < distances[i])), distance_zones[-1])
)

# 4. 过滤中心细胞
adata = adata[adata.obs.index != tumor_idx].copy()

# 5. 创建极坐标系
plt.figure(figsize=(15, 15))
ax = plt.subplot(111, polar=True)

# 6. 设置环形参数
num_zones = len(distance_zones)
zone_radii = np.linspace(0.1, 0.9, num_zones)
max_radius = 1.0  # 最外层半径

# 7. 获取细胞类型和分组
cell_types = sorted(adata.obs[cell_type_col].unique())
groups = sorted(adata.obs[group_col].unique())
theta_step = 2 * np.pi / len(cell_types)
cell_type_angles = {ct: i * theta_step for i, ct in enumerate(cell_types)}

# 8. 绘制距离环的网格线（灰色虚线）
for dist in distances:
    r_pos = 0.1 + (dist / distances[-1]) * 0.8
    ax.plot(np.linspace(0, 2*np.pi, 100), [r_pos]*100, 
            color='gray', linestyle='--', alpha=0.5)

# 9. 添加细胞类型分割线（黑色实线）[1,3](@ref)
for i, ct in enumerate(cell_types):
    angle = cell_type_angles[ct]
    # 从中心到最外层绘制分割线
    ax.plot([angle, angle], [0, max_radius], color='black', linewidth=1.5, alpha=0.7)

# 10. 中心位置添加细胞类型标签 [1,3](@ref)
ax.text(0, 0, tumor_cell_type, 
        ha='center', va='center', fontsize=14, weight='bold',
        bbox=dict(boxstyle="circle,pad=0.3", fc='white', ec='black', alpha=0.8))

# 11. 最外层添加细胞类型标签（200um外）[1,3](@ref)
label_radius = max_radius + 0.15  # 200um环外的位置
for ct, angle in cell_type_angles.items():
    if ct == tumor_cell_type: continue  # 跳过中心细胞
    label_angle = angle + theta_step/2  # 标签位于扇形中心
    x = label_radius * np.cos(label_angle)
    y = label_radius * np.sin(label_angle)
    ax.text(x, y, ct, ha='center', va='center', fontsize=10,
            rotation=np.degrees(label_angle) - 90,  # 径向旋转
            bbox=dict(boxstyle="round,pad=0.2", fc='white', ec='none', alpha=0.7))

# 12. 聚合数据（按类型/区域/分组）
agg_df = adata.obs.groupby(['celltype_final', 'distance_zone', group_col]).size().reset_index(name='count')

# 13. 绘制分组点（红/蓝）[1,3](@ref)
size_factor = 50  # 点大小缩放因子
for _, row in agg_df.iterrows():
    ct, zone, group = row['celltype_final'], row['distance_zone'], row[group_col]
    if ct == tumor_cell_type: continue
    
    zone_idx = distance_zones.index(zone)
    r = zone_radii[zone_idx]
    angle_center = cell_type_angles[ct] + theta_step/2
    
    # 分组偏移量（避免重叠）
    group_offset = 0.1 if group == 'disease' else -0.1
    
    # 绘制点（大小表示细胞数量）
    s = max(10, min(row['count'] * size_factor, 300))
    ax.scatter(angle_center + group_offset, r, s=s,
               color=group_colors[group], alpha=0.8, 
               edgecolors='w', linewidth=0.5)

# 14. 右下角添加图例 [9,10,11](@ref)
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', label='Disease',
               markerfacecolor='#FF0000', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='Normal',
               markerfacecolor='#0000FF', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='1 cell',
               markerfacecolor='gray', markersize=5),
    plt.Line2D([0], [0], marker='o', color='w', label='5 cells',
               markerfacecolor='gray', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='10 cells',
               markerfacecolor='gray', markersize=15)
]

ax.legend(handles=legend_elements, title='图例', 
          loc='lower right', bbox_to_anchor=(1.25, -0.05), 
          frameon=True, framealpha=0.9, ncol=2)

# 15. 美化图形
ax.set_theta_offset(np.pi/2)  # 0°位于顶部
ax.set_theta_direction(-1)     # 顺时针方向
ax.set_rlim(0, max_radius + 0.2)
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])

plt.tight_layout()
plt.savefig('空间分布分析图.png', dpi=300, bbox_inches='tight')
plt.show()
