import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import distance

# 1. 加载数据并确保名称唯一性
adata = sc.read_h5ad('integrated_data.h5ad')
adata.obs_names_make_unique()

# 2. 筛选目标细胞类型
target_celltypes = {
    "CD4+ TEM", "CD4+ Treg", "CD8+ Tex", "M2", "CTMCs", 
    "iDC", "POSTN+ FB", "THY1+ FB", "Kera Spinous", "Kera Inflammatory"
}
cell_mask = adata.obs['celltype_final'].isin(target_celltypes)
adata = adata[cell_mask].copy()

# 3. 配置参数
tumor_cell_type = 'POSTN+ FB'
cell_type_col = 'celltype_final'
group_col = 'group'
distances = [10, 50, 100, 150, 200]
group_colors = {'disease': '#FF0000', 'normal': '#0000FF'}

# 4. 定位中心细胞
tumor_idx = adata.obs[adata.obs[cell_type_col] == tumor_cell_type].index[0]
center_coords = adata.obsm['spatial'][adata.obs.index == tumor_idx][0]

# 5. 计算所有细胞到中心的距离
distances_to_center = distance.cdist([center_coords], adata.obsm['spatial'])[0]
adata.obs['distance_to_center'] = distances_to_center

# 6. 定义距离区域
distance_zones = [f"<{distances[0]}"] 
for i in range(1, len(distances)):
    distance_zones.append(f"{distances[i-1]}-{distances[i]}")
distance_zones.append(f">={distances[-1]}")

def assign_distance_zone(d):
    if d < distances[0]:
        return distance_zones[0]
    for i in range(1, len(distances)):
        if distances[i-1] <= d < distances[i]:
            return distance_zones[i]
    return distance_zones[-1]

adata.obs['distance_zone'] = adata.obs['distance_to_center'].apply(assign_distance_zone)

# 7. 过滤中心细胞
adata = adata[adata.obs.index != tumor_idx].copy()

# 8. 创建极坐标系
plt.figure(figsize=(15, 15))
ax = plt.subplot(111, polar=True)

# 9. 设置环形参数
max_radius = 1.0
zone_radii = [0.2, 0.35, 0.5, 0.65, 0.85]  # 手工调整距离环位置

# 10. 获取所有细胞类型（排除中心细胞）
all_cell_types = sorted(target_celltypes - {tumor_cell_type})
theta_step = 2 * np.pi / len(all_cell_types)

# 11. 固定细胞类型标签位置（按图片样式）
fixed_label_positions = {
    "CD8+ Tex": (1.15 * np.cos(0.2 * theta_step), 1.15 * np.sin(0.2 * theta_step), 0),
    "CD4+ Treg": (1.15 * np.cos(1.0 * theta_step), 1.15 * np.sin(1.0 * theta_step), 0),
    "CD4+ TEM": (1.15 * np.cos(1.8 * theta_step), 1.15 * np.sin(1.8 * theta_step), 0),
    "CTMCs": (1.15 * np.cos(2.6 * theta_step), 1.15 * np.sin(2.6 * theta_step), 0),
    "Kera Inflammatory": (1.15 * np.cos(3.4 * theta_step), 1.15 * np.sin(3.4 * theta_step), 0),
    "Kera Spinous": (1.15 * np.cos(4.2 * theta_step), 1.15 * np.sin(4.2 * theta_step), 0),
    "M2": (1.15 * np.cos(5.0 * theta_step), 1.15 * np.sin(5.0 * theta_step), 0),
    "iDC": (1.15 * np.cos(5.8 * theta_step), 1.15 * np.sin(5.8 * theta_step), 0),
    "THY1+ FB": (1.15 * np.cos(6.6 * theta_step), 1.15 * np.sin(6.6 * theta_step), 0),
}

# 12. 绘制所有固定位置的细胞类型标签
for ct, (x, y, rotation) in fixed_label_positions.items():
    # 计算智能旋转角度（避免文字倒置）
    degrees_angle = np.degrees(np.arctan2(y, x))
    if degrees_angle > 90 and degrees_angle < 270:
        rotation = degrees_angle - 180
    else:
        rotation = degrees_angle
    
    ax.text(
        x, y, ct, 
        ha='center', va='center', 
        fontsize=10, fontweight='bold',
        rotation=rotation,
        bbox=dict(boxstyle="round,pad=0.2", fc='white', ec='none', alpha=0.8)
    )

# 13. 绘制距离环并标注
for i, dist in enumerate(distances):
    # 在右侧标注距离
    ax.text(0, zone_radii[i], f"{dist}μm", 
            ha='left', va='center', fontsize=10, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.8))

# 14. 添加细胞类型分割线
for i in range(len(all_cell_types)):
    angle = i * theta_step
    ax.plot([angle, angle], [0, max_radius], color='black', linewidth=1, alpha=0.5)

# 15. 中心细胞标签
ax.text(0, 0, tumor_cell_type, 
        ha='center', va='center', fontsize=14, fontweight='bold',
        bbox=dict(boxstyle="circle,pad=0.3", fc='white', ec='black', alpha=0.8))

# 16. 聚合数据并绘制分组点
agg_df = adata.obs.groupby(['celltype_final', 'distance_zone', group_col]).size().reset_index(name='count')

# 根据图片排列手动设置细胞类型角度映射
cell_type_angles = {
    "CD8+ Tex": 0.2 * theta_step,
    "CD4+ Treg": 1.0 * theta_step,
    "CD4+ TEM": 1.8 * theta_step,
    "CTMCs": 2.6 * theta_step,
    "Kera Inflammatory": 3.4 * theta_step,
    "Kera Spinous": 4.2 * theta_step,
    "M2": 5.0 * theta_step,
    "iDC": 5.8 * theta_step,
    "THY1+ FB": 6.6 * theta_step
}

size_factor = 50
for _, row in agg_df.iterrows():
    ct, zone, group = row['celltype_final'], row['distance_zone'], row[group_col]
    
    # 跳过中心细胞类型
    if ct == tumor_cell_type: continue
    
    # 手动映射区域和半径
    zone_radius_map = {
        "<10": zone_radii[0],
        "10-50": zone_radii[1],
        "50-100": zone_radii[2],
        "100-150": zone_radii[3],
        "150-200": zone_radii[4],
        ">=200": 0.95
    }
    
    r = zone_radius_map.get(zone, 0.5)
    
    # 分组偏移
    group_offset = 0.1 if group == 'disease' else -0.1
    
    # 点大小（根据图片显示限制最大尺寸）
    s = min(max(20, row['count'] * size_factor), 300)
    
    # 绘制点
    ax.scatter(
        cell_type_angles[ct] + group_offset, 
        r, 
        s=s,
        color=group_colors[group], 
        alpha=0.8, 
        edgecolors='w', 
        linewidth=0.5
    )

# 17. 右下角图例（按图片样式）
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

# 精确匹配图片图例位置
ax.legend(
    handles=legend_elements, 
    title='Legend', 
    loc='lower right', 
    bbox_to_anchor=(1.15, 0.25),  # 匹配图片位置
    frameon=True, 
    framealpha=0.9, 
    ncol=1
)

# 18. 最终美化
ax.set_theta_offset(np.pi/2)
ax.set_theta_direction(-1)
ax.set_rlim(0, 1.15)
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
plt.title('Tumor Microenvironment Spatial Distribution', fontsize=16, pad=20)

plt.tight_layout()
plt.savefig('tumor_spatial_analysis_fixed.png', dpi=300, bbox_inches='tight')
plt.show()
