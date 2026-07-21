import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import distance

# 加载数据
adata = sc.read_h5ad('integrated_data.h5ad')
# 确保观测名称唯一
adata.obs_names_make_unique()

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
distance_zones = [f"<{distances[0]}"] 
for i in range(1, len(distances)):
    distance_zones.append(f"{distances[i-1]}-{distances[i]}")
distance_zones.append(f">={distances[-1]}")

# 修正：避免索引超出范围的函数
def assign_zone(d):
    if d < distances[0]:
        return distance_zones[0]
    for i in range(1, len(distances)):
        if distances[i-1] <= d < distances[i]:
            return distance_zones[i]
    return distance_zones[-1]

adata.obs['distance_zone'] = adata.obs['distance_to_center'].apply(assign_zone)

# 4. 过滤中心细胞
adata = adata[adata.obs.index != tumor_idx].copy()

# 5. 创建极坐标系
plt.figure(figsize=(15, 15))
ax = plt.subplot(111, polar=True)

# 6. 设置环形参数
num_zones = len(distance_zones)
zone_radii = np.linspace(0.1, 0.9, num_zones)
max_radius = 1.0  # 最外层半径

# 7. 获取所有细胞类型（排除中心细胞类型）
all_cell_types = sorted(target_celltypes - {tumor_cell_type})
theta_step = 2 * np.pi / len(all_cell_types)
cell_type_angles = {ct: i * theta_step for i, ct in enumerate(all_cell_types)}

# 8. 绘制距离环的网格线并标注距离
for i, dist in enumerate(distances):
    # 计算半径位置
    r_pos = 0.1 + (dist / distances[-1]) * 0.8
    
    # 绘制虚线环
    theta_line = np.linspace(0, 2*np.pi, 100)
    r_line = np.ones(100) * r_pos
    ax.plot(theta_line, r_line, color='gray', linestyle='--', alpha=0.7)
    
    # 添加距离标注（统一放在右侧）
    label_radius = 1.15  # 标签半径位置
    text_angle = 0  # 右侧位置，角度为0
    ax.text(text_angle, r_pos, f"{dist}μm", 
            ha='left', va='center', fontsize=10, fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.2'),
            transform=ax.transData._b)  # 使用数据坐标变换

# 9. 添加细胞类型分割线
for i, ct in enumerate(all_cell_types):
    angle = cell_type_angles[ct]
    # 从中心到最外层绘制分割线
    ax.plot([angle, angle], [0, max_radius], color='black', linewidth=1, alpha=0.5)

# 10. 中心位置添加细胞类型标签 
ax.text(0, 0, tumor_cell_type, 
        ha='center', va='center', fontsize=14, fontweight='bold',
        bbox=dict(boxstyle="circle,pad=0.3", fc='white', ec='black', alpha=0.8))

# 11. 确保所有细胞类型标签在最外层（统一位置）
label_radius = 1.15  # 标签半径位置
for ct, angle in cell_type_angles.items():
    # 计算标签位于扇形中心的角度
    label_angle = angle + theta_step/2
    
    # 计算坐标
    x = label_radius * np.cos(label_angle)
    y = label_radius * np.sin(label_angle)
    
    # 计算旋转角度（径向标签）
    rotation_angle = np.degrees(label_angle)
    # 调整角度使得文字方向正确
    if rotation_angle > 90 and rotation_angle < 270:
        rotation_angle += 180
    # 绘制标签
    ax.text(x, y, ct, 
            ha='center', va='center', 
            fontsize=10, fontweight='bold',
            rotation=rotation_angle,
            rotation_mode='anchor',
            bbox=dict(boxstyle="round,pad=0.2", fc='white', ec='none', alpha=0.8))

# 12. 聚合数据（按类型/区域/分组）
agg_df = adata.obs.groupby(['celltype_final', 'distance_zone', group_col]).size().reset_index(name='count')

# 13. 绘制分组点（红/蓝）
size_factor = 50  # 点大小缩放因子
for _, row in agg_df.iterrows():
    ct, zone, group = row['celltype_final'], row['distance_zone'], row[group_col]
    if ct == tumor_cell_type: 
        continue
    
    zone_idx = distance_zones.index(zone)
    r = zone_radii[zone_idx]
    angle_center = cell_type_angles[ct] + theta_step/2
    
    # 分组偏移量（避免重叠）
    group_offset = 0.1 if group == 'disease' else -0.1
    
    # 绘制点（大小表示细胞数量）
    s = max(20, min(row['count'] * size_factor, 500))  # 限制点大小范围
    ax.scatter(angle_center + group_offset, r, s=s,
               color=group_colors[group], alpha=0.8, 
               edgecolors='w', linewidth=0.5)

# 14. 右下角添加图例 
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

# 调整图例位置到右下角
ax.legend(handles=legend_elements, title='Legend', 
          loc='lower right', bbox_to_anchor=(1.3, 0.0), 
          frameon=True, framealpha=0.9, ncol=1)

# 15. 美化图形
ax.set_theta_offset(np.pi/2)  # 0°位于顶部
ax.set_theta_direction(-1)     # 顺时针方向
ax.set_rlim(0, max_radius + 0.2)
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
plt.title('Tumor Microenvironment Spatial Distribution', fontsize=16, pad=20)

plt.tight_layout()
plt.savefig('tumor_spatial_analysis.png', dpi=300, bbox_inches='tight')
plt.show()
