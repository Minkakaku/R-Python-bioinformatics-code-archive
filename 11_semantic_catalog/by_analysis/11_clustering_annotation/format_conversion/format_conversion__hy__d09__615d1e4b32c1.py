import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import distance

# 1. 加载数据并确保名称唯一性
adata = sc.read_h5ad('integrated_data.h5ad')
adata.obs_names_make_unique()

# 2. 筛选目标细胞类型 (按图片顺序排列)
target_celltypes = [
    "CD4+ TEM", "CD4+ Treg", "CD8+ Tex", "M2", 
    "iDC", "THY1+ FB", "Kera Spinous", 
    "Kera Inflammatory", "CTMCs", "Disease Stem","POSTN+ FB"
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
ax.set_ylim(0, 1.3)

# 7. 定义环形半径位置 
ring_radii = {
    "0-10": 0.15, "10-50": 0.35, 
    "50-100": 0.55, "100-150": 0.75, 
    "150-200": 0.95, "≥200": 1.15
}

# 8. 绘制距离环标注线和标签
for zone, r in ring_radii.items():
    # 绘制虚线环
    ax.plot(np.linspace(0, 2*np.pi, 100), [r]*100, 
            color='gray', linestyle='--', linewidth=0.7, alpha=0.7)
    
    # 添加距离标签（在右侧）
    ax.text(0, r, f"{zone.replace('-', 'μm')}", 
            ha='left', va='center', fontsize=10,
            bbox=dict(facecolor='white', alpha=0.8))

# 9. 设置细胞类型扇区 
n_types = len(target_celltypes)
theta_step = 2 * np.pi / n_types

# 10. 细胞类型标签位置（按图片样式调整）
label_positions = [
    (1, "CD4+ TEM", "center"),        # 左上角
    (2, "CD4+ Treg", "center"),        # 顶部
    (3, "CD8+ Tex", "center"),         # 右上角
    (4, "M2", "center"),               # 右侧
    (5, "iDC", "center"),              # 右下角
    (6, "THY1+ FB", "center"),         # 底部
    (7, "Kera Spinous", "center"),     # 左下角
    (8, "Kera Inflammatory", "center"),# 左侧
    (9, "CTMCs", "center"),            # 左上角内侧
    (10, "Disease Stem", "center")      # 中心外侧
]

label_radius = 1.25
for i, (pos, label, align) in enumerate(label_positions):
    angle = (i - 3.5) * theta_step  # 调整起始位置
    ax.text(angle, label_radius, label, ha='center', va='center', 
            fontsize=10, rotation=0, 
            bbox=dict(boxstyle="round,pad=0.2", fc='white', ec='none', alpha=0.8))

# 11. 扇区分割线
for i in range(n_types + 1):
    angle = i * theta_step
    ax.plot([angle, angle], [0, 1.15], 
            color='lightgray', linewidth=0.8)

# 12. 聚合数据并添加气泡值映射
    # 气泡大小映射：基于图中展示的值
for _, row in agg_df.iterrows():
    ct = row['celltype_final']
    zone = row['distance_zone']
    group = row[group_col]
    count = row['count']  # 直接获取细胞数量
    
    # 计算气泡位置（分组偏移）
    angle_idx = target_celltypes.index(ct)
    center_angle = angle_idx * theta_step + theta_step/2
    angle_offset = 0.1 if group == 'disease' else -0.1
    angle = center_angle + angle_offset
    
    # 计算气泡半径和大小
    r = ring_radii.get(zone, 0.9)
    bubble_size = 50 + (count / max_count) * 1500  # 尺寸与细胞数量线性相关
    
    # 绘制气泡（不添加文本标签）
    ax.scatter(
        angle, r, 
        s=bubble_size, 
        color=group_colors[group], 
        alpha=0.7,
        edgecolors='w',
        zorder=5
    )


# 13. 绘制气泡（按分组和位置）
#min_size = 50
#max_size = 2000

# 创建角度索引映射

# 14. 中心标注 - Fibroblasts
ax.text(0, 0, "Fibroblasts", ha='center', va='center', fontsize=14,
        bbox=dict(boxstyle="circle,pad=0.3", fc='white', ec='black', alpha=0.8))

# 15. 图例系统（右下角）
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', label='Disease',
               markerfacecolor='#FF0000', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='Normal',
               markerfacecolor='#0000FF', markersize=10),
    plt.Line2D([0], [0], marker='o', color='w', label='100',
               markerfacecolor='grey', markersize=5),
    plt.Line2D([0], [0], marker='o', color='w', label='500',
               markerfacecolor='grey', markersize=8),
    plt.Line2D([0], [0], marker='o', color='w', label='1000',
               markerfacecolor='grey', markersize=12),
    plt.Line2D([0], [0], marker='o', color='w', label='2000',
               markerfacecolor='grey', markersize=15)
]

ax.legend(handles=legend_elements, loc='lower right', 
          bbox_to_anchor=(1.15, 0.15), title="Legend",
          frameon=True, framealpha=0.9)

# 16. 美化输出
ax.set_xticks([])
ax.set_yticks([])
ax.spines['polar'].set_visible(False)
plt.title('Spatial Distribution of Cell Types', fontsize=16, pad=20)
plt.tight_layout()
plt.savefig('spatial_bubble_distribution_full.png', dpi=300, bbox_inches='tight')
plt.show()
