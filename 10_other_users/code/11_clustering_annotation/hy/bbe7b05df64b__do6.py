import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import distance
from matplotlib.path import Path
import matplotlib.patches as patches

# 1. 加载数据并确保名称唯一性
adata = sc.read_h5ad('integrated_data.h5ad')
adata.obs_names_make_unique()

# 2. 筛选目标细胞类型 (按图片顺序)
target_celltypes = [
    "CD4+ Tem", 
    "CD4+ Treg", 
    "CD8+ Tex", 
    "M2", 
    "iDC", 
    "POSTN+ FB", 
    "THY1+ FB", 
    "Kera Spinous", 
    "Kera Inflammatory", 
    "CTMCs",
    "Disease Stem"
]
cell_mask = adata.obs['celltype_final'].isin(target_celltypes)
adata = adata[cell_mask].copy()

# 3. 配置参数
tumor_cell_type = 'POSTN+ FB'  # 中心细胞类型（成纤维细胞）
cell_type_col = 'celltype_final'
group_col = 'group'
distances = [10, 50, 100, 150, 200]  # 距离环
group_colors = {'disease': '#FF0000', 'normal': '#0000FF'}  # 红/蓝配色

# 4. 定位中心细胞 - 随机选择一个成纤维细胞
tumor_cells = adata.obs[adata.obs[cell_type_col] == tumor_cell_type]
if len(tumor_cells) == 0:
    raise ValueError(f"未找到{tumor_cell_type}细胞")
center_cell_id = tumor_cells.index[0]  # 随机选择一个

# 5. 计算所有细胞到中心的距离
center_coords = adata.obsm['spatial'][adata.obs.index == center_cell_id][0]
all_coords = adata.obsm['spatial']
distances_to_center = distance.cdist([center_coords], all_coords)[0]
adata.obs['distance_to_center'] = distances_to_center

# 6. 定义距离区域
def assign_distance_zone(d):
    if d < 10:
        return "<10"
    elif 10 <= d < 50:
        return "10-50"
    elif 50 <= d < 100:
        return "50-100"
    elif 100 <= d < 150:
        return "100-150"
    elif 150 <= d < 200:
        return "150-200"
    else:
        return "≥200"

adata.obs['distance_zone'] = adata.obs['distance_to_center'].apply(assign_distance_zone)

# 7. 创建图表
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='polar')
ax.set_theta_offset(np.pi/2)  # 0度在顶部
ax.set_theta_direction(-1)    # 顺时针方向
ax.set_ylim(0, 1.2)           # 半径范围

# 8. 定义距离环位置
ring_radii = {
    "<10": 0.1,
    "10-50": 0.3,
    "50-100": 0.5,
    "100-150": 0.7,
    "150-200": 0.9,
    "≥200": 1.1
}

# 9. 绘制距离环标签（按图片样式）
dist_positions = {
    10: (0, 0.3, "10μm"),
    50: (0, 0.5, "50μm"),
    100: (0, 0.7, "100μm"),
    150: (0, 0.9, "150μm"),
    200: (0, 1.1, "200μm")
}

for dist, (x, r, text) in dist_positions.items():
    ax.text(x, r, text, ha='left', va='center', fontsize=10, 
           bbox=dict(facecolor='white', alpha=0.8))

# 10. 设置细胞类型扇区
n_types = len(target_celltypes)
theta_step = 2 * np.pi / n_types

# 11. 绘制细胞类型标签（按图片样式）
label_radius = 1.2  # 标签在最外层

label_positions = [
    (theta_step * 0.5, "CD4+ Tem"),
    (theta_step * 1.5, "CD4+ Treg"),
    (theta_step * 2.5, "CD8+ Tex"),
    (theta_step * 3.5, "M2"),
    (theta_step * 4.5, "iDC"),
    (theta_step * 5.5, "THY1+ FB"),
    (theta_step * 6.5, "Kera Spinous"),
    (theta_step * 7.5, "Kera Inflammatory"),
    (theta_step * 8.5, "CTMCs"),
    (theta_step * 9.5, "Disease Stem")
]

for angle, label in label_positions:
    # 计算文本旋转角度
    rotation_angle = np.degrees(angle)
    if rotation_angle > 90 and rotation_angle < 270:
        rotation_angle -= 180
    
    # 绘制标签
    ax.text(angle, label_radius, label, ha='center', va='center', 
            rotation=rotation_angle, fontsize=10, 
            bbox=dict(boxstyle="round,pad=0.2", fc='white', ec='none', alpha=0.8))

# 12. 绘制扇区分割线
for i in range(n_types + 1):
    ax.plot([i * theta_step, i * theta_step], [0, 1.15], color='gray', linewidth=0.5)

# 13. 聚合数据
agg_df = adata.obs.groupby(['celltype_final', 'distance_zone', group_col]).size().reset_index(name='count')

# 14. 绘制气泡（按分组）
size_factor = 0.01  # 大小缩放因子

# 固定气泡大小示例
bubble_size_map = {
    "CD4+ Tem": 100,
    "CD4+ Treg": 100,
    "CD8+ Tex": 100,
    "M2": 100,
    "iDC": 100,
    "THY1+ FB": 100,
    "Kera Spinous": 100,
    "Kera Inflammatory": 100,
    "CTMCs": 100,
}

for _, row in agg_df.iterrows():
    ct = row['celltype_final']
    if ct == tumor_cell_type: continue  # 跳过中心细胞
    
    # 匹配扇区角度
    angle_idx = target_celltypes.index(ct)
    angle = angle_idx * theta_step + theta_step/2
    r = ring_radii.get(row['distance_zone'], 0.9)
    
    # 使用预设气泡大小
    s = bubble_size_map.get(ct, 100)
    
    # 绘制气泡
    ax.scatter(angle, r, s=s, 
               color=group_colors[row[group_col]], 
               alpha=0.8, 
               edgecolors='w', 
               zorder=5)
    
    # 添加数量标签（按图片样式）
    ax.text(angle, r, str(s), ha='center', va='center', fontsize=8)

# 15. 添加中心细胞标签（Fibroblasts）
ax.text(0, 0, "POSTN+ FB", ha='center', va='center', fontsize=12,
        bbox=dict(boxstyle="circle,pad=0.3", fc='white', ec='black', alpha=0.8))

# 16. 隐藏不必要的元素
ax.set_xticks([])
ax.set_yticks([])
ax.spines['polar'].set_visible(False)

# 17. 添加图例（按图片右下角位置）
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', label='disease',
               markerfacecolor='#FF0000', markersize=8),
    plt.Line2D([0], [0], marker='o', color='w', label='normal',
               markerfacecolor='#0000FF', markersize=8),
    plt.Line2D([0], [0], marker='o', color='w', label='100',
               markerfacecolor='gray', markersize=5),
    plt.Line2D([0], [0], marker='o', color='w', label='500',
               markerfacecolor='gray', markersize=8),
    plt.Line2D([0], [0], marker='o', color='w', label='1000',
               markerfacecolor='gray', markersize=11),
    plt.Line2D([0], [0], marker='o', color='w', label='2000',
               markerfacecolor='gray', markersize=14)
]

# 在右下角添加图例
plt.legend(handles=legend_elements, title='Legend', loc='lower right', 
           bbox_to_anchor=(1.15, 0.15), frameon=True)

# 18. 保存并显示
plt.tight_layout()
plt.savefig('tumor_spatial_analysis_corrected.png', dpi=300, bbox_inches='tight')
plt.show()
