import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
import pandas as pd
from scipy.spatial import distance

# 加载数据
adata1 = sc.read_h5ad('integrated_data.h5ad')

# 定义目标细胞类型
target_celltypes = {
    "CD4+ TEM", "CD4+ Treg", "CD8+ Tex", "M2", "CTMCs", 
    "iDC", "POSTN+ FB", "THY1+ FB", "Kera Spinous", "Kera Inflammatory"
}

# 筛选目标细胞类型
cell_mask = adata1.obs['celltype_final'].isin(target_celltypes)
adata = adata1[cell_mask, :].copy()

# 配置参数
tumor_cell_type = 'POSTN+ FB'  # 设定为中心细胞的类型
cell_type_col = 'celltype_final'  # 细胞类型列名
group_col = 'group'  # 分组列名
distances = [10, 50, 100, 150, 200]  # 距离环设置
group_colors = {'disease': '#FF0000', 'normal': '#00FF00'}  # 分组颜色

# 1. 查找肿瘤细胞作为中心
tumor_indices = adata.obs[adata.obs[cell_type_col] == tumor_cell_type].index
if len(tumor_indices) == 0:
    raise ValueError(f"未找到标记为'{tumor_cell_type}'的细胞")

# 选择第一个肿瘤细胞作为中心点
center_cell_idx = tumor_indices[0]
center_coords = adata.obsm['spatial'][adata.obs.index == center_cell_idx][0]

# 2. 计算所有细胞到中心肿瘤细胞的距离
all_coords = adata.obsm['spatial']
distances_to_center = distance.cdist([center_coords], all_coords)[0]
adata.obs['distance_to_center'] = distances_to_center

# 3. 定义距离区域
distance_zones = [f"<{distances[0]}"] 
for i in range(1, len(distances)):
    distance_zones.append(f"{distances[i-1]}-{distances[i]}")
distance_zones.append(f">={distances[-1]}")

def get_distance_zone(dist):
    if dist < distances[0]:
        return distance_zones[0]
    for i in range(1, len(distances)):
        if dist < distances[i]:
            return distance_zones[i]
    return distance_zones[-1]

adata.obs['distance_zone'] = adata.obs['distance_to_center'].apply(get_distance_zone)

# 4. 过滤掉中心肿瘤细胞自身
adata = adata[adata.obs.index != center_cell_idx].copy()

# 5. 创建极坐标系图
plt.figure(figsize=(15, 15))
ax = plt.subplot(111, polar=True)

# 6. 设置极坐标参数
num_zones = len(distance_zones)
zone_radii = np.linspace(0.1, 0.9, num_zones)  # 每个区域的半径位置

# 7. 获取细胞类型和分组信息
cell_types = sorted(adata.obs[cell_type_col].unique())
groups = sorted(adata.obs[group_col].unique())

# 创建角度映射 - 每种细胞类型一个角度
theta_step = 2 * np.pi / len(cell_types)
cell_type_angles = {cell_type: i * theta_step for i, cell_type in enumerate(cell_types)}

# 8. 绘制距离环的网格线
distance_colors = plt.cm.viridis(np.linspace(0, 1, len(distances)))
for i, dist in enumerate(distances):
    # 计算半径位置（距离映射到半径）
    r_position = 0.1 + (dist / distances[-1]) * 0.8
    
    # 绘制虚线网格线
    theta_line = np.linspace(0, 2*np.pi, 100)
    r_line = np.ones(100) * r_position
    ax.plot(theta_line, r_line, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    # 添加距离标签
    ax.text(0, r_position + 0.03, f"{dist}μm", 
            ha='center', va='center', fontsize=9, 
            bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))

# 9. 聚合数据 - 计算每个区域/类型/分组的细胞数量
agg_df = adata.obs.groupby(
    ['celltype_final', 'distance_zone', group_col]
).size().reset_index(name='cell_count')

# 10. 绘制细胞点 - 大小表示细胞数量
scatter_objects = {'cell_types': {}, 'groups': {}}  # 用于图例
size_factor = 50  # 点大小缩放因子

for _, row in agg_df.iterrows():
    cell_type = row['celltype_final']
    zone = row['distance_zone']
    group = row[group_col]
    count = row['cell_count']
    
    # 获取区域索引
    zone_idx = distance_zones.index(zone)
    
    # 计算角度和半径
    theta = cell_type_angles[cell_type] + 0.5 * theta_step  # 扇形中心
    r = zone_radii[zone_idx]  # 环半径
    
    # 根据分组调整角度位置（避免重叠）
    group_offset = 0
    if group == 'disease':
        group_offset = -0.1  # 向左偏移
    elif group == 'normal':
        group_offset = 0.1   # 向右偏移
    
    # 绘制点（大小表示细胞数量）
    s = max(10, min(count * size_factor, 500))  # 限制点的大小范围
    scatter = ax.scatter(
        theta + group_offset, 
        r, 
        s=s,
        color=group_colors[group],
        alpha=0.8,
        edgecolors='w',
        linewidth=0.5
    )
    
    # 保存图例参考
    if cell_type not in scatter_objects['cell_types']:
        scatter_objects['cell_types'][cell_type] = scatter
    if group not in scatter_objects['groups']:
        scatter_objects['groups'][group] = scatter

# 11. 添加细胞类型标签（最外层）
label_radius = 1.05  # 标签位于环外
for cell_type, angle in cell_type_angles.items():
    x = label_radius * np.cos(angle + 0.5 * theta_step)
    y = label_radius * np.sin(angle + 0.5 * theta_step)
    
    ax.text(
        x, y, cell_type, 
        ha='center', va='center', fontsize=10,
        rotation=(angle + 0.5 * theta_step) * 180/np.pi - 90,  # 径向标签
        bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3')
    )

# 12. 设置极坐标显示参数
ax.set_theta_offset(np.pi/2)  # 0度在顶部
ax.set_theta_direction(-1)    # 顺时针方向
ax.set_rlim(0, 1.1)           # 半径范围
ax.grid(False)                # 关闭网格
ax.set_xticks([])             # 关闭角度刻度
ax.set_yticks([])             # 关闭径向刻度

# 13. 添加图例
# 分组图例
group_legend = [plt.Line2D([0], [0], marker='o', color='w', label=group,
                         markerfacecolor=group_colors[group], markersize=10)
              for group in scatter_objects['groups'].keys()]

# 细胞类型图例
cell_type_legend = [plt.Line2D([0], [0], marker='o', color='w', label=cell_type,
                             markerfacecolor='gray', markersize=8)
                  for cell_type in scatter_objects['cell_types'].keys()]

# 点大小图例（表示细胞数量）
size_legend = []
for count in [1, 5, 10, 20]:
    s = count * size_factor
    size_legend.append(plt.Line2D([0], [0], marker='o', color='w', label=f'{count} cells',
                               markerfacecolor='gray', markersize=np.sqrt(s)/2))

# 创建图例
ax.legend(handles=group_legend, title='分组', 
          loc='center left', bbox_to_anchor=(1.15, 0.8), 
          frameon=True, framealpha=0.8, title_fontsize=12)

ax.legend(handles=cell_type_legend, title='细胞类型', 
          loc='center left', bbox_to_anchor=(1.15, 0.6), 
          frameon=True, framealpha=0.8, title_fontsize=12)

ax.legend(handles=size_legend, title='细胞数量', 
          loc='center left', bbox_to_anchor=(1.15, 0.4), 
          frameon=True, framealpha=0.8, title_fontsize=12)

# 保持图例
plt.gca().add_artist(ax.legend(handles=group_legend, title='分组', 
                              loc='center left', bbox_to_anchor=(1.15, 0.8), 
                              frameon=True, framealpha=0.8, title_fontsize=12))
plt.gca().add_artist(ax.legend(handles=cell_type_legend, title='细胞类型', 
                              loc='center left', bbox_to_anchor=(1.15, 0.6), 
                              frameon=True, framealpha=0.8, title_fontsize=12))

plt.tight_layout()
plt.savefig('肿瘤细胞空间分布分析1.png', dpi=300, bbox_inches='tight')
plt.show()
