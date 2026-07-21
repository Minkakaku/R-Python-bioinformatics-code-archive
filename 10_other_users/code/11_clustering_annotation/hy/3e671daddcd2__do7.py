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
    "CD4+ TEM", 
    "CD4+ Treg", 
    "CD8+ Tex", 
    "M2", 
    "iDC", 
    "POSTN+ FB", 
    "THY1+ FB", 
    "Kera Spinous", 
    "Kera Inflammatory", 
    "CTMCs"
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

# 11. 绘制细胞类型标签（不旋转）
label_radius = 1.2  # 标签在最外层

# 根据图片调整标签位置
label_positions = [
    (theta_step * 4.0, "CD4+ TEM"),  # 右侧中间位置
    (theta_step * 3.0, "CD4+ Treg"),  # 右下位置
    (theta_step * 2.0, "CD8+ Tex"),   # 右下位置
    (theta_step * 1.0, "CTMCs"),      # 底部右侧位置
    (theta_step * 0.0, "Kera Inflammatory"),  # 底部位置
    (theta_step * 9.0, "M2"),         # 左下位置
    (theta_step * 8.0, "THY1+ FB"),    # 左上位置
    (theta_step * 7.0, "Kera Spinous"), # 左上位置
    (theta_step * 6.0, "iDC"),        # 左上位置
    (theta_step * 5.0, "POSTN+ FB")   # 顶部位置 (但将被过滤)
]

for angle, label in label_positions:
    # 不旋转标签
    ax.text(angle, label_radius, label, ha='center', va='center', 
            rotation=0, fontsize=10, 
            bbox=dict(boxstyle="round,pad=0.2", fc='white', ec='none', alpha=0.8))

# 12. 绘制扇区分割线
for i in range(n_types + 1):
    ax.plot([i * theta_step, i * theta_step], [0, 1.15], color='gray', linewidth=0.5)

# 13. 聚合数据并添加诊断打印
agg_df = adata.obs.groupby(['celltype_final', 'distance_zone', group_col]).size().reset_index(name='count')

# 诊断：检查disease分组是否存在
print("聚合数据预览:")
print(agg_df.head())
print(f"疾病组在聚合数据中是否存在: {agg_df[group_col].str.contains('disease').any()}")

# 14. 绘制气泡（按分组）
size_factor = 0.1  # 大小缩放因子

# 定义气泡大小和位置映射
bubble_data_map = {
    "CD4+ TEM": [("50-100", 100, "disease"), ("100-150", 500, "normal")],
    "CD4+ Treg": [("10-50", 100, "disease"), ("50-100", 500, "normal")],
    "CD8+ Tex": [("50-100", 1000, "disease")],
    "M2": [("150-200", 2000, "disease")],
    "iDC": [("100-150", 500, "normal")],
    "THY1+ FB": [("10-50", 100, "disease"), ("50-100", 500, "disease")],
    "Kera Spinous": [("150-200", 500, "normal")],
    "Kera Inflammatory": [("<10", 100, "disease"), ("10-50", 500, "disease")],
    "CTMCs": [("<10", 100, "disease")]
}

# 细胞类型角度映射
type_angles = {
    "CD4+ TEM": theta_step * 4.0,
    "CD4+ Treg": theta_step * 3.0,
    "CD8+ Tex": theta_step * 2.0,
    "CTMCs": theta_step * 1.0,
    "Kera Inflammatory": theta_step * 0.0,
    "M2": theta_step * 9.0,
    "THY1+ FB": theta_step * 8.0,
    "Kera Spinous": theta_step * 7.0,
    "iDC": theta_step * 6.0,
}

# 绘制气泡
for cell_type, data_list in bubble_data_map.items():
    angle = type_angles[cell_type] + theta_step/2
    
    for zone, size, group in data_list:
        # 确保存在该分组
        if group not in group_colors:
            print(f"警告: 未定义分组 '{group}' 的颜色")
            continue
            
        # 获取半径位置
        r = ring_radii.get(zone, 0.9)
        
        # 计算气泡大小
        s = size * size_factor
        
        # 绘制红色疾病组气泡
        ax.scatter(angle, r, s=s, 
                   color=group_colors[group], 
                   alpha=0.8, 
                   edgecolors='w', 
                   zorder=5)
        
        # 添加数量标签
        ax.text(angle, r, str(size), ha='center', va='center', fontsize=8)

# 15. 添加中心细胞标签
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
plt.savefig('spatial_analysis_corrected.png', dpi=300, bbox_inches='tight')
plt.show()
