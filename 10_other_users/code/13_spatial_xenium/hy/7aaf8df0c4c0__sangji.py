import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.patches import Polygon
import matplotlib.colors as mcolors
import os

# ================== 字体与CNS设置 ==================
font_manager.fontManager.addfont(
    '/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf'
)
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 16,
    'font.weight': 'normal',
    'pdf.fonttype': 42
})

# ================== 你的自定义颜色字典（删掉不存在的数据，仅保留附图颜色） ==================
color_dict = {
    'B cells':           '#1f77b4',  # 蓝色
    'Endothelials':      '#ff7f0e',  # 橙色
    'Epithelials':       '#2ca02c',  # 绿色
    'Fibroblasts':       '#d62728',  # 红色
    'Macrophages':       '#9467bd',  # 紫色
    'Neutrophils':       '#8c564b',  # 棕褐色
    'Plasmas':           '#e377c2',  # 粉红色
    'T cells':           '#bcbd22'   # 黄绿色
}

# ================== 数据加载与清洗 ==================
adata = sc.read('data_with_group_final.h5ad')

# 提取 N, A, C 占比数据
obs_df = adata.obs[['group_final', 'celltype']]
count_df = pd.crosstab(obs_df['group_final'], obs_df['celltype'])
percent_df = count_df.div(count_df.sum(axis=1), axis=0) * 100
percent_df = percent_df.reindex(['POST-NR', 'POST-R', 'PRE-NR', 'PRE-R'])

# 获取细胞类型列表（按数据的实际列排序）
celltypes = list(percent_df.columns)

# 将你的自定义颜色映射到当前存在的细胞类型上（没有在字典里的用灰色兜底）
celltype_colors = {ct: color_dict.get(ct, '#808080') for ct in celltypes}

# X轴映射和柱子宽度配置
x_map = {'POST-NR': 1, 'POST-R': 2, 'PRE-NR': 3, 'PRE-R': 4}  # 【已修复 1：PRE-R 写错的问题】
bar_width = 0.20  
half_width = bar_width / 2

# 准备画布
fig, ax = plt.subplots(figsize=(12, 8))
ax.set_xlim(0.4, 4.6)  # 【已修复 2：调整 x 轴范围，确保第4根柱子不超出画布】
ax.set_ylim(0, 105)

# ================== 计算并绘制淡色直连色块 ==================
group_finals = ['POST-NR', 'POST-R', 'PRE-NR', 'PRE-R']

# 预先计算每个状态、每种细胞类型的 y 轴底部和顶部坐标
y_ranges = {}
for group_final in group_finals:
    y_current = 0
    ranges = {}
    for ct in celltypes:
        val = percent_df.loc[group_final, ct]
        ranges[ct] = (y_current, y_current + val)
        y_current += val
    y_ranges[group_final] = ranges

# 画连接色块
for i in range(len(group_finals) - 1):
    group_final_left = group_finals[i]
    group_final_right = group_finals[i+1]
    
    x_start = x_map[group_final_left] + half_width
    x_end = x_map[group_final_right] - half_width
    
    for ct in celltypes:
        y_left_bottom, y_left_top = y_ranges[group_final_left][ct]
        y_right_bottom, y_right_top = y_ranges[group_final_right][ct]
        
        if y_left_top - y_left_bottom > 0 and y_right_top - y_right_bottom > 0:
            # 核心技巧：使用 to_rgba 将十六进制颜色加上透明度 (alpha=0.35)
            link_color = mcolors.to_rgba(celltype_colors[ct], alpha=0.35)
            
            link_poly = Polygon(
                [
                    (x_start, y_left_top),      # 左上
                    (x_end, y_right_top),       # 右上
                    (x_end, y_right_bottom),    # 右下
                    (x_start, y_left_bottom)    # 左下
                ],
                color=link_color, # 使用带透明度的淡色
                edgecolor=None
            )
            ax.add_patch(link_poly)

# ================== 绘制实色柱子（覆盖在上面） ==================
for group_final in group_finals:
    x_center = x_map[group_final]
    y_current = 0
    for ct in celltypes:
        val = percent_df.loc[group_final, ct]
        if val > 0:
            ax.bar(
                x_center, 
                val, 
                bottom=y_current, 
                width=bar_width, 
                color=celltype_colors[ct], # 原色、实心
                edgecolor='white', # 白色描边
                linewidth=0.8,
                alpha=1.0
            )
            y_current += val

# ================== 图表美化与保存 ==================
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(['POST-NR', 'POST-R', 'PRE-NR', 'PRE-R'], fontsize=18, fontweight='bold')
ax.set_ylabel('Fraction of cells (%)', fontsize=16)

# 去除网格
ax.grid(False)
# 去除顶部和右侧边框
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 生成图例（使用自定义颜色）
handles = [plt.Rectangle((0,0),1,1, color=celltype_colors[ct]) for ct in celltypes]
ax.legend(handles, celltypes, title='Annotation', bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

plt.tight_layout()

out_dir = '.'
plt.savefig(os.path.join(out_dir, 'alluvial_custom_colors.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'alluvial_custom_colors.pdf'), bbox_inches='tight')

plt.show()
