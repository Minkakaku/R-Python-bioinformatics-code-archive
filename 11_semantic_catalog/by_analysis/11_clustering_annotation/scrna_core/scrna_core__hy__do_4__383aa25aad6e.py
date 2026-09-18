import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import ticker
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable

# 创建包含数据的DataFrame
data = {
    "Symbol": ["SOCS3", "SOCS3", "SOCS3", "SOCS3", "SOCS1", "SOCS3", "SOCS3", "SOCS3", "SOCS3", 
               "STAT2", "STAT1", "JAK3", "STAT6", "TYK2", "JAK3", "IFNGR1", "IFNGR1"],
    "logfoldchanges": [-2.5291514, -2.0900233, -3.8919094, -2.829026, -0.835223, -5.2207017, 
                       -2.5210917, -3.9351864, -4.1098604, 1.5638305, 0.97596735, 1.3301098, 
                       0.5046044, 0.89681333, 1.0853907, 0.5783798, 0.86135674],
    "pct_nz_group": [0.20754290876242096, 0.31732840874115115, 0.48451829655861256, 
                     0.26964636542239684, 0.13015717092337917, 0.4542728635682159, 
                     0.3262987012987013, 0.42178217821782177, 0.3280828082808281, 
                     0.16230489379180849, 0.1491608966083793, 0.2642073778664008, 
                     0.235543369890329, 0.22504055449048813, 0.13641055891461437, 
                     0.22872732635304527, 0.1823297692740574],
    "group": ["Kerationocytes", "Stratum basal cells", "Fibroblasts", "T cells", "T cells", 
              "Mast cells", "Myeloid", "Macrophages", "Epithelials", "Fibroblasts", 
              "Fibroblasts", "T cells", "T cells", "Myeloid", "Myeloid", "Myeloid", "Macrophages"]
}

df = pd.DataFrame(data)

# 创建自定义颜色映射：负值为蓝色渐变，正值为红色渐变
colors_neg = ["#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"]
colors_pos = ["#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d"]
colors = colors_neg + colors_pos
cmap = LinearSegmentedColormap.from_list('custom_diverging', colors, N=256)

# 设置色标的范围（基于logfoldchanges的最大绝对值）
max_abs_value = np.abs(df['logfoldchanges']).max()
vmin, vmax = -max_abs_value, max_abs_value

# 设置点的大小范围
size_min, size_max = df['pct_nz_group'].min() * 800, df['pct_nz_group'].max() * 800

# 创建图表
plt.figure(figsize=(14, 10))  # 增大画布尺寸
ax = sns.scatterplot(
    data=df,
    x="Symbol",
    y="group",
    size="pct_nz_group",
    hue="logfoldchanges",
    palette=cmap,
    hue_norm=(vmin, vmax),
    sizes=(size_min, size_max),
    alpha=0.9,
    edgecolor="black",
    linewidth=0.8
)

# 添加标题和标签
plt.title("Cell Expression by Cell Type", fontsize=16, pad=16)
plt.xlabel("Genes", fontsize=14, labelpad=12)
plt.ylabel("Cell Types", fontsize=14, labelpad=12)
plt.xticks(rotation=45, ha='right')  # 旋转x轴标签
plt.grid(True, which='both', axis='y', linestyle='--', alpha=0.3)

# 创建分隔器用于颜色条布局
divider = make_axes_locatable(ax)

# 在右侧添加一个短的轴用于颜色图例
cax = divider.append_axes("right", size="5%", pad=0.1)

# 添加颜色条到这个小轴上
norm = plt.Normalize(vmin, vmax)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
colorbar = plt.colorbar(sm, cax=cax)
colorbar.ax.set_ylabel('logfoldchanges', fontsize=12, labelpad=10)
colorbar.set_ticks(np.arange(np.floor(vmin), np.ceil(vmax)+1))

# 关键修复：创建一个新的轴用于点大小图例
# 而不是尝试将同一个图例添加到多个坐标轴
size_ax = divider.append_axes("right", size="5%", pad=0.6)  # 在颜色条下方
size_ax.set_axis_off()  # 隐藏这个新轴的边框和刻度

# 创建点大小图例（只保留0.1, 0.25, 0.5）
target_sizes = [0.1, 0.25, 0.5]
size_handles = []
size_labels = []

for size in target_sizes:
    # 计算点大小（基于与数据中大小值的比例）
    marker_size = np.sqrt(size * 800)  # 与绘图使用相同的缩放
    
    # 创建图例元素
    size_handles.append(
        Line2D([0], [0], 
               marker='o', 
               color='w',
               markerfacecolor='gray',
               markersize=marker_size,
               markeredgecolor='black',
               markeredgewidth=0.8,
               linewidth=0,
               label=f'{size}')
    )
    size_labels.append(f'{size}')

# 在新轴上创建点大小图例
size_ax.legend(
    size_handles, 
    size_labels,
    loc='center',
    frameon=True,
    framealpha=1.0,
    title='pct_nz_group',
    title_fontsize=12
)

# 调整布局
plt.tight_layout(pad=2.0)

# 保存高质量图片
plt.savefig('gene_expression_dotplot.png', dpi=300, bbox_inches='tight')
plt.savefig('gene_expression_dotplot.pdf', bbox_inches='tight')

# 显示图表
plt.show()
