import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl

# 创建数据框（基于图片中的表格数据）
data = {
    'Symbol': ["SOCS1","SOCS3","SOCS3","SOCS3","SOCS3","SOCS3","SOCS3","SOCS3","SOCS3","STAT6","IFNGR1","IFNGR1","TYK2","STAT1","JAK3","JAK3","STAT2"],
    'logfoldchanges': [-0.835223,-5.2207017,-4.1098604,-3.9351864,-3.8919094,-2.829026,-2.5291514,-2.5210917,-2.0900233,0.5046044,0.5783798,0.86135674,0.89681333,0.97596735,1.0853907,1.3301098,1.5638305],
    'pct_nz_group': [0.13015717092337917,0.4542728635682159,0.3280828082808281,0.42178217821782177,0.48451829655861256,0.26964636542239684,0.20754290876242096,0.3262987012987013,0.31732840874115115,0.235543369890329,0.22872732635304527,0.1823297692740574,0.22504055449048813,0.1491608966083793,0.13641055891461437,0.2642073778664008,0.16230489379180849],
    'group': ["T cells","Mast cells","Epithelials","Macrophages","Fibroblasts","T cells","Kerationocytes","Myeloid","Stratum basal cells","T cells","Myeloid","Macrophages","Myeloid","Fibroblasts","Myeloid","T cells","Fibroblasts"]
}

df = pd.DataFrame(data)

# 创建自定义颜色映射：负值蓝色渐变，0为白色，正值红色渐变
colors = [(0, 'darkblue'), (0.4, 'lightblue'), (0.5, 'white'), (0.6, 'lightcoral'), (1, 'darkred')]
cmap = LinearSegmentedColormap.from_list('custom_red_blue', colors)

# 设置点的大小范围 (缩放因子基于pct_nz_group值)
min_size = 100
max_size = 700

# 计算点大小
size_norm = mpl.colors.Normalize(vmin=df['pct_nz_group'].min(), 
                                 vmax=df['pct_nz_group'].max())
sizes = min_size + (max_size - min_size) * size_norm(df['pct_nz_group'])

# 创建画布
fig, ax = plt.subplots(figsize=(10,6))

# 绘制点图
scatter = ax.scatter(
    x=df['Symbol'],
    y=df['group'],
    s=sizes,
    c=df['logfoldchanges'],
    cmap=cmap,
    vmin=-5,
    vmax=5,
    edgecolor='black',
    alpha=0.9
)
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)


# 设置轴标签
ax.set_xlabel('Symbol', fontsize=12, fontweight='bold')
ax.set_ylabel('Group', fontsize=12, fontweight='bold')
ax.set_title('Gene Expression Dot Plot', fontsize=14, fontweight='bold')
plt.xticks(rotation=45, ha='right')

# 1. 添加颜色图例（底部）
cax = fig.add_axes([0.83, 0.15, 0.03, 0.3])  # [left, bottom, width, height]
cbar = plt.colorbar(scatter, cax=cax, orientation='vertical')
cbar.set_label('logfoldchange', fontsize=10, labelpad=10)
cbar.set_ticks([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
cbar.ax.tick_params(labelsize=9)

# 2. 添加点大小图例（顶部）
size_legend_ax = fig.add_axes([0.83, 0.55, 0.03, 0.3])  # [left, bottom, width, height]

# 手动创建点大小图例（仅显示0.1, 0.25, 0.5）
legend_values = [0.1, 0.25, 0.5]
legend_sizes = min_size + (max_size - min_size) * size_norm(legend_values)

# 创建图例点
for size_val, size in zip(legend_values, legend_sizes):
    size_legend_ax.scatter([], [], s=size, c='black', label=f'{size_val:.2f}')

# 隐藏坐标轴
size_legend_ax.axis('off')

# 创建图例
size_legend = size_legend_ax.legend(
    loc='center',
    title='pct_nz_group',
    frameon=False,
    fontsize=9,
    title_fontsize=10,
    handletextpad=1.5,
    labelspacing=1.8  # 增加点与标签之间的间距
)

plt.tight_layout(rect=[0, 0, 0.8, 1])  # 为右侧图例留出空间
plt.show()
plt.savefig('gene_expression_dotplot.png', dpi=300, bbox_inches='tight')
plt.savefig('gene_expression_dotplot.pdf', bbox_inches='tight')
