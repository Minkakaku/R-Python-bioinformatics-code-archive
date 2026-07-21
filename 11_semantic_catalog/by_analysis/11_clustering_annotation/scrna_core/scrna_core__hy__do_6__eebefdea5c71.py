import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# 创建示例数据
data = {
    'Symbol': ['STAT2', 'JAK3', 'JAK3', 'STAT1', 'TYK2', 'IFNGR1', 'IFNGR1', 'STAT6', 
               'SOCS1', 'SOCS3', 'SOCS3', 'SOCS3', 'SOCS3', 'SOCS3', 'SOCS3', 'SOCS3', 'SOCS3'],
    'logfoldchanges': [1.5638305, 1.3301098, 1.0853907, 0.97596735, 0.89681333, 0.86135674, 
                       0.5783798, 0.5046044, -0.835223, -2.0900233, -2.5210917, -2.5291514, 
                       -2.829026, -3.8919094, -3.9351864, -4.1098604, -5.2207017],
    'pct_nz_group': [0.16230489379180849, 0.2642073778664008, 0.13641055891461437, 
                     0.1491608966083793, 0.22504055449048813, 0.1823297692740574, 
                     0.22872732635304527, 0.235543369890329, 0.13015717092337917, 
                     0.31732840874115115, 0.3262987012987013, 0.20754290876242096, 
                     0.26964636542239684, 0.48451829655861256, 0.42178217821782177, 
                     0.3280828082808281, 0.4542728635682159],
    'group': ['Fibroblasts', 'T cells', 'Myeloid', 'Fibroblasts', 'Myeloid', 'Macrophages', 
              'Myeloid', 'T cells', 'T cells', 'Stratum basal cells', 'Myeloid', 
              'Kerationocytes', 'T cells', 'Fibroblasts', 'Macrophages', 'Epithelials', 
              'Mast cells']
}

df = pd.DataFrame(data)

# 创建自定义颜色映射
colors = [(0, 'darkblue'), (0.4, 'lightblue'), (0.5, 'white'), (0.6, 'lightcoral'), (1, 'darkred')]
cmap = LinearSegmentedColormap.from_list('custom_red_blue', colors)

# 设置点的大小范围（调整为更合理的值）
min_size = 50
max_size = 300  # 减小最大点大小，避免无效值

# 创建归一化对象
size_norm = mpl.colors.Normalize(vmin=df['pct_nz_group'].min(), 
                                 vmax=df['pct_nz_group'].max())

# 计算点大小 - 确保所有值都大于0
sizes = np.clip(min_size + (max_size - min_size) * 
                (df['pct_nz_group'] - size_norm.vmin) / 
                (size_norm.vmax - size_norm.vmin + 1e-10), 
                min_size, max_size)  # 添加1e-10避免除以零

# 创建画布
plt.figure(figsize=(14, 10))

# 1. 主绘图区域
ax = plt.subplot2grid((10, 10), (0, 0), rowspan=8, colspan=7)

# 绘制散点图 - 使用有效的点大小
scatter = ax.scatter(
    x=df['Symbol'],
    y=df['group'],
    s=sizes,
    c=df['logfoldchanges'],
    cmap=cmap,
    vmin=-6,
    vmax=6,
    edgecolor='gray',
    alpha=0.8
)

# 设置标签和标题
ax.set_xlabel('Symbol', fontsize=12, fontweight='bold')
ax.set_ylabel('Group', fontsize=12, fontweight='bold')
ax.set_title('Gene Expression Dot Plot', fontsize=14, fontweight='bold')
plt.xticks(rotation=45, ha='right')

# 2. 点大小图例 - 右侧上方
ax_size = plt.subplot2grid((10, 10), (0, 7), rowspan=4, colspan=3)

# 创建点大小图例值 (使用NumPy数组)
size_values = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
legend_labels = [f"{val:.1f}" for val in size_values]
legend_sizes = min_size + (max_size - min_size) * (size_values - size_norm.vmin) / (size_norm.vmax - size_norm.vmin + 1e-10)
legend_sizes = np.clip(legend_sizes, min_size, max_size)  # 确保值有效

# 创建点大小图例
for i, size in enumerate(legend_sizes):
    ax_size.scatter([], [], s=size, color='gray', label=legend_labels[i])

# 设置图例位置和格式
ax_size.legend(title='pct_nz_group', 
               loc='center left',
               bbox_to_anchor=(0.5, 0.5),
               frameon=False,
               borderpad=1.5,
               title_fontsize=10,
               labelspacing=1.8)
ax_size.axis('off')
ax_size.set_title('Point Size Legend', fontsize=10, pad=10)

# 3. 颜色条 - 底部单独区域
cax = plt.subplot2grid((10, 10), (8, 0), rowspan=2, colspan=7)
cbar = plt.colorbar(scatter, cax=cax, orientation='horizontal')
cbar.set_label('logfoldchange', fontsize=10, labelpad=5)
cbar.set_ticks([-6, -4, -2, 0, 2, 4, 6])

# 调整布局
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.subplots_adjust(wspace=0.5, hspace=0.8)

plt.show()
plt.savefig('gene_expression_dotplot.png', dpi=300, bbox_inches='tight')
plt.savefig('gene_expression_dotplot.pdf', bbox_inches='tight')

