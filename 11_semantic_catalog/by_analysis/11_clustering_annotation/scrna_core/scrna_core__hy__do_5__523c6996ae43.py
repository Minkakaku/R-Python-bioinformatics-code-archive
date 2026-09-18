import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# 创建示例数据（根据您提供的表格）
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

# 创建自定义颜色映射：负值蓝色渐变，0为白色，正值红色渐变
colors = [(0, 'darkblue'), (0.4, 'lightblue'), (0.5, 'white'), (0.6, 'lightcoral'), (1, 'darkred')]
cmap = LinearSegmentedColormap.from_list('custom_red_blue', colors)

# 设置点的大小范围（根据pct_nz_group值）
min_size = 50
max_size = 1000
sizes = min_size + (max_size - min_size) * (df['pct_nz_group'] - df['pct_nz_group'].min()) / (df['pct_nz_group'].max() - df['pct_nz_group'].min())

# 创建图表
plt.figure(figsize=(12, 8))

# 绘制散点图
scatter = plt.scatter(
    x=df['Symbol'],
    y=df['group'],
    s=sizes,
    c=df['logfoldchanges'],
    cmap=cmap,
    vmin=-6,  # 根据数据范围设置
    vmax=6,   # 根据数据范围设置
    edgecolor='gray',
    alpha=0.8
)

# 添加颜色条
cbar = plt.colorbar(scatter, label='logfoldchange')
cbar.set_ticks([-6, -3, 0, 3, 6])  # 设置刻度位置

# 添加标签和标题
plt.xlabel('Symbol', fontsize=12, fontweight='bold')
plt.ylabel('Group', fontsize=12, fontweight='bold')
plt.title('Gene Expression Dot Plot', fontsize=14, fontweight='bold')

# 调整布局
plt.xticks(rotation=45, ha='right')  # 旋转x轴标签
plt.tight_layout()

# 添加图例说明点的大小
for size in [0.2, 0.4, 0.6]:
    plt.scatter([], [], s=min_size + (max_size - min_size)*size, 
                c='k', alpha=0.5, label=f'{size:.1f}')
plt.legend(title='pct_nz_group', loc='center left', bbox_to_anchor=(1.05, 0.5), 
           frameon=False, labelspacing=1.5)

plt.show()
plt.savefig('gene_expression_dotplot.png', dpi=300, bbox_inches='tight')
plt.savefig('gene_expression_dotplot.pdf', bbox_inches='tight')
