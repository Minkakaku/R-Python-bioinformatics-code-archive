import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# 从您提供的图片信息中提取精确数据 - 保持原始顺序
data = {
    'Pathway': [
        'Endocytosis',
        'Regulation Of Interleukin-6 Production',
        'Protein Phosphorylation',
        'Regulation Of Inflammatory Response',
        'Positive Regulation Of Leukocyte Cell-Cell Adhesion',
        'Regulation Of Tumor Necrosis Factor Production',
        'Receptor-Mediated Endocytosis',
        'Cellular Response To Cytokine Stimulus',
        'Receptor Internalization',
        'Cytokine-Mediated Signaling Pathway'
    ],
    'Gene_number': [10, 12, 10, 21, 12, 8, 12, 18, 9, 18],
    'P_value': [1.11e-06, 1.11e-06, 4.68e-07, 1.35e-07, 1.41e-08, 
                2.14e-08, 8.94e-09, 7.88e-09, 3.62e-09, 4.44e-10],
    'ratio': [0.063492063, 0.063492063, 0.087719298, 0.042, 0.094488189,
              0.181818182, 0.098360656, 0.058441558, 0.176470588, 0.070038911]
}

# 创建DataFrame - 保留原始顺序
df = pd.DataFrame(data)

# 计算负对数P值
df['-log10(P_value)'] = -np.log10(df['P_value'])

# 创建自定义从蓝到红的渐变色（低P值→红色，高P值→蓝色）
colors = [(0.10, 0.30, 0.85), (0.85, 0.10, 0.10)]  # 蓝 → 红
cmap = LinearSegmentedColormap.from_list('BluRed_Ramp', colors)

# 创建图形
plt.figure(figsize=(14, 10))
ax = plt.gca()

# 绘制散点图
scatter = ax.scatter(
    x='ratio', 
    y='Pathway',
    s=df['Gene_number'] * 30,  # 点大小基于基因数量
    c=df['-log10(P_value)'],   # 点颜色基于-log10(P值)
    cmap=cmap,
    edgecolor='k',             # 黑色边缘
    linewidth=0.8,
    alpha=0.85,
    zorder=3,                  # 确保点在网格上方
    data=df
)

# 添加颜色条（修复get_clim错误）
cbar = plt.colorbar(scatter, pad=0.02, shrink=0.3)
cbar.set_label('-log10(P value)', fontsize=12, labelpad=10)
cbar.ax.tick_params(labelsize=10)

# 添加"High"文本标记（使用正确的坐标）
vmin, vmax = cbar.mappable.get_clim()  # 获取颜色映射范围

# 设置图形标题和标签
plt.title('Pathway Enrichment Analysis', fontsize=16, fontweight='bold', pad=20)
plt.xlabel('Gene Ratio', fontsize=13, labelpad=15)
plt.ylabel('Pathway', fontsize=13, labelpad=15)
plt.grid(axis='x', linestyle='--', alpha=0.6, zorder=1)

# 设置刻度标签大小
plt.xticks(fontsize=11)
plt.yticks(fontsize=15,fontweight='bold')

# 创建分离良好的图例（修复重叠问题）
gene_sizes = [5,10,21]
legend_elements = []

# 创建间距合理的图例项
for size in gene_sizes:
    # 使用空白标签创建间距
    plt.scatter([], [], alpha=0, s=0, 
                label=f'\n ')  # 添加空行增加垂直间距
    
    # 添加实际点项
    scatter_point = plt.scatter([], [], s=size *20, 
                               edgecolors='k', color='gray', 
                               alpha=0.85, label=f'{size} genes')
    legend_elements.append(scatter_point)

# 创建图例（优化间距）
legend = ax.legend(handles=legend_elements, 
                   loc='upper right', 
                   bbox_to_anchor=(1.35, 1.00),  # 调整位置避免重叠
                   title='Gene Number', 
                   frameon=True, 
                   framealpha=0.9, 
                   title_fontsize=10,
                   handletextpad=4.0,  # 增加点与文本间距
                   borderpad=2.0,      # 增加图例外边距
                   labelspacing=1.2,   # 增加标签间距
                   columnspacing=1.5,  # 增加列间距
                   handlelength=1.0,   # 增加点长度
                   ncol=1)             # 单列布局更清晰

# 优化布局
plt.tight_layout()
plt.subplots_adjust(left=0.25, right=0.78, top=0.92)

plt.show()
plt.savefig("PATHWAY.pdf", dpi=300, bbox_inches='tight')
