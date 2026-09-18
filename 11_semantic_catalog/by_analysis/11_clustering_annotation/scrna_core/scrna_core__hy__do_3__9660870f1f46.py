import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker

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
plt.title("Gene Expression by Cell Type", fontsize=16, pad=16)
plt.xlabel("Genes", fontsize=14, labelpad=12)
plt.ylabel("Cell Types", fontsize=14, labelpad=12)
plt.xticks(rotation=45, ha='right')  # 旋转x轴标签
plt.grid(True, which='both', axis='y', linestyle='--', alpha=0.3)

# 获取图例句柄和标签
handles, labels = ax.get_legend_handles_labels()

# 移除默认生成的图例
if ax.legend_ is not None:
    ax.legend_.remove()

# 创建专门的点大小图例
# 获取唯一的点大小值并排序
unique_sizes = sorted(set(df['pct_nz_group']))
size_labels = [f"{size:.2f}" for size in unique_sizes]

# 创建对应的点大小句柄
size_handles = []
for size in unique_sizes:
    size_handles.append(
        plt.Line2D([0], [0], 
                  marker='o', 
                  color='w', 
                  markerfacecolor='gray',
                  markersize=np.sqrt(size * 800),  # 缩放点的大小
                  markeredgecolor='black')
    )

# 创建点大小图例，放在右下角
size_legend = ax.legend(
    size_handles,
    size_labels,
    labelspacing=1.5,
    borderpad=1.5,
    frameon=True,
    framealpha=0.9,
    title='pct_nz_group',
    loc='lower right',
    bbox_to_anchor=(1.05, 0.05)  # 放在右下角
)

# 添加颜色条（放在右上角）
norm = plt.Normalize(vmin, vmax)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cax = plt.colorbar(sm, ax=ax, 
                   location='right', 
                   pad=0.08, 
                   shrink=0.8,
                   aspect=30,
                   anchor=(0, 1.0))  # 固定在右上角

# 关键修复：使用正确的方法设置颜色条刻度
cax.ax.tick_params(labelsize=10)  # 设置刻度标签大小
cax.locator = ticker.MaxNLocator(nbins=6)  # 设置最大刻度数量为6
cax.update_ticks()  # 更新刻度

# 设置颜色条标签
cax.set_label('logfoldchanges', fontsize=12, labelpad=15)

# 调整布局
plt.tight_layout(pad=2.0)

# 保存高质量图片
plt.savefig('gene_expression_dotplot.png', dpi=300, bbox_inches='tight')
plt.savefig('gene_expression_dotplot.pdf', bbox_inches='tight')

# 显示图表
plt.show()
