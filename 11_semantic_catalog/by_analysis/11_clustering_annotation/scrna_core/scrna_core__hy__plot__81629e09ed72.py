import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# 1. 读取数据
df = pd.read_csv('data.csv')

# 2. 数据透视
df_pivot = df.pivot_table(
    index='names',
    columns='comparison',
    values='logfoldchanges',
    aggfunc='mean'
)

# 3. 指定列顺序
desired_order = ['A_vs_N', 'C_vs_A', 'C_vs_N']
df_pivot = df_pivot.reindex(columns=desired_order)

# 4. 填充缺失值为 0
df_pivot = df_pivot.fillna(0)

# 5. 计算颜色范围
max_abs = df_pivot.abs().max().max()
if max_abs < 0.1:
    max_abs = 2

# 6. 自定义颜色映射
colors = ["navy", "white", "#CD2626"]  # firebrick3 的十六进制
cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)

# 7. 绘制热图
plt.figure(figsize=(8, len(df_pivot.index) * 0.5 + 2))

ax = sns.heatmap(
    df_pivot,
    cmap=cmap,
    vmin=-max_abs,
    vmax=max_abs,
    center=0,
    annot=False,
    linewidths=0.5,
    linecolor='gray',
    cbar_kws={'label': 'log2 Fold Change'}
)

plt.title('Genes in different group', fontsize=14, pad=15)
plt.ylabel('Gene', fontsize=12)
plt.xlabel('Comparison', fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig('functional_evolution_heatmap.pdf', dpi=300, bbox_inches='tight')
plt.savefig('functional_evolution_heatmap.tiff', dpi=300, bbox_inches='tight')  # 删除了 compression

print("✅ 热图已保存")
plt.show()
