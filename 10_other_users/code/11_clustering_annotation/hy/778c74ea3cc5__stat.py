import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 设置中文字体（可选，避免中文显示问题）
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'DejaVu Sans']  # 根据系统选择
plt.rcParams['axes.unicode_minus'] = False

# 1. 读取h5ad文件（替换为你的文件路径）
adata = sc.read_h5ad("adata_celltype_stage.h5ad")  # 修改为实际路径

# 2. 从adata.obs中提取样本、阶段和细胞类型信息
df = adata.obs[['sample', 'stage', 'celltype']].copy()

# 3. 计算每个(sample, stage, celltype)组合的细胞数量
count_df = df.groupby(['sample', 'stage', 'celltype']).size().reset_index(name='count')

# 4. 计算每个(sample, stage)分组内各细胞类型的比例
# 先计算每个(sample, stage)分组的总细胞数
count_df['total_per_group'] = count_df.groupby(['sample', 'stage'])['count'].transform('sum')
# 再计算比例
count_df['proportion'] = count_df['count'] / count_df['total_per_group']

# 5. 检查数据（可选）
print("比例计算结果前5行：")
print(count_df.head())

# 6. 创建堆积条形图
# 方法：使用Seaborn的catplot，通过分面（facet）显示不同stage
# 调整数据格式：每个样本-阶段-细胞类型组合占一行
g = sns.catplot(
    data=count_df,
    kind='bar',
    x='sample',          # x轴为样本
    y='proportion',       # y轴为比例
    hue='celltype',       # 颜色区分细胞类型
    col='stage',          # 按阶段分面（每个阶段一个子图）
    dodge=False,          # 堆积图不需并置
    palette='Set3',       # 颜色方案
    legend_out=True,      # 图例放在外部
    height=5,             # 每个子图高度
    aspect=1.5            # 宽高比
)

# 7. 自定义图表
g.set_axis_labels("Sample", "Cell Type Proportion")
g.set_xticklabels(rotation=45)  # 旋转x轴标签
g.fig.suptitle("Cell Type Proportions by Sample and Stage", y=1.02)  # 总标题
plt.tight_layout()

# 8. 保存图表（可选）
plt.savefig("celltype_proportions.png", dpi=300, bbox_inches='tight')

# 9. 显示图表
plt.show()

# 可选：输出比例表格
proportion_table = count_df.pivot_table(
    index=['sample', 'stage'], 
    columns='celltype', 
    values='proportion', 
    fill_value=0
)
print("\n比例表格（样本×阶段×细胞类型）：")
print(proportion_table)
