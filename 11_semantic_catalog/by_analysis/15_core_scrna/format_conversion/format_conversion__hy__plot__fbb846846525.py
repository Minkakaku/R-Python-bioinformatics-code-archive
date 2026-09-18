import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap



# 1. 读取h5ad数据并加载聚类信息
adata = sc.read_h5ad("adata_celltype.h5ad")  # 替换为你的h5ad文件路径
print(adata)


target_genes = ["APOL2","GBP1","JAK1","JAK2","JAK3","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","WARS","IFNG","IFNGR1","IFNGR2"]  # 替换为目标基因
valid_genes = [gene for gene in target_genes if gene in adata.var_names]

n_genes = len(valid_genes)
n_cols = 3  # 每行3个子图
n_rows = (n_genes + n_cols - 1) // n_cols
fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols*5, n_rows*4))
axs = axs.flatten()

# 创建自定义颜色映射（灰白→红）
colors = ["#f0f0f0", "#ff0000"]  # 灰白色到红色
custom_cmap = LinearSegmentedColormap.from_list("custom_red", colors)

# 循环绘制每个基因的UMAP图
for i, gene in enumerate(valid_genes):
    sc.pl.umap(
        adata,
        color=gene,
        ax=axs[i],
        title=gene,
        show=False,
        color_map=custom_cmap,  # 颜色映射方案[4](@ref)
        frameon=False,
        legend_loc='on data'
    )
    
# 隐藏空白子图
for j in range(n_genes, len(axs)):
    axs[j].axis('off')

plt.tight_layout()
plt.savefig('gene_expression_umap.png', dpi=300)
plt.close()

if 'group' not in adata.obs.columns:
    raise KeyError("未找到分组信息 'group'，请检查adata.obs")

# 绘制分组的箱线图或小提琴图
fig, axs = plt.subplots(n_genes, 1, figsize=(8, 4*n_genes))

for i,gene in enumerate(valid_genes):
    png= str(gene)+"_violin_plot.png"
    # 箱线图（适合组间比较）
    # 绘制分组小提琴图，按 'group' 和 'celltype' 分组
    sc.pl.violin(
        adata,
        keys=gene,
        groupby='group',
        split=True,       # 分割不同组的小提琴图
        palette='muted',  # 颜色方案
        inner='stick',    # 内部显示数据点
        rotation=45,      # 标签旋转角度
        legend=False
    )
    plt.suptitle('Gene Expression by Group and Cell Type', y=1.02)  # 全局标题
    plt.tight_layout()
    plt.savefig(png, dpi=300)
    # 或使用箱线图（替换为以下代码）：
    # sc.pl.rank_genes_groups_violin(adata, groups=None, n_genes=1, gene_symbols=gene, groupby='group', ax=axs[i])

plt.tight_layout()
plt.savefig('gene_expression_by_group.png', dpi=300)
plt.close()
