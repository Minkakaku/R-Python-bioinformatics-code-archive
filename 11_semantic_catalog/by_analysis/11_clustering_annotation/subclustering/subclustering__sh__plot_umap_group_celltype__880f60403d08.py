"""
按组绘制UMAP图，每个分组子图中显示cell_type着色，
在整张图的左侧和底部添加贯穿全图的带箭头坐标轴。
输出 PNG 和 PDF。
"""

import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import matplotlib.transforms as mtransforms

# ============ 参数设置 ============
h5ad_path = r"/mnt/NexData_file/home/sh/Projects/FWSC20241287/05.Mission/result_recluster_Macrophages/recluster_Macrophages.h5ad"
out_dir = r"/mnt/NexData_file/home/sh/Projects/FWSC20241287/09.Mission/arrow_umap"
output_name = "dimplot_umap_group_celltype_Mac"

# ============ 读取数据 ============
adata = sc.read_h5ad(h5ad_path)
groups = adata.obs['group'].cat.categories.tolist()
n_groups = len(groups)

# ============ 绘图 ============
fig, axs = plt.subplots(nrows=1, ncols=n_groups, figsize=(n_groups * 5, 5))
if n_groups == 1:
    axs = [axs]

for i, grp in enumerate(groups):
    ax = axs[i]
    # 每个子图只显示该组的细胞
    sc.pl.umap(
        adata[adata.obs["group"] == grp],
        color="cell_type",
        legend_loc=None,       # 不在子图内显示图例
        title=grp,
        ax=ax,
        show=False,
        frameon=False,         # 去掉默认边框
        size=8,                # 点大小
    )
    ax.set_xlabel("")
    ax.set_ylabel("")

# 先做 tight_layout，确定最终布局
plt.tight_layout()

# ====== 在整张图上添加一组贯穿全图的带箭头坐标轴 ======
# 使用 fig.transFigure 坐标系（0,0 左下角，1,1 右上角）
# 箭头从左下角出发，x轴横跨整张图宽度，y轴纵跨整张图高度
arrow_color = 'black'
arrow_lw = 1.5
margin = 0.02  # 箭头起点离图边缘的距离（figure坐标比例）

# 起点
x0 = margin
y0 = margin

# x轴箭头：从左下角水平向右，跨整张图宽度
arrow_x = FancyArrowPatch(
    (x0, y0), (1.0 - margin, y0),
    arrowstyle='->', mutation_scale=15,
    color=arrow_color, linewidth=arrow_lw,
    clip_on=False,
    transform=fig.transFigure
)
fig.patches.append(arrow_x)

# y轴箭头：从左下角垂直向上，跨整张图高度
arrow_y = FancyArrowPatch(
    (x0, y0), (x0, 1.0 - margin),
    arrowstyle='->', mutation_scale=15,
    color=arrow_color, linewidth=arrow_lw,
    clip_on=False,
    transform=fig.transFigure
)
fig.patches.append(arrow_y)

# 添加标签（黑色、正常字体）
fig.text(
    0.5, margin * 0.3,
    "UMAP1", ha='center', va='top',
    fontsize=10, color=arrow_color,
    transform=fig.transFigure
)
fig.text(
    margin * 0.3, 0.5,
    "UMAP2", ha='right', va='center',
    fontsize=10, color=arrow_color,
    rotation=90,
    transform=fig.transFigure
)

# ============ 保存 ============
png_path = os.path.join(out_dir, f"{output_name}.png")
pdf_path = os.path.join(out_dir, f"{output_name}.pdf")
plt.savefig(png_path, facecolor='white', bbox_inches='tight', dpi=300, pad_inches=0.15)
plt.savefig(pdf_path, facecolor='white', bbox_inches='tight', pad_inches=0.15)
plt.close()
print(f"已保存: {png_path}")
print(f"已保存: {pdf_path}")
