"""
生成例图样式的热图：
- 左侧为热图主体（基因行×细胞列），不显示左侧y轴基因名
- 右侧用彩色方框展示每个细胞类型的 top 基因列表（允许重复基因）
- 底部显示细胞类型名称及对应颜色条
"""

import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import pandas as pd
import textwrap
import warnings
warnings.filterwarnings('ignore')

# ========== 参数设置 ==========
N_GENES = 12  # 每个细胞类型的 top 基因数
out_dir = os.path.dirname("/mnt/NexData_file/home/sh/Projects/FWSC20241287/14.Mission/heatmap_mac/")
h5ad_path = os.path.join("/mnt/NexData_file/home/sh/Projects/FWSC20241287/05.Mission/result_recluster_Macrophages/recluster_Macrophages.h5ad")

# ========== 读取数据 ==========
adata = sc.read_h5ad(h5ad_path)
adata.obs['cell_type'] = adata.obs['cell_type'].cat.add_categories('Remodeling Mac')
mask = adata.obs['cell_type'] == 'Remodeling Mac-Svep1'
adata.obs.loc[mask, 'cell_type'] = 'Remodeling Mac'
adata.obs['cell_type'] = adata.obs['cell_type'].cat.remove_unused_categories()

# 按 cell_type 计算差异基因
sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)

# ========== 获取每个细胞类型的 top 基因 ==========
cell_types = list(adata.obs['cell_type'].cat.categories)
top_genes_dict = {}  # {cell_type: [gene1, gene2, ...]}

for ct in cell_types:
    tmp = sc.get.rank_genes_groups_df(adata, group=ct)
    tmp = tmp[(tmp['logfoldchanges']>1) & (tmp['pvals_adj']<0.05) & (tmp['pct_nz_group'] >= 0.25)]
    genes = list(tmp.sort_values('logfoldchanges', ascending=False)['names'].head(N_GENES))
    top_genes_dict[ct] = genes

# 构建完整基因列表（允许重复，每组固定 N_GENES 行）
all_genes_rows = []  # [(gene_name, cell_type), ...]
for ct in cell_types:
    for g in top_genes_dict[ct]:
        all_genes_rows.append((g, ct))

gene_names_list = [row[0] for row in all_genes_rows]
n_total_rows = len(all_genes_rows)

# ========== 准备热图数据 ==========
# 按 cell_type 排序细胞
cell_order = []
cell_type_sizes = {}
for ct in cell_types:
    mask = adata.obs['cell_type'] == ct
    cells = adata.obs.index[mask].tolist()
    cell_order.extend(cells)
    cell_type_sizes[ct] = len(cells)

adata_ordered = adata[cell_order, :]

# 提取表达矩阵：行=基因（允许重复），列=细胞
available_genes = set(adata_ordered.var_names)
expr_matrix = np.zeros((n_total_rows, len(cell_order)))

for i, (gene, ct) in enumerate(all_genes_rows):
    if gene in available_genes:
        row_data = adata_ordered[:, gene].X
        if hasattr(row_data, 'toarray'):
            row_data = row_data.toarray()
        expr_matrix[i, :] = row_data.flatten()

# 按基因行标准化 (0-1缩放)
for i in range(n_total_rows):
    vmin = expr_matrix[i, :].min()
    vmax = expr_matrix[i, :].max()
    if vmax > vmin:
        expr_matrix[i, :] = (expr_matrix[i, :] - vmin) / (vmax - vmin)
    else:
        expr_matrix[i, :] = 0

# ========== 定义颜色方案 ==========
# 使用与例图类似的柔和颜色
pastel_colors = [
    '#1f77b4',
    '#ff7f0e',
    '#279e68',
    '#d62728',
    '#aa40fc',
    '#8c564b',
    '#e377c2',
    '#b5bd61'
]

ct_colors = {ct: pastel_colors[i % len(pastel_colors)] for i, ct in enumerate(cell_types)}

# ========== 绘图 ==========
fig_height = max(7, n_total_rows * 0.1 + 2)
fig = plt.figure(figsize=(12, fig_height))

# GridSpec 布局
gs = GridSpec(2, 2, width_ratios=[3, 1.5], height_ratios=[20, 1],
              wspace=0.01, hspace=0.08, figure=fig)

ax_heatmap = fig.add_subplot(gs[0, 0])
ax_genes = fig.add_subplot(gs[0, 1])
ax_colorbar_bottom = fig.add_subplot(gs[1, 0])

# ---------- 绘制热图 ----------
im = ax_heatmap.imshow(expr_matrix, aspect='auto', cmap='viridis',
                        interpolation='none', vmin=0, vmax=1)
# im = ax_heatmap.imshow(expr_matrix, aspect='auto', cmap='Reds',
                        # interpolation='none', vmin=0, vmax=1)
# 细胞类型垂直分隔线
pos = 0
for ct in cell_types:
    n_cells = cell_type_sizes[ct]
    pos += n_cells
    if pos < len(cell_order):
        ax_heatmap.axvline(x=pos - 0.5, color='black', linewidth=0.3, alpha=0.4)

# 基因组水平分隔线
row_pos = 0
for ct in cell_types:
    row_pos += N_GENES
    if row_pos < n_total_rows:
        ax_heatmap.axhline(y=row_pos - 0.5, color='black', linewidth=0.3, alpha=0.4)

# 右侧颜色条（标识基因所属细胞类型）
# 创建右侧的窄颜色条
right_bar_width = len(cell_order) * 0.012
for i, (gene, ct) in enumerate(all_genes_rows):
    color = ct_colors[ct]
    ax_heatmap.add_patch(mpatches.Rectangle(
        (len(cell_order) - 0.5, i - 0.5), right_bar_width, 1,
        facecolor=color, edgecolor='none', clip_on=False
    ))

# 隐藏y轴基因标签（例图不在左侧显示）
ax_heatmap.set_yticks([])
ax_heatmap.set_xticks([])
ax_heatmap.set_ylabel('')
ax_heatmap.set_xlabel('')

for spine in ax_heatmap.spines.values():
    spine.set_linewidth(0.5)
    spine.set_color('black')

# ---------- 右侧基因注释面板 ----------
ax_genes.set_xlim(0, 1)
ax_genes.set_ylim(n_total_rows, 0)
ax_genes.axis('off')

# 为每个细胞类型画彩色框+基因列表
row_pos = 0
for ct in cell_types:
    n_genes_ct = N_GENES
    color = ct_colors[ct]
    
    # 绘制圆角矩形框
    rect = mpatches.FancyBboxPatch(
        (0.01, row_pos + 0.15), 0.97, n_genes_ct - 0.3,
        boxstyle="round,pad=0.08",
        facecolor=color, alpha=0.3,
        edgecolor=color, linewidth=1.5
    )
    ax_genes.add_patch(rect)
    
    # 基因文本
    genes_text = ', '.join(top_genes_dict[ct])
    wrapped_text = textwrap.fill(genes_text, width=42)
    
    y_center = row_pos + n_genes_ct / 2
    
    ax_genes.text(0.5, y_center, wrapped_text,
                  ha='center', va='center',
                  fontsize=7.5, fontweight='normal',
                  color='#222222',
                  fontfamily='sans-serif',
                  linespacing=1.3)
    
    row_pos += n_genes_ct

# ---------- 底部细胞类型颜色条 ----------
ax_colorbar_bottom.set_xlim(0, len(cell_order))
ax_colorbar_bottom.set_ylim(0, 1)

pos = 0
for ct in cell_types:
    n_cells = cell_type_sizes[ct]
    color = ct_colors[ct]
    rect = mpatches.Rectangle((pos, 0), n_cells, 1,
                                facecolor=color, edgecolor='none', alpha=0.9)
    ax_colorbar_bottom.add_patch(rect)
    ax_colorbar_bottom.text(pos + n_cells/2, -0.5, ct,
                            ha='center', va='top',
                            fontsize=6.5, rotation=45,
                            fontfamily='sans-serif')
    pos += n_cells

ax_colorbar_bottom.set_xticks([])
ax_colorbar_bottom.set_yticks([])
for spine in ax_colorbar_bottom.spines.values():
    spine.set_visible(False)

ax_colorbar_bottom.text(len(cell_order)/2, -3.5, 'cell_type',
                        ha='center', va='top', fontsize=9,
                        fontfamily='sans-serif')

# 色标 - 动态定位紧贴热图左侧
fig.canvas.draw()  # 先渲染以获取正确的坐标
heatmap_bbox = ax_heatmap.get_position()
cbar_width = 0.012
cbar_height = heatmap_bbox.height * 0.4
cbar_x = heatmap_bbox.x0 - cbar_width - 0.04  # 紧贴热图左侧，留出刻度空间
cbar_y = heatmap_bbox.y0 + (heatmap_bbox.height - cbar_height) / 2
cbar_ax = fig.add_axes([cbar_x, cbar_y, cbar_width, cbar_height])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.ax.yaxis.set_ticks_position('left')
cbar.ax.yaxis.set_label_position('left')
cbar.ax.tick_params(labelsize=6.5)
cbar.outline.set_linewidth(0.5)

# ========== 保存 ==========
plt.savefig(os.path.join(out_dir, 'heatmap_styled.png'),
            dpi=200, bbox_inches='tight', pad_inches=0.2)
plt.savefig(os.path.join(out_dir, 'heatmap_styled.pdf'),
            bbox_inches='tight', pad_inches=0.2)
plt.close()
print(f"热图已保存: heatmap_styled.png / .pdf")
print(f"共 {len(cell_types)} 个细胞类型, 每组 {N_GENES} 个基因, 总计 {n_total_rows} 行")
