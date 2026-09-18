"""
绘制单细胞热图：每种细胞类型的 Top10 差异基因
右侧标注经典 marker 基因
风格参照 例图.png
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')
import textwrap

# ============================================================
# 1. 经典 Marker 基因定义
# ============================================================
# CLASSIC_MARKERS = {
#     'Macrophage 1': [
#         'Ly6c2','Spp1','Ccl12','Ms4a4c','Spon1','Ccl2','Phf11b','Scimp','Plac8','Treml2',
#     ],
#     'Macrophage 2': [
#         'Airn','Cstdc5','Krt14','Krt6a','Pf4','Krt6b','Ampd3','Krt5','Krt16','Cst6',
#     ],
#     'Macrophage 3': [
#         'C6','Chp2','Mgl2','Cd209d','Cd209g','Cd209f','Lyve1','Retnla','Sult1a1','Fcna',
#     ],
#     'Macrophage 4': [
#         'Cxcl3','Inhba','Il1b','Arg1','Ptgs2','AA467197','Il1a','Serpinb2','Nos2','Thbs1',
#     ],
#     'Macrophage 5': [
#         'Ly6a2','Chil3','Hp','Lrg1','Napsa','Wfdc21','Sorl1','Lcn2','Il1r2','Ifitm6',
#     ],
#     'Macrophage 6': [
#         'Cdk8','Gphn','Cmss1','Slfn4','Camk1d','Lars2','Retnlg','Gm19951','S100a9','mt-Nd5'
#     ]
# }

# # 细胞类型缩写
# CELL_TYPE_ABBREV = {
#     'Macrophage 1': 'Mac 1',
#     'Macrophage 2': 'Mac 2',
#     'Macrophage 3': 'Mac 3',
#     'Macrophage 4': 'Mac 4',
#     'Macrophage 5': 'Mac 5',
#     'Macrophage 6': 'Mac 6',
# }

# # 固定细胞类型顺序
# CELL_TYPE_ORDER = [
#     'Macrophage 1',
#     'Macrophage 2',
#     'Macrophage 3',
#     'Macrophage 4',
#     'Macrophage 5',
#     'Macrophage 6'
# ]

# # 为每种细胞类型分配颜色（参考例图风格）
# CELL_TYPE_COLORS = {
#     'Macrophage 1': '#1f77b4',
#     'Macrophage 2': '#ff7f0e',
#     'Macrophage 3': '#279e68',
#     'Macrophage 4': '#d62728',
#     'Macrophage 5': '#aa40fc',
#     'Macrophage 6': '#8c564b'
# }
CLASSIC_MARKERS = {
    'Fibroblast 1': [
        'Limch1','Smpd3','Efhd1','Dact2','Negr1','Sema3e','Car8','Igfbp5','Stmn4','Ppp2r2b',
    ],
    'Fibroblast 2': [
        'Prg4','Ifi206','Cxcl5','Mx1','Ifi209','Apol9b','Irf7','Ifit3b','Ccl5','Sh2d6',
    ],
    'Fibroblast 3': [
        'Angptl1','Col15a1','Hmcn2','Mfap4','Cilp','Gas6','Sparcl1','Cygb','Mgp','Cpxm2',
    ],
    'Fibroblast 4': [
        'Mirt1','Pik3ap1','Stfa2l1','Samsn1','Acod1','Ccl3','Ccr1','S100a9','Hcar2','Il36g',
    ],
    'Fibroblast 5': [
        'Prss35','Robo2','Megf6','Adgrl3','Grem1','Tnn','Kif26b','Brinp3','Fhod3','Col7a1',
    ]
}

# 细胞类型缩写
CELL_TYPE_ABBREV = {
    'Fibroblast 1': 'Fib 1',
    'Fibroblast 2': 'Fib 2',
    'Fibroblast 3': 'Fib 3',
    'Fibroblast 4': 'Fib 4',
    'Fibroblast 5': 'Fib 5'
}

# 固定细胞类型顺序
CELL_TYPE_ORDER = [
    'Fibroblast 1',
    'Fibroblast 2',
    'Fibroblast 3',
    'Fibroblast 4',
    'Fibroblast 5'
]

# 为每种细胞类型分配颜色（参考例图风格）
CELL_TYPE_COLORS = {
    'Fibroblast 1': '#1f77b4',
    'Fibroblast 2': '#ff7f0e',
    'Fibroblast 3': '#279e68',
    'Fibroblast 4': '#d62728',
    'Fibroblast 5': '#aa40fc'
}

# ============================================================
# 2. 加载数据并计算差异基因
# ============================================================
def load_and_compute(h5ad_path):
    print("Loading data...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"Data shape: {adata.shape}")

    need_compute = True
    if 'rank_genes_groups' in adata.uns:
        if adata.uns['rank_genes_groups']['params']['groupby'] == 'cell_type':
            need_compute = False
            print("Using existing rank_genes_groups (cell_type).")

    if need_compute:
        print("Computing rank_genes_groups by cell_type (wilcoxon)...")
        sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon',
                                use_raw=False, n_genes=50)
        print("Done.")

    return adata


# ============================================================
# 3. 提取 top10 基因并匹配经典 marker
# ============================================================
def get_top_genes_and_markers(adata, n_top=10):
    result = adata.uns['rank_genes_groups']
    top_genes_dict = {}
    marker_info_dict = {}

    for ct in CELL_TYPE_ORDER:
        if ct not in result['names'].dtype.names:
            print(f"Warning: {ct} not in rank_genes_groups, skipping.")
            continue
        names = result['names'][ct]
        filtered = [g for g in names if not g.startswith('ENSOCUG') and not g.startswith('RA_')]
        top10 = filtered[:n_top]
        top_genes_dict[ct] = top10

        classic = CLASSIC_MARKERS.get(ct, [])
        matched = [g for g in top10 if g in classic]
        marker_info_dict[ct] = matched
        print(f"  {ct}: top10 = {top10}")
        if matched:
            print(f"    -> Classic markers found: {matched}")

    return top_genes_dict, marker_info_dict


# ============================================================
# 4. 绘制热图（仿照例图风格）
# ============================================================
def plot_heatmap(adata, top_genes_dict, marker_info_dict, output_prefix):

    # --- 收集基因列表 ---
    gene_list = []
    gene_to_ct = {}
    ct_gene_ranges = {}
    idx = 0
    for ct in CELL_TYPE_ORDER:
        if ct not in top_genes_dict:
            continue
        start = idx
        for g in top_genes_dict[ct]:
            if g in adata.var_names and g not in gene_list:
                gene_list.append(g)
                gene_to_ct[g] = ct
                idx += 1
        end = idx
        ct_gene_ranges[ct] = (start, end)

    print(f"Plotting {len(gene_list)} genes...")

    # --- 按 cell_type 排序细胞 ---
    adata.obs['cell_type'] = pd.Categorical(
        adata.obs['cell_type'], categories=CELL_TYPE_ORDER, ordered=True
    )
    sorted_idx = adata.obs.sort_values('cell_type').index
    adata_sorted = adata[sorted_idx, gene_list].copy()

    sc.pp.scale(adata_sorted, max_value=3)

    if hasattr(adata_sorted.X, 'toarray'):
        expr_matrix = adata_sorted.X.toarray()
    else:
        expr_matrix = np.array(adata_sorted.X)

    # --- 降采样 ---
    n_cells_total = expr_matrix.shape[0]
    max_cells = 1000000
    if n_cells_total > max_cells:
        print(f"Downsampling from {n_cells_total} to ~{max_cells} cells...")
        cell_types_sorted = adata_sorted.obs['cell_type'].values
        sampled_indices = []
        for ct in CELL_TYPE_ORDER:
            ct_indices = np.where(cell_types_sorted == ct)[0]
            n_sample = max(10, int(max_cells * len(ct_indices) / n_cells_total))
            if len(ct_indices) > n_sample:
                sampled = np.sort(np.random.choice(ct_indices, n_sample, replace=False))
            else:
                sampled = ct_indices
            sampled_indices.append(sampled)
        sampled_indices = np.concatenate(sampled_indices)
        expr_matrix = expr_matrix[sampled_indices]
        cell_types_plot = cell_types_sorted[sampled_indices]
    else:
        cell_types_plot = adata_sorted.obs['cell_type'].values

    expr_matrix = np.clip(expr_matrix, -2, 2)
    expr_matrix = expr_matrix.T  # 行=基因, 列=细胞

    print(f"Final matrix shape: {expr_matrix.shape}")

    # ================================================================
    #  自定义 colormap：仿照例图 (浅紫/薰衣草 → 白 → 浅红 → 深红)
    #  底部 colorbar 显示 -2~2，低端蓝/高端红
    # ================================================================
    cmap_colors = [
        (0.0,  '#1A237E'),   # 深蓝
        (0.15, '#5C6BC0'),   # 蓝
        (0.30, '#B5C7E3'),   # 浅蓝
        (0.42, '#DDD6E8'),   # 薰衣草
        (0.50, '#F2EFF5'),   # 近白（微紫）
        (0.58, '#F5E0D8'),   # 很浅粉
        (0.70, '#E8A88C'),   # 浅橙红
        (0.82, '#C9604A'),   # 红
        (0.92, '#A0322A'),   # 深红
        (1.0,  '#7B1A1A'),   # 暗红
    ]
    cmap = LinearSegmentedColormap.from_list(
        'custom_heatmap',
        [(pos, col) for pos, col in cmap_colors],
        N=256
    )

    # ================================================================
    #  绑制图形
    # ================================================================
    n_genes = len(gene_list)
    fig_height = max(8, n_genes * 0.09 + 1.5)
    fig = plt.figure(figsize=(12, fig_height))

    # GridSpec: 热图 + 底部色带 + 右侧标注
    gs = gridspec.GridSpec(
        2, 2,
        width_ratios=[5, 3.0],
        height_ratios=[n_genes, 0.8],
        wspace=0.01, hspace=0.02
    )

    # ---- 主热图 ----
    ax_heat = fig.add_subplot(gs[0, 0])
    im = ax_heat.imshow(
        expr_matrix, aspect='auto', cmap=cmap,
        vmin=-2, vmax=2, interpolation='nearest'
    )
    ax_heat.set_xticks([])
    ax_heat.set_yticks([])
    ax_heat.set_ylabel('Genes', fontsize=12, labelpad=8)

    # ---- 右侧细胞类型色带（新增）----
    # 获取热图位置以精确定位右侧色带
    # heatmap_bbox = ax_heat.get_position()

    # # 设置右侧色带参数
    # right_bar_width = 0.5  # 色带宽度（可调整范围：0.05 到 10，根据实际效果选择，0.05-0.2通常合适）
    # right_bar_height = heatmap_bbox.height  # 色带高度与热图保持一致
    # right_bar_x = heatmap_bbox.x1 + 0.05  # 紧贴热图右侧（间隙可调）
    # right_bar_y = heatmap_bbox.y0  # 垂直对齐（顶部对齐）
    # ax_right_cbar = fig.add_axes([right_bar_x, right_bar_y, right_bar_width, right_bar_height])

    # # 准备颜色数据（垂直方向）
    # ct_colors_column = np.array([
    #     matplotlib.colors.to_rgb(CELL_TYPE_COLORS.get(ct, '#CCCCCC'))
    #     for ct in cell_types_plot
    # ]).reshape(-1, 1, 3)

    # ax_right_cbar.imshow(ct_colors_column, aspect='auto', interpolation='nearest')
    # ax_right_cbar.set_xticks([])
    # ax_right_cbar.set_yticks([])

    # 添加细胞类型分界线（细小的白线）
    cell_boundaries = []
    current_ct = cell_types_plot[0]
    for i, ct in enumerate(cell_types_plot):
        if ct != current_ct:
            cell_boundaries.append(i)
            current_ct = ct
    for b in cell_boundaries:
        ax_heat.axvline(x=b - 0.5, color='white', linewidth=0.5, alpha=0.8)

    # 添加基因组分界线
    for ct in CELL_TYPE_ORDER:
        if ct not in ct_gene_ranges:
            continue
        start, end = ct_gene_ranges[ct]
        if end < n_genes:
            ax_heat.axhline(y=end - 0.5, color='white', linewidth=0.5, alpha=0.8)

    # ---- 右侧标注区域 ----
    ax_annot = fig.add_subplot(gs[0, 1])
    ax_annot.set_xlim(0, 10)
    ax_annot.set_ylim(n_genes, 0)
    ax_annot.axis('off')

    for ct in CELL_TYPE_ORDER:
        if ct not in ct_gene_ranges or ct not in marker_info_dict:
            continue
        start, end = ct_gene_ranges[ct]
        if start == end:
            continue

        mid_y = (start + end) / 2
        color = CELL_TYPE_COLORS.get(ct, '#999999')
        abbrev = CELL_TYPE_ABBREV.get(ct, ct[:3].upper())
        markers_found = marker_info_dict[ct]

        # # 三角形标记 (◁)
        # tri = plt.Polygon(
        #     [[0.2, mid_y - 0.6], [0.2, mid_y + 0.6], [0.7, mid_y]],
        #     color=color, zorder=5, clip_on=False
        # )
        # ax_annot.add_patch(tri)

        # 标注文本
        if markers_found:
            # 有经典 marker: 显示 "缩写: marker1, marker2, ..."
            marker_str = ', '.join(markers_found)
            label = f"{abbrev}: {marker_str}"
            wrapped_label = '\n'.join(textwrap.wrap(label, width=len(label)//2))
            ax_annot.text(
                0.9, mid_y, wrapped_label,
                fontsize=12, fontweight='bold', color=color,
                va='center', ha='left',
                clip_on=False
            )
        else:
            # 无经典 marker: 只显示缩写
            ax_annot.text(
                0.9, mid_y, abbrev,
                fontsize=7.5, fontweight='normal', color=color,
                va='center', ha='left',
                clip_on=False
            )

    # ---- 底部细胞类型色带 ----
    ax_cbar_ct = fig.add_subplot(gs[1, 0])
    ct_colors_row = np.array([
        matplotlib.colors.to_rgb(CELL_TYPE_COLORS.get(ct, '#CCCCCC'))
        for ct in cell_types_plot
    ]).reshape(1, -1, 3)
    ax_cbar_ct.imshow(ct_colors_row, aspect='auto', interpolation='nearest')
    ax_cbar_ct.set_xticks([])
    ax_cbar_ct.set_yticks([])

    # ---- 底部右侧: colorbar + 图例 ----
    ax_legend = fig.add_subplot(gs[1, 1])
    ax_legend.axis('off')

    # Colorbar（放在底部右侧区域的中上位置）
    cbar_ax = fig.add_axes([0.12, 0.015, 0.28, 0.012])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks([-2, -1, 0, 1, 2])
    cbar.ax.tick_params(labelsize=6)
    cbar.outline.set_linewidth(0.5)

    # 图例（放在 colorbar 右侧）
    legend_elements = []
    for ct in CELL_TYPE_ORDER:
        if ct in ct_gene_ranges:
            legend_elements.append(
                mpatches.Patch(
                    facecolor=CELL_TYPE_COLORS.get(ct, '#999'),
                    edgecolor='none', label=ct
                )
            )
    fig.legend(
        handles=legend_elements,
        loc='lower right', bbox_to_anchor=(0.95, 0.001),
        ncol=3, fontsize=5.5, frameon=False,
        handlelength=1.0, handleheight=0.8,
        columnspacing=0.5, labelspacing=0.3
    )

    plt.subplots_adjust(left=0.04, right=0.97, top=0.99, bottom=0.045)

    # 保存
    for fmt in ['png', 'pdf']:
        out = f"{output_prefix}.{fmt}"
        fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"Saved: {out}")
    plt.close(fig)


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    import os

    # data_dir = r"/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/Macrophages/"
    # h5ad_path = os.path.join(data_dir, 'Macrophages.h5ad')
    data_dir = r"/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/Fibroblasts/"
    h5ad_path = os.path.join(data_dir, 'Fibroblasts.h5ad')
    output_prefix = os.path.join(data_dir, 'heatmap_top10_markers')

    np.random.seed(42)

    adata = load_and_compute(h5ad_path)
    top_genes_dict, marker_info_dict = get_top_genes_and_markers(adata, n_top=10)
    print(top_genes_dict)
    plot_heatmap(adata, CLASSIC_MARKERS, CLASSIC_MARKERS, output_prefix)

    print("\nAll done!")
