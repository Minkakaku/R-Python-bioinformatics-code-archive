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

# ============================================================
# 1. 经典 Marker 基因定义
# ============================================================
CLASSIC_MARKERS = {
    'Fibroblasts': [
        'COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'PDGFRA', 'VIM',
        'FBN1', 'COL14A1', 'CCDC80', 'AEBP1', 'SERPINF1', 'MGP',
        'PDGFRB', 'COL6A1', 'COL6A2', 'COL5A1', 'FBLN1', 'IGFBP6',
    ],
    'T cells': [
        'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B', 'TRAC', 'IL7R',
        'PTPRC', 'NKG7', 'GZMA', 'GZMB', 'GZMK', 'CCL5', 'CD2',
        'LCK', 'ITK', 'CD7', 'TRBC1', 'TRBC2',
    ],
    'Smooth muscle cells': [
        'ACTA2', 'MYH11', 'TAGLN', 'CNN1', 'TPM2', 'DES', 'CALD1',
        'LMOD1', 'TPM1', 'RGS5', 'MCAM', 'SORBS2', 'NEXN', 'MYOCD',
        'PRKG1', 'SYNPO2',
    ],
    'Neutrophils': [
        'S100A8', 'S100A9', 'CSF3R', 'FCGR3B', 'CXCR2', 'LYZ',
        'CXCL8', 'IL1B', 'MMP9', 'SELL', 'FPR1', 'NAMPT',
    ],
    'B cells': [
        'CD79A', 'CD79B', 'MS4A1', 'CD19', 'PAX5', 'BANK1', 'IGHM',
        'IGHD', 'HLA-DRA', 'CD22', 'BLK', 'FCER2', 'CR2',
        'REL', 'TCL1A',
    ],
    'Myeloids': [
        'CD14', 'CD68', 'CSF1R', 'AIF1', 'TYROBP', 'LYZ', 'HLA-DRA',
        'ITGAM', 'C1QA', 'C1QB', 'C1QC', 'CD163', 'MSR1', 'MARCO',
        'FCGR3A', 'FCER1G', 'CTSB', 'CTSS',
    ],
    'Erythroid cells': [
        'HBA1', 'HBA2', 'HBB', 'ALAS2', 'GYPA', 'GATA1', 'BPGM',
        'SNCA', 'EPB41', 'SLC4A1', 'ANK1', 'KLF1',
    ],
    'Endothelial cells': [
        'PECAM1', 'VWF', 'CDH5', 'FLT1', 'KDR', 'ENG', 'RAMP2',
        'CAV1', 'LDB2', 'PTPRB', 'CLDN5', 'ERG', 'EMCN', 'CAVIN2',
        'SELE', 'TM4SF1', 'TCF4',
    ],
    'Plasma cells': [
        'JCHAIN', 'MZB1', 'SDC1', 'XBP1', 'IGHG1', 'IGHA1', 'IGHM',
        'PRDM1', 'IRF4', 'SLAMF7',
    ],
    'Mast cells': [
        'KIT', 'CPA3', 'TPSAB1', 'TPSB2', 'MS4A2', 'HPGDS',
        'HDC', 'GATA2', 'FCER1A', 'CLU',
    ],
    'Schwann cells': [
        'S100B', 'PLP1', 'MPZ', 'SOX10', 'CDH19', 'MBP', 'CADM2',
        'SCN7A', 'NDRG1', 'PRX', 'EGR2', 'NCAM1',
    ],
}

# 细胞类型缩写
CELL_TYPE_ABBREV = {
    'Fibroblasts': 'FIB',
    'T cells': 'T',
    'Smooth muscle cells': 'SMC',
    'Neutrophils': 'NEU',
    'B cells': 'B',
    'Myeloids': 'MYE',
    'Erythroid cells': 'ERY',
    'Endothelial cells': 'ENDO',
    'Plasma cells': 'PLA',
    'Mast cells': 'MAST',
    'Schwann cells': 'SC',
}

# 固定细胞类型顺序
CELL_TYPE_ORDER = [
    'Endothelial cells',
    'Fibroblasts',
    'Smooth muscle cells',
    'T cells',
    'Myeloids',
    'Neutrophils',
    'B cells',
    'Plasma cells',
    'Mast cells',
    'Erythroid cells',
    'Schwann cells',
]

# 为每种细胞类型分配颜色（参考例图风格）
CELL_TYPE_COLORS = {
    'Endothelial cells': '#5B9BD5',
    'Fibroblasts': '#ED7D31',
    'Smooth muscle cells': '#A5A5A5',
    'T cells': '#C44E52',
    'Myeloids': '#70AD47',
    'Neutrophils': '#9B59B6',
    'B cells': '#E74C3C',
    'Plasma cells': '#E74C3C',#'#00B0F0',
    'Mast cells': '#FF6699',
    'Erythroid cells': '#2ECC71',
    'Schwann cells': '#8B4513',
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
                                pts=True)
        print("Done.")

    return adata


# ============================================================
# 3. 提取 top10 基因并匹配经典 marker
# ============================================================
def get_top_genes_and_markers(adata, n_top=10):
    result = sc.get.rank_genes_groups_df(adata, group=None,pval_cutoff=0.05)
    top_genes_dict = {}
    marker_info_dict = {}

    for ct in CELL_TYPE_ORDER:
        if ct not in result['group'].unique(): 
            print(f"Warning: {ct} not in rank_genes_groups, skipping.")
            continue
        result_ct = result[result['group'].isin([ct])].sort_values('logfoldchanges', ascending=False)
        names = result_ct['names'].tolist()

        filtered = [g for g in names if not g.startswith('ENSOCUG') and not g.startswith('RA_') and not g.startswith('RP')]
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
    max_cells = 5000
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
        (0.0,  '#08306B'),   # 深蓝
        (0.15, '#2171B5'),   # 蓝
        (0.3,  '#6BAED6'),   # 浅蓝
        (0.45, '#C6DBEF'),   # 很浅蓝
        (0.5,  '#F7F7F7'),   # 白
        (0.55, '#FDDBC7'),   # 很浅红
        (0.7,  '#F4A582'),   # 浅红
        (0.85, '#D6604D'),   # 红
        (1.0,  '#8B0000'),   # 深红
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

        # 三角形标记 (◁)
        tri = plt.Polygon(
            [[0.2, mid_y - 0.6], [0.2, mid_y + 0.6], [0.7, mid_y]],
            color=color, zorder=5, clip_on=False
        )
        ax_annot.add_patch(tri)

        # 标注文本
        if markers_found:
            # 有经典 marker: 显示 "缩写: marker1, marker2, ..."
            marker_str = ', '.join(markers_found)
            label = f"{abbrev}: {marker_str}"
            ax_annot.text(
                0.9, mid_y, label,
                fontsize=7.5, fontweight='bold', color=color,
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

    # # 图例（放在 colorbar 右侧）
    # legend_elements = []
    # for ct in CELL_TYPE_ORDER:
    #     if ct in ct_gene_ranges:
    #         legend_elements.append(
    #             mpatches.Patch(
    #                 facecolor=CELL_TYPE_COLORS.get(ct, '#999'),
    #                 edgecolor='none', label=ct
    #             )
    #         )
    # fig.legend(
    #     handles=legend_elements,
    #     loc='lower right', bbox_to_anchor=(0.95, 0.001),
    #     ncol=3, fontsize=5.5, frameon=False,
    #     handlelength=1.0, handleheight=0.8,
    #     columnspacing=0.5, labelspacing=0.3
    # )

    # plt.subplots_adjust(left=0.04, right=0.97, top=0.99, bottom=0.045)

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

    data_dir = r'/mnt/NexData_file/home/zwl/Project/SC2601160389/售后1/'
    h5ad_path = os.path.join(data_dir, 'result.h5ad')
    output_prefix = os.path.join('/mnt/NexData_file/home/zwl/Project/SC2601160389/售后2/', 'heatmap/heatmap_top10_markers')

    np.random.seed(42)

    adata = load_and_compute(h5ad_path)
    top_genes_dict, marker_info_dict = get_top_genes_and_markers(adata, n_top=10)
    plot_heatmap(adata, top_genes_dict, marker_info_dict, output_prefix)

    print("\nAll done!")
