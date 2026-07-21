#!/usr/bin/env python3
"""
整合分析流水线：细胞组成分析 + ROE分析 + DEG分析(scanpy rank_genes_groups) + 富集分析
完整版（优化版）：
1. 箱线图Y轴从0开始，添加ns标注
2. 不使用FDR校正，直接使用p_value
3. 火山图优化：去掉基因标签，图例在右上角，小字体无边框，0在中间
4. 富集分析：KEGG + GO_BP + GO_CC + GO_MF，气泡图X轴大幅放宽、Y轴间距拉开、尺寸动态适配
5. 柱状图Y轴间距拉开，上下调图形高度一致，长名称自动扩宽
"""

import os
import sys
import warnings
import logging
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import font_manager
from matplotlib.lines import Line2D
import seaborn as sns
from scipy import stats
from scipy.stats import mannwhitneyu, chi2_contingency
from itertools import combinations
import statsmodels.stats.multitest as smm

# 减少Matplotlib的字体日志输出
logging.getLogger('fontTools').setLevel(logging.WARNING)
logging.getLogger('matplotlib.font_manager').setLevel(logging.WARNING)

warnings.filterwarnings('ignore')

# ============================================================
# 字体设置
# ============================================================

ARIAL_PATH = '/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf'
if os.path.exists(ARIAL_PATH):
    font_manager.fontManager.addfont(ARIAL_PATH)
else:
    alt_paths = [
        '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf',
        '/System/Library/Fonts/Arial.ttf',
        'C:/Windows/Fonts/arial.ttf'
    ]
    for alt in alt_paths:
        if os.path.exists(alt):
            font_manager.fontManager.addfont(alt)
            break
    else:
        print("Warning: Arial font not found, using default sans-serif")

# ============================================================
# 全局字体大小配置
# ============================================================

FONT_TITLE = 65
FONT_AXIS_LABEL = 60
FONT_TICK = 50
FONT_STAT = 48
FONT_LEGEND = 48
FONT_ANNOTATION = 45
FONT_SMALL = 36  # 更小的字体用于图例
FONT_COLORBAR = 48
FONT_SUPTITLE = 65
FONT_HEATMAP_CELL = 30
FONT_HEATMAP_AXIS = 45

# ============================================================
# matplotlib全局设置
# ============================================================

plt.rcParams.update({
    'font.family': 'arial',
    'font.weight': 'heavy',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'font.size': FONT_TICK,
    'axes.titlesize': FONT_TITLE,
    'axes.labelsize': FONT_AXIS_LABEL,
    'xtick.labelsize': FONT_TICK,
    'ytick.labelsize': FONT_TICK,
    'legend.fontsize': FONT_LEGEND,
    'figure.titlesize': FONT_TITLE,
    'lines.linewidth': 3.0,
    'axes.linewidth': 3.0,
    'xtick.major.width': 2.5,
    'ytick.major.width': 2.5,
    'xtick.major.size': 10,
    'ytick.major.size': 10,
    'xtick.major.pad': 12,
    'ytick.major.pad': 12,
    'axes.labelpad': 20,
    'axes.titlepad': 30,
})

# ============================================================
# 配色方案
# ============================================================

COLORS_GROUPS = {
    'PRE_pCR': '#D62728',
    'PRE_nonpCR': '#FF7F0E',
    'POST_pCR': '#1F77B4',
    'POST_nonpCR': '#2CA02C',
}

COLOR_UP = '#D62728'
COLOR_DOWN = '#1F77B4'
COLOR_NS = '#BDBDBD'
COLOR_PALETTE = ['#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
                 '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']

# ============================================================
# 工具函数
# ============================================================

def safe_name(name):
    return str(name).replace(' ', '_').replace('/', '_').replace('\\', '_').replace('(', '').replace(')', '').replace('+', 'plus')

def calculate_roe(observed_table):
    row_sums = observed_table.sum(axis=1).values.reshape(-1, 1)
    col_sums = observed_table.sum(axis=0).values.reshape(1, -1)
    total = observed_table.sum().sum()
    expected = (row_sums @ col_sums) / total
    roe = observed_table / expected
    return roe

def add_significance_annotation_v2(ax, x1, x2, y, p_value, height_offset=0.03):
    """
    添加显著性标注 - 改进版
    如果p_value >= 0.05，显示 'ns'
    """
    if p_value < 0.001:
        sig_text = '****'
    elif p_value < 0.01:
        sig_text = '***'
    elif p_value < 0.05:
        sig_text = '**'
    else:
        sig_text = 'ns'  # 不显著时显示 ns

    y_max = y + height_offset
    ax.plot([x1, x1, x2, x2], [y, y_max, y_max, y], lw=3.0, color='black')
    ax.text((x1 + x2) / 2, y_max + height_offset * 0.3, sig_text, ha='center', va='bottom',
            fontsize=FONT_ANNOTATION, fontweight='heavy', color='black')
    return y_max + height_offset * 1.5

# ============================================================
# 第一部分：细胞比例分析与箱线图（Y轴从0开始，不使用FDR，添加ns）
# ============================================================

def calculate_cell_proportions(adata, group_col='group_sub', celltype_col='singler_best', sample_col='sample'):
    logger = logging.getLogger('Pipeline')
    obs = adata.obs[[sample_col, group_col, celltype_col]].copy()
    obs = obs.dropna(subset=[celltype_col, group_col])
    for col in [sample_col, group_col, celltype_col]:
        obs[col] = obs[col].astype(str).astype('category')
    counts = obs.groupby([sample_col, group_col, celltype_col], observed=False).size().reset_index(name='count')
    group_totals = counts.groupby([sample_col, group_col], observed=False)['count'].transform('sum')
    counts['proportion'] = counts['count'] / group_totals
    counts['percentage'] = counts['proportion'] * 100
    logger.info(f"细胞比例计算完成，共 {counts[celltype_col].nunique()} 种细胞类型")
    return counts, obs

def statistical_test_proportions_no_fdr(prop_df, group_col='group_sub', celltype_col='singler_best', p_threshold=0.05):
    """
    统计检验 - 不使用FDR校正，直接使用p_value
    """
    logger = logging.getLogger('Pipeline')
    groups = sorted(prop_df[group_col].unique())
    celltypes = sorted(prop_df[celltype_col].unique())
    results = []
    for ct in celltypes:
        sub = prop_df[prop_df[celltype_col] == ct]
        for i, g1 in enumerate(groups):
            for g2 in groups[i+1:]:
                vals1 = sub[sub[group_col] == g1]['proportion'].values
                vals2 = sub[sub[group_col] == g2]['proportion'].values
                if len(vals1) >= 3 and len(vals2) >= 3:
                    stat, p = mannwhitneyu(vals1, vals2, alternative='two-sided')
                    results.append({
                        'celltype': ct, 'group1': g1, 'group2': g2,
                        'p_value': p, 'mean_diff': vals1.mean() - vals2.mean(),
                        'mean1': vals1.mean(), 'mean2': vals2.mean()
                    })
    result_df = pd.DataFrame(results)
    sig_celltypes = set()
    if len(result_df) > 0:
        sig_rows = result_df[result_df['p_value'] < p_threshold]
        sig_celltypes = set(sig_rows['celltype'].unique())
    logger.info(f"统计检验完成，{len(sig_celltypes)} 种细胞类型有显著差异 (p < {p_threshold})")
    return result_df, sig_celltypes

def plot_celltype_boxplot_cns_fixed_v2(prop_df, ct, group_col='group_sub', output_dir='.',
                                        stat_results=None, show_significance=True):
    """
    CNS级别箱线图 - 改进版
    改进：
    1. 添加 'ns' 标注
    2. 显著性标注样式调整
    """
    sub = prop_df[prop_df['singler_best'] == ct].copy()
    if sub.shape[0] == 0:
        return None, None

    groups = sorted(sub[group_col].unique())
    n_groups = len(groups)
    sub['percentage'] = sub['proportion'] * 100

    # ============================================================
    # 箱线图尺寸
    # ============================================================
    fig_width = max(14, n_groups * 6)
    fig_height = 14

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # 准备数据
    plot_data = []
    positions = []
    colors = []
    for i, g in enumerate(groups):
        vals = sub[sub[group_col] == g]['percentage'].values
        vals_clean = vals[~np.isnan(vals) & ~np.isinf(vals)]
        plot_data.append(vals_clean)
        positions.append(i)
        colors.append(COLORS_GROUPS.get(g, COLOR_PALETTE[i % len(COLOR_PALETTE)]))

    # ============================================================
    # 绘制箱线图
    # ============================================================
    bp = ax.boxplot(plot_data, positions=positions, widths=0.6,
                     patch_artist=True, showfliers=False,
                     boxprops=dict(linewidth=3.5, edgecolor='black'),
                     whiskerprops=dict(linewidth=3.0, color='black'),
                     capprops=dict(linewidth=3.0, color='black'),
                     medianprops=dict(linewidth=4.0, color='darkred'))

    # 填充箱子颜色
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # 添加散点
    for i, vals in enumerate(plot_data):
        if len(vals) > 0:
            jitter = np.random.normal(0, 0.08, size=len(vals))
            ax.scatter(np.full(len(vals), positions[i]) + jitter, vals,
                      color='black', s=200, alpha=0.6, zorder=3, linewidths=2.0,
                      edgecolors='white')

    # ============================================================
    # 显著性标注 - 使用p_value，包含ns
    # ============================================================
    if show_significance and stat_results is not None:
        # 获取该细胞类型的所有比较结果
        ct_results = stat_results[stat_results['celltype'] == ct]
        if len(ct_results) > 0:
            # 计算合适的Y轴位置
            all_vals = []
            for vals in plot_data:
                if len(vals) > 0:
                    all_vals.extend(vals)
            
            if len(all_vals) > 0:
                y_max_data = max(all_vals)
            else:
                y_max_data = 1.0
            
            # 显著性标注放在数据上方
            y_pos = y_max_data * 1.15
            height_offset = y_max_data * 0.08
            
            for _, row in ct_results.iterrows():
                try:
                    x1 = groups.index(row['group1'])
                    x2 = groups.index(row['group2'])
                    y_pos = add_significance_annotation_v2(ax, x1, x2, y_pos, row['p_value'], height_offset)
                except ValueError:
                    continue

    # ============================================================
    # 美化 - Y轴从0开始
    # ============================================================
    ax.set_xticks(positions)
    ax.set_xticklabels(groups, rotation=30, ha='right', fontsize=FONT_TICK, fontweight='heavy')
    ax.set_ylabel('Cell Proportion', fontsize=FONT_AXIS_LABEL, fontweight='heavy')
    ax.set_title(f'{ct}', fontsize=FONT_TITLE, fontweight='heavy', pad=35)

    # Y轴从0开始
    all_valid_vals = []
    for vals in plot_data:
        if len(vals) > 0:
            vals_clean = vals[~np.isnan(vals) & ~np.isinf(vals)]
            if len(vals_clean) > 0:
                all_valid_vals.extend(vals_clean)
    
    if len(all_valid_vals) > 0:
        y_max = max(all_valid_vals)
        if np.isnan(y_max) or np.isinf(y_max) or y_max <= 0:
            y_max = 10.0
        ax.set_ylim(0, y_max * 1.4)
    else:
        ax.set_ylim(0, 20)
    
    ax.set_xlim(-0.8, n_groups - 0.2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(4.0)
    ax.spines['bottom'].set_linewidth(4.0)

    ax.yaxis.grid(True, linestyle='--', alpha=0.3, linewidth=2.5)
    ax.set_axisbelow(True)

    ax.tick_params(axis='both', which='major', width=4.0, length=15, pad=18)
    ax.tick_params(axis='y', which='major', pad=18)

    plt.tight_layout(pad=4.0)

    safe_ct = safe_name(ct)
    fig_path_png = os.path.join(output_dir, f'boxplot_{safe_ct}_v2.png')
    fig_path_pdf = os.path.join(output_dir, f'boxplot_{safe_ct}_v2.pdf')
    fig.savefig(fig_path_png, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    fig.savefig(fig_path_pdf, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)

    return fig_path_png, ct

def plot_all_boxplots_fixed_v2(prop_df, output_dir, group_col='group_sub',
                                stat_results=None, sig_celltypes=None):
    """
    绘制所有细胞类型的箱线图（改进版）
    """
    logger = logging.getLogger('Pipeline')
    all_dir = os.path.join(output_dir, 'boxplots_all')
    sig_dir = os.path.join(output_dir, 'sig')
    os.makedirs(all_dir, exist_ok=True)
    os.makedirs(sig_dir, exist_ok=True)

    celltypes = sorted(prop_df['singler_best'].unique())
    all_paths = []
    sig_paths = []

    for ct in celltypes:
        is_sig = sig_celltypes is not None and ct in sig_celltypes
        save_dir = sig_dir if is_sig else all_dir
        fig_path, _ = plot_celltype_boxplot_cns_fixed_v2(
            prop_df, ct, group_col=group_col, output_dir=save_dir,
            stat_results=stat_results, show_significance=is_sig
        )
        if fig_path:
            all_paths.append(fig_path)
            if is_sig:
                sig_paths.append(fig_path)

    logger.info(f"箱线图绘制完成: 总计 {len(all_paths)} 张, 显著 {len(sig_paths)} 张 (sig文件夹)")

    if sig_celltypes and len(sig_celltypes) > 0:
        plot_summary_boxplot_grid_fixed_v2(prop_df, sig_celltypes, output_dir, group_col, stat_results)

    return all_paths, sig_paths

def plot_summary_boxplot_grid_fixed_v2(prop_df, sig_celltypes, output_dir, group_col, stat_results):
    """
    汇总箱线图网格（改进版）
    """
    n_cells = len(sig_celltypes)
    if n_cells == 0:
        return

    n_cols = min(2, n_cells)
    n_rows = (n_cells + n_cols - 1) // n_cols

    subplot_width = 14
    subplot_height = 12

    fig_width = n_cols * subplot_width + 3
    fig_height = n_rows * subplot_height + 3

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
    if n_cells == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    groups = sorted(prop_df[group_col].unique())

    for idx, ct in enumerate(sorted(sig_celltypes)):
        ax = axes[idx]
        sub = prop_df[prop_df['singler_best'] == ct].copy()
        sub['percentage'] = sub['proportion'] * 100

        plot_data = []
        positions = []
        colors_list = []
        for i, g in enumerate(groups):
            vals = sub[sub[group_col] == g]['percentage'].values
            vals_clean = vals[~np.isnan(vals) & ~np.isinf(vals)]
            plot_data.append(vals_clean)
            positions.append(i)
            colors_list.append(COLORS_GROUPS.get(g, COLOR_PALETTE[i % len(COLOR_PALETTE)]))

        # 绘制箱线图
        bp = ax.boxplot(plot_data, positions=positions, widths=0.6,
                         patch_artist=True, showfliers=False,
                         boxprops=dict(linewidth=3.0, edgecolor='black'),
                         whiskerprops=dict(linewidth=3.0, color='black'),
                         capprops=dict(linewidth=3.0, color='black'),
                         medianprops=dict(linewidth=4.0, color='darkred'))

        for patch, color in zip(bp['boxes'], colors_list):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        # 添加散点
        for i, vals in enumerate(plot_data):
            if len(vals) > 0:
                jitter = np.random.normal(0, 0.06, size=len(vals))
                ax.scatter(np.full(len(vals), positions[i]) + jitter, vals,
                          color='black', s=120, alpha=0.6, zorder=3, linewidths=1.5,
                          edgecolors='white')

        # 显著性标注
        ct_results = stat_results[stat_results['celltype'] == ct]
        if len(ct_results) > 0:
            all_valid_vals = []
            for vals in plot_data:
                if len(vals) > 0:
                    all_valid_vals.extend(vals)
            
            if len(all_valid_vals) > 0:
                y_max = max(all_valid_vals)
            else:
                y_max = 1.0
            
            height_offset = y_max * 0.10
            y_pos = y_max * 1.08
            for _, row in ct_results.iterrows():
                try:
                    x1 = groups.index(row['group1'])
                    x2 = groups.index(row['group2'])
                    y_pos = add_significance_annotation_v2(ax, x1, x2, y_pos, row['p_value'], height_offset)
                except ValueError:
                    continue

        ax.set_xticks(positions)
        ax.set_xticklabels(groups, rotation=30, ha='right', fontsize=FONT_SMALL + 4, fontweight='heavy')
        ax.set_title(f'{ct}', fontsize=FONT_SMALL + 8, fontweight='heavy', pad=20)
        ax.set_ylabel('Cell Proportion', fontsize=FONT_SMALL + 2, fontweight='heavy')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Y轴从0开始
        all_valid_vals = []
        for vals in plot_data:
            if len(vals) > 0:
                vals_clean = vals[~np.isnan(vals) & ~np.isinf(vals)]
                if len(vals_clean) > 0:
                    all_valid_vals.extend(vals_clean)
        
        if len(all_valid_vals) > 0:
            y_max = max(all_valid_vals)
            if np.isnan(y_max) or np.isinf(y_max) or y_max <= 0:
                y_max = 10.0
            ax.set_ylim(0, y_max * 1.35)
        else:
            ax.set_ylim(0, 10)
        
        ax.set_xlim(-0.8, len(groups) - 0.2)

    for idx in range(n_cells, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle('Significant Cell Type Proportions Across Groups',
                 fontsize=FONT_SUPTITLE, fontweight='heavy', y=1.02)
    plt.tight_layout(pad=5.0)

    summary_path = os.path.join(output_dir, 'summary_significant_boxplots.png')
    summary_pdf = os.path.join(output_dir, 'summary_significant_boxplots.pdf')
    fig.savefig(summary_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    fig.savefig(summary_pdf, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)

    logging.getLogger('Pipeline').info(f"汇总图已保存: {summary_path}")

# ============================================================
# 第二部分：ROE分析与热图
# ============================================================

def perform_roe_analysis(obs_df, group_col='group_sub', celltype_col='singler_best'):
    logger = logging.getLogger('Pipeline')
    crosstab = pd.crosstab(obs_df[celltype_col], obs_df[group_col])
    roe = calculate_roe(crosstab)
    groups = sorted(obs_df[group_col].unique())
    comparisons = list(combinations(groups, 2))

    roe_results = []
    log2fc_data = {}

    for g1, g2 in comparisons:
        comp_name = f"{g1}\nvs\n{g2}"
        comp_name_simple = f"{g1}_vs_{g2}"
        log2fc = np.log2(roe[g1] / roe[g2])
        log2fc_data[comp_name] = log2fc
        for ct in roe.index:
            roe_results.append({
                'celltype': ct, 'comparison': comp_name_simple,
                'ROE_group1': roe.loc[ct, g1], 'ROE_group2': roe.loc[ct, g2],
                'Log2FC': log2fc[ct] if not np.isinf(log2fc[ct]) else np.nan,
                'Enriched_in': g1 if log2fc[ct] > 0 else (g2 if log2fc[ct] < 0 else 'None')
            })

    roe_df = pd.DataFrame(roe_results)
    log2fc_df = pd.DataFrame(log2fc_data)
    log2fc_df.index.name = 'celltype'
    log2fc_df = log2fc_df.reset_index()
    logger.info(f"ROE分析完成，共 {len(comparisons)} 组比较")
    return roe_df, log2fc_df, crosstab

def plot_roe_heatmap_cns(log2fc_df, output_dir, title="ROE Analysis",
                          max_abs=None, top_n=None, sig_celltypes=None):
    logger = logging.getLogger('Pipeline')
    plot_df = log2fc_df.set_index('celltype')

    if sig_celltypes is not None and len(sig_celltypes) > 0:
        plot_df = plot_df[plot_df.index.isin(sig_celltypes)]
    elif top_n is not None and top_n < len(plot_df):
        mean_abs = plot_df.abs().mean(axis=1)
        top_cells = mean_abs.nlargest(top_n).index
        plot_df = plot_df.loc[top_cells]

    if len(plot_df) == 0:
        logger.warning("没有数据用于热图")
        return None

    plot_df['mean_abs'] = plot_df.abs().mean(axis=1)
    plot_df = plot_df.sort_values('mean_abs', ascending=False)
    plot_df = plot_df.drop(columns='mean_abs')

    display_cols = [c.replace('\nvs\n', ' vs ') for c in plot_df.columns]

    if max_abs is None:
        max_abs = plot_df.abs().max().max()
        max_abs = max(max_abs, 0.5)

    n_rows = len(plot_df)
    n_cols = len(plot_df.columns)
    fig_width = max(12, n_cols * 6)
    fig_height = max(10, n_rows * 1.5)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    im = ax.imshow(plot_df.values, cmap=cmap, aspect='auto',
                    vmin=-max_abs, vmax=max_abs)

    for i in range(n_rows):
        for j in range(n_cols):
            val = plot_df.values[i, j]
            if not np.isnan(val):
                text_color = 'white' if abs(val) > max_abs * 0.5 else 'black'
                ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                       fontsize=FONT_HEATMAP_CELL, color=text_color, fontweight='heavy')

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(display_cols, rotation=45, ha='right',
                       fontsize=FONT_HEATMAP_AXIS, fontweight='heavy')
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(plot_df.index, fontsize=FONT_HEATMAP_AXIS, fontweight='heavy')

    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Log2(Fold Change)', fontsize=FONT_COLORBAR, fontweight='heavy', labelpad=20)
    cbar.ax.tick_params(labelsize=FONT_SMALL, width=3.0, length=8)
    cbar.outline.set_linewidth(3.0)

    ax.set_title(title, fontsize=FONT_TITLE, fontweight='heavy', pad=30)
    ax.set_xlabel('Comparison', fontsize=FONT_AXIS_LABEL, fontweight='heavy', labelpad=20)
    ax.set_ylabel('Cell Type', fontsize=FONT_AXIS_LABEL, fontweight='heavy', labelpad=20)

    for spine in ax.spines.values():
        spine.set_linewidth(3.0)

    plt.tight_layout(pad=4.0)

    heatmap_path = os.path.join(output_dir, 'roe_heatmap.png')
    heatmap_path_pdf = os.path.join(output_dir, 'roe_heatmap.pdf')
    fig.savefig(heatmap_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    fig.savefig(heatmap_path_pdf, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)

    logger.info(f"ROE热图已保存: {heatmap_path}")
    return heatmap_path

# ============================================================
# 第三部分：DEG分析 - 使用scanpy rank_genes_groups
# 优化火山图：去掉基因标签，图例在右上角，小字体无边框，0在中间
# ============================================================

def perform_dea_scanpy(adata, celltype, group1, group2,
                       celltype_col='singler_best', group_col='group_sub',
                       min_cells=10, method='wilcoxon'):
    """
    使用scanpy的rank_genes_groups进行差异表达分析
    """
    logger = logging.getLogger('Pipeline')

    # 筛选特定细胞类型
    celltype_data = adata[adata.obs[celltype_col] == celltype].copy()
    group1_data = celltype_data[celltype_data.obs[group_col] == group1].copy()
    group2_data = celltype_data[celltype_data.obs[group_col] == group2].copy()

    if group1_data.n_obs < min_cells or group2_data.n_obs < min_cells:
        logger.warning(f"{celltype} 细胞数不足 (组1:{group1_data.n_obs}, 组2:{group2_data.n_obs})，跳过DEG分析")
        return None

    # 合并数据
    dea_data = sc.concat([group1_data, group2_data])

    # 修复raw丢失的问题
    if dea_data.raw is None:
        dea_data.raw = dea_data.copy()

    # 使用原始计数
    dea_data.X = dea_data.raw.X
    dea_data.obs['comparison_group'] = pd.Categorical(
        [group1] * group1_data.n_obs + [group2] * group2_data.n_obs
    )

    # 标准化
    sc.pp.normalize_total(dea_data, target_sum=1e4)
    sc.pp.log1p(dea_data)

    # rank_genes_groups
    sc.tl.rank_genes_groups(dea_data, groupby='comparison_group', method=method, use_raw=False)

    # 获取结果
    de_results = sc.get.rank_genes_groups_df(dea_data, group=group1)
    de_results['celltype'] = celltype
    de_results['comparison'] = f"{group1}_vs_{group2}"
    de_results['-log10_padj'] = -np.log10(de_results['pvals_adj'].clip(lower=1e-300))

    logger.info(f"DEG完成: {celltype} ({group1} vs {group2})，共 {len(de_results)} 个基因")
    return de_results

def plot_volcano_cns_optimized(de_results, output_dir, comparison_name, celltype_name,
                               logfc_threshold=0.25, pval_threshold=0.05):
    """
    CNS级别火山图 - 优化版
    优化：
    1. 去掉基因标签
    2. 图例放在右上角，小字体无边框
    3. 整体放大
    4. X轴0放在中间
    """
    if de_results is None or len(de_results) == 0:
        return None

    df = de_results.copy()

    df['category'] = 'NS'
    up_mask = (df['pvals_adj'] < pval_threshold) & (df['logfoldchanges'] > logfc_threshold)
    down_mask = (df['pvals_adj'] < pval_threshold) & (df['logfoldchanges'] < -logfc_threshold)
    df.loc[up_mask, 'category'] = 'Up'
    df.loc[down_mask, 'category'] = 'Down'

    n_up = up_mask.sum()
    n_down = down_mask.sum()

    # ============================================================
    # 整体放大：增大图表尺寸
    # ============================================================
    fig, ax = plt.subplots(figsize=(22, 18))

    colors = {'NS': COLOR_NS, 'Up': COLOR_UP, 'Down': COLOR_DOWN}
    sizes = {'NS': 30, 'Up': 80, 'Down': 80}
    alphas = {'NS': 0.3, 'Up': 0.8, 'Down': 0.8}
    zorders = {'NS': 1, 'Up': 3, 'Down': 2}

    for cat in ['NS', 'Down', 'Up']:
        subset = df[df['category'] == cat]
        if len(subset) > 0:
            ax.scatter(subset['logfoldchanges'], subset['-log10_padj'],
                      c=colors[cat], s=sizes[cat], alpha=alphas[cat],
                      zorder=zorders[cat], label=f'{cat} ({len(subset)})',
                      edgecolors='none')

    # ============================================================
    # 阈值线
    # ============================================================
    ax.axhline(-np.log10(pval_threshold), linestyle='--', color='grey',
               linewidth=3.0, alpha=0.6)
    ax.axvline(logfc_threshold, linestyle='--', color='grey',
               linewidth=3.0, alpha=0.6)
    ax.axvline(-logfc_threshold, linestyle='--', color='grey',
               linewidth=3.0, alpha=0.6)
    
    # ============================================================
    # 关键：在X轴0位置画一条垂直基准线
    # ============================================================
    ax.axvline(0, linestyle='-', color='black', linewidth=2.0, alpha=0.5)

    ax.set_xlabel('Log2 Fold Change', fontsize=FONT_AXIS_LABEL, fontweight='heavy', labelpad=25)
    ax.set_ylabel('-Log10 Adjusted P-value', fontsize=FONT_AXIS_LABEL, fontweight='heavy', labelpad=25)
    ax.set_title(f'{celltype_name}\n{comparison_name}',
                 fontsize=FONT_TITLE, fontweight='heavy', pad=35)

    # ============================================================
    # 图例：右上角，小字体，无边框
    # ============================================================
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_UP,
               markersize=15, label=f'Up ({n_up})'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_DOWN,
               markersize=15, label=f'Down ({n_down})'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_NS,
               markersize=12, label=f'NS ({len(df) - n_up - n_down})')
    ]
    
    ax.legend(handles=legend_elements, loc='upper right', 
             frameon=False, fontsize=FONT_SMALL,
             borderaxespad=1.0, handletextpad=0.5)

    # ============================================================
    # 设置X轴范围，确保0在中间
    # ============================================================
    # 计算最大绝对值
    max_fc = max(abs(df['logfoldchanges'].max()), abs(df['logfoldchanges'].min()))
    # 加上一点边距
    x_limit = max_fc * 1.1
    ax.set_xlim(-x_limit, x_limit)
    
    # 确保0在X轴中间
    ax.set_xlim(-x_limit, x_limit)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(4.0)
    ax.spines['bottom'].set_linewidth(4.0)
    ax.grid(True, linestyle='--', alpha=0.2, linewidth=2.5)

    ax.tick_params(axis='both', which='major', width=4.0, length=15, pad=18)

    plt.tight_layout(pad=4.0)

    safe_ct = safe_name(celltype_name)
    safe_comp = safe_name(comparison_name)
    volcano_path = os.path.join(output_dir, f'volcano_{safe_comp}_{safe_ct}.png')
    volcano_path_pdf = os.path.join(output_dir, f'volcano_{safe_comp}_{safe_ct}.pdf')

    fig.savefig(volcano_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    fig.savefig(volcano_path_pdf, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)

    return volcano_path

# ============================================================
# 第四部分：富集分析 - KEGG + GO_BP + GO_CC + GO_MF
# 气泡图优化：X轴更宽，气泡大小动态调整，Y轴间距拉开
# 柱状图优化：Y轴间距拉开，上下调图形高度一致
# ============================================================

def perform_enrichment_analysis(deg_df, output_dir, organism='human'):
    """
    对DEG结果进行富集分析，KEGG + GO_BP + GO_CC + GO_MF
    """
    logger = logging.getLogger('Pipeline')
    
    try:
        import gseapy as gp
    except ImportError:
        logger.warning("gseapy未安装，跳过富集分析")
        logger.warning("请安装: pip install gseapy")
        return None
    
    if deg_df is None or len(deg_df) == 0:
        logger.warning("没有DEG数据，跳过富集分析")
        return None
    
    # 提取上调基因和下调基因
    up_genes = deg_df[
        (deg_df['pvals_adj'] < 0.05) & 
        (deg_df['logfoldchanges'] > 0.25)
    ]['names'].tolist()
    
    down_genes = deg_df[
        (deg_df['pvals_adj'] < 0.05) & 
        (deg_df['logfoldchanges'] < -0.25)
    ]['names'].tolist()
    
    logger.info(f"富集分析: 上调基因 {len(up_genes)} 个，下调基因 {len(down_genes)} 个")
    
    # 设置基因集
    if organism == 'human':
        kegg_set = 'KEGG_2021_Human'
        go_bp_set = 'GO_Biological_Process_2021'
        go_cc_set = 'GO_Cellular_Component_2021'
        go_mf_set = 'GO_Molecular_Function_2021'
    else:
        kegg_set = 'KEGG_2021_Mouse'
        go_bp_set = 'GO_Biological_Process_2021'
        go_cc_set = 'GO_Cellular_Component_2021'
        go_mf_set = 'GO_Molecular_Function_2021'
    
    # 存储所有富集结果
    all_enrichment_results = []
    
    # 对上调基因进行富集
    if len(up_genes) >= 10:
        logger.info("对上调基因进行富集分析")
        
        # KEGG富集
        try:
            kegg_enr = gp.enrichr(gene_list=up_genes, gene_sets=kegg_set, 
                                 organism=organism, outdir=None)
            if kegg_enr.results is not None and len(kegg_enr.results) > 0:
                kegg_enr.results['gene_set'] = 'KEGG'
                kegg_enr.results['direction'] = 'up'
                all_enrichment_results.append(kegg_enr.results)
        except Exception as e:
            logger.warning(f"KEGG富集分析失败: {e}")
        
        # GO BP富集
        try:
            go_bp_enr = gp.enrichr(gene_list=up_genes, gene_sets=go_bp_set, 
                                  organism=organism, outdir=None)
            if go_bp_enr.results is not None and len(go_bp_enr.results) > 0:
                go_bp_enr.results['gene_set'] = 'GO_BP'
                go_bp_enr.results['direction'] = 'up'
                all_enrichment_results.append(go_bp_enr.results)
        except Exception as e:
            logger.warning(f"GO BP富集分析失败: {e}")
        
        # GO CC富集
        try:
            go_cc_enr = gp.enrichr(gene_list=up_genes, gene_sets=go_cc_set, 
                                  organism=organism, outdir=None)
            if go_cc_enr.results is not None and len(go_cc_enr.results) > 0:
                go_cc_enr.results['gene_set'] = 'GO_CC'
                go_cc_enr.results['direction'] = 'up'
                all_enrichment_results.append(go_cc_enr.results)
        except Exception as e:
            logger.warning(f"GO CC富集分析失败: {e}")
        
        # GO MF富集
        try:
            go_mf_enr = gp.enrichr(gene_list=up_genes, gene_sets=go_mf_set, 
                                  organism=organism, outdir=None)
            if go_mf_enr.results is not None and len(go_mf_enr.results) > 0:
                go_mf_enr.results['gene_set'] = 'GO_MF'
                go_mf_enr.results['direction'] = 'up'
                all_enrichment_results.append(go_mf_enr.results)
        except Exception as e:
            logger.warning(f"GO MF富集分析失败: {e}")
    
    # 对下调基因进行富集
    if len(down_genes) >= 10:
        logger.info("对下调基因进行富集分析")
        
        # KEGG富集
        try:
            kegg_enr = gp.enrichr(gene_list=down_genes, gene_sets=kegg_set, 
                                 organism=organism, outdir=None)
            if kegg_enr.results is not None and len(kegg_enr.results) > 0:
                kegg_enr.results['gene_set'] = 'KEGG'
                kegg_enr.results['direction'] = 'down'
                all_enrichment_results.append(kegg_enr.results)
        except Exception as e:
            logger.warning(f"KEGG富集分析失败: {e}")
        
        # GO BP富集
        try:
            go_bp_enr = gp.enrichr(gene_list=down_genes, gene_sets=go_bp_set, 
                                  organism=organism, outdir=None)
            if go_bp_enr.results is not None and len(go_bp_enr.results) > 0:
                go_bp_enr.results['gene_set'] = 'GO_BP'
                go_bp_enr.results['direction'] = 'down'
                all_enrichment_results.append(go_bp_enr.results)
        except Exception as e:
            logger.warning(f"GO BP富集分析失败: {e}")
        
        # GO CC富集
        try:
            go_cc_enr = gp.enrichr(gene_list=down_genes, gene_sets=go_cc_set, 
                                  organism=organism, outdir=None)
            if go_cc_enr.results is not None and len(go_cc_enr.results) > 0:
                go_cc_enr.results['gene_set'] = 'GO_CC'
                go_cc_enr.results['direction'] = 'down'
                all_enrichment_results.append(go_cc_enr.results)
        except Exception as e:
            logger.warning(f"GO CC富集分析失败: {e}")
        
        # GO MF富集
        try:
            go_mf_enr = gp.enrichr(gene_list=down_genes, gene_sets=go_mf_set, 
                                  organism=organism, outdir=None)
            if go_mf_enr.results is not None and len(go_mf_enr.results) > 0:
                go_mf_enr.results['gene_set'] = 'GO_MF'
                go_mf_enr.results['direction'] = 'down'
                all_enrichment_results.append(go_mf_enr.results)
        except Exception as e:
            logger.warning(f"GO MF富集分析失败: {e}")
    
    # 保存所有结果
    if all_enrichment_results:
        # 合并所有结果
        all_enrichment = pd.concat(all_enrichment_results, ignore_index=True)
        all_enrichment.to_csv(os.path.join(output_dir, 'enrichment_results.csv'), index=False)
        
        # 分别提取KEGG、GO_BP、GO_CC、GO_MF结果
        kegg_results = all_enrichment[all_enrichment['gene_set'] == 'KEGG']
        go_bp_results = all_enrichment[all_enrichment['gene_set'] == 'GO_BP']
        go_cc_results = all_enrichment[all_enrichment['gene_set'] == 'GO_CC']
        go_mf_results = all_enrichment[all_enrichment['gene_set'] == 'GO_MF']
        
        # 分别保存
        if len(kegg_results) > 0:
            kegg_results.to_csv(os.path.join(output_dir, 'kegg_enrichment_results.csv'), index=False)
            # 绘制KEGG气泡图
            plot_bubble_dynamic(kegg_results, output_dir, 'KEGG')
        
        if len(go_bp_results) > 0:
            go_bp_results.to_csv(os.path.join(output_dir, 'go_bp_enrichment_results.csv'), index=False)
            # 绘制GO_BP气泡图
            plot_bubble_dynamic(go_bp_results, output_dir, 'GO_BP')
        
        if len(go_cc_results) > 0:
            go_cc_results.to_csv(os.path.join(output_dir, 'go_cc_enrichment_results.csv'), index=False)
            # 绘制GO_CC气泡图
            plot_bubble_dynamic(go_cc_results, output_dir, 'GO_CC')
        
        if len(go_mf_results) > 0:
            go_mf_results.to_csv(os.path.join(output_dir, 'go_mf_enrichment_results.csv'), index=False)
            # 绘制GO_MF气泡图
            plot_bubble_dynamic(go_mf_results, output_dir, 'GO_MF')
        
        # 绘制柱状图（所有通路）
        plot_enrichment_barplot(all_enrichment, output_dir)
        
        logger.info(f"富集分析完成，保存了 KEGG、GO_BP、GO_CC、GO_MF 结果")
    else:
        logger.warning("没有找到显著的富集结果")
    
    return all_enrichment_results

def plot_enrichment_barplot(enrichment_df, output_dir, top_n=10):
    """
    绘制富集分析柱状图 - 优化版
    优化：
    1. Y轴间距拉开，防止标签重叠
    2. 上下调图形高度保持一致（同一条目数时）
    3. 通路名很长时动态放大图形宽度
    """
    if enrichment_df is None or len(enrichment_df) == 0:
        return
    
    # 筛选显著结果
    sig_enrich = enrichment_df[enrichment_df['Adjusted P-value'] < 0.05]
    if len(sig_enrich) == 0:
        return
    
    # 分别处理上调、下调
    for direction in ['up', 'down']:
        subset = sig_enrich[sig_enrich['direction'] == direction]
        if len(subset) == 0:
            continue
        
        # 取前top_n个
        subset = subset.nsmallest(top_n, 'Adjusted P-value')
        n_pathways = len(subset)
        
        # ============================================================
        # 动态高度：根据条目数 + 固定边距
        # 使用 y_spacing_factor 增加垂直间距
        # ============================================================
        y_spacing_factor = 1.5  # 拉开间距
        fig_height = max(10, n_pathways * y_spacing_factor + 4)
        
        # ============================================================
        # 动态宽度：考虑通路名称长度
        # ============================================================
        max_term_len = max(len(str(t)) for t in subset['Term'])
        fig_width = max(20, max_term_len * 0.25 + 12)  # 字符数 * 系数
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        # Y坐标拉开间距
        y_pos = np.arange(len(subset)) * y_spacing_factor
        
        bars = ax.barh(y_pos, -np.log10(subset['Adjusted P-value'].values),
                       color='#D62728' if direction == 'up' else '#1F77B4')
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(subset['Term'], fontsize=FONT_SMALL, fontweight='heavy')
        ax.set_xlabel('-Log10 Adjusted P-value', fontsize=FONT_AXIS_LABEL, fontweight='heavy')
        ax.set_title(f'{direction.upper()} Regulated Genes - Pathway Enrichment',
                     fontsize=FONT_TITLE, fontweight='heavy', pad=30)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(3.0)
        ax.spines['bottom'].set_linewidth(3.0)
        
        ax.grid(True, axis='x', linestyle='--', alpha=0.3, linewidth=2.5)
        ax.set_axisbelow(True)
        
        # 设置 Y 轴范围留出额外空间
        ax.set_ylim(y_pos[-1] + y_spacing_factor * 0.8, y_pos[0] - y_spacing_factor * 0.8)
        
        plt.tight_layout(pad=4.0)
        
        save_path = os.path.join(output_dir, f'enrichment_barplot_{direction}.png')
        save_pdf = os.path.join(output_dir, f'enrichment_barplot_{direction}.pdf')
        fig.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        fig.savefig(save_pdf, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        plt.close(fig)

def plot_bubble_dynamic(enrichment_df, output_dir, category_name, top_n=15):
    """
    绘制富集气泡图 - 优化版
    优化：
    1. X轴范围更宽（Gene Ratio 系数增大至 3.5）
    2. 气泡大小和Y轴间距动态调整，防止重叠
    3. 通路名很长时动态放大图形宽度和高度
    """
    if enrichment_df is None or len(enrichment_df) == 0:
        return
    
    # 筛选显著结果
    sig_enrich = enrichment_df[enrichment_df['Adjusted P-value'] < 0.05]
    if len(sig_enrich) == 0:
        logger = logging.getLogger('Pipeline')
        logger.warning(f"没有显著的{category_name}富集结果")
        return
    
    for direction in ['up', 'down']:
        subset = sig_enrich[sig_enrich['direction'] == direction]
        if len(subset) == 0:
            continue
        
        subset = subset.nsmallest(top_n, 'Adjusted P-value')
        n_pathways = len(subset)
        
        subset['Gene_Ratio'] = subset['Overlap'].apply(lambda x: float(x.split('/')[0]) / float(x.split('/')[1]))
        
        # ============================================================
        # 动态宽度和高度
        # ============================================================
        max_term_len = max(len(str(t)) for t in subset['Term'])
        # 宽度：基础宽度 + 每个通路的增量 + 最长名字的额外宽度
        fig_width = max(22, n_pathways * 1.2 + 14 + max_term_len * 0.18)
        fig_width = min(fig_width, 50)  # 上限防止过于夸张
        
        # 高度：条目间距增大
        y_spacing = 1.8
        fig_height = max(12, n_pathways * y_spacing + 6)
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        # Y坐标拉开间距
        y_positions = np.arange(len(subset)) * y_spacing
        
        # 气泡大小
        log_p_values = -np.log10(subset['Adjusted P-value'].values)
        min_size = 80
        max_size = 500
        if log_p_values.max() > log_p_values.min():
            sizes = min_size + (log_p_values - log_p_values.min()) / (log_p_values.max() - log_p_values.min()) * (max_size - min_size)
        else:
            sizes = np.full(len(log_p_values), min_size + (max_size - min_size) / 2)
        
        scatter = ax.scatter(
            subset['Gene_Ratio'],
            y_positions,
            s=sizes,
            c=log_p_values,
            cmap='Reds' if direction == 'up' else 'Blues',
            alpha=0.8,
            edgecolors='black',
            linewidth=1.5
        )
        
        ax.set_yticks(y_positions)
        ax.set_yticklabels(subset['Term'], fontsize=FONT_SMALL, fontweight='heavy')
        ax.set_xlabel('Gene Ratio', fontsize=FONT_AXIS_LABEL, fontweight='heavy')
        
        if category_name == 'KEGG':
            ylabel = 'KEGG Pathway'
        elif category_name == 'GO_BP':
            ylabel = 'GO Biological Process'
        elif category_name == 'GO_CC':
            ylabel = 'GO Cellular Component'
        elif category_name == 'GO_MF':
            ylabel = 'GO Molecular Function'
        else:
            ylabel = 'Pathway'
        ax.set_ylabel(ylabel, fontsize=FONT_AXIS_LABEL, fontweight='heavy')
        ax.set_title(f'{direction.upper()} Regulated Genes - {category_name} Enrichment',
                     fontsize=FONT_TITLE, fontweight='heavy', pad=30)
        
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('-Log10 Adjusted P-value', fontsize=FONT_COLORBAR, fontweight='heavy', labelpad=20)
        cbar.ax.tick_params(labelsize=FONT_SMALL, width=3.0, length=8)
        cbar.outline.set_linewidth(3.0)
        
        # ============================================================
        # X轴范围极大放宽
        # ============================================================
        max_gene_ratio = subset['Gene_Ratio'].max()
        x_lim_factor = 3.5
        x_max = max(max_gene_ratio * x_lim_factor, 0.3)
        ax.set_xlim(0, x_max)
        
        # Y轴范围留出边距
        ax.set_ylim(y_positions[-1] + y_spacing * 0.8, y_positions[0] - y_spacing * 0.8)
        
        ax.grid(True, axis='x', linestyle='--', alpha=0.3, linewidth=2.5)
        ax.set_axisbelow(True)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(3.0)
        ax.spines['bottom'].set_linewidth(3.0)
        
        plt.tight_layout(pad=4.0)
        
        category_lower = category_name.lower()
        save_path = os.path.join(output_dir, f'{category_lower}_bubble_{direction}.png')
        save_pdf = os.path.join(output_dir, f'{category_lower}_bubble_{direction}.pdf')
        fig.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        fig.savefig(save_pdf, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        plt.close(fig)
        
        logger = logging.getLogger('Pipeline')
        logger.info(f"{category_name}气泡图已保存: {save_path}")

# ============================================================
# 主流水线（包含富集分析）
# ============================================================

def run_integrated_pipeline_with_enrichment(adata_path, output_base_dir='./Integrated_Analysis_Results',
                                             group_col='group_sub', celltype_col='singler_best',
                                             sample_col='sample', run_deg=True, run_enrichment=True):
    """
    运行完整的整合分析流水线 - 包含富集分析
    """
    os.makedirs(output_base_dir, exist_ok=True)
    log_file = os.path.join(output_base_dir, 'pipeline.log')

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger('Pipeline')

    logger.info("="*60)
    logger.info("Integrated Analysis Pipeline Started (with Enrichment)")
    logger.info("Includes: Cell Proportion + ROE + DEG (scanpy) + Enrichment")
    logger.info("FIX 1: Boxplots Y-axis starts at 0, with ns annotation")
    logger.info("FIX 2: No FDR correction, direct p-value")
    logger.info("FIX 3: Volcano plots optimized (no gene labels, legend on top-right, 0 in center)")
    logger.info("FIX 4: Enrichment: KEGG + GO_BP + GO_CC + GO_MF with optimized bubble plots")
    logger.info("FIX 5: Bubble plot X-axis extended, Y spacing increased, dynamic size for long terms")
    logger.info("FIX 6: Bar plot Y spacing unified, long term width adjusted")
    logger.info("="*60)

    # 加载数据
    logger.info(f"Loading data from: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    logger.info(f"Data loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    required_cols = [sample_col, group_col, celltype_col]
    for col in required_cols:
        if col not in adata.obs.columns:
            raise ValueError(f"Required column '{col}' not found in adata.obs")

    # 保存原始计数
    if adata.raw is None:
        adata.raw = adata.copy()
        logger.info("保存原始计数到adata.raw")

    # ========================================
    # Part 1: 细胞比例分析
    # ========================================
    logger.info("\n" + "="*40)
    logger.info("Part 1: Cell Proportion Analysis")
    logger.info("="*40)

    prop_output_dir = os.path.join(output_base_dir, '01_Cell_Proportions')
    os.makedirs(prop_output_dir, exist_ok=True)

    prop_df, obs_df = calculate_cell_proportions(adata, group_col, celltype_col, sample_col)
    prop_df.to_csv(os.path.join(prop_output_dir, 'cell_proportions.csv'), index=False)

    stat_results, sig_celltypes = statistical_test_proportions_no_fdr(prop_df, group_col, celltype_col)
    stat_results.to_csv(os.path.join(prop_output_dir, 'statistical_test_results.csv'), index=False)

    with open(os.path.join(prop_output_dir, 'significant_celltypes.txt'), 'w') as f:
        for ct in sorted(sig_celltypes):
            f.write(f"{ct}\n")

    logger.info(f"显著差异细胞类型 ({len(sig_celltypes)}): {sorted(sig_celltypes)}")

    logger.info("Plotting boxplots...")
    all_boxplots, sig_boxplots = plot_all_boxplots_fixed_v2(
        prop_df, prop_output_dir, group_col, stat_results, sig_celltypes
    )

    # ========================================
    # Part 2: ROE分析
    # ========================================
    logger.info("\n" + "="*40)
    logger.info("Part 2: ROE Analysis")
    logger.info("="*40)

    roe_output_dir = os.path.join(output_base_dir, '02_ROE_Analysis')
    os.makedirs(roe_output_dir, exist_ok=True)

    roe_df, log2fc_df, crosstab = perform_roe_analysis(obs_df, group_col, celltype_col)
    roe_df.to_csv(os.path.join(roe_output_dir, 'roe_results.csv'), index=False)
    log2fc_df.to_csv(os.path.join(roe_output_dir, 'log2fc_results.csv'), index=False)
    crosstab.to_csv(os.path.join(roe_output_dir, 'crosstab.csv'))

    plot_roe_heatmap_cns(log2fc_df, roe_output_dir, title="Cell Type Enrichment (ROE Analysis)")

    if len(sig_celltypes) > 0:
        plot_roe_heatmap_cns(log2fc_df, roe_output_dir,
                              title="Significant Cell Types - ROE Analysis",
                              sig_celltypes=sig_celltypes)

    # ========================================
    # Part 3: DEG分析 - 为所有细胞类型
    # ========================================
    if run_deg:
        logger.info("\n" + "="*40)
        logger.info("Part 3: DEG Analysis (scanpy rank_genes_groups)")
        logger.info("Running DEG for ALL cell types")
        logger.info("="*40)

        deg_output_dir = os.path.join(output_base_dir, '03_DEG_Analysis')
        os.makedirs(deg_output_dir, exist_ok=True)

        groups = sorted(adata.obs[group_col].unique())
        comparisons = list(combinations(groups, 2))

        all_deg_results = []
        
        # 获取所有细胞类型
        all_celltypes = sorted(adata.obs[celltype_col].unique())
        logger.info(f"Running DEG for all {len(all_celltypes)} cell types")

        for celltype in all_celltypes:
            ct_output_dir = os.path.join(deg_output_dir, safe_name(celltype))
            os.makedirs(ct_output_dir, exist_ok=True)

            ct_logger = logging.getLogger(f'DEG_{celltype}')
            ct_logger.setLevel(logging.INFO)
            fh = logging.FileHandler(os.path.join(ct_output_dir, 'deg_analysis.log'), mode='w')
            fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            ct_logger.addHandler(fh)

            for g1, g2 in comparisons:
                comp_name = f"{g1}_vs_{g2}"
                comp_dir = os.path.join(ct_output_dir, comp_name)
                os.makedirs(comp_dir, exist_ok=True)

                ct_logger.info(f"DEG: {celltype} - {comp_name}")

                de_results = perform_dea_scanpy(
                    adata, celltype, g1, g2,
                    celltype_col=celltype_col,
                    group_col=group_col
                )

                if de_results is not None and len(de_results) > 0:
                    de_results.to_csv(os.path.join(comp_dir, 'deg_results.csv'), index=False)
                    
                    # 使用优化后的火山图
                    plot_volcano_cns_optimized(de_results, comp_dir, comp_name, celltype)

                    up_genes = de_results[
                        (de_results['pvals_adj'] < 0.05) &
                        (de_results['logfoldchanges'] > 0.25)
                    ]['names'].tolist()

                    down_genes = de_results[
                        (de_results['pvals_adj'] < 0.05) &
                        (de_results['logfoldchanges'] < -0.25)
                    ]['names'].tolist()

                    with open(os.path.join(comp_dir, 'up_genes.txt'), 'w') as f:
                        f.write('\n'.join(up_genes))
                    with open(os.path.join(comp_dir, 'down_genes.txt'), 'w') as f:
                        f.write('\n'.join(down_genes))

                    all_deg_results.append(de_results)
                    ct_logger.info(f"  Up: {len(up_genes)}, Down: {len(down_genes)}")
                    
                    # ========================================
                    # Part 4: 富集分析
                    # ========================================
                    if run_enrichment and (len(up_genes) >= 10 or len(down_genes) >= 10):
                        enrichment_dir = os.path.join(comp_dir, 'Enrichment')
                        os.makedirs(enrichment_dir, exist_ok=True)
                        
                        ct_logger.info(f"Running enrichment analysis for {celltype} - {comp_name}")
                        perform_enrichment_analysis(de_results, enrichment_dir)
                        
                else:
                    ct_logger.warning(f"  No DEG results (insufficient cells or other error)")

            ct_logger.removeHandler(fh)

        if all_deg_results:
            all_deg_df = pd.concat(all_deg_results, ignore_index=True)
            all_deg_df.to_csv(os.path.join(deg_output_dir, 'all_deg_results.csv'), index=False)
            logger.info(f"所有DEG结果汇总: {len(all_deg_df)} 条记录，来自 {len(all_celltypes)} 个细胞类型")
            
            # 生成DEG汇总统计
            deg_summary = all_deg_df.groupby(['celltype', 'comparison']).agg({
                'pvals_adj': lambda x: (x < 0.05).sum(),
                'logfoldchanges': lambda x: ((x > 0.25) | (x < -0.25)).sum()
            }).reset_index()
            deg_summary.columns = ['celltype', 'comparison', 'sig_genes', 'fc_genes']
            deg_summary.to_csv(os.path.join(deg_output_dir, 'deg_summary.csv'), index=False)
        else:
            logger.warning("No DEG results generated")

    # ========================================
    # 汇总报告
    # ========================================
    logger.info("\n" + "="*40)
    logger.info("Generating Summary Report")
    logger.info("="*40)

    summary_path = os.path.join(output_base_dir, 'analysis_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("="*60 + "\n")
        f.write("Integrated Analysis Pipeline - Summary Report\n")
        f.write("="*60 + "\n\n")
        f.write(f"Date: {pd.Timestamp.now()}\n")
        f.write(f"Input: {adata_path}\n")
        f.write(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}\n")
        f.write(f"Group column: {group_col}\n")
        f.write(f"Cell type column: {celltype_col}\n")
        f.write(f"Sample column: {sample_col}\n\n")
        f.write(f"Groups: {sorted(adata.obs[group_col].unique())}\n")
        f.write(f"Cell types: {adata.obs[celltype_col].nunique()}\n\n")
        f.write("Font: Arial, PDF fonttype: 42\n")
        f.write(f"Title: {FONT_TITLE}pt, Axis: {FONT_AXIS_LABEL}pt, Tick: {FONT_TICK}pt\n\n")
        f.write(f"Significant cell types (proportion): {len(sig_celltypes)}\n")
        for ct in sorted(sig_celltypes):
            f.write(f"  - {ct}\n")
        f.write("\nFIXES:\n")
        f.write("1. Boxplots: Y-axis starts at 0, with ns annotation\n")
        f.write("2. Statistics: No FDR correction, direct p-value\n")
        f.write("3. Volcano: No gene labels, legend on top-right, 0 in center\n")
        f.write("4. Enrichment: KEGG + GO_BP + GO_CC + GO_MF with optimized bubble plots\n")
        f.write("5. Bubble plot: X-axis extended, Y spacing increased, dynamic size for long terms\n")
        f.write("6. Bar plot: Y spacing unified, long term width adjusted\n")
        f.write("\nOutput:\n")
        f.write("01_Cell_Proportions/\n")
        f.write("  ├── boxplots_all/\n")
        f.write("  ├── sig/\n")
        f.write("  └── summary_significant_boxplots.png/pdf\n")
        f.write("02_ROE_Analysis/\n")
        f.write("  └── roe_heatmap.png/pdf\n")
        f.write("03_DEG_Analysis/\n")
        f.write("  └── {ALL celltypes}/{comparison}/\n")
        f.write("      ├── deg_results.csv\n")
        f.write("      ├── volcano_*.png/pdf (optimized)\n")
        f.write("      ├── up_genes.txt\n")
        f.write("      ├── down_genes.txt\n")
        f.write("      └── Enrichment/\n")
        f.write("          ├── enrichment_results.csv\n")
        f.write("          ├── kegg_enrichment_results.csv\n")
        f.write("          ├── go_bp_enrichment_results.csv\n")
        f.write("          ├── go_cc_enrichment_results.csv\n")
        f.write("          ├── go_mf_enrichment_results.csv\n")
        f.write("          ├── kegg_bubble_{up,down}.png/pdf\n")
        f.write("          ├── go_bp_bubble_{up,down}.png/pdf\n")
        f.write("          ├── go_cc_bubble_{up,down}.png/pdf\n")
        f.write("          └── go_mf_bubble_{up,down}.png/pdf\n")

    logger.info(f"\nSummary: {summary_path}")
    logger.info("\n" + "="*60)
    logger.info("Pipeline completed successfully!")
    logger.info(f"Results: {output_base_dir}")
    logger.info("="*60)

    return output_base_dir

# ============================================================
# 主程序入口
# ============================================================

if __name__ == "__main__":
    ADATA_PATH = "M2_C3_All_K_Rtrajectory_PostR_vs_PreR.h5ad"
    OUTPUT_DIR = "./Integrated_Analysis_Results_Final"
    GROUP_COL = "group_sub"
    CELLTYPE_COL = "singler_best"
    SAMPLE_COL = "sample"
    RUN_DEG = True
    RUN_ENRICHMENT = True

    if not os.path.exists(ADATA_PATH):
        print(f"Error: Input file not found: {ADATA_PATH}")
        sys.exit(1)

    run_integrated_pipeline_with_enrichment(
        adata_path=ADATA_PATH,
        output_base_dir=OUTPUT_DIR,
        group_col=GROUP_COL,
        celltype_col=CELLTYPE_COL,
        sample_col=SAMPLE_COL,
        run_deg=RUN_DEG,
        run_enrichment=RUN_ENRICHMENT
    )
