import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import mannwhitneyu
# 不再需要 multipletests
import scanpy as sc
from typing import Dict
import os
import warnings
warnings.filterwarnings('ignore')

# ========== 全局绘图设置（Nature 风格） ==========
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 8
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 0.6
plt.rcParams['ytick.major.width'] = 0.6
plt.rcParams['xtick.major.size'] = 3
plt.rcParams['ytick.major.size'] = 3
plt.rcParams['lines.markersize'] = 3
sns.set_context("paper", font_scale=1.2, rc={"lines.linewidth": 0.8})

PALETTE_CANCER = {'MPC': '#E64B35', 'TNBC': '#4DBBD5'}

# ================== 数据准备 ==================
def prepare_proportions(adata: sc.AnnData,
                        sample_col: str = 'sample',
                        cancer_col: str = 'cancer_type',
                        celltype_col: str = 'annotation') -> pd.DataFrame:
    obs = adata.obs[[sample_col, cancer_col, celltype_col]].copy()
    obs = obs.dropna(subset=[celltype_col])
    count = obs.groupby([sample_col, cancer_col, celltype_col]).size().reset_index(name='count')
    total_per_sample = count.groupby([sample_col, cancer_col])['count'].sum().reset_index(name='total')
    prop = count.merge(total_per_sample, on=[sample_col, cancer_col])
    prop['proportion'] = prop['count'] / prop['total']
    prop = prop.rename(columns={celltype_col: 'cell_type', cancer_col: 'cancer_type', sample_col: 'sample'})
    return prop[['sample', 'cancer_type', 'cell_type', 'proportion']]

# ================== 统计检验（原始 p 值，无校正） ==================
def compute_significance(prop_df: pd.DataFrame,
                         group_col: str = 'cancer_type',
                         value_col: str = 'proportion',
                         celltype_col: str = 'cell_type') -> Dict[str, Dict]:
    groups = prop_df[group_col].unique()
    if len(groups) != 2:
        raise ValueError("分组变量必须有两个水平")
    g1, g2 = groups[0], groups[1]

    results = {}
    celltypes = prop_df[celltype_col].unique()

    for ct in celltypes:
        vals1 = prop_df.loc[(prop_df[celltype_col] == ct) & (prop_df[group_col] == g1), value_col]
        vals2 = prop_df.loc[(prop_df[celltype_col] == ct) & (prop_df[group_col] == g2), value_col]
        if len(vals1) < 3 or len(vals2) < 3:
            results[ct] = {'p_val': np.nan, 'signif': 'ns', 'U': np.nan}
            continue
        stat, p = mannwhitneyu(vals1, vals2, alternative='two-sided')
        results[ct] = {'p_val': p, 'U': stat}
        # 直接使用原始 p 值，不进行 FDR 校正
        if p < 0.001:
            results[ct]['signif'] = '***'
        elif p < 0.01:
            results[ct]['signif'] = '**'
        elif p < 0.05:
            results[ct]['signif'] = '*'
        else:
            results[ct]['signif'] = 'ns'
    return results

# ================== 原图 A：小提琴 + 显著性（所有细胞类型） ==================
def plot_violin_significance(prop_df, signif_dict, save_path, max_cols=6):
    celltypes = sorted(prop_df['cell_type'].unique())
    n_types = len(celltypes)
    n_cols = min(max_cols, n_types)
    n_rows = int(np.ceil(n_types / n_cols))
    figsize = (2.8 * n_cols, 2.5 * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    for idx, ct in enumerate(celltypes):
        ax = axes[idx]
        sub = prop_df[prop_df['cell_type'] == ct]
        sns.violinplot(data=sub, x='cancer_type', y='proportion',
                       palette=PALETTE_CANCER, inner='box', linewidth=0.6, scale='width', ax=ax, cut=0)
        sns.stripplot(data=sub, x='cancer_type', y='proportion',
                      color='black', alpha=0.5, size=3, jitter=True, ax=ax)
        ax.set_title(ct, fontsize=9, pad=4)
        ax.set_xlabel('')
        ax.set_ylabel('Proportion' if idx % n_cols == 0 else '')
        ax.set_ylim(-0.02, 1.02)
        signif = signif_dict.get(ct, {}).get('signif', 'ns')
        if signif != 'ns':
            y_max = sub['proportion'].max() * 1.05 if len(sub) else 1
            h = y_max * 0.05
            x1, x2 = 0, 1
            ax.plot([x1, x1, x2, x2], [y_max, y_max+h, y_max+h, y_max], lw=0.6, c='black')
            ax.text((x1+x2)/2, y_max+h, signif, ha='center', va='bottom', fontsize=10)
    for j in range(idx+1, len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout(pad=1.2)
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"[Fig A] 已保存: {save_path}")

# ================== 原图 B：堆叠柱状图 ==================
def plot_stacked_bar(prop_df, save_path, figsize=(12,6)):
    wide = prop_df.pivot_table(index=['sample','cancer_type'], columns='cell_type', values='proportion', fill_value=0).reset_index()
    wide = wide.sort_values(['cancer_type','sample'])
    sample_order = wide['sample'].tolist()
    celltypes = sorted(prop_df['cell_type'].unique())
    palette_cell = sns.color_palette("husl", len(celltypes))
    col_dict = dict(zip(celltypes, palette_cell))

    fig, ax = plt.subplots(figsize=figsize)
    bottom = np.zeros(len(wide))
    for ct in celltypes:
        vals = wide[ct].values
        ax.bar(range(len(wide)), vals, bottom=bottom, color=col_dict[ct], label=ct, width=0.8)
        bottom += vals
    ax.set_xticks(range(len(wide)))
    labels = [f"{row['sample']}\n({row['cancer_type']})" for _, row in wide.iterrows()]
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel('Proportion')
    ax.legend(title='Cell type', bbox_to_anchor=(1.01,1), loc='upper left', fontsize=7)
    cancer_types = wide['cancer_type'].unique()
    for ct in cancer_types:
        idxs = wide[wide['cancer_type']==ct].index
        if len(idxs)==0: continue
        xmin, xmax = min(idxs)-0.4, max(idxs)+0.4
        ax.axvspan(xmin, xmax, facecolor=PALETTE_CANCER.get(ct,'#888'), alpha=0.1, zorder=0)
        ax.text((xmin+xmax)/2, 1.05, ct, ha='center', fontsize=9, color=PALETTE_CANCER.get(ct,'#888'), weight='bold')
    ax.set_xlim(-0.5, len(wide)-0.5)
    ax.set_ylim(0,1.15)
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"[Fig B] 已保存: {save_path}")

# ================== 原图 C：平均占比差异热图 ==================
def plot_heatmap_diff(prop_df, save_path, figsize=(6,4)):
    mean_prop = prop_df.groupby(['cell_type','cancer_type'])['proportion'].mean().reset_index()
    pivot = mean_prop.pivot(index='cell_type', columns='cancer_type', values='proportion')
    for ct in ['MPC','TNBC']:
        if ct not in pivot.columns: pivot[ct]=0
    pivot = pivot[['MPC','TNBC']]
    pivot['diff (MPC-TNBC)'] = pivot['MPC'] - pivot['TNBC']
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdBu_r', center=0, linewidths=0.5, ax=ax, cbar_kws={'label':'Proportion / Difference'})
    ax.set_title('Mean cell-type proportion and difference')
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"[Fig C] 已保存: {save_path}")

# ================== 原图 D：UMAP 分面（若存在） ==================
def plot_umap_facet(adata, celltype_col='annotation', cancer_col='cancer_type', basis='X_umap', save_path='FigD_UMAP.pdf'):
    if basis not in adata.obsm_keys():
        print("[跳过] UMAP 坐标不存在")
        return
    celltypes = sorted(adata.obs[celltype_col].dropna().unique())
    if len(celltypes) > 12:
        top = adata.obs[celltype_col].value_counts().index[:12]
        celltypes = sorted(top)
    cancer_types = sorted(adata.obs[cancer_col].dropna().unique())
    n_rows, n_cols = len(cancer_types), len(celltypes)
    figsize = (2.2*n_cols, 2.2*n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    umap_coords = adata.obsm[basis]
    for i, ct_cancer in enumerate(cancer_types):
        for j, ct_cell in enumerate(celltypes):
            ax = axes[i][j]
            idx = (adata.obs[cancer_col]==ct_cancer) & (adata.obs[celltype_col]==ct_cell)
            ax.scatter(umap_coords[:,0], umap_coords[:,1], s=0.5, c='lightgrey', rasterized=True)
            ax.scatter(umap_coords[idx,0], umap_coords[idx,1], s=1, c=PALETTE_CANCER[ct_cancer], rasterized=True)
            ax.set_title(f'{ct_cancer} | {ct_cell}', fontsize=7, pad=2)
            ax.set_xticks([]); ax.set_yticks([])
            if j==0: ax.set_ylabel(ct_cancer, fontsize=8)
            if i==0: ax.set_xlabel(ct_cell, fontsize=7, rotation=0, ha='center')
    plt.tight_layout(pad=0.8)
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"[Fig D] 已保存: {save_path}")

# ================== 新增图 E：显著细胞类型小提琴图（使用原始 p 值） ==================
def plot_significant_violin(prop_df, signif_dict, save_path, max_cols=4):
    # 筛选显著（signif != 'ns'）的细胞类型
    sig_types = [ct for ct, info in signif_dict.items() if info.get('signif', 'ns') != 'ns']
    if not sig_types:
        print("没有显著差异的细胞类型，跳过图 E")
        return
    sig_types = sorted(sig_types)
    n = len(sig_types)
    n_cols = min(max_cols, n)
    n_rows = int(np.ceil(n / n_cols))
    figsize = (3.2 * n_cols, 2.8 * n_rows)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    for idx, ct in enumerate(sig_types):
        ax = axes[idx]
        sub = prop_df[prop_df['cell_type'] == ct]
        sns.violinplot(data=sub, x='cancer_type', y='proportion',
                       palette=PALETTE_CANCER, inner='box', linewidth=0.8, scale='width', ax=ax, cut=0)
        sns.stripplot(data=sub, x='cancer_type', y='proportion',
                      color='black', alpha=0.6, size=4, jitter=True, ax=ax)
        ax.set_title(ct, fontsize=10, pad=6)
        ax.set_xlabel('')
        ax.set_ylabel('Proportion' if idx % n_cols == 0 else '')
        p_val = signif_dict[ct]['p_val']
        signif_mark = signif_dict[ct]['signif']
        p_text = f"p = {p_val:.2e}" if p_val < 0.001 else f"p = {p_val:.3f}"
        ax.text(0.5, 0.95, f"{signif_mark} {p_text}", transform=ax.transAxes,
                ha='center', va='top', fontsize=9, weight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
        ax.set_ylim(-0.02, max(sub['proportion'].max()*1.1, 0.2))
    for j in range(idx+1, len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout(pad=1.5)
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"[Fig E] 显著细胞类型小提琴图已保存: {save_path}")

# ================== 新增图 F：差异火山图（使用原始 p 值） ==================
def plot_volcano_diff(prop_df, signif_dict, save_path, figsize=(8,6)):
    means = prop_df.groupby(['cell_type','cancer_type'])['proportion'].mean().unstack(fill_value=0)
    diff = means['MPC'] - means['TNBC']
    celltypes = diff.index.tolist()
    neg_log10_p = []
    sig_labels = []
    for ct in celltypes:
        p_val = signif_dict[ct]['p_val']
        if np.isnan(p_val) or p_val == 0:
            neg_log10_p.append(np.nan)
            sig_labels.append('ns')
        else:
            neg_log10_p.append(-np.log10(p_val))
            sig_labels.append(signif_dict[ct]['signif'])
    df_volcano = pd.DataFrame({'cell_type': celltypes, 'diff': diff.values,
                               '-log10(p)': neg_log10_p, 'signif': sig_labels})
    significant = df_volcano['signif'] != 'ns'

    fig, ax = plt.subplots(figsize=figsize)
    ns = df_volcano[~significant]
    ax.scatter(ns['diff'], ns['-log10(p)'], color='grey', alpha=0.5, s=60, edgecolor='white')
    sig = df_volcano[significant]
    colors = ['#E64B35' if x > 0 else '#4DBBD5' for x in sig['diff']]
    ax.scatter(sig['diff'], sig['-log10(p)'], color=colors, s=100, edgecolor='black', linewidth=0.5)
    for _, row in sig.iterrows():
        ax.annotate(row['cell_type'], (row['diff'], row['-log10(p)']),
                    textcoords="offset points", xytext=(5,5), ha='left', fontsize=8,
                    arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.axhline(-np.log10(0.05), color='grey', linewidth=0.6, linestyle=':', label='p=0.05')
    ax.set_xlabel('Mean proportion difference (MPC − TNBC)', fontsize=11)
    ax.set_ylabel('−log₁₀ (p value)', fontsize=11)
    ax.legend()
    sns.despine(ax=ax)
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"[Fig F] 差异火山图已保存: {save_path}")

# ================== 统计表格导出（无校正 p 值列） ==================
def export_statistics_table(prop_df, signif_dict, save_path):
    stats = prop_df.groupby(['cell_type','cancer_type'])['proportion'].agg(['mean','std']).unstack()
    stats.columns = ['_'.join(col).strip() for col in stats.columns.values]
    for ct in ['MPC','TNBC']:
        for stat in ['mean','std']:
            col_name = f'{stat}_{ct}'
            if col_name not in stats.columns:
                stats[col_name] = np.nan
    stats['mean_diff'] = stats['mean_MPC'] - stats['mean_TNBC']

    stats['U_statistic'] = np.nan
    stats['p_value'] = np.nan
    stats['significance'] = ''
    for ct in stats.index:
        if ct in signif_dict:
            info = signif_dict[ct]
            stats.loc[ct, 'U_statistic'] = info.get('U', np.nan)
            stats.loc[ct, 'p_value'] = info.get('p_val', np.nan)
            stats.loc[ct, 'significance'] = info.get('signif', '')
    cols = ['mean_MPC','std_MPC','mean_TNBC','std_TNBC','mean_diff',
            'U_statistic','p_value','significance']
    stats = stats[cols]
    stats.to_csv(save_path, float_format='%.6f')
    print(f"[表格] 统计结果已导出: {save_path}")

# ================== 主流程 ==================
def main(adata_path='final_all.h5ad',
         sample_col='patient',          # 维度为 patient
         cancer_col='cancer_type',
         celltype_col='annotation',
         save_dir='./figures',
         run_umap=True):
    os.makedirs(save_dir, exist_ok=True)

    adata = sc.read(adata_path)
    print(f"原始细胞数: {adata.n_obs}")
    adata = adata[adata.obs[cancer_col].isin(['MPC','TNBC'])].copy()
    print(f"MPC/TNBC 细胞数: {adata.n_obs}")

    prop_df = prepare_proportions(adata, sample_col=sample_col, cancer_col=cancer_col, celltype_col=celltype_col)
    signif_dict = compute_significance(prop_df)
    print("显著性检验完成（基于患者，原始 p 值，无校正）")

    # 输出显著差异的细胞类型
    sig_list = [(ct, info['p_val'], info['signif']) for ct, info in signif_dict.items() if info['signif'] != 'ns']
    if sig_list:
        print("\n原始 p < 0.05 的细胞类型：")
        for ct, p, s in sig_list:
            print(f"  {ct}: p={p:.4f} {s}")
    else:
        print("\n未发现原始 p < 0.05 的细胞类型。")

    # 所有图及表格输出
    plot_violin_significance(prop_df, signif_dict,
                            save_path=os.path.join(save_dir, 'FigA_all_violin.pdf'))
    plot_stacked_bar(prop_df,
                     save_path=os.path.join(save_dir, 'FigB_stacked_bar.pdf'))
    plot_heatmap_diff(prop_df,
                      save_path=os.path.join(save_dir, 'FigC_heatmap_diff.pdf'))
    if run_umap:
        plot_umap_facet(adata, celltype_col=celltype_col, cancer_col=cancer_col,
                        save_path=os.path.join(save_dir, 'FigD_UMAP_facet.pdf'))

    plot_significant_violin(prop_df, signif_dict,
                            save_path=os.path.join(save_dir, 'FigE_significant_violin.pdf'))
    plot_volcano_diff(prop_df, signif_dict,
                      save_path=os.path.join(save_dir, 'FigF_volcano_diff.pdf'))
    export_statistics_table(prop_df, signif_dict,
                            save_path=os.path.join(save_dir, 'celltype_proportion_statistics.csv'))

    print(f"所有结果已保存至: {save_dir}")

if __name__ == "__main__":
    main()
