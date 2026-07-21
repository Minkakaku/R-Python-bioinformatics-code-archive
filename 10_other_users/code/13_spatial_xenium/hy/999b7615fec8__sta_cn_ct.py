"""
Cell-type proportion comparison across CNs and groups
======================================================
- Input: adata_with_CNs_10.h5ad
- Columns used: sample, CNs_10, group_final, celltype
- Four contrasts:
    'PRE-R' vs 'PRE-NR'
    'PRE-R' vs 'POST-R'
    'PRE-NR' vs 'POST-NR'
    'POST-R' vs 'POST-NR'
- For each contrast and each CN, produce:
    * Individual boxplots for every cell type
    * A combined plot of cell types with significant differences (raw p < 0.05)
- Output: hierarchical folder structure with publication‑quality figures.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from typing import List, Tuple

# ================== 字体设置 ==================
from matplotlib import font_manager

font_manager.fontManager.addfont('/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf')

plt.rcParams['font.family'] = 'arial'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ============================================================
# 0. Load data and define parameters
# ============================================================
adata = sc.read('adata_with_CNs_10.h5ad')

CONTRASTS = [
    ('PRE-R', 'PRE-NR'),
    ('PRE-R', 'POST-R'),
    ('PRE-NR', 'POST-NR'),
    ('POST-R', 'POST-NR'),
]

# CNs sorted numerically
cn_series = adata.obs['CNs_10']
cn_list = sorted(cn_series.unique(), key=lambda x: int(x.replace('CN', '')))
celltype_order = sorted(adata.obs['celltype'].dropna().unique())

# ============================================================
# 1. Compute per‑sample, per‑CN cell‑type proportions
# ============================================================
def compute_proportions(adata, cn_list, celltype_order):
    rows = []
    for (sample, cn), sub in adata.obs.groupby(['sample', 'CNs_10']):
        counts = sub['celltype'].value_counts()
        total = counts.sum()
        for ct in celltype_order:
            prop = counts.get(ct, 0) / total if total > 0 else np.nan
            rows.append({
                'sample': sample,
                'group_final': sub['group_final'].iloc[0],
                'CN': cn,
                'celltype': ct,
                'proportion': prop,
            })
    return pd.DataFrame(rows)

prop_df = compute_proportions(adata, cn_list, celltype_order)

# ============================================================
# 2. Statistical testing (Mann–Whitney U) – NO FDR correction
# ============================================================
pvalue_records = []

for cn in cn_list:
    cn_data = prop_df[prop_df['CN'] == cn]
    for (g1, g2) in CONTRASTS:
        contrast_name = f"{g1}_vs_{g2}"
        for ct in celltype_order:
            vals_g1 = cn_data[(cn_data['group_final'] == g1) & (cn_data['celltype'] == ct)]['proportion'].dropna()
            vals_g2 = cn_data[(cn_data['group_final'] == g2) & (cn_data['celltype'] == ct)]['proportion'].dropna()
            if len(vals_g1) < 3 or len(vals_g2) < 3:
                pval = np.nan
            else:
                _, pval = mannwhitneyu(vals_g1, vals_g2, alternative='two-sided')
            pvalue_records.append((cn, contrast_name, ct, pval))

pval_df = pd.DataFrame(pvalue_records, columns=['CN', 'contrast', 'celltype', 'pvalue'])
pval_df['pvalue'] = pd.to_numeric(pval_df['pvalue'], errors='coerce')
# No FDR correction – keep raw p‑values

# Helper for stars
def pval_to_stars(p):
    if pd.isna(p):
        return ''
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'

# ============================================================
# 3. Plotting functions
# ============================================================
sns.set_style('ticks')
plt.rcParams['font.size'] = 7
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['ytick.major.width'] = 0.5

# Color palette for groups (colourblind‑friendly)
group_palette = {'PRE-R': '#D62728', 'PRE-NR': '#1F77B4',
                 'POST-R': '#D62728', 'POST-NR': '#1F77B4'}

def draw_single_celltype_boxplot(ax, data, g1, g2, ct, pval):
    """Draw a boxplot for a single cell type on given axes."""
    sub = data[data['celltype'] == ct]
    sns.boxplot(
        data=sub, x='group_final', y='proportion',
        order=[g1, g2], palette=group_palette,
        ax=ax, linewidth=0.5, fliersize=2, width=0.5,
        showmeans=False
    )
    # Overlay jittered points
    sns.stripplot(
        data=sub, x='group_final', y='proportion',
        order=[g1, g2], color='black', size=2, alpha=0.7,
        ax=ax, jitter=True
    )
    # Significance annotation
    stars = pval_to_stars(pval)
    if stars and stars != 'ns':
        y_max = sub['proportion'].max()
        y_pos = y_max * 1.15
        ax.text(0.5, y_pos, stars, ha='center', va='bottom', fontsize=9, fontweight='bold')
    # Styling
    ax.set_title(ct, fontsize=10, fontweight='bold')
    ax.set_ylabel('Proportion', fontsize=9)
    ax.set_xlabel('')
    ax.tick_params(labelsize=8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

def draw_combined_significant(ax, data, g1, g2, sig_cts, pval_dict):
    """Draw boxplots for significant cell types only."""
    sub = data[data['celltype'].isin(sig_cts)]
    if sub.empty:
        ax.text(0.5, 0.5, 'No significant cell types', ha='center', va='center',
                transform=ax.transAxes, fontsize=10)
        ax.axis('off')
        return
    sns.boxplot(
        data=sub, x='celltype', y='proportion', hue='group_final',
        order=sig_cts, hue_order=[g1, g2], palette=group_palette,
        ax=ax, linewidth=0.5, fliersize=2, width=0.6,
        showmeans=False
    )
    # Significance stars
    for i, ct in enumerate(sig_cts):
        p = pval_dict.get(ct, np.nan)
        stars = pval_to_stars(p)
        if stars and stars != 'ns':
            ct_data = sub[sub['celltype'] == ct]
            y_max = ct_data['proportion'].max() if len(ct_data) > 0 else 0
            y_pos = y_max * 1.1
            ax.text(i, y_pos, stars, ha='center', va='bottom', fontsize=8, fontweight='bold')
    ax.set_title(f'{cn} significant cell types', fontsize=10, fontweight='bold')
    ax.set_ylabel('Proportion', fontsize=9)
    ax.set_xlabel('')
    ax.tick_params(axis='x', rotation=45, labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.legend(frameon=False, fontsize=8, title='')
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

# ============================================================
# 4. Main generation loop
# ============================================================
base_out_dir = './CT_enrichment_comparison'
os.makedirs(base_out_dir, exist_ok=True)

for (g1, g2) in CONTRASTS:
    contrast_name = f"{g1}_vs_{g2}"
    contrast_dir = os.path.join(base_out_dir, contrast_name)
    os.makedirs(contrast_dir, exist_ok=True)
    
    for cn in cn_list:
        cn_dir = os.path.join(contrast_dir, cn)
        os.makedirs(cn_dir, exist_ok=True)
        
        # Subset data for this contrast and CN
        cn_data = prop_df[(prop_df['CN'] == cn) &
                          (prop_df['group_final'].isin([g1, g2]))].copy()
        
        # ---- Plot individual cell types ----
        for ct in celltype_order:
            fig, ax = plt.subplots(figsize=(3, 3))
            # Get raw p‑value for this combination
            pval_row = pval_df.loc[
                (pval_df['CN'] == cn) &
                (pval_df['contrast'] == contrast_name) &
                (pval_df['celltype'] == ct), 'pvalue'
            ]
            pval = pval_row.values[0] if len(pval_row) > 0 else np.nan
            draw_single_celltype_boxplot(ax, cn_data, g1, g2, ct, pval)
            fig.tight_layout()
            fig.savefig(os.path.join(cn_dir, f'{ct}.pdf'), dpi=300, bbox_inches='tight')
            fig.savefig(os.path.join(cn_dir, f'{ct}.png'), dpi=300, bbox_inches='tight')
            plt.close(fig)
        
        # ---- Plot combined significant cell types (raw p < 0.05) ----
        sig_mask = (pval_df['CN'] == cn) & (pval_df['contrast'] == contrast_name)
        sig_sub = pval_df.loc[sig_mask & (pval_df['pvalue'] < 0.05)]
        sig_cts = sig_sub['celltype'].tolist()
        pval_dict = dict(zip(sig_sub['celltype'], sig_sub['pvalue']))
        
        fig, ax = plt.subplots(figsize=(max(4, len(sig_cts)*1.2), 3.5))
        draw_combined_significant(ax, cn_data, g1, g2, sig_cts, pval_dict)
        fig.tight_layout()
        fig.savefig(os.path.join(cn_dir, 'significant_combined.pdf'), dpi=300, bbox_inches='tight')
        fig.savefig(os.path.join(cn_dir, 'significant_combined.png'), dpi=300, bbox_inches='tight')
        plt.close(fig)

print(f"All figures saved under: {base_out_dir}")
