#!/usr/bin/env python3
"""
CN Proportion AUROC Analysis (Python version)
-----------------------------------------------
Calculate per-sample cell proportion in each CN, then compare A vs N and C vs N
using directional AUROC and Wilcoxon test.
Generate heatmaps for CNs_10 and CNs_15 with stars for AUROC > 0.75.
"""

import os, re, sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from matplotlib.colors import LinearSegmentedColormap

# ========== Configuration ==========
H5AD_PATH = "All_filter_state_samp_new_with_CNs10-15.h5ad"
OUT_DIR = "./auroc_results"
os.makedirs(OUT_DIR, exist_ok=True)

CN_COLS = ["CNs_10", "CNs_15"]
SAMPLE_COL = "sample_new"
STATE_COL = "state"
REF_GROUP = "N"                # baseline
POS_GROUPS = ["A", "C"]        # groups to compare against baseline
AUC_STAR_CUTOFF = 0.75
P_ADJUST_METHOD = "holm"       # multiple testing correction across CNs
FIG_DPI = 200
# ===================================

def sort_cn_labels(labels):
    """Sort CN labels by trailing numeric part (e.g. CN0 < CN1 < CN10)."""
    def key(x):
        s = str(x)
        m = re.search(r'(\d+)$', s)
        if m:
            return (0, int(m.group(1)))
        return (1, s)
    return sorted(labels, key=key)

def compute_sample_cn_fraction(obs_df, cn_col):
    """
    For each sample, compute the fraction of cells in each CN.
    Returns a DataFrame (samples x CNs) with proportions.
    """
    sample_totals = obs_df.groupby(SAMPLE_COL).size()
    cn_counts = obs_df.groupby([SAMPLE_COL, cn_col]).size().unstack(fill_value=0)
    cn_fraction = cn_counts.div(sample_totals, axis=0)
    # ensure all CNs are present (even if some samples miss a CN)
    all_cns = sort_cn_labels(obs_df[cn_col].dropna().unique())
    cn_fraction = cn_fraction.reindex(columns=all_cns, fill_value=0.0)
    return cn_fraction

def directional_auc(pos_vals, neg_vals):
    """
    Compute directional AUROC (Mann-Whitney U / (n1*n2)).
    >0.5 means pos group tends to have larger values.
    """
    pos = np.asarray(pos_vals, dtype=float)
    neg = np.asarray(neg_vals, dtype=float)
    pos = pos[np.isfinite(pos)]
    neg = neg[np.isfinite(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    # Combine and rank
    all_vals = np.concatenate([pos, neg])
    ranks = np.argsort(np.argsort(all_vals)) + 1  # simple ranking (no ties handling? better use scipy)
    # More accurate: use scipy.stats.rankdata with average ties
    from scipy.stats import rankdata
    ranks = rankdata(all_vals, method='average')
    U = np.sum(ranks[:len(pos)]) - len(pos) * (len(pos) + 1) / 2
    return U / (len(pos) * len(neg))

def compute_auroc_table(cn_fraction, state_series, pos_group, ref_group=REF_GROUP):
    """
    For a single positive group (e.g. 'A') vs reference ('N'), compute per-CN:
    - n_pos, n_ref
    - mean_pos, mean_ref
    - median_pos, median_ref
    - AUROC
    - p_value (Wilcoxon)
    Returns a DataFrame.
    """
    # restrict to needed states
    needed_states = [ref_group, pos_group]
    mask = state_series.isin(needed_states)
    state_sub = state_series.loc[mask]
    frac_sub = cn_fraction.loc[mask.index.intersection(cn_fraction.index)]

    results = []
    for cn in frac_sub.columns:
        pos_samples = state_sub[state_sub == pos_group].index
        ref_samples = state_sub[state_sub == ref_group].index
        pos_vals = frac_sub.loc[pos_samples.intersection(frac_sub.index), cn].values
        ref_vals = frac_sub.loc[ref_samples.intersection(frac_sub.index), cn].values

        pos_vals = pos_vals[np.isfinite(pos_vals)]
        ref_vals = ref_vals[np.isfinite(ref_vals)]

        n_pos = len(pos_vals)
        n_ref = len(ref_vals)
        if n_pos < 2 or n_ref < 2:
            p_val = np.nan
        else:
            try:
                _, p_val = mannwhitneyu(pos_vals, ref_vals, alternative='two-sided')
            except ValueError:
                p_val = np.nan

        auroc = directional_auc(pos_vals, ref_vals)

        results.append({
            'CN': cn,
            'positive_group': pos_group,
            'reference_group': ref_group,
            'n_positive': n_pos,
            'n_reference': n_ref,
            'mean_positive': np.mean(pos_vals) if n_pos > 0 else np.nan,
            'mean_reference': np.mean(ref_vals) if n_ref > 0 else np.nan,
            'median_positive': np.median(pos_vals) if n_pos > 0 else np.nan,
            'median_reference': np.median(ref_vals) if n_ref > 0 else np.nan,
            'auroc': auroc,
            'p_value': p_val
        })

    df = pd.DataFrame(results)
    # Multiple testing correction across CNs (within this comparison)
    if len(df) > 1:
        from statsmodels.stats.multitest import multipletests
        pvals = df['p_value'].values
        mask = ~np.isnan(pvals)
        if np.any(mask):
            corrected = np.full_like(pvals, np.nan, dtype=float)
            _, corrected_vals, _, _ = multipletests(pvals[mask], method=P_ADJUST_METHOD)
            corrected[mask] = corrected_vals
            df['p_adj'] = corrected
        else:
            df['p_adj'] = np.nan
    else:
        df['p_adj'] = df['p_value']
    return df

def draw_auroc_heatmap(auroc_df, cn_col, pos_group, out_prefix):
    """
    Draw a heatmap with one row (single comparison) and CN columns.
    AUROC is the fill value; stars added where AUROC > AUC_STAR_CUTOFF.
    """
    # Prepare matrix
    cn_order = sort_cn_labels(auroc_df['CN'].unique())
    auroc_series = auroc_df.set_index('CN')['auroc']
    auroc_series = auroc_series.reindex(cn_order)
    mat = auroc_series.values.reshape(1, -1)
    star_mat = np.where(np.isfinite(mat) & (mat > AUC_STAR_CUTOFF), "*", "").reshape(1, -1)

    # Colormap: white to red (or blue-white-red for AUC 0-1)
    # Use white (0.5) to red (1) since we expect positive direction.
    # For consistency with original R script we can use a diverging colormap with blue for low and red for high.
    # But user previously requested "white to red". Here we apply a red gradient for values >0.5,
    # and blue for <0.5? Since we only care about >0.5, white-to-red is okay.
    cmap = LinearSegmentedColormap.from_list('white_red', ['white', 'darkred'])
    vmin, vmax = 0.0, 1.0

    fig, ax = plt.subplots(figsize=(max(6, len(cn_order)*0.5), 1.5))
    sns.heatmap(mat, annot=star_mat, fmt='', cmap=cmap,
                vmin=vmin, vmax=vmax, linewidths=0.5, ax=ax,
                cbar_kws={"label": "AUROC"},
                xticklabels=cn_order, yticklabels=[f"{pos_group} vs {REF_GROUP}"],
                annot_kws={'size': 16, 'weight': 'bold'})
    ax.set_title(f"AUROC: {pos_group} vs {REF_GROUP} ({cn_col})")
    ax.set_ylabel("")
    plt.tight_layout()
    plt.savefig(out_prefix + ".png", dpi=FIG_DPI, bbox_inches='tight')
    plt.savefig(out_prefix + ".pdf", bbox_inches='tight')
    plt.close()
    print(f"Saved heatmap: {out_prefix}.png/.pdf")

def main():
    print(f"Loading {H5AD_PATH}...")
    adata = sc.read_h5ad(H5AD_PATH)

    # Check required columns
    required = [SAMPLE_COL, STATE_COL] + CN_COLS
    for col in required:
        assert col in adata.obs.columns, f"Missing column: {col}"

    obs_df = adata.obs[[SAMPLE_COL, STATE_COL] + CN_COLS].dropna()
    print(f"Total cells: {len(obs_df)}")

    state_series = obs_df.groupby(SAMPLE_COL)[STATE_COL].first()  # one state per sample

    for cn_col in CN_COLS:
        print(f"\n=== Processing {cn_col} ===")
        cn_frac = compute_sample_cn_fraction(obs_df, cn_col)
        # Align samples
        common_samples = cn_frac.index.intersection(state_series.index)
        cn_frac = cn_frac.loc[common_samples]
        state_series = state_series.loc[common_samples]

        for pos_group in POS_GROUPS:
            print(f"  Comparison: {pos_group} vs {REF_GROUP}")
            auroc_df = compute_auroc_table(cn_frac, state_series, pos_group, REF_GROUP)

            # Save CSV
            csv_file = os.path.join(OUT_DIR, f"{cn_col}_{pos_group}_vs_{REF_GROUP}_AUROC.csv")
            auroc_df.to_csv(csv_file, index=False)
            print(f"  Saved AUROC table to {csv_file}")

            # Heatmap
            heatmap_prefix = os.path.join(OUT_DIR, f"{cn_col}_{pos_group}_vs_{REF_GROUP}_AUROC_heatmap")
            draw_auroc_heatmap(auroc_df, cn_col, pos_group, heatmap_prefix)

    print("\nAll done. Results are in:", os.path.abspath(OUT_DIR))

if __name__ == "__main__":
    main()
