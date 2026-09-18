#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integrated analysis of cell type differences between inner and outer mucus zones
Statistical analysis across all samples
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import warnings
from pathlib import Path
from scipy import stats
from scipy.stats import mannwhitneyu, kruskal, ranksums, ttest_ind, ttest_rel, wilcoxon, shapiro, levene
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
warnings.filterwarnings('ignore')

def load_and_prepare_data(h5ad_path):
    """
    Load h5ad data and prepare for analysis
    
    Parameters:
    h5ad_path: h5ad file path
    
    Returns:
    adata: AnnData object
    """
    print(f"Loading h5ad file: {h5ad_path}")
    
    try:
        adata = sc.read_h5ad(h5ad_path)
    except Exception as e:
        print(f"Failed to read file: {e}")
        return None
    
    # Check required columns
    required_cols = ['mucus_zone', 'sample', 'annotation']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}")
        print(f"Available columns: {adata.obs.columns.tolist()}")
        return None
    
    print(f"Data shape: {adata.shape}")
    print(f"Number of cells: {adata.n_obs}")
    print(f"Number of genes: {adata.n_vars}")
    
    # Print basic information
    print("\n=== Basic Information ===")
    print(f"Number of groups: {adata.obs['mucus_zone'].nunique()}")
    print(f"Number of samples: {adata.obs['sample'].nunique()}")
    print(f"Number of cell types: {adata.obs['annotation'].nunique()}")
    
    # Print group information
    print(f"\nGroup distribution:")
    group_counts = adata.obs['mucus_zone'].value_counts()
    for group, count in group_counts.items():
        percentage = count / adata.n_obs * 100
        print(f"  {group}: {count} cells ({percentage:.2f}%)")
    
    return adata

def calculate_group_level_proportions(adata, min_cells_per_sample_zone=10):
    """
    Calculate cell type proportions at group level (inner vs outer)
    
    Parameters:
    adata: AnnData object
    min_cells_per_sample_zone: minimum cells per sample per zone
    
    Returns:
    group_prop_df: DataFrame with group-level proportions
    sample_zone_prop_df: DataFrame with sample-zone level proportions
    """
    print(f"\n=== Calculating group-level proportions ===")
    
    # Get all samples and cell types
    all_samples = adata.obs['sample'].unique()
    all_celltypes = adata.obs['annotation'].unique()
    
    print(f"Number of samples: {len(all_samples)}")
    print(f"Number of cell types: {len(all_celltypes)}")
    
    # Store proportions for each sample-zone combination
    sample_zone_data = []
    
    for sample in all_samples:
        sample_mask = adata.obs['sample'] == sample
        sample_data = adata[sample_mask]
        
        for zone in ['inner', 'outer']:
            zone_mask = sample_data.obs['mucus_zone'] == zone
            zone_data = sample_data[zone_mask]
            
            if len(zone_data) < min_cells_per_sample_zone:
                print(f"  Sample {sample}, zone {zone}: {len(zone_data)} cells (< {min_cells_per_sample_zone}), skipping")
                continue
            
            # Calculate proportions
            celltype_counts = zone_data.obs['annotation'].value_counts()
            total_cells = len(zone_data)
            
            row_data = {
                'sample': sample,
                'mucus_zone': zone,
                'total_cells': total_cells
            }
            
            for celltype in all_celltypes:
                count = celltype_counts.get(celltype, 0)
                proportion = count / total_cells * 100
                row_data[celltype] = proportion
            
            sample_zone_data.append(row_data)
    
    # Create sample-zone level DataFrame
    sample_zone_prop_df = pd.DataFrame(sample_zone_data)
    
    if sample_zone_prop_df.empty:
        print("Error: No valid sample-zone combinations")
        return None, None
    
    # Create group-level proportions (mean across samples)
    group_prop_data = []
    
    for zone in sample_zone_prop_df['mucus_zone'].unique():
        zone_data = sample_zone_prop_df[sample_zone_prop_df['mucus_zone'] == zone]
        
        row_data = {'mucus_zone': zone, 'n_samples': len(zone_data)}
        
        for celltype in all_celltypes:
            if celltype in zone_data.columns:
                mean_prop = zone_data[celltype].mean()
                std_prop = zone_data[celltype].std()
                sem_prop = zone_data[celltype].sem()
                
                row_data[f'{celltype}_mean'] = mean_prop
                row_data[f'{celltype}_std'] = std_prop
                row_data[f'{celltype}_sem'] = sem_prop
                row_data[f'{celltype}_n'] = len(zone_data[celltype].dropna())
        
        group_prop_data.append(row_data)
    
    group_prop_df = pd.DataFrame(group_prop_data)
    
    print(f"Sample-zone combinations analyzed: {len(sample_zone_prop_df)}")
    print(f"Inner samples: {len(sample_zone_prop_df[sample_zone_prop_df['mucus_zone'] == 'inner'])}")
    print(f"Outer samples: {len(sample_zone_prop_df[sample_zone_prop_df['mucus_zone'] == 'outer'])}")
    
    return group_prop_df, sample_zone_prop_df

def perform_group_comparison_tests(sample_zone_prop_df, celltype, test_type='mannwhitney'):
    """
    Perform statistical tests comparing inner vs outer for a cell type
    
    Parameters:
    sample_zone_prop_df: DataFrame with sample-zone level proportions
    celltype: cell type to analyze
    test_type: statistical test to use
    
    Returns:
    Dictionary with test results
    """
    # Extract data for the cell type
    inner_data = sample_zone_prop_df[
        (sample_zone_prop_df['mucus_zone'] == 'inner') & 
        (sample_zone_prop_df[celltype].notna())
    ][celltype].values
    
    outer_data = sample_zone_prop_df[
        (sample_zone_prop_df['mucus_zone'] == 'outer') & 
        (sample_zone_prop_df[celltype].notna())
    ][celltype].values
    
    if len(inner_data) < 3 or len(outer_data) < 3:
        return {
            'annotation': celltype,
            'n_inner': len(inner_data),
            'n_outer': len(outer_data),
            'test_performed': False,
            'reason': f'Insufficient data: inner={len(inner_data)}, outer={len(outer_data)}'
        }
    
    # Calculate basic statistics
    inner_mean = np.mean(inner_data)
    inner_median = np.median(inner_data)
    inner_std = np.std(inner_data)
    inner_sem = inner_std / np.sqrt(len(inner_data))
    
    outer_mean = np.mean(outer_data)
    outer_median = np.median(outer_data)
    outer_std = np.std(outer_data)
    outer_sem = outer_std / np.sqrt(len(outer_data))
    
    # Calculate fold change
    fold_change = inner_mean / outer_mean if outer_mean > 0 else np.inf
    log2_fold_change = np.log2(fold_change) if fold_change != np.inf else np.inf
    
    # Perform statistical test
    p_value = np.nan
    test_name = ""
    statistic = np.nan
    
    try:
        if test_type == 'mannwhitney':
            # Mann-Whitney U test (non-parametric, independent samples)
            statistic, p_value = mannwhitneyu(inner_data, outer_data, alternative='two-sided')
            test_name = "Mann-Whitney U test"
        
        elif test_type == 'ttest_ind':
            # Independent t-test (parametric)
            # Check normality assumption
            _, p_norm_inner = shapiro(inner_data) if len(inner_data) < 5000 else (np.nan, 0.05)
            _, p_norm_outer = shapiro(outer_data) if len(outer_data) < 5000 else (np.nan, 0.05)
            
            # Check equal variance
            _, p_var = levene(inner_data, outer_data)
            
            if p_norm_inner > 0.05 and p_norm_outer > 0.05:
                if p_var > 0.05:
                    # Equal variance
                    statistic, p_value = ttest_ind(inner_data, outer_data, equal_var=True)
                    test_name = "Independent t-test (equal variance)"
                else:
                    # Unequal variance
                    statistic, p_value = ttest_ind(inner_data, outer_data, equal_var=False)
                    test_name = "Welch's t-test (unequal variance)"
            else:
                # Not normal, fall back to Mann-Whitney
                statistic, p_value = mannwhitneyu(inner_data, outer_data, alternative='two-sided')
                test_name = "Mann-Whitney U test (non-normal data)"
        
        elif test_type == 'ranksums':
            # Wilcoxon rank-sum test
            statistic, p_value = ranksums(inner_data, outer_data)
            test_name = "Wilcoxon rank-sum test"
        
        else:
            raise ValueError(f"Unknown test type: {test_type}")
    
    except Exception as e:
        print(f"  Error performing test for {celltype}: {e}")
        return {
            'annotation': celltype,
            'n_inner': len(inner_data),
            'n_outer': len(outer_data),
            'test_performed': False,
            'reason': f'Test error: {str(e)}'
        }
    
    # Calculate effect size
    if test_name.startswith("Mann-Whitney") or test_name.startswith("Wilcoxon"):
        # Common language effect size
        n1, n2 = len(inner_data), len(outer_data)
        u_stat = statistic if test_name.startswith("Mann-Whitney") else None
        
        if u_stat is not None:
            # For Mann-Whitney U
            effect_size = u_stat / (n1 * n2)  # Proportion of pairs where inner > outer
        else:
            effect_size = np.nan
    
    else:
        # Cohen's d for t-tests
        pooled_std = np.sqrt(((len(inner_data)-1)*inner_std**2 + (len(outer_data)-1)*outer_std**2) / 
                           (len(inner_data) + len(outer_data) - 2))
        effect_size = (inner_mean - outer_mean) / pooled_std if pooled_std > 0 else np.nan
    
    return {
        'annotation': celltype,
        'n_inner': len(inner_data),
        'n_outer': len(outer_data),
        'test_performed': True,
        'test_name': test_name,
        'statistic': statistic,
        'p_value': p_value,
        'p_value_adj': np.nan,  # Will be filled after multiple testing correction
        'significant': p_value < 0.05,
        'inner_mean': inner_mean,
        'inner_median': inner_median,
        'inner_std': inner_std,
        'inner_sem': inner_sem,
        'outer_mean': outer_mean,
        'outer_median': outer_median,
        'outer_std': outer_std,
        'outer_sem': outer_sem,
        'fold_change': fold_change,
        'log2_fold_change': log2_fold_change,
        'effect_size': effect_size,
        'enrichment': 'inner' if inner_mean > outer_mean else 'outer'
    }

def perform_paired_analysis(sample_zone_prop_df, celltype):
    """
    Perform paired analysis for samples that have both inner and outer data
    
    Parameters:
    sample_zone_prop_df: DataFrame with sample-zone level proportions
    celltype: cell type to analyze
    
    Returns:
    Dictionary with paired test results
    """
    # Create pivot table: sample vs zone
    pivot_df = sample_zone_prop_df.pivot_table(
        index='sample', 
        columns='mucus_zone', 
        values=celltype
    )
    
    # Keep only samples with both inner and outer data
    paired_df = pivot_df.dropna()
    
    if len(paired_df) < 3:
        return {
            'annotation': celltype,
            'n_pairs': len(paired_df),
            'test_performed': False,
            'reason': f'Insufficient paired samples: {len(paired_df)}'
        }
    
    inner_paired = paired_df['inner'].values
    outer_paired = paired_df['outer'].values
    
    # Calculate differences
    differences = inner_paired - outer_paired
    mean_diff = np.mean(differences)
    std_diff = np.std(differences)
    
    # Paired t-test
    t_stat, p_value_ttest = ttest_rel(inner_paired, outer_paired)
    
    # Wilcoxon signed-rank test
    try:
        t_stat_wilcoxon, p_value_wilcoxon = wilcoxon(inner_paired, outer_paired)
    except:
        p_value_wilcoxon = np.nan
    
    # Calculate paired effect size
    cohen_d = mean_diff / np.std(differences) if np.std(differences) > 0 else np.nan
    
    return {
        'annotation': celltype,
        'n_pairs': len(paired_df),
        'test_performed': True,
        'test_name': 'Paired analysis',
        'p_value_ttest': p_value_ttest,
        'p_value_wilcoxon': p_value_wilcoxon,
        'significant_ttest': p_value_ttest < 0.05,
        'significant_wilcoxon': p_value_wilcoxon < 0.05 if not np.isnan(p_value_wilcoxon) else False,
        'mean_difference': mean_diff,
        'std_difference': std_diff,
        'cohen_d': cohen_d,
        'inner_mean_paired': np.mean(inner_paired),
        'outer_mean_paired': np.mean(outer_paired)
    }

def analyze_all_samples_comparison(adata, output_dir="integrated_group_comparison", 
                                   min_mean_proportion=0.1, test_type='mannwhitney',
                                   perform_paired_analysis_flag=True):
    """
    Integrated analysis comparing inner vs outer across all samples
    
    Parameters:
    adata: AnnData object
    output_dir: output directory
    min_mean_proportion: minimum mean proportion threshold
    test_type: statistical test to use
    perform_paired_analysis_flag: whether to perform paired analysis
    """
    print(f"\n=== Integrated analysis: Inner vs Outer across all samples ===")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Calculate proportions
    group_prop_df, sample_zone_prop_df = calculate_group_level_proportions(adata)
    
    if group_prop_df is None or sample_zone_prop_df is None:
        print("Error: Could not calculate proportions")
        return None
    
    # Save proportion data
    group_prop_path = output_path / "group_level_proportions.csv"
    group_prop_df.to_csv(group_prop_path, encoding='utf-8', index=False)
    
    sample_zone_path = output_path / "sample_zone_proportions.csv"
    sample_zone_prop_df.to_csv(sample_zone_path, encoding='utf-8', index=False)
    
    print(f"\nProportion data saved:")
    print(f"  Group-level: {group_prop_path}")
    print(f"  Sample-zone level: {sample_zone_path}")
    
    # Get cell types (exclude metadata columns)
    metadata_cols = ['sample', 'mucus_zone', 'total_cells']
    all_celltypes = [col for col in sample_zone_prop_df.columns if col not in metadata_cols]
    
    # Filter cell types by mean proportion
    filtered_celltypes = []
    
    for celltype in all_celltypes:
        # Calculate overall mean proportion
        overall_mean = sample_zone_prop_df[celltype].mean()
        
        if overall_mean >= min_mean_proportion:
            filtered_celltypes.append(celltype)
        else:
            print(f"  Filtered out {celltype}: mean proportion {overall_mean:.2f}% < {min_mean_proportion}%")
    
    print(f"\nAnalyzing {len(filtered_celltypes)} cell types "
          f"(filtered out {len(all_celltypes) - len(filtered_celltypes)} rare types)")
    
    # Perform statistical tests
    print(f"\nPerforming statistical tests ({test_type})...")
    
    test_results = []
    paired_results = []
    
    for i, celltype in enumerate(filtered_celltypes, 1):
        print(f"  {i}/{len(filtered_celltypes)}: {celltype}")
        
        # Unpaired comparison
        result = perform_group_comparison_tests(sample_zone_prop_df, celltype, test_type)
        test_results.append(result)
        
        # Paired comparison
        if perform_paired_analysis_flag:
            paired_result = perform_paired_analysis(sample_zone_prop_df, celltype)
            paired_results.append(paired_result)
    
    # Convert to DataFrames
    test_df = pd.DataFrame(test_results)
    paired_df = pd.DataFrame(paired_results) if paired_results else pd.DataFrame()
    
    # Apply multiple testing correction
    if not test_df.empty and 'p_value' in test_df.columns:
        valid_pvals = test_df['p_value'].dropna()
        if len(valid_pvals) > 0:
            _, pvals_adj, _, _ = multipletests(valid_pvals, method='fdr_bh')
            test_df.loc[valid_pvals.index, 'p_value_adj'] = pvals_adj
            test_df['significant_adj'] = test_df['p_value_adj'] < 0.05
    
    # Save results
    test_results_path = output_path / "statistical_test_results.csv"
    test_df.to_csv(test_results_path, encoding='utf-8', index=False)
    print(f"Statistical test results saved: {test_results_path}")
    
    if not paired_df.empty:
        paired_results_path = output_path / "paired_analysis_results.csv"
        paired_df.to_csv(paired_results_path, encoding='utf-8', index=False)
        print(f"Paired analysis results saved: {paired_results_path}")
    
    # Create visualizations
    print(f"\nCreating visualizations...")
    
    # 1. Volcano plot
    create_volcano_plot(test_df, output_path)
    
    # 2. Bar plot of significant cell types
    create_significant_barplot(test_df, output_path)
    
    # 3. Heatmap of proportions
    create_proportion_heatmap(test_df, sample_zone_prop_df, output_path)
    
    # 4. Box plots for top significant cell types
    create_top_celltype_boxplots(test_df, sample_zone_prop_df, output_path)
    
    # 5. Correlation between inner and outer
    create_correlation_plot(sample_zone_prop_df, test_df, output_path)
    
    # 6. Paired comparison plots
    if not paired_df.empty:
        create_paired_comparison_plots(paired_df, sample_zone_prop_df, output_path)
    
    # Create summary report
    create_summary_report(test_df, paired_df, sample_zone_prop_df, output_path)
    
    return test_df, paired_df, sample_zone_prop_df

def create_volcano_plot(test_df, output_dir, pval_threshold=0.05, 
                       fc_threshold=1.5, top_n=20, figsize=(12, 8)):
    """
    Create volcano plot showing -log10(p) vs log2(fold change)
    """
    if test_df.empty or 'log2_fold_change' not in test_df.columns or 'p_value' not in test_df.columns:
        print("  Cannot create volcano plot: missing data")
        return
    
    # Prepare data
    plot_df = test_df.copy()
    plot_df = plot_df.dropna(subset=['log2_fold_change', 'p_value'])
    
    if plot_df.empty:
        print("  Cannot create volcano plot: no valid data points")
        return
    
    # Calculate -log10(p)
    plot_df['neg_log10_p'] = -np.log10(plot_df['p_value'])
    
    # Mark significance
    plot_df['significant'] = (plot_df['p_value'] < pval_threshold) & \
                            (abs(plot_df['log2_fold_change']) > np.log2(fc_threshold))
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create scatter plot
    scatter = ax.scatter(plot_df['log2_fold_change'], plot_df['neg_log10_p'],
                        c=plot_df['significant'].map({True: 'red', False: 'gray'}),
                        alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
    
    # Add threshold lines
    ax.axhline(y=-np.log10(pval_threshold), color='red', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(x=np.log2(fc_threshold), color='blue', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(x=-np.log2(fc_threshold), color='blue', linestyle='--', alpha=0.5, linewidth=1)
    
    # Label top N significant points
    significant_df = plot_df[plot_df['significant']].copy()
    if len(significant_df) > 0:
        # Sort by significance
        significant_df = significant_df.sort_values('neg_log10_p', ascending=False)
        
        for i, (idx, row) in enumerate(significant_df.head(top_n).iterrows()):
            ax.annotate(row['annotation'], 
                       (row['log2_fold_change'], row['neg_log10_p']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.8,
                       bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.5))
    
    # Set labels and title
    ax.set_xlabel('log2(Fold Change) [Inner/Outer]', fontsize=12)
    ax.set_ylabel('-log10(p-value)', fontsize=12)
    ax.set_title('Volcano Plot: Inner vs Outer Cell Type Differences\n' +
                f'({len(plot_df)} cell types, {len(significant_df)} significant)', 
                fontsize=14, fontweight='bold')
    
    ax.grid(True, alpha=0.3)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', alpha=0.6, label=f'Significant (p<{pval_threshold}, |FC|>{fc_threshold})'),
        Patch(facecolor='gray', alpha=0.6, label='Not significant')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    
    # Save figure
    png_path = output_dir / "volcano_plot.png"
    pdf_path = output_dir / "volcano_plot.pdf"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  Volcano plot saved: {png_path}")

def create_significant_barplot(test_df, output_dir, top_n=25, figsize=(14, 8)):
    """
    Create bar plot showing top significant cell types
    """
    if test_df.empty or 'significant' not in test_df.columns:
        print("  Cannot create significant bar plot: missing data")
        return
    
    # Filter significant results
    significant_df = test_df[test_df['significant'] == True].copy()
    
    if significant_df.empty:
        print("  No significant results to plot")
        return
    
    # Sort by fold change
    significant_df = significant_df.sort_values('fold_change', ascending=False)
    
    # Take top N
    plot_df = significant_df.head(min(top_n, len(significant_df)))
    
    # Create figure
    fig, axes = plt.subplots(2, 1, figsize=figsize, gridspec_kw={'height_ratios': [3, 1]})
    
    # Top plot: Mean proportions
    ax1 = axes[0]
    
    x = np.arange(len(plot_df))
    width = 0.35
    
    # Get data
    inner_means = plot_df['inner_mean'].values
    outer_means = plot_df['outer_mean'].values
    inner_sem = plot_df['inner_sem'].values
    outer_sem = plot_df['outer_sem'].values
    
    # Create bars with error bars
    rects1 = ax1.bar(x - width/2, inner_means, width, yerr=inner_sem,
                     label='Inner', color='lightblue', alpha=0.8,
                     capsize=3, error_kw={'elinewidth': 1, 'capthick': 1})
    
    rects2 = ax1.bar(x + width/2, outer_means, width, yerr=outer_sem,
                     label='Outer', color='lightcoral', alpha=0.8,
                     capsize=3, error_kw={'elinewidth': 1, 'capthick': 1})
    
    # Add p-value annotations
    for i, (idx, row) in enumerate(plot_df.iterrows()):
        p_val = row['p_value']
        
        # Format p-value
        if p_val < 0.0001:
            p_text = 'p < 0.0001'
        elif p_val < 0.001:
            p_text = f'p = {p_val:.4f}'
        elif p_val < 0.01:
            p_text = f'p = {p_val:.3f}'
        else:
            p_text = f'p = {p_val:.2f}'
        
        # Add star for significance
        star_count = 0
        if p_val < 0.001:
            stars = '***'
        elif p_val < 0.01:
            stars = '**'
        elif p_val < 0.05:
            stars = '*'
        else:
            stars = ''
        
        # Position annotation
        max_val = max(row['inner_mean'], row['outer_mean'])
        ax1.text(i, max_val + 1, stars, ha='center', va='bottom', 
                fontsize=10, color='red', fontweight='bold')
        
        # Add p-value text below
        ax1.text(i, -1, p_text, ha='center', va='top', 
                fontsize=7, rotation=45, alpha=0.7)
    
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Mean Proportion (%)')
    ax1.set_title(f'Top {len(plot_df)} Significantly Different Cell Types\n' +
                  'Mean ± SEM', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(plot_df['annotation'], rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Bottom plot: Fold changes
    ax2 = axes[1]
    
    # Get fold changes
    fold_changes = plot_df['fold_change'].values
    log2_fc = plot_df['log2_fold_change'].values
    
    # Color by direction
    colors = ['green' if fc > 1 else 'red' for fc in fold_changes]
    
    # Create bar plot
    rects3 = ax2.bar(x, log2_fc, color=colors, alpha=0.7)
    
    # Add horizontal lines
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.5)
    ax2.axhline(y=np.log2(1.5), color='blue', linestyle='--', linewidth=0.5, alpha=0.5)
    ax2.axhline(y=-np.log2(1.5), color='blue', linestyle='--', linewidth=0.5, alpha=0.5)
    
    ax2.set_xlabel('Cell Type')
    ax2.set_ylabel('log2(Fold Change)')
    ax2.set_title('Fold Change (Inner/Outer)', fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels(plot_df['annotation'], rotation=45, ha='right')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add fold change values
    for rect, fc, logfc in zip(rects3, fold_changes, log2_fc):
        height = rect.get_height()
        if abs(height) > 0.1:
            ax2.text(rect.get_x() + rect.get_width()/2., 
                    height + (0.1 if height > 0 else -0.3),
                    f'{fc:.2f}x', ha='center', va='center' if height > 0 else 'top', 
                    fontsize=7, fontweight='bold')
    
    plt.tight_layout()
    
    # Save figure
    png_path = output_dir / "significant_celltypes_barplot.png"
    pdf_path = output_dir / "significant_celltypes_barplot.pdf"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  Significant cell types bar plot saved: {png_path}")

def create_proportion_heatmap(test_df, sample_zone_prop_df, output_dir, 
                            top_n=30, figsize=(16, 12)):
    """
    Create heatmap of cell type proportions
    """
    if test_df.empty or sample_zone_prop_df.empty:
        print("  Cannot create proportion heatmap: missing data")
        return
    
    # Get top significant cell types
    significant_df = test_df[test_df['significant'] == True].copy()
    
    if significant_df.empty:
        print("  No significant results for heatmap")
        return
    
    # Sort by fold change
    significant_df = significant_df.sort_values('fold_change', ascending=False)
    
    # Take top N
    top_celltypes = significant_df.head(min(top_n, len(significant_df)))['annotation'].tolist()
    
    if not top_celltypes:
        return
    
    # Prepare data for heatmap
    heatmap_data = []
    sample_order = []
    
    # Separate inner and outer
    for zone in ['inner', 'outer']:
        zone_data = sample_zone_prop_df[sample_zone_prop_df['mucus_zone'] == zone]
        
        for _, sample_row in zone_data.iterrows():
            sample_name = sample_row['sample']
            row_data = []
            
            for celltype in top_celltypes:
                if celltype in sample_row:
                    row_data.append(sample_row[celltype])
                else:
                    row_data.append(0)
            
            heatmap_data.append(row_data)
            sample_order.append(f"{sample_name} ({zone})")
    
    heatmap_array = np.array(heatmap_data)
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=figsize, gridspec_kw={'width_ratios': [4, 1]})
    
    # Main heatmap
    ax1 = axes[0]
    
    # Normalize for better visualization
    heatmap_norm = (heatmap_array - heatmap_array.min()) / (heatmap_array.max() - heatmap_array.min() + 1e-10)
    
    im = ax1.imshow(heatmap_norm, cmap='YlOrRd', aspect='auto')
    
    # Set ticks
    ax1.set_yticks(range(len(sample_order)))
    ax1.set_yticklabels(sample_order, fontsize=8)
    ax1.set_xticks(range(len(top_celltypes)))
    ax1.set_xticklabels(top_celltypes, rotation=45, ha='right', fontsize=9)
    
    ax1.set_title('Cell Type Proportions Across Samples\n' +
                 f'Top {len(top_celltypes)} significant cell types', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label('Normalized Proportion', fontsize=12)
    
    # Add zone separator
    n_inner = len(sample_zone_prop_df[sample_zone_prop_df['mucus_zone'] == 'inner'])
    if n_inner > 0 and n_inner < len(sample_order):
        ax1.axhline(y=n_inner - 0.5, color='black', linestyle='-', linewidth=2)
        ax1.text(0, n_inner/2, 'Inner', ha='center', va='center', 
                transform=ax1.get_yaxis_transform(), fontsize=10, 
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
        ax1.text(0, n_inner + (len(sample_order)-n_inner)/2, 'Outer', ha='center', va='center', 
                transform=ax1.get_yaxis_transform(), fontsize=10,
                bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.7))
    
    # Side panel: Mean proportions and statistics
    ax2 = axes[1]
    ax2.axis('off')
    
    # Create summary text
    summary_text = "Summary of Top Cell Types:\n\n"
    
    for i, celltype in enumerate(top_celltypes[:10]):  # Show first 10
        celltype_data = test_df[test_df['annotation'] == celltype]
        if not celltype_data.empty:
            row = celltype_data.iloc[0]
            summary_text += f"{i+1}. {celltype}:\n"
            summary_text += f"   Inner: {row['inner_mean']:.2f}%\n"
            summary_text += f"   Outer: {row['outer_mean']:.2f}%\n"
            summary_text += f"   FC: {row['fold_change']:.2f}x\n"
            summary_text += f"   p = {row['p_value']:.4f}\n\n"
    
    ax2.text(0.05, 0.95, summary_text, transform=ax2.transAxes,
             fontsize=8, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.7))
    
    plt.tight_layout()
    
    # Save figure
    png_path = output_dir / "proportion_heatmap.png"
    pdf_path = output_dir / "proportion_heatmap.pdf"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  Proportion heatmap saved: {png_path}")

def create_top_celltype_boxplots(test_df, sample_zone_prop_df, output_dir, 
                               top_n=12, figsize=(16, 12)):
    """
    Create box plots for top significant cell types
    """
    if test_df.empty or sample_zone_prop_df.empty:
        print("  Cannot create box plots: missing data")
        return
    
    # Get top significant cell types
    significant_df = test_df[test_df['significant'] == True].copy()
    
    if significant_df.empty:
        print("  No significant results for box plots")
        return
    
    # Sort by p-value
    significant_df = significant_df.sort_values('p_value')
    
    # Take top N
    top_celltypes = significant_df.head(min(top_n, len(significant_df)))['annotation'].tolist()
    
    if not top_celltypes:
        return
    
    # Determine grid layout
    n_cols = 4
    n_rows = (len(top_celltypes) + n_cols - 1) // n_cols
    
    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    # Flatten axes array
    axes_flat = axes.flatten()
    
    for idx, celltype in enumerate(top_celltypes):
        if idx >= len(axes_flat):
            break
        
        ax = axes_flat[idx]
        
        # Prepare data
        inner_data = sample_zone_prop_df[
            (sample_zone_prop_df['mucus_zone'] == 'inner') & 
            (sample_zone_prop_df[celltype].notna())
        ][celltype].values
        
        outer_data = sample_zone_prop_df[
            (sample_zone_prop_df['mucus_zone'] == 'outer') & 
            (sample_zone_prop_df[celltype].notna())
        ][celltype].values
        
        # Create box plot
        box_data = [inner_data, outer_data]
        box = ax.boxplot(box_data, labels=['Inner', 'Outer'], patch_artist=True)
        
        # Set colors
        colors = ['lightblue', 'lightcoral']
        for patch, color in zip(box['boxes'], colors):
            patch.set_facecolor(color)
        
        # Add individual points
        for i, data in enumerate(box_data, 1):
            x = np.random.normal(i, 0.04, size=len(data))
            ax.scatter(x, data, alpha=0.6, color='black', s=20, edgecolor='white')
        
        # Add significance annotation
        celltype_data = test_df[test_df['annotation'] == celltype]
        if not celltype_data.empty:
            p_val = celltype_data.iloc[0]['p_value']
            fc = celltype_data.iloc[0]['fold_change']
            
            # Format p-value
            if p_val < 0.001:
                p_text = 'p < 0.001'
            elif p_val < 0.01:
                p_text = f'p = {p_val:.3f}'
            else:
                p_text = f'p = {p_val:.2f}'
            
            # Add stars
            if p_val < 0.001:
                stars = '***'
            elif p_val < 0.01:
                stars = '**'
            elif p_val < 0.05:
                stars = '*'
            else:
                stars = ''
            
            # Add annotation
            y_max = max(np.max(inner_data) if len(inner_data) > 0 else 0,
                       np.max(outer_data) if len(outer_data) > 0 else 0)
            
            ax.text(1.5, y_max * 0.9, stars, ha='center', fontsize=12, 
                   color='red', fontweight='bold')
            ax.text(1.5, y_max * 0.8, p_text, ha='center', fontsize=8, 
                   bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))
        
        # Set title and labels
        ax.set_title(f'{celltype}\nFC = {fc:.2f}x' if not celltype_data.empty else celltype, 
                    fontsize=10, fontweight='bold')
        ax.set_ylabel('Proportion (%)', fontsize=8)
        ax.tick_params(axis='both', labelsize=8)
        ax.grid(True, alpha=0.3)
    
    # Hide unused subplots
    for idx in range(len(top_celltypes), len(axes_flat)):
        axes_flat[idx].axis('off')
    
    plt.suptitle(f'Top {len(top_celltypes)} Significantly Different Cell Types\nBox Plots with Individual Samples', 
                fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    # Save figure
    png_path = output_dir / "top_celltypes_boxplots.png"
    pdf_path = output_dir / "top_celltypes_boxplots.pdf"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  Top cell types box plots saved: {png_path}")

def create_correlation_plot(sample_zone_prop_df, test_df, output_dir, 
                          top_n=20, figsize=(14, 10)):
    """
    Create correlation plot between inner and outer proportions
    """
    if sample_zone_prop_df.empty:
        print("  Cannot create correlation plot: missing data")
        return
    
    # Get paired samples
    pivot_df = sample_zone_prop_df.pivot_table(
        index='sample', 
        columns='mucus_zone', 
        values=[col for col in sample_zone_prop_df.columns if col not in ['sample', 'mucus_zone', 'total_cells']],
        aggfunc='first'
    )
    
    # This is a MultiIndex, let's simplify
    inner_cols = [col[1] for col in pivot_df.columns if col[0] == 'inner']
    outer_cols = [col[1] for col in pivot_df.columns if col[0] == 'outer']
    
    # Find common cell types
    common_celltypes = list(set(inner_cols) & set(outer_cols))
    
    if not common_celltypes:
        print("  No common cell types for correlation plot")
        return
    
    # Get top significant cell types
    if not test_df.empty:
        significant_celltypes = test_df[test_df['significant'] == True]['annotation'].tolist()
        plot_celltypes = [ct for ct in significant_celltypes if ct in common_celltypes][:top_n]
    else:
        plot_celltypes = common_celltypes[:top_n]
    
    if not plot_celltypes:
        plot_celltypes = common_celltypes[:top_n]
    
    # Determine grid layout
    n_cols = 5
    n_rows = (len(plot_celltypes) + n_cols - 1) // n_cols
    
    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    # Flatten axes array
    axes_flat = axes.flatten()
    
    for idx, celltype in enumerate(plot_celltypes):
        if idx >= len(axes_flat):
            break
        
        ax = axes_flat[idx]
        
        # Get data
        inner_vals = []
        outer_vals = []
        
        for sample in pivot_df.index:
            inner_val = pivot_df.loc[sample, ('inner', celltype)] if ('inner', celltype) in pivot_df.columns else np.nan
            outer_val = pivot_df.loc[sample, ('outer', celltype)] if ('outer', celltype) in pivot_df.columns else np.nan
            
            if not np.isnan(inner_val) and not np.isnan(outer_val):
                inner_vals.append(inner_val)
                outer_vals.append(outer_val)
        
        if len(inner_vals) < 3:
            ax.text(0.5, 0.5, 'Insufficient\ndata', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=8)
            ax.set_title(celltype, fontsize=9)
            continue
        
        # Create scatter plot
        ax.scatter(inner_vals, outer_vals, alpha=0.6, s=30, color='blue', edgecolor='black')
        
        # Add diagonal line (y = x)
        max_val = max(max(inner_vals), max(outer_vals))
        ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.5, linewidth=1)
        
        # Calculate correlation
        r, p = stats.pearsonr(inner_vals, outer_vals)
        
        # Add correlation line
        if len(inner_vals) > 1:
            z = np.polyfit(inner_vals, outer_vals, 1)
            p_fit = np.poly1d(z)
            ax.plot(inner_vals, p_fit(inner_vals), "b-", alpha=0.5, linewidth=1)
        
        # Set labels and title
        ax.set_xlabel('Inner %', fontsize=8)
        ax.set_ylabel('Outer %', fontsize=8)
        ax.set_title(f'{celltype}\nr = {r:.2f}, p = {p:.3f}', fontsize=9, fontweight='bold')
        
        # Set equal aspect ratio
        ax.set_aspect('equal', adjustable='box')
        ax.grid(True, alpha=0.3)
    
    # Hide unused subplots
    for idx in range(len(plot_celltypes), len(axes_flat)):
        axes_flat[idx].axis('off')
    
    plt.suptitle('Correlation Between Inner and Outer Proportions\n' +
                'For Each Cell Type Across Samples', 
                fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    # Save figure
    png_path = output_dir / "inner_outer_correlation.png"
    pdf_path = output_dir / "inner_outer_correlation.pdf"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  Correlation plot saved: {png_path}")

def create_paired_comparison_plots(paired_df, sample_zone_prop_df, output_dir, 
                                 top_n=15, figsize=(14, 10)):
    """
    Create plots for paired analysis results
    """
    if paired_df.empty:
        print("  Cannot create paired comparison plots: missing data")
        return
    
    # Filter for cell types with significant paired results
    significant_paired = paired_df[
        (paired_df['test_performed']) & 
        (paired_df['significant_wilcoxon'] == True)
    ].copy()
    
    if significant_paired.empty:
        print("  No significant paired results")
        return
    
    # Sort by mean difference
    significant_paired = significant_paired.sort_values('mean_difference', key=abs, ascending=False)
    
    # Take top N
    top_paired = significant_paired.head(min(top_n, len(significant_paired)))
    
    # Create figure
    fig, axes = plt.subplots(2, 1, figsize=figsize, gridspec_kw={'height_ratios': [3, 1]})
    
    # Top plot: Paired differences
    ax1 = axes[0]
    
    x = np.arange(len(top_paired))
    width = 0.6
    
    # Get mean differences
    mean_diffs = top_paired['mean_difference'].values
    std_diffs = top_paired['std_difference'].values
    
    # Color by direction
    colors = ['green' if diff > 0 else 'red' for diff in mean_diffs]
    
    # Create bar plot
    bars = ax1.bar(x, mean_diffs, width, yerr=std_diffs, 
                  color=colors, alpha=0.7, capsize=5,
                  error_kw={'elinewidth': 1, 'capthick': 1})
    
    # Add horizontal line at 0
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    # Add significance stars
    for i, (idx, row) in enumerate(top_paired.iterrows()):
        p_val = row['p_value_wilcoxon']
        
        if p_val < 0.001:
            stars = '***'
        elif p_val < 0.01:
            stars = '**'
        elif p_val < 0.05:
            stars = '*'
        else:
            stars = ''
        
        if stars:
            height = row['mean_difference']
            ax1.text(i, height + (0.1 if height > 0 else -0.5), stars, 
                    ha='center', va='center' if height > 0 else 'top',
                    fontsize=10, color='red', fontweight='bold')
        
        # Add p-value
        p_text = f'p={p_val:.3f}' if p_val >= 0.001 else 'p<0.001'
        ax1.text(i, -0.5 if height > 0 else 0.5, p_text, ha='center', 
                va='top' if height > 0 else 'bottom', fontsize=7, rotation=45)
    
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Mean Difference (Inner - Outer) %')
    ax1.set_title('Paired Analysis: Inner - Outer Differences\n' +
                 f'Top {len(top_paired)} cell types with significant paired differences', 
                 fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(top_paired['annotation'], rotation=45, ha='right')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Bottom plot: Number of pairs
    ax2 = axes[1]
    
    n_pairs = top_paired['n_pairs'].values
    
    bars2 = ax2.bar(x, n_pairs, width, color='skyblue', alpha=0.7)
    
    ax2.set_xlabel('Cell Type')
    ax2.set_ylabel('Number of Paired Samples')
    ax2.set_title('Number of Samples with Both Inner and Outer Data', fontsize=10)
    ax2.set_xticks(x)
    ax2.set_xticklabels(top_paired['annotation'], rotation=45, ha='right')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add count labels
    for bar, count in zip(bars2, n_pairs):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{int(count)}', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    
    # Save figure
    png_path = output_dir / "paired_analysis_results.png"
    pdf_path = output_dir / "paired_analysis_results.pdf"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  Paired analysis plot saved: {png_path}")

def create_summary_report(test_df, paired_df, sample_zone_prop_df, output_dir):
    """
    Create comprehensive summary report
    """
    print(f"\n=== Creating summary report ===")
    
    # Create report lines
    report_lines = []
    report_lines.append("="*80)
    report_lines.append("INTEGRATED ANALYSIS REPORT: Inner vs Outer Mucus Zone Comparison")
    report_lines.append("="*80)
    report_lines.append(f"\nAnalysis date: {pd.Timestamp.now()}")
    
    # Basic statistics
    report_lines.append("\n" + "-"*80)
    report_lines.append("BASIC STATISTICS")
    report_lines.append("-"*80)
    
    n_samples = sample_zone_prop_df['sample'].nunique()
    n_inner_samples = len(sample_zone_prop_df[sample_zone_prop_df['mucus_zone'] == 'inner'])
    n_outer_samples = len(sample_zone_prop_df[sample_zone_prop_df['mucus_zone'] == 'outer'])
    
    report_lines.append(f"Total unique samples: {n_samples}")
    report_lines.append(f"Samples with inner data: {n_inner_samples}")
    report_lines.append(f"Samples with outer data: {n_outer_samples}")
    
    if not paired_df.empty:
        avg_pairs = paired_df['n_pairs'].mean() if 'n_pairs' in paired_df.columns else 0
        report_lines.append(f"Average paired samples per cell type: {avg_pairs:.1f}")
    
    # Cell type statistics
    n_celltypes_analyzed = len(test_df[test_df['test_performed'] == True])
    n_significant = len(test_df[(test_df['test_performed']) & (test_df['significant'] == True)])
    
    report_lines.append(f"\nCell types analyzed: {n_celltypes_analyzed}")
    report_lines.append(f"Significantly different cell types: {n_significant}")
    report_lines.append(f"Percentage significant: {n_significant/n_celltypes_analyzed*100:.1f}%")
    
    # Significant cell types breakdown
    if n_significant > 0:
        inner_enriched = len(test_df[(test_df['significant']) & (test_df['enrichment'] == 'inner')])
        outer_enriched = len(test_df[(test_df['significant']) & (test_df['enrichment'] == 'outer')])
        
        report_lines.append(f"  Enriched in inner zone: {inner_enriched}")
        report_lines.append(f"  Enriched in outer zone: {outer_enriched}")
    
    # Top significant cell types
    if n_significant > 0:
        report_lines.append("\n" + "-"*80)
        report_lines.append("TOP SIGNIFICANT CELL TYPES")
        report_lines.append("-"*80)
        
        # Sort by p-value
        significant_df = test_df[test_df['significant'] == True].copy()
        significant_df = significant_df.sort_values('p_value')
        
        for i, (idx, row) in enumerate(significant_df.head(10).iterrows(), 1):
            report_lines.append(f"{i}. {row['annotation']}:")
            report_lines.append(f"   Inner: {row['inner_mean']:.2f}%, Outer: {row['outer_mean']:.2f}%")
            report_lines.append(f"   Fold change: {row['fold_change']:.2f}x (log2FC: {row['log2_fold_change']:.2f})")
            report_lines.append(f"   p-value: {row['p_value']:.4e}, n_inner: {row['n_inner']}, n_outer: {row['n_outer']}")
    
    # Paired analysis summary
    if not paired_df.empty and 'significant_wilcoxon' in paired_df.columns:
        n_sig_paired = len(paired_df[paired_df['significant_wilcoxon'] == True])
        
        report_lines.append("\n" + "-"*80)
        report_lines.append("PAIRED ANALYSIS SUMMARY")
        report_lines.append("-"*80)
        report_lines.append(f"Cell types with significant paired differences: {n_sig_paired}")
        
        if n_sig_paired > 0:
            sig_paired = paired_df[paired_df['significant_wilcoxon'] == True].copy()
            sig_paired = sig_paired.sort_values('mean_difference', key=abs, ascending=False)
            
            for i, (idx, row) in enumerate(sig_paired.head(5).iterrows(), 1):
                report_lines.append(f"{i}. {row['annotation']}:")
                report_lines.append(f"   Mean difference: {row['mean_difference']:.2f}%")
                report_lines.append(f"   p-value (Wilcoxon): {row['p_value_wilcoxon']:.4e}")
                report_lines.append(f"   n_pairs: {row['n_pairs']}")
    
    # Methodological details
    report_lines.append("\n" + "-"*80)
    report_lines.append("METHODOLOGICAL DETAILS")
    report_lines.append("-"*80)
    
    if not test_df.empty:
        test_name = test_df['test_name'].iloc[0] if 'test_name' in test_df.columns else "Unknown"
        report_lines.append(f"Statistical test used: {test_name}")
    
    report_lines.append(f"Multiple testing correction: Benjamini-Hochberg FDR")
    report_lines.append(f"Significance threshold: p < 0.05")
    
    # Save report
    report_path = output_dir / "integrated_analysis_report.txt"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"Summary report saved: {report_path}")
    
    # Also create a shorter CSV summary
    summary_data = {
        'statistic': [
            'n_samples', 'n_inner_samples', 'n_outer_samples',
            'n_celltypes_analyzed', 'n_significant', 'percent_significant',
            'n_inner_enriched', 'n_outer_enriched'
        ],
        'value': [
            n_samples, n_inner_samples, n_outer_samples,
            n_celltypes_analyzed, n_significant, n_significant/n_celltypes_analyzed*100 if n_celltypes_analyzed > 0 else 0,
            inner_enriched if 'inner_enriched' in locals() else 0,
            outer_enriched if 'outer_enriched' in locals() else 0
        ]
    }
    
    summary_df = pd.DataFrame(summary_data)
    summary_path = output_dir / "analysis_summary.csv"
    summary_df.to_csv(summary_path, index=False, encoding='utf-8')
    
    print(f"Analysis summary saved: {summary_path}")

def main():
    """Main function"""
    # File path configuration
    h5ad_path = "extracted_samples.h5ad"  # Your h5ad file path
    output_dir = "integrated_inner_outer_comparison"  # Output directory
    
    print("="*60)
    print("Integrated Analysis: Inner vs Outer Cell Type Comparison")
    print("="*60)
    
    try:
        # 1. Load data
        adata = load_and_prepare_data(h5ad_path)
        if adata is None:
            sys.exit(1)
        
        # Check if we have inner and outer zones
        zones = adata.obs['mucus_zone'].unique()
        if 'inner' not in zones or 'outer' not in zones:
            print(f"\nError: Data must contain both 'inner' and 'outer' zones")
            print(f"Available zones: {zones}")
            sys.exit(1)
        
        # 2. Perform integrated analysis
        print(f"\nStarting integrated analysis...")
        test_df, paired_df, sample_zone_prop_df = analyze_all_samples_comparison(
            adata, output_dir, 
            min_mean_proportion=0.1,  # Minimum 0.1% mean proportion
            test_type='mannwhitney',  # Mann-Whitney U test (non-parametric)
            perform_paired_analysis_flag=True
        )
        
        if test_df is None:
            print("\nAnalysis failed.")
            sys.exit(1)
        
        # 3. Print final summary
        print(f"\n" + "="*60)
        print("ANALYSIS COMPLETE!")
        print("="*60)
        
        n_significant = len(test_df[test_df['significant'] == True]) if not test_df.empty else 0
        n_analyzed = len(test_df[test_df['test_performed'] == True]) if not test_df.empty else 0
        
        print(f"\nResults Summary:")
        print(f"  Cell types analyzed: {n_analyzed}")
        print(f"  Significantly different: {n_significant}")
        
        if n_significant > 0:
            # Count enrichment directions
            inner_enriched = len(test_df[(test_df['significant']) & (test_df['enrichment'] == 'inner')])
            outer_enriched = len(test_df[(test_df['significant']) & (test_df['enrichment'] == 'outer')])
            
            print(f"    Enriched in inner zone: {inner_enriched}")
            print(f"    Enriched in outer zone: {outer_enriched}")
            
            # Show top 5
            top_significant = test_df[test_df['significant'] == True].sort_values('p_value').head(5)
            print(f"\n  Top 5 most significant cell types:")
            for i, (idx, row) in enumerate(top_significant.iterrows(), 1):
                print(f"    {i}. {row['annotation']}: p={row['p_value']:.4e}, FC={row['fold_change']:.2f}x")
        
        print(f"\nOutput directory: {output_dir}")
        print(f"\nGenerated files:")
        
        # List generated files
        output_path = Path(output_dir)
        file_count = 0
        for file_path in sorted(output_path.rglob("*")):
            if file_path.is_file():
                file_size = file_path.stat().st_size / 1024  # Size in KB
                print(f"  {file_path.relative_to(output_path.parent)} ({file_size:.1f} KB)")
                file_count += 1
        
        print(f"\nTotal files generated: {file_count}")
        print(f"\nAnalysis complete!")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
