#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sample-wise cell type comparison between inner and outer mucus zones
Statistical analysis of cell type differences for each sample
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
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import chi2_contingency, fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests
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
    required_cols = ['mucus_zone', 'sample', 'CNs_10']
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
    print(f"Number of cell types: {adata.obs['CNs_10'].nunique()}")
    
    # Print group information
    print(f"\nGroup distribution:")
    group_counts = adata.obs['mucus_zone'].value_counts()
    for group, count in group_counts.items():
        percentage = count / adata.n_obs * 100
        print(f"  {group}: {count} cells ({percentage:.2f}%)")
    
    # Print sample information
    print(f"\nSample distribution:")
    sample_counts = adata.obs['sample'].value_counts()
    for sample, count in sample_counts.head(20).items():
        percentage = count / adata.n_obs * 100
        group = adata.obs[adata.obs['sample'] == sample]['mucus_zone'].iloc[0] if not adata.obs[adata.obs['sample'] == sample]['mucus_zone'].empty else "Unknown"
        print(f"  {sample}: {count} cells ({percentage:.2f}%), group: {group}")
    
    return adata

def analyze_sample_inner_outer_comparison(adata, output_dir="sample_wise_comparison", 
                                          min_cells_per_zone=10, 
                                          pval_threshold=0.05,
                                          create_plots=True):
    """
    Analyze inner vs outer cell type composition for each sample
    
    Parameters:
    adata: AnnData object
    output_dir: output directory
    min_cells_per_zone: minimum cells required in each zone
    pval_threshold: p-value threshold for significance
    create_plots: whether to create visualizations
    """
    print(f"\n=== Analyzing inner vs outer comparison for each sample ===")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Get unique samples
    all_samples = adata.obs['sample'].unique()
    all_celltypes = adata.obs['CNs_10'].unique()
    
    print(f"Total samples: {len(all_samples)}")
    print(f"Total cell types: {len(all_celltypes)}")
    
    # Create directories for plots
    if create_plots:
        barplots_dir = output_path / "barplots"
        barplots_dir.mkdir(exist_ok=True)
        heatmaps_dir = output_path / "heatmaps"
        heatmaps_dir.mkdir(exist_ok=True)
        stacked_bar_dir = output_path / "stacked_bars"
        stacked_bar_dir.mkdir(exist_ok=True)
    
    # Dictionary to store results
    all_results = []
    summary_stats = []
    sample_comparison_details = []
    
    # Analyze each sample
    for i, sample in enumerate(all_samples, 1):
        print(f"\nProcessing sample {i}/{len(all_samples)}: {sample}")
        
        # Filter data for this sample
        sample_mask = adata.obs['sample'] == sample
        sample_data = adata[sample_mask]
        
        if len(sample_data) == 0:
            print(f"  No data for sample {sample}, skipping...")
            continue
        
        # Check if sample has both inner and outer zones
        zones_in_sample = sample_data.obs['mucus_zone'].unique()
        
        if len(zones_in_sample) < 2:
            print(f"  Sample {sample} has only {len(zones_in_sample)} zone(s): {zones_in_sample}, skipping comparison...")
            continue
        
        # Separate inner and outer data
        inner_mask = sample_data.obs['mucus_zone'] == 'inner'
        outer_mask = sample_data.obs['mucus_zone'] == 'outer'
        
        inner_data = sample_data[inner_mask] if inner_mask.any() else None
        outer_data = sample_data[outer_mask] if outer_mask.any() else None
        
        if inner_data is None or outer_data is None:
            print(f"  Sample {sample} missing inner or outer data, skipping...")
            continue
        
        inner_cells = len(inner_data)
        outer_cells = len(outer_data)
        
        print(f"  Inner cells: {inner_cells}, Outer cells: {outer_cells}")
        
        if inner_cells < min_cells_per_zone or outer_cells < min_cells_per_zone:
            print(f"  Sample {sample} has insufficient cells (< {min_cells_per_zone}) in one zone, skipping...")
            continue
        
        # Calculate cell type proportions
        inner_counts = inner_data.obs['CNs_10'].value_counts()
        outer_counts = outer_data.obs['CNs_10'].value_counts()
        
        inner_proportions = (inner_counts / inner_cells * 100).to_dict()
        outer_proportions = (outer_counts / outer_cells * 100).to_dict()
        
        # Get all cell types present in either zone
        all_celltypes_in_sample = set(inner_counts.index.tolist() + outer_counts.index.tolist())
        
        # Prepare data for statistical tests
        test_results = []
        
        for celltype in all_celltypes_in_sample:
            inner_count = inner_counts.get(celltype, 0)
            outer_count = outer_counts.get(celltype, 0)
            
            # Create 2x2 contingency table
            # Row 0: inner, Row 1: outer
            # Column 0: celltype present, Column 1: other cell types
            other_inner = inner_cells - inner_count
            other_outer = outer_cells - outer_count
            
            contingency_table = np.array([[inner_count, other_inner],
                                          [outer_count, other_outer]])
            
            # Fisher's exact test for cell type enrichment
            try:
                if inner_count + outer_count > 0:
                    oddsratio, p_value_fisher = fisher_exact(contingency_table, alternative='two-sided')
                    
                    # Chi-squared test
                    chi2, p_value_chi2, dof, expected = chi2_contingency(contingency_table, correction=True)
                    
                    # Calculate fold change
                    inner_prop = (inner_count + 1) / (inner_cells + 1)  # Add pseudocount
                    outer_prop = (outer_count + 1) / (outer_cells + 1)  # Add pseudocount
                    
                    if outer_prop > 0:
                        fold_change = inner_prop / outer_prop
                        log2_fold_change = np.log2(fold_change)
                    else:
                        fold_change = np.inf
                        log2_fold_change = np.inf
                    
                    # Determine enrichment direction
                    if inner_prop > outer_prop:
                        enrichment = "inner"
                    elif inner_prop < outer_prop:
                        enrichment = "outer"
                    else:
                        enrichment = "equal"
                    
                    test_results.append({
                        'sample': sample,
                        'CNs_10': celltype,
                        'inner_count': inner_count,
                        'outer_count': outer_count,
                        'inner_proportion': inner_proportions.get(celltype, 0),
                        'outer_proportion': outer_proportions.get(celltype, 0),
                        'fold_change': fold_change,
                        'log2_fold_change': log2_fold_change,
                        'p_value_fisher': p_value_fisher,
                        'p_value_chi2': p_value_chi2,
                        'odds_ratio': oddsratio,
                        'enrichment': enrichment,
                        'total_cells': inner_cells + outer_cells
                    })
                    
                    # Save detailed comparison
                    sample_comparison_details.append({
                        'sample': sample,
                        'celltype': celltype,
                        'inner_cells': inner_count,
                        'outer_cells': outer_count,
                        'inner_total': inner_cells,
                        'outer_total': outer_cells,
                        'inner_prop': inner_proportions.get(celltype, 0),
                        'outer_prop': outer_proportions.get(celltype, 0),
                        'fold_change': fold_change,
                        'p_value': p_value_fisher
                    })
            except Exception as e:
                print(f"    Error testing {celltype}: {e}")
        
        if not test_results:
            print(f"  No valid test results for sample {sample}")
            continue
        
        # Convert to DataFrame
        results_df = pd.DataFrame(test_results)
        
        # Multiple testing correction
        results_df['p_value_fisher_adj'] = multipletests(results_df['p_value_fisher'], 
                                                         method='fdr_bh')[1]
        results_df['p_value_chi2_adj'] = multipletests(results_df['p_value_chi2'], 
                                                      method='fdr_bh')[1]
        
        # Mark significant results
        results_df['significant_fisher'] = results_df['p_value_fisher_adj'] < pval_threshold
        results_df['significant_chi2'] = results_df['p_value_chi2_adj'] < pval_threshold
        
        # Save per-sample results
        sample_results_path = output_path / f"{sample}_celltype_comparison.csv"
        results_df.to_csv(sample_results_path, index=False, encoding='utf-8')
        
        # Add to all results
        all_results.append(results_df)
        
        # Create visualizations
        if create_plots and len(results_df) > 0:
            # 1. Bar plot comparing proportions
            create_comparison_barplot(results_df, sample, barplots_dir, pval_threshold)
            
            # 2. Heatmap of fold changes
            create_fold_change_heatmap(results_df, sample, heatmaps_dir, pval_threshold)
            
            # 3. Stacked bar plot
            create_stacked_barplot(inner_data, outer_data, sample, stacked_bar_dir)
        
        # Create summary statistics
        n_significant_fisher = results_df['significant_fisher'].sum()
        n_significant_chi2 = results_df['significant_chi2'].sum()
        
        # Find most enriched cell types
        if len(results_df) > 0:
            inner_enriched = results_df[(results_df['significant_fisher']) & 
                                        (results_df['enrichment'] == 'inner')]
            outer_enriched = results_df[(results_df['significant_fisher']) & 
                                        (results_df['enrichment'] == 'outer')]
            
            summary_stats.append({
                'sample': sample,
                'inner_cells': inner_cells,
                'outer_cells': outer_cells,
                'total_cells': inner_cells + outer_cells,
                'n_celltypes': len(all_celltypes_in_sample),
                'n_significant_fisher': n_significant_fisher,
                'n_significant_chi2': n_significant_chi2,
                'n_inner_enriched': len(inner_enriched),
                'n_outer_enriched': len(outer_enriched),
                'most_enriched_inner': inner_enriched.nlargest(1, 'fold_change')['CNs_10'].iloc[0] if len(inner_enriched) > 0 else 'None',
                'most_enriched_outer': outer_enriched.nlargest(1, 'fold_change')['CNs_10'].iloc[0] if len(outer_enriched) > 0 else 'None',
                'max_fold_change': results_df['fold_change'].max() if not results_df.empty else 0
            })
        
        print(f"  Analyzed {len(results_df)} cell types, {n_significant_fisher} significantly different (Fisher's exact)")
    
    # Combine all results
    if all_results:
        combined_results = pd.concat(all_results, ignore_index=True)
        
        # Save combined results
        combined_path = output_path / "all_samples_combined_results.csv"
        combined_results.to_csv(combined_path, index=False, encoding='utf-8')
        print(f"\nCombined results saved: {combined_path}")
        
        # Create summary statistics
        summary_df = pd.DataFrame(summary_stats)
        summary_path = output_path / "sample_summary_statistics.csv"
        summary_df.to_csv(summary_path, index=False, encoding='utf-8')
        print(f"Summary statistics saved: {summary_path}")
        
        # Save detailed comparisons
        details_df = pd.DataFrame(sample_comparison_details)
        details_path = output_path / "detailed_comparisons.csv"
        details_df.to_csv(details_path, index=False, encoding='utf-8')
        print(f"Detailed comparisons saved: {details_path}")
        
        # Create summary plots
        create_summary_plots(summary_df, combined_results, output_path, pval_threshold)
        
        # Create meta-analysis across all samples
        perform_meta_analysis(combined_results, output_path, pval_threshold)
        
        return combined_results, summary_df, details_df
    else:
        print("\nNo valid comparisons found. Check if samples have both inner and outer zones.")
        return None, None, None

def create_comparison_barplot(results_df, sample_name, output_dir, pval_threshold=0.05, 
                             top_n=20, figsize=(14, 8)):
    """
    Create bar plot comparing inner vs outer proportions
    """
    # Sort by fold change
    results_sorted = results_df.sort_values('fold_change', ascending=False)
    
    # Take top N (both directions)
    top_results = results_sorted.head(top_n)
    
    if len(top_results) == 0:
        return
    
    # Create figure
    fig, axes = plt.subplots(2, 1, figsize=figsize, 
                            gridspec_kw={'height_ratios': [3, 1]})
    
    # Top subplot: Bar plot
    ax1 = axes[0]
    
    x = np.arange(len(top_results))
    width = 0.35
    
    # Get proportions
    inner_props = top_results['inner_proportion'].values
    outer_props = top_results['outer_proportion'].values
    
    # Create bars
    rects1 = ax1.bar(x - width/2, inner_props, width, label='Inner', color='lightblue', alpha=0.8)
    rects2 = ax1.bar(x + width/2, outer_props, width, label='Outer', color='lightcoral', alpha=0.8)
    
    # Add significance markers
    for i, (_, row) in enumerate(top_results.iterrows()):
        if row['significant_fisher']:
            # Get max value for positioning
            max_val = max(row['inner_proportion'], row['outer_proportion'])
            ax1.text(i, max_val + 0.5, '*', ha='center', va='bottom', fontsize=12, color='red')
    
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Proportion (%)')
    ax1.set_title(f'Sample {sample_name}: Inner vs Outer Cell Type Proportions\n' +
                  f'Top {len(top_results)} cell types by fold change', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(top_results['CNs_10'], rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for rect in rects1:
        height = rect.get_height()
        if height > 0.5:  # Only label if > 0.5%
            ax1.text(rect.get_x() + rect.get_width()/2., height + 0.1,
                    f'{height:.1f}', ha='center', va='bottom', fontsize=8)
    
    for rect in rects2:
        height = rect.get_height()
        if height > 0.5:  # Only label if > 0.5%
            ax1.text(rect.get_x() + rect.get_width()/2., height + 0.1,
                    f'{height:.1f}', ha='center', va='bottom', fontsize=8)
    
    # Bottom subplot: Fold change
    ax2 = axes[1]
    
    # Get fold changes
    fold_changes = top_results['log2_fold_change'].values
    
    # Create bar plot for fold change
    colors = ['lightgreen' if fc > 0 else 'lightgray' for fc in fold_changes]
    rects3 = ax2.bar(x, fold_changes, color=colors, alpha=0.7)
    
    # Add horizontal line at 0
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    ax2.set_xlabel('Cell Type')
    ax2.set_ylabel('log2(Fold Change)')
    ax2.set_title('Log2 Fold Change (Inner/Outer)', fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels(top_results['CNs_10'], rotation=45, ha='right')
    ax2.grid(True, alpha=0.3)
    
    # Add fold change values
    for rect, fc in zip(rects3, fold_changes):
        height = rect.get_height()
        if abs(height) > 0.1:  # Only label if > 0.1
            ax2.text(rect.get_x() + rect.get_width()/2., 
                    height + (0.1 if height > 0 else -0.3),
                    f'{fc:.2f}', ha='center', va='center' if height > 0 else 'top', 
                    fontsize=8)
    
    plt.tight_layout()
    
    # Save figure
    safe_sample_name = str(sample_name).replace('/', '_').replace(' ', '_').replace('\\', '_')
    png_path = output_dir / f"{safe_sample_name}_comparison.png"
    pdf_path = output_dir / f"{safe_sample_name}_comparison.pdf"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)

def create_fold_change_heatmap(results_df, sample_name, output_dir, pval_threshold=0.05,
                              top_n=30, figsize=(12, 10)):
    """
    Create heatmap of fold changes
    """
    # Sort by significance and fold change
    results_sorted = results_df.sort_values(['significant_fisher', 'fold_change'], 
                                           ascending=[False, False])
    
    # Take top N
    top_results = results_sorted.head(top_n)
    
    if len(top_results) == 0:
        return
    
    # Prepare data for heatmap
    heatmap_data = []
    annotations = []
    p_values = []
    
    for _, row in top_results.iterrows():
        heatmap_data.append([row['log2_fold_change']])
        annotations.append(row['CNs_10'])
        p_values.append(row['p_value_fisher_adj'])
    
    heatmap_matrix = np.array(heatmap_data)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap
    im = ax.imshow(heatmap_matrix, cmap='coolwarm', aspect='auto',
                  vmin=-max(abs(heatmap_matrix.min()), abs(heatmap_matrix.max())),
                  vmax=max(abs(heatmap_matrix.min()), abs(heatmap_matrix.max())))
    
    # Set ticks
    ax.set_yticks(range(len(annotations)))
    ax.set_yticklabels(annotations, fontsize=10)
    ax.set_xticks([])
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('log2(Fold Change)\nInner/Outer', fontsize=12)
    
    # Add p-value annotation
    for i, p_val in enumerate(p_values):
        if p_val < pval_threshold:
            ax.text(0, i, '★', ha='center', va='center', 
                   color='gold', fontsize=12, fontweight='bold')
    
    # Set title
    ax.set_title(f'Sample {sample_name}: Cell Type Enrichment (Inner vs Outer)\n' +
                 f'★ = p < {pval_threshold} (FDR corrected)', 
                 fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    
    # Save figure
    safe_sample_name = str(sample_name).replace('/', '_').replace(' ', '_').replace('\\', '_')
    png_path = output_dir / f"{safe_sample_name}_fold_change_heatmap.png"
    pdf_path = output_dir / f"{safe_sample_name}_fold_change_heatmap.pdf"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)

def create_stacked_barplot(inner_data, outer_data, sample_name, output_dir, 
                          top_n=15, figsize=(12, 8)):
    """
    Create stacked bar plot comparing cell type composition
    """
    # Get cell type proportions
    inner_counts = inner_data.obs['CNs_10'].value_counts()
    outer_counts = outer_data.obs['CNs_10'].value_counts()
    
    inner_props = inner_counts / len(inner_data) * 100
    outer_props = outer_counts / len(outer_data) * 100
    
    # Combine and get top N cell types
    all_celltypes = set(inner_counts.index.tolist() + outer_counts.index.tolist())
    celltype_sums = {}
    
    for ct in all_celltypes:
        celltype_sums[ct] = inner_props.get(ct, 0) + outer_props.get(ct, 0)
    
    # Get top N cell types
    top_celltypes = sorted(celltype_sums.items(), key=lambda x: x[1], reverse=True)[:top_n]
    top_celltypes = [ct[0] for ct in top_celltypes]
    
    # Prepare data for plotting
    inner_top = [inner_props.get(ct, 0) for ct in top_celltypes]
    outer_top = [outer_props.get(ct, 0) for ct in top_celltypes]
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    # Plot 1: Stacked bar
    x = np.arange(len(top_celltypes))
    
    ax1 = axes[0]
    ax1.bar(x, inner_top, label='Inner', color='lightblue', alpha=0.8)
    ax1.bar(x, outer_top, bottom=inner_top, label='Outer', color='lightcoral', alpha=0.8)
    
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Proportion (%)')
    ax1.set_title(f'Sample {sample_name}: Cell Type Composition\nStacked Bar Plot', 
                 fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(top_celltypes, rotation=45, ha='right')
    ax1.legend()
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Plot 2: Side-by-side bar
    ax2 = axes[1]
    width = 0.35
    
    rects1 = ax2.bar(x - width/2, inner_top, width, label='Inner', color='lightblue', alpha=0.8)
    rects2 = ax2.bar(x + width/2, outer_top, width, label='Outer', color='lightcoral', alpha=0.8)
    
    ax2.set_xlabel('Cell Type')
    ax2.set_ylabel('Proportion (%)')
    ax2.set_title('Side-by-Side Comparison', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(top_celltypes, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    # Save figure
    safe_sample_name = str(sample_name).replace('/', '_').replace(' ', '_').replace('\\', '_')
    png_path = output_dir / f"{safe_sample_name}_composition.png"
    pdf_path = output_dir / f"{safe_sample_name}_composition.pdf"
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)

def create_summary_plots(summary_df, combined_results, output_dir, pval_threshold=0.05):
    """
    Create summary plots across all samples
    """
    if summary_df.empty or combined_results.empty:
        return
    
    # 1. Bar plot: Number of significant cell types per sample
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Significant cell types per sample
    ax1 = axes[0, 0]
    samples = summary_df['sample']
    n_significant = summary_df['n_significant_fisher']
    
    bars = ax1.bar(range(len(samples)), n_significant, color='skyblue')
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Number of Significant Cell Types')
    ax1.set_title(f'Significant Cell Types per Sample\n(p < {pval_threshold}, FDR corrected)', 
                 fontsize=12, fontweight='bold')
    ax1.set_xticks(range(len(samples)))
    ax1.set_xticklabels(samples, rotation=45, ha='right')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                f'{int(height)}', ha='center', va='bottom', fontsize=8)
    
    # Plot 2: Inner vs Outer cell counts
    ax2 = axes[0, 1]
    x = np.arange(len(summary_df))
    width = 0.35
    
    bars_inner = ax2.bar(x - width/2, summary_df['inner_cells'], width, 
                        label='Inner', color='lightblue')
    bars_outer = ax2.bar(x + width/2, summary_df['outer_cells'], width, 
                        label='Outer', color='lightcoral')
    
    ax2.set_xlabel('Sample')
    ax2.set_ylabel('Number of Cells')
    ax2.set_title('Cell Counts per Sample', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(summary_df['sample'], rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Enrichment pattern heatmap
    ax3 = axes[1, 0]
    
    # Get most consistently enriched cell types
    enrichment_counts = combined_results.groupby(['CNs_10', 'enrichment']).size().unstack(fill_value=0)
    if 'inner' in enrichment_counts.columns and 'outer' in enrichment_counts.columns:
        # Calculate net enrichment score
        enrichment_counts['net_score'] = enrichment_counts.get('inner', 0) - enrichment_counts.get('outer', 0)
        top_enriched = enrichment_counts.nlargest(20, 'net_score').index.tolist()
        bottom_enriched = enrichment_counts.nsmallest(20, 'net_score').index.tolist()
        
        # Combine and get unique
        top_celltypes = list(set(top_enriched + bottom_enriched))[:20]
        
        # Create heatmap data
        heatmap_data = []
        for ct in top_celltypes:
            row = []
            for sample in summary_df['sample'].head(10):  # Top 10 samples
                sample_data = combined_results[(combined_results['sample'] == sample) & 
                                              (combined_results['CNs_10'] == ct)]
                if not sample_data.empty:
                    if sample_data.iloc[0]['significant_fisher']:
                        if sample_data.iloc[0]['enrichment'] == 'inner':
                            row.append(1)  # Enriched in inner
                        else:
                            row.append(-1)  # Enriched in outer
                    else:
                        row.append(0)  # Not significant
                else:
                    row.append(0)  # Not present
            
            heatmap_data.append(row)
        
        if heatmap_data:
            im = ax3.imshow(heatmap_data, cmap='coolwarm', aspect='auto',
                          vmin=-1, vmax=1)
            
            ax3.set_yticks(range(len(top_celltypes)))
            ax3.set_yticklabels(top_celltypes, fontsize=8)
            ax3.set_xticks(range(min(10, len(summary_df))))
            ax3.set_xticklabels(summary_df['sample'].head(10), rotation=45, ha='right', fontsize=8)
            
            cbar = plt.colorbar(im, ax=ax3, fraction=0.046, pad=0.04)
            cbar.set_ticks([-1, 0, 1])
            cbar.set_ticklabels(['Outer', 'NS', 'Inner'])
            
            ax3.set_title('Cell Type Enrichment Pattern\n(Top 20 most variable)', 
                         fontsize=12, fontweight='bold')
    
    # Plot 4: Fold change distribution
    ax4 = axes[1, 1]
    
    # Filter significant results
    significant_results = combined_results[combined_results['significant_fisher']]
    
    if not significant_results.empty:
        inner_enriched = significant_results[significant_results['enrichment'] == 'inner']['log2_fold_change']
        outer_enriched = significant_results[significant_results['enrichment'] == 'outer']['log2_fold_change']
        
        bins = np.linspace(-5, 5, 51)
        ax4.hist(inner_enriched, bins=bins, alpha=0.5, label='Inner Enriched', color='blue')
        ax4.hist(outer_enriched, bins=bins, alpha=0.5, label='Outer Enriched', color='red')
        
        ax4.axvline(x=0, color='black', linestyle='--', linewidth=0.5)
        ax4.set_xlabel('log2(Fold Change)')
        ax4.set_ylabel('Frequency')
        ax4.set_title('Distribution of Significant Fold Changes', fontsize=12, fontweight='bold')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save summary plot
    summary_png = output_dir / "summary_statistics.png"
    summary_pdf = output_dir / "summary_statistics.pdf"
    
    plt.savefig(summary_png, dpi=300, bbox_inches='tight')
    plt.savefig(summary_pdf, bbox_inches='tight')
    plt.close(fig)
    
    # 2. Meta-analysis summary
    create_meta_analysis_plot(combined_results, output_dir)

def create_meta_analysis_plot(combined_results, output_dir):
    """
    Create meta-analysis summary plot
    """
    if combined_results.empty:
        return
    
    # Aggregate results by cell type
    celltype_stats = []
    
    for celltype in combined_results['CNs_10'].unique():
        ct_data = combined_results[combined_results['CNs_10'] == celltype]
        
        if len(ct_data) > 1:
            # Perform random-effects meta-analysis of fold changes
            log_fcs = ct_data['log2_fold_change'].values
            weights = 1 / (ct_data['p_value_fisher'] + 1e-10)  # Weight by inverse p-value
            
            weighted_mean = np.average(log_fcs, weights=weights)
            n_samples = len(ct_data)
            n_sig = ct_data['significant_fisher'].sum()
            mean_inner = ct_data['inner_proportion'].mean()
            mean_outer = ct_data['outer_proportion'].mean()
            
            celltype_stats.append({
                'CNs_10': celltype,
                'n_samples': n_samples,
                'n_significant': n_sig,
                'weighted_mean_fc': weighted_mean,
                'mean_inner_prop': mean_inner,
                'mean_outer_prop': mean_outer,
                'mean_fold_change': np.mean(ct_data['fold_change']),
                'median_fold_change': np.median(ct_data['fold_change'])
            })
    
    if not celltype_stats:
        return
    
    stats_df = pd.DataFrame(celltype_stats)
    stats_df = stats_df.sort_values('n_significant', ascending=False)
    
    # Create meta-analysis plot
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    
    # Plot 1: Number of significant samples per cell type
    ax1 = axes[0]
    top_stats = stats_df.head(20)
    
    y_pos = np.arange(len(top_stats))
    ax1.barh(y_pos, top_stats['n_significant'], color='lightgreen')
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(top_stats['CNs_10'])
    ax1.set_xlabel('Number of Significant Samples')
    ax1.set_title('Most Consistently Different Cell Types\n(Number of samples with p < 0.05)', 
                 fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='x')
    
    # Plot 2: Fold change consistency
    ax2 = axes[1]
    
    # Filter for cell types with at least 3 significant samples
    consistent_celltypes = stats_df[stats_df['n_significant'] >= 3]
    
    if len(consistent_celltypes) > 0:
        top_consistent = consistent_celltypes.head(15)
        
        x = np.arange(len(top_consistent))
        width = 0.35
        
        bars1 = ax2.bar(x - width/2, top_consistent['mean_inner_prop'], width, 
                       label='Mean Inner %', color='lightblue')
        bars2 = ax2.bar(x + width/2, top_consistent['mean_outer_prop'], width, 
                       label='Mean Outer %', color='lightcoral')
        
        ax2.set_xlabel('Cell Type')
        ax2.set_ylabel('Mean Proportion (%)')
        ax2.set_title('Average Proportions for Consistently Different Cell Types\n(≥3 significant samples)', 
                     fontsize=12, fontweight='bold')
        ax2.set_xticks(x)
        ax2.set_xticklabels(top_consistent['CNs_10'], rotation=45, ha='right')
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    # Save meta-analysis plot
    meta_png = output_dir / "meta_analysis_summary.png"
    meta_pdf = output_dir / "meta_analysis_summary.pdf"
    
    plt.savefig(meta_png, dpi=300, bbox_inches='tight')
    plt.savefig(meta_pdf, bbox_inches='tight')
    plt.close(fig)
    
    # Save meta-analysis results
    meta_results_path = output_dir / "meta_analysis_results.csv"
    stats_df.to_csv(meta_results_path, index=False, encoding='utf-8')

def perform_meta_analysis(combined_results, output_dir, pval_threshold=0.05):
    """
    Perform meta-analysis across all samples
    """
    print("\n=== Performing meta-analysis across all samples ===")
    
    if combined_results.empty:
        print("No data for meta-analysis")
        return
    
    # Prepare meta-analysis data
    meta_data = []
    
    for celltype in combined_results['CNs_10'].unique():
        ct_data = combined_results[combined_results['CNs_10'] == celltype]
        
        if len(ct_data) > 1:
            # Calculate summary statistics
            n_samples = len(ct_data)
            n_sig = ct_data['significant_fisher'].sum()
            sig_percentage = n_sig / n_samples * 100
            
            # Calculate weighted average fold change
            weights = 1 / (ct_data['p_value_fisher'] + 1e-10)
            weighted_fc = np.average(ct_data['log2_fold_change'], weights=weights)
            
            # Calculate consistency score
            # Positive if enriched in same direction in most samples
            inner_enriched = len(ct_data[ct_data['enrichment'] == 'inner'])
            outer_enriched = len(ct_data[ct_data['enrichment'] == 'outer'])
            
            if inner_enriched > outer_enriched:
                consistency = 'inner'
                consistency_score = inner_enriched / n_samples
            elif outer_enriched > inner_enriched:
                consistency = 'outer'
                consistency_score = outer_enriched / n_samples
            else:
                consistency = 'mixed'
                consistency_score = 0.5
            
            meta_data.append({
                'CNs_10': celltype,
                'n_samples': n_samples,
                'n_significant': n_sig,
                'percent_significant': sig_percentage,
                'mean_inner_prop': ct_data['inner_proportion'].mean(),
                'mean_outer_prop': ct_data['outer_proportion'].mean(),
                'mean_fold_change': np.mean(ct_data['fold_change']),
                'median_fold_change': np.median(ct_data['fold_change']),
                'weighted_mean_log2fc': weighted_fc,
                'consistent_enrichment': consistency,
                'consistency_score': consistency_score,
                'min_p_value': ct_data['p_value_fisher'].min(),
                'max_fold_change': ct_data['fold_change'].max(),
                'min_fold_change': ct_data['fold_change'].min()
            })
    
    if not meta_data:
        print("No valid data for meta-analysis")
        return
    
    # Create DataFrame
    meta_df = pd.DataFrame(meta_data)
    
    # Sort by consistency and significance
    meta_df = meta_df.sort_values(['consistency_score', 'percent_significant'], ascending=False)
    
    # Save meta-analysis results
    meta_path = output_dir / "celltype_meta_analysis.csv"
    meta_df.to_csv(meta_path, index=False, encoding='utf-8')
    print(f"Meta-analysis results saved: {meta_path}")
    
    # Identify consistently different cell types
    consistent_celltypes = meta_df[
        (meta_df['n_significant'] >= 3) & 
        (meta_df['consistency_score'] >= 0.7)  # At least 70% consistent
    ]
    
    print(f"\nConsistently different cell types (≥3 significant samples, ≥70% consistent):")
    print(f"Total: {len(consistent_celltypes)} cell types")
    
    for _, row in consistent_celltypes.iterrows():
        print(f"  {row['CNs_10']}: {row['n_significant']}/{row['n_samples']} samples significant, "
              f"enriched in {row['consistent_enrichment']} ({row['consistency_score']:.1%})")
    
    # Create meta-analysis report
    report_lines = []
    report_lines.append("="*60)
    report_lines.append("META-ANALYSIS REPORT: Inner vs Outer Mucus Zone Comparison")
    report_lines.append("="*60)
    report_lines.append(f"\nAnalysis time: {pd.Timestamp.now()}")
    report_lines.append(f"Total samples analyzed: {combined_results['sample'].nunique()}")
    report_lines.append(f"Total cell types analyzed: {len(meta_df)}")
    report_lines.append(f"\nConsistency threshold: ≥3 significant samples, ≥70% agreement")
    report_lines.append(f"Cell types meeting threshold: {len(consistent_celltypes)}")
    
    if len(consistent_celltypes) > 0:
        report_lines.append("\n" + "-"*60)
        report_lines.append("CONSISTENTLY ENRICHED IN INNER ZONE:")
        report_lines.append("-"*60)
        
        inner_enriched = consistent_celltypes[consistent_celltypes['consistent_enrichment'] == 'inner']
        for _, row in inner_enriched.iterrows():
            report_lines.append(f"{row['CNs_10']}: {row['n_significant']}/{row['n_samples']} samples "
                              f"({row['percent_significant']:.1f}%), "
                              f"mean fold change: {row['mean_fold_change']:.2f}")
        
        report_lines.append("\n" + "-"*60)
        report_lines.append("CONSISTENTLY ENRICHED IN OUTER ZONE:")
        report_lines.append("-"*60)
        
        outer_enriched = consistent_celltypes[consistent_celltypes['consistent_enrichment'] == 'outer']
        for _, row in outer_enriched.iterrows():
            report_lines.append(f"{row['CNs_10']}: {row['n_significant']}/{row['n_samples']} samples "
                              f"({row['percent_significant']:.1f}%), "
                              f"mean fold change: {row['mean_fold_change']:.2f}")
    
    # Save report
    report_path = output_dir / "meta_analysis_report.txt"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"\nMeta-analysis report saved: {report_path}")
    
    return meta_df

def main():
    """Main function"""
    # File path configuration
    h5ad_path = "extracted_samples.h5ad"  # Your h5ad file path
    output_dir = "sample_wise_inner_outer_comparison_CN"  # Output directory
    
    print("="*60)
    print("Sample-wise Inner vs Outer Cell Type Comparison Tool")
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
        
        # 2. Analyze each sample
        print(f"\nStarting analysis of {adata.obs['sample'].nunique()} samples...")
        combined_results, summary_df, details_df = analyze_sample_inner_outer_comparison(
            adata, output_dir, min_cells_per_zone=10, pval_threshold=0.05, create_plots=True
        )
        
        if combined_results is None:
            print("\nNo valid comparisons could be performed.")
            print("Possible reasons:")
            print("1. Samples don't have both inner and outer zones")
            print("2. Too few cells in one or both zones")
            print("3. No cell types meet the analysis criteria")
            sys.exit(1)
        
        # 3. Print final summary
        print(f"\n" + "="*60)
        print("ANALYSIS COMPLETE!")
        print("="*60)
        
        n_samples_analyzed = summary_df['sample'].nunique() if summary_df is not None else 0
        n_comparisons = len(combined_results) if combined_results is not None else 0
        
        print(f"\nSummary:")
        print(f"  Samples analyzed: {n_samples_analyzed}")
        print(f"  Total comparisons: {n_comparisons}")
        
        if summary_df is not None:
            total_significant = summary_df['n_significant_fisher'].sum()
            avg_significant = summary_df['n_significant_fisher'].mean()
            print(f"  Total significant findings: {total_significant}")
            print(f"  Average per sample: {avg_significant:.1f}")
        
        print(f"\nOutput directory: {output_dir}")
        print(f"\nGenerated files:")
        
        # List generated files
        output_path = Path(output_dir)
        for file_path in sorted(output_path.rglob("*")):
            if file_path.is_file():
                file_size = file_path.stat().st_size / 1024  # Size in KB
                print(f"  {file_path.relative_to(output_path.parent)} ({file_size:.1f} KB)")
        
        # Show most consistently different cell types
        if combined_results is not None:
            print(f"\nTop 5 most frequently different cell types:")
            
            # Count how many samples each cell type is significant in
            sig_counts = combined_results[combined_results['significant_fisher']].groupby('CNs_10').size()
            if len(sig_counts) > 0:
                for i, (celltype, count) in enumerate(sig_counts.nlargest(5).items(), 1):
                    total_samples = combined_results[combined_results['CNs_10'] == celltype]['sample'].nunique()
                    print(f"  {i}. {celltype}: {count}/{total_samples} samples")
            
        print(f"\nAnalysis complete!")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
