#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cell type group comparison analysis tool
Statistical analysis of cell type differences between groups
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
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import kruskal, mannwhitneyu
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
    
    # Print celltype information
    print(f"\nCell type distribution (top 20):")
    celltype_counts = adata.obs['annotation'].value_counts()
    for i, (celltype, count) in enumerate(celltype_counts.head(20).items()):
        percentage = count / adata.n_obs * 100
        print(f"  {i+1}. {celltype}: {count} cells ({percentage:.2f}%)")
    
    return adata

def calculate_celltype_proportions(adata, min_cells_per_sample=10):
    """
    Calculate cell type proportions for each sample
    
    Parameters:
    adata: AnnData object
    min_cells_per_sample: minimum cell count threshold per sample
    
    Returns:
    prop_df: proportion DataFrame
    sample_info: sample information DataFrame
    """
    print(f"\n=== Calculating cell type proportions ===")
    
    # Get all samples
    all_samples = adata.obs['sample'].unique()
    all_celltypes = adata.obs['annotation'].unique()
    
    print(f"Number of samples: {len(all_samples)}")
    print(f"Number of cell types: {len(all_celltypes)}")
    
    # Calculate cell type counts for each sample
    celltype_counts_list = []
    sample_info_list = []
    
    for sample in all_samples:
        # Get data for current sample
        sample_mask = adata.obs['sample'] == sample
        sample_data = adata[sample_mask]
        
        # Skip samples with too few cells
        if len(sample_data) < min_cells_per_sample:
            print(f"  Skipping sample {sample}: only {len(sample_data)} cells (< {min_cells_per_sample})")
            continue
        
        # Get group information for sample
        sample_group = sample_data.obs['mucus_zone'].iloc[0] if len(sample_data.obs['mucus_zone'].unique()) == 1 else sample_data.obs['mucus_zone'].iloc[0]
        
        # Count each cell type
        celltype_counts = sample_data.obs['annotation'].value_counts()
        
        # Calculate proportions
        total_cells = len(sample_data)
        
        # Create row data
        row_data = {'sample': sample, 'mucus_zone': sample_group, 'total_cells': total_cells}
        for celltype in all_celltypes:
            count = celltype_counts.get(celltype, 0)
            proportion = count / total_cells * 100 if total_cells > 0 else 0
            row_data[celltype] = proportion
        
        celltype_counts_list.append(row_data)
        
        # Save sample information
        sample_info_list.append({
            'sample': sample,
            'mucus_zone': sample_group,
            'total_cells': total_cells
        })
    
    # Create DataFrames
    prop_df = pd.DataFrame(celltype_counts_list)
    sample_info_df = pd.DataFrame(sample_info_list)
    
    # Set index
    prop_df.set_index('sample', inplace=True)
    
    print(f"Valid samples: {len(prop_df)}")
    print(f"Samples excluded: {len(all_samples) - len(prop_df)}")
    
    return prop_df, sample_info_df

def perform_statistical_tests(prop_df, celltype, group_col='mucus_zone'):
    """
    Perform statistical tests for specified cell type
    
    Parameters:
    prop_df: proportion DataFrame
    celltype: cell type name
    group_col: group column name
    
    Returns:
    Statistical results dictionary
    """
    # Extract data for current cell type
    data = prop_df[[celltype, group_col]].copy()
    data.columns = ['proportion', 'mucus_zone']
    
    # Remove missing values
    data = data.dropna()
    
    # Get groups
    groups = data['mucus_zone'].unique()
    
    if len(groups) < 2:
        return {
            'annotation': celltype,
            'n_groups': len(groups),
            'test_performed': False,
            'message': f'Only {len(groups)} group(s), cannot perform group comparison'
        }
    
    # Extract data for each group
    group_data = []
    for group in groups:
        group_values = data[data['mucus_zone'] == group]['proportion'].values
        group_data.append(group_values)
    
    # Perform statistical tests
    if len(groups) == 2:
        # Two-group comparison: Mann-Whitney U test
        group1_data = group_data[0]
        group2_data = group_data[1]
        
        if len(group1_data) > 0 and len(group2_data) > 0:
            stat, p_value = mannwhitneyu(group1_data, group2_data, alternative='two-sided')
            test_name = "Mann-Whitney U test"
        else:
            return {
                'annotation': celltype,
                'n_groups': len(groups),
                'test_performed': False,
                'message': 'Insufficient data in some groups'
            }
    
    else:
        # Multi-group comparison: Kruskal-Wallis test
        if all(len(g) > 0 for g in group_data):
            stat, p_value = kruskal(*group_data)
            test_name = "Kruskal-Wallis test"
            
            # If significant, perform post-hoc test
            if p_value < 0.05:
                # Prepare data for post-hoc test
                posthoc_data = []
                for group in groups:
                    group_values = data[data['mucus_zone'] == group]['proportion'].values
                    for value in group_values:
                        posthoc_data.append([value, group])
                
                posthoc_df = pd.DataFrame(posthoc_data, columns=['proportion', 'mucus_zone'])
                
                # Perform Tukey HSD test
                tukey = pairwise_tukeyhsd(endog=posthoc_df['proportion'],
                                         groups=posthoc_df['mucus_zone'],
                                         alpha=0.05)
                posthoc_results = str(tukey.summary())
            else:
                posthoc_results = None
        else:
            return {
                'annotation': celltype,
                'n_groups': len(groups),
                'test_performed': False,
                'message': 'Insufficient data in some groups'
            }
    
    # Calculate basic statistics
    summary_stats = {}
    for group in groups:
        group_values = data[data['mucus_zone'] == group]['proportion'].values
        if len(group_values) > 0:
            summary_stats[group] = {
                'mean': np.mean(group_values),
                'median': np.median(group_values),
                'std': np.std(group_values),
                'min': np.min(group_values),
                'max': np.max(group_values),
                'n_samples': len(group_values)
            }
        else:
            summary_stats[group] = {
                'mean': np.nan,
                'median': np.nan,
                'std': np.nan,
                'min': np.nan,
                'max': np.nan,
                'n_samples': 0
            }
    
    result = {
        'annotation': celltype,
        'n_groups': len(groups),
        'test_performed': True,
        'test_name': test_name,
        'statistic': stat,
        'p_value': p_value,
        'significant': p_value < 0.05,
        'summary_stats': summary_stats
    }
    
    if len(groups) > 2 and p_value < 0.05 and 'posthoc_results' in locals():
        result['posthoc_results'] = posthoc_results
    
    return result

def create_boxplot(prop_df, celltype, output_dir, stats_result=None, figsize=(10, 6)):
    """
    Create boxplot for specified cell type
    
    Parameters:
    prop_df: proportion DataFrame
    celltype: cell type name
    output_dir: output directory
    stats_result: statistical results
    figsize: figure size
    """
    # Prepare data
    plot_data = prop_df[[celltype, 'mucus_zone']].copy()
    plot_data.columns = ['proportion', 'mucus_zone']
    
    # Remove missing values
    plot_data = plot_data.dropna()
    
    if len(plot_data) == 0:
        print(f"  {celltype}: No valid data, skipping plot")
        return None
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(figsize[0] * 1.5, figsize[1]), 
                             gridspec_kw={'width_ratios': [3, 1]})
    
    # Main plot: Boxplot
    ax_main = axes[0]
    
    # Create boxplot
    box_plot = ax_main.boxplot([plot_data[plot_data['mucus_zone'] == group]['proportion'].values 
                                for group in plot_data['mucus_zone'].unique()],
                               labels=plot_data['mucus_zone'].unique(),
                               patch_artist=True)
    
    # Set colors for boxes
    colors = plt.cm.Set3(np.linspace(0, 1, len(plot_data['mucus_zone'].unique())))
    for patch, color in zip(box_plot['boxes'], colors):
        patch.set_facecolor(color)
    
    # Add individual data points
    for i, group in enumerate(plot_data['mucus_zone'].unique(), 1):
        group_data = plot_data[plot_data['mucus_zone'] == group]['proportion'].values
        x_pos = np.random.normal(i, 0.04, size=len(group_data))
        ax_main.scatter(x_pos, group_data, alpha=0.6, color='black', s=30, edgecolor='white')
    
    # Set title and labels
    title = f'{celltype} Distribution Across Groups\n'
    if stats_result and stats_result.get('test_performed', False):
        p_value = stats_result['p_value']
        significance = 'Significant' if stats_result['significant'] else 'Not significant'
        title += f'{stats_result["test_name"]}: p = {p_value:.4e} ({significance})'
    
    ax_main.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax_main.set_xlabel('Group', fontsize=12)
    ax_main.set_ylabel('Cell Proportion (%)', fontsize=12)
    ax_main.tick_params(axis='x', rotation=45)
    
    # Add statistical significance annotation
    if stats_result and stats_result.get('test_performed', False) and stats_result['significant']:
        # Find max y value for annotation positioning
        y_max = plot_data['proportion'].max()
        y_range = plot_data['proportion'].max() - plot_data['proportion'].min()
        
        # Add star annotation
        ax_main.text(0.5, 0.95, f'p = {stats_result["p_value"]:.4e} ***',
                    transform=ax_main.transAxes, fontsize=12,
                    ha='center', va='top',
                    bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))
    
    # Side plot: Statistical summary
    ax_stats = axes[1]
    ax_stats.axis('off')
    
    if stats_result and 'summary_stats' in stats_result:
        summary_text = f"Statistical Summary:\n"
        summary_text += f"Cell Type: {celltype}\n"
        summary_text += f"Test: {stats_result.get('test_name', 'N/A')}\n"
        summary_text += f"p-value: {stats_result.get('p_value', 'N/A'):.4e}\n"
        summary_text += f"Significant: {'Yes' if stats_result.get('significant', False) else 'No'}\n\n"
        
        for group, stats_dict in stats_result['summary_stats'].items():
            summary_text += f"Group {group}:\n"
            summary_text += f"  Samples: {stats_dict['n_samples']}\n"
            summary_text += f"  Mean: {stats_dict['mean']:.2f}%\n"
            summary_text += f"  Median: {stats_dict['median']:.2f}%\n"
            summary_text += f"  Std: {stats_dict['std']:.2f}\n"
        
        ax_stats.text(0.05, 0.95, summary_text, transform=ax_stats.transAxes,
                     fontsize=9, verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.7))
    
    plt.tight_layout()
    
    # Save figure
    safe_celltype = str(celltype).replace('/', '_').replace(' ', '_').replace('\\', '_')
    
    # Save as PNG
    png_path = output_dir / f'{safe_celltype}_boxplot.png'
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    
    # Save as PDF
    pdf_path = output_dir / f'{safe_celltype}_boxplot.pdf'
    plt.savefig(pdf_path, bbox_inches='tight')
    
    plt.close(fig)
    
    return png_path, pdf_path

def create_violinplot(prop_df, celltype, output_dir, stats_result=None, figsize=(10, 6)):
    """
    Create violin plot for specified cell type
    
    Parameters:
    prop_df: proportion DataFrame
    celltype: cell type name
    output_dir: output directory
    stats_result: statistical results
    figsize: figure size
    """
    # Prepare data
    plot_data = prop_df[[celltype, 'mucus_zone']].copy()
    plot_data.columns = ['proportion', 'mucus_zone']
    
    # Remove missing values
    plot_data = plot_data.dropna()
    
    if len(plot_data) == 0:
        return None
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create violin plot
    violin_parts = ax.violinplot([plot_data[plot_data['mucus_zone'] == group]['proportion'].values 
                                  for group in plot_data['mucus_zone'].unique()],
                                 showmeans=True, showmedians=True)
    
    # Add individual data points
    for i, group in enumerate(plot_data['mucus_zone'].unique(), 1):
        group_data = plot_data[plot_data['mucus_zone'] == group]['proportion'].values
        x_pos = np.random.normal(i, 0.04, size=len(group_data))
        ax.scatter(x_pos, group_data, alpha=0.6, color='black', s=30, edgecolor='white')
    
    # Set x-tick labels
    ax.set_xticks(range(1, len(plot_data['mucus_zone'].unique()) + 1))
    ax.set_xticklabels(plot_data['mucus_zone'].unique())
    
    # Set title and labels
    title = f'{celltype} Distribution (Violin Plot)\n'
    if stats_result and stats_result.get('test_performed', False):
        p_value = stats_result['p_value']
        significance = 'Significant' if stats_result['significant'] else 'Not significant'
        title += f'{stats_result["test_name"]}: p = {p_value:.4e} ({significance})'
    
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Group', fontsize=12)
    ax.set_ylabel('Cell Proportion (%)', fontsize=12)
    
    plt.tight_layout()
    
    # Save figure
    safe_celltype = str(celltype).replace('/', '_').replace(' ', '_').replace('\\', '_')
    
    # Save as PNG
    png_path = output_dir / f'{safe_celltype}_violinplot.png'
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    
    # Save as PDF
    pdf_path = output_dir / f'{safe_celltype}_violinplot.pdf'
    plt.savefig(pdf_path, bbox_inches='tight')
    
    plt.close(fig)
    
    return png_path, pdf_path

def analyze_celltype_differences(adata, output_dir="celltype_group_comparison", 
                                 min_proportion=0.1, min_samples_per_group=3,
                                 create_violin_plots=False):
    """
    Analyze cell type differences between groups
    
    Parameters:
    adata: AnnData object
    output_dir: output directory
    min_proportion: minimum average proportion threshold
    min_samples_per_group: minimum samples per group
    create_violin_plots: whether to create violin plots
    """
    print(f"\n=== Analyzing cell type differences between groups ===")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Calculate cell type proportions
    prop_df, sample_info_df = calculate_celltype_proportions(adata)
    
    if prop_df.empty:
        print("Error: No valid sample data")
        return None
    
    # Save proportion data
    prop_output_path = output_path / "celltype_proportions.csv"
    prop_df.to_csv(prop_output_path, encoding='utf-8')
    print(f"Cell type proportion data saved: {prop_output_path}")
    
    # Save sample information
    sample_info_path = output_path / "sample_info.csv"
    sample_info_df.to_csv(sample_info_path, encoding='utf-8', index=False)
    print(f"Sample information saved: {sample_info_path}")
    
    # Get all cell types
    all_celltypes = [col for col in prop_df.columns if col not in ['mucus_zone', 'total_cells']]
    
    # Filter rare cell types
    filtered_celltypes = []
    for celltype in all_celltypes:
        # Calculate average proportion
        mean_proportion = prop_df[celltype].mean()
        
        # Check if each group has enough samples
        groups_with_enough_samples = []
        for group in prop_df['mucus_zone'].unique():
            group_mask = prop_df['mucus_zone'] == group
            group_samples = group_mask.sum()
            
            if group_samples >= min_samples_per_group:
                groups_with_enough_samples.append(group)
        
        # Keep if average proportion above threshold and at least two groups have enough samples
        if mean_proportion >= min_proportion and len(groups_with_enough_samples) >= 2:
            filtered_celltypes.append(celltype)
        else:
            print(f"  Skipping {celltype}: mean proportion {mean_proportion:.2f}% < {min_proportion}% "
                  f"or less than {min_samples_per_group} samples/group")
    
    print(f"\nAnalyzing {len(filtered_celltypes)} cell types (filtered out {len(all_celltypes) - len(filtered_celltypes)} rare types)")
    
    # Create plot directories
    boxplots_dir = output_path / "boxplots"
    boxplots_dir.mkdir(exist_ok=True)
    
    if create_violin_plots:
        violinplots_dir = output_path / "violinplots"
        violinplots_dir.mkdir(exist_ok=True)
    
    # Perform statistical tests and create plots for each cell type
    stats_results = []
    plot_files = []
    
    for i, celltype in enumerate(filtered_celltypes, 1):
        print(f"\nProcessing cell type {i}/{len(filtered_celltypes)}: {celltype}")
        
        # Perform statistical test
        stats_result = perform_statistical_tests(prop_df, celltype)
        stats_results.append(stats_result)
        
        # Create boxplot
        if stats_result.get('test_performed', False):
            png_path, pdf_path = create_boxplot(prop_df, celltype, boxplots_dir, stats_result)
            if png_path and pdf_path:
                plot_files.append({
                    'annotation': celltype,
                    'plot_type': 'boxplot',
                    'png_path': png_path,
                    'pdf_path': pdf_path,
                    'p_value': stats_result['p_value'],
                    'significant': stats_result['significant']
                })
        
        # Create violin plot if requested
        if create_violin_plots and stats_result.get('test_performed', False):
            v_png_path, v_pdf_path = create_violinplot(prop_df, celltype, violinplots_dir, stats_result)
            if v_png_path and v_pdf_path:
                plot_files.append({
                    'annotation': celltype,
                    'plot_type': 'violinplot',
                    'png_path': v_png_path,
                    'pdf_path': v_pdf_path,
                    'p_value': stats_result['p_value'],
                    'significant': stats_result['significant']
                })
    
    # Save statistical results
    stats_df = pd.DataFrame(stats_results)
    stats_path = output_path / "statistical_results.csv"
    stats_df.to_csv(stats_path, encoding='utf-8', index=False)
    print(f"\nStatistical results saved: {stats_path}")
    
    # Create summary report
    create_summary_report(stats_results, prop_df, output_path)
    
    # Create summary plot for significant cell types
    create_significant_summary_plot(stats_results, prop_df, output_path)
    
    return stats_results, prop_df, plot_files

def create_summary_report(stats_results, prop_df, output_dir):
    """
    Create summary report
    
    Parameters:
    stats_results: statistical results list
    prop_df: proportion DataFrame
    output_dir: output directory
    """
    print(f"\n=== Creating summary report ===")
    
    # Extract significant results
    significant_results = [r for r in stats_results if r.get('significant', False)]
    non_significant_results = [r for r in stats_results if not r.get('significant', True)]
    
    # Create report
    report_lines = []
    report_lines.append("="*60)
    report_lines.append("Cell Type Group Comparison Analysis Report")
    report_lines.append("="*60)
    report_lines.append(f"\nAnalysis time: {pd.Timestamp.now()}")
    report_lines.append(f"Total cell types analyzed: {len(stats_results)}")
    report_lines.append(f"Significantly different cell types: {len(significant_results)}")
    report_lines.append(f"Non-significant cell types: {len(non_significant_results)}")
    
    if significant_results:
        report_lines.append("\n" + "="*60)
        report_lines.append("Significantly Different Cell Types:")
        report_lines.append("="*60)
        
        # Sort by p-value
        significant_results_sorted = sorted(significant_results, key=lambda x: x.get('p_value', 1))
        
        for i, result in enumerate(significant_results_sorted, 1):
            celltype = result['annotation']
            p_value = result.get('p_value', 'N/A')
            test_name = result.get('test_name', 'N/A')
            
            report_lines.append(f"\n{i}. {celltype}")
            report_lines.append(f"   Test: {test_name}")
            report_lines.append(f"   p-value: {p_value:.4e}")
            
            # Add group comparisons
            if 'summary_stats' in result:
                report_lines.append("   Group statistics:")
                for group, stats in result['summary_stats'].items():
                    report_lines.append(f"     Group {group}: mean={stats['mean']:.2f}%, "
                                      f"samples={stats['n_samples']}")
    
    # Save report
    report_path = output_dir / "analysis_report.txt"
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"Analysis report saved: {report_path}")
    
    # Create CSV of significant results
    if significant_results:
        sig_df = pd.DataFrame([
            {
                'annotation': r['annotation'],
                'test': r.get('test_name', ''),
                'p_value': r.get('p_value', ''),
                'significant': r.get('significant', False)
            }
            for r in significant_results
        ])
        sig_path = output_dir / "significant_results.csv"
        sig_df.to_csv(sig_path, encoding='utf-8', index=False)
        print(f"Significant results list saved: {sig_path}")

def create_significant_summary_plot(stats_results, prop_df, output_dir, 
                                    top_n=20, figsize=(14, 10)):
    """
    Create summary plot for significant cell types
    
    Parameters:
    stats_results: statistical results list
    prop_df: proportion DataFrame
    output_dir: output directory
    top_n: show top N most significant
    figsize: figure size
    """
    # Extract significant results
    significant_results = [r for r in stats_results if r.get('significant', False)]
    
    if not significant_results:
        print("No significantly different cell types, skipping summary plot")
        return
    
    # Sort by p-value, take top N
    significant_results_sorted = sorted(significant_results, 
                                       key=lambda x: x.get('p_value', 1))
    top_results = significant_results_sorted[:min(top_n, len(significant_results_sorted))]
    
    # Prepare data
    plot_data_list = []
    
    for result in top_results:
        celltype = result['annotation']
        p_value = result.get('p_value', 1)
        
        # Get data for this cell type
        celltype_data = prop_df[[celltype, 'mucus_zone']].copy()
        celltype_data.columns = ['proportion', 'mucus_zone']
        celltype_data['annotation'] = celltype
        celltype_data['log10_p_value'] = -np.log10(p_value) if p_value > 0 else 20
        
        plot_data_list.append(celltype_data)
    
    if not plot_data_list:
        return
    
    plot_data = pd.concat(plot_data_list, ignore_index=True)
    
    # Create figure
    fig, axes = plt.subplots(2, 1, figsize=figsize, 
                            gridspec_kw={'height_ratios': [3, 1]})
    
    # Top subplot: Heatmap of mean proportions
    ax_heatmap = axes[0]
    
    # Calculate mean proportion for each cell type in each group
    heatmap_data = []
    for result in top_results:
        celltype = result['annotation']
        row = {'annotation': celltype, 'log10_p_value': -np.log10(result.get('p_value', 1))}
        
        for group in prop_df['mucus_zone'].unique():
            group_mask = prop_df['mucus_zone'] == group
            mean_prop = prop_df.loc[group_mask, celltype].mean()
            row[group] = mean_prop
        
        heatmap_data.append(row)
    
    heatmap_df = pd.DataFrame(heatmap_data)
    heatmap_df.set_index('annotation', inplace=True)
    
    # Create heatmap
    im = ax_heatmap.imshow(heatmap_df.drop('log10_p_value', axis=1).values, 
                          aspect='auto', cmap='YlOrRd')
    
    # Set axes
    ax_heatmap.set_yticks(range(len(heatmap_df)))
    ax_heatmap.set_yticklabels(heatmap_df.index, fontsize=10)
    ax_heatmap.set_xticks(range(len(heatmap_df.columns)-1))
    ax_heatmap.set_xticklabels([col for col in heatmap_df.columns if col != 'log10_p_value'], 
                              fontsize=10, rotation=45, ha='right')
    
    ax_heatmap.set_title(f'Top {len(top_results)} Significantly Different Cell Types - Mean Proportions (%)', 
                        fontsize=14, fontweight='bold', pad=20)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax_heatmap, fraction=0.046, pad=0.04)
    cbar.set_label('Mean Cell Proportion (%)', fontsize=12)
    
    # Bottom subplot: p-value information
    ax_pvalues = axes[1]
    ax_pvalues.axis('off')
    
    # Create p-value table
    pvalue_text = "Cell Type Significance Ranking:\n"
    for i, result in enumerate(top_results, 1):
        celltype = result['annotation']
        p_value = result.get('p_value', 1)
        test_name = result.get('test_name', 'N/A')
        
        pvalue_text += f"{i}. {celltype}: p = {p_value:.2e} ({test_name})\n"
    
    ax_pvalues.text(0.05, 0.95, pvalue_text, transform=ax_pvalues.transAxes,
                   fontsize=9, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    plt.tight_layout()
    
    # Save figure
    summary_png = output_dir / "significant_celltypes_summary.png"
    summary_pdf = output_dir / "significant_celltypes_summary.pdf"
    
    plt.savefig(summary_png, dpi=300, bbox_inches='tight')
    plt.savefig(summary_pdf, bbox_inches='tight')
    plt.close()
    
    print(f"Significant cell types summary plot saved: {summary_png}")

def main():
    """Main function"""
    # File path configuration
    h5ad_path = "extracted_samples.h5ad"  # Replace with your h5ad file path
    output_dir = "celltype_group_analysis"  # Output directory
    create_violin_plots = False  # Set to True to also create violin plots
    
    print("="*60)
    print("Cell Type Group Comparison Analysis Tool")
    print("="*60)
    
    try:
        # 1. Load data
        adata = load_and_prepare_data(h5ad_path)
        if adata is None:
            sys.exit(1)
        
        # 2. Analyze cell type differences
        stats_results, prop_df, plot_files = analyze_celltype_differences(
            adata, output_dir, min_proportion=0.1, min_samples_per_group=3,
            create_violin_plots=create_violin_plots
        )
        
        if stats_results is None:
            sys.exit(1)
        
        # 3. Print final summary
        print(f"\n" + "="*60)
        print("Analysis Complete!")
        print("="*60)
        print(f"Input file: {h5ad_path}")
        print(f"Output directory: {output_dir}")
        print(f"Cell types analyzed: {len(stats_results)}")
        print(f"Boxplots generated: {sum(1 for p in plot_files if p['plot_type'] == 'boxplot')}")
        
        if create_violin_plots:
            print(f"Violin plots generated: {sum(1 for p in plot_files if p['plot_type'] == 'violinplot')}")
        
        # Count significant results
        significant_count = sum(1 for r in stats_results if r.get('significant', False))
        print(f"Significantly different cell types: {significant_count}")
        
        # Show file list
        print(f"\nGenerated files:")
        output_path = Path(output_dir)
        for file_type in ['csv', 'txt', 'png', 'pdf']:
            for file_path in output_path.rglob(f"*.{file_type}"):
                if file_path.is_file():
                    relative_path = file_path.relative_to(output_path.parent)
                    print(f"  {relative_path}")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
