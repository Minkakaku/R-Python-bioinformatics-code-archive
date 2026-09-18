import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy.sparse import issparse
from scipy.spatial import KDTree
from scipy.stats import pearsonr, spearmanr, kruskal
from scipy.ndimage import gaussian_filter1d
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# ================= 配置 =================
input_file = "All_filter_state_samp_new_with_CNs10-15.h5ad"
output_dir = "distance_expression_analysis"
genes_to_analyze = ["MYC", "AXIN2", "RNF43"]
target_celltype = "WNT+ ISC"
source_celltype = "matri-prod CAF"
max_distance = 250
bin_width = 25  # 改为25um
state_column = 'state'  # state列名，请根据实际数据调整
# ========================================

os.makedirs(output_dir, exist_ok=True)

print(f"Reading {input_file}...")
adata = sc.read_h5ad(input_file)
print(f"Data shape: {adata.shape}")

# ---- 检查细胞类型 ----
if 'celltype' not in adata.obs.columns:
    raise KeyError("'celltype' column not found in adata.obs.")

if target_celltype not in adata.obs['celltype'].values:
    available = adata.obs['celltype'].unique().tolist()
    raise ValueError(f"'{target_celltype}' not found. Available: {available[:10]}...")

if source_celltype not in adata.obs['celltype'].values:
    available = adata.obs['celltype'].unique().tolist()
    raise ValueError(f"'{source_celltype}' not found. Available: {available[:10]}...")

# ---- 检查state列 ----
if state_column not in adata.obs.columns:
    print(f"⚠ Warning: '{state_column}' not found in adata.obs.")
    print(f"Available columns: {adata.obs.columns.tolist()}")
    print("Will skip state-specific analysis.")
    has_state = False
else:
    has_state = True
    print(f"✓ Found state column: '{state_column}'")
    print(f"  States: {adata.obs[state_column].unique().tolist()}")

# ---- 自动识别样本列和空间坐标 ----
sample_col = None
for col in ['sample', 'Sample', 'patient', 'Patient', 'batch', 'Batch', 'orig.ident', 'slice', 'Slide']:
    if col in adata.obs.columns:
        sample_col = col
        break
if sample_col is None:
    print("No sample column found, treating all cells as one sample.")
    adata.obs['sample'] = 'All_Samples'
    sample_col = 'sample'
print(f"\nSample column detected: '{sample_col}'")

spatial_key = None
for key in ['spatial', 'X_spatial', 'coordinates', 'obsm/spatial', 'X_umap']:
    if key in adata.obsm:
        spatial_key = key
        break
if spatial_key is None:
    for key in adata.obsm.keys():
        if adata.obsm[key].shape[1] >= 2:
            spatial_key = key
            break
if spatial_key is None:
    raise RuntimeError("No 2D spatial coordinates found in adata.obsm.")
print(f"Spatial coordinates: '{spatial_key}'")

# ---- 准备表达矩阵 ----
if adata.raw is not None:
    print("\nUsing adata.raw as expression matrix...")
    expr_adata = adata.raw.to_adata()
    expr_adata.obsm = adata.obsm.copy()
    expr_adata.obs = adata.obs.copy()
else:
    print("No raw found, normalizing + log1p...")
    expr_adata = adata.copy()
    sc.pp.normalize_total(expr_adata, target_sum=1e4)
    sc.pp.log1p(expr_adata)

# ---- 检查基因存在性 ----
valid_genes = {}
for g in genes_to_analyze:
    if g in expr_adata.var_names:
        valid_genes[g] = g
        print(f"✓ Found '{g}'")
    else:
        match = [v for v in expr_adata.var_names if v.upper() == g.upper()]
        if match:
            valid_genes[g] = match[0]
            print(f"✓ Case-insensitive match: '{g}' -> '{match[0]}'")
        else:
            print(f"✗ Gene '{g}' not found!")
            valid_genes[g] = None

missing = [g for g in genes_to_analyze if valid_genes.get(g) is None]
if missing:
    raise SystemExit(f"Missing required genes: {missing}")

# ---- 初始化结果存储 ----
all_results = []

print(f"\nProcessing samples...")
samples = expr_adata.obs[sample_col].unique()
print(f"Found {len(samples)} sample(s)")

for sample in samples:
    print(f"\n--- Sample: {sample} ---")
    
    sample_mask = expr_adata.obs[sample_col] == sample
    sample_data = expr_adata[sample_mask].copy()
    
    target_mask = sample_data.obs['celltype'] == target_celltype
    source_mask = sample_data.obs['celltype'] == source_celltype
    
    n_target = target_mask.sum()
    n_source = source_mask.sum()
    print(f"  {target_celltype}: {n_target} cells")
    print(f"  {source_celltype}: {n_source} cells")
    
    if n_target == 0:
        print(f"  ⚠ No target cells found in this sample, skipping...")
        continue
    if n_source == 0:
        print(f"  ⚠ No source cells found in this sample, skipping...")
        continue
    
    target_coords = sample_data[target_mask].obsm[spatial_key][:, :2]
    source_coords = sample_data[source_mask].obsm[spatial_key][:, :2]
    
    kdtree = KDTree(source_coords)
    distances, indices = kdtree.query(target_coords, k=1)
    
    print(f"  Distance range: {distances.min():.2f} - {distances.max():.2f}")
    
    target_expr = sample_data[target_mask]
    gene_exprs = {}
    for g in genes_to_analyze:
        gene_name = valid_genes[g]
        X = target_expr[:, gene_name].X
        if issparse(X):
            gene_exprs[g] = X.toarray().flatten()
        else:
            gene_exprs[g] = np.asarray(X).flatten()
    
    # 获取state信息
    if has_state:
        target_states = sample_data[target_mask].obs[state_column].values
    
    for i, dist in enumerate(distances):
        result_dict = {
            'sample': sample,
            'distance': dist
        }
        for g in genes_to_analyze:
            result_dict[f'{g}_expr'] = gene_exprs[g][i]
        if has_state:
            result_dict['state'] = target_states[i]
        all_results.append(result_dict)

if len(all_results) == 0:
    raise RuntimeError("No valid cell pairs found in any sample!")

results_df = pd.DataFrame(all_results)
print(f"\nTotal analyzed cells: {len(results_df)}")
print(f"Distance statistics:")
print(f"  Mean: {results_df['distance'].mean():.2f}")
print(f"  Median: {results_df['distance'].median():.2f}")
print(f"  Range: {results_df['distance'].min():.2f} - {results_df['distance'].max():.2f}")

# ---- 过滤距离范围 ----
filtered_df = results_df[results_df['distance'] <= max_distance].copy()
print(f"\nCells within {max_distance} units: {len(filtered_df)} ({len(filtered_df)/len(results_df)*100:.1f}%)")

if len(filtered_df) < 10:
    print("⚠ Warning: Very few cells within the specified distance range!")

# ---- 创建距离分箱（25um宽度） ----
bins = np.arange(0, max_distance + bin_width, bin_width)
filtered_df['distance_bin'] = pd.cut(filtered_df['distance'], bins=bins, right=False)

print(f"\nNumber of bins: {len(bins)-1} (width = {bin_width}μm)")

# ================================================================
# 1. 总体分箱图（每个基因单独画）
# ================================================================
print("\n" + "="*60)
print("Generating individual binned plots...")
print("="*60)

gene_colors = ['#E41A1C', '#4DAF4A', '#FFD700']

# 存储总体分箱相关性
overall_binned_correlations = []

for idx, gene in enumerate(genes_to_analyze):
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    
    expr_col = f'{gene}_expr'
    color = gene_colors[idx % len(gene_colors)]
    
    # 计算每个分箱的统计量
    bin_stats = filtered_df.groupby('distance_bin', observed=False).agg(
        mean_expr=(expr_col, 'mean'),
        sem_expr=(expr_col, 'sem'),
        count=(expr_col, 'count')
    ).reset_index()
    
    bin_stats['bin_center'] = bin_stats['distance_bin'].apply(lambda x: x.mid)
    bin_stats_filtered = bin_stats[bin_stats['count'] >= 3].copy()
    
    if len(bin_stats_filtered) >= 3:
        # 绘制误差棒
        ax.errorbar(bin_stats_filtered['bin_center'], 
                    bin_stats_filtered['mean_expr'],
                    yerr=bin_stats_filtered['sem_expr'], 
                    fmt='o-', 
                    color=color,
                    capsize=5, 
                    capthick=2, 
                    linewidth=2.5, 
                    markersize=8,
                    markerfacecolor=color,
                    markeredgecolor='white',
                    markeredgewidth=1.5,
                    label='Mean ± SEM',
                    zorder=3)
        
        # 计算分箱级别的相关性
        bin_r, bin_p = pearsonr(bin_stats_filtered['bin_center'], 
                                bin_stats_filtered['mean_expr'])
        
        overall_binned_correlations.append({
            'Gene': gene,
            'Binned_Pearson_r': bin_r,
            'Binned_Pearson_p': bin_p,
            'n_bins': len(bin_stats_filtered),
            'n_cells': len(filtered_df)
        })
        
        # 在标题下方显示统计信息
        ax.set_title(f'{gene}', fontsize=16, fontweight='bold', pad=30)
        ax.text(0.5, 0.92, f"r = {bin_r:.4f}, p = {bin_p:.2e}", 
                transform=ax.transAxes, fontsize=11, ha='center', va='bottom',
                fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='white', 
                         alpha=0.9, edgecolor='gray', linewidth=1))
    else:
        ax.set_title(f'{gene}', fontsize=16, fontweight='bold', pad=20)
        ax.text(0.5, 0.5, f'Insufficient data\n({len(bin_stats_filtered)} bins with ≥3 cells)',
                transform=ax.transAxes, ha='center', va='center', fontsize=12, color='red')
    
    ax.set_xlabel(f'Distance to nearest {source_celltype} (μm)', fontsize=13, fontweight='bold')
    ax.set_ylabel(f'Mean {gene} expression', fontsize=13, fontweight='bold')
    ax.set_xlim(0, max_distance)
    ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='lower right', fontsize=9, framealpha=0.9)
    
    plt.tight_layout()
    
    fig_path = os.path.join(output_dir, f"binned_{gene}_expression.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    print(f"  ✓ {gene} binned plot saved to: {fig_path}")

# 保存总体分箱相关性
if overall_binned_correlations:
    overall_corr_df = pd.DataFrame(overall_binned_correlations)
    overall_corr_df.to_csv(os.path.join(output_dir, "overall_binned_correlations.csv"), index=False)
    print(f"  ✓ Overall binned correlations saved to: {os.path.join(output_dir, 'overall_binned_correlations.csv')}")

# ================================================================
# 2. State分组的曲线图（每个基因一个图，多条曲线）
# ================================================================
if has_state:
    print("\n" + "="*60)
    print("Generating state-specific binned plots...")
    print("="*60)
    
    states = filtered_df['state'].unique()
    n_states = len(states)
    
    # 定义state颜色映射：N=绿色，A=蓝色，C=红色
    state_color_map = {
        'N': '#2E8B57',  # 海洋绿
        'A': '#1E88E5',  # 亮蓝
        'C': '#E53935'   # 亮红
    }
    
    # 如果还有其他state，使用备用颜色
    default_colors = ['#FF8F00', '#8E24AA', '#00897B', '#6D4C41', '#546E7A']
    state_markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
    
    # 为每个state分配颜色
    state_colors = {}
    state_marker_map = {}
    for i, state in enumerate(states):
        if state in state_color_map:
            state_colors[state] = state_color_map[state]
        else:
            # 如果state不在预定义中，使用备用颜色
            state_colors[state] = default_colors[i % len(default_colors)]
        state_marker_map[state] = state_markers[i % len(state_markers)]
    
    # 存储state分组的相关性
    state_binned_correlations = []
    
    for idx, gene in enumerate(genes_to_analyze):
        fig, ax = plt.subplots(1, 1, figsize=(11, 8))
        expr_col = f'{gene}_expr'
        
        state_stats_list = []
        state_corr_list = []
        
        # 先计算所有state的数据，确定y轴范围
        all_means = []
        
        for state in states:
            state_data = filtered_df[filtered_df['state'] == state]
            
            if len(state_data) < 10:
                print(f"  ⚠ State '{state}' has only {len(state_data)} cells for {gene}, skipping...")
                continue
            
            # 计算每个分箱的统计量
            state_data_copy = state_data.copy()
            state_data_copy['distance_bin'] = pd.cut(state_data_copy['distance'], bins=bins, right=False)
            bin_stats = state_data_copy.groupby('distance_bin', observed=False).agg(
                mean_expr=(expr_col, 'mean'),
                sem_expr=(expr_col, 'sem'),
                count=(expr_col, 'count')
            ).reset_index()
            
            bin_stats['bin_center'] = bin_stats['distance_bin'].apply(lambda x: x.mid)
            bin_stats_filtered = bin_stats[bin_stats['count'] >= 3].copy()
            
            if len(bin_stats_filtered) >= 3:
                all_means.extend(bin_stats_filtered['mean_expr'].values)
                
                color = state_colors[state]
                marker = state_marker_map[state]
                
                ax.errorbar(bin_stats_filtered['bin_center'], 
                            bin_stats_filtered['mean_expr'],
                            yerr=bin_stats_filtered['sem_expr'], 
                            fmt=f'{marker}-', 
                            color=color,
                            capsize=4, 
                            capthick=1.5, 
                            linewidth=2.5, 
                            markersize=9,
                            markerfacecolor=color,
                            markeredgecolor='white',
                            markeredgewidth=1.5,
                            label=f'{state}',
                            zorder=3,
                            elinewidth=1.5)
                
                # 计算该state的分箱相关性
                state_r, state_p = pearsonr(bin_stats_filtered['bin_center'], 
                                           bin_stats_filtered['mean_expr'])
                
                state_corr_list.append({
                    'Gene': gene,
                    'State': state,
                    'Binned_Pearson_r': state_r,
                    'Binned_Pearson_p': state_p,
                    'n_bins': len(bin_stats_filtered),
                    'n_cells': len(state_data)
                })
                
                state_stats_list.append({
                    'state': state,
                    'n_cells': len(state_data),
                    'n_bins': len(bin_stats_filtered),
                    'r': state_r,
                    'p': state_p,
                    'color': color
                })
        
        # 保存state相关性
        if state_corr_list:
            state_binned_correlations.extend(state_corr_list)
        
        if len(state_stats_list) >= 2:
            # 设置标题
            ax.set_title(f'{gene}', fontsize=17, fontweight='bold', pad=35)
            
            # 在图的中上部显示相关性信息（y轴0.82位置，水平居中）
            corr_lines = []
            for s in state_stats_list:
                corr_lines.append(f"{s['state']}: r = {s['r']:.4f}, p = {s['p']:.2e}")
            
            corr_text = "\n".join(corr_lines)
            
            # 放在中上部，不遮挡曲线
            if all_means:
                y_max = max(all_means) if all_means else 1
                y_range = y_max - (min(all_means) if all_means else 0)
                if y_range < 0.01:
                    y_range = 0.5
                y_pos = 0.82
            else:
                y_pos = 0.82
            
            props = dict(boxstyle='round', facecolor='white', alpha=0.92, 
                        edgecolor='gray', linewidth=1.5)
            
            ax.text(0.5, y_pos, corr_text, transform=ax.transAxes, 
                    fontsize=10, ha='center', va='center',
                    fontweight='bold',
                    bbox=props,
                    linespacing=1.6)
            
            ax.set_xlabel(f'Distance to nearest {source_celltype} (μm)', fontsize=13, fontweight='bold')
            ax.set_ylabel(f'Mean {gene} expression', fontsize=13, fontweight='bold')
            ax.set_xlim(0, max_distance)
            
            # 设置y轴范围，留出顶部空间给相关性信息
            if all_means:
                y_min = min(all_means) if all_means else 0
                y_max = max(all_means) if all_means else 1
                y_padding = (y_max - y_min) * 0.25 if y_max > y_min else 0.3
                ax.set_ylim(y_min - y_padding * 0.1, y_max + y_padding)
            
            ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            
            # 图例放在右上角
            ax.legend(loc='upper right', fontsize=10, framealpha=0.9, 
                     bbox_to_anchor=(1.0, 1.0))
        else:
            ax.text(0.5, 0.5, f'Insufficient data for state comparison\n(need at least 2 states with ≥3 bins)',
                    transform=ax.transAxes, ha='center', va='center', fontsize=12, color='red')
            ax.set_title(f'{gene} - Insufficient Data', fontsize=14, fontweight='bold')
            ax.set_xlabel(f'Distance to nearest {source_celltype} (μm)', fontsize=13, fontweight='bold')
            ax.set_ylabel(f'Mean {gene} expression', fontsize=13, fontweight='bold')
        
        plt.tight_layout()
        
        fig_path = os.path.join(output_dir, f"binned_{gene}_by_state.png")
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
        plt.close()
        print(f"  ✓ {gene} state-specific plot saved to: {fig_path}")
    
    # 保存state分组相关性
    if state_binned_correlations:
        state_corr_df = pd.DataFrame(state_binned_correlations)
        state_corr_df.to_csv(os.path.join(output_dir, "state_binned_correlations.csv"), index=False)
        print(f"  ✓ State binned correlations saved to: {os.path.join(output_dir, 'state_binned_correlations.csv')}")

# ================================================================
# 3. 原始数据的相关性（不分组）
# ================================================================
print("\n" + "="*60)
print("Calculating overall correlations...")
print("="*60)

overall_correlations = []
for gene in genes_to_analyze:
    expr_col = f'{gene}_expr'
    pearson_r, pearson_p = pearsonr(filtered_df['distance'], filtered_df[expr_col])
    spearman_r, spearman_p = spearmanr(filtered_df['distance'], filtered_df[expr_col])
    
    overall_correlations.append({
        'Gene': gene,
        'Pearson_r': pearson_r,
        'Pearson_p': pearson_p,
        'Spearman_r': spearman_r,
        'Spearman_p': spearman_p,
        'n_cells': len(filtered_df)
    })
    
    print(f"  {gene}:")
    print(f"    Pearson r = {pearson_r:.4f} (p = {pearson_p:.2e})")
    print(f"    Spearman r = {spearman_r:.4f} (p = {spearman_p:.2e})")

overall_corr_df = pd.DataFrame(overall_correlations)
overall_corr_df.to_csv(os.path.join(output_dir, "overall_correlations.csv"), index=False)

# ================================================================
# 4. State分组的原始数据相关性
# ================================================================
if has_state:
    print("\n" + "="*60)
    print("Calculating state-specific correlations...")
    print("="*60)
    
    state_correlations = []
    
    for state in filtered_df['state'].unique():
        state_data = filtered_df[filtered_df['state'] == state]
        
        if len(state_data) < 10:
            print(f"  ⚠ State '{state}' has only {len(state_data)} cells, skipping...")
            continue
        
        for gene in genes_to_analyze:
            expr_col = f'{gene}_expr'
            pearson_r, pearson_p = pearsonr(state_data['distance'], state_data[expr_col])
            spearman_r, spearman_p = spearmanr(state_data['distance'], state_data[expr_col])
            
            state_correlations.append({
                'State': state,
                'Gene': gene,
                'Pearson_r': pearson_r,
                'Pearson_p': pearson_p,
                'Spearman_r': spearman_r,
                'Spearman_p': spearman_p,
                'n_cells': len(state_data)
            })
    
    if state_correlations:
        state_corr_df = pd.DataFrame(state_correlations)
        state_corr_df.to_csv(os.path.join(output_dir, "state_correlations.csv"), index=False)
        print(f"  ✓ State correlations saved to: {os.path.join(output_dir, 'state_correlations.csv')}")
        
        # 打印结果
        print("\nState-specific correlations:")
        for state in state_corr_df['State'].unique():
            print(f"\n  State: {state}")
            state_subset = state_corr_df[state_corr_df['State'] == state]
            for _, row in state_subset.iterrows():
                print(f"    {row['Gene']}: r = {row['Pearson_r']:.4f} (p = {row['Pearson_p']:.2e})")

# ================================================================
# 5. 保存完整结果
# ================================================================
print("\n" + "="*60)
print("Saving complete results...")
print("="*60)

# 保存详细数据
results_df.to_csv(os.path.join(output_dir, "distance_expression_data.csv"), index=False)
print(f"  ✓ Detailed data saved to: {os.path.join(output_dir, 'distance_expression_data.csv')}")

# 保存过滤后的数据（含state信息）
filtered_df.to_csv(os.path.join(output_dir, "filtered_distance_expression_data.csv"), index=False)
print(f"  ✓ Filtered data saved to: {os.path.join(output_dir, 'filtered_distance_expression_data.csv')}")

# ---- 创建综合统计摘要 ----
summary_stats = []
for gene in genes_to_analyze:
    expr_col = f'{gene}_expr'
    pearson_r, pearson_p = pearsonr(filtered_df['distance'], filtered_df[expr_col])
    spearman_r, spearman_p = spearmanr(filtered_df['distance'], filtered_df[expr_col])
    
    # 查找分箱相关性
    bin_r = np.nan
    bin_p = np.nan
    for bc in overall_binned_correlations:
        if bc['Gene'] == gene:
            bin_r = bc['Binned_Pearson_r']
            bin_p = bc['Binned_Pearson_p']
            break
    
    summary_stats.append({
        'Gene': gene,
        'N_cells': len(filtered_df),
        'Pearson_r': pearson_r,
        'Pearson_p': pearson_p,
        'Spearman_r': spearman_r,
        'Spearman_p': spearman_p,
        'Binned_Pearson_r': bin_r,
        'Binned_Pearson_p': bin_p,
        'Mean_expression': filtered_df[expr_col].mean(),
        'Median_distance': filtered_df['distance'].median()
    })

summary_df = pd.DataFrame(summary_stats)
summary_df.to_csv(os.path.join(output_dir, "correlation_summary.csv"), index=False)
print(f"  ✓ Summary statistics saved to: {os.path.join(output_dir, 'correlation_summary.csv')}")

# ---- 打印最终总结 ----
print("\n" + "="*60)
print("KEY RESULTS SUMMARY")
print("="*60)
for _, row in summary_df.iterrows():
    print(f"\n{row['Gene']}:")
    print(f"  Cells analyzed: {int(row['N_cells']):,}")
    print(f"  Pearson r = {row['Pearson_r']:.4f} (p = {row['Pearson_p']:.2e})")
    print(f"  Spearman r = {row['Spearman_r']:.4f} (p = {row['Spearman_p']:.2e})")
    if not np.isnan(row['Binned_Pearson_r']):
        print(f"  Binned r = {row['Binned_Pearson_r']:.4f} (p = {row['Binned_Pearson_p']:.2e})")
    print(f"  Mean expression: {row['Mean_expression']:.4f}")

print(f"\n✅ All done! Results saved to: {output_dir}")
