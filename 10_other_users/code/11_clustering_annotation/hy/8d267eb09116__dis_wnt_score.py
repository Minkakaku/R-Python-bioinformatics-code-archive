import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy.sparse import issparse
from scipy.spatial import KDTree
from scipy.stats import pearsonr, spearmanr, kendalltau
import warnings
warnings.filterwarnings('ignore')

# ================= 配置 =================
input_file = "All_filter_state_samp_new_with_CNs10-15.h5ad"
output_dir = "wnt_score"  # 输出到wnt_score文件夹
genes_used = ["MYC", "AXIN2", "RNF43"]
target_celltype = "WNT+ ISC"
source_celltype = "matri-prod CAF"
max_distance = 250
bin_width = 25
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

# ---- 准备表达矩阵（使用原始counts） ----
if adata.raw is not None:
    print("\nUsing adata.raw for raw counts...")
    raw_adata = adata.raw.to_adata()
    raw_adata.obsm = adata.obsm.copy()
    raw_adata.obs = adata.obs.copy()
else:
    print("No raw found, using adata.X (assuming raw counts)...")
    raw_adata = adata.copy()

# ---- 检查基因存在性 ----
valid_genes = {}
for g in genes_used:
    if g in raw_adata.var_names:
        valid_genes[g] = g
        print(f"✓ Found '{g}'")
    else:
        match = [v for v in raw_adata.var_names if v.upper() == g.upper()]
        if match:
            valid_genes[g] = match[0]
            print(f"✓ Case-insensitive match: '{g}' -> '{match[0]}'")
        else:
            print(f"✗ Gene '{g}' not found!")
            valid_genes[g] = None

missing = [g for g in genes_used if valid_genes.get(g) is None]
if missing:
    raise SystemExit(f"Missing required genes: {missing}")

# ---- 初始化结果存储 ----
all_results = []

print(f"\nProcessing samples...")
samples = raw_adata.obs[sample_col].unique()
print(f"Found {len(samples)} sample(s)")

for sample in samples:
    print(f"\n--- Sample: {sample} ---")
    
    sample_mask = raw_adata.obs[sample_col] == sample
    sample_data = raw_adata[sample_mask].copy()
    
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
    
    # 获取目标细胞的原始counts
    target_raw = sample_data[target_mask]
    raw_counts = {}
    for g in genes_used:
        gene_name = valid_genes[g]
        X = target_raw[:, gene_name].X
        if issparse(X):
            raw_counts[g] = X.toarray().flatten()
        else:
            raw_counts[g] = np.asarray(X).flatten()
        # 确保非负
        raw_counts[g] = np.maximum(raw_counts[g], 0)
    
    # 获取state信息
    if has_state:
        target_states = sample_data[target_mask].obs[state_column].values
    
    for i, dist in enumerate(distances):
        # ---- 计算WNT突破评分 ----
        # WNT突破评分 = MYC / ((AXIN2 + RNF43) / 2)
        myc = raw_counts['MYC'][i]
        axin2 = raw_counts['AXIN2'][i]
        rnf43 = raw_counts['RNF43'][i]
        
        # 分母：AXIN2和RNF43的平均值
        denominator = (axin2 + rnf43) / 2
        
        # 计算WNT突破评分（避免除零）
        if denominator > 0:
            wnt_breakthrough_score = myc / denominator
        else:
            wnt_breakthrough_score = np.nan  # 如果分母为0，设为NaN
        
        result_dict = {
            'sample': sample,
            'distance': dist,
            'MYC_raw': myc,
            'AXIN2_raw': axin2,
            'RNF43_raw': rnf43,
            'WNT_breakthrough_score': wnt_breakthrough_score  # MYC / ((AXIN2 + RNF43) / 2)
        }
        
        if has_state:
            result_dict['state'] = target_states[i]
        all_results.append(result_dict)

if len(all_results) == 0:
    raise RuntimeError("No valid cell pairs found in any sample!")

results_df = pd.DataFrame(all_results)
print(f"\nTotal analyzed cells: {len(results_df)}")

# 移除NaN值（分母为0的细胞）
valid_mask = ~results_df['WNT_breakthrough_score'].isna()
results_df = results_df[valid_mask].copy()
print(f"Cells with valid WNT breakthrough score: {len(results_df)} (removed {valid_mask.sum() - len(results_df)} cells with denominator=0)")

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
# 1. 总体分箱图
# ================================================================
print("\n" + "="*60)
print("Generating binned plots for WNT breakthrough score...")
print("="*60)

# ---- 1.1 WNT突破评分分箱图 ----
fig, ax = plt.subplots(1, 1, figsize=(10, 7))

expr_col = 'WNT_breakthrough_score'
color = '#E74C3C'  # 红色

# 计算每个分箱的统计量
bin_stats = filtered_df.groupby('distance_bin', observed=False).agg(
    mean_expr=(expr_col, 'mean'),
    sem_expr=(expr_col, 'sem'),
    count=(expr_col, 'count')
).reset_index()

bin_stats['bin_center'] = bin_stats['distance_bin'].apply(lambda x: x.mid)
bin_stats_filtered = bin_stats[bin_stats['count'] >= 3].copy()

if len(bin_stats_filtered) >= 3:
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
    spearman_r, spearman_p = spearmanr(bin_stats_filtered['bin_center'], 
                                       bin_stats_filtered['mean_expr'])
    kendall_tau, kendall_p = kendalltau(bin_stats_filtered['bin_center'], 
                                        bin_stats_filtered['mean_expr'])
    
    ax.set_title('WNT Breakthrough Score', fontsize=16, fontweight='bold', pad=30)
    ax.text(0.5, 0.92, 
            f"Pearson r = {bin_r:.4f}, p = {bin_p:.2e}\n"
            f"Spearman r = {spearman_r:.4f}, p = {spearman_p:.2e}\n"
            f"Kendall τ = {kendall_tau:.4f}, p = {kendall_p:.2e}", 
            transform=ax.transAxes, fontsize=10, ha='center', va='bottom',
            fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', 
                     alpha=0.9, edgecolor='gray', linewidth=1))
    
    # 保存分箱相关性
    binned_corr_data = {
        'Feature': 'WNT Breakthrough Score',
        'Binned_Pearson_r': bin_r,
        'Binned_Pearson_p': bin_p,
        'Binned_Spearman_r': spearman_r,
        'Binned_Spearman_p': spearman_p,
        'Binned_Kendall_tau': kendall_tau,
        'Binned_Kendall_p': kendall_p,
        'n_bins': len(bin_stats_filtered),
        'n_cells': len(filtered_df)
    }
    binned_corr_df = pd.DataFrame([binned_corr_data])
    binned_corr_df.to_csv(os.path.join(output_dir, "binned_correlations.csv"), index=False)
else:
    ax.set_title('WNT Breakthrough Score', fontsize=16, fontweight='bold', pad=20)
    ax.text(0.5, 0.5, f'Insufficient data\n({len(bin_stats_filtered)} bins with ≥3 cells)',
            transform=ax.transAxes, ha='center', va='center', fontsize=12, color='red')

ax.set_xlabel(f'Distance to nearest {source_celltype} (μm)', fontsize=13, fontweight='bold')
ax.set_ylabel('Mean WNT Breakthrough Score\n(MYC / ((AXIN2 + RNF43) / 2))', fontsize=13, fontweight='bold')
ax.set_xlim(0, max_distance)
ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(loc='lower right', fontsize=9, framealpha=0.9)

plt.tight_layout()

fig_path = os.path.join(output_dir, "binned_WNT_breakthrough_score.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
plt.close()
print(f"  ✓ WNT breakthrough score binned plot saved to: {fig_path}")

# ---- 1.2 单个基因的分箱图（用于参考） ----
gene_colors = ['#E41A1C', '#4DAF4A', '#FFD700']
gene_labels = ['MYC', 'AXIN2', 'RNF43']

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for idx, gene in enumerate(['MYC', 'AXIN2', 'RNF43']):
    ax = axes[idx]
    expr_col = f'{gene}_raw'
    color = gene_colors[idx]
    
    bin_stats = filtered_df.groupby('distance_bin', observed=False).agg(
        mean_expr=(expr_col, 'mean'),
        sem_expr=(expr_col, 'sem'),
        count=(expr_col, 'count')
    ).reset_index()
    
    bin_stats['bin_center'] = bin_stats['distance_bin'].apply(lambda x: x.mid)
    bin_stats_filtered = bin_stats[bin_stats['count'] >= 3].copy()
    
    if len(bin_stats_filtered) >= 3:
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
        
        bin_r, bin_p = pearsonr(bin_stats_filtered['bin_center'], 
                                bin_stats_filtered['mean_expr'])
        
        ax.set_title(f'{gene} (raw counts)', fontsize=14, fontweight='bold', pad=25)
        ax.text(0.5, 0.92, f"r = {bin_r:.4f}, p = {bin_p:.2e}", 
                transform=ax.transAxes, fontsize=10, ha='center', va='bottom',
                fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='white', 
                         alpha=0.9, edgecolor='gray', linewidth=1))
    else:
        ax.set_title(f'{gene} (raw counts)', fontsize=14, fontweight='bold', pad=20)
        ax.text(0.5, 0.5, f'Insufficient data',
                transform=ax.transAxes, ha='center', va='center', fontsize=12, color='red')
    
    ax.set_xlabel(f'Distance to nearest {source_celltype} (μm)', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'Mean {gene} raw counts', fontsize=12, fontweight='bold')
    ax.set_xlim(0, max_distance)
    ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='lower right', fontsize=9, framealpha=0.9)

plt.tight_layout()
fig_path = os.path.join(output_dir, "binned_individual_genes.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
plt.close()
print(f"  ✓ Individual genes binned plot saved to: {fig_path}")

# ================================================================
# 2. State分组的曲线图
# ================================================================
if has_state:
    print("\n" + "="*60)
    print("Generating state-specific plots...")
    print("="*60)
    
    states = filtered_df['state'].unique()
    
    # 定义state颜色映射
    state_color_map = {
        'N': '#2E8B57',
        'A': '#1E88E5',
        'C': '#E53935'
    }
    default_colors = ['#FF8F00', '#8E24AA', '#00897B', '#6D4C41', '#546E7A']
    state_markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h']
    
    state_colors = {}
    state_marker_map = {}
    for i, state in enumerate(states):
        if state in state_color_map:
            state_colors[state] = state_color_map[state]
        else:
            state_colors[state] = default_colors[i % len(default_colors)]
        state_marker_map[state] = state_markers[i % len(state_markers)]
    
    # ---- 2.1 WNT突破评分 ----
    fig, ax = plt.subplots(1, 1, figsize=(11, 8))
    expr_col = 'WNT_breakthrough_score'
    
    state_stats_list = []
    all_means = []
    
    for state in states:
        state_data = filtered_df[filtered_df['state'] == state]
        if len(state_data) < 10:
            print(f"  ⚠ State '{state}' has only {len(state_data)} cells, skipping...")
            continue
        
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
            
            state_r, state_p = pearsonr(bin_stats_filtered['bin_center'], 
                                       bin_stats_filtered['mean_expr'])
            spearman_r, spearman_p = spearmanr(bin_stats_filtered['bin_center'], 
                                               bin_stats_filtered['mean_expr'])
            kendall_tau, kendall_p = kendalltau(bin_stats_filtered['bin_center'], 
                                                bin_stats_filtered['mean_expr'])
            
            state_stats_list.append({
                'state': state,
                'pearson_r': state_r,
                'pearson_p': state_p,
                'spearman_r': spearman_r,
                'spearman_p': spearman_p,
                'kendall_tau': kendall_tau,
                'kendall_p': kendall_p,
                'color': color
            })
    
    if len(state_stats_list) >= 2:
        ax.set_title('WNT Breakthrough Score by State', fontsize=17, fontweight='bold', pad=35)
        
        # 显示统计信息
        corr_lines = []
        for s in state_stats_list:
            corr_lines.append(f"{s['state']}: r = {s['pearson_r']:.4f}, p = {s['pearson_p']:.2e}")
        corr_text = "\n".join(corr_lines)
        
        kendall_lines = [f"{s['state']}: τ = {s['kendall_tau']:.4f}, p = {s['kendall_p']:.2e}" for s in state_stats_list]
        kendall_text = "\n".join(kendall_lines)
        
        full_text = f"Pearson:\n{corr_text}\n\nMann-Kendall (τ):\n{kendall_text}"
        
        y_pos = 0.82
        props = dict(boxstyle='round', facecolor='white', alpha=0.92, 
                    edgecolor='gray', linewidth=1.5)
        ax.text(0.5, y_pos, full_text, transform=ax.transAxes, 
                fontsize=9, ha='center', va='center',
                fontweight='bold', bbox=props, linespacing=1.5)
        
        ax.set_xlabel(f'Distance to nearest {source_celltype} (μm)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Mean WNT Breakthrough Score\n(MYC / ((AXIN2 + RNF43) / 2))', fontsize=13, fontweight='bold')
        ax.set_xlim(0, max_distance)
        
        if all_means:
            y_min = min(all_means) if all_means else 0
            y_max = max(all_means) if all_means else 1
            y_padding = (y_max - y_min) * 0.25 if y_max > y_min else 0.3
            ax.set_ylim(y_min - y_padding * 0.1, y_max + y_padding)
        
        ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.legend(loc='upper right', fontsize=10, framealpha=0.9, bbox_to_anchor=(1.0, 1.0))
    else:
        ax.text(0.5, 0.5, f'Insufficient data for state comparison',
                transform=ax.transAxes, ha='center', va='center', fontsize=12, color='red')
        ax.set_title('WNT Breakthrough Score - Insufficient Data', fontsize=14, fontweight='bold')
        ax.set_xlabel(f'Distance to nearest {source_celltype} (μm)', fontsize=13, fontweight='bold')
        ax.set_ylabel('Mean WNT Breakthrough Score', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    fig_path = os.path.join(output_dir, "binned_WNT_breakthrough_score_by_state.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()
    print(f"  ✓ WNT breakthrough score state-specific plot saved to: {fig_path}")

# ================================================================
# 3. 整体相关性分析
# ================================================================
print("\n" + "="*60)
print("Calculating overall correlations...")
print("="*60)

# WNT突破评分的相关性
score = filtered_df['WNT_breakthrough_score'].values
distances = filtered_df['distance'].values

pearson_r, pearson_p = pearsonr(distances, score)
spearman_r, spearman_p = spearmanr(distances, score)
kendall_tau, kendall_p = kendalltau(distances, score)

print(f"\nWNT Breakthrough Score:")
print(f"  Pearson r = {pearson_r:.4f} (p = {pearson_p:.2e})")
print(f"  Spearman r = {spearman_r:.4f} (p = {spearman_p:.2e})")
print(f"  Kendall τ = {kendall_tau:.4f} (p = {kendall_p:.2e})")

# 单个基因的相关性
gene_correlations = []
for gene in ['MYC', 'AXIN2', 'RNF43']:
    expr = filtered_df[f'{gene}_raw'].values
    p_r, p_p = pearsonr(distances, expr)
    s_r, s_p = spearmanr(distances, expr)
    k_t, k_p = kendalltau(distances, expr)
    
    gene_correlations.append({
        'Gene': gene,
        'Pearson_r': p_r,
        'Pearson_p': p_p,
        'Spearman_r': s_r,
        'Spearman_p': s_p,
        'Kendall_tau': k_t,
        'Kendall_p': k_p
    })
    print(f"\n{gene} (raw counts):")
    print(f"  Pearson r = {p_r:.4f} (p = {p_p:.2e})")
    print(f"  Spearman r = {s_r:.4f} (p = {s_p:.2e})")
    print(f"  Kendall τ = {k_t:.4f} (p = {k_p:.2e})")

# ================================================================
# 4. 保存结果
# ================================================================
print("\n" + "="*60)
print("Saving results...")
print("="*60)

# 保存详细数据
results_df.to_csv(os.path.join(output_dir, "wnt_breakthrough_data.csv"), index=False)
print(f"  ✓ Detailed data saved to: {os.path.join(output_dir, 'wnt_breakthrough_data.csv')}")

filtered_df.to_csv(os.path.join(output_dir, "filtered_wnt_breakthrough_data.csv"), index=False)
print(f"  ✓ Filtered data saved to: {os.path.join(output_dir, 'filtered_wnt_breakthrough_data.csv')}")

# 保存相关性摘要
summary_data = {
    'Feature': ['WNT Breakthrough Score'] + [f'{g} (raw)' for g in ['MYC', 'AXIN2', 'RNF43']],
    'Pearson_r': [pearson_r] + [c['Pearson_r'] for c in gene_correlations],
    'Pearson_p': [pearson_p] + [c['Pearson_p'] for c in gene_correlations],
    'Spearman_r': [spearman_r] + [c['Spearman_r'] for c in gene_correlations],
    'Spearman_p': [spearman_p] + [c['Spearman_p'] for c in gene_correlations],
    'Kendall_tau': [kendall_tau] + [c['Kendall_tau'] for c in gene_correlations],
    'Kendall_p': [kendall_p] + [c['Kendall_p'] for c in gene_correlations],
    'N_cells': [len(filtered_df)] * 4
}

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(os.path.join(output_dir, "correlation_summary.csv"), index=False)
print(f"  ✓ Correlation summary saved to: {os.path.join(output_dir, 'correlation_summary.csv')}")

# ---- 打印最终总结 ----
print("\n" + "="*60)
print("KEY RESULTS SUMMARY - WNT Breakthrough Score")
print("="*60)
print(f"\nWNT Breakthrough Score = MYC / ((AXIN2 + RNF43) / 2)")
print(f"\nCells analyzed: {len(filtered_df):,}")
print(f"\nOverall correlations with distance to {source_celltype}:")
print(f"  Pearson r = {pearson_r:.4f} (p = {pearson_p:.2e})")
print(f"  Spearman r = {spearman_r:.4f} (p = {spearman_p:.2e})")
print(f"  Kendall τ = {kendall_tau:.4f} (p = {kendall_p:.2e})")

print(f"\n✅ All done! Results saved to: {output_dir}")
