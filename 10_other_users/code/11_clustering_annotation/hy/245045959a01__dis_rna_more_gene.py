import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import os
from scipy.sparse import issparse
from scipy.spatial import KDTree
from scipy.stats import spearmanr, kruskal, kendalltau
from scipy.ndimage import gaussian_filter1d
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# ================= 字体设置 =================
# 如果字体路径不对，请注释掉下面这行，改用系统默认
font_manager.fontManager.addfont('/home/hy/software/miniconda3/envs/sprna/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf')

# ================== 全局设置 ==================
plt.rcParams['font.family'] = 'arial'
plt.rcParams['font.weight'] = 'heavy'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ================= 配置 =================
input_file = "All_filter_state_samp_new_with_CNs10-15.h5ad"
output_dir = "distance_expression_analysis_more_genes"
genes_to_analyze = ["CCND1", "CCND2", "CDK4", "LGR5", "SOX9", "TCF7", "ASCL2", "NKD1", "DKK1", "SP5", "c-JUN", "VEGFA", "MMP7", "FOSL1"]
target_celltype = "WNT+ ISC"
source_celltype = "matri-prod CAF"
max_distance = 250
bin_width = 25
state_column = 'state'
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
    raw_counts_adata = expr_adata.copy()
else:
    print("No raw found, using normalized data...")
    expr_adata = adata.copy()
    sc.pp.normalize_total(expr_adata, target_sum=1e4)
    sc.pp.log1p(expr_adata)
    raw_counts_adata = adata.copy()

# ---- 检查基因存在性（容错处理） ----
valid_genes = {}
found_any = False
for g in genes_to_analyze:
    # 1. 精确匹配
    if g in expr_adata.var_names:
        valid_genes[g] = g
        print(f"✓ Found '{g}'")
        found_any = True
        continue
    
    # 2. 大小写+符号容错匹配（去除 - 和 _ 后比较）
    g_clean = g.replace('-', '').replace('_', '').upper()
    match = None
    for v in expr_adata.var_names:
        v_clean = v.replace('-', '').replace('_', '').upper()
        if v_clean == g_clean:
            match = v
            break
    
    if match:
        valid_genes[g] = match
        print(f"✓ Fuzzy match: '{g}' -> '{match}'")
        found_any = True
    else:
        print(f"✗ Gene '{g}' not found in data, will skip this gene.")
        valid_genes[g] = None

# 过滤掉未找到的基因，只保留有效基因
final_genes = []
for g in genes_to_analyze:
    if valid_genes.get(g) is not None:
        final_genes.append(g)

if not final_genes:
    raise SystemExit("No valid genes found in the data! Check gene names or your dataset.")

print(f"\nWill analyze {len(final_genes)} genes: {final_genes}")

# ---- 初始化结果存储 ----
all_results = []

print(f"\nProcessing samples...")
samples = expr_adata.obs[sample_col].unique()
print(f"Found {len(samples)} sample(s)")

for sample in samples:
    print(f"\n--- Sample: {sample} ---")
    
    sample_mask = expr_adata.obs[sample_col] == sample
    sample_data = expr_adata[sample_mask].copy()
    raw_sample_data = raw_counts_adata[sample_mask].copy()
    
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
    
    # 获取目标细胞的基因表达
    target_expr = sample_data[target_mask]
    gene_exprs = {}
    for g in final_genes:
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
        for g in final_genes:
            result_dict[f'{g}_expr'] = gene_exprs[g][i]
        
        if has_state:
            result_dict['state'] = target_states[i]
        all_results.append(result_dict)

if len(all_results) == 0:
    raise RuntimeError("No valid cell pairs found in any sample!")

results_df = pd.DataFrame(all_results)
print(f"\nTotal analyzed cells: {len(results_df)}")

# ---- 过滤距离范围 ----
filtered_df = results_df[results_df['distance'] <= max_distance].copy()
print(f"\nCells within {max_distance} units: {len(filtered_df)} ({len(filtered_df)/len(results_df)*100:.1f}%)")

if len(filtered_df) < 10:
    print("⚠ Warning: Very few cells within the specified distance range!")

# ---- 创建距离分箱 ----
bins = np.arange(0, max_distance + bin_width, bin_width)
filtered_df['distance_bin'] = pd.cut(filtered_df['distance'], bins=bins, right=False)

print(f"\nNumber of bins: {len(bins)-1} (width = {bin_width}μm)")

# ================================================================
# 1. 总体分箱图（仅单基因）—— 修改图例位置
# ================================================================
print("\n" + "="*60)
print("Generating individual binned plots...")
print("="*60)

# 为最终基因动态分配颜色
gene_colors = ['#E41A1C', '#4DAF4A', '#FFD700', '#9B59B6', '#FF6B6B', '#1E88E5', '#8E24AA', '#00897B', '#FF8F00', '#6D4C41', '#D32F2F', '#388E3C', '#1976D2', '#F57C00']
overall_binned_correlations = []

for idx, gene in enumerate(final_genes):
    fig, ax = plt.subplots(1, 1, figsize=(16, 10))
    
    expr_col = f'{gene}_expr'
    color = gene_colors[idx % len(gene_colors)]
    
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
                    capsize=12, 
                    capthick=6, 
                    linewidth=9,
                    markersize=32,
                    markerfacecolor=color,
                    markeredgecolor='white',
                    markeredgewidth=6,
                    label='Mean ± SEM',
                    zorder=3,
                    elinewidth=6)
        
        kendall_tau, kendall_p = kendalltau(bin_stats_filtered['bin_center'], 
                                            bin_stats_filtered['mean_expr'])
        
        overall_binned_correlations.append({
            'Feature': gene,
            'Type': 'Single Gene',
            'Binned_Kendall_tau': kendall_tau,
            'Binned_Kendall_p': kendall_p,
            'n_bins': len(bin_stats_filtered),
            'n_cells': len(filtered_df)
        })
        
        ax.set_title(f'{gene}', fontsize=62, fontweight='heavy', pad=50)
        
        # ---------- 统计文本 ----------
        ax.text(1.02, 0.95, f"Mann-Kendall (τ)\nτ = {kendall_tau:.4f}\np = {kendall_p:.2e}", 
                transform=ax.transAxes, fontsize=42, ha='left', va='top',
                fontweight='heavy', linespacing=2.2)
        
        # ---------- 将图例标签移到统计文本下方 ----------
        # 删除原有的图例（注释掉 ax.legend）
        # 手动添加一个带颜色圆点的说明
        ax.text(1.02, 0.5, f"●  Mean ± SEM", transform=ax.transAxes, fontsize=42,
                fontweight='heavy', color=color, va='top')
    else:
        ax.set_title(f'{gene}', fontsize=62, fontweight='heavy', pad=20)
        ax.text(0.5, 0.5, f'Insufficient data\n({len(bin_stats_filtered)} bins with ≥3 cells)',
                transform=ax.transAxes, ha='center', va='center', fontsize=42, color='red',
                fontweight='heavy')
    
    ax.set_xlabel('Distance', fontsize=55, fontweight='heavy')
    ax.set_ylabel('Expression', fontsize=55, fontweight='heavy')
    ax.set_xlim(0, max_distance)
    ax.tick_params(labelsize=46, width=5, length=14)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('heavy')
    ax.grid(True, alpha=0.25, linestyle='--', linewidth=2.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(6)
    ax.spines['bottom'].set_linewidth(6)
    
    # 不再使用 ax.legend()
    # ax.legend(loc='lower right', ...)   # 已移除
    
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
# 2. State分组的曲线图（仅单基因）—— 保持不变
# ================================================================
if has_state:
    print("\n" + "="*60)
    print("Generating state-specific binned plots...")
    print("="*60)
    
    states = filtered_df['state'].unique()
    n_states = len(states)
    
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
    
    state_binned_correlations = []
    
    for idx, gene in enumerate(final_genes):
        fig, ax = plt.subplots(1, 1, figsize=(22, 12))
        expr_col = f'{gene}_expr'
        
        state_stats_list = []
        state_corr_list = []
        all_means = []
        
        for state in states:
            state_data = filtered_df[filtered_df['state'] == state]
            if len(state_data) < 10:
                print(f"  ⚠ State '{state}' has only {len(state_data)} cells for {gene}, skipping...")
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
                            capsize=14, 
                            capthick=6, 
                            linewidth=11,
                            markersize=40,
                            markerfacecolor=color,
                            markeredgecolor='white',
                            markeredgewidth=6,
                            label=f'{state}',
                            zorder=3,
                            elinewidth=7)
                
                kendall_tau, kendall_p = kendalltau(bin_stats_filtered['bin_center'], 
                                                    bin_stats_filtered['mean_expr'])
                
                state_corr_list.append({
                    'Feature': gene,
                    'Type': 'Single Gene',
                    'State': state,
                    'Binned_Kendall_tau': kendall_tau,
                    'Binned_Kendall_p': kendall_p,
                    'n_bins': len(bin_stats_filtered),
                    'n_cells': len(state_data)
                })
                
                state_stats_list.append({
                    'state': state,
                    'tau': kendall_tau,
                    'tau_p': kendall_p,
                    'color': color
                })
        
        if state_corr_list:
            state_binned_correlations.extend(state_corr_list)
        
        if len(state_stats_list) >= 2:
            ax.set_title(f'{gene}', fontsize=65, fontweight='heavy', pad=60)
            
            kendall_lines = [f"{s['state']}: τ = {s['tau']:.4f}, p = {s['tau_p']:.2e}" for s in state_stats_list]
            kendall_text = "Mann-Kendall (τ):\n" + "\n".join(kendall_lines)
            
            ax.text(1.02, 0.95, kendall_text, transform=ax.transAxes, 
                    fontsize=48, ha='left', va='top',
                    fontweight='heavy', linespacing=2.5)
            
            ax.set_xlabel('Distance', fontsize=60, fontweight='heavy')
            ax.set_ylabel('Expression', fontsize=60, fontweight='heavy')
            ax.set_xlim(0, max_distance)
            ax.tick_params(labelsize=50, width=6, length=16)
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontweight('heavy')
            
            if all_means:
                y_min = min(all_means) if all_means else 0
                y_max = max(all_means) if all_means else 1
                y_padding = (y_max - y_min) * 0.25 if y_max > y_min else 0.3
                ax.set_ylim(y_min - y_padding * 0.1, y_max + y_padding)
            
            ax.grid(True, alpha=0.25, linestyle='--', linewidth=2.5)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(6)
            ax.spines['bottom'].set_linewidth(6)
            
            # State图图例保留（显示各State名称）
            ax.legend(loc='center left', fontsize=48, framealpha=0.9,
                     bbox_to_anchor=(1.02, 0.02), markerscale=2.5,
                     prop={'weight': 'heavy', 'size': 48},
                     handletextpad=0.5, labelspacing=0.8)
        else:
            ax.text(0.5, 0.5, f'Insufficient data for state comparison',
                    transform=ax.transAxes, ha='center', va='center', fontsize=48, color='red',
                    fontweight='heavy')
            ax.set_title(f'{gene} - Insufficient Data', fontsize=55, fontweight='heavy')
            ax.set_xlabel('Distance', fontsize=60, fontweight='heavy')
            ax.set_ylabel('Expression', fontsize=60, fontweight='heavy')
            ax.tick_params(labelsize=50)
        
        plt.tight_layout()
        fig_path = os.path.join(output_dir, f"binned_{gene}_by_state.png")
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.savefig(fig_path.replace('.png', '.pdf'), bbox_inches='tight')
        plt.close()
        print(f"  ✓ {gene} state-specific plot saved to: {fig_path}")
    
    if state_binned_correlations:
        state_corr_df = pd.DataFrame(state_binned_correlations)
        state_corr_df.to_csv(os.path.join(output_dir, "state_binned_correlations.csv"), index=False)
        print(f"  ✓ State binned correlations saved to: {os.path.join(output_dir, 'state_binned_correlations.csv')}")

# ================================================================
# 3. 原始数据的相关性（Spearman和Kendall）—— 不变
# ================================================================
print("\n" + "="*60)
print("Calculating overall correlations (Spearman & Kendall)...")
print("="*60)

overall_correlations = []

all_features = []
for gene in final_genes:
    all_features.append({'name': gene, 'col': f'{gene}_expr', 'type': 'Single Gene'})

for feature in all_features:
    expr_col = feature['col']
    feature_name = feature['name']
    feature_type = feature['type']
    
    spearman_r, spearman_p = spearmanr(filtered_df['distance'], filtered_df[expr_col])
    kendall_tau, kendall_p = kendalltau(filtered_df['distance'], filtered_df[expr_col])
    
    overall_correlations.append({
        'Feature': feature_name,
        'Type': feature_type,
        'Spearman_r': spearman_r,
        'Spearman_p': spearman_p,
        'Kendall_tau': kendall_tau,
        'Kendall_p': kendall_p,
        'n_cells': len(filtered_df)
    })
    
    print(f"  {feature_name} ({feature_type}):")
    print(f"    Spearman r = {spearman_r:.4f} (p = {spearman_p:.2e})")
    print(f"    Kendall τ = {kendall_tau:.4f} (p = {kendall_p:.2e})")

overall_corr_df = pd.DataFrame(overall_correlations)
overall_corr_df.to_csv(os.path.join(output_dir, "overall_correlations.csv"), index=False)
print(f"  ✓ Overall correlations saved to: {os.path.join(output_dir, 'overall_correlations.csv')}")

# ================================================================
# 4. State分组的原始数据相关性—— 不变
# ================================================================
if has_state:
    print("\n" + "="*60)
    print("Calculating state-specific correlations (Spearman & Kendall)...")
    print("="*60)
    
    state_correlations = []
    
    for state in filtered_df['state'].unique():
        state_data = filtered_df[filtered_df['state'] == state]
        
        if len(state_data) < 10:
            print(f"  ⚠ State '{state}' has only {len(state_data)} cells, skipping...")
            continue
        
        for feature in all_features:
            expr_col = feature['col']
            feature_name = feature['name']
            feature_type = feature['type']
            
            spearman_r, spearman_p = spearmanr(state_data['distance'], state_data[expr_col])
            kendall_tau, kendall_p = kendalltau(state_data['distance'], state_data[expr_col])
            
            state_correlations.append({
                'State': state,
                'Feature': feature_name,
                'Type': feature_type,
                'Spearman_r': spearman_r,
                'Spearman_p': spearman_p,
                'Kendall_tau': kendall_tau,
                'Kendall_p': kendall_p,
                'n_cells': len(state_data)
            })
    
    if state_correlations:
        state_corr_df = pd.DataFrame(state_correlations)
        state_corr_df.to_csv(os.path.join(output_dir, "state_correlations.csv"), index=False)
        print(f"  ✓ State correlations saved to: {os.path.join(output_dir, 'state_correlations.csv')}")
        
        print("\nState-specific correlations (Kendall τ):")
        for state in state_corr_df['State'].unique():
            print(f"\n  State: {state}")
            state_subset = state_corr_df[state_corr_df['State'] == state]
            for _, row in state_subset.iterrows():
                print(f"    {row['Feature']}: τ = {row['Kendall_tau']:.4f} (p = {row['Kendall_p']:.2e})")

# ================================================================
# 5. 保存完整结果—— 不变
# ================================================================
print("\n" + "="*60)
print("Saving complete results...")
print("="*60)

results_df.to_csv(os.path.join(output_dir, "distance_expression_data.csv"), index=False)
print(f"  ✓ Detailed data saved to: {os.path.join(output_dir, 'distance_expression_data.csv')}")

filtered_df.to_csv(os.path.join(output_dir, "filtered_distance_expression_data.csv"), index=False)
print(f"  ✓ Filtered data saved to: {os.path.join(output_dir, 'filtered_distance_expression_data.csv')}")

# 创建综合统计摘要
summary_stats = []
for feature in all_features:
    feature_name = feature['name']
    feature_type = feature['type']
    expr_col = feature['col']
    
    spearman_r, spearman_p = spearmanr(filtered_df['distance'], filtered_df[expr_col])
    kendall_tau, kendall_p = kendalltau(filtered_df['distance'], filtered_df[expr_col])
    
    bin_tau = np.nan
    bin_tau_p = np.nan
    for bc in overall_binned_correlations:
        if bc['Feature'] == feature_name:
            bin_tau = bc.get('Binned_Kendall_tau', np.nan)
            bin_tau_p = bc.get('Binned_Kendall_p', np.nan)
            break
    
    summary_stats.append({
        'Feature': feature_name,
        'Type': feature_type,
        'N_cells': len(filtered_df),
        'Spearman_r': spearman_r,
        'Spearman_p': spearman_p,
        'Kendall_tau': kendall_tau,
        'Kendall_p': kendall_p,
        'Binned_Kendall_tau': bin_tau,
        'Binned_Kendall_p': bin_tau_p,
        'Mean_expression': filtered_df[expr_col].mean(),
        'Median_distance': filtered_df['distance'].median()
    })

summary_df = pd.DataFrame(summary_stats)
summary_df.to_csv(os.path.join(output_dir, "correlation_summary.csv"), index=False)
print(f"  ✓ Summary statistics saved to: {os.path.join(output_dir, 'correlation_summary.csv')}")

# ---- 打印最终总结 ----
print("\n" + "="*60)
print("KEY RESULTS SUMMARY (Mann-Kendall trend test)")
print("="*60)
for _, row in summary_df.iterrows():
    print(f"\n{row['Feature']} ({row['Type']}):")
    print(f"  Cells analyzed: {int(row['N_cells']):,}")
    print(f"  Spearman r = {row['Spearman_r']:.4f} (p = {row['Spearman_p']:.2e})")
    print(f"  Kendall τ = {row['Kendall_tau']:.4f} (p = {row['Kendall_p']:.2e})")
    if not np.isnan(row['Binned_Kendall_tau']):
        print(f"  Binned Kendall τ = {row['Binned_Kendall_tau']:.4f} (p = {row['Binned_Kendall_p']:.2e})")
    print(f"  Mean expression: {row['Mean_expression']:.4f}")

print(f"\n✅ All done! Results saved to: {output_dir}")
