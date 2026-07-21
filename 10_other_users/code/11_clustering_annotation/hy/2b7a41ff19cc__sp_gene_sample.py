import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.sparse import issparse
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings('ignore')

# ================= 配置 =================
input_file = "All_filter_state_samp_new_with_CNs10-15.h5ad"
output_dir = "gene_spatial_plots_by_sample"
genes_to_plot = ["CTNNB1", "MYC", "AXIN2", "RNF43"]
# ========================================

os.makedirs(output_dir, exist_ok=True)

print(f"Reading {input_file}...")
adata = sc.read_h5ad(input_file)
print(f"Data shape: {adata.shape}")

# 使用 raw 或归一化后的表达矩阵
if adata.raw is not None:
    print("Using adata.raw as expression matrix...")
    expr_adata = adata.raw.to_adata()
    expr_adata.obsm = adata.obsm.copy()
    expr_adata.obs = adata.obs.copy()
else:
    print("No raw found, normalizing + log1p...")
    expr_adata = adata.copy()
    sc.pp.normalize_total(expr_adata, target_sum=1e4)
    sc.pp.log1p(expr_adata)

# 基因存在性检查
valid_genes = []
for g in genes_to_plot:
    if g in expr_adata.var_names:
        valid_genes.append(g)
        print(f"✓ Found '{g}'")
    else:
        match = [v for v in expr_adata.var_names if v.upper() == g.upper()]
        if match:
            valid_genes.append(match[0])
            print(f"✓ Case-insensitive match: '{g}' -> '{match[0]}'")
        else:
            print(f"✗ Gene '{g}' not found!")
if not valid_genes:
    raise SystemExit("No target genes found.")

# 自动识别样本列
sample_col = None
for col in ['sample', 'Sample', 'patient', 'Patient', 'batch', 'Batch', 'orig.ident', 'slice', 'Slide']:
    if col in expr_adata.obs.columns:
        sample_col = col
        break
if sample_col is None:
    print("No sample column found, treating all cells as one sample.")
    expr_adata.obs['sample'] = 'All_Samples'
    sample_col = 'sample'
print(f"Sample column detected: '{sample_col}'")

# 自动识别空间坐标
spatial_key = None
for key in ['spatial', 'X_spatial', 'coordinates', 'obsm/spatial', 'X_umap']:
    if key in expr_adata.obsm:
        spatial_key = key
        break
if spatial_key is None:
    for key in expr_adata.obsm.keys():
        if expr_adata.obsm[key].shape[1] >= 2:
            spatial_key = key
            break
if spatial_key is None:
    raise RuntimeError("No 2D spatial coordinates found in adata.obsm.")
print(f"Spatial coordinates: '{spatial_key}'")

# 获取所有实际存在的样本名
samples_in_data = sorted(expr_adata.obs[sample_col].unique())
print(f"Found {len(samples_in_data)} sample(s): {samples_in_data}")

# 输出每个样本的细胞总数
print("\n=== Cell counts per sample ===")
sample_counts = expr_adata.obs[sample_col].value_counts().sort_index()
for sample, cnt in sample_counts.items():
    print(f"{sample}: {cnt} cells")
sample_counts.to_csv(os.path.join(output_dir, "sample_cell_counts.csv"))
print(f"Cell count table saved to {os.path.join(output_dir, 'sample_cell_counts.csv')}\n")

# 颜色映射
cmap = LinearSegmentedColormap.from_list('white_to_red',
    ['#F5F5F5', '#FFCCCC', '#FF9999', '#FF6666', '#FF3333', '#FF0000', '#CC0000'], N=256)

def get_expression(adata_sub, gene):
    if gene not in adata_sub.var_names:
        return None
    if hasattr(adata_sub, 'layers') and 'counts' in adata_sub.layers:
        X = adata_sub[:, gene].layers['counts']
    else:
        X = adata_sub[:, gene].X
    return X.toarray().flatten() if issparse(X) else np.asarray(X).flatten()

# 核心循环：基因 → 样本
for gene in valid_genes:
    for sample in samples_in_data:
        mask = expr_adata.obs[sample_col] == sample
        sub = expr_adata[mask].copy()
        if sub.n_obs == 0:
            continue

        expr = get_expression(sub, gene)
        if expr is None:
            continue

        coords = sub.obsm[spatial_key][:, :2]
        cell_count = sub.n_obs

        pts = 2 if cell_count > 10000 else 3 if cell_count > 5000 else 4 if cell_count > 2000 else 5

        vmin, vmax = np.percentile(expr, 1), np.percentile(expr, 99)
        if vmax - vmin < 0.001:
            vmax = vmin + 1

        plot_dir = os.path.join(output_dir, "individual_plots", gene, sample)
        os.makedirs(plot_dir, exist_ok=True)

        fig, ax = plt.subplots(figsize=(8, 8))
        scatter = ax.scatter(coords[:, 0], coords[:, 1], c=expr, cmap=cmap,
                             s=pts, edgecolors='none', alpha=0.7, vmin=vmin, vmax=vmax)
        ax.set_title(f"{gene} - {sample}", fontsize=16, fontweight='bold')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("X"); ax.set_ylabel("Y")

        x_range = coords[:, 0].max() - coords[:, 0].min()
        y_range = coords[:, 1].max() - coords[:, 1].min()
        pad = 0.05
        if x_range > 0 and y_range > 0:
            ax.set_xlim(coords[:, 0].min() - pad*x_range, coords[:, 0].max() + pad*x_range)
            ax.set_ylim(coords[:, 1].min() - pad*y_range, coords[:, 1].max() + pad*y_range)

        mean_expr = np.mean(expr)
        pos_rate = np.sum(expr > 0) / len(expr) * 100
        ax.text(0.02, 0.98,
                f"Cells: {cell_count:,}\nMean: {mean_expr:.3f}\nPositive: {pos_rate:.1f}%",
                transform=ax.transAxes, fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

        plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04, label='Expression')
        plt.tight_layout()

        file_base = f"{gene}_{sample}_spatial"
        plt.savefig(os.path.join(plot_dir, f"{file_base}.png"), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(plot_dir, f"{file_base}.pdf"), bbox_inches='tight')
        plt.close()
        print(f"Saved: {file_base}.png / .pdf")

print(f"\nDone! All plots saved to: {output_dir}")
