import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy.sparse import issparse
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# ================= 配置 =================
input_file = "All_filter_state_samp_new_with_CNs10-15.h5ad"
output_dir = "gene_classification_spatial_85"
genes_to_use = ["MYC", "AXIN2", "RNF43"]

# 需要保留的 5 种细胞类型
target_celltypes = [
    "WNT+ ISC",
    "IGF2+ prolif stem tumor Epi",
    "CXCL5+ Inf stem tumor Epi",
    "Hypo Stem tumor Epi",
    "ESC-like Epi"
]
# ========================================

os.makedirs(output_dir, exist_ok=True)

print(f"Reading {input_file}...")
adata = sc.read_h5ad(input_file)
print(f"Data shape: {adata.shape}")

# 检查 celltype 列是否存在
if 'celltype' not in adata.obs.columns:
    raise KeyError("'celltype' column not found in adata.obs. Please check your data.")

# ---- 标记目标细胞和其他细胞 ----
adata.obs['is_target'] = adata.obs['celltype'].isin(target_celltypes)
print(f"\nTarget cell types: {target_celltypes}")
print(f"Number of target cells: {adata.obs['is_target'].sum()}")
print(f"Number of other cells: {(~adata.obs['is_target']).sum()}")

# 统计其他细胞类型的组成（供参考）
other_types = adata.obs.loc[~adata.obs['is_target'], 'celltype'].value_counts()
print("\nOther cell types encountered:")
for ct, cnt in other_types.items():
    print(f"  {ct}: {cnt}")

# ---- 准备表达矩阵（log1p归一化） ----
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

# 确保 is_target 列存在
expr_adata.obs['is_target'] = adata.obs['is_target']

# ---- 检查基因存在性 ----
valid_genes = {}
for g in genes_to_use:
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

missing = [g for g in genes_to_use if valid_genes.get(g) is None]
if missing:
    raise SystemExit(f"Missing required genes: {missing}")

# ---- 仅使用目标细胞计算 75% 分位数阈值 ----
target_mask = expr_adata.obs['is_target'].values
thresholds = {}
for g in genes_to_use:
    gene = valid_genes[g]
    # 提取目标细胞的表达值
    X_all = expr_adata[:, gene].X
    if issparse(X_all):
        expr_all = X_all.toarray().flatten()
    else:
        expr_all = np.asarray(X_all).flatten()
    expr_target = expr_all[target_mask]
    if len(expr_target) == 0:
        raise RuntimeError(f"No target cells found for gene {g}, cannot compute threshold.")
    thr = np.percentile(expr_target, 85)
    thresholds[g] = thr
    print(f"Threshold for {g} (75% quantile, based on target cells only): {thr:.4f}")

# ---- 对每个细胞进行分类（仅目标细胞参与红/绿/黄/灰逻辑） ----
gene_myc = valid_genes["MYC"]
gene_axin2 = valid_genes["AXIN2"]
gene_rnf43 = valid_genes["RNF43"]

X_myc = expr_adata[:, gene_myc].X
X_axin2 = expr_adata[:, gene_axin2].X
X_rnf43 = expr_adata[:, gene_rnf43].X
if issparse(X_myc):
    myc_vals = X_myc.toarray().flatten()
    axin2_vals = X_axin2.toarray().flatten()
    rnf43_vals = X_rnf43.toarray().flatten()
else:
    myc_vals = np.asarray(X_myc).flatten()
    axin2_vals = np.asarray(X_axin2).flatten()
    rnf43_vals = np.asarray(X_rnf43).flatten()

thr_myc = thresholds["MYC"]
thr_axin2 = thresholds["AXIN2"]
thr_rnf43 = thresholds["RNF43"]

myc_high = myc_vals > thr_myc
axin2_high = axin2_vals > thr_axin2
rnf43_high = rnf43_vals > thr_rnf43

# 初始化分类标签（全部设为其他）
labels = np.full(len(myc_vals), "Other (lightblue)", dtype=object)

# 仅对目标细胞进行红/绿/黄/灰分类
target_idx = np.where(target_mask)[0]
myc_high_target = myc_high[target_idx]
axin2_high_target = axin2_high[target_idx]
rnf43_high_target = rnf43_high[target_idx]

# 黄色：MYC高且(AXIN2高或RNF43高)
yellow_mask_target = myc_high_target & (axin2_high_target | rnf43_high_target)
# 红色：仅MYC高
red_mask_target = myc_high_target & ~(axin2_high_target | rnf43_high_target)
# 绿色：仅AXIN2或RNF43高
green_mask_target = (axin2_high_target | rnf43_high_target) & ~myc_high_target
# 灰色：其余（均低）
gray_mask_target = ~(yellow_mask_target | red_mask_target | green_mask_target)

labels[target_idx[yellow_mask_target]] = "Yellow (MYC+ & AXIN2/RNF43+)"
labels[target_idx[red_mask_target]] = "Red (MYC high only)"
labels[target_idx[green_mask_target]] = "Green (AXIN2/RNF43 high only)"
labels[target_idx[gray_mask_target]] = "Gray (all low)"

expr_adata.obs["classification"] = labels

# ---- 整体统计（包括其他细胞） ----
print("\n=== Classification summary (all cells) ===")
class_counts = pd.Series(labels).value_counts()
for cls, cnt in class_counts.items():
    print(f"{cls}: {cnt} cells ({cnt/len(labels)*100:.1f}%)")

# 单独统计目标细胞中的红/绿/黄/灰比例
target_labels = labels[target_mask]
print("\n=== Classification within target cells only ===")
target_class_counts = pd.Series(target_labels).value_counts()
for cls, cnt in target_class_counts.items():
    print(f"{cls}: {cnt} cells ({cnt/len(target_labels)*100:.1f}%)")

# ---- 自动识别样本列和空间坐标 ----
sample_col = None
for col in ['sample', 'Sample', 'patient', 'Patient', 'batch', 'Batch', 'orig.ident', 'slice', 'Slide']:
    if col in expr_adata.obs.columns:
        sample_col = col
        break
if sample_col is None:
    print("No sample column found, treating all cells as one sample.")
    expr_adata.obs['sample'] = 'All_Samples'
    sample_col = 'sample'
print(f"\nSample column detected: '{sample_col}'")

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

samples_in_data = sorted(expr_adata.obs[sample_col].unique())
print(f"Found {len(samples_in_data)} sample(s): {samples_in_data}")

# ---- 颜色映射（增加浅蓝色） ----
color_map = {
    "Red (MYC high only)": "#E41A1C",
    "Green (AXIN2/RNF43 high only)": "#4DAF4A",
    "Yellow (MYC+ & AXIN2/RNF43+)": "#FFD700",
    "Gray (all low)": "#999999",
    "Other (lightblue)": "#ADD8E6"   # 浅蓝色
}

# ---- 为每个样本绘图（点大小固定为 1） ----
for sample in samples_in_data:
    mask = expr_adata.obs[sample_col] == sample
    sub = expr_adata[mask].copy()
    if sub.n_obs == 0:
        continue

    sub_labels = sub.obs["classification"]
    coords = sub.obsm[spatial_key][:, :2]
    cell_count = sub.n_obs

    # 点大小固定为 1
    pts = 1

    fig, ax = plt.subplots(figsize=(10, 10))
    # 按类别分别绘制以获得正确图例
    for cls, color in color_map.items():
        cls_mask = sub_labels == cls
        if np.sum(cls_mask) == 0:
            continue
        ax.scatter(coords[cls_mask, 0], coords[cls_mask, 1],
                   c=color, s=pts, edgecolors='none', alpha=0.7, label=cls)

    ax.set_title(f"Cell classification - {sample}", fontsize=16, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("X"); ax.set_ylabel("Y")

    # 坐标轴范围留白
    x_range = coords[:, 0].max() - coords[:, 0].min()
    y_range = coords[:, 1].max() - coords[:, 1].min()
    pad = 0.05
    if x_range > 0 and y_range > 0:
        ax.set_xlim(coords[:, 0].min() - pad*x_range, coords[:, 0].max() + pad*x_range)
        ax.set_ylim(coords[:, 1].min() - pad*y_range, coords[:, 1].max() + pad*y_range)

    # 统计信息文本框
    cnt = Counter(sub_labels)
    stats_text = f"Sample: {sample}\nTotal cells: {cell_count:,}\n"
    # 按固定顺序显示
    for cls in ["Yellow (MYC+ & AXIN2/RNF43+)", "Red (MYC high only)", 
                "Green (AXIN2/RNF43 high only)", "Gray (all low)", "Other (lightblue)"]:
        stats_text += f"{cls}: {cnt.get(cls, 0)} cells\n"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    ax.legend(loc='lower right', markerscale=2, fontsize=9)
    plt.tight_layout()

    plot_dir = os.path.join(output_dir, "individual_plots", sample)
    os.makedirs(plot_dir, exist_ok=True)
    file_base = f"classification_{sample}"
    plt.savefig(os.path.join(plot_dir, f"{file_base}.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(plot_dir, f"{file_base}.pdf"), bbox_inches='tight')
    plt.close()
    print(f"Saved: {file_base}.png / .pdf")

# 保存阈值信息（仅三个基因）
threshold_df = pd.DataFrame(list(thresholds.items()), columns=["Gene", "75%_threshold"])
threshold_df.to_csv(os.path.join(output_dir, "gene_thresholds.csv"), index=False)
print(f"\nGene thresholds saved to {os.path.join(output_dir, 'gene_thresholds.csv')}")

print(f"\nAll done! Results saved to: {output_dir}")
