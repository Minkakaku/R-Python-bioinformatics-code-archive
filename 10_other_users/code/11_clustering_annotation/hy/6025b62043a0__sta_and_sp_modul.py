import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy.sparse import issparse
from collections import Counter
from scipy.stats import chi2_contingency, kruskal
from itertools import combinations
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# ================= 配置 =================
input_file = "All_filter_state_samp_new_with_CNs10-15.h5ad"
output_dir = "gene_classification_spatial_75_new"
genes_to_use = ["MYC", "AXIN2", "RNF43"]

# 需要保留的 5 种细胞类型
target_celltypes = [
    "WNT+ ISC",
    "IGF2+ prolif stem tumor Epi",
    "CXCL5+ Inf stem tumor Epi",
    "Hypo Stem tumor Epi",
    "ESC-like Epi"
]

# 分组列：用于比较不同组间红/绿/黄细胞比例的差异
group_col = 'state'
# ========================================

os.makedirs(output_dir, exist_ok=True)

print(f"Reading {input_file}...")
adata = sc.read_h5ad(input_file)
print(f"Data shape: {adata.shape}")

# 检查必需列是否存在
if 'celltype' not in adata.obs.columns:
    raise KeyError("'celltype' column not found in adata.obs. Please check your data.")
if group_col not in adata.obs.columns:
    raise KeyError(f"'{group_col}' column not found in adata.obs. Please check your data.")

# ---- 标记目标细胞和其他细胞 ----
adata.obs['is_target'] = adata.obs['celltype'].isin(target_celltypes)
print(f"\nTarget cell types: {target_celltypes}")
print(f"Number of target cells: {adata.obs['is_target'].sum()}")
print(f"Number of other cells: {(~adata.obs['is_target']).sum()}")

# 统计其他细胞类型的组成
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
    X_all = expr_adata[:, gene].X
    if issparse(X_all):
        expr_all = X_all.toarray().flatten()
    else:
        expr_all = np.asarray(X_all).flatten()
    expr_target = expr_all[target_mask]
    if len(expr_target) == 0:
        raise RuntimeError(f"No target cells found for gene {g}, cannot compute threshold.")
    thr = np.percentile(expr_target, 75)
    thresholds[g] = thr
    print(f"Threshold for {g} (75% quantile, based on target cells only): {thr:.4f}")

# ---- 对每个细胞进行分类 ----
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

labels = np.full(len(myc_vals), "Other (lightblue)", dtype=object)

target_idx = np.where(target_mask)[0]
myc_high_target = myc_high[target_idx]
axin2_high_target = axin2_high[target_idx]
rnf43_high_target = rnf43_high[target_idx]

yellow_mask_target = myc_high_target & (axin2_high_target | rnf43_high_target)
red_mask_target = myc_high_target & ~(axin2_high_target | rnf43_high_target)
green_mask_target = (axin2_high_target | rnf43_high_target) & ~myc_high_target
gray_mask_target = ~(yellow_mask_target | red_mask_target | green_mask_target)

labels[target_idx[yellow_mask_target]] = "Yellow (MYC+ & AXIN2/RNF43+)"
labels[target_idx[red_mask_target]] = "Red (MYC high only)"
labels[target_idx[green_mask_target]] = "Green (AXIN2/RNF43 high only)"
labels[target_idx[gray_mask_target]] = "Gray (all low)"

expr_adata.obs["classification"] = labels

# ---- 整体统计 ----
print("\n=== Classification summary (all cells) ===")
class_counts = pd.Series(labels).value_counts()
for cls, cnt in class_counts.items():
    print(f"{cls}: {cnt} cells ({cnt/len(labels)*100:.1f}%)")

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

# ---- 颜色映射 ----
color_map = {
    "Red (MYC high only)": "#E41A1C",
    "Green (AXIN2/RNF43 high only)": "#4DAF4A",
    "Yellow (MYC+ & AXIN2/RNF43+)": "#FFD700",
    "Gray (all low)": "#999999",
    "Other (lightblue)": "#ADD8E6"
}

# ---- 为每个样本绘图 ----
for sample in samples_in_data:
    mask = expr_adata.obs[sample_col] == sample
    sub = expr_adata[mask].copy()
    if sub.n_obs == 0:
        continue

    sub_labels = sub.obs["classification"]
    coords = sub.obsm[spatial_key][:, :2]
    cell_count = sub.n_obs
    pts = 1

    fig, ax = plt.subplots(figsize=(10, 10))
    for cls, color in color_map.items():
        cls_mask = sub_labels == cls
        if np.sum(cls_mask) == 0:
            continue
        ax.scatter(coords[cls_mask, 0], coords[cls_mask, 1],
                   c=color, s=pts, edgecolors='none', alpha=0.7, label=cls)

    ax.set_title(f"Cell classification - {sample}", fontsize=16, fontweight='bold')
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("X"); ax.set_ylabel("Y")

    x_range = coords[:, 0].max() - coords[:, 0].min()
    y_range = coords[:, 1].max() - coords[:, 1].min()
    pad = 0.05
    if x_range > 0 and y_range > 0:
        ax.set_xlim(coords[:, 0].min() - pad*x_range, coords[:, 0].max() + pad*x_range)
        ax.set_ylim(coords[:, 1].min() - pad*y_range, coords[:, 1].max() + pad*y_range)

    cnt = Counter(sub_labels)
    stats_text = f"Sample: {sample}\nTotal cells: {cell_count:,}\n"
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

# 保存阈值信息
threshold_df = pd.DataFrame(list(thresholds.items()), columns=["Gene", "75%_threshold"])
threshold_df.to_csv(os.path.join(output_dir, "gene_thresholds.csv"), index=False)
print(f"\nGene thresholds saved to {os.path.join(output_dir, 'gene_thresholds.csv')}")

# ==================== 红/绿/黄细胞在 state 组间的差异比较（细胞水平） ====================
print(f"\n========== Group comparison for Red/Green/Yellow cells by '{group_col}' ==========")

groups = expr_adata.obs[group_col].astype(str).unique()
if len(groups) < 2:
    print("Only one group found. Cannot perform statistical comparison.")
else:
    # 筛选属于红/绿/黄的细胞
    rgy_mask = expr_adata.obs["classification"].isin([
        "Red (MYC high only)",
        "Green (AXIN2/RNF43 high only)",
        "Yellow (MYC+ & AXIN2/RNF43+)"
    ])
    rgy_df = expr_adata.obs.loc[rgy_mask, [group_col, "classification"]].copy()
    rgy_df.columns = ["group", "class"]
    print(f"\nTotal Red/Green/Yellow cells across groups: {len(rgy_df)}")

    # 构建列联表: 行=group, 列=class
    crosstab = pd.crosstab(rgy_df["group"], rgy_df["class"])
    for cls in ["Red (MYC high only)", "Green (AXIN2/RNF43 high only)", "Yellow (MYC+ & AXIN2/RNF43+)"]:
        if cls not in crosstab.columns:
            crosstab[cls] = 0
    crosstab = crosstab[["Red (MYC high only)", "Green (AXIN2/RNF43 high only)", "Yellow (MYC+ & AXIN2/RNF43+)"]]

    print("\n--- Contingency table (counts) ---")
    print(crosstab)

    prop_table = crosstab.div(crosstab.sum(axis=1), axis=0)
    print("\n--- Proportion table ---")
    print(prop_table.round(4))

    crosstab.to_csv(os.path.join(output_dir, "rgy_group_counts.csv"))
    prop_table.to_csv(os.path.join(output_dir, "rgy_group_proportions.csv"))
    print(f"\nContingency tables saved to {output_dir}")

    try:
        chi2, p_overall, dof, expected = chi2_contingency(crosstab)
        print(f"\nOverall chi-square test: chi2 = {chi2:.4f}, p = {p_overall:.4e}, dof = {dof}")
    except ValueError as e:
        print(f"Chi-square test failed: {e}")
        chi2, p_overall = np.nan, np.nan

    group_list = crosstab.index.tolist()
    pair_results = []
    for g1, g2 in combinations(group_list, 2):
        sub_table = crosstab.loc[[g1, g2]]
        try:
            _, p_pair, _, _ = chi2_contingency(sub_table)
        except ValueError:
            p_pair = np.nan
        pair_results.append((g1, g2, p_pair))

    if pair_results:
        print("\n--- Pairwise group comparisons (chi-square test) ---")
        for g1, g2, p_val in pair_results:
            print(f"{g1} vs {g2}: p = {p_val:.4e}")

    stats_summary = pd.DataFrame({
        "Test": ["Overall chi-square"],
        "Statistic": [chi2],
        "p_value": [p_overall]
    })
    if pair_results:
        pair_df = pd.DataFrame(pair_results, columns=["Group1", "Group2", "p_value"])
        stats_summary = pd.concat([stats_summary, pair_df], ignore_index=True)
    stats_summary.to_csv(os.path.join(output_dir, "rgy_group_test_results.csv"), index=False)
    print(f"Statistical results saved to {os.path.join(output_dir, 'rgy_group_test_results.csv')}")

    # 堆叠条形图
    fig, ax = plt.subplots(figsize=(max(6, len(group_list)*1.2), 5))
    colors = [color_map["Red (MYC high only)"],
              color_map["Green (AXIN2/RNF43 high only)"],
              color_map["Yellow (MYC+ & AXIN2/RNF43+)"]]
    prop_table.plot(kind="bar", stacked=True, ax=ax, color=colors, edgecolor="white", width=0.7)
    ax.set_title(f"Proportion of Red/Green/Yellow cells by {group_col}", fontsize=14, fontweight='bold')
    ax.set_xlabel(group_col)
    ax.set_ylabel("Proportion")
    ax.legend(title="Classification", bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.tight_layout()
    fig_path = os.path.join(output_dir, "rgy_group_proportions_barplot.png")
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.savefig(fig_path.replace(".png", ".pdf"), bbox_inches='tight')
    plt.close()
    print(f"Bar plot saved to {fig_path}")

# ==================== 新增：样本-状态组合的箱线图（样本维度统计） ====================
print(f"\n========== Sample-state combination boxplots ==========")

# 创建样本-状态组合列
expr_adata.obs['sample_state'] = (
    expr_adata.obs[sample_col].astype(str) + "_" + expr_adata.obs[group_col].astype(str)
)

# 仅保留目标细胞，因为红/绿/黄分类只针对目标细胞
target_adata = expr_adata[expr_adata.obs['is_target']].copy()

# 按 sample_state 统计分类计数
count_df = target_adata.obs.groupby('sample_state')['classification'].value_counts().unstack(fill_value=0)
# 确保包含红绿黄三类
for cls in ['Red (MYC high only)', 'Green (AXIN2/RNF43 high only)', 'Yellow (MYC+ & AXIN2/RNF43+)']:
    if cls not in count_df.columns:
        count_df[cls] = 0

# 总目标细胞数
total_target = count_df.sum(axis=1)
# 计算各分类在目标细胞中的比例
prop_df = count_df.div(total_target, axis=0)

# 获取每个 sample_state 对应的 state 分组（安全方法：从原数据提取映射）
state_map = expr_adata.obs[['sample_state', group_col]].drop_duplicates().set_index('sample_state')[group_col].astype(str)
prop_df['state'] = state_map.loc[prop_df.index]

# 可选：过滤掉目标细胞数过少的组合（例如<10），可根据需要调整
min_cells = 10
prop_df = prop_df.loc[total_target >= min_cells].copy()
if len(prop_df) < len(total_target):
    removed = len(total_target) - len(prop_df)
    print(f"Filtered out {removed} sample-state combinations with fewer than {min_cells} target cells.")

print(f"\nNumber of valid sample-state combinations: {len(prop_df)}")

# 绘制三个分类的箱线图（按 state 分组）
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
class_names = ['Red (MYC high only)', 'Green (AXIN2/RNF43 high only)', 'Yellow (MYC+ & AXIN2/RNF43+)']
for i, cls in enumerate(class_names):
    ax = axes[i]
    # 箱线图 + 散点
    sns.boxplot(x='state', y=cls, data=prop_df.reset_index(), ax=ax, palette='Set2', width=0.5)
    sns.stripplot(x='state', y=cls, data=prop_df.reset_index(), ax=ax, color='black', alpha=0.7, size=6)
    ax.set_title(f'{cls.split(" ")[0]} cells proportion', fontsize=12)
    ax.set_ylabel('Proportion in target cells')
    ax.set_xlabel(group_col)
    ax.tick_params(axis='x', rotation=45)
plt.tight_layout()
fig_path = os.path.join(output_dir, "sample_state_boxplots.png")
plt.savefig(fig_path, dpi=300)
plt.savefig(fig_path.replace('.png', '.pdf'))
plt.close()
print(f"Boxplots saved to {fig_path}")

# 保存比例表
prop_df.to_csv(os.path.join(output_dir, "sample_state_proportions.csv"))
print(f"Sample-state proportions saved to {os.path.join(output_dir, 'sample_state_proportions.csv')}")

# 对每个分类进行 Kruskal-Wallis 检验（样本级别）
print("\n--- Kruskal-Wallis test for each category (sample-state level) ---")
for cls in class_names:
    groups_data = [group[cls].values for name, group in prop_df.groupby('state')]
    if len(groups_data) >= 2:
        stat, p = kruskal(*groups_data)
        print(f"{cls}: statistic = {stat:.4f}, p = {p:.4e}")
    else:
        print(f"{cls}: not enough groups for test")

print(f"\nAll done! Results saved to: {output_dir}")
