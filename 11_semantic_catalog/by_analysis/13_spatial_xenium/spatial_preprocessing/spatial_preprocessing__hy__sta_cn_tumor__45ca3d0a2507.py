import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, ttest_ind, ks_2samp
from statsmodels.stats.multitest import multipletests
import os

# ========================= 配置 =========================
h5ad_path = "All_filter_with_region_CNs10.h5ad"
output_dir = "./cn_comparison_plots"
os.makedirs(output_dir, exist_ok=True)

group1 = "PRE-R"
group2 = "PRE-NR"

# ========================= 加载并转换数据 =========================
print("加载数据...")
adata = sc.read_h5ad(h5ad_path)
obs = adata.obs.copy()

if "CNs_10" not in obs.columns:
    raise KeyError("数据中无 CNs_10 列！")

# 转换为哑变量，列名 CN0 ~ CN9
cn_dummies = pd.get_dummies(obs["CNs_10"])
cn_cols = [f"CN{i}" for i in range(10)]
for c in cn_cols:
    if c not in cn_dummies.columns:
        cn_dummies[c] = 0
cn_dummies = cn_dummies[cn_cols]
obs = pd.concat([obs, cn_dummies], axis=1)

# ========================= 筛选肿瘤区域 + 目标分组 =========================
tumor_mask = obs["region"] == "tumor"
group_mask = obs["group_final"].isin([group1, group2])
obs_filt = obs[tumor_mask & group_mask].copy()
print(f"Tumor + 目标组细胞数: {obs_filt.shape[0]}")

# ========================= 计算每个样本各 CN 占比 =========================
grouped = obs_filt.groupby("sample")
total_cells = grouped.size()
cn_counts = grouped[cn_cols].sum()
cn_ratio = cn_counts.div(total_cells, axis=0)

# 添加分组信息
sample_group = obs_filt[["sample", "group_final"]].drop_duplicates().set_index("sample")
cn_ratio = cn_ratio.join(sample_group)

print(f"参与比较的样本数: {cn_ratio.shape[0]} ({group1}: {(cn_ratio['group_final']==group1).sum()}, {group2}: {(cn_ratio['group_final']==group2).sum()})")

# ========================= 输出占比矩阵 =========================
ratio_matrix_path = os.path.join(output_dir, "cn_ratio_matrix.csv")
cn_ratio.to_csv(ratio_matrix_path)
print(f"占比矩阵已保存至: {ratio_matrix_path}")

# ========================= 执行多种统计检验 =========================
results = []
pvals_mw = []  # 用于FDR校正的Mann-Whitney p值

for cn in cn_cols:
    vals1 = cn_ratio[cn_ratio["group_final"] == group1][cn].dropna().values
    vals2 = cn_ratio[cn_ratio["group_final"] == group2][cn].dropna().values

    if len(vals1) == 0 or len(vals2) == 0:
        print(f"跳过 {cn}: 某组无样本")
        continue

    res_entry = {"CN": cn,
                 f"{group1}_n": len(vals1),
                 f"{group2}_n": len(vals2),
                 f"{group1}_mean": np.mean(vals1),
                 f"{group2}_mean": np.mean(vals2),
                 f"{group1}_median": np.median(vals1),
                 f"{group2}_median": np.median(vals2)}

    # 1. Mann-Whitney U (双侧)
    u_stat, mw_p = mannwhitneyu(vals1, vals2, alternative="two-sided")
    res_entry["MW_U"] = u_stat
    res_entry["MW_p"] = mw_p
    pvals_mw.append(mw_p)

    # 2. Welch's t-test (双侧)
    t_stat, t_p = ttest_ind(vals1, vals2, equal_var=False, alternative="two-sided")
    res_entry["Welch_t"] = t_stat
    res_entry["Welch_p"] = t_p

    # 3. Kolmogorov-Smirnov (双侧)
    ks_stat, ks_p = ks_2samp(vals1, vals2, alternative="two-sided")
    res_entry["KS_stat"] = ks_stat
    res_entry["KS_p"] = ks_p

    results.append(res_entry)

res_df = pd.DataFrame(results)

# ========================= 多重检验校正（对MW和Welch和KS的p值） =========================
# 这里以MW的p值为例做FDR校正，也可对其他检验做
if pvals_mw:
    reject, qvals, _, _ = multipletests(pvals_mw, method="fdr_bh")
    res_df["MW_q"] = qvals
    res_df["MW_FDR_sig"] = ["***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns" for q in qvals]
else:
    res_df["MW_q"] = np.nan
    res_df["MW_FDR_sig"] = "ns"

# 保存完整统计结果
result_csv_path = os.path.join(output_dir, "cn_comparison_multi_test.csv")
res_df.to_csv(result_csv_path, index=False, float_format="%.6e")
print(f"多检验结果表格已保存至: {result_csv_path}")

# ========================= 绘图（仍使用MW原始p值标注星号） =========================
sns.set_style("whitegrid")
plt.rcParams["figure.dpi"] = 150

for i, row in res_df.iterrows():
    cn = row["CN"]
    # 提取绘图数据
    plot_df = cn_ratio[[cn, "group_final"]].copy()
    plot_df.columns = ["ratio", "group"]

    # 使用原始MW p值来标记星号（你也可以换成MW_q）
    p_val = row["MW_p"]
    if p_val < 0.001:
        sig_label = "***"
    elif p_val < 0.01:
        sig_label = "**"
    elif p_val < 0.05:
        sig_label = "*"
    else:
        sig_label = "ns"

    fig, ax = plt.subplots(figsize=(5, 5))
    palette = {group1: "#E64B35", group2: "#4DBBD5"}
    order = [group1, group2]

    sns.boxplot(x="group", y="ratio", data=plot_df, order=order,
                palette=palette, width=0.5, linewidth=1.2, ax=ax)
    sns.stripplot(x="group", y="ratio", data=plot_df, order=order,
                  color="black", size=6, alpha=0.7, jitter=True, ax=ax)

    if sig_label != "ns":
        y_max = plot_df["ratio"].max()
        y_range = plot_df["ratio"].max() - plot_df["ratio"].min()
        if y_range == 0:
            y_range = 0.01
        y_line = y_max + 0.05 * y_range
        y_star = y_line + 0.02 * y_range
        ax.plot([0, 0, 1, 1], [y_line, y_star, y_star, y_line],
                lw=1.5, c="black")
        ax.text(0.5, y_star, sig_label, ha="center", va="bottom",
                fontsize=14, fontweight="bold")

    ax.set_title(f"{cn} 占比 (Tumor)", fontsize=14)
    ax.set_ylabel(f"{cn} 阳性细胞占比")
    ax.set_xlabel("Group")
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"{cn}_comparison.png")
    fig.savefig(save_path)
    plt.close(fig)
    print(f"已保存图片: {save_path}")

print("所有分析完成！")
