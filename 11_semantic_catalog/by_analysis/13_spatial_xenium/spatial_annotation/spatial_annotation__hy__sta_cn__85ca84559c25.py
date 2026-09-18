#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
针对不同 CNs_10 子集重复进行细胞类型比例分析
子集1: CN2, CN7, CN8 合并
子集2: CN7, CN8 合并
每个子集内比较 MPC vs TNBC 组间差异（原始 P 值）
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import os, warnings
warnings.filterwarnings("ignore")

# ==================== 全局参数 ====================
H5AD_PATH = "filtered_with_cancer_type.h5ad"
CELL_TYPES = [
    "CD8+_Tem", "M2", "MKI67+_T", "B_follicular",
    "CD4+_Tex", "LYVE1+_TAM", "CD8+_Tcm", "NK",
    "IFN-T", "B_germinal_center", "EMT_M2", "Treg",
    "CD8+_Tex", "CD4+_Tem", "B_memory", "TAM"
]
GROUP_COL = "cancer_type"          # MPC vs TNBC
SAMPLE_COL = "sample"
ANNO_COL = "annotation"
CN_COL = "CNs_10"            # 拷贝数亚群列

# CNS 风格设置
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.8,
})
sns.set_style("white")
GROUP_COLORS = {"MPC": "#0072B2", "TNBC": "#D55E00"}

# ==================== 分析函数（复用） ====================
def run_analysis(adata_subset, output_dir):
    """在给定细胞子集上运行比例分析并绘图"""
    os.makedirs(output_dir, exist_ok=True)
    obs = adata_subset.obs.copy()
    # 确保类型为字符串
    obs[SAMPLE_COL] = obs[SAMPLE_COL].astype(str)
    obs[GROUP_COL] = obs[GROUP_COL].astype(str)
    obs[ANNO_COL] = obs[ANNO_COL].astype(str)

    # 1. 计算比例
    total_per_sample = obs.groupby(SAMPLE_COL).size().reset_index(name="total_cells")
    obs_target = obs[obs[ANNO_COL].isin(CELL_TYPES)].copy()
    if obs_target.empty:
        print("  警告：该子集中没有目标细胞类型，跳过分析。")
        return
    ct_counts = obs_target.groupby([SAMPLE_COL, ANNO_COL]).size().reset_index(name="cell_count")
    prop_df = ct_counts.merge(total_per_sample, on=SAMPLE_COL)
    prop_df["proportion"] = prop_df["cell_count"] / prop_df["total_cells"]

    # 添加分组
    sample_group = obs[[SAMPLE_COL, GROUP_COL]].drop_duplicates().set_index(SAMPLE_COL)
    prop_df[GROUP_COL] = prop_df[SAMPLE_COL].map(sample_group[GROUP_COL])

    # 总比例
    total_prop = prop_df.groupby([SAMPLE_COL, GROUP_COL])["proportion"].sum().reset_index()
    total_prop = total_prop.rename(columns={"proportion": "total_proportion"})

    # 只保留 MPC 和 TNBC
    prop_df = prop_df[prop_df[GROUP_COL].isin(["MPC", "TNBC"])]
    total_prop = total_prop[total_prop[GROUP_COL].isin(["MPC", "TNBC"])]
    if prop_df.empty:
        print("  警告：过滤 MPC/TNBC 后无数据，跳过分析。")
        return

    # 2. 总体比例图
    print(f"  正在绘制总体比例图...")
    fig, ax = plt.subplots(figsize=(3.5, 4.5))
    sns.boxplot(
        data=total_prop, x=GROUP_COL, y="total_proportion",
        order=["MPC", "TNBC"], palette=GROUP_COLORS, width=0.5,
        linewidth=1, ax=ax
    )
    sns.stripplot(
        data=total_prop, x=GROUP_COL, y="total_proportion",
        order=["MPC", "TNBC"], color="black", size=5, alpha=0.4, ax=ax
    )
    mpc_vals = total_prop[total_prop[GROUP_COL]=="MPC"]["total_proportion"].values
    tnbc_vals = total_prop[total_prop[GROUP_COL]=="TNBC"]["total_proportion"].values
    if len(mpc_vals)>=3 and len(tnbc_vals)>=3:
        _, p = mannwhitneyu(mpc_vals, tnbc_vals, alternative='two-sided')
        sig = ""
        if p < 0.001: sig = "***"
        elif p < 0.01: sig = "**"
        elif p < 0.05: sig = "*"
        if sig:
            y_max = total_prop["total_proportion"].max()
            y_min, y_max_ax = ax.get_ylim()
            y_line = y_max * 1.05
            ax.plot([0, 0, 1, 1], [y_line, y_line + 0.01*y_max, y_line + 0.01*y_max, y_line],
                    lw=1.2, c='black')
            ax.text(0.5, y_line + 0.015*y_max, sig, ha='center', va='bottom', fontsize=14, fontweight='bold')
            ax.set_ylim(y_min, y_line + 0.15*y_max)
        print(f"  总体比例 P = {p:.4f} {sig}")

    ax.set_ylabel("Proportion of selected cell types")
    ax.set_xlabel("")
    ax.set_title("Total selected cell type proportion")
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "overall_proportion_boxplot.pdf"))
    plt.close()

    # 3. 多细胞类型联合箱线图
    print(f"  正在绘制联合箱线图...")
    unique_ct = sorted(prop_df[ANNO_COL].unique())
    if not unique_ct:
        print("  警告：无细胞类型可绘制，跳过。")
        return

    p_dict = {}
    for ct in unique_ct:
        sub = prop_df[prop_df[ANNO_COL] == ct]
        mpc = sub[sub[GROUP_COL]=="MPC"]["proportion"].values
        tnbc = sub[sub[GROUP_COL]=="TNBC"]["proportion"].values
        if len(mpc)>=3 and len(tnbc)>=3:
            _, p = mannwhitneyu(mpc, tnbc, alternative='two-sided')
            p_dict[ct] = p
        else:
            p_dict[ct] = 1.0

    fig, ax = plt.subplots(figsize=(len(unique_ct)*0.8 + 2, 5))
    sns.boxplot(
        data=prop_df, x=ANNO_COL, y="proportion", hue=GROUP_COL,
        order=unique_ct, hue_order=["MPC", "TNBC"], palette=GROUP_COLORS,
        width=0.6, linewidth=0.8, fliersize=0, ax=ax
    )
    sns.stripplot(
        data=prop_df, x=ANNO_COL, y="proportion", hue=GROUP_COL,
        order=unique_ct, hue_order=["MPC", "TNBC"],
        color="black", dodge=True, size=3, alpha=0.5, ax=ax
    )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], ["MPC", "TNBC"], frameon=False, title="")

    y_max_all = prop_df["proportion"].max()
    if pd.isna(y_max_all):
        y_max_all = 0.0   # 安全回退
    current_y = y_max_all * 1.05
    y_step = y_max_all * 0.06

    for i, ct in enumerate(unique_ct):
        p = p_dict[ct]
        sig = ""
        if p < 0.001: sig = "***"
        elif p < 0.01: sig = "**"
        elif p < 0.05: sig = "*"
        if sig:
            x1, x2 = i - 0.2, i + 0.2
            ax.plot([x1, x1, x2, x2], [current_y, current_y + y_step/3, current_y + y_step/3, current_y],
                    lw=1.2, c='black', clip_on=False)
            ax.text((x1 + x2)/2, current_y + y_step/2, sig, ha='center', va='bottom',
                    fontsize=12, fontweight='bold')
            current_y += y_step
        print(f"  {ct}: P = {p:.4f} {sig}")

    # 防止 NaN 导致 set_ylim 失败
    top_limit = current_y + y_step
    if pd.isna(top_limit) or np.isinf(top_limit):
        top_limit = None
    if top_limit is not None:
        ax.set_ylim(top=top_limit)
    ax.set_ylabel("Proportion")
    ax.set_xlabel("")
    ax.set_title("Cell type proportions in MPC vs TNBC")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "celltype_proportion_multipanel_boxplot.pdf"))
    plt.close()
    print(f"  结果已保存至: {output_dir}\n")


# ==================== 主流程 ====================
print("读取数据...")
adata = sc.read_h5ad(H5AD_PATH)
if CN_COL not in adata.obs.columns:
    raise KeyError(f"obs 中缺少列: {CN_COL}")

# 定义两个子集
subsets = {
    "CN2CN7CN8": ["CN2", "CN7", "CN8"],
    "CN7CN8": ["CN7", "CN8"]
}

for folder_name, cn_values in subsets.items():
    print(f"\n===== 分析子集: {folder_name} (CNs = {cn_values}) =====")
    mask = adata.obs[CN_COL].isin(cn_values)
    adata_sub = adata[mask].copy()
    print(f"  细胞数: {adata_sub.n_obs}")
    if adata_sub.n_obs == 0:
        print("  警告: 该子集无细胞，跳过。")
        continue
    run_analysis(adata_sub, output_dir=folder_name)

print("\n全部分析完成！")
