#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
统计指定细胞类型的样本比例，比较 MPC vs TNBC
输出：1) 总体比例箱线图  2) 各细胞类型联合箱线图（原始 P 值显著性标注）
绘图风格：CNS 主刊级别
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
import os, warnings
warnings.filterwarnings("ignore")

# ==================== 参数设置 ====================
H5AD_PATH = "MPCvsTNBC_all_sample_ct.h5ad"
CELL_TYPES = [
    "CD8+_Tem", "M2", "MKI67+_T", "B_follicular",
    "CD4+_Tex", "LYVE1+_TAM", "CD8+_Tcm", "NK",
    "IFN-T", "B_germinal_center", "EMT_M2", "Treg",
    "CD8+_Tex", "CD4+_Tem", "B_memory", "TAM"
]
GROUP_COL = "group"
SAMPLE_COL = "sample"
ANNO_COL = "annotation"
OUTPUT_DIR = "./celltype_proportion_plots"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# CNS 级图形设置
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
GROUP_COLORS = {"MPC": "#0072B2", "TNBC": "#D55E00"}  # 色盲友好配色

# ==================== 1. 数据预处理 ====================
print("1. 读取数据并计算比例...")
adata = sc.read_h5ad(H5AD_PATH)
obs = adata.obs.copy()
for col in [SAMPLE_COL, GROUP_COL, ANNO_COL]:
    if col not in obs.columns:
        raise KeyError(f"obs 中缺少列: {col}")

obs[SAMPLE_COL] = obs[SAMPLE_COL].astype(str)
obs[GROUP_COL] = obs[GROUP_COL].astype(str)
obs[ANNO_COL] = obs[ANNO_COL].astype(str)

# 每个样本的总细胞数
total_per_sample = obs.groupby(SAMPLE_COL).size().reset_index(name="total_cells")
# 筛选目标细胞类型
obs_target = obs[obs[ANNO_COL].isin(CELL_TYPES)].copy()
# 每个样本每个细胞类型的数目
ct_counts = obs_target.groupby([SAMPLE_COL, ANNO_COL]).size().reset_index(name="cell_count")
# 合并总细胞数并计算比例
prop_df = ct_counts.merge(total_per_sample, on=SAMPLE_COL)
prop_df["proportion"] = prop_df["cell_count"] / prop_df["total_cells"]

# 添加分组信息
sample_group = obs[[SAMPLE_COL, GROUP_COL]].drop_duplicates().set_index(SAMPLE_COL)
prop_df[GROUP_COL] = prop_df[SAMPLE_COL].map(sample_group[GROUP_COL])

# 计算每个样本的总目标细胞比例（所有 CELL_TYPES 之和）
total_prop = prop_df.groupby([SAMPLE_COL, GROUP_COL])["proportion"].sum().reset_index()
total_prop = total_prop.rename(columns={"proportion": "total_proportion"})

# 只保留 MPC 和 TNBC
prop_df = prop_df[prop_df[GROUP_COL].isin(["MPC", "TNBC"])]
total_prop = total_prop[total_prop[GROUP_COL].isin(["MPC", "TNBC"])]

# ==================== 2. 总体比例箱线图 ====================
print("2. 绘制总体比例图...")
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
# 添加显著性
mpc_vals = total_prop[total_prop[GROUP_COL]=="MPC"]["total_proportion"].values
tnbc_vals = total_prop[total_prop[GROUP_COL]=="TNBC"]["total_proportion"].values
if len(mpc_vals)>=3 and len(tnbc_vals)>=3:
    stat, p = mannwhitneyu(mpc_vals, tnbc_vals, alternative='two-sided')
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
    # 打印 p 值
    print(f"总体比例: P = {p:.4f} {'*' if p<0.05 else ''}")

ax.set_ylabel("Proportion of selected cell types")
ax.set_xlabel("")
ax.set_title("Total selected cell type proportion")
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "overall_proportion_boxplot.pdf"))
plt.close()

# ==================== 3. 多细胞类型联合箱线图 ====================
print("3. 绘制各细胞类型联合箱线图...")
unique_ct = sorted(prop_df[ANNO_COL].unique())
n_ct = len(unique_ct)

# 计算每个细胞类型的 P 值
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

# 图形尺寸
fig, ax = plt.subplots(figsize=(n_ct*0.8 + 2, 5))

# 绘制箱线图：每个细胞类型两个箱子
sns.boxplot(
    data=prop_df, x=ANNO_COL, y="proportion", hue=GROUP_COL,
    order=unique_ct, hue_order=["MPC", "TNBC"], palette=GROUP_COLORS,
    width=0.6, linewidth=0.8, fliersize=0, ax=ax
)
# 添加散点（抖动）
sns.stripplot(
    data=prop_df, x=ANNO_COL, y="proportion", hue=GROUP_COL,
    order=unique_ct, hue_order=["MPC", "TNBC"],
    color="black", dodge=True, size=3, alpha=0.5, ax=ax
)

# 图例只保留两个颜色框
handles, labels = ax.get_legend_handles_labels()
# 剔除散点图例（重复的），只保留箱线图图例
n = len(unique_ct)*2  # 散点图例在后面
ax.legend(handles[:2], ["MPC", "TNBC"], frameon=False, title="")

# 标注显著性
y_max_all = prop_df["proportion"].max()
y_offsets = {}  # 记录每个细胞类型的标注线高度
current_y = y_max_all * 1.05
y_step = y_max_all * 0.06

for i, ct in enumerate(unique_ct):
    p = p_dict[ct]
    sig = ""
    if p < 0.001: sig = "***"
    elif p < 0.01: sig = "**"
    elif p < 0.05: sig = "*"
    if sig:
        # 在组间画线
        x1, x2 = i - 0.2, i + 0.2
        ax.plot([x1, x1, x2, x2], [current_y, current_y + y_step/3, current_y + y_step/3, current_y],
                lw=1.2, c='black', clip_on=False)
        ax.text((x1 + x2)/2, current_y + y_step/2, sig, ha='center', va='bottom', fontsize=12, fontweight='bold')
        current_y += y_step
    # 打印 p 值信息
    print(f"{ct}: P = {p:.4f} {sig}")

# 调整 y 轴上限以容纳标注
ax.set_ylim(top=current_y + y_step)
ax.set_ylabel("Proportion")
ax.set_xlabel("")
ax.set_title("Cell type proportions in MPC vs TNBC")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "celltype_proportion_multipanel_boxplot.pdf"))
plt.close()

print(f"所有图形已保存至: {OUTPUT_DIR}")
