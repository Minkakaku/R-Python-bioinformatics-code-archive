#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
按样本维度统计 Tumor 区域细胞类型组成差异：
  比较1：PRE-R vs PRE-NR
  比较2：PRE-NR vs POST-NR
每个细胞类型绘制单独箱线图，差异显著时标注 p 值和显著性星号。
"""

import os
import re
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import scanpy as sc

warnings.filterwarnings('ignore')

# ==================== 配置 ====================
H5AD_PATH = "All_filter_with_region_CNs10.h5ad"
OUTPUT_DIR = "tumor_composition_analysis"
FDR_THRESHOLD = 0.1             # 多重检验校正后的显著性阈值
MIN_CELLS_PER_SAMPLE = 20       # 样本最低 Tumor 细胞数
MIN_SAMPLES_PER_GROUP = 2       # 每组最少样本数

# 两个独立的比较任务
COMPARISONS = [
    (["PRE-R", "PRE-NR"], "PRE-R_vs_PRE-NR"),
    (["PRE-NR", "POST-NR"], "PRE-NR_vs_POST-NR"),
]

# ==================== 函数 ====================
def extract_comparison_data(adata, groups, region='tumor'):
    """
    从 adata 中提取指定 groups 在 tumor 区域的细胞，返回仅含两个组的 AnnData 子集。
    """
    mask = (adata.obs['region'] == region) & (adata.obs['group_final'].isin(groups))
    return adata[mask].copy()

def compute_proportions(adata_sub):
    """
    计算每个样本的细胞类型比例（仅基于该样本的 Tumor 细胞）。
    返回：DataFrame，行为样本，列为细胞类型，值为比例（0~1）。
    """
    # 统计样本 × 细胞类型 的细胞数
    count_df = adata_sub.obs.groupby(['sample', 'celltype']).size().unstack(fill_value=0)
    # 只保留 Tumor 细胞总数 >= MIN_CELLS_PER_SAMPLE 的样本
    sample_totals = count_df.sum(axis=1)
    valid_samples = sample_totals[sample_totals >= MIN_CELLS_PER_SAMPLE].index
    count_df = count_df.loc[valid_samples]
    # 转换为比例
    prop_df = count_df.div(count_df.sum(axis=1), axis=0)
    return prop_df

def add_group_info(prop_df, adata_sub):
    """将样本对应的分组信息添加到比例矩阵"""
    sample_group = adata_sub.obs[['sample', 'group_final']].drop_duplicates().set_index('sample')
    prop_df = prop_df.join(sample_group)
    return prop_df

def run_differential_test(prop_df_with_group):
    """
    对比例矩阵中的每个细胞类型进行 Mann-Whitney U 检验，
    返回包含 p 值、FDR、中位数、log2FC 的结果 DataFrame。
    """
    df = prop_df_with_group.copy()
    groups = df['group_final'].unique()
    if len(groups) != 2:
        raise ValueError(f"需要恰好两个组，但得到: {groups}")
    g1, g2 = groups

    mask1 = df['group_final'] == g1
    mask2 = df['group_final'] == g2
    data1 = df[mask1].drop(columns='group_final')
    data2 = df[mask2].drop(columns='group_final')

    if data1.shape[0] < MIN_SAMPLES_PER_GROUP or data2.shape[0] < MIN_SAMPLES_PER_GROUP:
        return None

    results = []
    for ct in data1.columns:
        v1 = data1[ct].values
        v2 = data2[ct].values
        if np.array_equal(v1, v2):
            pval = 1.0
        else:
            try:
                _, pval = mannwhitneyu(v1, v2, alternative='two-sided')
            except ValueError:
                pval = 1.0

        med1 = np.median(v1) if len(v1) > 0 else 0
        med2 = np.median(v2) if len(v2) > 0 else 0
        pseudo = 1e-6
        log2fc = np.log2((med1 + pseudo) / (med2 + pseudo))

        results.append({
            'celltype': ct,
            f'median_{g1}': med1,
            f'median_{g2}': med2,
            'log2FC': log2fc,
            'p_value': pval
        })

    res_df = pd.DataFrame(results)
    # 多重比较校正
    if len(res_df) > 0:
        _, fdr, _, _ = multipletests(res_df['p_value'], method='fdr_bh')
        res_df['FDR'] = fdr
        res_df['significant'] = res_df['FDR'] < FDR_THRESHOLD
    return res_df

def plot_one_celltype(prop_df, group_col, ct, res_row, comp_name, out_dir):
    """绘制单个细胞类型的箱线图，仅包含该任务的两个组"""
    plot_data = prop_df[[ct, group_col]].copy()
    groups_present = sorted(plot_data[group_col].unique())  # 保证顺序一致
    if len(groups_present) != 2:
        return

    fig, ax = plt.subplots(figsize=(5, 6))
    sns.boxplot(data=plot_data, x=group_col, y=ct, order=groups_present, palette='Set2', ax=ax)
    sns.stripplot(data=plot_data, x=group_col, y=ct, order=groups_present,
                  color='black', size=5, alpha=0.6, ax=ax)

    # 标注显著性（如果存在且显著）
    if res_row is not None and res_row['significant']:
        pval = res_row['p_value']
        star_str = ''
        if pval < 0.001:
            star_str = '***'
        elif pval < 0.01:
            star_str = '**'
        elif pval < 0.05:
            star_str = '*'

        # 确定横线的 y 位置
        all_vals = plot_data[ct].values
        y_max = all_vals.max()
        y_range = all_vals.max() - all_vals.min() if all_vals.max() > all_vals.min() else 0.1
        y_line = y_max + y_range * 0.1
        h = y_range * 0.03

        # 画出连接两组的横线
        x1, x2 = 0, 1
        ax.plot([x1, x1, x2, x2], [y_line, y_line + h, y_line + h, y_line],
                lw=1.5, c='black')
        # 在横线上方写 p 值和星号
        ax.text((x1 + x2) / 2, y_line + h + y_range * 0.02,
                f"p = {pval:.4f} {star_str}",
                ha='center', va='bottom', fontsize=11)

    ax.set_ylabel('Proportion')
    ax.set_title(f'{ct}\n{comp_name}')
    plt.tight_layout()
    # 生成合法文件名
    safe_ct = re.sub(r'[\\/*?:"<>| ]', '_', ct)
    fpath = os.path.join(out_dir, f"{safe_ct}.png")
    fig.savefig(fpath, dpi=300, bbox_inches='tight')
    plt.close(fig)

# ==================== 主流程 ====================
def main():
    print(f"读取数据: {H5AD_PATH}")
    adata = sc.read_h5ad(H5AD_PATH)

    # 检查必要列
    for col in ['region', 'group_final', 'sample', 'celltype']:
        if col not in adata.obs.columns:
            raise KeyError(f"obs 中缺少列: {col}")

    print(f"原始细胞数: {adata.n_obs}")
    tumor_mask = adata.obs['region'] == 'tumor'
    print(f"Tumor 区域细胞数: {tumor_mask.sum()}")

    # 创建输出根目录
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for groups_list, comp_name in COMPARISONS:
        print(f"\n{'='*60}")
        print(f"处理比较: {comp_name}")
        # 1. 提取仅含这两个组的 Tumor 数据
        adata_comp = extract_comparison_data(adata, groups_list)
        print(f"  用于比较的细胞总数: {adata_comp.n_obs}")
        if adata_comp.n_obs == 0:
            print("  无细胞，跳过")
            continue

        # 2. 计算样本-细胞类型比例
        prop_df = compute_proportions(adata_comp)
        if prop_df.empty:
            print("  无满足细胞数要求的样本，跳过")
            continue
        prop_df = add_group_info(prop_df, adata_comp)
        print(f"  样本数: {prop_df.shape[0]}, 细胞类型数: {prop_df.shape[1]-1}")
        print(f"  各组样本数: {prop_df['group_final'].value_counts().to_dict()}")

        # 3. 差异检验
        res_df = run_differential_test(prop_df)
        if res_df is None:
            print("  统计检验失败（可能样本数不足）")
            continue
        # 保存检验结果
        comp_out_dir = os.path.join(OUTPUT_DIR, comp_name)
        os.makedirs(comp_out_dir, exist_ok=True)
        res_df.to_csv(os.path.join(comp_out_dir, "differential_abundance.csv"), index=False)
        print(f"  显著细胞类型数 (FDR<{FDR_THRESHOLD}): {res_df['significant'].sum()}")

        # 4. 为每个细胞类型绘制箱线图
        celltypes = prop_df.columns.drop('group_final')
        for ct in celltypes:
            row = res_df[res_df['celltype'] == ct]
            if not row.empty:
                row_data = row.iloc[0]
            else:
                row_data = None
            plot_one_celltype(prop_df, 'group_final', ct, row_data, comp_name, comp_out_dir)

        print(f"  箱线图已保存至: {comp_out_dir}")

    print("\n所有分析完成！")

if __name__ == "__main__":
    main()
