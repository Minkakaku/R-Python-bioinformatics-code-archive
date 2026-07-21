#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
按样本维度统计 Tumor vs Normal 区域细胞类型组成差异：
  比较同一批样本在 tumor 与 normal 两个区域内的细胞类型比例，
  使用 Mann-Whitney U 检验（非配对），多重检验校正后标注显著性。
  输出结果至 ct_comp 目录。
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
OUTPUT_DIR = "ct_comp"                # 输出目录
FDR_THRESHOLD = 0.1                   # 多重检验校正后的显著性阈值
MIN_CELLS_PER_SAMPLE_REGION = 20      # 每个样本在单个区域的最低细胞数
MIN_SAMPLES_PER_REGION = 2            # 每个区域最少样本数

# ==================== 函数 ====================
def compute_region_proportions(adata):
    """
    分别统计每个样本在 tumor / normal 区域的细胞类型比例。
    返回: DataFrame，索引为 (sample, region)，列为细胞类型 + region 列
    """
    # 只保留 tumor 和 normal 区域
    valid_regions = ['tumor', 'normal']
    mask = adata.obs['region'].isin(valid_regions)
    sub = adata[mask].copy()

    # 统计 sample × region × celltype 的细胞数
    count_df = sub.obs.groupby(['sample', 'region', 'celltype']).size().unstack(fill_value=0)
    # 过滤每个 (sample, region) 组合的总细胞数 >= MIN_CELLS_PER_SAMPLE_REGION
    totals = count_df.sum(axis=1)
    valid_idx = totals[totals >= MIN_CELLS_PER_SAMPLE_REGION].index
    count_df = count_df.loc[valid_idx]

    if count_df.empty:
        return pd.DataFrame()

    # 转换为比例
    prop_df = count_df.div(count_df.sum(axis=1), axis=0)

    # 将多级索引转为列，得到 sample, region, 各细胞类型比例
    prop_df = prop_df.reset_index()  # 现在包含 sample, region 两列 + 细胞类型列
    return prop_df

def run_differential_test(prop_df):
    """
    对比例矩阵进行 tumor vs normal 的 Mann-Whitney U 检验，
    返回结果 DataFrame。
    """
    df = prop_df.copy()
    regions = df['region'].unique()
    # 确保只有两个区域
    if set(regions) != {'tumor', 'normal'}:
        raise ValueError(f"区域列必须仅包含 tumor 和 normal，但现为: {regions}")

    # 分别提取两组数据
    tumor_mask = df['region'] == 'tumor'
    normal_mask = df['region'] == 'normal'
    data_tumor = df[tumor_mask].drop(columns=['sample', 'region'])
    data_normal = df[normal_mask].drop(columns=['sample', 'region'])

    if data_tumor.shape[0] < MIN_SAMPLES_PER_REGION or data_normal.shape[0] < MIN_SAMPLES_PER_REGION:
        print("  任一组样本数不足，跳过统计检验")
        return None

    results = []
    for ct in data_tumor.columns:
        v1 = data_tumor[ct].values
        v2 = data_normal[ct].values
        if np.array_equal(v1, v2):
            pval = 1.0
        else:
            try:
                _, pval = mannwhitneyu(v1, v2, alternative='two-sided')
            except ValueError:
                pval = 1.0

        med_t = np.median(v1) if len(v1) > 0 else 0
        med_n = np.median(v2) if len(v2) > 0 else 0
        pseudo = 1e-6
        log2fc = np.log2((med_t + pseudo) / (med_n + pseudo))

        results.append({
            'celltype': ct,
            'median_tumor': med_t,
            'median_normal': med_n,
            'log2FC': log2fc,
            'p_value': pval
        })

    res_df = pd.DataFrame(results)
    if len(res_df) > 0:
        _, fdr, _, _ = multipletests(res_df['p_value'], method='fdr_bh')
        res_df['FDR'] = fdr
        res_df['significant'] = res_df['FDR'] < FDR_THRESHOLD
    return res_df

def plot_one_celltype(prop_df, ct, res_row, out_dir):
    """绘制单个细胞类型的 tumor vs normal 箱线图 + 显著性标注"""
    plot_data = prop_df[['region', ct]].copy()
    # 确保两组都存在
    if plot_data['region'].nunique() != 2:
        return

    fig, ax = plt.subplots(figsize=(5, 6))
    sns.boxplot(data=plot_data, x='region', y=ct, order=['normal', 'tumor'],
                palette='Set2', ax=ax)
    sns.stripplot(data=plot_data, x='region', y=ct, order=['normal', 'tumor'],
                  color='black', size=5, alpha=0.6, ax=ax)

    if res_row is not None and res_row['significant']:
        pval = res_row['p_value']
        star_str = ''
        if pval < 0.001:
            star_str = '***'
        elif pval < 0.01:
            star_str = '**'
        elif pval < 0.05:
            star_str = '*'

        all_vals = plot_data[ct].values
        y_max = all_vals.max()
        y_range = all_vals.max() - all_vals.min() if all_vals.max() > all_vals.min() else 0.1
        y_line = y_max + y_range * 0.1
        h = y_range * 0.03

        x1, x2 = 0, 1
        ax.plot([x1, x1, x2, x2], [y_line, y_line + h, y_line + h, y_line],
                lw=1.5, c='black')
        ax.text((x1 + x2) / 2, y_line + h + y_range * 0.02,
                f"p = {pval:.4f} {star_str}",
                ha='center', va='bottom', fontsize=11)

    ax.set_ylabel('Proportion')
    ax.set_title(f'{ct}\nTumor vs Normal')
    plt.tight_layout()

    safe_ct = re.sub(r'[\\/*?:"<>| ]', '_', ct)
    fpath = os.path.join(out_dir, f"{safe_ct}.png")
    fig.savefig(fpath, dpi=300, bbox_inches='tight')
    plt.close(fig)

# ==================== 主流程 ====================
def main():
    print(f"读取数据: {H5AD_PATH}")
    adata = sc.read_h5ad(H5AD_PATH)

    # 检查必要列
    for col in ['region', 'sample', 'celltype']:
        if col not in adata.obs.columns:
            raise KeyError(f"obs 中缺少列: {col}")

    print(f"原始细胞数: {adata.n_obs}")
    print(f"Tumor 细胞数: {(adata.obs['region']=='tumor').sum()}")
    print(f"Normal 细胞数: {(adata.obs['region']=='normal').sum()}")

    # 计算 tumor / normal 区域比例
    prop_df = compute_region_proportions(adata)
    if prop_df.empty:
        print("没有满足细胞数要求的样本-区域组合，退出。")
        return

    print(f"满足阈值的 (sample, region) 组合数: {prop_df.shape[0]}")
    print(f"各区域样本数:\n{prop_df['region'].value_counts().to_string()}")

    # 创建输出目录
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 差异检验
    res_df = run_differential_test(prop_df)
    if res_df is None:
        print("统计检验未执行。")
        return

    # 保存检验结果
    res_df.to_csv(os.path.join(OUTPUT_DIR, "differential_abundance.csv"), index=False)
    sig_count = res_df['significant'].sum()
    print(f"显著细胞类型数 (FDR<{FDR_THRESHOLD}): {sig_count}")

    # 为每个细胞类型绘制箱线图
    celltypes = prop_df.columns.drop(['sample', 'region'])
    for ct in celltypes:
        row = res_df[res_df['celltype'] == ct]
        row_data = row.iloc[0] if not row.empty else None
        plot_one_celltype(prop_df, ct, row_data, OUTPUT_DIR)

    print(f"完成！箱线图与统计表已保存至 {OUTPUT_DIR}/")

if __name__ == "__main__":
    main()
