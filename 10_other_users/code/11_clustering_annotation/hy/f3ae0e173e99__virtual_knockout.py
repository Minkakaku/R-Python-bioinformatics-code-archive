#!/usr/bin/env python3
"""
基于 scTenifold 的虚拟敲除脚本 - 最终版
使用 virtual_knockout 并传入完整的质控、网络构建、张量分解、流形对齐参数
数据读取部分直接调用 scTenifold 的 anndata_to_dataframe
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.sparse import issparse
from scTenifold import virtual_knockout
from scTenifold.core._networks import anndata_to_dataframe

def load_adata(adata_file):
    """加载 AnnData 数据"""
    adata = sc.read_h5ad(adata_file)
    if not adata.obs_names.is_unique:
        adata.obs_names_make_unique()
    print(f"AnnData 加载完成: {adata.n_obs} 细胞, {adata.n_vars} 基因")
    return adata

def extract_celltype(adata, celltype_col, target_celltype):
    """
    提取指定细胞类型的表达矩阵 (基因×细胞)
    使用 scTenifold 的 anndata_to_dataframe 进行转换
    """
    if celltype_col not in adata.obs.columns:
        raise ValueError(f"列 '{celltype_col}' 不在 adata.obs 中")

    adata_sub = adata[adata.obs[celltype_col] == target_celltype].copy()
    print(f"提取到 {adata_sub.n_obs} 个 '{target_celltype}' 细胞")

    if adata_sub.n_obs == 0:
        raise ValueError(f"细胞类型 '{target_celltype}' 不存在")

    # 直接使用 scTenifold 的转换函数，返回基因×细胞的 DataFrame
    expr_df = anndata_to_dataframe(adata_sub)
    print(f"表达矩阵维度 (基因×细胞): {expr_df.shape}")
    return expr_df

def visualize_results(result_df, output_prefix):
    """
    可视化虚拟敲除结果
    注意: 传入的 result_df 应该已经按 p-value 升序排列
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 1. 差异调控基因条形图 (Top 20，即 p-value 最小的20个)
    top20 = result_df.head(20)
    axes[0].barh(top20['Gene'][::-1], top20['Distance'][::-1], color='steelblue')
    axes[0].set_xlabel('Distance (差异调控距离)')
    axes[0].set_title('Top 20 显著差异调控基因 (按 p-value 排序)')

    # 2. 差异调控基因的 p-value 分布
    p_safe = result_df['p-value'].clip(lower=1e-300)  # 避免 log(0)
    axes[1].scatter(range(len(result_df)), -np.log10(p_safe),
                    alpha=0.6, s=10)
    axes[1].axhline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
    axes[1].set_xlabel('基因排名 (按 p-value 升序)')
    axes[1].set_ylabel('-log10(p-value)')
    axes[1].set_title('差异调控显著性分布')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_visualization.png", dpi=150)
    plt.show()
    print(f"可视化已保存至 {output_prefix}_visualization.png")

def main():
    parser = argparse.ArgumentParser(description="scTenifold 虚拟敲除脚本 (完整版)")
    parser.add_argument("--adata", type=str, required=True, help="AnnData 文件 (.h5ad)")
    parser.add_argument("--cell-type-col", type=str, default="celltype", help="adata.obs 中细胞类型列名")
    parser.add_argument("--cell-type", type=str, required=True, help="目标细胞类型")
    parser.add_argument("--gene", type=str, required=True, help="敲除的基因名")
    parser.add_argument("--output", type=str, default="ko_result", help="输出前缀")
    args = parser.parse_args()

    # 1. 加载数据
    adata = load_adata(args.adata)

    # 2. 提取指定细胞类型
    expr_df = extract_celltype(adata, args.cell-type-col, args.cell-type)

    # 检查基因是否存在
    if args.gene not in expr_df.index:
        raise ValueError(f"基因 '{args.gene}' 不在该细胞类型的表达矩阵中")

    # 3. 执行虚拟敲除（完整参数）
    print(f"执行虚拟敲除: 基因 '{args.gene}'")

    # ----- 参数设置（可根据数据调整）-----
    qc_kws = {
        "min_lib_size": 1,
        "plot": False,
        "max_mito_ratio": 1.0,   # 若已提前过滤，可设为1.0关闭；否则建议0.2
        "min_percent": 0,
        "min_exp_avg": 0,
        "min_exp_sum": 0
    }
    network_kws = {
        "n_nets": 2,             # 网络构建的子采样次数
        "n_samp_cells": 5000,      # 每次采样的细胞数（可根据细胞总数调整）
        "q": 0                   # 分位数阈值
    }
    td_kws = {
        "K": 2,                  # 张量分解的秩
        "max_iter": 20,
        "init": "random"
    }
    ma_kws = {
        "d": 2                   # 流形对齐维度
    }

    result = virtual_knockout(
        data=expr_df,
        ko_genes=[args.gene],
        qc_kws=qc_kws,
        network_kws=network_kws,
        ko_method="default",
        td_kws=td_kws,
        ma_kws=ma_kws
    )

    # 4. 按 p-value 排序并保存
    result = result.sort_values("p-value")  # 升序：显著的在前
    result.to_csv(f"{args.output}_knockout_results.csv", index=False)
    print(f"结果已保存至 {args.output}_knockout_results.csv")
    print("前5个显著差异基因:")
    print(result.head())

    # 5. 可视化结果
    visualize_results(result, args.output)

if __name__ == "__main__":
    main()
