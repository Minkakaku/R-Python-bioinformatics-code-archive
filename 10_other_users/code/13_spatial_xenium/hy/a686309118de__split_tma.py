#!/usr/bin/env python3
"""
根据 sample 文件拆分 h5ad 数据。
每个 .sample 文件包含该 TMA 对应的样本名（每行一个），
脚本将在 adata.obs["sample"] 中匹配，并输出对应的子集 h5ad。
"""

import argparse
import os
import sys
import scanpy as sc
import pandas as pd

def load_sample_list(sample_file):
    """从 .sample 文件中读取样本名列表，忽略空行和首尾空白"""
    samples = []
    with open(sample_file, "r") as f:
        for line in f:
            s = line.strip()
            if s:  # 忽略空行
                samples.append(s)
    return samples

def main():
    parser = argparse.ArgumentParser(description="按 sample 文件拆分 h5ad")
    parser.add_argument("--h5ad", required=True, help="输入的 h5ad 文件，例如 All_filter.h5ad")
    parser.add_argument("--sample_files", nargs="+", required=True,
                        help="一个或多个 .sample 文件，例如 TMA1.sample TMA2.sample TMA.sample")
    parser.add_argument("--output_dir", default=".", help="输出目录，默认为当前目录")
    args = parser.parse_args()

    # 读取 h5ad
    print(f"读取 {args.h5ad} ...")
    adata = sc.read_h5ad(args.h5ad)

    if "sample" not in adata.obs.columns:
        print("错误：h5ad 的 obs 中没有 'sample' 列！")
        sys.exit(1)

    # 确保输出目录存在
    os.makedirs(args.output_dir, exist_ok=True)

    # 对每个 sample 文件进行拆分
    for sf in args.sample_files:
        if not os.path.exists(sf):
            print(f"警告：{sf} 不存在，跳过")
            continue

        print(f"处理 {sf} ...")
        sample_names = load_sample_list(sf)

        # 匹配样本
        mask = adata.obs["sample"].isin(sample_names)
        sub_adata = adata[mask].copy()

        n_total = adata.n_obs
        n_sub = sub_adata.n_obs
        print(f"  样本列表包含 {len(sample_names)} 个样本，匹配到 {n_sub} / {n_total} 个细胞")

        if n_sub == 0:
            print(f"  警告：未匹配到任何细胞，仍会输出空 annoata 对象")

        # 输出文件名：将 .sample 替换为 .h5ad
        basename = os.path.basename(sf)
        if basename.endswith(".sample"):
            out_name = basename[:-7] + ".h5ad"  # 去掉 .sample
        else:
            out_name = basename + ".h5ad"

        out_path = os.path.join(args.output_dir, out_name)
        sub_adata.write(out_path)
        print(f"  已保存至: {out_path}")

    print("完成！")

if __name__ == "__main__":
    main()
