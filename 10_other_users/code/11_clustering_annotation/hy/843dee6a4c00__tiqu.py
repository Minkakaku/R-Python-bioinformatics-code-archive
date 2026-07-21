#!/usr/bin/env python3
"""
从 All_filter_state_samp_new_with_CNs10-15.h5ad 中提取指定样本。
"""

import anndata as ad

# 输入文件
input_file = "All_filter_state_samp_new_with_CNs10-15.h5ad"
output_file = "data_4_band.h5ad"

# 要提取的样本列表
samples_to_keep = [
    "009D-1",
    "SR202500003B",
    "SR202500005D",
    "SR202500005F-1",
    "SR202500009D-1",
    "SR202500009D-2",
    "SR2025000016J",
    "SR2025000016K",
]

# 读取数据
print(f"正在读取 {input_file}...")
adata = ad.read_h5ad(input_file)
print(f"总细胞数: {adata.n_obs}")

# 检查样本列是否存在
if "sample" not in adata.obs.columns:
    raise KeyError("obs 中缺少 'sample' 列，请检查数据。")

# 筛选存在的样本
existing_samples = set(adata.obs["sample"].unique())
missing = [s for s in samples_to_keep if s not in existing_samples]
if missing:
    print(f"警告: 以下样本在原数据中不存在，将被忽略: {missing}")

# 进行子集提取
mask = adata.obs["sample"].isin(samples_to_keep)
adata_sub = adata[mask].copy()

print(f"提取到的细胞数: {adata_sub.n_obs}")
print("样本分布:")
print(adata_sub.obs["sample"].value_counts())

# 保存
adata_sub.write_h5ad(output_file)
print(f"已保存为 {output_file}")
