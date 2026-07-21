#!/usr/bin/env python3
"""
filter_by_group_and_CN5.py

根据 group 列中的子字符串筛选细胞，并进一步只保留 CNs_10 == 'CN5' 的细胞。
添加 cellchat_group 列。
"""

import anndata as ad
import numpy as np

# 输入文件
input_h5ad = "merged_with_CNs_Her2_low.h5ad"

# 读取数据
adata = ad.read_h5ad(input_h5ad)
print(f"读取 {input_h5ad}: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")

# 检查必需的列
required_cols = ["group", "CNs_10"]
for col in required_cols:
    if col not in adata.obs.columns:
        raise KeyError(f"obs 中缺少 '{col}' 列，可用列: {', '.join(adata.obs.columns)}")

# 打印 group 和 CNs_10 的部分值以便确认
print("group 唯一值:", adata.obs["group"].unique().tolist())
print("CNs_10 唯一值:", adata.obs["CNs_10"].unique().tolist())

# ---------- 提取1：pre_pCR + pre_NpCR ----------
# 条件：group 包含 _pCR_P 或 _NpCR_P，且 CNs_10 == "CN5"
mask1_group = adata.obs["group"].str.contains("_pCR_P", na=False) | adata.obs["group"].str.contains("_NpCR_P", na=False)
mask1_cn = adata.obs["CNs_10"] == "CN5"
mask1 = mask1_group & mask1_cn

adata1 = adata[mask1].copy()
# 添加 cellchat_group
adata1.obs["cellchat_group"] = np.where(
    adata1.obs["group"].str.contains("_pCR_P", na=False),
    "pre_pCR",
    "pre_NpCR"  # 剩余为 _NpCR_P
)
print(f"提取1（pre_pCR + pre_NpCR，且 CNs_10=CN5）得到 {adata1.n_obs} 个细胞")
adata1.write("pre_pCR_NpCR.h5ad")
print("已保存 pre_pCR_NpCR.h5ad")

# ---------- 提取2：post_NpCR + pre_NpCR ----------
# 条件：group 包含 _NpCR_A 或 _NpCR_P，且 CNs_10 == "CN5"
mask2_group = adata.obs["group"].str.contains("_NpCR_A", na=False) | adata.obs["group"].str.contains("_NpCR_P", na=False)
mask2_cn = adata.obs["CNs_10"] == "CN5"
mask2 = mask2_group & mask2_cn

adata2 = adata[mask2].copy()
# 添加 cellchat_group
adata2.obs["cellchat_group"] = np.where(
    adata2.obs["group"].str.contains("_NpCR_A", na=False),
    "post_NpCR",
    "pre_NpCR"  # 对应 _NpCR_P
)
print(f"提取2（post_NpCR + pre_NpCR，且 CNs_10=CN5）得到 {adata2.n_obs} 个细胞")
adata2.write("all_pCR_NpCR.h5ad")
print("已保存 all_pCR_NpCR.h5ad")
