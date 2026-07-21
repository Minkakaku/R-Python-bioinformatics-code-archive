#!/usr/bin/env python3
"""
绘制指定基因在特定细胞类型中的平均表达热图。
使用 raw 矩阵（已标准化和缩放），颜色：白→红。
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ---------- 参数 ----------
H5AD_PATH = "MPCvsTNBC_all_sample_ct.h5ad"
OUT_DIR = "./heatmap_expression"
os.makedirs(OUT_DIR, exist_ok=True)

# 目标基因（新增 LAMC1, HAS2, VEGFA）
TARGET_GENES = [
    "ANGPTL4","MKI67","S100A1","E2F1","CDK6","HIF1A","SOX9","COL2A1", "SDC2", "SNAI2","COL11A1","SQSTM1", "FGFR2", "TGFB2", "NTRK2", "BNIP3",
    "LAMC1", "HAS2", "VEGFA"
]

TARGET_CELLTYPES = [
    "FGFR2+_Basal_Epithelium",
    "Stress-responsive_glandular_Epithelium",
    "EMT_Epithelium",
    "Proliferative_Epithelium"
]

OUT_PREFIX = "target_gene_expression_heatmap"

# ---------- 1. 读取数据 ----------
print(f"Reading {H5AD_PATH} ...")
adata = sc.read_h5ad(H5AD_PATH)

if adata.raw is None:
    raise ValueError("AnnData object does not contain 'raw'. Please use raw matrix.")
print("Using raw matrix (normalized & scaled data).")

if 'annotation' not in adata.obs.columns:
    raise KeyError("'annotation' column missing in adata.obs.")

# 筛选细胞类型
cell_mask = adata.obs['annotation'].isin(TARGET_CELLTYPES)
adata_sub = adata[cell_mask, :].copy()
print(f"Filtered to {adata_sub.n_obs} cells.")

# 检查目标基因是否在 raw.var_names 中
raw_var = adata.raw.var_names
missing = [g for g in TARGET_GENES if g not in raw_var]
if missing:
    print(f"Warning: Genes not found in raw: {missing}")
    TARGET_GENES = [g for g in TARGET_GENES if g in raw_var]
if not TARGET_GENES:
    raise ValueError("No target genes found in raw data.")

# 从 raw 中提取对应细胞的表达矩阵
cell_barcodes = adata_sub.obs_names
raw_sub = adata.raw[cell_barcodes, TARGET_GENES]  # 返回 AnnData
expr_matrix = raw_sub.X

# 转为稠密矩阵（若数据量不大）
if hasattr(expr_matrix, "toarray"):
    expr_matrix = expr_matrix.toarray()

expr_df = pd.DataFrame(expr_matrix, index=cell_barcodes, columns=TARGET_GENES)
expr_df['celltype'] = adata_sub.obs['annotation'].values

# ---------- 2. 计算每种细胞类型的平均表达 ----------
mean_expr = expr_df.groupby('celltype')[TARGET_GENES].mean()
mean_expr = mean_expr.reindex(TARGET_CELLTYPES)
print("Average expression per cell type:")
print(mean_expr.round(3))

# ---------- 3. 绘制热图（白色→红色）----------
# 自定义 colormap：白色到红色
white_red_cmap = sns.light_palette("red", as_cmap=True)

plt.figure(figsize=(len(TARGET_GENES)*0.8 + 3, len(TARGET_CELLTYPES)*0.6 + 1))
ax = sns.heatmap(
    mean_expr,
    annot=True,
    fmt='.2f',
    cmap=white_red_cmap,
    linewidths=0.5,
    cbar_kws={'label': 'Mean expression (scaled)'}
)
ax.set_title('Average expression of target genes in selected cell types')
ax.set_xlabel('Gene')
ax.set_ylabel('Cell Type')
plt.tight_layout()

out_png = os.path.join(OUT_DIR, f"{OUT_PREFIX}.png")
out_pdf = os.path.join(OUT_DIR, f"{OUT_PREFIX}.pdf")
plt.savefig(out_png, dpi=300, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
plt.close()
print(f"Heatmap saved: {out_png}, {out_pdf}")

out_csv = os.path.join(OUT_DIR, f"{OUT_PREFIX}_mean_expr.csv")
mean_expr.to_csv(out_csv)
print(f"Mean expression table saved: {out_csv}")
