#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CD8+ T细胞功能状态基因集评分（更新版，适配 Xenium 5k）
输入: All_filter.h5ad (原始计数矩阵)
输出: All_filter_score.h5ad (添加评分列)
"""

import scanpy as sc
import warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# ---------- 1. 加载数据 ----------
print(">>> 读取 All_filter.h5ad ...")
adata = sc.read_h5ad("All_filter.h5ad")
print(f"数据维度: {adata.n_obs} 细胞/斑点 × {adata.n_vars} 基因")

# ---------- 2. 最新基因集定义 ----------
gene_sets_raw = {
    "Effector_CD8T_5k": [
        "GZMB", "PRF1", "IFNG", "GZMA", "GZMM",
        "CST7", "XCL1", "TBX21", "CX3CR1", "ZNF683"
    ],
    "Proliferating_CD8T_5k": [
        "MKI67", "PCNA", "TOP2A", "STMN1", "MCM5",
        "CCNB2", "CDK1", "AURKB", "BIRC5"
    ],
    "TerminalExhausted_CD8T_5k": [
        "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT",
        "TOX", "EOMES", "BATF", "PRDM1", "CD38", "CXCL13"
    ],
    "IFN_Response_5k": [
        "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1",
        "OASL", "OAS3", "RSAD2", "STAT1", "STAT2",
        "IRF7", "DDX58"
    ],
    "Hypoxia_Stress_5k": [
        "HIF1A", "BNIP3", "EGLN3", "NR4A1", "EGR1"
    ]
}

# ---------- 3. 自动过滤缺失基因 ----------
print("\n>>> 检查基因覆盖并过滤：")
filtered_sets = {}
for name, genes in gene_sets_raw.items():
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        print(f"⚠️ {name} 缺失 {len(missing)} 个基因: {missing}")
    if len(present) < 3:
        print(f"❌ {name} 剩余基因不足3个，跳过评分。")
        continue
    filtered_sets[name] = present
    print(f"✅ {name} 将使用 {len(present)} 个基因: {present}")

if not filtered_sets:
    raise RuntimeError("没有可用的基因集，请检查基因名称！")

# ---------- 4. 标准化原始计数 ----------
print("\n>>> 预处理：normalize_total (1e4) -> log1p")
adata.raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ---------- 5. 计算评分 ----------
print("\n>>> 计算基因集评分 ...")
for name, genes in filtered_sets.items():
    ctrl = min(50, len(genes))      # 对照基因数不超过基因集大小
    print(f"  计算 {name} (ctrl_size={ctrl})")
    sc.tl.score_genes(
        adata,
        gene_list=genes,
        score_name=f"{name}_score",
        ctrl_size=ctrl
    )

new_cols = [f"{n}_score" for n in filtered_sets.keys()]
print(f"\n已在 .obs 中添加: {new_cols}")

# ---------- 6. 保存 ----------
output_file = "All_filter_score.h5ad"
print(f"\n>>> 保存到 {output_file} ...")
adata.write(output_file)
print("完成！")
