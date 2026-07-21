#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1. 将 all_with_group.csv 中的通路评分按 cell_id + sample 匹配添加到 h5ad 的 obs 中，保存为新 h5ad。
2. 提取所有上皮细胞（基于 obs 中的 celltype 列）。
3. 对上皮细胞表达矩阵进行标准化（normalize_total + log1p）。
4. 对每个通路-基因对，计算标准化后基因表达量与对应通路评分的 Spearman 相关系数，输出 cor.csv。
"""

import pandas as pd
import numpy as np
import scanpy as sc
from scipy.stats import spearmanr
from scipy.sparse import issparse
import warnings
warnings.filterwarnings('ignore')

# ======================== 路径配置 ========================
H5AD_INPUT = "All_filter_state_samp_new_with_CNs10-15.h5ad"
CSV_SCORE = "all_with_group.csv"
GENE_MAP_FILE = "gene_func.csv"          # 通路-基因映射文件，应包含 Function, Pathway, Gene 三列
H5AD_OUTPUT = "All_filter_state_samp_new_with_CNs10-15_score.h5ad"
COR_OUTPUT = "cor.csv"

# ======================== 上皮细胞类型列表 ========================
EPI_TYPES = [
    "CXCL5+ Inf stem tumor Epi",
    "IGF2+ prolif stem tumor Epi",
    "Hypo Stem tumor Epi",
    "WNT+ ISC",
    "ESC-like Epi",
    "Prol Metab Epi",
    "EMT-like Inv Ent",
    "Inf Hyp Epi",
    "Early stress Epi",
    "Metab stress Epi",
    "Mat Ent",
    "MUC5B+ CGC",
    "MUC5AC+ MGC",
    "PCSK1+ AGC",
    "Hypo Metab Epi"
]

# ======================== 1. 读取数据 ========================
print("1. 正在读取 h5ad 文件...")
adata = sc.read_h5ad(H5AD_INPUT)
print(f"h5ad 形状: {adata.shape}")

print("2. 正在读取通路评分 CSV...")
df_score = pd.read_csv(CSV_SCORE)
print(f"CSV 评分数据形状: {df_score.shape}")

# 评分列名列表
score_cols = [
    "DNA_repair_score", "EMT_score", "angiogenesis_score",
    "cellular_metabolism_score", "cycling_score",
    "dedifferentiation_score", "hypoxia_score",
    "inflammation_score", "oncogenic_TFs_score", "stemness_score"
]

# ======================== 2. 匹配评分并添加到 adata.obs ========================
# 确保 cell_id 和 sample 类型一致
df_score['cell_id'] = df_score['cell_id'].astype(str)
df_score['sample'] = df_score['sample'].astype(str)
adata.obs['cell_id'] = adata.obs['cell_id'].astype(str)
adata.obs['sample'] = adata.obs['sample'].astype(str)

# 构建待合并的评分子表，重命名 group 列避免冲突
df_merge = df_score[['cell_id', 'sample'] + score_cols + ['group']].copy()
df_merge.rename(columns={'group': 'group_from_csv'}, inplace=True)

# ---- 关键修复1：去除重复的 (cell_id, sample) ----
dup_mask = df_merge.duplicated(subset=['cell_id', 'sample'])
if dup_mask.any():
    print(f"警告: 评分数据中有 {dup_mask.sum()} 个重复的 (cell_id, sample) 记录，将保留第一个。")
    df_merge = df_merge.drop_duplicates(subset=['cell_id', 'sample'], keep='first')

# ---- 关键修复2：确保评分列是浮点数，group列是字符串 ----
for col in score_cols:
    df_merge[col] = pd.to_numeric(df_merge[col], errors='coerce')
df_merge['group_from_csv'] = df_merge['group_from_csv'].astype(str)

# 左连接：将评分添加到 adata.obs（不改变行数）
print("3. 正在按 cell_id + sample 匹配评分...")
adata.obs = adata.obs.merge(df_merge, on=['cell_id', 'sample'], how='left')

# 检查匹配率
matched = adata.obs[score_cols].notna().all(axis=1).sum()
print(f"成功匹配到所有评分的细胞数: {matched} / {adata.n_obs}")

# ---- 关键修复3：再次确保合并后的评分列是数值类型，避免写入h5ad时报错 ----
for col in score_cols:
    if col in adata.obs.columns:
        adata.obs[col] = pd.to_numeric(adata.obs[col], errors='coerce').astype(float)

if 'group_from_csv' in adata.obs.columns:
    adata.obs['group_from_csv'] = adata.obs['group_from_csv'].astype(str)

# ======================== 3. 保存添加评分后的 h5ad ========================
print("4. 正在保存新 h5ad...")
adata.write(H5AD_OUTPUT)
print(f"已保存至: {H5AD_OUTPUT}")

# ======================== 4. 提取上皮细胞 ========================
print("5. 正在提取上皮细胞...")
if 'celltype' not in adata.obs.columns:
    raise KeyError("h5ad 的 obs 中缺少 'celltype' 列，无法筛选上皮细胞。")

mask_epi = adata.obs['celltype'].isin(EPI_TYPES)
adata_epi = adata[mask_epi].copy()
print(f"上皮细胞数量: {adata_epi.n_obs}")

# ======================== 5. 标准化上皮细胞表达矩阵 ========================
print("6. 标准化上皮细胞表达矩阵 (normalize_total + log1p)...")
# 假设原始矩阵为计数矩阵；若已标准化可注释掉下面两行
sc.pp.normalize_total(adata_epi, target_sum=1e4)   # 每细胞归一化到 10,000
sc.pp.log1p(adata_epi)                             # log(1 + x) 转换

# ======================== 6. 读取通路-基因映射 ========================
print("7. 正在读取通路-基因映射文件...")

# 尝试自动检测分隔符（制表符或逗号）
try:
    # 首先尝试制表符
    df_map = pd.read_csv(GENE_MAP_FILE, sep='\t')
    if df_map.shape[1] < 2 or 'Function' not in df_map.columns:
        raise ValueError("制表符读取失败，列缺失")
except:
    try:
        # 其次尝试多个空格
        df_map = pd.read_csv(GENE_MAP_FILE, sep='\s+', engine='python')
    except:
        # 默认逗号
        df_map = pd.read_csv(GENE_MAP_FILE)

# 统一列名：如果没有 'Function', 'Pathway', 'Gene'，尝试根据位置赋值
if df_map.shape[1] == 3 and not all(col in df_map.columns for col in ['Function', 'Pathway', 'Gene']):
    df_map.columns = ['Function', 'Pathway', 'Gene']
elif 'Function' not in df_map.columns:
    raise KeyError(f"映射文件中缺少 'Function' 列，当前列名: {df_map.columns.tolist()}。请检查文件格式。")

print(f"映射文件形状: {df_map.shape}")

# 构建从 Function (通路简称) 到评分列名的对应关系
func_to_score = {col.replace('_score', ''): col for col in score_cols}
# 筛选出在评分列中存在的通路
valid_func = set(func_to_score.keys()) & set(df_map['Function'].unique())
print(f"有效通路数: {len(valid_func)}")

# ======================== 7. 计算相关性 ========================
print("8. 正在计算基因表达与通路评分的相关性（Spearman）...")

# 获取上皮细胞的表达矩阵（转为稠密数组，注意内存；若过大可采用稀疏批量处理）
X = adata_epi.X
if issparse(X):
    X = X.toarray()

# 基因列表
genes = adata_epi.var_names.tolist()
gene_to_idx = {g: i for i, g in enumerate(genes)}

results = []
for func in valid_func:
    score_col = func_to_score[func]
    if score_col not in adata_epi.obs.columns:
        print(f"警告: 上皮细胞 obs 中无列 {score_col}，跳过通路 {func}")
        continue

    scores = adata_epi.obs[score_col].values.astype(float)
    # 剔除评分为 NaN 的细胞
    valid_mask = ~np.isnan(scores)
    if valid_mask.sum() < 5:
        print(f"警告: {func} 的有效样本数不足5，跳过")
        continue

    scores_valid = scores[valid_mask]
    X_valid = X[valid_mask, :]

    # 获取该通路对应的基因
    genes_in_func = df_map[df_map['Function'] == func]['Gene'].unique()
    for gene in genes_in_func:
        if gene not in gene_to_idx:
            continue
        gidx = gene_to_idx[gene]
        expr = X_valid[:, gidx]

        # 排除基因表达恒定（方差为零）的情况
        if np.max(expr) == np.min(expr):
            continue

        try:
            rho, pval = spearmanr(expr, scores_valid)
        except Exception as e:
            print(f"计算 {func} - {gene} 时出错: {e}")
            continue

        results.append({
            'Function': func,
            'Gene': gene,
            'Spearman_r': rho,
            'p_value': pval
        })

# ======================== 8. 输出结果 ========================
print(f"9. 共计算 {len(results)} 对相关性，正在保存至 {COR_OUTPUT}...")
df_cor = pd.DataFrame(results)
df_cor.to_csv(COR_OUTPUT, index=False)
print("全部完成！")
