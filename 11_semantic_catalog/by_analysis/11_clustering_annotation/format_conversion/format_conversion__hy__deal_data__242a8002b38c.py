import scanpy as sc
import pandas as pd
import numpy as np

# 读取h5ad数据
adata = sc.read('All_filter_state_samp_new_with_CNs10-15.h5ad')

# 定义细胞类型到三大类的映射规则
epithelial_celltypes = [
    'CXCL5+ Inf stem tumor Epi',
    'IGF2+ prolif stem tumor Epi',
    'Hypo Stem tumor Epi',
    'WNT+ ISC',
    'ESC-like Epi',
    'Prol Metab Epi',
    'EMT-like Inv Ent',
    'Inf Hyp Epi',
    'Early stress Epi',
    'Metab stress Epi',
    'Mat Ent',
    'MUC5B+ CGC',
    'MUC5AC+ MGC',
    'PCSK1+ AGC',
    'Hypo Metab Epi'
]

immune_celltypes = [
    'CTSL+ PLAUR+ Macro',
    'STAB1+ Macro',
    'CD4+ T',
    'γδ T',
    'CD8+ T',
    'Treg',
    'MKI67+ T',
    'DC',
    'GCB',
    'MBC',
    'Mast',
    'Plasmas',
    'Neutrophils'
]

stromal_celltypes = [
    'Angiogenic EC',
    'Venule EC',
    'Lymphatic EC',
    'Arteriole EC',
    'ADAMDEC1+ Homeo-Immuno Fib',
    'Act Strom CAF',
    'matri-prod CAF',
    'RSPO3+ GREM1+ Crypt Fib',
    'WNT5B+ VEGFA+ Villus Fib',
    'Myofib',
    'apCAF',
    'Com Fib',
    'CXCL12+ Fib',
    'Pericyte-like Fib'
]

# 创建新列并分配类别
adata.obs['celltype_zq'] = np.nan  # 初始化为NaN

# 假设原细胞类型存储在adata.obs['celltype']中
# 如果不是，请将下面的'celltype'替换为实际的列名
celltype_col = 'celltype'  # 请根据实际情况修改

# 检查原细胞类型列是否存在
if celltype_col not in adata.obs.columns:
    print(f"警告: 列 '{celltype_col}' 不存在!")
    print("可用的列有:", list(adata.obs.columns))
    # 尝试找到可能的细胞类型列
    possible_cols = [col for col in adata.obs.columns if 'cell' in col.lower() or 'type' in col.lower()]
    if possible_cols:
        print(f"可能的细胞类型列: {possible_cols}")
else:
    # 分配上皮细胞类型
    for ct in epithelial_celltypes:
        mask = adata.obs[celltype_col] == ct
        adata.obs.loc[mask, 'celltype_zq'] = 'Epithelial'
    
    # 分配免疫细胞类型
    for ct in immune_celltypes:
        mask = adata.obs[celltype_col] == ct
        adata.obs.loc[mask, 'celltype_zq'] = 'Immune'
    
    # 分配基质细胞类型
    for ct in stromal_celltypes:
        mask = adata.obs[celltype_col] == ct
        adata.obs.loc[mask, 'celltype_zq'] = 'Stromal'
    
    # 检查是否有未分配的细胞类型
    unassigned = adata.obs['celltype_zq'].isna().sum()
    if unassigned > 0:
        print(f"警告: 有 {unassigned} 个细胞未分配类别")
        # 显示未分配的细胞类型
        unassigned_types = adata.obs.loc[adata.obs['celltype_zq'].isna(), celltype_col].unique()
        print(f"未分配的细胞类型: {list(unassigned_types)}")
    
    # 统计各类别数量
    print("细胞类型分布:")
    print(adata.obs['celltype_zq'].value_counts(dropna=False))
    
    # 保存新数据
    output_file = 'All_filter_state_samp_new_ctzq_with_CNs10-15.h5ad'
    adata.write(output_file)
    print(f"数据已保存为: {output_file}")
