import scanpy as sc
import numpy as np

# 1. 读取数据
adata = sc.read_h5ad('All_filter_state_new_samp_new.h5ad')  # 请替换为您的实际文件名

print(f"原始数据: {adata.shape[0]} 个细胞")

# 2. 提取目标细胞类型
target_celltypes = [
    "CXCL5+ Inf stem tumor Epi",
    "Prol Metab Epi", 
    "EMT-like Inv Ent",
    "ESC-like Epi",
    "WNT+ ISC",
    "MUC5B+ CGC",
    "Hypo Stem tumor Epi",
    "Hypo Metab Epi",
    "Inf Hyp Epi",
    "Early stress Epi",
    "Metab stress Epi",
    "MUC5AC+ MGC",
    "Mat Ent",
    "PCSK1+ AGC",
    "IGF2+ prolif stem tumor Epi"
]

# 3. 筛选细胞
mask = adata.obs['celltype'].isin(target_celltypes)
adata_epi = adata[mask, :].copy()

print(f"提取后数据: {adata_epi.shape[0]} 个细胞")

# 4. 保存为Epi_new.h5ad
adata_epi.write('Epi_new.h5ad')
print("已保存为 Epi_new.h5ad")
