import scanpy as sc
import numpy as np

adata = sc.read_h5ad("integrated_data.h5ad") 

# 定义指定的细胞类型列表
target_celltypes = [
    "CD4+ CD247+T", "CD4+ Tem", "CD4+ Treg", "CD8+ Tex", "Disease Stem",
    "FAMS", "Kera Basal", "Kera Granular", "Kera Inflammatory", "Kera Spinous",
    "LAMP3+ DC", "M2", "MARCO+ Mø", "Normal Stem", "SBC Proliferation",
    "SBC Spinous", "SIGLEC1+ Mø", "SOCS3+ CTMCs", "SOCS3+ FB", "SOCS3+ M1",
    "THY1+ FB", "cDC1", "cDC2", "iDC", "\u200bCTMCs","POSTN+ FB"
]

# 创建原始数据的副本
new_adata = adata.copy()

# 初始化ct列为"others"
new_adata.obs['ct'] = "others"

# 遍历目标细胞类型
for ct in target_celltypes:
    # 检查细胞类型是否存在
    if ct in new_adata.obs['celltype_final'].values:
        # 创建布尔掩码选择该细胞类型
        mask = new_adata.obs['celltype_final'] == ct
        
        # 将符合条件的细胞设置为目标细胞类型名称
        new_adata.obs.loc[mask, 'ct'] = ct
    else:
        print(f"警告: 细胞类型 '{ct}' 在数据中未找到，已跳过")

# 输出细胞类型统计
print("ct列细胞类型分布:")
print(new_adata.obs['ct'].value_counts())

# 保存新数据到h5ad文件
output_path = "modified_data_with_ct.h5ad"
new_adata.write_h5ad(output_path)
print(f"\n新数据已保存至: {output_path}")
