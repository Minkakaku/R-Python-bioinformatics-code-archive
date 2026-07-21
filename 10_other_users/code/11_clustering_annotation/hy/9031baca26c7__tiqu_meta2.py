import scanpy as sc
import pandas as pd
import os

# 读取h5ad文件
adata = sc.read('cell_type_Adipocytes.h5ad')

# 查看数据基本信息
print(f"数据形状: {adata.shape}")
print(f"观测变量(obs)列名: {list(adata.obs.columns)}")
print(f"前5个细胞的obs信息:")
print(adata.obs.head())

# 检查是否存在所需的列
if 'sample' not in adata.obs.columns:
    print("警告: 未找到'sample'列，以下是可用的列:")
    print(list(adata.obs.columns))
    # 尝试寻找可能的样本列
    sample_candidates = [col for col in adata.obs.columns if 'sample' in col.lower() or 'group' in col.lower() or 'condition' in col.lower()]
    if sample_candidates:
        print(f"可能的样本列: {sample_candidates}")
        # 使用第一个候选列
        sample_col = sample_candidates[0]
    else:
        # 如果没有找到，使用默认值
        print("使用默认值'unknown'作为样本")
        adata.obs['sample'] = 'unknown'
        sample_col = 'sample'
else:
    sample_col = 'sample'

if 'celltype' not in adata.obs.columns:
    print("警告: 未找到'celltype'列，以下是可用的列:")
    print(list(adata.obs.columns))
    # 尝试寻找可能的细胞类型列
    celltype_candidates = [col for col in adata.obs.columns if 'cell' in col.lower() and 'type' in col.lower()]
    if celltype_candidates:
        print(f"可能的细胞类型列: {celltype_candidates}")
        # 使用第一个候选列
        celltype_col = celltype_candidates[0]
    else:
        # 如果没有找到，创建一个默认列
        print("使用默认值'unknown_celltype'作为细胞类型")
        adata.obs['celltype'] = 'unknown_celltype'
        celltype_col = 'celltype'
else:
    celltype_col = 'cell_type'

# 创建包含所需信息的DataFrame
meta_data = pd.DataFrame({
    'cell_index': adata.obs.index,  # 细胞的index
    'sample': adata.obs[sample_col],
    'celltype': adata.obs[celltype_col]
})

# 保存为CSV文件
output_file = 'Adipocytes_meta.csv'
meta_data.to_csv(output_file, index=False)

print(f"\n提取的信息:")
print(f"总细胞数: {meta_data.shape[0]}")
print(f"样本数量: {meta_data['sample'].nunique()}")
print(f"细胞类型数量: {meta_data['celltype'].nunique()}")
print(f"\n样本分布:")
print(meta_data['sample'].value_counts())
print(f"\n细胞类型分布:")
print(meta_data['celltype'].value_counts())
print(f"\n数据已保存为: {output_file}")

# 可选：显示前几行数据
print(f"\nCSV文件前5行:")
print(meta_data.head())
