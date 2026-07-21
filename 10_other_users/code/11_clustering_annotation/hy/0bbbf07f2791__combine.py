import scanpy as sc
import anndata as ad
import numpy as np
import os

def merge_h5ad_files(file_list, output_file="merged_data.h5ad"):
    """
    合并多个h5ad文件并保留原始数据
    
    参数:
    file_list: 要合并的h5ad文件路径列表
    output_file: 输出的合并文件路径
    """
    
    adata_list = []
    
    print("开始读取h5ad文件...")
    for i, file_path in enumerate(file_list):
        print(f"读取文件 {i+1}/{len(file_list)}: {os.path.basename(file_path)}")
        
        # 读取h5ad文件
        adata = sc.read_h5ad(file_path)
        
        # 确保每个数据集都有原始计数数据
        if adata.raw is None:
            # 如果没有raw数据，将当前数据设为raw
            adata.raw = adata
        else:
            # 确保raw数据是正确的格式
            adata.raw = adata.raw.to_adata() if hasattr(adata.raw, 'to_adata') else adata.raw
        
        # 为每个细胞添加来源标识
        adata.obs['batch'] = os.path.basename(file_path).replace('.h5ad', '')
        
        adata_list.append(adata)
    
    print("所有文件读取完成，开始合并...")
    
    if len(adata_list) == 0:
        raise ValueError("没有找到可合并的数据集")
    
    # 方法1: 尝试直接合并（适用于相同基因集）
    try:
        print("尝试基于共同基因进行合并...")
        
        # 找出所有数据集的共同基因
        common_genes = set(adata_list[0].var_names)
        for adata in adata_list[1:]:
            common_genes = common_genes.intersection(set(adata.var_names))
        
        print(f"发现 {len(common_genes)} 个共同基因")
        
        # 如果共同基因数量合理，使用共同基因子集进行合并
        if len(common_genes) > 1000:  # 要有1000个共同基因
            # 对每个数据集取共同基因子集
            adata_subset_list = []
            for adata in adata_list:
                adata_subset = adata[:, list(common_genes)].copy()
                adata_subset_list.append(adata_subset)
            
            # 合并数据
            merged_adata = ad.concat(adata_subset_list, join='inner', label='batch_key')
            print("使用共同基因合并成功")
            
        else:
            # 共同基因太少，使用外连接合并
            print("共同基因较少，使用外连接合并...")
            merged_adata = ad.concat(adata_list, join='outer', label='batch_key')
            
    except Exception as e:
        print(f"方法1失败: {e}, 尝试外连接合并...")
        # 方法2: 外连接合并（处理不同基因集）
        merged_adata = ad.concat(adata_list, join='outer', label='batch_key')
    
    # 处理合并后的raw数据
    print("处理原始计数数据...")
    raw_data_list = []
    for adata in adata_list:
        raw_adata = adata.raw.to_adata() if hasattr(adata.raw, 'to_adata') else adata.raw
        # 确保只取共同基因（如果之前用了共同基因）
        if 'common_genes' in locals() and len(common_genes) > 1000:
            common_genes_list = list(common_genes.intersection(set(raw_adata.var_names)))
            raw_adata = raw_adata[:, common_genes_list]
        raw_data_list.append(raw_adata)
    
    # 合并raw数据
    try:
        merged_raw_adata = ad.concat(raw_data_list, join='inner' if 'common_genes' in locals() and len(common_genes) > 1000 else 'outer')
        merged_adata.raw = merged_raw_adata
        print("原始计数数据合并成功")
    except Exception as e:
        print(f"原始数据合并失败: {e}, 将使用当前数据作为原始数据")
        merged_adata.raw = merged_adata
    
    # 添加基本的质量控制指标
    print("计算质量控制指标...")
    sc.pp.calculate_qc_metrics(merged_adata, inplace=True)
    
    # 保存合并后的数据
    print(f"保存合并数据到: {output_file}")
    merged_adata.write_h5ad(output_file)
    
    # 输出合并结果统计
    print("\n=== 合并完成 ===")
    print(f"总细胞数: {merged_adata.n_obs}")
    print(f"总基因数: {merged_adata.n_var}")
    print(f"批次信息: {merged_adata.obs['batch'].value_counts().to_dict()}")
    
    return merged_adata

# 您的文件列表
file_list = [
    'B.h5ad',
    'celltype_Pericytes.h5ad', 
    'celltype_Plasmas.h5ad',
    'Epithelials.h5ad',
    'Fibro.h5ad',
    'Macrophages.h5ad',
    'T.h5ad'
]

# 确保文件存在
existing_files = [f for f in file_list if os.path.exists(f)]
missing_files = set(file_list) - set(existing_files)

if missing_files:
    print(f"警告: 以下文件不存在: {missing_files}")

if existing_files:
    # 执行合并
    merged_data = merge_h5ad_files(existing_files, "merged_single_cell_data.h5ad")
    
    print("\n合并成功！您可以使用以下代码验证结果:")
    print("""
# 验证合并结果
import scanpy as sc
adata = sc.read_h5ad("merged_single_cell_data.h5ad")
print(f"数据集形状: {adata.shape}")
print("批次分布:")
print(adata.obs['batch'].value_counts())

# 访问原始计数数据
if adata.raw is not None:
    print(f"原始数据形状: {adata.raw.shape}")
    """)
else:
    print("错误: 没有找到可用的h5ad文件")
