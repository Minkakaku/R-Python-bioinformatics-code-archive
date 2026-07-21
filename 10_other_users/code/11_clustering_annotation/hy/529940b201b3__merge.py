import scanpy as sc
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

def merge_epithelial_h5ad_files(file_list, output_path='merged_epithelial_cells.h5ad'):
    """
    合并多个上皮细胞类型的h5ad文件，保留raw矩阵
    
    Args:
        file_list: h5ad文件路径列表
        output_path: 合并后输出文件路径
    """
    print(f"开始合并 {len(file_list)} 个h5ad文件...")
    
    adata_list = []
    raw_shapes = []
    
    # 1. 读取所有h5ad文件
    for i, file_path in enumerate(file_list):
        print(f"\n读取文件 {i+1}: {file_path}")
        
        if not os.path.exists(file_path):
            print(f"  错误: 文件不存在")
            continue
            
        try:
            adata = sc.read_h5ad(file_path)
            print(f"  细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
            
            # 检查是否有raw矩阵
            if adata.raw is not None:
                raw_adata = adata.raw.to_adata()
                print(f"  raw矩阵形状: {raw_adata.shape}")
                raw_shapes.append(raw_adata.shape)
            else:
                print(f"  警告: 文件没有raw矩阵")
            
            # 添加来源信息
            source_name = os.path.basename(file_path).replace('.h5ad', '')
            adata.obs['source_file'] = source_name
            
            adata_list.append(adata)
            
        except Exception as e:
            print(f"  读取失败: {e}")
    
    if len(adata_list) == 0:
        print("没有成功读取任何文件")
        return None
    
    print(f"\n成功读取 {len(adata_list)} 个文件")
    
    # 2. 检查基因是否一致
    print("\n检查基因一致性...")
    
    # 检查主数据的基因
    all_var_names = [adata.var_names for adata in adata_list]
    common_genes = set(all_var_names[0])
    
    for i in range(1, len(all_var_names)):
        common_genes = common_genes.intersection(set(all_var_names[i]))
    
    print(f"  主数据共同基因数: {len(common_genes)}")
    
    # 检查raw数据的基因
    raw_genes_list = []
    for adata in adata_list:
        if adata.raw is not None:
            raw_adata = adata.raw.to_adata()
            raw_genes_list.append(set(raw_adata.var_names))
    
    if raw_genes_list:
        common_raw_genes = set.intersection(*raw_genes_list)
        print(f"  raw数据共同基因数: {len(common_raw_genes)}")
    else:
        common_raw_genes = set()
        print(f"  没有找到raw矩阵")
    
    # 3. 对齐基因（如果需要）
    aligned_adata_list = []
    
    for adata in adata_list:
        # 确保主数据使用共同基因
        if len(common_genes) > 0 and len(common_genes) < len(adata.var_names):
            print(f"  对齐基因: 从 {len(adata.var_names)} 到 {len(common_genes)} 个基因")
            adata_aligned = adata[:, list(common_genes)].copy()
        else:
            adata_aligned = adata.copy()
        
        aligned_adata_list.append(adata_aligned)
    
    # 4. 合并主数据
    print(f"\n合并主数据...")
    
    # 使用scanpy的concat函数合并
    merged_adata = sc.concat(
        aligned_adata_list,
        join='outer',  # 使用outer join保留所有列
        label='source_file',  # 添加来源标签
        index_unique='-'  # 使索引唯一
    )
    
    print(f"  合并后主数据形状: {merged_adata.shape}")
    print(f"  细胞来源分布:")
    print(merged_adata.obs['source_file'].value_counts())
    
    # 5. 合并raw矩阵
    print(f"\n合并raw矩阵...")
    
    # 检查哪些adata有raw矩阵
    raw_adata_list = []
    for adata in adata_list:
        if adata.raw is not None:
            raw_adata = adata.raw.to_adata()
            
            # 确保使用共同的raw基因
            if len(common_raw_genes) > 0 and len(common_raw_genes) < len(raw_adata.var_names):
                raw_adata = raw_adata[:, list(common_raw_genes)].copy()
            
            raw_adata_list.append(raw_adata)
    
    if raw_adata_list:
        # 合并raw数据
        merged_raw_adata = sc.concat(
            raw_adata_list,
            join='outer',
            index_unique='-'
        )
        
        print(f"  合并后raw矩阵形状: {merged_raw_adata.shape}")
        
        # 设置raw属性
        merged_adata.raw = merged_raw_adata
    else:
        print("  警告: 没有可合并的raw矩阵")
    
    # 6. 保存合并后的数据
    print(f"\n保存合并数据到: {output_path}")
    merged_adata.write(output_path)
    
    # 7. 验证保存的数据
    print(f"\n验证保存的数据...")
    try:
        loaded_data = sc.read_h5ad(output_path)
        print(f"  加载成功: 主数据形状={loaded_data.shape}")
        
        if loaded_data.raw is not None:
            loaded_raw = loaded_data.raw.to_adata()
            print(f"  raw矩阵形状={loaded_raw.shape}")
            
            # 检查细胞是否匹配
            if loaded_data.n_obs == loaded_raw.n_obs:
                print(f"  ✓ 主数据和raw数据细胞数匹配")
            else:
                print(f"  ✗ 主数据和raw数据细胞数不匹配")
        else:
            print(f"  ✗ 保存的文件没有raw矩阵")
            
    except Exception as e:
        print(f"  验证失败: {e}")
    
    return merged_adata

def main():
    # 定义要合并的文件列表
    file_list = [
        "celltype_Epithelials.h5ad",
        "celltype_ILIRN+_Epithelials.h5ad", 
        "celltype_Secretory_Epithelials.h5ad",
        "celltype_Stem_Epithelials.h5ad"
    ]
    
    # 检查文件是否存在
    existing_files = []
    for file_path in file_list:
        if os.path.exists(file_path):
            existing_files.append(file_path)
        else:
            print(f"警告: 文件不存在: {file_path}")
    
    if len(existing_files) == 0:
        print("错误: 没有找到任何输入文件")
        return
    
    print(f"找到 {len(existing_files)} 个文件:")
    for file_path in existing_files:
        print(f"  - {file_path}")
    
    # 合并文件
    merged_data = merge_epithelial_h5ad_files(
        file_list=existing_files,
        output_path='merged_all_epithelial_cells.h5ad'
    )
    
    if merged_data is not None:
        print(f"\n合并完成!")
        print(f"总细胞数: {merged_data.n_obs}")
        print(f"总基因数: {merged_data.n_vars}")
        print(f"输出文件: merged_all_epithelial_cells.h5ad")
        
        # 打印详细信息
        print(f"\n详细信息:")
        print(f"观察值列: {merged_data.obs.columns.tolist()}")
        
        if 'annotation' in merged_data.obs.columns:
            print(f"\n细胞类型分布:")
            print(merged_data.obs['annotation'].value_counts())
        
        print(f"\n文件来源分布:")
        print(merged_data.obs['source_file'].value_counts())

if __name__ == "__main__":
    main()
