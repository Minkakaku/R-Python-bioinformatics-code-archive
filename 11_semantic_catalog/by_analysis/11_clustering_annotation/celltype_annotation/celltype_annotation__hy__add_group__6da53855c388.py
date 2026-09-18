import scanpy as sc
import pandas as pd
import numpy as np

def add_sample_annotations_v2(adata_path, output_path=None):
    """
    根据新的样本信息表为AnnData对象添加TIME、PCR、TREATMENT、Patient四列
    
    参数:
    adata_path: 合并后的h5ad文件路径
    output_path: 输出文件路径（可选）
    """
    
    # 读取合并后的数据
    print("读取合并后的数据...")
    adata = sc.read_h5ad(adata_path)
    print(f"原始数据形状: {adata.shape}")
    
    # 检查sample列是否存在
    if 'sample' not in adata.obs.columns:
        raise ValueError("obs中不存在'sample'列，请检查数据")
    
    print(f"样本列中包含的唯一样本数: {adata.obs['sample'].nunique()}")
    print(f"样本示例: {adata.obs['sample'].head().tolist()}")
    
    # 创建新的样本信息映射表（四列）
    sample_info_dict = {
        'TIME': {
            'POST_11-1': 'PRE', 'POST_11-2': 'PRE', 'POST_11-3': 'PRE',
            'POST_3': 'POST', 'POST_5': 'POST', 'POST_6': 'POST', 'POST_8': 'POST',
            'PRE_10-1': 'PRE', 'PRE_10-2': 'PRE', 'PRE_11-1': 'PRE', 'PRE_11-2': 'PRE',
            'PRE_12-1': 'PRE', 'PRE_12-2': 'PRE', 'PRE_13-2': 'PRE', 'PRE_14-1': 'PRE',
            'PRE_14-2': 'PRE', 'PRE_14-3': 'PRE', 'PRE_3-1': 'PRE', 'PRE_3-2': 'PRE',
            'PRE_5-1': 'PRE', 'PRE_5-2': 'PRE', 'PRE_6-1': 'PRE', 'PRE_6-2': 'PRE', 'PRE_8': 'PRE'
        },
        'PCR': {
            'POST_11-1': 'Pcr', 'POST_11-2': 'Pcr', 'POST_11-3': 'Pcr',
            'POST_3': 'non-pCR', 'POST_5': 'Pcr', 'POST_6': 'Pcr', 'POST_8': 'Pcr',
            'PRE_10-1': 'non-pCR', 'PRE_10-2': 'non-pCR', 'PRE_11-1': 'Pcr', 'PRE_11-2': 'Pcr',
            'PRE_12-1': 'non-pCR', 'PRE_12-2': 'non-pCR', 'PRE_13-2': 'Pcr', 'PRE_14-1': 'Pcr',
            'PRE_14-2': 'Pcr', 'PRE_14-3': 'Pcr', 'PRE_3-1': 'non-pCR', 'PRE_3-2': 'non-pCR',
            'PRE_5-1': 'Pcr', 'PRE_5-2': 'Pcr', 'PRE_6-1': 'Pcr', 'PRE_6-2': 'Pcr', 'PRE_8': 'Pcr'
        },
        'TREATMENT': {
            'POST_11-1': 'B', 'POST_11-2': 'B', 'POST_11-3': 'B',
            'POST_3': 'A', 'POST_5': 'A', 'POST_6': 'A', 'POST_8': 'A',
            'PRE_10-1': 'A', 'PRE_10-2': 'A', 'PRE_11-1': 'B', 'PRE_11-2': 'B',
            'PRE_12-1': 'B', 'PRE_12-2': 'B', 'PRE_13-2': 'B', 'PRE_14-1': 'B',
            'PRE_14-2': 'B', 'PRE_14-3': 'B', 'PRE_3-1': 'A', 'PRE_3-2': 'A',
            'PRE_5-1': 'A', 'PRE_5-2': 'A', 'PRE_6-1': 'A', 'PRE_6-2': 'A', 'PRE_8': 'A'
        },
        'Patient': {
            'POST_11-1': 'POST_11', 'POST_11-2': 'POST_11', 'POST_11-3': 'POST_11',
            'POST_3': 'POST_3', 'POST_5': 'POST_5', 'POST_6': 'POST_6', 'POST_8': 'POST_8',
            'PRE_10-1': 'PRE_10', 'PRE_10-2': 'PRE_10', 'PRE_11-1': 'PRE_11', 'PRE_11-2': 'PRE_11',
            'PRE_12-1': 'PRE_12', 'PRE_12-2': 'PRE_12', 'PRE_13-2': 'PRE_13', 'PRE_14-1': 'PRE_14',
            'PRE_14-2': 'PRE_14', 'PRE_14-3': 'PRE_14', 'PRE_3-1': 'PRE_3', 'PRE_3-2': 'PRE_3',
            'PRE_5-1': 'PRE_5', 'PRE_5-2': 'PRE_5', 'PRE_6-1': 'PRE_6', 'PRE_6-2': 'PRE_6', 'PRE_8': 'PRE_8'
        }
    }
    
    # 为每个样本添加新列
    print("\n添加样本注释信息...")
    
    # 使用map函数进行映射[6,7](@ref)
    for col_name, mapping_dict in sample_info_dict.items():
        adata.obs[col_name] = adata.obs['sample'].map(mapping_dict)
        print(f"添加列: {col_name}")
    
    # 检查映射结果
    missing_samples = adata.obs[adata.obs['TIME'].isna()]['sample'].unique()
    if len(missing_samples) > 0:
        print(f"\n警告: 以下样本在映射表中未找到，对应注释将为NaN: {list(missing_samples)}")
    else:
        print("✅ 所有样本都成功映射")
    
    # 统计各列的分布[8](@ref)
    print("\n各列数值分布:")
    for col in ['TIME', 'PCR', 'TREATMENT', 'Patient']:
        if col in adata.obs.columns:
            value_counts = adata.obs[col].value_counts()
            print(f"\n{col}列分布:")
            for value, count in value_counts.items():
                print(f"  {value}: {count} 个细胞")
    
    # 保存结果
    if output_path is None:
        output_path = adata_path.replace('.h5ad', '_with_annotations_v2.h5ad')
    
    print(f"\n保存结果到: {output_path}")
    adata.write_h5ad(output_path)
    
    # 显示最终的数据信息
    print("\n" + "="*60)
    print("最终数据信息:")
    print(f"数据形状: {adata.shape}")
    print(f"观察列(obs)包含: {list(adata.obs.columns)}")
    print(f"样本数: {adata.n_obs}")
    print(f"基因数: {adata.n_vars}")
    
    # 显示各列的独特值数量
    print(f"\n各列独特值数量:")
    for col in ['TIME', 'PCR', 'TREATMENT', 'Patient']:
        if col in adata.obs.columns:
            unique_count = adata.obs[col].nunique()
            print(f"  {col}: {unique_count} 个独特值")
    print("="*60)
    
    return adata

# 使用示例
if __name__ == "__main__":
    # 设置文件路径
    input_file = "merged_data_with_raw.h5ad"  # 您合并后的h5ad文件
    output_file = "merged_data_with_annotations_v2.h5ad"  # 输出文件
    
    try:
        # 执行添加注释的操作
        annotated_adata = add_sample_annotations_v2(
            adata_path=input_file,
            output_path=output_file
        )
        
        print("\n 样本注释添加完成！")
        print(f"输出文件: {output_file}")
        
        # 验证结果
        print("\n验证前10个细胞的注释信息:")
        print(annotated_adata.obs[['sample', 'TIME', 'PCR', 'TREATMENT', 'Patient']].head(10))
        
        # 显示各患者的样本分布
        print("\n各患者的样本分布:")
        patient_sample_counts = annotated_adata.obs.groupby('Patient')['sample'].nunique()
        print(patient_sample_counts)
        
    except Exception as e:
        print(f"\n❌ 处理过程中出现错误: {e}")
