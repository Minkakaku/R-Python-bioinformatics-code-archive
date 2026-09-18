import scanpy as sc
import pandas as pd
import numpy as np
import os
from pathlib import Path

def add_patient_from_table(h5ad_path, mapping_file_path, output_path=None, 
                          sample_col='sample', sep='\t', header=True):
    """
    从映射表中读取sample-patient对应关系，并添加到h5ad文件的obs中
    
    参数:
    h5ad_path: h5ad文件路径
    mapping_file_path: 包含映射关系的文本文件路径
    output_path: 输出h5ad文件路径，如果为None则覆盖原文件
    sample_col: h5ad文件中样本列的名称
    sep: 映射文件的分隔符
    header: 映射文件是否有表头
    """
    
    # 读取h5ad文件
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"原始数据形状: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    
    # 检查sample列是否存在
    if sample_col not in adata.obs.columns:
        raise ValueError(f"在adata.obs中未找到'{sample_col}'列。可用的列有: {list(adata.obs.columns)}")
    
    print(f"数据中包含的样本数量: {adata.obs[sample_col].nunique()}")
    print("样本示例:", list(adata.obs[sample_col].unique())[:10])
    
    # 读取映射文件
    print(f"\n正在读取映射文件: {mapping_file_path}")
    
    # 根据是否有表头决定读取方式
    if header:
        mapping_df = pd.read_csv(mapping_file_path, sep=sep)
        # 假设第一列是sample，第二列是patient
        if len(mapping_df.columns) < 2:
            raise ValueError(f"映射文件至少需要2列，但只找到{len(mapping_df.columns)}列")
        
        sample_col_mapping = mapping_df.columns[0]
        patient_col = mapping_df.columns[1]
    else:
        mapping_df = pd.read_csv(mapping_file_path, sep=sep, header=None)
        if len(mapping_df.columns) < 2:
            raise ValueError(f"映射文件至少需要2列，但只找到{len(mapping_df.columns)}列")
        
        sample_col_mapping = 'sample'
        patient_col = 'patient'
        mapping_df.columns = [sample_col_mapping, patient_col] + [f'col_{i}' for i in range(2, len(mapping_df.columns))]
    
    print(f"映射文件包含 {len(mapping_df)} 行映射关系")
    print(f"映射文件列名: {list(mapping_df.columns)}")
    
    # 创建sample到patient的映射字典
    sample_to_patient = dict(zip(mapping_df[sample_col_mapping], mapping_df[patient_col]))
    
    print(f"\n映射关系示例:")
    for i, (sample, patient) in enumerate(sample_to_patient.items()):
        if i < 5:  # 只显示前5个示例
            print(f"  {sample} -> {patient}")
        else:
            break
    
    # 检查是否有重复的sample
    duplicate_samples = mapping_df[sample_col_mapping].duplicated().sum()
    if duplicate_samples > 0:
        print(f"警告: 映射文件中有 {duplicate_samples} 个重复的sample")
        # 显示重复的sample
        duplicates = mapping_df[mapping_df[sample_col_mapping].duplicated(keep=False)]
        print("重复的sample:")
        for sample in duplicates[sample_col_mapping].unique():
            print(f"  {sample}: {list(duplicates[duplicates[sample_col_mapping] == sample][patient_col])}")
        
        # 使用第一个出现的映射
        sample_to_patient = mapping_df.drop_duplicates(subset=[sample_col_mapping], keep='first')
        sample_to_patient = dict(zip(sample_to_patient[sample_col_mapping], sample_to_patient[patient_col]))
        print("已使用每个sample的第一个映射关系")
    
    # 检查是否已有patient列
    if 'patient' in adata.obs.columns:
        print("警告: 已存在patient列，将覆盖...")
        # 可以选择备份原列
        # adata.obs['patient_old'] = adata.obs['patient']
    
    # 根据sample列添加patient列
    print(f"\n正在添加patient列...")
    adata.obs['patient'] = adata.obs[sample_col].map(sample_to_patient)
    
    # 统计匹配情况
    matched_count = adata.obs['patient'].notna().sum()
    unmatched_count = adata.obs['patient'].isna().sum()
    unmatched_samples = adata.obs[adata.obs['patient'].isna()][sample_col].unique()
    
    print(f"\n=== 匹配结果 ===")
    print(f"成功匹配: {matched_count} 个细胞 ({matched_count/adata.n_obs*100:.2f}%)")
    print(f"未匹配: {unmatched_count} 个细胞 ({unmatched_count/adata.n_obs*100:.2f}%)")
    
    if len(unmatched_samples) > 0:
        print(f"未匹配的样本数量: {len(unmatched_samples)}")
        print("未匹配的样本示例:")
        for sample in unmatched_samples[:10]:  # 只显示前10个
            print(f"  - {sample}")
        if len(unmatched_samples) > 10:
            print(f"  ... 还有 {len(unmatched_samples) - 10} 个未显示")
        
        # 处理未匹配的样本
        choice = input("\n请选择处理未匹配样本的方式:\n1. 将patient设置为sample值\n2. 设置为'Unknown'\n3. 保留为NaN\n请输入选择(1/2/3): ")
        
        if choice == '1':
            adata.obs['patient'] = adata.obs['patient'].fillna(adata.obs[sample_col])
            print("已将未匹配样本的patient设置为sample值")
        elif choice == '2':
            adata.obs['patient'] = adata.obs['patient'].fillna('Unknown')
            print("已将未匹配样本的patient设置为'Unknown'")
        else:
            print("未匹配样本的patient将保留为NaN")
    
    # 显示最终的映射关系
    print(f"\n=== 最终的样本-patient映射 ===")
    sample_patient_final = adata.obs[[sample_col, 'patient']].drop_duplicates().sort_values(sample_col)
    
    for idx, row in sample_patient_final.iterrows():
        print(f"  {row[sample_col]} -> {row['patient']}")
    
    # 保存结果
    if output_path is None:
        # 在原文件名后添加 _with_patient
        base_name = os.path.splitext(h5ad_path)[0]
        output_path = f"{base_name}_with_patient.h5ad"
    
    print(f"\n正在保存结果到: {output_path}")
    adata.write(output_path)
    
    # 显示统计信息
    print(f"\n=== 最终统计 ===")
    print(f"总细胞数: {adata.n_obs}")
    print(f"总基因数: {adata.n_vars}")
    print(f"唯一patient数: {adata.obs['patient'].nunique()}")
    print(f"唯一sample数: {adata.obs[sample_col].nunique()}")
    
    # 显示patient分布
    if adata.obs['patient'].notna().any():
        print(f"\npatient分布 (前10个):")
        patient_counts = adata.obs['patient'].value_counts().head(10)
        for patient, count in patient_counts.items():
            percentage = count / adata.n_obs * 100
            print(f"  {patient}: {count} 个细胞 ({percentage:.2f}%)")
    
    return adata

# 使用示例
if __name__ == "__main__":
    # 文件路径
    h5ad_path = "adata_with_CNs10_TMA_5.h5ad"  # 请替换为您的h5ad文件路径
    mapping_file_path = "releate.txt"  # 映射文件路径
    
    # 输出文件路径（可选）
    output_path = None  # 如果为None，则自动生成
    
    try:
        # 检查文件是否存在
        if not os.path.exists(h5ad_path):
            raise FileNotFoundError(f"h5ad文件不存在: {h5ad_path}")
        
        if not os.path.exists(mapping_file_path):
            raise FileNotFoundError(f"映射文件不存在: {mapping_file_path}")
        
        # 执行添加patient列的操作
        adata = add_patient_from_table(
            h5ad_path=h5ad_path,
            mapping_file_path=mapping_file_path,
            output_path=output_path,
            sample_col='sample',  # 如果您的样本列名不同，请修改此项
            sep='\t',             # 如果您的文件使用其他分隔符，请修改此项
            header=True           # 如果您的映射文件没有表头，请设置为False
        )
        
        print(f"\n处理完成！")
        
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()
