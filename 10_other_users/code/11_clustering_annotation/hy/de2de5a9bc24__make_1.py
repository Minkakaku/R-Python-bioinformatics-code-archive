import scanpy as sc
import pandas as pd
import os
from pathlib import Path

def extract_samples_and_save(h5ad_path, target_samples, output_path):
    """
    从h5ad文件中提取指定样本并保存
    
    参数:
    h5ad_path: h5ad文件路径
    target_samples: 要提取的样本列表
    output_path: 输出文件路径
    """
    # 读取h5ad文件
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"原始数据: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    print(f"原始样本数量: {adata.obs['sample'].nunique()}")
    
    # 检查sample列是否存在
    if 'sample' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'sample'列。可用的列有: " + ", ".join(adata.obs.columns))
    
    # 检查哪些样本存在于数据中
    existing_samples = [sample for sample in target_samples if sample in adata.obs['sample'].values]
    missing_samples = [sample for sample in target_samples if sample not in adata.obs['sample'].values]
    
    print(f"目标样本数: {len(target_samples)}")
    print(f"找到的样本数: {len(existing_samples)}")
    print(f"缺失的样本数: {len(missing_samples)}")
    
    if missing_samples:
        print(f"未找到的样本: {missing_samples}")
    
    if len(existing_samples) == 0:
        print("错误: 没有找到任何目标样本！")
        print("数据集中存在的样本:", list(adata.obs['sample'].unique())[:20])
        return None
    
    # 提取目标样本
    adata_subset = adata[adata.obs['sample'].isin(existing_samples), :].copy()
    
    print(f"提取后的细胞数: {adata_subset.n_obs}")
    print(f"提取后的基因数: {adata_subset.n_vars}")
    print("\n样本分布统计:")
    sample_stats = adata_subset.obs['sample'].value_counts()
    print(sample_stats)
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 保存为新的h5ad文件
    adata_subset.write(output_path)
    print(f"数据已保存为: {output_path}")
    print("-" * 50)
    
    return adata_subset

def process_sample_files(sample_list_file, h5ad_path):
    """
    处理sample文件列表，为每个样本提取对应的数据
    
    参数:
    sample_list_file: 包含样本标识符的文件路径
    h5ad_path: 原始h5ad文件路径
    """
    # 读取sample文件
    with open(sample_list_file, 'r') as f:
        sample_ids = [line.strip() for line in f if line.strip()]
    
    print(f"从 {sample_list_file} 中读取到 {len(sample_ids)} 个样本标识符")
    print("=" * 80)
    
    # 读取原始h5ad文件一次
    print("正在读取原始h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"原始数据: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    print(f"原始样本数量: {adata.obs['sample'].nunique()}")
    print("=" * 80)
    
    # 处理每个样本标识符
    for sample_id in sample_ids:
        print(f"\n处理样本: {sample_id}")
        
        # 检查对应的.sample文件是否存在
        sample_file = f"{sample_id}.sample"
        if not os.path.exists(sample_file):
            print(f"警告: 找不到样本文件 {sample_file}，跳过此样本")
            continue
        
        # 读取.sample文件中的目标样本列表
        with open(sample_file, 'r') as f:
            target_samples = [line.strip() for line in f if line.strip()]
        
        print(f"从 {sample_file} 中读取到 {len(target_samples)} 个目标样本: {target_samples}")
        
        # 检查哪些目标样本存在于数据中
        existing_samples = [sample for sample in target_samples if sample in adata.obs['sample'].values]
        missing_samples = [sample for sample in target_samples if sample not in adata.obs['sample'].values]
        
        if len(existing_samples) == 0:
            print(f"错误: 在数据中未找到 {sample_id} 的任何目标样本，跳过")
            if missing_samples:
                print(f"未找到的样本: {missing_samples}")
            continue
        
        if missing_samples:
            print(f"警告: 部分样本在数据中不存在: {missing_samples}")
        
        # 提取目标样本
        adata_subset = adata[adata.obs['sample'].isin(existing_samples), :].copy()
        
        print(f"提取后的细胞数: {adata_subset.n_obs}")
        print(f"提取后的样本分布:")
        sample_stats = adata_subset.obs['sample'].value_counts()
        for sample_name, count in sample_stats.items():
            print(f"  {sample_name}: {count} 个细胞")
        
        # 创建输出目录和文件路径
        output_dir = sample_id
        output_path = os.path.join(output_dir, f"{sample_id}.h5ad")
        
        # 确保输出目录存在
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # 保存为新的h5ad文件
        adata_subset.write(output_path)
        print(f"数据已保存为: {output_path}")
        print("-" * 50)

# 执行提取
if __name__ == "__main__":
    # 文件路径配置
    h5ad_path = "merged_single_cell_data_with_group_new.h5ad"  # 统一的h5ad文件
    sample_list_file = "sample"  # 包含样本标识符的文件
    
    # 检查文件是否存在
    if not os.path.exists(h5ad_path):
        print(f"错误: 找不到h5ad文件: {h5ad_path}")
    elif not os.path.exists(sample_list_file):
        print(f"错误: 找不到sample文件: {sample_list_file}")
    else:
        # 处理所有样本
        process_sample_files(sample_list_file, h5ad_path)
        print("\n处理完成!")
