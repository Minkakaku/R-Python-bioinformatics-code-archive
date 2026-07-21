import scanpy as sc
import pandas as pd
import numpy as np
import os

# 1. 读取h5ad文件
input_file = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/state/All_filter_with_state_new.h5ad"
print(f"读取h5ad文件: {input_file}")
adata = sc.read_h5ad(input_file)

# 2. 打印adata的基本信息
print("\nadata基本信息:")
print(adata)
print(f"\n数据形状: {adata.shape}")
print(f"观测变量(obs)列: {list(adata.obs.columns)}")

# 3. 检查是否有'sample'列
if 'sample' in adata.obs.columns:
    print(f"\n样本列信息:")
    print(f"样本列数据类型: {adata.obs['sample'].dtype}")
    print(f"样本列唯一值数量: {adata.obs['sample'].nunique()}")
    
    # 获取所有唯一的样本名称
    unique_samples = adata.obs['sample'].unique()
    print(f"前10个样本名称: {unique_samples[:10]}")
    
    # 检查特定样本是否存在
    target_sample = 'SR202500016K'
    
    # 方法1: 使用in运算符检查
    if target_sample in adata.obs['sample'].values:
        print(f"\n样本 '{target_sample}' 存在!")
        
        # 方法2: 获取该样本的细胞数量
        sample_mask = adata.obs['sample'] == target_sample
        sample_cells_count = sample_mask.sum()
        print(f"样本 '{target_sample}' 有 {sample_cells_count} 个细胞")
        
        # 提取该样本的数据
        sample_data = adata[adata.obs['sample'] == target_sample, :]
        print(f"提取的数据形状: {sample_data.shape}")
        
        # 如果有state列，查看分布
        if 'state' in adata.obs.columns:
            state_counts = sample_data.obs['state'].value_counts(dropna=False)
            print(f"\n样本 '{target_sample}' 的state分布:")
            for state, count in state_counts.items():
                percentage = count / sample_data.n_obs * 100
                state_display = 'NaN' if pd.isna(state) else state
                print(f"  {state_display}: {count} ({percentage:.2f}%)")
    else:
        print(f"\n样本 '{target_sample}' 不存在!")
        print(f"最相似的样本名称:")
        
        # 查找相似的样本名称
        all_samples = adata.obs['sample'].unique()
        similar_samples = [s for s in all_samples if str(s).startswith('SR202500016')]
        if similar_samples:
            print(f"  以'SK202500016'开头的样本: {similar_samples}")
        else:
            print(f"  没有找到相似的样本")
    
    # 显示所有样本的统计信息
    print(f"\n所有样本的细胞数量:")
    sample_counts = adata.obs['sample'].value_counts()
    for sample, count in sample_counts.head(20).items():  # 只显示前20个
        print(f"  {sample}: {count} 个细胞")
    
    if len(sample_counts) > 20:
        print(f"  ... 还有 {len(sample_counts) - 20} 个样本未显示")
    
else:
    print("错误: adata.obs中没有'sample'列")

# 4. 如果有state列，查看整体分布
if 'state' in adata.obs.columns:
    print(f"\nstate列信息:")
    print(f"state列数据类型: {adata.obs['state'].dtype}")
    print(f"state列唯一值: {adata.obs['state'].unique()}")
    
    state_counts = adata.obs['state'].value_counts(dropna=False)
    print(f"state分布:")
    for state, count in state_counts.items():
        percentage = count / adata.n_obs * 100
        state_display = 'NaN' if pd.isna(state) else state
        print(f"  {state_display}: {count} ({percentage:.2f}%)")

# 5. 创建完整的样本列表
if 'sample' in adata.obs.columns:
    output_dir = "./sample_analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    # 保存所有样本列表
    all_samples = adata.obs['sample'].unique()
    samples_file = os.path.join(output_dir, "all_samples.txt")
    
    with open(samples_file, 'w') as f:
        f.write("所有样本列表:\n")
        f.write("=" * 40 + "\n")
        for sample in sorted(all_samples):
            f.write(f"{sample}\n")
    
    print(f"\n所有样本列表已保存到: {samples_file}")
    
    # 保存详细的样本统计
    stats_file = os.path.join(output_dir, "sample_statistics.csv")
    sample_stats = adata.obs.groupby('sample').size().reset_index()
    sample_stats.columns = ['sample', 'cell_count']
    sample_stats = sample_stats.sort_values('cell_count', ascending=False)
    sample_stats.to_csv(stats_file, index=False)
    
    print(f"样本统计已保存到: {stats_file}")
    
    # 如果有state列，保存每个样本的state分布
    if 'state' in adata.obs.columns:
        state_by_sample = adata.obs.groupby(['sample', 'state']).size().unstack(fill_value=0)
        state_by_sample_file = os.path.join(output_dir, "state_by_sample.csv")
        state_by_sample.to_csv(state_by_sample_file)
        print(f"各样本state分布已保存到: {state_by_sample_file}")

print("\n分析完成!")
