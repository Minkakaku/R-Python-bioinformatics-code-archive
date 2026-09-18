import scanpy as sc
import pandas as pd
import numpy as np
import random
import json
import os
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

def stratified_sample_h5ad(input_h5ad, output_h5ad, subset_size=20000, random_seed=42):
    """
    从h5ad数据中按照sample和celltype的比例抽取指定数量的细胞子集
    
    参数:
    - input_h5ad: 输入h5ad文件路径
    - output_h5ad: 输出h5ad文件路径
    - subset_size: 要抽取的细胞数
    - random_seed: 随机种子
    """
    
    # 设置随机种子保证结果可重复
    random.seed(random_seed)
    np.random.seed(random_seed)
    
    print(f"正在读取数据: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    
    # 检查必需列
    required_cols = ['sample', 'celltype']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"缺少必需列: {missing_cols}")
    
    total_cells = adata.shape[0]
    print(f"原始数据: {total_cells} 个细胞")
    print(f"样本数: {adata.obs['sample'].nunique()}")
    print(f"细胞类型数: {adata.obs['celltype'].nunique()}")
    
    # 创建DataFrame方便处理
    df = pd.DataFrame({
        'sample': adata.obs['sample'],
        'celltype': adata.obs['celltype']
    }, index=adata.obs.index)
    
    # 分析原始分布
    print("\n原始数据分布分析:")
    print("="*50)
    
    # sample分布
    sample_counts = df['sample'].value_counts()
    sample_proportions = sample_counts / total_cells
    print(f"\nsample分布 (前10个):")
    for i, (sample, count) in enumerate(sample_counts.head(10).items()):
        prop = sample_proportions[sample] * 100
        print(f"  {sample}: {count} cells ({prop:.2f}%)")
    
    # celltype分布
    celltype_counts = df['celltype'].value_counts()
    celltype_proportions = celltype_counts / total_cells
    print(f"\ncelltype分布 (前10个):")
    for i, (celltype, count) in enumerate(celltype_counts.head(10).items()):
        prop = celltype_proportions[celltype] * 100
        print(f"  {celltype}: {count} cells ({prop:.2f}%)")
    
    # 组合分布
    combined_counts = df.groupby(['sample', 'celltype']).size()
    combined_proportions = combined_counts / total_cells
    print(f"\n前20个sample-celltype组合:")
    for i, ((sample, celltype), count) in enumerate(combined_counts.sort_values(ascending=False).head(20).items()):
        prop = combined_proportions[(sample, celltype)] * 100
        print(f"  {i+1:2d}. {sample} - {celltype}: {count} cells ({prop:.3f}%)")
    
    # 计算每个组合应该抽取的数量
    print(f"\n计算每个sample-celltype组合的抽样数量...")
    target_counts = {}
    
    # 计算每个组合的目标数量（向下取整）
    for (sample, celltype), count in combined_counts.items():
        proportion = count / total_cells
        target = int(np.floor(proportion * subset_size))
        if target > 0:
            target_counts[(sample, celltype)] = min(target, count)
    
    # 计算当前分配的总数
    allocated = sum(target_counts.values())
    remaining = subset_size - allocated
    print(f"已分配: {allocated} 个细胞, 剩余: {remaining} 个细胞需要分配")
    
    # 将剩余细胞分配给比例较大的组合
    if remaining > 0:
        # 按比例排序
        sorted_combinations = sorted(combined_counts.items(), key=lambda x: x[1], reverse=True)
        
        for (sample, celltype), count in sorted_combinations:
            if remaining <= 0:
                break
                
            # 这个组合已分配的数量
            allocated_in_group = target_counts.get((sample, celltype), 0)
            
            # 这个组合还能分配多少
            available = count - allocated_in_group
            if available > 0:
                # 分配1个细胞
                target_counts[(sample, celltype)] = allocated_in_group + 1
                remaining -= 1
    
    # 如果还有剩余，分配给随机组合
    if remaining > 0:
        all_combinations = list(combined_counts.items())
        for _ in range(remaining):
            (sample, celltype), count = random.choice(all_combinations)
            target_counts[(sample, celltype)] = target_counts.get((sample, celltype), 0) + 1
    
    # 验证总数
    final_total = sum(target_counts.values())
    if final_total != subset_size:
        # 调整最大的组合
        max_key = max(target_counts, key=target_counts.get)
        target_counts[max_key] += (subset_size - final_total)
    
    print(f"最终分配总数: {sum(target_counts.values())} 个细胞")
    
    # 执行抽样
    print(f"\n开始抽样...")
    selected_indices = []
    
    for (sample, celltype), target_count in target_counts.items():
        # 获取这个组合的所有细胞索引
        group_indices = df[(df['sample'] == sample) & (df['celltype'] == celltype)].index.tolist()
        
        if len(group_indices) >= target_count:
            # 随机抽取
            selected = random.sample(group_indices, target_count)
        else:
            # 如果不够，全部选取
            selected = group_indices
            print(f"  警告: {sample}-{celltype} 只有 {len(group_indices)} 个细胞，少于目标 {target_count}")
        
        selected_indices.extend(selected)
    
    # 确保数量正确
    if len(selected_indices) > subset_size:
        selected_indices = random.sample(selected_indices, subset_size)
    elif len(selected_indices) < subset_size:
        # 补充随机细胞
        remaining_needed = subset_size - len(selected_indices)
        all_indices = df.index.tolist()
        remaining_indices = [idx for idx in all_indices if idx not in selected_indices]
        additional = random.sample(remaining_indices, min(remaining_needed, len(remaining_indices)))
        selected_indices.extend(additional)
    
    print(f"实际抽取细胞数: {len(selected_indices)}")
    
    # 创建子集
    adata_subset = adata[selected_indices, :].copy()
    
    # 保存子集
    adata_subset.write_h5ad(output_h5ad)
    print(f"\n子集已保存到: {output_h5ad}")
    print(f"子集形状: {adata_subset.shape}")
    
    # 分析子集分布
    print(f"\n子集数据分布分析:")
    print("="*50)
    
    # sample分布
    subset_sample_counts = adata_subset.obs['sample'].value_counts()
    subset_sample_proportions = subset_sample_counts / len(adata_subset)
    
    print(f"\n子集sample分布 (前10个):")
    for i, (sample, count) in enumerate(subset_sample_counts.head(10).items()):
        orig_prop = sample_proportions.get(sample, 0) * 100
        subset_prop = subset_sample_proportions[sample] * 100
        diff = subset_prop - orig_prop
        print(f"  {sample}: {count} cells ({subset_prop:.2f}%, 原始: {orig_prop:.2f}%, 差异: {diff:+.2f}%)")
    
    # celltype分布
    subset_celltype_counts = adata_subset.obs['celltype'].value_counts()
    subset_celltype_proportions = subset_celltype_counts / len(adata_subset)
    
    print(f"\n子集celltype分布 (前10个):")
    for i, (celltype, count) in enumerate(subset_celltype_counts.head(10).items()):
        orig_prop = celltype_proportions.get(celltype, 0) * 100
        subset_prop = subset_celltype_proportions[celltype] * 100
        diff = subset_prop - orig_prop
        print(f"  {celltype}: {count} cells ({subset_prop:.2f}%, 原始: {orig_prop:.2f}%, 差异: {diff:+.2f}%)")
    
    # 保存抽样信息
    sampling_info = {
        'input_file': input_h5ad,
        'output_file': output_h5ad,
        'total_cells_original': int(total_cells),
        'subset_size_requested': subset_size,
        'subset_size_actual': len(adata_subset),
        'random_seed': random_seed,
        'original_sample_distribution': {k: float(v) for k, v in sample_proportions.items()},
        'subset_sample_distribution': {k: float(v) for k, v in subset_sample_proportions.items()},
        'original_celltype_distribution': {k: float(v) for k, v in celltype_proportions.items()},
        'subset_celltype_distribution': {k: float(v) for k, v in subset_celltype_proportions.items()},
        'target_counts_per_group': {f"{k[0]}_{k[1]}": int(v) for k, v in target_counts.items()},
        'selected_indices': selected_indices
    }
    
    info_file = output_h5ad.replace('.h5ad', '_sampling_info.json')
    with open(info_file, 'w') as f:
        json.dump(sampling_info, f, indent=2)
    
    print(f"\n抽样信息已保存到: {info_file}")
    
    return adata_subset, sampling_info

# 主函数
if __name__ == "__main__":
    # 配置参数
    input_file = "data_sp_ep.h5ad"  # 你的输入文件
    output_file = "data_sp_ep_subset_20000.h5ad"  # 输出文件
    subset_size = 20000
    random_seed = 42
    
    try:
        # 运行抽样
        subset_data, sampling_info = stratified_sample_h5ad(
            input_h5ad=input_file,
            output_h5ad=output_file,
            subset_size=subset_size,
            random_seed=random_seed
        )
        
        print(f"\n抽样完成!")
        print(f"原始数据: {sampling_info['total_cells_original']} 个细胞")
        print(f"抽样子集: {sampling_info['subset_size_actual']} 个细胞")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
