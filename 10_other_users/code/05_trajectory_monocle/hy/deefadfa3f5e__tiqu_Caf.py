import scanpy as sc
import pandas as pd
import numpy as np

# 1. 读取h5ad文件
adata = sc.read_h5ad('adata_with_band_by_sample.h5ad')

print(f"原始数据形状: {adata.shape}")
print(f"细胞类型: {adata.obs['celltype'].unique()}")

# 检查是否有raw
if adata.raw is not None:
    print(f"原始数据有raw矩阵，形状: {adata.raw.shape}")
else:
    print("警告: 原始数据没有raw矩阵")

# 2. 提取指定的细胞类型
target_celltypes = [
    'myCAF', 
    'iCAF', 
    'Matrix Fibroblasts', 
    'Matrix remodeling CAFs', 
    'cycling CAF', 
    'EMT-like CAFs'
]

# 过滤出目标细胞类型
adata_filtered = adata[adata.obs['celltype'].isin(target_celltypes)].copy()

print(f"\n过滤后数据形状: {adata_filtered.shape}")
print(f"各细胞类型数量:")
print(adata_filtered.obs['celltype'].value_counts())

# 3. 按样本+细胞类型进行分层采样
total_cells_to_sample = 30000

# 获取所有样本
samples = adata_filtered.obs['sample'].unique()
print(f"\n共有 {len(samples)} 个样本")

# 存储采样后的索引
sampled_indices = []

# 对每个样本进行处理
for sample in samples:
    # 获取该样本的所有细胞
    sample_mask = adata_filtered.obs['sample'] == sample
    sample_adata = adata_filtered[sample_mask]
    
    # 计算该样本中每种细胞类型的比例
    celltype_counts = sample_adata.obs['celltype'].value_counts()
    celltype_props = celltype_counts / celltype_counts.sum()
    
    # 计算该样本应该抽取的总细胞数（按样本总细胞数比例）
    sample_total_cells = len(sample_adata)
    
    # 按样本总细胞数比例分配采样数量
    sample_size = int(round((sample_total_cells / adata_filtered.shape[0]) * total_cells_to_sample))
    
    # 如果某个样本计算得到的采样数为0，至少分配1个（如果该样本有细胞的话）
    if sample_size == 0 and sample_total_cells > 0:
        sample_size = 1
    
    print(f"\n样本 {sample}: 总细胞数={sample_total_cells}, 计划采样={sample_size}")
    
    # 如果该样本细胞数少于计划采样数，全部取用
    if sample_total_cells <= sample_size:
        sampled_indices.extend(sample_adata.obs.index.tolist())
        print(f"  -> 样本细胞数不足，全部取用 ({sample_total_cells} 个)")
        continue
    
    # 按细胞类型比例在该样本内进行采样
    sampled_indices_sample = []
    for celltype, prop in celltype_props.items():
        # 计算该细胞类型应该抽取的数量
        n_cells_for_type = int(round(prop * sample_size))
        
        # 获取该样本中该细胞类型的所有细胞索引
        type_mask = (adata_filtered.obs['sample'] == sample) & (adata_filtered.obs['celltype'] == celltype)
        type_indices = adata_filtered[type_mask].obs.index.tolist()
        
        # 如果该类型细胞数少于计划数，全部取用
        if len(type_indices) <= n_cells_for_type:
            sampled_indices_sample.extend(type_indices)
            print(f"  {celltype}: 计划{n_cells_for_type}个, 实际{len(type_indices)}个 (不足)")
        else:
            # 随机抽取
            sampled_indices_sample.extend(np.random.choice(type_indices, size=n_cells_for_type, replace=False))
            print(f"  {celltype}: 抽取{n_cells_for_type}个")
    
    # 修正由于四舍五入导致的采样数量偏差
    # 如果采样数量少于计划，从该样本剩余细胞中随机补充
    if len(sampled_indices_sample) < sample_size:
        # 获取该样本中还未被采样的细胞
        sampled_set = set(sampled_indices_sample)
        all_sample_indices = set(sample_adata.obs.index.tolist())
        remaining_indices = list(all_sample_indices - sampled_set)
        
        if remaining_indices:
            n_need = sample_size - len(sampled_indices_sample)
            n_need = min(n_need, len(remaining_indices))
            sampled_indices_sample.extend(np.random.choice(remaining_indices, size=n_need, replace=False))
            print(f"  -> 补充抽取 {n_need} 个细胞")
    
    # 如果采样数量多于计划，随机删除多余的
    elif len(sampled_indices_sample) > sample_size:
        n_remove = len(sampled_indices_sample) - sample_size
        sampled_indices_sample = np.random.choice(sampled_indices_sample, size=sample_size, replace=False).tolist()
        print(f"  -> 随机删除 {n_remove} 个细胞")
    
    sampled_indices.extend(sampled_indices_sample)
    print(f"  该样本最终采样: {len(sampled_indices_sample)} 个细胞")

# 创建采样后的数据
adata_sampled = adata_filtered[sampled_indices].copy()

print(f"\n{'='*50}")
print(f"采样完成！")
print(f"采样后数据形状: {adata_sampled.shape}")
print(f"\n各样本采样后细胞数量:")
print(adata_sampled.obs['sample'].value_counts().sort_index())

print(f"\n各细胞类型采样后数量:")
print(adata_sampled.obs['celltype'].value_counts())

# 4. 保存结果，同时保留raw矩阵
print(f"\n{'='*50}")
print("保存数据...")

# 方法1: 直接保存（会保留raw矩阵）
output_file = 'adata_sampled_30000.h5ad'
adata_sampled.write(output_file)
print(f"数据已保存到: {output_file}")

# 验证保存的文件是否包含raw矩阵
print("\n验证保存的文件...")
adata_loaded = sc.read_h5ad(output_file)
if adata_loaded.raw is not None:
    print(f"✓ 成功! raw矩阵已保留，形状: {adata_loaded.raw.shape}")
else:
    print("✗ 警告: raw矩阵未保留")

# 5. 额外保存一份带有详细统计信息的报告
print(f"\n生成采样报告...")
report_file = 'sampling_report.txt'
with open(report_file, 'w') as f:
    f.write("="*60 + "\n")
    f.write("采样报告\n")
    f.write("="*60 + "\n\n")
    
    f.write(f"原始数据总细胞数: {adata.shape[0]}\n")
    f.write(f"过滤后细胞数: {adata_filtered.shape[0]}\n")
    f.write(f"目标细胞类型: {', '.join(target_celltypes)}\n")
    f.write(f"采样后细胞数: {adata_sampled.shape[0]}\n\n")
    
    f.write("各细胞类型数量变化:\n")
    f.write("-"*40 + "\n")
    original_counts = adata_filtered.obs['celltype'].value_counts()
    sampled_counts = adata_sampled.obs['celltype'].value_counts()
    for ct in target_celltypes:
        orig = original_counts.get(ct, 0)
        samp = sampled_counts.get(ct, 0)
        f.write(f"{ct:25s}: 原始 {orig:6d} -> 采样 {samp:6d} ({(samp/orig*100):.1f}%)\n")
    
    f.write("\n各样本采样情况:\n")
    f.write("-"*40 + "\n")
    sample_counts_orig = adata_filtered.obs['sample'].value_counts()
    sample_counts_samp = adata_sampled.obs['sample'].value_counts()
    for sample in sorted(sample_counts_orig.index):
        orig = sample_counts_orig.get(sample, 0)
        samp = sample_counts_samp.get(sample, 0)
        f.write(f"{sample:20s}: 原始 {orig:6d} -> 采样 {samp:6d} ({(samp/orig*100):.1f}%)\n")
    
    f.write("\n" + "="*60 + "\n")
    
print(f"采样报告已保存到: {report_file}")

print("\n所有操作完成！")
