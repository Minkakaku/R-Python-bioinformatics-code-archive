#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scanpy as sc
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

# 定义要合并的文件列表
file_list = [
    "adata_with_CNs10_TMA-1_with_patient_MZ.h5ad",
    "adata_with_CNs10_TMA-3_with_patient_MZ.h5ad", 
    "adata_with_CNs10_TMA-4_with_patient_MZ.h5ad",
    "adata_with_CNs10_TMA-5_with_patient_MZ.h5ad",
    "adata_with_CNs10_TMA_with_patient_MZ.h5ad"
]

# 输出文件
output_file = "merged_TMA_data.h5ad"

print("="*60)
print("合并5个TMA h5ad文件")
print("="*60)

# 检查文件是否存在
existing_files = []
for f in file_list:
    if os.path.exists(f):
        existing_files.append(f)
        print(f"找到文件: {f}")
    else:
        print(f"警告: 未找到文件: {f}")

if len(existing_files) == 0:
    print("错误: 没有找到任何文件")
    exit(1)

print(f"\n将合并 {len(existing_files)} 个文件")

# 读取所有文件
adata_list = []
for i, file_path in enumerate(existing_files):
    print(f"\n读取文件 {i+1}: {file_path}")
    try:
        adata = sc.read_h5ad(file_path)
        print(f"  细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
        
        # 检查是否有raw矩阵
        if adata.raw is not None:
            print(f"  有raw矩阵: 是, 形状: {adata.raw.shape}")
        else:
            print(f"  有raw矩阵: 否")
        
        # 检查obs索引是否有重复
        if adata.obs.index.duplicated().any():
            print(f"  警告: 当前文件的obs索引有重复值")
            print(f"  重复的cell_id: {adata.obs.index[adata.obs.index.duplicated()].unique()[:5]}")
            print(f"  将在合并时自动处理重复索引")
        
        # 添加批次信息
        batch_name = os.path.splitext(os.path.basename(file_path))[0]
        adata.obs['batch'] = batch_name
        
        adata_list.append(adata)
    except Exception as e:
        print(f"  读取失败: {e}")

if len(adata_list) == 0:
    print("错误: 没有成功读取任何文件")
    exit(1)

# 合并文件
print(f"\n开始合并 {len(adata_list)} 个文件...")

# 检查基因是否一致
first_var = adata_list[0].var.index.tolist()
all_same = True
for i, adata in enumerate(adata_list[1:], 1):
    if not adata.var.index.equals(pd.Index(first_var)):
        all_same = False
        print(f"警告: 文件 {i+1} 的基因与第一个文件不一致")
        break

if all_same:
    print("所有文件的基因一致，使用内连接合并")
    join_method = 'inner'
else:
    print("文件的基因不一致，使用外连接合并（缺失值用0填充）")
    join_method = 'outer'

# 执行合并
try:
    merged_adata = sc.concat(
        adata_list,
        join=join_method,
        label='batch',
        keys=[os.path.splitext(os.path.basename(f))[0] for f in existing_files],
        index_unique='-'  # 只在有重复时才添加后缀
    )
    
    print(f"\n主数据合并完成!")
    print(f"合并后数据:")
    print(f"  细胞数: {merged_adata.n_obs}")
    print(f"  基因数: {merged_adata.n_vars}")
    print(f"  obs形状: {merged_adata.obs.shape}")
    print(f"  var形状: {merged_adata.var.shape}")
    
    # 检查合并后的索引是否有重复
    if merged_adata.obs.index.duplicated().any():
        print(f"  警告: 合并后仍有重复索引")
        duplicate_count = merged_adata.obs.index.duplicated().sum()
        print(f"  重复索引数量: {duplicate_count}")
    else:
        print(f"  合并后索引唯一: 是")
    
    # 处理raw矩阵
    print(f"\n处理raw矩阵...")
    
    # 检查所有文件是否都有raw矩阵
    all_have_raw = all([adata.raw is not None for adata in adata_list])
    
    if all_have_raw:
        print("所有文件都有raw矩阵，将合并raw矩阵")
        
        # 分别提取每个文件的raw矩阵
        raw_adata_list = []
        for i, adata in enumerate(adata_list):
            # 从raw对象创建AnnData
            raw_adata = sc.AnnData(
                X=adata.raw.X,
                var=adata.raw.var,
                obs=adata.obs[[]].copy()  # 创建一个空的obs DataFrame，但保持正确的索引
            )
            
            # 添加批次信息
            batch_name = adata.obs['batch'].iloc[0] if 'batch' in adata.obs.columns else f"batch_{i+1}"
            raw_adata.obs['batch'] = batch_name
            
            raw_adata_list.append(raw_adata)
        
        # 合并raw矩阵
        merged_raw = sc.concat(
            raw_adata_list,
            join=join_method,
            label='batch',
            keys=[os.path.splitext(os.path.basename(f))[0] for f in existing_files],
            index_unique='-'
        )
        
        # 设置回merged_adata
        merged_adata.raw = merged_raw
        print(f"  raw矩阵合并完成:")
        print(f"    raw矩阵形状: {merged_adata.raw.shape}")
    else:
        print("不是所有文件都有raw矩阵，跳过raw矩阵合并")
        print(f"  有raw矩阵的文件数: {sum([1 for adata in adata_list if adata.raw is not None])}/{len(adata_list)}")
    
    # 检查批次分布
    if 'batch' in merged_adata.obs.columns:
        print(f"\n批次分布:")
        batch_counts = merged_adata.obs['batch'].value_counts()
        for batch, count in batch_counts.items():
            print(f"  {batch}: {count} 细胞 ({count/merged_adata.n_obs*100:.1f}%)")
    
    # 保存合并后的文件
    print(f"\n保存到: {output_file}")
    merged_adata.write_h5ad(output_file)
    print("保存成功!")
    
    # 验证
    print(f"\n验证保存的文件:")
    verify_adata = sc.read_h5ad(output_file)
    print(f"  验证通过: 成功读取 {verify_adata.n_obs} 个细胞")
    
    if verify_adata.raw is not None:
        print(f"  raw矩阵存在: 是, 形状: {verify_adata.raw.shape}")
    else:
        print(f"  raw矩阵存在: 否")
    
    print("\n" + "="*60)
    print("合并完成!")
    print("="*60)
    
except Exception as e:
    print(f"合并过程中出错: {e}")
    import traceback
    traceback.print_exc()
