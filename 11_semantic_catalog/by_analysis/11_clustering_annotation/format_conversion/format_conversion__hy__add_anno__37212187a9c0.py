import scanpy as sc
import pandas as pd
import numpy as np
import sys
from pathlib import Path

def main():
    # 文件路径
    h5ad_file = "adata_filtered.h5ad"
    meta_file = "meta.data"
    output_file = "all_raw.h5ad"
    
    # 检查文件是否存在
    for file in [h5ad_file, meta_file]:
        if not Path(file).exists():
            print(f"错误: 文件 {file} 不存在！")
            sys.exit(1)
    
    print("1. 读取h5ad数据...")
    adata = sc.read_h5ad(h5ad_file)
    print(f"   h5ad数据形状: {adata.shape}")
    print(f"   h5ad obs列: {adata.obs.columns.tolist()}")
    
    # 检查h5ad中是否有必要的列
    required_h5ad_cols = ['cell_id', 'sample']
    for col in required_h5ad_cols:
        if col not in adata.obs.columns:
            print(f"错误: h5ad数据中缺少列 '{col}'")
            print(f"可用列: {adata.obs.columns.tolist()}")
            sys.exit(1)
    
    print("\n2. 读取meta.data文件...")
    # 尝试不同格式读取meta数据
    try:
        # 先尝试制表符分隔
        meta_df = pd.read_csv(meta_file, sep='\t')
    except:
        try:
            # 尝试逗号分隔
            meta_df = pd.read_csv(meta_file, sep=',')
        except:
            try:
                # 尝试空格分隔
                meta_df = pd.read_csv(meta_file, sep='\s+')
            except:
                print("错误: 无法读取meta.data文件，请确保是TSV、CSV或空格分隔格式")
                sys.exit(1)
    
    print(f"   meta.data数据形状: {meta_df.shape}")
    print(f"   meta.data列: {meta_df.columns.tolist()}")
    
    # 检查meta中是否有必要的列
    required_meta_cols = ['cell_id', 'sample', 'celltype']
    for col in required_meta_cols:
        if col not in meta_df.columns:
            print(f"错误: meta.data中缺少列 '{col}'")
            print(f"可用列: {meta_df.columns.tolist()}")
            sys.exit(1)
    
    # 确保列名匹配
    meta_df = meta_df[required_meta_cols].copy()
    
    print("\n3. 数据预处理...")
    # 去除可能的空格
    for col in ['cell_id', 'sample']:
        adata.obs[col] = adata.obs[col].astype(str).str.strip()
        meta_df[col] = meta_df[col].astype(str).str.strip()
    
    # 创建组合键用于匹配 (cell_id + sample)
    adata.obs['match_key'] = adata.obs['cell_id'].astype(str) + '_' + adata.obs['sample'].astype(str)
    meta_df['match_key'] = meta_df['cell_id'].astype(str) + '_' + meta_df['sample'].astype(str)
    
    print(f"   h5ad中细胞总数: {adata.shape[0]}")
    print(f"   meta.data中细胞总数: {len(meta_df)}")
    
    # 检查重复的匹配键
    h5ad_duplicates = adata.obs['match_key'].duplicated().sum()
    meta_duplicates = meta_df['match_key'].duplicated().sum()
    
    if h5ad_duplicates > 0:
        print(f"  警告: h5ad数据中有 {h5ad_duplicates} 个重复的cell_id+sample组合")
        # 显示前几个重复的匹配键
        duplicates = adata.obs[adata.obs['match_key'].duplicated(keep=False)]['match_key'].unique()[:5]
        print(f"  重复示例: {duplicates.tolist()}")
    
    if meta_duplicates > 0:
        print(f"  警告: meta.data中有 {meta_duplicates} 个重复的cell_id+sample组合")
        # 显示前几个重复的匹配键
        duplicates = meta_df[meta_df['match_key'].duplicated(keep=False)]['match_key'].unique()[:5]
        print(f"  重复示例: {duplicates.tolist()}")
        # 去除meta中的重复，保留第一个
        meta_df = meta_df.drop_duplicates(subset='match_key', keep='first')
        print(f"  去重后meta.data中细胞数: {len(meta_df)}")
    
    print("\n4. 匹配细胞...")
    # 找到匹配的细胞
    matched_keys = set(adata.obs['match_key']).intersection(set(meta_df['match_key']))
    
    print(f"   匹配到的细胞数量: {len(matched_keys)}")
    print(f"   匹配率: {len(matched_keys)/len(meta_df)*100:.2f}% (基于meta.data)")
    print(f"   匹配率: {len(matched_keys)/adata.shape[0]*100:.2f}% (基于h5ad)")
    
    if len(matched_keys) == 0:
        print("错误: 没有匹配到任何细胞，请检查数据")
        # 检查可能的匹配问题
        print("检查前5个细胞的匹配键:")
        print("h5ad前5个:", adata.obs['match_key'].head().tolist())
        print("meta.data前5个:", meta_df['match_key'].head().tolist())
        sys.exit(1)
    
    # 创建匹配索引
    matched_indices = adata.obs['match_key'].isin(matched_keys)
    matched_adata = adata[matched_indices, :].copy()
    
    print(f"   从h5ad中提取的细胞数: {matched_adata.shape[0]}")
    
    # 创建匹配映射字典
    celltype_map = dict(zip(meta_df['match_key'], meta_df['celltype']))
    
    # 添加celltype信息
    matched_adata.obs['celltype'] = matched_adata.obs['match_key'].map(celltype_map)
    
    # 检查是否有缺失的celltype
    missing_celltype = matched_adata.obs['celltype'].isnull().sum()
    if missing_celltype > 0:
        print(f"  警告: 有 {missing_celltype} 个细胞没有匹配到celltype信息")
    
    # 清理临时列
    if 'match_key' in matched_adata.obs.columns:
        matched_adata.obs = matched_adata.obs.drop('match_key', axis=1)
    
    # 检查最终列
    print(f"\n5. 结果数据结构:")
    print(f"   细胞数: {matched_adata.shape[0]}, 基因数: {matched_adata.shape[1]}")
    print(f"   obs列: {sorted(matched_adata.obs.columns.tolist())}")
    
    # 统计celltype分布
    if 'celltype' in matched_adata.obs.columns:
        celltype_counts = matched_adata.obs['celltype'].value_counts()
        print(f"\n   Celltype分布:")
        for i, (ctype, count) in enumerate(celltype_counts.items()):
            percentage = count / len(matched_adata) * 100
            print(f"     {i+1:2d}. {ctype:25s}: {count:5d} 个细胞 ({percentage:6.2f}%)")
    
    # 统计sample分布
    if 'sample' in matched_adata.obs.columns:
        sample_counts = matched_adata.obs['sample'].value_counts()
        print(f"\n   Sample分布:")
        for i, (sample, count) in enumerate(sample_counts.items()):
            percentage = count / len(matched_adata) * 100
            print(f"     {i+1:2d}. {sample:20s}: {count:5d} 个细胞 ({percentage:6.2f}%)")
    
    # 6. 保存结果
    print(f"\n6. 保存结果到 {output_file}...")
    matched_adata.write_h5ad(output_file, compression='gzip')
    print(f"   文件已保存")
    
    # 7. 保存匹配信息
    info_file = "matching_info.txt"
    with open(info_file, 'w') as f:
        f.write("细胞匹配信息\n")
        f.write("=" * 50 + "\n")
        f.write(f"输入h5ad文件: {h5ad_file}\n")
        f.write(f"输入meta文件: {meta_file}\n")
        f.write(f"输出文件: {output_file}\n")
        f.write(f"原始h5ad细胞数: {adata.shape[0]}\n")
        f.write(f"原始meta.data细胞数: {len(meta_df)}\n")
        f.write(f"匹配到的细胞数: {len(matched_keys)}\n")
        f.write(f"匹配率(基于meta.data): {len(matched_keys)/len(meta_df)*100:.2f}%\n")
        f.write(f"匹配率(基于h5ad): {len(matched_keys)/adata.shape[0]*100:.2f}%\n")
        f.write(f"最终细胞数: {matched_adata.shape[0]}\n")
        
        if 'celltype' in matched_adata.obs.columns:
            f.write("\nCelltype分布:\n")
            for ctype, count in matched_adata.obs['celltype'].value_counts().items():
                f.write(f"  {ctype}: {count} 个细胞\n")
        
        if 'sample' in matched_adata.obs.columns:
            f.write("\nSample分布:\n")
            for sample, count in matched_adata.obs['sample'].value_counts().items():
                f.write(f"  {sample}: {count} 个细胞\n")
    
    print(f"\n7. 匹配信息已保存到: {info_file}")
    
    # 8. 保存匹配的细胞列表
    cells_file = "matched_cells.txt"
    matched_cells = matched_adata.obs[['cell_id', 'sample', 'celltype']].copy()
    matched_cells.to_csv(cells_file, sep='\t', index=True)
    print(f"   匹配的细胞列表已保存到: {cells_file}")

if __name__ == "__main__":
    main()
