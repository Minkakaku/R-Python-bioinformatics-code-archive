import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys
from pathlib import Path

def main():
    # 定义要处理的h5ad文件列表
    h5ad_files = [
        "B_cells.h5ad",
        "Endothelials.h5ad", 
        "Epithelials.h5ad",
        "Fibroblasts.h5ad",
        "Macrophages.h5ad",
        "Neuron-like.h5ad",
        "Neutrophils.h5ad",
        "Plasmas.h5ad",
        "T_cells.h5ad"
    ]
    
    # 检查文件是否存在
    missing_files = []
    for file in h5ad_files:
        if not Path(file).exists():
            missing_files.append(file)
    
    if missing_files:
        print(f"错误: 以下文件不存在: {missing_files}")
        sys.exit(1)
    
    print("开始从多个h5ad文件中提取元数据...")
    print("="*60)
    
    # 存储所有提取的数据
    all_meta_data = []
    
    # 定义要提取的列
    columns_to_extract = ['cell_id', 'sample', 'group', 'celltype']
    
    for i, h5ad_file in enumerate(h5ad_files, 1):
        print(f"处理文件 {i}/{len(h5ad_files)}: {h5ad_file}")
        
        try:
            # 读取h5ad文件
            adata = sc.read_h5ad(h5ad_file)
            print(f"  数据形状: {adata.shape}")
            print(f"  obs列: {adata.obs.columns.tolist()}")
            
            # 检查所需的列是否存在
            existing_columns = [col for col in columns_to_extract if col in adata.obs.columns]
            missing_columns = [col for col in columns_to_extract if col not in adata.obs.columns]
            
            if missing_columns:
                print(f"  警告: 缺少以下列: {missing_columns}")
                print(f"  将用NaN填充缺失的列")
            
            # 提取存在的列
            if existing_columns:
                # 创建数据副本
                meta_df = adata.obs[existing_columns].copy()
                
                # 添加缺失的列并用NaN填充
                for col in missing_columns:
                    meta_df[col] = np.nan
                
                # 确保列的顺序一致
                meta_df = meta_df[columns_to_extract]
                
                # 添加来源文件信息
                meta_df['source_file'] = h5ad_file
                
                all_meta_data.append(meta_df)
                print(f"  成功提取 {len(meta_df)} 个细胞的元数据")
            else:
                print(f"  错误: 文件中没有找到任何所需的列")
                continue
                
        except Exception as e:
            print(f"  处理文件 {h5ad_file} 时出错: {e}")
            continue
    
    if not all_meta_data:
        print("错误: 没有成功提取任何数据")
        sys.exit(1)
    
    # 合并所有数据
    print(f"\n合并所有提取的数据...")
    combined_meta = pd.concat(all_meta_data, ignore_index=True)
    
    print(f"合并后的数据形状: {combined_meta.shape}")
    print(f"总细胞数: {len(combined_meta)}")
    
    # 检查列
    print(f"\n数据列: {combined_meta.columns.tolist()}")
    
    # 检查重复的cell_id
    duplicate_cell_ids = combined_meta['cell_id'].duplicated().sum()
    if duplicate_cell_ids > 0:
        print(f"警告: 发现 {duplicate_cell_ids} 个重复的cell_id")
        # 显示前几个重复的cell_id
        duplicates = combined_meta[combined_meta['cell_id'].duplicated(keep=False)]
        print(f"重复的cell_id示例: {duplicates['cell_id'].unique()[:10].tolist()}")
    else:
        print("cell_id列中没有重复值")
    
    # 统计每个来源文件的细胞数
    print(f"\n各文件贡献的细胞数:")
    file_counts = combined_meta['source_file'].value_counts()
    for file, count in file_counts.items():
        print(f"  {file}: {count} 个细胞 ({count/len(combined_meta):.2%})")
    
    # 统计celltype分布
    if 'celltype' in combined_meta.columns:
        print(f"\nCelltype分布:")
        celltype_counts = combined_meta['celltype'].value_counts(dropna=False)
        for ctype, count in celltype_counts.items():
            percentage = count / len(combined_meta) * 100
            ctype_name = 'NaN' if pd.isna(ctype) else ctype
            print(f"  {ctype_name:20s}: {count:6d} 个细胞 ({percentage:6.2f}%)")
    
    # 统计group分布
    if 'group' in combined_meta.columns:
        print(f"\nGroup分布:")
        group_counts = combined_meta['group'].value_counts(dropna=False)
        for group, count in group_counts.items():
            percentage = count / len(combined_meta) * 100
            group_name = 'NaN' if pd.isna(group) else group
            print(f"  {group_name:20s}: {count:6d} 个细胞 ({percentage:6.2f}%)")
    
    # 统计sample分布
    if 'sample' in combined_meta.columns:
        print(f"\nSample分布:")
        sample_counts = combined_meta['sample'].value_counts(dropna=False)
        for sample, count in sample_counts.items():
            percentage = count / len(combined_meta) * 100
            sample_name = 'NaN' if pd.isna(sample) else sample
            print(f"  {sample_name:20s}: {count:6d} 个细胞 ({percentage:6.2f}%)")
    
    # 保存合并的元数据
    output_file = "meta.data"
    print(f"\n保存合并的元数据到 {output_file}...")
    
    # 保存为制表符分隔的文件
    combined_meta.to_csv(output_file, sep='\t', index=False, encoding='utf-8-sig')
    print(f"元数据已保存到: {output_file}")
    
    # 同时保存为CSV格式
    csv_file = "meta.csv"
    combined_meta.to_csv(csv_file, index=False, encoding='utf-8-sig')
    print(f"元数据(CSV格式)已保存到: {csv_file}")
    
    # 保存统计信息
    stats_file = "meta_statistics.txt"
    with open(stats_file, 'w', encoding='utf-8') as f:
        f.write("元数据统计信息\n")
        f.write("="*60 + "\n\n")
        
        f.write(f"总细胞数: {len(combined_meta)}\n")
        f.write(f"数据列: {', '.join(combined_meta.columns.tolist())}\n")
        f.write(f"来源文件数: {len(file_counts)}\n\n")
        
        f.write("各文件贡献的细胞数:\n")
        f.write("-"*40 + "\n")
        for file, count in file_counts.items():
            f.write(f"{file}: {count} 个细胞 ({count/len(combined_meta):.2%})\n")
        
        if 'celltype' in combined_meta.columns:
            f.write(f"\nCelltype分布:\n")
            f.write("-"*40 + "\n")
            for ctype, count in celltype_counts.items():
                percentage = count / len(combined_meta) * 100
                ctype_name = 'NaN' if pd.isna(ctype) else ctype
                f.write(f"{ctype_name}: {count} 个细胞 ({percentage:.2f}%)\n")
        
        if 'group' in combined_meta.columns:
            f.write(f"\nGroup分布:\n")
            f.write("-"*40 + "\n")
            for group, count in group_counts.items():
                percentage = count / len(combined_meta) * 100
                group_name = 'NaN' if pd.isna(group) else group
                f.write(f"{group_name}: {count} 个细胞 ({percentage:.2f}%)\n")
        
        if 'sample' in combined_meta.columns:
            f.write(f"\nSample分布:\n")
            f.write("-"*40 + "\n")
            for sample, count in sample_counts.items():
                percentage = count / len(combined_meta) * 100
                sample_name = 'NaN' if pd.isna(sample) else sample
                f.write(f"{sample_name}: {count} 个细胞 ({percentage:.2f}%)\n")
        
        f.write(f"\n前5行数据预览:\n")
        f.write("-"*40 + "\n")
        f.write(combined_meta.head().to_string())
        f.write("\n")
    
    print(f"统计信息已保存到: {stats_file}")
    
    # 创建各列的缺失值统计
    missing_stats = combined_meta.isnull().sum()
    missing_percentage = (missing_stats / len(combined_meta) * 100).round(2)
    
    missing_df = pd.DataFrame({
        'column': missing_stats.index,
        'missing_count': missing_stats.values,
        'missing_percentage': missing_percentage.values
    })
    
    missing_file = "missing_values_statistics.csv"
    missing_df.to_csv(missing_file, index=False, encoding='utf-8-sig')
    print(f"缺失值统计已保存到: {missing_file}")
    
    print(f"\n缺失值统计:")
    for col in columns_to_extract:
        if col in combined_meta.columns:
            missing = combined_meta[col].isnull().sum()
            percentage = missing / len(combined_meta) * 100
            print(f"  {col:15s}: {missing:6d} 个缺失值 ({percentage:6.2f}%)")
    
    print("\n处理完成!")

if __name__ == "__main__":
    main()
