import scanpy as sc
import pandas as pd
import numpy as np
import sys
from pathlib import Path

def main():
    # 文件路径
    input_file = "All_raw.h5ad"
    output_file = "All_filter.h5ad"
    
    # 检查输入文件是否存在
    if not Path(input_file).exists():
        print(f"错误: 文件 {input_file} 不存在！")
        sys.exit(1)
    
    print("1. 读取h5ad数据...")
    adata = sc.read_h5ad(input_file)
    print(f"   原始数据形状: {adata.shape}")
    print(f"   细胞数: {adata.shape[0]}, 基因数: {adata.shape[1]}")
    
    # 检查是否有celltype列
    if 'celltype' not in adata.obs.columns:
        print("错误: 数据中缺少 'celltype' 列")
        print(f"可用列: {adata.obs.columns.tolist()}")
        sys.exit(1)
    
    # 检查是否有sample列
    if 'sample' not in adata.obs.columns:
        print("错误: 数据中缺少 'sample' 列，无法进行group_final和patient_final映射")
        print(f"可用列: {adata.obs.columns.tolist()}")
        sys.exit(1)
    
    # 2. 添加group_final和patient_final列
    print("\n2. 添加group_final和patient_final列...")
    
    # 定义样本到group_final的映射关系
    sample_to_group_final = {
        "PRE_5-1": "PRE-R", "PRE_5-2": "PRE-R", "PRE_6-1": "PRE-R", 
        "PRE_6-2": "PRE-R", "PRE_6-3": "PRE-R", "PRE_8": "PRE-R", 
        "PRE_15-1": "PRE-R", "PRE_15-2": "PRE-R", "PRE_15-3": "PRE-R", 
        "PRE_15-4": "PRE-R", "PRE_16": "PRE-R", "PRE_17-1": "PRE-R", 
        "PRE_17-2": "PRE-R", "PRE_17-3": "PRE-R", "PRE_18-1": "PRE-R", 
        "PRE_18-2": "PRE-R", "PRE_19": "PRE-R", "PRE_3-1": "PRE-NR", 
        "PRE_3-2": "PRE-NR", "PRE_10-1": "PRE-NR", "PRE_10-2": "PRE-NR", 
        "PRE_20-1": "PRE-NR", "PRE_20-2": "PRE-NR", "PRE_20-3": "PRE-NR", 
        "PRE_20-4": "PRE-NR", "PRE_21": "PRE-NR", "PRE_22-1": "PRE-NR", 
        "PRE_22-2": "PRE-NR", "PRE_23-1": "PRE-NR", "PRE_23-2": "PRE-NR", 
        "PRE_25-1": "PRE-NR", "PRE_25-2": "PRE-NR", "PRE_25-3": "PRE-NR", 
        "PRE_25-4": "PRE-NR", "PRE_26-1": "PRE-NR", "PRE_26-3": "PRE-NR", 
        "PRE_27-1": "PRE-NR", "PRE_27-2": "PRE-NR", "POST_5": "POST-R", 
        "POST_6": "POST-R", "POST_8": "POST-R", "POST_15": "POST-R", 
        "POST_16": "POST-R", "POST_18": "POST-R", "POST_19": "POST-R", 
        "POST_3": "POST-NR", "POST_20": "POST-NR", "POST_21": "POST-NR", 
        "POST_22-1": "POST-NR", "POST_25": "POST-NR", "POST_26": "POST-NR", 
        "POST_27": "POST-NR"
    }
    
    # 定义样本到patient_final的映射关系
    sample_to_patient_final = {
        "POST_3": "P3", "POST_8": "P8", "PRE_8": "P8", "POST_6": "P6", 
        "PRE_3-1": "P3", "PRE_3-2": "P3", "PRE_6-1": "P6", "PRE_6-2": "P6", 
        "PRE_5-1": "P5", "PRE_5-2": "P5", "PRE_10-1": "P10", "PRE_10-2": "P10", 
        "POST_5": "P5", "POST_16": "P16", "POST_19": "P19", "POST_22-1": "P22", 
        "POST_22-2": "P22", "PRE_15-1": "P15", "PRE_15-2": "P15", "PRE_15-3": "P15", 
        "PRE_15-4": "P15", "PRE_16": "P16", "PRE_17-1": "P17", "PRE_17-2": "P17", 
        "PRE_17-3": "P17", "PRE_20-1": "P20", "PRE_20-2": "P20", "PRE_21": "P21", 
        "PRE_22-1": "P22", "PRE_22-2": "P22", "PRE_23-1": "P23", "PRE_23-2": "P23", 
        "PRE_24-1": "P24", "PRE_24-2": "P24", "PRE_26-1": "P26", "PRE_26-2": "P26", 
        "PRE_26-3": "P26", "POST_15": "P15", "POST_18": "P18", "POST_20": "P20", 
        "POST_21": "P21", "POST_25": "P25", "POST_26": "P26", "POST_27": "P27", 
        "PRE_18-1": "P18", "PRE_18-2": "P18", "PRE_19": "P19", "PRE_25-1": "P25", 
        "PRE_25-2": "P25", "PRE_25-3": "P25", "PRE_25-4": "P25", "PRE_27-1": "P27", 
        "PRE_27-2": "P27", "PRE_6-3": "P3", "PRE_20-3": "P20", "PRE_20-4": "P20"
    }
    
    # 应用映射
    adata.obs['group_final'] = adata.obs['sample'].map(sample_to_group_final)
    adata.obs['patient_final'] = adata.obs['sample'].map(sample_to_patient_final)
    
    # 检查未映射的样本
    unmapped_group = adata.obs[adata.obs['group_final'].isna()]['sample'].unique()
    unmapped_patient = adata.obs[adata.obs['patient_final'].isna()]['sample'].unique()
    
    if len(unmapped_group) > 0:
        print(f"   警告: 以下样本未找到group_final映射关系: {unmapped_group.tolist()}")
        print(f"   未映射的细胞数: {(adata.obs['group_final'].isna()).sum()}")
    
    if len(unmapped_patient) > 0:
        print(f"   警告: 以下样本未找到patient_final映射关系: {unmapped_patient.tolist()}")
        print(f"   未映射的细胞数: {(adata.obs['patient_final'].isna()).sum()}")
    
    # 统计group_final分布
    if 'group_final' in adata.obs.columns:
        group_final_counts = adata.obs['group_final'].value_counts(dropna=False)
        print(f"\n   Group_final分布:")
        for i, (group, count) in enumerate(group_final_counts.items()):
            group_name = 'NaN' if pd.isna(group) else group
            print(f"     {i+1:2d}. {group_name:20s}: {count:6d} 个细胞 ({count/adata.shape[0]:.2%})")
    
    # 统计patient_final分布
    if 'patient_final' in adata.obs.columns:
        patient_final_counts = adata.obs['patient_final'].value_counts(dropna=False)
        print(f"\n   Patient_final分布:")
        for i, (patient, count) in enumerate(patient_final_counts.items()):
            patient_name = 'NaN' if pd.isna(patient) else patient
            print(f"     {i+1:2d}. {patient_name:20s}: {count:6d} 个细胞 ({count/adata.shape[0]:.2%})")
    
    # 3. 修改celltype列
    print("\n3. 修改celltype列...")
    print(f"   原始celltype类型: {adata.obs['celltype'].nunique()} 种")
    print(f"   原始celltype分布:")
    celltype_counts = adata.obs['celltype'].value_counts()
    for i, (ctype, count) in enumerate(celltype_counts.items()):
        print(f"     {i+1:2d}. {ctype:30s}: {count:6d} 个细胞 ({count/adata.shape[0]:.2%})")
    
    # 定义需要修改的映射关系
    celltype_mapping = {
        'M2-like TRM': 'M2',
        'DC': 'LAMP3+DC',
        'cDC1': 'cDC'
    }
    
    # 备份原始celltype
    # adata.obs['celltype_old'] = adata.obs['celltype'].copy()
    
    # 应用映射
    adata.obs['celltype'] = adata.obs['celltype'].map(lambda x: celltype_mapping.get(x, x))
    
    print(f"   修改后celltype类型: {adata.obs['celltype'].nunique()} 种")
    print(f"   修改后的celltype分布:")
    celltype_counts_after = adata.obs['celltype'].value_counts()
    for i, (ctype, count) in enumerate(celltype_counts_after.items()):
        print(f"     {i+1:2d}. {ctype:30s}: {count:6d} 个细胞 ({count/adata.shape[0]:.2%})")
    
    # 4. 剔除指定celltype类型的细胞
    print("\n4. 剔除指定celltype类型的细胞...")
    
    # 定义需要剔除的细胞类型列表
    remove_celltypes = [
        'B-CAF',
        'B-Epi', 
        'B-Myeloid',
        'B-T',
        'Endo-CAF',
        'Endo-Plasma',
        'Macro-Epi',
        'Macro-CAF',
        'Macro-Mkl67+ Epi',
        'Schwann ike',
        'T-CAF',
        'T-Endo',
        'T-Macro',
        'T-Epi',
        'Endo like'
    ]
    
    # 检查要剔除的细胞类型是否存在
    existing_remove_types = [ctype for ctype in remove_celltypes if ctype in adata.obs['celltype'].unique()]
    non_existing_types = [ctype for ctype in remove_celltypes if ctype not in adata.obs['celltype'].unique()]
    
    if existing_remove_types:
        print(f"   将剔除以下细胞类型: {existing_remove_types}")
        
        # 统计要剔除的细胞数量
        remove_mask = adata.obs['celltype'].isin(existing_remove_types)
        cells_to_remove = remove_mask.sum()
        
        print(f"   将剔除的细胞数: {cells_to_remove} ({cells_to_remove/adata.shape[0]:.2%})")
        print(f"   剔除前细胞数: {adata.shape[0]}")
        
        # 剔除这些细胞
        adata_filtered = adata[~remove_mask, :].copy()
        
        print(f"   剔除后细胞数: {adata_filtered.shape[0]}")
        print(f"   保留的细胞比例: {adata_filtered.shape[0]/adata.shape[0]:.2%}")
        
        # 统计每个被剔除类型的细胞数
        print(f"\n   被剔除的细胞类型分布:")
        for ctype in existing_remove_types:
            count = (adata.obs['celltype'] == ctype).sum()
            if count > 0:
                print(f"     {ctype:30s}: {count:6d} 个细胞 ({count/adata.shape[0]:.2%})")
    
    if non_existing_types:
        print(f"\n   注意: 以下指定的剔除类型在数据中不存在: {non_existing_types}")
    
    # 如果没有要剔除的类型，则使用原始数据
    if not existing_remove_types:
        print("   没有找到需要剔除的细胞类型，使用原始数据")
        adata_filtered = adata.copy()
    
    # 5. 检查过滤后的数据
    print(f"\n5. 过滤后数据统计:")
    print(f"   细胞数: {adata_filtered.shape[0]}, 基因数: {adata_filtered.shape[1]}")
    print(f"   celltype类型: {adata_filtered.obs['celltype'].nunique()} 种")
    print(f"   过滤后celltype分布:")
    
    filtered_celltype_counts = adata_filtered.obs['celltype'].value_counts()
    for i, (ctype, count) in enumerate(filtered_celltype_counts.items()):
        print(f"     {i+1:2d}. {ctype:30s}: {count:6d} 个细胞 ({count/adata_filtered.shape[0]:.2%})")
    
    # 统计过滤后的group_final分布
    if 'group_final' in adata_filtered.obs.columns:
        filtered_group_counts = adata_filtered.obs['group_final'].value_counts(dropna=False)
        print(f"\n   过滤后group_final分布:")
        for i, (group, count) in enumerate(filtered_group_counts.items()):
            group_name = 'NaN' if pd.isna(group) else group
            print(f"     {i+1:2d}. {group_name:20s}: {count:6d} 个细胞 ({count/adata_filtered.shape[0]:.2%})")
    
    # 统计过滤后的patient_final分布
    if 'patient_final' in adata_filtered.obs.columns:
        filtered_patient_counts = adata_filtered.obs['patient_final'].value_counts(dropna=False)
        print(f"\n   过滤后patient_final分布:")
        for i, (patient, count) in enumerate(filtered_patient_counts.items()):
            patient_name = 'NaN' if pd.isna(patient) else patient
            print(f"     {i+1:2d}. {patient_name:20s}: {count:6d} 个细胞 ({count/adata_filtered.shape[0]:.2%})")
    
    # 6. 保存过滤后的数据
    print(f"\n6. 保存数据到 {output_file}...")
    adata_filtered.write_h5ad(output_file, compression='gzip')
    print(f"   文件已保存")
    
    # 7. 保存处理信息
    info_file = "filtering_info.txt"
    with open(info_file, 'w') as f:
        f.write("细胞过滤信息\n")
        f.write("=" * 50 + "\n")
        f.write(f"输入文件: {input_file}\n")
        f.write(f"输出文件: {output_file}\n")
        f.write(f"原始细胞数: {adata.shape[0]}\n")
        f.write(f"过滤后细胞数: {adata_filtered.shape[0]}\n")
        f.write(f"剔除细胞数: {adata.shape[0] - adata_filtered.shape[0]}\n")
        f.write(f"保留比例: {adata_filtered.shape[0]/adata.shape[0]:.2%}\n\n")
        
        f.write("celltype修改映射:\n")
        for old, new in celltype_mapping.items():
            f.write(f"  {old} -> {new}\n")
        
        f.write("\n剔除的细胞类型:\n")
        for ctype in remove_celltypes:
            f.write(f"  {ctype}\n")
        
        f.write("\n原始celltype分布:\n")
        for ctype, count in celltype_counts.items():
            f.write(f"  {ctype}: {count} 个细胞 ({count/adata.shape[0]:.2%})\n")
        
        f.write("\n修改后celltype分布:\n")
        for ctype, count in celltype_counts_after.items():
            f.write(f"  {ctype}: {count} 个细胞 ({count/adata.shape[0]:.2%})\n")
        
        f.write("\n过滤后celltype分布:\n")
        for ctype, count in filtered_celltype_counts.items():
            f.write(f"  {ctype}: {count} 个细胞 ({count/adata_filtered.shape[0]:.2%})\n")
        
        f.write("\n被剔除的细胞类型统计:\n")
        for ctype in existing_remove_types:
            count = (adata.obs['celltype'] == ctype).sum()
            if count > 0:
                f.write(f"  {ctype}: {count} 个细胞 ({count/adata.shape[0]:.2%})\n")
        
        f.write("\n原始group_final分布:\n")
        for group, count in adata.obs['group_final'].value_counts(dropna=False).items():
            group_name = 'NaN' if pd.isna(group) else group
            f.write(f"  {group_name}: {count} 个细胞 ({count/adata.shape[0]:.2%})\n")
        
        f.write("\n原始patient_final分布:\n")
        for patient, count in adata.obs['patient_final'].value_counts(dropna=False).items():
            patient_name = 'NaN' if pd.isna(patient) else patient
            f.write(f"  {patient_name}: {count} 个细胞 ({count/adata.shape[0]:.2%})\n")
        
        f.write("\n过滤后group_final分布:\n")
        for group, count in adata_filtered.obs['group_final'].value_counts(dropna=False).items():
            group_name = 'NaN' if pd.isna(group) else group
            f.write(f"  {group_name}: {count} 个细胞 ({count/adata_filtered.shape[0]:.2%})\n")
        
        f.write("\n过滤后patient_final分布:\n")
        for patient, count in adata_filtered.obs['patient_final'].value_counts(dropna=False).items():
            patient_name = 'NaN' if pd.isna(patient) else patient
            f.write(f"  {patient_name}: {count} 个细胞 ({count/adata_filtered.shape[0]:.2%})\n")
    
    print(f"\n7. 处理信息已保存到: {info_file}")
    
    # 8. 创建细胞类型统计表
    stats_df = pd.DataFrame({
        'celltype': filtered_celltype_counts.index,
        'count': filtered_celltype_counts.values,
        'percentage': (filtered_celltype_counts.values / adata_filtered.shape[0] * 100).round(2)
    })
    stats_df = stats_df.sort_values('count', ascending=False)
    stats_df.to_csv("filtered_celltype_statistics.csv", index=False, encoding='utf-8-sig')
    print(f"   细胞类型统计已保存到: filtered_celltype_statistics.csv")
    
    # 9. 创建分组统计表
    if 'group_final' in adata_filtered.obs.columns:
        group_stats = adata_filtered.obs['group_final'].value_counts(dropna=False)
        group_stats_df = pd.DataFrame({
            'group_final': group_stats.index,
            'count': group_stats.values,
            'percentage': (group_stats.values / adata_filtered.shape[0] * 100).round(2)
        })
        group_stats_df = group_stats_df.sort_values('count', ascending=False)
        group_stats_df.to_csv("filtered_group_final_statistics.csv", index=False, encoding='utf-8-sig')
        print(f"   group_final统计已保存到: filtered_group_final_statistics.csv")
    
    if 'patient_final' in adata_filtered.obs.columns:
        patient_stats = adata_filtered.obs['patient_final'].value_counts(dropna=False)
        patient_stats_df = pd.DataFrame({
            'patient_final': patient_stats.index,
            'count': patient_stats.values,
            'percentage': (patient_stats.values / adata_filtered.shape[0] * 100).round(2)
        })
        patient_stats_df = patient_stats_df.sort_values('count', ascending=False)
        patient_stats_df.to_csv("filtered_patient_final_statistics.csv", index=False, encoding='utf-8-sig')
        print(f"   patient_final统计已保存到: filtered_patient_final_statistics.csv")
    
    print("\n处理完成!")

if __name__ == "__main__":
    main()
