import scanpy as sc
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def add_group_info_and_save(h5ad_path, before_patients, after_patients, 
                           before_group_name='MPC', after_group_name='TNBC',
                           patient_col='patient', output_h5ad="MPCvsTNBC_all_sample_ct.h5ad"):
    """
    为h5ad数据添加分组信息并保存
    
    参数:
    h5ad_path: h5ad文件路径
    before_patients: MPC组patient列表
    after_patients: TNBC组patient列表
    before_group_name: MPC组的名称
    after_group_name: TNBC组的名称
    patient_col: patient列名
    output_h5ad: 输出h5ad文件名
    """
    
    # 读取h5ad文件
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"原始数据形状: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    
    # 检查必要的列是否存在
    if patient_col not in adata.obs.columns:
        raise ValueError(f"在adata.obs中未找到'{patient_col}'列")
    
    # 合并所有目标patient
    all_patients = before_patients + after_patients
    existing_patients = [p for p in all_patients if p in adata.obs[patient_col].values]
    missing_patients = [p for p in all_patients if p not in adata.obs[patient_col].values]
    
    print(f"\n=== Patient匹配情况 ===")
    print(f"目标patient数量: {len(all_patients)}")
    print(f"找到的patient数量: {len(existing_patients)}")
    print(f"未找到的patient: {missing_patients}")
    
    if len(existing_patients) < 2:
        print("警告: 找到的patient数量较少!")
    
    # 提取目标patient数据
    adata_subset = adata[adata.obs[patient_col].isin(existing_patients), :].copy()
    
    print(f"\n=== 数据提取结果 ===")
    print(f"提取后的细胞数: {adata_subset.n_obs}")
    
    # 添加分组信息
    adata_subset.obs['group'] = 'unknown'
    adata_subset.obs.loc[adata_subset.obs[patient_col].isin(before_patients), 'group'] = before_group_name
    adata_subset.obs.loc[adata_subset.obs[patient_col].isin(after_patients), 'group'] = after_group_name
    
    # 统计分组情况
    group_counts = adata_subset.obs['group'].value_counts()
    print(f"\n=== 分组统计 ===")
    for group, count in group_counts.items():
        print(f"{group}组: {count} 个细胞")
    
    # 统计每个patient的细胞数量
    patient_counts = adata_subset.obs[patient_col].value_counts()
    print(f"\n=== 各Patient细胞数量统计 ===")
    for patient in existing_patients:
        count = patient_counts[patient]
        group = before_group_name if patient in before_patients else after_group_name
        print(f"Patient {patient} ({group}组): {count} 个细胞")
    
    # 保存数据
    if adata_subset.raw is None:
        adata_subset.raw = adata_subset
    
    adata_subset.write(output_h5ad, compression='gzip')
    print(f"\n数据已保存为: {output_h5ad}")
    
    # 创建分组统计表
    create_group_summary(adata_subset, before_group_name, after_group_name, patient_col)
    
    return adata_subset

def create_group_summary(adata, before_group_name, after_group_name, patient_col):
    """
    创建分组统计表
    """
    print(f"\n=== 分组详细统计 ===")
    
    # 获取分组信息
    groups = adata.obs['group'].unique()
    
    for group in groups:
        group_data = adata[adata.obs['group'] == group]
        patients = group_data.obs[patient_col].unique()
        
        print(f"\n{group}组 ({len(patients)}个patient):")
        for patient in patients:
            patient_data = group_data[group_data.obs[patient_col] == patient]
            print(f"  Patient {patient}: {len(patient_data)} 个细胞")
    
    # 保存分组信息到CSV
    summary_df = adata.obs[['group', patient_col]].drop_duplicates()
    summary_df = summary_df.sort_values(['group', patient_col])
    summary_df.to_csv("group_assignment.csv", index=False)
    print(f"\n分组信息已保存为: group_assignment.csv")
    
    # 创建细胞类型统计
    if 'annotation' in adata.obs.columns:
        celltype_stats = adata.obs.groupby(['group', 'annotation']).size().reset_index(name='count')
        celltype_stats = celltype_stats.pivot_table(index='annotation', columns='group', values='count', fill_value=0)
        celltype_stats.to_csv("celltype_counts_by_group.csv")
        print(f"细胞类型计数已保存为: celltype_counts_by_group.csv")
    
    return summary_df

# 使用示例
if __name__ == "__main__":
    # 文件路径
    h5ad_path = "adata_celltype_RCTD_removed_with_patient.h5ad"  # 请替换为您的h5ad文件路径
    
    # 定义分组
    before_patients = ["Z2038176-15","Q2324341-1","Z1530765-11","F1706967-11","F1805049-2","Z1915172-7","Z1954693-9","ZF2207170-8","ZF2202916-8","ZF2208790-1","Z2262137-5","Z2271577-12","ZF2304103-8","Z2323818-8","Z2336500-6","Z2345459-7","Z2332034-5","Z2417954-8","Q2236500-1","Q2406275-2","Z2101519-1","Z2123354-1","Z1642635-14","Z2224197-10","Z1736838-15","Z1944219-14","Z1430967-9","Z2057500-35","Z2171766-3","Z2353385-10","Z1635989-12","Z1814675-2","Z1826868-6","Z1952876-9","Z1802235-10"]
    after_patients = ["Z2279427-2","H2510098-1","H2510098-2","Q2243285-1","Z2370717-2","Q2337205-2","Q2233699-2","Q2239242-1","Q2319383-1","Q2322393-1","Q2337350-1"]
    
    try:
        adata_with_groups = add_group_info_and_save(
            h5ad_path, 
            before_patients, 
            after_patients,
            before_group_name='MPC',
            after_group_name='TNBC',
            output_h5ad="MPCvsTNBC_all_sample_ct.h5ad"
        )
        
        # 打印最终结果摘要
        print(f"\n=== 最终结果摘要 ===")
        print(f"处理后的数据形状: {adata_with_groups.n_obs} 个细胞, {adata_with_groups.n_vars} 个基因")
        
        # 检查分组是否正确
        group_counts = adata_with_groups.obs['group'].value_counts()
        print(f"\n分组细胞数量:")
        for group, count in group_counts.items():
            print(f"  {group}: {count} 个细胞 ({count/adata_with_groups.n_obs*100:.1f}%)")
        
        # 检查patient数量
        patient_counts_by_group = adata_with_groups.obs.groupby('group')['patient'].nunique()
        print(f"\n分组patient数量:")
        for group, count in patient_counts_by_group.items():
            print(f"  {group}: {count} 个patient")
        
        print(f"\n处理完成!")
        
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()
