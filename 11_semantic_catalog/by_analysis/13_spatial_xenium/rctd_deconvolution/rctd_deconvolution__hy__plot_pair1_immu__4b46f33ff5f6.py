import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pathlib import Path

def process_paired_samples_with_specific_celltypes(h5ad_path, patient1, patient2, 
                                                 target_celltypes, annotation_col='annotation', 
                                                 patient_col='patient'):
    """
    处理配对样本数据：提取、统计并可视化指定细胞类型
    """
    # 读取h5ad文件
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"原始数据形状: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    
    # 检查必要的列是否存在
    if patient_col not in adata.obs.columns:
        raise ValueError(f"在adata.obs中未找到'{patient_col}'列")
    
    if annotation_col not in adata.obs.columns:
        raise ValueError(f"在adata.obs中未找到'{annotation_col}'列")
    
    # 检查目标patient是否存在
    target_patients = [patient1, patient2]
    existing_patients = [patient for patient in target_patients if patient in adata.obs[patient_col].values]
    missing_patients = [patient for patient in target_patients if patient not in adata.obs[patient_col].values]
    
    print(f"\n=== Patient匹配情况 ===")
    print(f"目标patient: {target_patients}")
    print(f"找到的patient: {existing_patients}")
    
    if len(existing_patients) < 2:
        print("错误: 没有找到足够的patient!")
        return None
    
    # 提取目标patient数据
    adata_subset = adata[adata.obs[patient_col].isin(existing_patients), :].copy()
    
    print(f"提取后的细胞数: {adata_subset.n_obs}")
    
    # 统计每个patient的细胞数量
    patient_counts = adata_subset.obs[patient_col].value_counts()
    print(f"\n=== 细胞数量统计 ===")
    for patient in existing_patients:
        count = patient_counts[patient]
        print(f"Patient {patient}: {count} 个细胞")
    
    # 检查指定的细胞类型是否存在
    all_celltypes = adata_subset.obs[annotation_col].unique()
    existing_celltypes = [ct for ct in target_celltypes if ct in all_celltypes]
    
    print(f"找到的指定细胞类型: {len(existing_celltypes)} 种")
    
    # 修复关键错误：先添加'Others'类别
    if 'celltype_grouped' not in adata_subset.obs.columns:
        # 确保annotation列是分类类型
        if not pd.api.types.is_categorical_dtype(adata_subset.obs[annotation_col]):
            adata_subset.obs[annotation_col] = adata_subset.obs[annotation_col].astype('category')
        
        # 创建新的分类列，包含所有现有类别加上'Others'
        all_categories = list(adata_subset.obs[annotation_col].cat.categories) + ['Others']
        
        # 创建分组列
        adata_subset.obs['celltype_grouped'] = adata_subset.obs[annotation_col].copy()
        
        # 添加'Others'类别
        adata_subset.obs['celltype_grouped'] = adata_subset.obs['celltype_grouped'].cat.add_categories(['Others'])
    
    # 现在安全地设置'Others'值
    other_mask = ~adata_subset.obs[annotation_col].isin(existing_celltypes)
    adata_subset.obs.loc[other_mask, 'celltype_grouped'] = 'Others'
    
    # 计算细胞类型分布
    celltype_counts = pd.crosstab(adata_subset.obs[patient_col], adata_subset.obs['celltype_grouped'])
    
    # 确保patient顺序
    celltype_counts = celltype_counts.reindex([patient1, patient2])
    
    # 转换为百分比
    celltype_percent = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100
    
    print(f"\n细胞类型分布(百分比):")
    print(celltype_percent)
    
    # 设置图形样式
    plt.style.use('default')
    plt.rcParams['font.family'] = 'DejaVu Sans'
    
    # 创建图形
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # 准备颜色
    n_celltypes = len(celltype_percent.columns)
    colors = plt.cm.Set3(np.linspace(0, 1, n_celltypes))
    
    # 创建堆积图
    bottom = np.zeros(len(celltype_percent))
    bar_width = 0.7
    
    for i, celltype in enumerate(celltype_percent.columns):
        values = celltype_percent[celltype].values
        bars = ax.bar(range(len(celltype_percent)), values, bottom=bottom, 
                      width=bar_width, label=celltype, color=colors[i], 
                      alpha=0.8, edgecolor='white', linewidth=1)
        bottom += values
        
        # 添加百分比标签
        for j, (val, bar) in enumerate(zip(values, bars)):
            if val > 3:
                height = bar.get_height()
                y_pos = bar.get_y() + height / 2
                ax.text(bar.get_x() + bar.get_width()/2, y_pos, 
                        f'{val:.1f}%', ha='center', va='center', 
                        fontsize=9, fontweight='bold')
    
    # 设置图表属性
    ax.set_title('Target cell types', fontsize=16, fontweight='bold', pad=20)
    ax.set_ylabel('Proportion (%)', fontsize=14)
    ax.set_xlabel('Patient', fontsize=14)
    ax.set_xticks(range(len(celltype_percent)))
    ax.set_xticklabels(celltype_percent.index, rotation=0)
    ax.grid(axis='y', alpha=0.3, linestyle='--', zorder=0)
    ax.set_ylim(0, 100)
    
    # 添加图例
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
              borderaxespad=0., fontsize=10, title='Cell Type')
    
    plt.tight_layout()
    
    # 保存图片
    output_filename = f"Immu_celltype_comparison_{patient1}_vs_{patient2}.png"
    output_filename1 = f"Immu_celltype_comparison_{patient1}_vs_{patient2}.pdf"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.savefig(output_filename1, dpi=300, bbox_inches='tight')
    print(f"\n图表已保存为: {output_filename}")
    
    # 显示详细统计信息
    print(f"\n=== 详细统计信息 ===")
    for patient in existing_patients:
        patient_data = adata_subset[adata_subset.obs[patient_col] == patient]
        celltypes = patient_data.obs['celltype_grouped'].value_counts()
        total_cells = len(patient_data)
        
        print(f"\nPatient {patient}:")
        print(f"总细胞数: {total_cells}")
        print("细胞类型分布:")
        for celltype, count in celltypes.items():
            percentage = (count / total_cells) * 100
            print(f"  {celltype}: {count} 个细胞 ({percentage:.2f}%)")
    
    plt.show()
    
    # 保存提取后的数据
    output_h5ad = "paired_1_immu.h5ad"
    
    # 确保保留raw矩阵
    if adata_subset.raw is None:
        adata_subset.raw = adata_subset
    
    adata_subset.write(output_h5ad)
    print(f"\n提取后的数据已保存为: {output_h5ad}")
    
    return adata_subset, celltype_percent

# 使用示例
if __name__ == "__main__":
    # 文件路径和patient名称
    h5ad_path = "adata_celltype_RCTD_removed_with_patient.h5ad"  # 替换为您的文件路径
    patient1 = "Z1944219-14"
    patient2 = "Z2489020-13"
    
    # 指定要分析的细胞类型
    target_celltypes = ["CD8+_Tem","M2","MKI67+_T","B_follicular","CD4+_Tex","LYVE1+_TAM","CD8+_Tcm","NK","IFN-T","B_germinal_center","EMT_M2","Treg","CD8+_Tex","CD4+_Tem","B_memory","TAM"]
    
    try:
        result = process_paired_samples_with_specific_celltypes(
            h5ad_path, patient1, patient2, target_celltypes
        )
        
        if result is not None:
            print(f"\n处理完成!")
            
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()
