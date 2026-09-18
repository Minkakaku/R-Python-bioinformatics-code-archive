import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pathlib import Path

def process_paired_samples(h5ad_path, patient1, patient2, annotation_col='annotation', patient_col='patient'):
    """
    处理配对样本数据：提取、统计并可视化
    
    参数:
    h5ad_path: h5ad文件路径
    patient1: 第一个patient名称
    patient2: 第二个patient名称
    annotation_col: 细胞类型所在的列名
    patient_col: patient信息所在的列名
    """
    
    # 读取h5ad文件
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"原始数据形状: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    
    # 检查必要的列是否存在
    if patient_col not in adata.obs.columns:
        raise ValueError(f"在adata.obs中未找到'{patient_col}'列。可用的列有: {list(adata.obs.columns)}")
    
    if annotation_col not in adata.obs.columns:
        raise ValueError(f"在adata.obs中未找到'{annotation_col}'列。可用的列有: {list(adata.obs.columns)}")
    
    # 检查目标patient是否存在
    target_patients = [patient1, patient2]
    existing_patients = [patient for patient in target_patients if patient in adata.obs[patient_col].values]
    missing_patients = [patient for patient in target_patients if patient not in adata.obs[patient_col].values]
    
    print(f"\n=== Patient匹配情况 ===")
    print(f"目标patient: {target_patients}")
    print(f"找到的patient: {existing_patients}")
    print(f"未找到的patient: {missing_patients}")
    
    if len(existing_patients) < 2:
        print("错误: 没有找到足够的patient!")
        print("数据集中存在的patient示例:", list(adata.obs[patient_col].unique())[:10])
        return None
    
    # 提取目标patient数据
    adata_subset = adata[adata.obs[patient_col].isin(existing_patients), :].copy()
    
    print(f"\n=== 提取结果 ===")
    print(f"提取后的细胞数: {adata_subset.n_obs}")
    
    # 统计每个patient的细胞数量
    patient_counts = adata_subset.obs[patient_col].value_counts()
    print(f"\n=== 细胞数量统计 ===")
    for patient in existing_patients:
        count = patient_counts[patient]
        print(f"Patient {patient}: {count} 个细胞")
    
    # 计算细胞类型分布
    celltype_counts = pd.crosstab(adata_subset.obs[patient_col], adata_subset.obs[annotation_col])
    
    # 确保patient顺序：Q2236500-1在前，Z2315540-8在后
    celltype_counts = celltype_counts.reindex([patient1, patient2])
    
    # 转换为百分比
    celltype_percent = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100
    
    print(f"\n细胞类型分布(百分比):")
    print(celltype_percent)
    
    # 设置图形样式
    plt.style.use('default')
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.size'] = 12
    
    # 创建图形
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # 1. 绘制竖向堆积柱状图（百分比）
    n_celltypes = len(celltype_percent.columns)
    colors = plt.cm.Set3(np.linspace(0, 1, n_celltypes))
    
    # 创建堆积图
    bottom = np.zeros(len(celltype_percent))
    bar_width = 0.6
    
    for i, celltype in enumerate(celltype_percent.columns):
        values = celltype_percent[celltype].values
        bars = ax1.bar(range(len(celltype_percent)), values, bottom=bottom, 
                      width=bar_width, label=celltype, color=colors[i], 
                      alpha=0.8, edgecolor='white', linewidth=1)
        bottom += values
        
        # 在柱子上添加百分比标签（只显示大于5%的）
        for j, (val, bar) in enumerate(zip(values, bars)):
            if val > 5:  # 只在占比大于5%时显示标签
                height = bar.get_height()
                y_pos = bar.get_y() + height / 2
                ax1.text(bar.get_x() + bar.get_width()/2, y_pos, 
                        f'{val:.1f}%', ha='center', va='center', 
                        fontsize=9, fontweight='bold')
    
    # 设置图表属性
    ax1.set_title('细胞类型占比比较', fontsize=16, fontweight='bold', pad=20)
    ax1.set_ylabel('细胞类型占比 (%)', fontsize=14)
    ax1.set_xlabel('Patient', fontsize=14)
    ax1.set_xticks(range(len(celltype_percent)))
    ax1.set_xticklabels(celltype_percent.index, rotation=0)
    ax1.grid(axis='y', alpha=0.3, linestyle='--', zorder=0)
    ax1.set_ylim(0, 100)
    
    # 添加图例
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
              borderaxespad=0., fontsize=10, title='细胞类型', 
              title_fontsize=11)
    
    # 2. 绘制横向堆积图（更好地显示细胞类型标签）
    celltype_percent_transposed = celltype_percent.T
    
    left = np.zeros(len(celltype_percent_transposed.columns))
    for i, celltype in enumerate(celltype_percent_transposed.index):
        values = celltype_percent_transposed.loc[celltype].values
        bars = ax2.barh(range(len(celltype_percent_transposed.columns)), values, left=left, 
                       label=celltype, color=colors[i], alpha=0.8, edgecolor='white', linewidth=0.5)
        left += values
        
        # 添加百分比文本
        for j, (val, bar) in enumerate(zip(values, bars)):
            if val > 5:  # 只在占比大于5%时显示
                width = bar.get_width()
                x_pos = bar.get_x() + width / 2
                ax2.text(x_pos, bar.get_y() + bar.get_height()/2, f'{val:.1f}%', 
                        ha='center', va='center', fontsize=9, fontweight='bold')
    
    ax2.set_title('细胞类型占比横向堆积图', fontsize=16, fontweight='bold', pad=20)
    ax2.set_xlabel('百分比 (%)', fontsize=14)
    ax2.set_ylabel('Patient', fontsize=14)
    ax2.set_yticks(range(len(celltype_percent_transposed.columns)))
    ax2.set_yticklabels(celltype_percent_transposed.columns)
    ax2.grid(axis='x', alpha=0.3, linestyle='--', zorder=0)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    output_filename = f"celltype_comparison_{patient1}_vs_{patient2}.png"
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\n图表已保存为: {output_filename}")
    
    # 显示详细统计信息
    print(f"\n=== 详细统计信息 ===")
    for patient in existing_patients:
        patient_data = adata_subset[adata_subset.obs[patient_col] == patient]
        celltypes = patient_data.obs[annotation_col].value_counts()
        total_cells = len(patient_data)
        
        print(f"\nPatient {patient}:")
        print(f"总细胞数: {total_cells}")
        print("细胞类型分布:")
        for celltype, count in celltypes.items():
            percentage = (count / total_cells) * 100
            print(f"  {celltype}: {count} 个细胞 ({percentage:.2f}%)")
    
    plt.show()
    
    # 保存提取后的数据（包含raw矩阵）
    output_h5ad = "paired_1.h5ad"
    
    # 确保保留raw矩阵数据
    if adata_subset.raw is not None:
        print(f"保留raw矩阵数据，形状: {adata_subset.raw.X.shape}")
    else:
        print("警告: 未找到raw矩阵数据，将使用当前表达矩阵")
        # 将当前数据设置为raw
        adata_subset.raw = adata_subset
    
    # 保存h5ad文件
    adata_subset.write(output_h5ad)
    print(f"\n提取后的数据已保存为: {output_h5ad}")
    
    return adata_subset, celltype_percent

# 使用示例
if __name__ == "__main__":
    # 文件路径和patient名称
    h5ad_path = "adata_celltype_RCTD_removed_with_patient.h5ad"  # 请替换为您的h5ad文件路径
    patient1 = "Q2236500-1"  # 将显示在前（左侧）
    patient2 = "Z2315540-8"  # 将显示在后（右侧）
    
    try:
        # 执行处理流程
        result = process_paired_samples(h5ad_path, patient1, patient2)
        
        if result is not None:
            adata_subset, celltype_percent = result
            
            # 显示最终总结
            print(f"\n{'='*50}")
            print("处理完成总结:")
            print(f"提取的细胞总数: {adata_subset.n_obs}")
            print(f"包含的基因数: {adata_subset.n_vars}")
            print(f"唯一的细胞类型数: {adata_subset.obs['annotation'].nunique()}")
            print(f"输出文件: paired_1.h5ad")
            print(f"可视化文件: celltype_comparison_{patient1}_vs_{patient2}.png")
            
    except Exception as e:
        print(f"发生错误: {str(e)}")
        import traceback
        traceback.print_exc()
