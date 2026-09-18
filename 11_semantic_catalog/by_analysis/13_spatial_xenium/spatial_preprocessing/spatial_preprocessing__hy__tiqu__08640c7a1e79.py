#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
h5ad数据处理脚本
根据患者列表提取样本并添加癌症类型标签
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys
import warnings
warnings.filterwarnings('ignore')

# 患者列表定义
short_patients = ["Z1802235-10","Z1736838-15","Q2324341-1"]
long_patients = ["Z1954693-9","Z1952876-9","Z1915172-7","Z1826868-6","F1805049-2","Z1814675-2","F1706967-11","Z1642635-14","Z1635989-12","Z1530765-11"]

def load_and_process_data(h5ad_path, output_path=None):
    """
    加载h5ad数据，根据患者列表提取样本并添加癌症类型标签
    
    参数:
    h5ad_path: 输入h5ad文件路径
    output_path: 输出h5ad文件路径，如果为None则不保存
    
    返回:
    处理后的AnnData对象
    """
    print(f"正在加载h5ad文件: {h5ad_path}")
    
    try:
        adata = sc.read_h5ad(h5ad_path)
    except Exception as e:
        print(f"加载文件失败: {e}")
        return None
    
    print(f"原始数据形状: {adata.shape}")
    print(f"细胞数量: {adata.n_obs}")
    print(f"基因数量: {adata.n_vars}")
    
    # 检查patient列是否存在
    if 'patient' not in adata.obs.columns:
        print("错误: obs中不存在'patient'列")
        print(f"可用的列: {adata.obs.columns.tolist()}")
        return None
    
    # 统计原始数据中的患者分布
    print(f"\n原始数据中患者数量: {adata.obs['patient'].nunique()}")
    
    # 创建患者列表集合
    mpc_set = set(short_patients)
    tnbc_set = set(long_patients)
    
    # 检查患者列表是否有重叠
    overlap = mpc_set.intersection(tnbc_set)
    if overlap:
        print(f"警告: 患者列表有重叠: {overlap}")
    
    # 筛选属于指定列表的患者
    selected_patients = mpc_set.union(tnbc_set)
    
    # 创建筛选掩码
    patient_mask = adata.obs['patient'].isin(selected_patients)
    
    print(f"\n目标患者总数: {len(selected_patients)}")
    print(f"  short患者: {len(mpc_set)}")
    print(f"  long患者: {len(tnbc_set)}")
    
    # 统计在数据中实际找到的患者
    found_patients = set(adata.obs['patient'][patient_mask].unique())
    found_mpc = found_patients.intersection(mpc_set)
    found_tnbc = found_patients.intersection(tnbc_set)
    
    print(f"\n在数据中找到的目标患者: {len(found_patients)}")
    print(f"  找到的short患者: {len(found_mpc)}")
    print(f"  找到的long患者: {len(found_tnbc)}")
    
    # 检查未找到的患者
    not_found = selected_patients - found_patients
    if not_found:
        print(f"  未找到的患者 ({len(not_found)}个): {sorted(list(not_found))[:10]}{'...' if len(not_found) > 10 else ''}")
    
    # 如果没有任何目标患者，报错
    if len(found_patients) == 0:
        print("错误: 在数据中没有找到任何目标患者")
        return None
    
    # 提取目标患者的数据
    adata_filtered = adata[patient_mask].copy()
    print(f"\n筛选后数据形状: {adata_filtered.shape}")
    print(f"筛选后细胞数量: {adata_filtered.n_obs}")
    
    # 添加OS列
    print("\n添加OS列...")
    
    # 创建新的OS列
    OS = []
    
    for patient in adata_filtered.obs['patient']:
        if patient in mpc_set:
            OS.append('short')
        elif patient in tnbc_set:
            OS.append('long')
        else:
            OS.append('Unknown')
    
    adata_filtered.obs['OS'] = OS
    
    # 统计癌症类型分布
    print("\n癌症类型分布:")
    cancer_counts = adata_filtered.obs['OS'].value_counts()
    for cancer, count in cancer_counts.items():
        percentage = count / adata_filtered.n_obs * 100
        print(f"  {cancer}: {count} 个细胞 ({percentage:.2f}%)")
    
    # 按患者统计
    print("\n按患者统计:")
    patient_stats = adata_filtered.obs.groupby(['patient', 'OS']).size().reset_index(name='cell_count')
    
    for cancer in ['short', 'long']:
        cancer_patients = patient_stats[patient_stats['OS'] == cancer]
        if not cancer_patients.empty:
            print(f"\n  {cancer}患者 ({len(cancer_patients)}个):")
            for _, row in cancer_patients.iterrows():
                print(f"    {row['patient']}: {row['cell_count']} 个细胞")
    
    # 保存处理后的数据
    if output_path:
        print(f"\n保存处理后的数据到: {output_path}")
        adata_filtered.write_h5ad(output_path)
        print("保存完成")
    
    return adata_filtered

def create_summary_statistics(adata, output_dir="output"):
    """
    创建汇总统计信息
    
    参数:
    adata: 处理后的AnnData对象
    output_dir: 输出目录
    """
    print(f"\n=== 创建汇总统计 ===")
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 癌症类型统计
    cancer_summary = adata.obs['OS'].value_counts().reset_index()
    cancer_summary.columns = ['OS', 'cell_count']
    cancer_summary['percentage'] = cancer_summary['cell_count'] / len(adata) * 100
    
    # 2. 患者统计
    patient_summary = adata.obs.groupby(['OS', 'patient']).size().reset_index(name='cell_count')
    
    # 3. 样本统计（如果有sample列）
    if 'sample' in adata.obs.columns:
        sample_summary = adata.obs.groupby(['OS', 'sample']).size().reset_index(name='cell_count')
        sample_summary_path = os.path.join(output_dir, "sample_summary.csv")
        sample_summary.to_csv(sample_summary_path, index=False, encoding='utf-8-sig')
        print(f"样本统计已保存: {sample_summary_path}")
    
    # 保存统计文件
    cancer_summary_path = os.path.join(output_dir, "OS_summary.csv")
    patient_summary_path = os.path.join(output_dir, "patient_summary.csv")
    
    cancer_summary.to_csv(cancer_summary_path, index=False, encoding='utf-8-sig')
    patient_summary.to_csv(patient_summary_path, index=False, encoding='utf-8-sig')
    
    print(f"癌症类型统计已保存: {cancer_summary_path}")
    print(f"患者统计已保存: {patient_summary_path}")
    
    # 打印汇总信息
    print(f"\n汇总统计:")
    print(f"  总细胞数: {len(adata)}")
    print(f"  癌症类型: {adata.obs['OS'].nunique()} 种")
    print(f"  患者数量: {adata.obs['patient'].nunique()}")
    
    if 'sample' in adata.obs.columns:
        print(f"  样本数量: {adata.obs['sample'].nunique()}")
    
    return cancer_summary, patient_summary

def visualize_data(adata, output_dir="output"):
    """
    创建数据可视化
    
    参数:
    adata: 处理后的AnnData对象
    output_dir: 输出目录
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        print(f"\n=== 创建可视化图表 ===")
        
        # 创建可视化目录
        viz_dir = os.path.join(output_dir, "visualizations")
        os.makedirs(viz_dir, exist_ok=True)
        
        # 1. 癌症类型细胞数量条形图
        plt.figure(figsize=(8, 6))
        
        cancer_counts = adata.obs['OS'].value_counts()
        bars = plt.bar(range(len(cancer_counts)), cancer_counts.values, 
                       color=['skyblue', 'lightcoral', 'lightgray'])
        
        plt.title('Cell Count by OS', fontsize=14, pad=20)
        plt.xlabel('OS', fontsize=12)
        plt.ylabel('Number of Cells', fontsize=12)
        plt.xticks(range(len(cancer_counts)), cancer_counts.index)
        
        # 在柱子上添加数值
        for i, bar in enumerate(bars):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{int(height)}', ha='center', va='bottom', fontsize=10)
        
        plt.tight_layout()
        bar_chart_path = os.path.join(viz_dir, "OS_cell_count.png")
        plt.savefig(bar_chart_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  癌症类型细胞数量图已保存: {bar_chart_path}")
        
        # 2. 患者细胞数量点图
        plt.figure(figsize=(12, 8))
        
        # 按癌症类型分组
        patient_stats = adata.obs.groupby(['OS', 'patient']).size().reset_index(name='cell_count')
        
        # 创建点图
        for i, OS in enumerate(patient_stats['OS'].unique()):
            cancer_data = patient_stats[patient_stats['OS'] == OS]
            plt.scatter([i] * len(cancer_data), cancer_data['cell_count'], 
                       alpha=0.7, s=100, label=OS)
            
            # 添加中位数线
            median_val = cancer_data['cell_count'].median()
            plt.hlines(median_val, i-0.2, i+0.2, colors='red', linewidth=2, 
                      label=f'{OS} median' if i == 0 else "")
        
        plt.title('Cell Count per Patient by OS', fontsize=14, pad=20)
        plt.xlabel('OS', fontsize=12)
        plt.ylabel('Number of Cells per Patient', fontsize=12)
        plt.xticks(range(len(patient_stats['OS'].unique())), 
                  patient_stats['OS'].unique())
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        patient_chart_path = os.path.join(viz_dir, "patient_cell_count.png")
        plt.savefig(patient_chart_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  患者细胞数量点图已保存: {patient_chart_path}")
        
        print(f"  所有可视化图表已保存到: {viz_dir}")
        
    except Exception as e:
        print(f"  创建可视化图表时出错: {e}")
        import traceback
        traceback.print_exc()

def main():
    """主函数"""
    # 文件路径配置
    h5ad_path = "CN_add.h5ad"  # 替换为您的h5ad文件路径
    output_h5ad_path = "filtered_with_OS.h5ad"  # 输出h5ad文件路径
    output_dir = "output"  # 输出目录
    
    print("="*60)
    print("h5ad数据处理工具 - 根据患者列表提取样本并添加癌症类型")
    print("="*60)
    
    try:
        # 1. 加载和处理数据
        adata_filtered = load_and_process_data(h5ad_path, output_h5ad_path)
        
        if adata_filtered is None:
            print("数据处理失败，程序退出")
            sys.exit(1)
        
        # 2. 创建汇总统计
        cancer_summary, patient_summary = create_summary_statistics(adata_filtered, output_dir)
        
        # 3. 创建可视化图表
        visualize_data(adata_filtered, output_dir)
        
        # 4. 打印最终摘要
        print(f"\n" + "="*60)
        print("处理完成!")
        print("="*60)
        print(f"输入文件: {h5ad_path}")
        print(f"输出h5ad文件: {output_h5ad_path}")
        print(f"输出目录: {output_dir}")
        
        # 显示文件列表
        print(f"\n生成的文件:")
        for root, dirs, files in os.walk(output_dir):
            for file in files:
                file_path = os.path.join(root, file)
                relative_path = os.path.relpath(file_path, output_dir)
                print(f"  {relative_path}")
    
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
