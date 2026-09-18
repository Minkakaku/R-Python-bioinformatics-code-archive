#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

def analyze_mucus_zone_distribution(adata):
    """
    分析mucus_zone分布
    任务1: 统计每个sample/pateint中inner细胞的比例
    任务2: 统计每个sample/pateint中inner细胞annotation类型分布
    """
    
    print("="*60)
    print("开始分析mucus_zone分布")
    print("="*60)
    
    # 检查必要的列是否存在
    required_cols = ['sample', 'patient', 'mucus_zone', 'annotation']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        raise ValueError(f"缺少必要的列: {missing_cols}")
    
    # 创建输出目录
    output_dir = "mucus_zone_analysis_results"
    os.makedirs(output_dir, exist_ok=True)
    
    # 任务1: 统计每个sample/pateint中inner细胞的比例
    print("\n1. 统计每个sample/pateint中inner细胞的比例")
    
    # 1.1 sample级别的inner细胞比例
    print("\n   sample级别统计:")
    sample_stats = calculate_inner_proportion_by_group(adata, group_col='sample')
    print(f"    分析完成，共统计 {len(sample_stats)} 个sample")
    
    # 1.2 patient级别的inner细胞比例
    print("\n   patient级别统计:")
    patient_stats = calculate_inner_proportion_by_group(adata, group_col='patient')
    print(f"    分析完成，共统计 {len(patient_stats)} 个patient")
    
    # 任务2: 统计每个sample/pateint中inner细胞的annotation类型分布
    print("\n2. 统计每个sample/pateint中inner细胞的annotation类型分布")
    
    # 2.1 sample级别的inner细胞annotation分布
    print("\n   sample级别inner细胞annotation分布:")
    sample_inner_annotation_dist = calculate_annotation_distribution_inner(
        adata, group_col='sample'
    )
    
    # 2.2 patient级别的inner细胞annotation分布
    print("\n   patient级别inner细胞annotation分布:")
    patient_inner_annotation_dist = calculate_annotation_distribution_inner(
        adata, group_col='patient'
    )
    
    # 保存结果
    print("\n3. 保存分析结果")
    save_results(
        sample_stats, 
        patient_stats, 
        sample_inner_annotation_dist, 
        patient_inner_annotation_dist, 
        output_dir
    )
    
    # 可视化结果
    print("\n4. 生成可视化图表")
    create_visualizations(
        sample_stats, 
        patient_stats, 
        sample_inner_annotation_dist, 
        patient_inner_annotation_dist, 
        output_dir
    )
    
    print("\n" + "="*60)
    print("分析完成!")
    print("="*60)
    print(f"所有结果已保存到: {output_dir}/")
    
    return {
        'sample_stats': sample_stats,
        'patient_stats': patient_stats,
        'sample_annotation_dist': sample_inner_annotation_dist,
        'patient_annotation_dist': patient_inner_annotation_dist
    }

def calculate_inner_proportion_by_group(adata, group_col='sample'):
    """
    计算每个分组中inner细胞的比例
    
    参数:
    adata: AnnData对象
    group_col: 分组列名 ('sample' 或 'patient')
    """
    print(f"    按 {group_col} 分组统计...")
    
    # 获取分组列表
    groups = adata.obs[group_col].unique()
    print(f"    发现 {len(groups)} 个{group_col}")
    
    # 统计结果列表
    results = []
    
    for group in groups:
        # 筛选当前分组的细胞
        group_mask = adata.obs[group_col] == group
        group_cells = adata.obs[group_mask]
        
        # 统计总细胞数
        total_cells = len(group_cells)
        
        # 统计inner细胞数
        inner_cells = group_cells[group_cells['mucus_zone'] == 'inner']
        inner_count = len(inner_cells)
        
        # 计算比例
        if total_cells > 0:
            inner_proportion = inner_count / total_cells * 100
        else:
            inner_proportion = 0
        
        # 添加到结果
        results.append({
            group_col: group,
            'total_cells': total_cells,
            'inner_cells': inner_count,
            'inner_proportion(%)': inner_proportion
        })
    
    # 转换为DataFrame
    results_df = pd.DataFrame(results)
    
    # 按inner比例排序
    results_df = results_df.sort_values('inner_proportion(%)', ascending=False)
    
    return results_df

def calculate_annotation_distribution_inner(adata, group_col='sample'):
    """
    计算每个分组中inner细胞的annotation类型分布
    
    参数:
    adata: AnnData对象
    group_col: 分组列名 ('sample' 或 'patient')
    """
    print(f"    按 {group_col} 分组统计inner细胞的annotation分布...")
    
    # 筛选inner细胞
    inner_mask = adata.obs['mucus_zone'] == 'inner'
    inner_adata = adata[inner_mask]
    
    if len(inner_adata) == 0:
        print(f"    警告: 没有找到inner细胞")
        return pd.DataFrame()
    
    # 获取所有annotation类型
    all_annotations = inner_adata.obs['annotation'].unique()
    print(f"    发现 {len(all_annotations)} 种annotation类型")
    
    # 获取分组列表
    groups = inner_adata.obs[group_col].unique()
    
    # 创建交叉表
    crosstab = pd.crosstab(
        inner_adata.obs[group_col], 
        inner_adata.obs['annotation']
    )
    
    # 计算比例
    proportion_crosstab = crosstab.div(crosstab.sum(axis=1), axis=0) * 100
    
    return proportion_crosstab

def save_results(sample_stats, patient_stats, 
                 sample_annotation_dist, patient_annotation_dist, output_dir):
    """
    保存所有分析结果
    """
    # 保存sample级别的inner比例
    sample_stats_path = os.path.join(output_dir, "sample_inner_proportion.tsv")
    sample_stats.to_csv(sample_stats_path, sep='\t', index=False)
    print(f"    已保存: {sample_stats_path}")
    
    # 保存patient级别的inner比例
    patient_stats_path = os.path.join(output_dir, "patient_inner_proportion.tsv")
    patient_stats.to_csv(patient_stats_path, sep='\t', index=False)
    print(f"    已保存: {patient_stats_path}")
    
    # 保存sample级别的inner细胞annotation分布
    if not sample_annotation_dist.empty:
        sample_annotation_path = os.path.join(output_dir, "sample_inner_annotation_distribution.tsv")
        sample_annotation_dist.to_csv(sample_annotation_path, sep='\t')
        print(f"    已保存: {sample_annotation_path}")
    
    # 保存patient级别的inner细胞annotation分布
    if not patient_annotation_dist.empty:
        patient_annotation_path = os.path.join(output_dir, "patient_inner_annotation_distribution.tsv")
        patient_annotation_dist.to_csv(patient_annotation_path, sep='\t')
        print(f"    已保存: {patient_annotation_path}")
    
    # 创建汇总报告
    create_summary_report(sample_stats, patient_stats, 
                          sample_annotation_dist, patient_annotation_dist, output_dir)

def create_summary_report(sample_stats, patient_stats, 
                          sample_annotation_dist, patient_annotation_dist, output_dir):
    """
    创建汇总报告
    """
    report_path = os.path.join(output_dir, "analysis_summary.txt")
    
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("MUCUS_ZONE 分析结果汇总报告\n")
        f.write("="*80 + "\n\n")
        
        # 1. 样本概述
        f.write("1. 数据概述\n")
        f.write("-"*40 + "\n")
        f.write(f"样本总数: {len(sample_stats)}\n")
        f.write(f"病人总数: {len(patient_stats)}\n")
        
        # 2. inner细胞比例统计
        f.write("\n2. Inner细胞比例统计\n")
        f.write("-"*40 + "\n")
        
        # sample级别
        f.write("\nSample级别:\n")
        f.write(f"  平均inner比例: {sample_stats['inner_proportion(%)'].mean():.2f}%\n")
        f.write(f"  inner比例范围: {sample_stats['inner_proportion(%)'].min():.2f}% - {sample_stats['inner_proportion(%)'].max():.2f}%\n")
        f.write(f"  inner比例中位数: {sample_stats['inner_proportion(%)'].median():.2f}%\n")
        
        # 统计inner比例较高的样本
        high_inner_samples = sample_stats[sample_stats['inner_proportion(%)'] > 50]
        f.write(f"  inner比例>50%的样本数: {len(high_inner_samples)}\n")
        
        if len(high_inner_samples) > 0:
            f.write("  高inner比例样本 (前10个):\n")
            for i, (_, row) in enumerate(high_inner_samples.head(10).iterrows()):
                f.write(f"    {row['sample']}: {row['inner_proportion(%)']:.2f}% "
                       f"({row['inner_cells']}/{row['total_cells']})\n")
        
        # patient级别
        f.write("\nPatient级别:\n")
        f.write(f"  平均inner比例: {patient_stats['inner_proportion(%)'].mean():.2f}%\n")
        f.write(f"  inner比例范围: {patient_stats['inner_proportion(%)'].min():.2f}% - {patient_stats['inner_proportion(%)'].max():.2f}%\n")
        f.write(f"  inner比例中位数: {patient_stats['inner_proportion(%)'].median():.2f}%\n")
        
        # 3. annotation分布统计
        f.write("\n3. Inner细胞Annotation类型分布\n")
        f.write("-"*40 + "\n")
        
        if not sample_annotation_dist.empty:
            f.write("\nSample级别的主要annotation类型:\n")
            # 计算所有样本中平均比例最高的annotation类型
            mean_annotation_proportions = sample_annotation_dist.mean(axis=0)
            top_annotations = mean_annotation_proportions.nlargest(10)
            
            for i, (ann, prop) in enumerate(top_annotations.items()):
                f.write(f"  {i+1}. {ann}: {prop:.2f}%\n")
        
        if not patient_annotation_dist.empty:
            f.write("\nPatient级别的主要annotation类型:\n")
            mean_annotation_proportions_patient = patient_annotation_dist.mean(axis=0)
            top_annotations_patient = mean_annotation_proportions_patient.nlargest(10)
            
            for i, (ann, prop) in enumerate(top_annotations_patient.items()):
                f.write(f"  {i+1}. {ann}: {prop:.2f}%\n")
        
        # 4. 文件列表
        f.write("\n4. 输出文件\n")
        f.write("-"*40 + "\n")
        f.write("sample_inner_proportion.tsv - sample级别inner比例\n")
        f.write("patient_inner_proportion.tsv - patient级别inner比例\n")
        
        if not sample_annotation_dist.empty:
            f.write("sample_inner_annotation_distribution.tsv - sample级别annotation分布\n")
        
        if not patient_annotation_dist.empty:
            f.write("patient_inner_annotation_distribution.tsv - patient级别annotation分布\n")
        
        f.write(f"可视化图表保存到: {output_dir}/plots/\n")
    
    print(f"    已保存: {report_path}")

def create_visualizations(sample_stats, patient_stats, 
                          sample_annotation_dist, patient_annotation_dist, output_dir):
    """
    创建可视化图表
    """
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    # 设置图形风格
    plt.style.use('default')
    sns.set_palette("Set2")
    
    print("    生成可视化图表...")
    
    # 1. sample级别的inner比例条形图
    fig1, ax1 = plt.subplots(figsize=(14, 8))
    
    # 只显示前30个sample（如果太多）
    display_samples = sample_stats.head(30)
    
    bars = ax1.bar(range(len(display_samples)), 
                  display_samples['inner_proportion(%)'], 
                  color='steelblue', alpha=0.8)
    
    ax1.set_xlabel('Sample', fontsize=12)
    ax1.set_ylabel('Inner Cell Proportion (%)', fontsize=12)
    ax1.set_title('Inner Cell Proportion by Sample (Top 30)', fontsize=14, fontweight='bold')
    ax1.set_xticks(range(len(display_samples)))
    ax1.set_xticklabels(display_samples['sample'], rotation=45, ha='right', fontsize=8)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # 添加数值标签
    for i, bar in enumerate(bars):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "sample_inner_proportion.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(plots_dir, "sample_inner_proportion.pdf"), bbox_inches='tight')
    plt.close(fig1)
    print("      已生成: sample_inner_proportion.png/pdf")
    
    # 2. patient级别的inner比例条形图
    fig2, ax2 = plt.subplots(figsize=(12, 6))
    
    bars2 = ax2.bar(range(len(patient_stats)), 
                   patient_stats['inner_proportion(%)'], 
                   color='lightcoral', alpha=0.8)
    
    ax2.set_xlabel('Patient', fontsize=12)
    ax2.set_ylabel('Inner Cell Proportion (%)', fontsize=12)
    ax2.set_title('Inner Cell Proportion by Patient', fontsize=14, fontweight='bold')
    ax2.set_xticks(range(len(patient_stats)))
    ax2.set_xticklabels(patient_stats['patient'], rotation=45, ha='right')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # 添加数值标签
    for i, bar in enumerate(bars2):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "patient_inner_proportion.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(plots_dir, "patient_inner_proportion.pdf"), bbox_inches='tight')
    plt.close(fig2)
    print("      已生成: patient_inner_proportion.png/pdf")
    
    # 3. inner比例分布直方图
    fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(14, 6))
    
    # sample级别分布
    ax3a.hist(sample_stats['inner_proportion(%)'], bins=20, alpha=0.7, 
             color='steelblue', edgecolor='black')
    ax3a.axvline(sample_stats['inner_proportion(%)'].mean(), color='red', 
                linestyle='--', linewidth=2, label=f'Mean: {sample_stats["inner_proportion(%)"].mean():.2f}%')
    ax3a.set_xlabel('Inner Proportion (%)')
    ax3a.set_ylabel('Number of Samples')
    ax3a.set_title('Distribution of Inner Proportion (Sample Level)')
    ax3a.legend()
    ax3a.grid(True, alpha=0.3)
    
    # patient级别分布
    ax3b.hist(patient_stats['inner_proportion(%)'], bins=20, alpha=0.7, 
             color='lightcoral', edgecolor='black')
    ax3b.axvline(patient_stats['inner_proportion(%)'].mean(), color='red', 
                linestyle='--', linewidth=2, label=f'Mean: {patient_stats["inner_proportion(%)"].mean():.2f}%')
    ax3b.set_xlabel('Inner Proportion (%)')
    ax3b.set_ylabel('Number of Patients')
    ax3b.set_title('Distribution of Inner Proportion (Patient Level)')
    ax3b.legend()
    ax3b.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "inner_proportion_distribution.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(plots_dir, "inner_proportion_distribution.pdf"), bbox_inches='tight')
    plt.close(fig3)
    print("      已生成: inner_proportion_distribution.png/pdf")
    
    # 4. sample级别annotation热图（如果数据不是太多）
    if not sample_annotation_dist.empty and len(sample_annotation_dist) <= 50:
        # 只显示前20个最丰富的annotation类型
        top_annotations = sample_annotation_dist.mean(axis=0).nlargest(20).index
        heatmap_data = sample_annotation_dist[top_annotations]
        
        fig4, ax4 = plt.subplots(figsize=(16, 10))
        
        sns.heatmap(heatmap_data.T, cmap='YlOrRd', ax=ax4, 
                   cbar_kws={'label': 'Proportion (%)'}, fmt='.1f')
        
        ax4.set_xlabel('Sample', fontsize=12)
        ax4.set_ylabel('Annotation Type', fontsize=12)
        ax4.set_title('Annotation Distribution in Inner Cells (Sample Level)', 
                     fontsize=14, fontweight='bold')
        ax4.set_xticklabels(ax4.get_xticklabels(), rotation=45, ha='right')
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "sample_inner_annotation_heatmap.png"), 
                   dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(plots_dir, "sample_inner_annotation_heatmap.pdf"), 
                   bbox_inches='tight')
        plt.close(fig4)
        print("      已生成: sample_inner_annotation_heatmap.png/pdf")
    
    # 5. patient级别annotation热图
    if not patient_annotation_dist.empty:
        # 只显示前20个最丰富的annotation类型
        top_annotations_patient = patient_annotation_dist.mean(axis=0).nlargest(20).index
        heatmap_data_patient = patient_annotation_dist[top_annotations_patient]
        
        fig5, ax5 = plt.subplots(figsize=(12, 8))
        
        sns.heatmap(heatmap_data_patient.T, cmap='YlOrRd', ax=ax5, 
                   cbar_kws={'label': 'Proportion (%)'}, fmt='.1f')
        
        ax5.set_xlabel('Patient', fontsize=12)
        ax5.set_ylabel('Annotation Type', fontsize=12)
        ax5.set_title('Annotation Distribution in Inner Cells (Patient Level)', 
                     fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "patient_inner_annotation_heatmap.png"), 
                   dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(plots_dir, "patient_inner_annotation_heatmap.pdf"), 
                   bbox_inches='tight')
        plt.close(fig5)
        print("      已生成: patient_inner_annotation_heatmap.png/pdf")
    
    # 6. 最丰富的annotation类型条形图
    if not patient_annotation_dist.empty:
        # 计算patient级别的平均annotation分布
        mean_annotation_proportions = patient_annotation_dist.mean(axis=0)
        top_20_annotations = mean_annotation_proportions.nlargest(20)
        
        fig6, ax6 = plt.subplots(figsize=(12, 8))
        
        colors = plt.cm.Set2(np.arange(len(top_20_annotations)) / len(top_20_annotations))
        bars6 = ax6.barh(range(len(top_20_annotations)), 
                        top_20_annotations.values, 
                        color=colors, alpha=0.8)
        
        ax6.set_xlabel('Average Proportion in Inner Cells (%)', fontsize=12)
        ax6.set_ylabel('Annotation Type', fontsize=12)
        ax6.set_title('Top 20 Annotation Types in Inner Cells (Patient Level)', 
                     fontsize=14, fontweight='bold')
        ax6.set_yticks(range(len(top_20_annotations)))
        ax6.set_yticklabels(top_20_annotations.index)
        ax6.invert_yaxis()  # 最高的在顶部
        ax6.grid(True, alpha=0.3, axis='x')
        
        # 添加数值标签
        for i, (bar, value) in enumerate(zip(bars6, top_20_annotations.values)):
            width = bar.get_width()
            ax6.text(width + 0.5, bar.get_y() + bar.get_height()/2,
                    f'{value:.2f}%', va='center', fontsize=9)
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "top_annotation_types.png"), 
                   dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(plots_dir, "top_annotation_types.pdf"), 
                   bbox_inches='tight')
        plt.close(fig6)
        print("      已生成: top_annotation_types.png/pdf")

def main():
    """
    主函数
    """
    # 配置参数
    h5ad_file = "merged_TMA_data.h5ad"
    
    print("="*60)
    print("Mucus Zone 分析脚本")
    print("="*60)
    
    try:
        # 1. 读取数据
        print("正在读取h5ad文件...")
        adata = sc.read_h5ad(h5ad_file)
        print(f"数据形状: {adata.shape}")
        print(f"观察到的变量: {adata.obs.columns.tolist()}")
        
        # 检查必要的列
        required_cols = ['sample', 'patient', 'mucus_zone', 'annotation']
        for col in required_cols:
            if col in adata.obs.columns:
                print(f"  {col}: 存在")
            else:
                print(f"  {col}: 不存在 - 错误!")
        
        # 显示数据概览
        print(f"\n数据概览:")
        print(f"  总细胞数: {adata.n_obs}")
        print(f"  mucus_zone分布:")
        if 'mucus_zone' in adata.obs.columns:
            print(adata.obs['mucus_zone'].value_counts())
        print(f"  sample数量: {adata.obs['sample'].nunique()}")
        print(f"  patient数量: {adata.obs['patient'].nunique()}")
        if 'annotation' in adata.obs.columns:
            print(f"  annotation类型数量: {adata.obs['annotation'].nunique()}")
        
        # 2. 执行分析
        results = analyze_mucus_zone_distribution(adata)
        
        # 3. 显示关键结果摘要
        print("\n" + "="*60)
        print("关键结果摘要")
        print("="*60)
        
        sample_stats = results['sample_stats']
        patient_stats = results['patient_stats']
        
        print(f"\nSample级别统计:")
        print(f"  平均inner比例: {sample_stats['inner_proportion(%)'].mean():.2f}%")
        print(f"  最高inner比例: {sample_stats['inner_proportion(%)'].max():.2f}% "
              f"({sample_stats.loc[sample_stats['inner_proportion(%)'].idxmax(), 'sample']})")
        print(f"  最低inner比例: {sample_stats['inner_proportion(%)'].min():.2f}% "
              f"({sample_stats.loc[sample_stats['inner_proportion(%)'].idxmin(), 'sample']})")
        
        print(f"\nPatient级别统计:")
        print(f"  平均inner比例: {patient_stats['inner_proportion(%)'].mean():.2f}%")
        print(f"  最高inner比例: {patient_stats['inner_proportion(%)'].max():.2f}% "
              f"({patient_stats.loc[patient_stats['inner_proportion(%)'].idxmax(), 'patient']})")
        print(f"  最低inner比例: {patient_stats['inner_proportion(%)'].min():.2f}% "
              f"({patient_stats.loc[patient_stats['inner_proportion(%)'].idxmin(), 'patient']})")
        
        # 显示inner比例>50%的样本
        high_inner_samples = sample_stats[sample_stats['inner_proportion(%)'] > 50]
        if len(high_inner_samples) > 0:
            print(f"\nInner比例>50%的样本 ({len(high_inner_samples)}个):")
            for _, row in high_inner_samples.head(5).iterrows():
                print(f"  {row['sample']}: {row['inner_proportion(%)']:.2f}%")
            if len(high_inner_samples) > 5:
                print(f"  ... 还有{len(high_inner_samples)-5}个")
        
        print("\n分析完成! 请查看输出目录获取详细结果。")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
