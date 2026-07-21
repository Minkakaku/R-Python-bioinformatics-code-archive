#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CN分布统计工具 - 生成Y/N格式表格
统计每个CN值中细胞类型在样本中的分布
"""

import pandas as pd
import numpy as np
import os
import sys
import warnings
from pathlib import Path
warnings.filterwarnings('ignore')

def load_data(filepath):
    """
    加载CSV数据文件
    
    参数:
    filepath: 文件路径
    
    返回:
    DataFrame对象
    """
    print(f"正在读取数据文件: {filepath}")
    
    # 尝试不同的分隔符
    df = pd.read_csv(filepath, sep='\t')
    
    # 检查必要的列
    required_cols = ['CNs_10', 'sample', 'annotation']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        print(f"错误: 文件中缺少必要的列: {missing_cols}")
        print(f"文件中可用的列: {df.columns.tolist()}")
        return None
    
    print(f"数据形状: {df.shape}")
    print(f"总细胞数: {len(df)}")
    
    return df

def clean_zwsp_in_dataframe(df):
    """
    清理DataFrame中的零宽空格字符
    
    参数:
    df: DataFrame对象
    
    返回:
    清理后的DataFrame
    """
    print("清理数据中的零宽空格字符...")
    
    # 零宽空格字符列表
    zwsp_chars = ['\u200b', '\u200c', '\u200d', '\u200e', '\u200f', '\ufeff']
    
    # 清理列名
    original_columns = df.columns.tolist()
    cleaned_columns = []
    
    for col in original_columns:
        if isinstance(col, str):
            cleaned = col
            for char in zwsp_chars:
                cleaned = cleaned.replace(char, '')
            cleaned_columns.append(cleaned)
        else:
            cleaned_columns.append(col)
    
    df.columns = cleaned_columns
    
    # 清理字符串列中的数据
    str_columns = df.select_dtypes(include=['object']).columns.tolist()
    cleaned_count = 0
    
    for col in str_columns:
        # 检查该列是否包含零宽空格
        if df[col].dtype == 'object':
            for char in zwsp_chars:
                if df[col].astype(str).str.contains(char, regex=False).any():
                    df[col] = df[col].astype(str).str.replace(char, '', regex=False)
                    cleaned_count += 1
    
    if cleaned_count > 0:
        print(f"清理了 {cleaned_count} 个列中的零宽空格字符")
    else:
        print("未发现零宽空格字符")
    
    return df

def generate_cn_distribution_tables(df, output_dir="cn_distribution_yn"):
    """
    生成每个CN的细胞类型-样本分布表（Y/N格式）
    
    参数:
    df: 包含数据的DataFrame
    output_dir: 输出目录
    
    返回:
    统计结果字典
    """
    print("\n" + "="*60)
    print("生成CN分布表 (Y/N格式)")
    print("="*60)
    
    # 创建输出目录
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # 获取所有唯一的CN值
    cn_values = df['CNs_10'].dropna().unique()
    cn_values = sorted(cn_values)
    
    print(f"发现 {len(cn_values)} 个不同的CN值:")
    for i, cn in enumerate(cn_values, 1):
        print(f"  {i}. {cn}")
    
    # 获取所有样本
    all_samples = sorted(df['sample'].dropna().unique())
    print(f"\n发现 {len(all_samples)} 个样本")
    
    # 获取所有细胞类型
    all_celltypes = sorted(df['annotation'].dropna().unique())
    print(f"发现 {len(all_celltypes)} 种细胞类型")
    
    # 为每个CN值生成分布表
    cn_stats = {}
    
    for cn_value in cn_values:
        print(f"\n=== 处理 CN: {cn_value} ===")
        
        # 筛选当前CN的数据
        cn_mask = df['CNs_10'] == cn_value
        cn_df = df[cn_mask].copy()
        
        if len(cn_df) == 0:
            print(f"  CN {cn_value} 没有数据，跳过")
            continue
        
        print(f"  CN {cn_value} 的细胞总数: {len(cn_df)}")
        
        # 创建细胞类型-样本分布矩阵
        # 初始化一个全为'N'的矩阵
        distribution_matrix = pd.DataFrame(
            'N', 
            index=all_celltypes, 
            columns=all_samples
        )
        
        # 获取当前CN中每个细胞类型在哪些样本中出现
        for celltype in all_celltypes:
            # 筛选当前细胞类型的数据
            celltype_mask = cn_df['annotation'] == celltype
            celltype_samples = cn_df.loc[celltype_mask, 'sample'].unique()
            
            # 在这些样本中标记为'Y'
            for sample in celltype_samples:
                if sample in distribution_matrix.columns:
                    distribution_matrix.at[celltype, sample] = 'Y'
        
        # 添加统计列
        distribution_matrix['出现样本数'] = (distribution_matrix == 'Y').sum(axis=1)
        distribution_matrix['总样本数'] = len(all_samples)
        distribution_matrix['出现比例'] = distribution_matrix['出现样本数'] / len(all_samples)
        
        # 按出现样本数排序
        distribution_matrix = distribution_matrix.sort_values(
            '出现样本数', 
            ascending=False
        )
        
        # 保存到CSV文件
        safe_cn_name = str(cn_value).replace('/', '_').replace(' ', '_')
        output_filename = f"CN_{safe_cn_name}_distribution_YN.csv"
        output_filepath = output_path / output_filename
        
        distribution_matrix.to_csv(output_filepath, encoding='utf-8-sig')
        
        print(f"  已保存: {output_filepath}")
        print(f"  表格形状: {distribution_matrix.shape}")
        
        # 显示统计信息
        celltypes_with_data = (distribution_matrix['出现样本数'] > 0).sum()
        print(f"  有数据的细胞类型数: {celltypes_with_data}")
        print(f"  无数据的细胞类型数: {len(all_celltypes) - celltypes_with_data}")
        
        # 显示前几行
        print("\n  表格预览 (前5行):")
        preview_rows = min(5, len(distribution_matrix))
        preview_cols = min(5, len(distribution_matrix.columns))
        preview = distribution_matrix.iloc[:preview_rows, :preview_cols]
        print(preview)
        
        # 记录统计信息
        cn_stats[cn_value] = {
            'total_cells': len(cn_df),
            'celltypes_with_data': celltypes_with_data,
            'total_celltypes': len(all_celltypes),
            'total_samples': len(all_samples),
            'output_file': output_filepath
        }
    
    return cn_stats, all_celltypes, all_samples

def create_summary_statistics(df, cn_stats, all_celltypes, all_samples, output_dir="cn_distribution_yn"):
    """
    创建汇总统计信息
    
    参数:
    df: 原始数据
    cn_stats: CN统计信息字典
    all_celltypes: 所有细胞类型列表
    all_samples: 所有样本列表
    output_dir: 输出目录
    """
    print(f"\n=== 创建汇总统计 ===")
    
    # 创建汇总目录
    summary_dir = Path(output_dir) / "summary"
    summary_dir.mkdir(exist_ok=True)
    
    # 1. CN总体统计
    cn_summary = df['CNs_10'].value_counts().reset_index()
    cn_summary.columns = ['CN', '细胞数量']
    cn_summary['百分比'] = cn_summary['细胞数量'] / len(df) * 100
    
    # 添加详细统计
    detailed_stats = []
    for cn_value, stats in cn_stats.items():
        detailed_stats.append({
            'CN': cn_value,
            '细胞数量': stats['total_cells'],
            '细胞类型数': stats['celltypes_with_data'],
            '总细胞类型数': stats['total_celltypes'],
            '样本数': stats['total_samples']
        })
    
    detailed_df = pd.DataFrame(detailed_stats)
    cn_summary = pd.merge(cn_summary, detailed_df, on='CN', how='left')
    
    cn_summary_path = summary_dir / "CN_overall_summary.csv"
    cn_summary.to_csv(cn_summary_path, encoding='utf-8-sig', index=False)
    print(f"CN总体统计已保存: {cn_summary_path}")
    
    # 2. 样本-CN交叉表
    sample_cn_crosstab = pd.crosstab(df['sample'], df['CNs_10'])
    sample_cn_crosstab_path = summary_dir / "sample_CN_crosstab.csv"
    sample_cn_crosstab.to_csv(sample_cn_crosstab_path, encoding='utf-8-sig')
    print(f"样本-CN交叉表已保存: {sample_cn_crosstab_path}")
    
    # 3. 细胞类型出现频率汇总
    celltype_summary = []
    for celltype in all_celltypes:
        # 统计该细胞类型出现在哪些CN中
        celltype_data = df[df['annotation'] == celltype]
        cn_counts = celltype_data['CNs_10'].value_counts()
        
        for cn_value, count in cn_counts.items():
            celltype_summary.append({
                'celltype': celltype,
                'CN': cn_value,
                '细胞数量': count
            })
    
    celltype_summary_df = pd.DataFrame(celltype_summary)
    celltype_summary_path = summary_dir / "celltype_CN_summary.csv"
    celltype_summary_df.to_csv(celltype_summary_path, encoding='utf-8-sig', index=False)
    print(f"细胞类型-CN统计已保存: {celltype_summary_path}")
    
    # 4. 创建热图格式的汇总表
    # 创建一个矩阵，显示每个CN中哪些细胞类型有数据
    heatmap_data = []
    for celltype in all_celltypes[:50]:  # 只显示前50种细胞类型，避免太大
        row = {'celltype': celltype}
        for cn_value in cn_stats.keys():
            # 检查该细胞类型在当前CN中是否有数据
            mask = (df['CNs_10'] == cn_value) & (df['annotation'] == celltype)
            has_data = df[mask].shape[0] > 0
            row[cn_value] = 'Y' if has_data else 'N'
        heatmap_data.append(row)
    
    if heatmap_data:
        heatmap_df = pd.DataFrame(heatmap_data)
        heatmap_df.set_index('celltype', inplace=True)
        heatmap_path = summary_dir / "celltype_CN_heatmap_summary.csv"
        heatmap_df.to_csv(heatmap_path, encoding='utf-8-sig')
        print(f"细胞类型-CN热图汇总已保存: {heatmap_path}")
    
    # 打印简要统计
    print(f"\n汇总统计:")
    print(f"  总细胞数: {len(df)}")
    print(f"  样本数量: {len(all_samples)}")
    print(f"  细胞类型数量: {len(all_celltypes)}")
    print(f"  CN值数量: {len(cn_stats)}")
    
    return cn_summary, sample_cn_crosstab, celltype_summary_df

def create_visualizations(df, cn_summary, output_dir="cn_distribution_yn"):
    """
    创建可视化图表
    
    参数:
    df: 原始数据
    cn_summary: CN统计信息
    output_dir: 输出目录
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        print(f"\n=== 创建可视化图表 ===")
        
        # 设置中文字体
        plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False
        
        # 创建可视化目录
        viz_dir = Path(output_dir) / "visualizations"
        viz_dir.mkdir(exist_ok=True)
        
        # 1. CN分布条形图
        plt.figure(figsize=(12, 8))
        
        # 按细胞数量排序
        cn_summary_sorted = cn_summary.sort_values('细胞数量', ascending=False)
        
        bars = plt.bar(range(len(cn_summary_sorted)), cn_summary_sorted['细胞数量'], 
                       color='skyblue', edgecolor='black')
        
        plt.title('CN分布统计', fontsize=16, pad=20)
        plt.xlabel('CN', fontsize=14)
        plt.ylabel('细胞数量', fontsize=14)
        plt.xticks(range(len(cn_summary_sorted)), cn_summary_sorted['CN'], rotation=45, ha='right')
        
        # 在每个柱子上添加数值
        for i, bar in enumerate(bars):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{int(height)}', ha='center', va='bottom', fontsize=9)
        
        plt.tight_layout()
        bar_chart_path = viz_dir / "CN_distribution_bar.png"
        plt.savefig(bar_chart_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  CN分布条形图已保存: {bar_chart_path}")
        
        # 2. 样本-CN热图
        if len(df['sample'].unique()) > 1 and len(df['CNs_10'].unique()) > 1:
            plt.figure(figsize=(12, 8))
            
            # 创建样本-CN交叉表
            sample_cn_crosstab = pd.crosstab(df['sample'], df['CNs_10'])
            
            # 绘制热图
            sns.heatmap(sample_cn_crosstab, annot=True, fmt='d', cmap='YlOrRd', 
                       cbar_kws={'label': '细胞数量'})
            
            plt.title('各样本中CN分布热图', fontsize=16, pad=20)
            plt.xlabel('CN', fontsize=14)
            plt.ylabel('样本', fontsize=14)
            plt.tight_layout()
            
            heatmap_path = viz_dir / "sample_CN_heatmap.png"
            plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  样本-CN热图已保存: {heatmap_path}")
        
        # 3. 前20种细胞类型的CN分布堆叠条形图
        top_celltypes = df['annotation'].value_counts().head(20).index
        celltype_cn_data = df[df['annotation'].isin(top_celltypes)]
        
        if len(celltype_cn_data) > 0:
            plt.figure(figsize=(14, 10))
            
            # 创建交叉表
            celltype_cn_crosstab = pd.crosstab(celltype_cn_data['annotation'], 
                                               celltype_cn_data['CNs_10'])
            
            # 绘制堆叠条形图
            celltype_cn_crosstab.plot(kind='bar', stacked=True, figsize=(14, 10), 
                                     colormap='tab20c')
            
            plt.title('前20种细胞类型的CN分布（堆叠条形图）', fontsize=16, pad=20)
            plt.xlabel('细胞类型', fontsize=14)
            plt.ylabel('细胞数量', fontsize=14)
            plt.xticks(rotation=45, ha='right')
            plt.legend(title='CN', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            
            stacked_bar_path = viz_dir / "celltype_CN_stacked_bar.png"
            plt.savefig(stacked_bar_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  细胞类型-CN堆叠条形图已保存: {stacked_bar_path}")
        
        print(f"  所有可视化图表已保存到: {viz_dir}")
        
    except Exception as e:
        print(f"  创建可视化图表时出错: {e}")
        import traceback
        traceback.print_exc()

def main():
    """主函数"""
    # 文件路径配置
    csv_path = "CNs_10_all_cell_features.csv"  # 替换为您的CSV文件路径
    output_dir = "cn_distribution_yn"  # 输出目录
    
    print("=" * 60)
    print("CN分布统计工具 (Y/N格式)")
    print("=" * 60)
    
    try:
        # 1. 加载数据
        df = load_data(csv_path)
        if df is None:
            print("无法加载数据，程序退出")
            sys.exit(1)
        
        # 2. 清理零宽空格字符
        df = clean_zwsp_in_dataframe(df)
        
        # 3. 生成CN分布表
        cn_stats, all_celltypes, all_samples = generate_cn_distribution_tables(df, output_dir)
        
        if not cn_stats:
            print("没有生成任何CN分布表")
            sys.exit(1)
        
        # 4. 创建汇总统计
        cn_summary, sample_cn_crosstab, celltype_summary = create_summary_statistics(
            df, cn_stats, all_celltypes, all_samples, output_dir
        )
        
        # 5. 创建可视化图表
        create_visualizations(df, cn_summary, output_dir)
        
        # 打印最终摘要
        print(f"\n" + "="*60)
        print("处理完成!")
        print("="*60)
        print(f"输入文件: {csv_path}")
        print(f"输出目录: {output_dir}")
        print(f"处理的CN数量: {len(cn_stats)}")
        print(f"生成的分布表: {len(cn_stats)} 个")
        print(f"样本数量: {len(all_samples)}")
        print(f"细胞类型数量: {len(all_celltypes)}")
        
        # 显示文件列表
        print(f"\n生成的文件:")
        output_path = Path(output_dir)
        for file_path in output_path.rglob("*.csv"):
            if file_path.is_file():
                relative_path = file_path.relative_to(output_path)
                print(f"  {relative_path}")
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
