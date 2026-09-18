#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scanpy as sc
import pandas as pd
import argparse
import os
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def add_mucus_zone_annotation(h5ad_file, tma_inner_file, output_file):
    """
    为h5ad文件添加mucus_zone注释
    
    参数:
    h5ad_file: 输入的h5ad文件路径
    tma_inner_file: TMA_inner文件路径
    output_file: 输出文件路径
    """
    
    print(f"1. 读取h5ad文件: {h5ad_file}")
    # 读取h5ad文件
    adata = sc.read_h5ad(h5ad_file)
    
    print(f"   obs形状: {adata.obs.shape}")
    print(f"   var形状: {adata.var.shape}")
    print(f"   样本数: {adata.n_obs}, 基因数: {adata.n_vars}")
    print(f"   obs列名: {list(adata.obs.columns)}")
    
    # 检查cell_id列是否存在
    if 'cell_id' not in adata.obs.columns:
        print("错误: h5ad的obs中没有cell_id列!")
        print(f"   obs中的列: {list(adata.obs.columns)}")
        return None
    
    print(f"\n   检查cell_id列:")
    print(f"   cell_id数据类型: {adata.obs['cell_id'].dtype}")
    print(f"   前5个cell_id值: {adata.obs['cell_id'].head().tolist()}")
    print(f"   cell_id是否有重复: {adata.obs['cell_id'].duplicated().any()}")
    if adata.obs['cell_id'].duplicated().any():
        print(f"   警告: cell_id列有重复值!")
        print(f"   重复值示例: {adata.obs['cell_id'][adata.obs['cell_id'].duplicated()].head().tolist()}")
    
    print(f"\n2. 读取TMA_inner文件: {tma_inner_file}")
    # 读取TMA_inner文件，无表头
    try:
        tma_inner = pd.read_csv(tma_inner_file, header=None, names=['cell_id'])
    except Exception as e:
        print(f"   读取TMA_inner文件失败: {e}")
        # 尝试用空格分隔
        try:
            tma_inner = pd.read_csv(tma_inner_file, header=None, sep='\s+', names=['cell_id'])
        except Exception as e2:
            print(f"   用空格分隔也失败: {e2}")
            return None
    
    print(f"   TMA_inner中的cell_id数量: {len(tma_inner)}")
    print(f"   前5个cell_id: {tma_inner['cell_id'].head().tolist()}")
    print(f"   TMA_inner中cell_id数据类型: {tma_inner['cell_id'].dtype}")
    
    # 确保TMA_inner中的cell_id是字符串类型
    tma_inner['cell_id'] = tma_inner['cell_id'].astype(str)
    
    # 将cell_id转换为集合
    inner_cell_ids = set(tma_inner['cell_id'].tolist())
    print(f"   转换为集合后数量: {len(inner_cell_ids)}")
    
    print(f"\n3. 确保adata.obs中的cell_id是字符串类型")
    adata.obs['cell_id'] = adata.obs['cell_id'].astype(str)
    
    print(f"\n4. 添加mucus_zone列")
    # 初始化mucus_zone列为'outer'
    adata.obs['mucus_zone'] = 'outer'
    
    # 统计inner细胞数量
    inner_mask = adata.obs['cell_id'].isin(inner_cell_ids)
    adata.obs.loc[inner_mask, 'mucus_zone'] = 'inner'
    
    inner_cells_found = inner_mask.sum()
    inner_cells_not_found = len(inner_cell_ids) - inner_cells_found
    
    print(f"   找到的inner细胞数: {inner_cells_found}")
    print(f"   未找到的inner细胞数: {inner_cells_not_found}")
    
    if inner_cells_not_found > 0:
        print(f"   警告: {inner_cells_not_found}个TMA_inner中的cell_id在adata中未找到")
        # 显示前5个未找到的cell_id
        not_found_ids = inner_cell_ids - set(adata.obs['cell_id'])
        not_found_examples = list(not_found_ids)[:5]
        print(f"   示例未找到的cell_id: {not_found_examples}")
    
    # 统计mucus_zone分布
    zone_counts = adata.obs['mucus_zone'].value_counts()
    print(f"\n   mucus_zone分布:")
    for zone, count in zone_counts.items():
        print(f"     {zone}: {count}个细胞 ({count/len(adata.obs)*100:.1f}%)")
    
    # 检查cell_id列和mucus_zone列的前几个值
    print(f"\n   检查前5个细胞的标注结果:")
    sample_df = adata.obs[['cell_id', 'mucus_zone']].head()
    for idx, row in sample_df.iterrows():
        print(f"     {row['cell_id']}: {row['mucus_zone']}")
    
    # 检查是否有patient列，如果有，显示patient与mucus_zone的交叉表
    if 'patient' in adata.obs.columns:
        print(f"\n5. Patient与mucus_zone分布:")
        cross_tab = pd.crosstab(adata.obs['patient'], adata.obs['mucus_zone'])
        print(cross_tab)
        
        # 计算每个patient中inner的比例
        if 'inner' in cross_tab.columns:
            cross_tab['inner_percent'] = cross_tab['inner'] / (cross_tab['inner'] + cross_tab.get('outer', 0)) * 100
            print(f"\n   每个patient中inner细胞比例:")
            print(cross_tab[['inner', 'inner_percent']].round(2))
    
    print(f"\n6. 保存到新文件: {output_file}")
    # 保存到新文件
    adata.write_h5ad(output_file)
    
    print(f"   文件保存成功!")
    
    # 验证保存的文件
    print(f"\n7. 验证保存的文件:")
    adata_verify = sc.read_h5ad(output_file)
    print(f"   读取验证: 成功")
    print(f"   mucus_zone列存在: {'mucus_zone' in adata_verify.obs.columns}")
    if 'mucus_zone' in adata_verify.obs.columns:
        verify_counts = adata_verify.obs['mucus_zone'].value_counts()
        print(f"   验证mucus_zone分布:")
        for zone, count in verify_counts.items():
            print(f"     {zone}: {count}个细胞")
    
    return adata

def visualize_results(adata, output_dir="results"):
    """
    可视化mucus_zone分布
    
    参数:
    adata: AnnData对象
    output_dir: 输出目录
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        os.makedirs(output_dir, exist_ok=True)
        
        print(f"\n8. 生成可视化图表:")
        
        # 1. mucus_zone分布饼图
        fig1, ax1 = plt.subplots(figsize=(8, 6))
        zone_counts = adata.obs['mucus_zone'].value_counts()
        colors = {'inner': 'lightcoral', 'outer': 'steelblue'}
        
        ax1.pie(zone_counts.values, labels=zone_counts.index, 
               autopct='%1.1f%%', colors=[colors.get(z, 'gray') for z in zone_counts.index])
        ax1.set_title('Distribution of mucus_zone', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/mucus_zone_distribution.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_dir}/mucus_zone_distribution.pdf", dpi=300, bbox_inches='tight')
        plt.close(fig1)
        print(f"   图表1保存: {output_dir}/mucus_zone_distribution.png/pdf")
        
        # 2. 如果存在UMAP坐标，可视化mucus_zone
        if 'X_umap' in adata.obsm:
            fig2, ax2 = plt.subplots(figsize=(10, 8))
            
            # 创建颜色映射
            scatter = ax2.scatter(
                adata.obsm['X_umap'][:, 0],
                adata.obsm['X_umap'][:, 1],
                c=[colors.get(z, 'gray') for z in adata.obs['mucus_zone']],
                s=10, alpha=0.7, edgecolors='none'
            )
            
            ax2.set_xlabel('UMAP1', fontsize=12)
            ax2.set_ylabel('UMAP2', fontsize=12)
            ax2.set_title('mucus_zone on UMAP', fontsize=14, fontweight='bold')
            
            # 创建图例
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor=colors['inner'], label='inner'),
                              Patch(facecolor=colors['outer'], label='outer')]
            ax2.legend(handles=legend_elements, loc='best')
            
            plt.tight_layout()
            plt.savefig(f"{output_dir}/mucus_zone_umap.png", dpi=300, bbox_inches='tight')
            plt.savefig(f"{output_dir}/mucus_zone_umap.pdf", dpi=300, bbox_inches='tight')
            plt.close(fig2)
            print(f"   图表2保存: {output_dir}/mucus_zone_umap.png/pdf")
        
        # 3. 如果存在patient列，创建分组可视化
        if 'patient' in adata.obs.columns:
            # 统计每个patient的mucus_zone分布
            patient_zone = pd.crosstab(adata.obs['patient'], adata.obs['mucus_zone'])
            
            fig3, ax3 = plt.subplots(figsize=(12, 8))
            
            if 'inner' in patient_zone.columns and 'outer' in patient_zone.columns:
                x = range(len(patient_zone))
                width = 0.35
                
                ax3.bar([i - width/2 for i in x], patient_zone['inner'], width, 
                       label='inner', color='lightcoral', alpha=0.8)
                ax3.bar([i + width/2 for i in x], patient_zone['outer'], width, 
                       label='outer', color='steelblue', alpha=0.8)
                
                ax3.set_xlabel('Patient', fontsize=12)
                ax3.set_ylabel('Cell Count', fontsize=12)
                ax3.set_title('mucus_zone Distribution by Patient', fontsize=14, fontweight='bold')
                ax3.set_xticks(x)
                ax3.set_xticklabels(patient_zone.index, rotation=45, ha='right')
                ax3.legend()
                ax3.grid(True, alpha=0.3, axis='y')
                
                plt.tight_layout()
                plt.savefig(f"{output_dir}/mucus_zone_by_patient.png", dpi=300, bbox_inches='tight')
                plt.savefig(f"{output_dir}/mucus_zone_by_patient.pdf", dpi=300, bbox_inches='tight')
                plt.close(fig3)
                print(f"   图表3保存: {output_dir}/mucus_zone_by_patient.png/pdf")
        
        print(f"   所有图表已保存到: {output_dir}/")
        
    except ImportError as e:
        print(f"   注意: matplotlib或seaborn未安装，跳过可视化: {e}")
    except Exception as e:
        print(f"   可视化过程中出错: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='为h5ad文件添加mucus_zone注释')
    parser.add_argument('h5ad_file', help='输入的h5ad文件路径')
    parser.add_argument('tma_inner_file', help='TMA_inner文件路径（只有一列cell_id，无表头）')
    parser.add_argument('--output_file', default='adata_with_CNs10_TMA_with_patient_MZ.h5ad',
                       help='输出文件路径，默认为adata_with_CNs10_TMA_with_patient_MZ.h5ad')
    parser.add_argument('--visualize', action='store_true', 
                       help='生成可视化图表')
    parser.add_argument('--output_dir', default='mucus_zone_results',
                       help='可视化输出目录，默认为mucus_zone_results')
    
    args = parser.parse_args()
    
    print("="*60)
    print("添加mucus_zone注释脚本")
    print("="*60)
    
    # 检查输入文件是否存在
    if not os.path.exists(args.h5ad_file):
        print(f"错误: h5ad文件不存在: {args.h5ad_file}")
        exit(1)
    
    if not os.path.exists(args.tma_inner_file):
        print(f"错误: TMA_inner文件不存在: {args.tma_inner_file}")
        exit(1)
    
    # 运行主函数
    adata = add_mucus_zone_annotation(args.h5ad_file, args.tma_inner_file, args.output_file)
    
    if adata is not None:
        print(f"\n" + "="*60)
        print("处理完成!")
        print("="*60)
        print(f"输入文件: {args.h5ad_file}")
        print(f"TMA_inner文件: {args.tma_inner_file}")
        print(f"输出文件: {args.output_file}")
        print(f"总细胞数: {adata.n_obs}")
        
        if 'mucus_zone' in adata.obs.columns:
            zone_counts = adata.obs['mucus_zone'].value_counts()
            print(f"mucus_zone分布:")
            for zone, count in zone_counts.items():
                print(f"  {zone}: {count}个细胞 ({count/adata.n_obs*100:.1f}%)")
        
        # 如果需要可视化
        if args.visualize:
            visualize_results(adata, args.output_dir)
    else:
        print("\n处理失败!")
