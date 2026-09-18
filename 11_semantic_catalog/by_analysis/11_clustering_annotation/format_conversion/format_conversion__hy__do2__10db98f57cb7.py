import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体和绘图样式
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

def load_and_process_data(file_path):
    """加载并处理单细胞数据"""
    print("正在加载数据...")
    adata = sc.read_h5ad(file_path)
    print(f"数据加载成功！细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
    
    # 检查必要的列是否存在
    required_cols = ['sample', 'celltype', 'TIME']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        print(f"警告: 缺少以下列: {missing_cols}")
        print(f"可用列: {list(adata.obs.columns)}")
        # 尝试使用可能的替代列名
        if 'sample' in adata.obs.columns and 'sample' not in adata.obs.columns:
            adata.obs['sample'] = adata.obs['sample']
        if 'cell_type' in adata.obs.columns and 'celltype' not in adata.obs.columns:
            adata.obs['celltype'] = adata.obs['cell_type']
        if 'treatment' in adata.obs.columns and 'TIME' not in adata.obs.columns:
            adata.obs['TIME'] = adata.obs['treatment']
    
    return adata

def calculate_celltype_proportions(adata, sample_col='sample', 
                                 celltype_col='celltype', treatment_col='TIME'):
    """计算每种细胞类型在每个样本中的比例"""
    print("正在计算细胞类型比例...")
    
    # 计算每个样本中每种细胞类型的细胞数量
    celltype_counts = pd.crosstab(adata.obs[sample_col], adata.obs[celltype_col])
    
    # 计算每个样本的细胞总数
    total_cells_per_sample = celltype_counts.sum(axis=1)
    
    # 计算每种细胞类型在每个样本中的比例
    celltype_proportions = celltype_counts.div(total_cells_per_sample, axis=0)
    
    # 重置索引并转换为长格式
    celltype_proportions  = celltype_proportions.reset_index()
    celltype_proportions.rename(columns={'index': sample_col}, inplace=True)
    
    # 转换为长格式
    proportion_long = pd.melt(
        celltype_proportions, 
        id_vars=[sample_col], 
        value_name='Proportion', 
        var_name='CellType'
    )
    
    # 添加治疗组信息
    sample_meta = adata.obs[[sample_col, treatment_col]].drop_duplicates()
    proportion_long = proportion_long.merge(sample_meta, on=sample_col)
    
    print(f"计算完成！共 {len(proportion_long)} 个样本-细胞类型组合")
    print(f"细胞类型: {proportion_long['CellType'].unique().tolist()}")
    print(f"治疗组: {proportion_long[treatment_col].unique().tolist()}")
    
    return proportion_long

def plot_individual_celltype_boxplots(proportion_data, treatment_col='TIME', 
                                    output_dir='group_b', figsize=(10, 6)):
    """为每个细胞类型绘制单独的分组比较箱线图[1,5](@ref)"""
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    celltypes = proportion_data['CellType'].unique()
    print(f"开始为 {len(celltypes)} 种细胞类型绘制单独的箱线图...")
    
    # 为每个细胞类型创建单独的箱线图[5](@ref)
    for celltype in celltypes:
        print(f"处理细胞类型: {celltype}")
        
        # 创建图形
        plt.figure(figsize=figsize)
        
        # 筛选当前细胞类型的数据
        celltype_data = proportion_data[proportion_data['CellType'] == celltype]
        
        # 绘制箱线图[1,5](@ref)
        ax = sns.boxplot(data=celltype_data, x=treatment_col, y='Proportion', 
                        palette="Set2", showfliers=True, width=0.6)
        
        # 添加数据点（散点）[7](@ref)
        sns.stripplot(data=celltype_data, x=treatment_col, y='Proportion', 
                     color='black', size=4, jitter=0.2, alpha=0.6)
        
        # 设置标题和标签
        plt.title(f'{celltype} - Cell Proportion by Treatment Group', 
                 fontsize=14, fontweight='bold', pad=15)
        plt.xlabel('Treatment Group', fontsize=12)
        plt.ylabel('Cell Proportion', fontsize=12)
        
        # 美化图形
        plt.xticks(rotation=45, ha='right')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # 保存图形[1](@ref)
        filename = f"{celltype.replace(' ', '_').replace('/', '_')}_boxplot.png"
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.close()
        
        print(f"已保存: {filepath}")
    
    print(f"所有细胞类型的箱线图已保存到 {output_dir} 文件夹")

def create_summary_report(proportion_data, treatment_col='TIME', output_dir='group_b'):
    """创建细胞类型比例的汇总报告"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 计算汇总统计
    summary_stats = proportion_data.groupby(['CellType', treatment_col])['Proportion'].agg([
        'count', 'mean', 'std', 'min', 'max', 'median'
    ]).round(4)
    
    # 保存汇总统计
    summary_stats.to_csv(os.path.join(output_dir, 'celltype_proportions_summary.csv'))
    
    # 创建汇总可视化
    plt.figure(figsize=(12, 8))
    
    # 按中位数排序细胞类型
    median_order = proportion_data.groupby('CellType')['Proportion'].median().sort_values(ascending=False).index
    
    sns.boxplot(data=proportion_data, x='CellType', y='Proportion', hue=treatment_col,
                palette="Set2", order=median_order)
    plt.xticks(rotation=45, ha='right')
    plt.title('All Cell Types - Proportion Summary', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'all_celltypes_summary.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    return summary_stats

def main():
    """主函数"""
    # 文件路径
    file_path = 'extracted_specific_patients.h5ad'
    
    try:
        # 1. 加载数据
        adata = load_and_process_data(file_path)
        
        # 2. 计算细胞类型比例
        proportion_data = calculate_celltype_proportions(adata)
        
        if proportion_data is None or len(proportion_data) == 0:
            print("错误: 无法计算细胞类型比例")
            return
        
        # 3. 为每个细胞类型绘制单独的箱线图（保存到group_b文件夹）
        print("\n为每个细胞类型绘制单独的箱线图...")
        plot_individual_celltype_boxplots(proportion_data, output_dir='group_b')
        
        # 4. 创建汇总报告
        print("创建汇总报告...")
        summary_stats = create_summary_report(proportion_data, output_dir='group_b')
        
        # 5. 保存原始数据
        proportion_data.to_csv(os.path.join('group_b', 'celltype_proportions_detailed.csv'), index=False)
        
        print("\n" + "="*50)
        print("分析完成！")
        print("生成的文件（保存在 group_b 文件夹）:")
        print("- 每个细胞类型的单独箱线图 PNG 文件")
        print("- celltype_proportions_summary.csv: 汇总统计")
        print("- all_celltypes_summary.png: 所有细胞类型汇总图")
        print("- celltype_proportions_detailed.csv: 详细比例数据")
        print("="*50)
        
        # 显示汇总统计
        print("\n细胞类型比例摘要:")
        print(summary_stats)
        
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
