import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

def create_output_folder(folder_name='plot'):
    """
    创建输出文件夹
    """
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"创建文件夹: {folder_name}")
    return folder_name

def load_data(adata_path, relate_path):
    """
    加载数据
    """
    print("加载数据...")
    
    # 加载h5ad数据
    adata = sc.read_h5ad(adata_path)
    print(f"h5ad数据加载完成: {adata.shape}")
    
    # 加载相关文件
    relate_df = pd.read_csv(relate_path,sep="\t")
    print(f"相关文件加载完成: {relate_df.shape}")
    print(f"相关文件前几行:\n{relate_df.head()}")
    
    return adata, relate_df

def calculate_celltype_proportions(adata, patient_col='patient', annotation_col='annotation'):
    """
    计算每个患者的细胞类型比例
    """
    print("\n计算细胞类型比例...")
    
    # 确保所需列存在
    if patient_col not in adata.obs.columns:
        raise ValueError(f"adata.obs中未找到'{patient_col}'列")
    if annotation_col not in adata.obs.columns:
        raise ValueError(f"adata.obs中未找到'{annotation_col}'列")
    
    # 计算每个患者的细胞总数
    patient_cell_counts = adata.obs[patient_col].value_counts()
    print(f"发现 {len(patient_cell_counts)} 个患者")
    
    # 计算每个患者中每种细胞类型的数量
    celltype_counts = pd.crosstab(adata.obs[patient_col], adata.obs[annotation_col])
    
    # 计算比例
    celltype_proportions = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100
    
    # 添加总计行
    celltype_proportions['total_cells'] = patient_cell_counts
    
    print(f"细胞类型比例矩阵形状: {celltype_proportions.shape}")
    print(f"细胞类型: {celltype_proportions.columns.tolist()}")
    print(f"前5个患者:\n{celltype_proportions.head()}")
    
    return celltype_proportions

def merge_and_calculate_correlation(celltype_proportions, relate_df, patient_col='patient'):
    """
    合并数据并计算相关性
    """
    print("\n合并数据并计算相关性...")
    
    # 检查celltype_proportions的索引
    print(f"celltype_proportions索引名称: {celltype_proportions.index.name}")
    print(f"celltype_proportions前5行索引:\n{celltype_proportions.index[:5]}")
    
    # 检查relate_df的列
    print(f"relate_df列: {relate_df.columns.tolist()}")
    print(f"relate_df前5行:\n{relate_df.head()}")
    
    # 重置celltype_proportions的索引，将患者从索引转换为列
    if celltype_proportions.index.name == patient_col or celltype_proportions.index.name is None:
        celltype_proportions_reset = celltype_proportions.reset_index()
        # 重命名索引列为患者列名
        if 'index' in celltype_proportions_reset.columns:
            celltype_proportions_reset = celltype_proportions_reset.rename(columns={'index': patient_col})
    else:
        celltype_proportions_reset = celltype_proportions.copy()
        # 如果患者列已经是列而不是索引，确保列名正确
        if patient_col not in celltype_proportions_reset.columns and celltype_proportions_reset.index.name is not None:
            celltype_proportions_reset[patient_col] = celltype_proportions_reset.index
    
    print(f"重置索引后的celltype_proportions列: {celltype_proportions_reset.columns.tolist()}")
    
    # 合并数据
    merged_df = pd.merge(celltype_proportions_reset, relate_df, on=patient_col, how='inner')
    
    print(f"合并后数据形状: {merged_df.shape}")
    print(f"合并后患者数: {merged_df[patient_col].nunique()}")
    print(f"合并后数据前5行:\n{merged_df.head()}")
    
    # 计算相关性
    correlation_results = []
    
    # 获取所有细胞类型列（排除患者列和prop列）
    celltype_cols = [col for col in merged_df.columns 
                     if col not in [patient_col, 'prop', 'total_cells']]
    
    print(f"将计算相关性的细胞类型列: {celltype_cols}")
    
    for celltype in celltype_cols:
        # 删除缺失值
        valid_data = merged_df[[celltype, 'prop']].dropna()
        
        if len(valid_data) >= 3:  # 至少需要3个点计算相关性
            # 计算Pearson相关系数
            corr_coef, p_value = stats.pearsonr(valid_data[celltype], valid_data['prop'])
            
            # 计算样本数
            n_samples = len(valid_data)
            
            correlation_results.append({
                'celltype': celltype,
                'correlation': corr_coef,
                'p_value': p_value,
                'n_samples': n_samples
            })
    
    # 转换为DataFrame
    correlation_df = pd.DataFrame(correlation_results)
    
    if len(correlation_df) > 0:
        # 按相关性绝对值排序
        correlation_df['abs_correlation'] = correlation_df['correlation'].abs()
        correlation_df = correlation_df.sort_values('abs_correlation', ascending=False)
    
    print(f"计算了 {len(correlation_df)} 种细胞类型的相关性")
    
    return merged_df, correlation_df

def create_correlation_plots(merged_df, correlation_df, save_dir='plot'):
    """
    创建相关性可视化图形
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 检查是否有相关性结果
    if len(correlation_df) == 0:
        print("警告: 没有计算到任何相关性，跳过绘图")
        return None, None, None
    
    # 获取细胞类型列
    celltype_cols = correlation_df['celltype'].tolist()
    
    # 1. 创建相关性柱状图
    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # 左图：相关性柱状图
    bar_colors = ['red' if x < 0 else 'blue' for x in correlation_df['correlation']]
    bars = ax1.bar(range(len(correlation_df)), correlation_df['correlation'], color=bar_colors)
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Correlation Coefficient')
    ax1.set_title('Correlation with Clinical Parameter')
    ax1.set_xticks(range(len(correlation_df)))
    ax1.set_xticklabels(correlation_df['celltype'], rotation=45, ha='right')
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    # 在柱状图上添加数值
    for i, (bar, corr, p_val) in enumerate(zip(bars, correlation_df['correlation'], correlation_df['p_value'])):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + (0.01 if height >= 0 else -0.02),
                f'{corr:.2f}', ha='center', va='bottom' if height >= 0 else 'top', fontsize=8)
        
        # 添加显著性标记
        if p_val < 0.001:
            sig_symbol = '***'
        elif p_val < 0.01:
            sig_symbol = '**'
        elif p_val < 0.05:
            sig_symbol = '*'
        else:
            sig_symbol = ''
        
        if sig_symbol:
            ax1.text(bar.get_x() + bar.get_width()/2., height + (0.03 if height >= 0 else -0.04),
                    sig_symbol, ha='center', va='bottom' if height >= 0 else 'top', fontsize=10, fontweight='bold')
    
    # 右图：p值热图
    p_value_matrix = correlation_df[['celltype', 'p_value']].copy()
    p_value_matrix['log10_p'] = -np.log10(p_value_matrix['p_value'])
    p_value_matrix = p_value_matrix.sort_values('log10_p', ascending=False)
    
    # 创建颜色映射
    colors = ['red' if p < 0.05 else 'gray' for p in correlation_df['p_value']]
    ax2.barh(range(len(p_value_matrix)), p_value_matrix['log10_p'], color=colors)
    ax2.set_xlabel('-log10(p-value)')
    ax2.set_ylabel('Cell Type')
    ax2.set_title('Statistical Significance (-log10 p-value)')
    ax2.set_yticks(range(len(p_value_matrix)))
    ax2.set_yticklabels(p_value_matrix['celltype'])
    ax2.axvline(x=-np.log10(0.05), color='red', linestyle='--', linewidth=1, alpha=0.7, label='p=0.05')
    ax2.axvline(x=-np.log10(0.01), color='orange', linestyle='--', linewidth=1, alpha=0.7, label='p=0.01')
    ax2.legend()
    
    plt.suptitle('Cell Type Proportion Correlation Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # 保存柱状图
    bar_plot_png = os.path.join(save_dir, 'correlation_bar_plot.png')
    bar_plot_pdf = os.path.join(save_dir, 'correlation_bar_plot.pdf')
    fig1.savefig(bar_plot_png, dpi=300, bbox_inches='tight')
    fig1.savefig(bar_plot_pdf, bbox_inches='tight')
    print(f"保存相关性柱状图: {bar_plot_png}")
    print(f"保存相关性柱状图: {bar_plot_pdf}")
    
    plt.show()
    plt.close(fig1)
    
    # 2. 创建散点图矩阵（显示前4个最相关的细胞类型）
    n_top = min(4, len(celltype_cols))
    top_celltypes = correlation_df['celltype'].head(n_top).tolist()
    
    fig2, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    for i, celltype in enumerate(top_celltypes):
        if i >= len(axes):
            break
            
        ax = axes[i]
        data = merged_df[[celltype, 'prop', 'patient']].dropna()
        
        if len(data) > 0:
            # 计算回归线
            x = data[celltype]
            y = data['prop']
            
            # 计算回归参数
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            line_x = np.array([x.min(), x.max()])
            line_y = intercept + slope * line_x
            
            # 绘制散点图
            scatter = ax.scatter(x, y, alpha=0.6, edgecolors='w', linewidth=0.5)
            ax.plot(line_x, line_y, 'r-', linewidth=2, 
                   label=f'r={r_value:.2f}, p={p_value:.3f}')
            
            # 添加患者标签
            for idx, row in data.iterrows():
                ax.annotate(row['patient'], (row[celltype], row['prop']), 
                           fontsize=6, alpha=0.7)
            
            ax.set_xlabel(f'{celltype} Proportion (%)')
            ax.set_ylabel('Clinical Parameter')
            ax.set_title(f'{celltype}\nr={r_value:.2f}, p={p_value:.3f}')
            ax.legend()
            ax.grid(True, alpha=0.3)
    
    # 隐藏多余的子图
    for i in range(len(top_celltypes), len(axes)):
        axes[i].axis('off')
    
    plt.suptitle('Top Correlated Cell Types - Scatter Plots', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # 保存散点图
    scatter_plot_png = os.path.join(save_dir, 'correlation_scatter_plot.png')
    scatter_plot_pdf = os.path.join(save_dir, 'correlation_scatter_plot.pdf')
    fig2.savefig(scatter_plot_png, dpi=300, bbox_inches='tight')
    fig2.savefig(scatter_plot_pdf, bbox_inches='tight')
    print(f"保存相关性散点图: {scatter_plot_png}")
    print(f"保存相关性散点图: {scatter_plot_pdf}")
    
    plt.show()
    plt.close(fig2)
    
    # 3. 创建热图
    fig3, ax3 = plt.subplots(figsize=(12, 8))
    
    # 准备热图数据
    heatmap_data = merged_df.copy()
    
    # 标准化数据
    from sklearn.preprocessing import StandardScaler
    
    # 选择要显示在热图中的列
    display_celltypes = correlation_df['celltype'].head(10).tolist()  # 只显示前10个
    heatmap_cols = display_celltypes + ['prop']
    
    if 'total_cells' in heatmap_data.columns:
        heatmap_cols = ['total_cells'] + heatmap_cols
    
    # 创建热图矩阵
    heatmap_matrix = heatmap_data[heatmap_cols].set_index(heatmap_data['patient'])
    
    # 标准化
    scaler = StandardScaler()
    heatmap_matrix_scaled = pd.DataFrame(
        scaler.fit_transform(heatmap_matrix),
        columns=heatmap_matrix.columns,
        index=heatmap_matrix.index
    )
    
    # 创建热图
    sns.heatmap(heatmap_matrix_scaled, cmap='RdBu_r', center=0, 
                ax=ax3, cbar_kws={'label': 'Z-score'})
    
    ax3.set_title('Patient-Level Cell Type Proportions and Clinical Parameter (Z-scored)', 
                  fontsize=14, fontweight='bold')
    ax3.set_xlabel('Variables')
    ax3.set_ylabel('Patients')
    
    plt.tight_layout()
    
    # 保存热图
    heatmap_png = os.path.join(save_dir, 'correlation_heatmap.png')
    heatmap_pdf = os.path.join(save_dir, 'correlation_heatmap.pdf')
    fig3.savefig(heatmap_png, dpi=300, bbox_inches='tight')
    fig3.savefig(heatmap_pdf, bbox_inches='tight')
    print(f"保存相关性热图: {heatmap_png}")
    print(f"保存相关性热图: {heatmap_pdf}")
    
    plt.show()
    plt.close(fig3)
    
    return fig1, fig2, fig3

def save_results(correlation_df, merged_df, save_dir='plot'):
    """
    保存结果到文件
    """
    save_dir = create_output_folder(save_dir)
    
    # 保存相关性结果
    correlation_file = os.path.join(save_dir, 'correlation_results.csv')
    correlation_df.to_csv(correlation_file, index=False)
    print(f"保存相关性结果: {correlation_file}")
    
    # 保存合并后的数据
    merged_file = os.path.join(save_dir, 'merged_data.csv')
    merged_df.to_csv(merged_file, index=False)
    print(f"保存合并数据: {merged_file}")
    
    # 保存统计摘要
    summary_file = os.path.join(save_dir, 'correlation_analysis_summary.txt')
    
    with open(summary_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("      Cell Type Proportion Correlation Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("1. DATA OVERVIEW\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total patients: {merged_df['patient'].nunique()}\n")
        f.write(f"Total cell types analyzed: {len(correlation_df)}\n")
        f.write(f"Clinical parameter column: 'prop'\n\n")
        
        f.write("2. CORRELATION RESULTS\n")
        f.write("-" * 40 + "\n")
        
        # 显著相关的细胞类型
        if len(correlation_df) > 0:
            significant_celltypes = correlation_df[correlation_df['p_value'] < 0.05]
            f.write(f"Significantly correlated cell types (p < 0.05): {len(significant_celltypes)}\n\n")
            
            f.write("Top 10 correlations:\n")
            f.write("-" * 30 + "\n")
            for i, row in correlation_df.head(10).iterrows():
                f.write(f"{row['celltype']}:\n")
                f.write(f"  Correlation: {row['correlation']:.3f}\n")
                f.write(f"  p-value: {row['p_value']:.4f}\n")
                f.write(f"  n_samples: {row['n_samples']}\n")
                f.write("\n")
            
            f.write("3. STATISTICAL SUMMARY\n")
            f.write("-" * 40 + "\n")
            f.write(f"Mean absolute correlation: {correlation_df['abs_correlation'].mean():.3f}\n")
            f.write(f"Median absolute correlation: {correlation_df['abs_correlation'].median():.3f}\n")
            f.write(f"Maximum correlation: {correlation_df['correlation'].max():.3f}\n")
            f.write(f"Minimum correlation: {correlation_df['correlation'].min():.3f}\n")
            f.write(f"Number with p < 0.05: {len(significant_celltypes)}\n")
            f.write(f"Number with p < 0.01: {len(correlation_df[correlation_df['p_value'] < 0.01])}\n")
            f.write(f"Number with p < 0.001: {len(correlation_df[correlation_df['p_value'] < 0.001])}\n")
        else:
            f.write("No correlation results available.\n")
        
        f.write("\n" + "=" * 60 + "\n")
        f.write("Analysis completed successfully!\n")
        f.write("=" * 60 + "\n")
    
    print(f"保存统计摘要: {summary_file}")
    
    return correlation_file, merged_file, summary_file

def main():
    """
    主函数
    """
    # 文件路径
    adata_path = 'adata_with_band_by_sample_improved.h5ad'
    relate_path = 'releate.csv'  # 假设您的相关文件是CSV格式
    
    # 检查文件是否存在
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"未找到h5ad文件: {adata_path}")
    if not os.path.exists(relate_path):
        raise FileNotFoundError(f"未找到相关文件: {relate_path}")
    
    # 1. 加载数据
    adata, relate_df = load_data(adata_path, relate_path)
    
    # 2. 计算细胞类型比例
    celltype_proportions = calculate_celltype_proportions(
        adata, 
        patient_col='patient', 
        annotation_col='annotation'
    )
    
    # 3. 合并数据并计算相关性
    merged_df, correlation_df = merge_and_calculate_correlation(
        celltype_proportions, 
        relate_df, 
        patient_col='patient'
    )
    
    # 4. 创建可视化图形
    print("\n创建可视化图形...")
    create_correlation_plots(merged_df, correlation_df, save_dir='plot')
    
    # 5. 保存结果
    print("\n保存结果...")
    save_results(correlation_df, merged_df, save_dir='plot')
    
    # 6. 显示关键结果
    print("\n" + "=" * 60)
    print("关键结果摘要")
    print("=" * 60)
    
    if len(correlation_df) > 0:
        # 显示相关性最强的5个细胞类型
        print("\n相关性最强的5个细胞类型:")
        print("-" * 40)
        for i, row in correlation_df.head(5).iterrows():
            sig_symbol = ""
            if row['p_value'] < 0.001:
                sig_symbol = "***"
            elif row['p_value'] < 0.01:
                sig_symbol = "**"
            elif row['p_value'] < 0.05:
                sig_symbol = "*"
            
            print(f"{row['celltype']}: r = {row['correlation']:.3f} {sig_symbol} (p = {row['p_value']:.4f})")
        
        # 显示相关性最弱的5个细胞类型
        print("\n相关性最弱的5个细胞类型:")
        print("-" * 40)
        for i, row in correlation_df.tail(5).iterrows():
            sig_symbol = ""
            if row['p_value'] < 0.001:
                sig_symbol = "***"
            elif row['p_value'] < 0.01:
                sig_symbol = "**"
            elif row['p_value'] < 0.05:
                sig_symbol = "*"
            
            print(f"{row['celltype']}: r = {row['correlation']:.3f} {sig_symbol} (p = {row['p_value']:.4f})")
    else:
        print("\n没有计算到相关性结果")
    
    print("\n分析完成！所有结果已保存到 'plot' 文件夹。")
    
    return adata, relate_df, merged_df, correlation_df

# 运行主函数
if __name__ == "__main__":
    adata, relate_df, merged_df, correlation_df = main()
