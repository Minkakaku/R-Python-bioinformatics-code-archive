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
    
    return celltype_proportions

def merge_and_calculate_correlation(celltype_proportions, relate_df, patient_col='patient'):
    """
    合并数据并计算相关性
    """
    print("\n合并数据并计算相关性...")
    
    # 重置celltype_proportions的索引，将患者从索引转换为列
    celltype_proportions_reset = celltype_proportions.reset_index()
    # 重命名索引列为患者列名
    if 'index' in celltype_proportions_reset.columns:
        celltype_proportions_reset = celltype_proportions_reset.rename(columns={'index': patient_col})
    
    print(f"重置索引后的celltype_proportions列: {celltype_proportions_reset.columns.tolist()}")
    
    # 合并数据
    merged_df = pd.merge(celltype_proportions_reset, relate_df, on=patient_col, how='inner')
    
    print(f"合并后数据形状: {merged_df.shape}")
    print(f"合并后患者数: {merged_df[patient_col].nunique()}")
    
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
            try:
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
            except Exception as e:
                print(f"  警告: 计算细胞类型 {celltype} 的相关性时出错: {e}")
                print(f"  细胞类型 {celltype} 的数据范围: {valid_data[celltype].min():.2f} 到 {valid_data[celltype].max():.2f}")
    
    # 转换为DataFrame
    correlation_df = pd.DataFrame(correlation_results)
    
    if len(correlation_df) > 0:
        # 按相关性绝对值排序
        correlation_df['abs_correlation'] = correlation_df['correlation'].abs()
        correlation_df = correlation_df.sort_values('abs_correlation', ascending=False)
    
    print(f"计算了 {len(correlation_df)} 种细胞类型的相关性")
    
    return merged_df, correlation_df

def create_correlation_summary_plots(merged_df, correlation_df, save_dir='plot'):
    """
    创建相关性摘要图（柱状图、p值图、热图）
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 检查是否有相关性结果
    if len(correlation_df) == 0:
        print("警告: 没有计算到任何相关性，跳过摘要绘图")
        return None, None
    
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
    
    # 2. 创建热图
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
    
    return fig1, fig3

def create_individual_scatter_plots(merged_df, correlation_df, save_dir='plot', max_per_page=9):
    """
    为每个细胞类型创建单独的散点图（拟合曲线）
    每页显示多个子图
    
    参数:
    ----------
    max_per_page : int, default=9
        每页显示的最大子图数量
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 检查是否有相关性结果
    if len(correlation_df) == 0:
        print("警告: 没有计算到任何相关性，跳过散点图绘制")
        return
    
    # 获取所有细胞类型
    celltypes = correlation_df['celltype'].tolist()
    n_celltypes = len(celltypes)
    
    print(f"将为 {n_celltypes} 个细胞类型创建散点图")
    
    # 计算需要多少页
    n_pages = (n_celltypes + max_per_page - 1) // max_per_page
    
    # 为每个页面创建图形
    for page in range(n_pages):
        # 计算当前页的细胞类型索引范围
        start_idx = page * max_per_page
        end_idx = min((page + 1) * max_per_page, n_celltypes)
        page_celltypes = celltypes[start_idx:end_idx]
        
        # 计算子图的行列数
        n_subplots = len(page_celltypes)
        n_cols = min(3, n_subplots)
        n_rows = (n_subplots + n_cols - 1) // n_cols
        
        # 创建图形
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
        
        # 如果只有一个子图，确保axes是数组
        if n_subplots == 1:
            axes = np.array([axes])
        
        # 展平axes数组以便遍历
        axes_flat = axes.flatten() if n_rows * n_cols > 1 else axes
        
        # 为当前页的每个细胞类型创建子图
        for idx, (celltype, ax) in enumerate(zip(page_celltypes, axes_flat[:n_subplots])):
            data = merged_df[[celltype, 'prop', 'patient']].dropna()
            
            if len(data) > 0:
                # 获取细胞类型相关统计
                celltype_stats = correlation_df[correlation_df['celltype'] == celltype].iloc[0]
                corr_coef = celltype_stats['correlation']
                p_value = celltype_stats['p_value']
                
                # 提取x和y值
                x = data[celltype]
                y = data['prop']
                
                # 检查x值是否有变化
                if np.std(x) == 0 or len(np.unique(x)) == 1:
                    # x值没有变化，无法计算回归线
                    print(f"  警告: 细胞类型 {celltype} 的x值完全相同，跳过回归线绘制")
                    
                    # 绘制散点图（不添加样本标签）
                    ax.scatter(x, y, alpha=0.6, edgecolors='w', linewidth=0.5, s=50)
                    
                    # 设置标题
                    sig_symbol = ""
                    if p_value < 0.001:
                        sig_symbol = "***"
                    elif p_value < 0.01:
                        sig_symbol = "**"
                    elif p_value < 0.05:
                        sig_symbol = "*"
                    
                    title = f'{celltype}\nr={corr_coef:.3f}{sig_symbol} (p={p_value:.3g})'
                    ax.set_title(title, fontsize=10)
                    
                    # 设置坐标轴标签
                    ax.set_xlabel(f'{celltype} Proportion (%)')
                    ax.set_ylabel('Clinical Parameter')
                    
                    # 添加网格
                    ax.grid(True, alpha=0.3, linestyle='--')
                    
                    # 添加文本标注样本数
                    ax.text(0.05, 0.95, f'n={len(data)}\n(No variance in x)', 
                           transform=ax.transAxes, fontsize=9, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                else:
                    # 可以计算回归线
                    try:
                        # 计算回归参数
                        slope, intercept, r_value, p_val_reg, std_err = stats.linregress(x, y)
                        line_x = np.array([x.min(), x.max()])
                        line_y = intercept + slope * line_x
                        
                        # 绘制散点图
                        ax.scatter(x, y, alpha=0.6, edgecolors='w', linewidth=0.5, s=50)
                        
                        # 绘制拟合线
                        ax.plot(line_x, line_y, 'r-', linewidth=2)
                        
                        # 设置标题
                        sig_symbol = ""
                        if p_value < 0.001:
                            sig_symbol = "***"
                        elif p_value < 0.01:
                            sig_symbol = "**"
                        elif p_value < 0.05:
                            sig_symbol = "*"
                        
                        title = f'{celltype}\nr={corr_coef:.3f}{sig_symbol} (p={p_value:.3g})'
                        ax.set_title(title, fontsize=10)
                        
                        # 设置坐标轴标签
                        ax.set_xlabel(f'{celltype} Proportion (%)')
                        ax.set_ylabel('Clinical Parameter')
                        
                        # 添加网格
                        ax.grid(True, alpha=0.3, linestyle='--')
                        
                        # 添加文本标注样本数
                        ax.text(0.05, 0.95, f'n={len(data)}', transform=ax.transAxes,
                               fontsize=9, verticalalignment='top',
                               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                        
                    except Exception as e:
                        print(f"  警告: 计算细胞类型 {celltype} 的回归线时出错: {e}")
                        # 只绘制散点图
                        ax.scatter(x, y, alpha=0.6, edgecolors='w', linewidth=0.5, s=50)
                        ax.set_title(f'{celltype} (Regression failed)', fontsize=10)
                        ax.set_xlabel(f'{celltype} Proportion (%)')
                        ax.set_ylabel('Clinical Parameter')
                        ax.grid(True, alpha=0.3, linestyle='--')
        
        # 隐藏未使用的子图
        for idx in range(n_subplots, len(axes_flat)):
            axes_flat[idx].axis('off')
        
        plt.suptitle(f'Cell Type Proportion vs Clinical Parameter (Page {page+1}/{n_pages})', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # 为suptitle留出空间
        
        # 保存当前页
        scatter_page_png = os.path.join(save_dir, f'scatter_plots_page_{page+1}.png')
        scatter_page_pdf = os.path.join(save_dir, f'scatter_plots_page_{page+1}.pdf')
        fig.savefig(scatter_page_png, dpi=300, bbox_inches='tight')
        fig.savefig(scatter_page_pdf, bbox_inches='tight')
        print(f"保存散点图第{page+1}页: {scatter_page_png}")
        print(f"保存散点图第{page+1}页: {scatter_page_pdf}")
        
        plt.show()
        plt.close(fig)
    
    # 另外，为相关性最强的几个细胞类型单独创建大图
    n_top = min(6, len(celltypes))
    top_celltypes = celltypes[:n_top]
    
    print(f"\n为前{n_top}个最相关的细胞类型创建单独的大图...")
    
    for i, celltype in enumerate(top_celltypes):
        fig, ax = plt.subplots(figsize=(8, 6))
        
        data = merged_df[[celltype, 'prop', 'patient']].dropna()
        
        if len(data) > 0:
            # 获取细胞类型相关统计
            celltype_stats = correlation_df[correlation_df['celltype'] == celltype].iloc[0]
            corr_coef = celltype_stats['correlation']
            p_value = celltype_stats['p_value']
            
            # 提取x和y值
            x = data[celltype]
            y = data['prop']
            
            # 检查x值是否有变化
            if np.std(x) == 0 or len(np.unique(x)) == 1:
                # x值没有变化，无法计算回归线
                print(f"  警告: 细胞类型 {celltype} 的x值完全相同，跳过回归线绘制")
                
                # 绘制散点图
                ax.scatter(x, y, alpha=0.7, edgecolors='k', linewidth=0.5, s=80, color='steelblue')
                
                # 设置标题
                sig_symbol = ""
                if p_value < 0.001:
                    sig_symbol = "***"
                elif p_value < 0.01:
                    sig_symbol = "**"
                elif p_value < 0.05:
                    sig_symbol = "*"
                
                title = f'{celltype} vs Clinical Parameter\nr={corr_coef:.3f}{sig_symbol} (p={p_value:.3g})'
                ax.set_title(title, fontsize=14, fontweight='bold')
                
                # 设置坐标轴标签
                ax.set_xlabel(f'{celltype} Proportion (%)', fontsize=12)
                ax.set_ylabel('Clinical Parameter', fontsize=12)
                
                # 添加网格
                ax.grid(True, alpha=0.3, linestyle='--')
                
                # 添加统计信息文本
                stats_text = f'r = {corr_coef:.3f}\np = {p_value:.3g}\nn = {len(data)}\n(No variance in x)'
                ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                       fontsize=11, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            else:
                try:
                    # 计算回归参数
                    slope, intercept, r_value, p_val_reg, std_err = stats.linregress(x, y)
                    line_x = np.array([x.min(), x.max()])
                    line_y = intercept + slope * line_x
                    
                    # 绘制散点图
                    ax.scatter(x, y, alpha=0.7, edgecolors='k', linewidth=0.5, s=80, color='steelblue')
                    
                    # 绘制拟合线
                    ax.plot(line_x, line_y, 'r-', linewidth=3, label='Linear fit')
                    
                    # 设置标题
                    sig_symbol = ""
                    if p_value < 0.001:
                        sig_symbol = "***"
                    elif p_value < 0.01:
                        sig_symbol = "**"
                    elif p_value < 0.05:
                        sig_symbol = "*"
                    
                    title = f'{celltype} vs Clinical Parameter\nr={corr_coef:.3f}{sig_symbol} (p={p_value:.3g})'
                    ax.set_title(title, fontsize=14, fontweight='bold')
                    
                    # 设置坐标轴标签
                    ax.set_xlabel(f'{celltype} Proportion (%)', fontsize=12)
                    ax.set_ylabel('Clinical Parameter', fontsize=12)
                    
                    # 添加网格
                    ax.grid(True, alpha=0.3, linestyle='--')
                    
                    # 添加统计信息文本
                    stats_text = f'r = {corr_coef:.3f}\np = {p_value:.3g}\nn = {len(data)}'
                    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                           fontsize=11, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
                    
                    # 添加图例
                    ax.legend(loc='upper right')
                except Exception as e:
                    print(f"  警告: 计算细胞类型 {celltype} 的回归线时出错: {e}")
                    # 只绘制散点图
                    ax.scatter(x, y, alpha=0.7, edgecolors='k', linewidth=0.5, s=80, color='steelblue')
                    ax.set_title(f'{celltype} vs Clinical Parameter\nr={corr_coef:.3f}', fontsize=14, fontweight='bold')
                    ax.set_xlabel(f'{celltype} Proportion (%)', fontsize=12)
                    ax.set_ylabel('Clinical Parameter', fontsize=12)
                    ax.grid(True, alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        
        # 保存单独的大图
        scatter_single_png = os.path.join(save_dir, f'scatter_{celltype}.png')
        scatter_single_pdf = os.path.join(save_dir, f'scatter_{celltype}.pdf')
        fig.savefig(scatter_single_png, dpi=300, bbox_inches='tight')
        fig.savefig(scatter_single_pdf, bbox_inches='tight')
        print(f"保存单独散点图: {scatter_single_png}")
        
        plt.show()
        plt.close(fig)

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
    
    # 4. 创建相关性摘要图
    print("\n创建相关性摘要图...")
    create_correlation_summary_plots(merged_df, correlation_df, save_dir='plot')
    
    # 5. 为每个细胞类型创建单独的散点图
    print("\n为每个细胞类型创建散点图...")
    create_individual_scatter_plots(merged_df, correlation_df, save_dir='plot', max_per_page=9)
    
    # 6. 保存结果
    print("\n保存结果...")
    save_results(correlation_df, merged_df, save_dir='plot')
    
    # 7. 显示关键结果
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
