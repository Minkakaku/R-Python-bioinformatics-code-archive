import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.stats.multitest as multi
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
    relate_df = pd.read_csv(relate_path)
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

def merge_and_calculate_correlation(celltype_proportions, relate_df, patient_col='patient',
                                   p_cutoff=0.05, r_cutoff=0.3, apply_fdr=False):
    """
    合并数据并计算相关性，应用cutoff
    
    参数:
    ----------
    p_cutoff : float, default=0.05
        显著性阈值
    r_cutoff : float, default=0.3
        相关系数绝对值阈值
    apply_fdr : bool, default=False
        是否应用FDR校正
    """
    print("\n合并数据并计算相关性...")
    print(f"参数设置: p_cutoff={p_cutoff}, r_cutoff={r_cutoff}, apply_fdr={apply_fdr}")
    
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
                
                # 计算95%置信区间
                if n_samples >= 3:
                    # Fisher z变换
                    z = np.arctanh(corr_coef)
                    se = 1.0 / np.sqrt(n_samples - 3)
                    z_lower = z - 1.96 * se
                    z_upper = z + 1.96 * se
                    ci_lower = np.tanh(z_lower)
                    ci_upper = np.tanh(z_upper)
                else:
                    ci_lower, ci_upper = np.nan, np.nan
                
                correlation_results.append({
                    'celltype': celltype,
                    'correlation': corr_coef,
                    'p_value': p_value,
                    'n_samples': n_samples,
                    'ci_lower': ci_lower,
                    'ci_upper': ci_upper,
                    'mean_proportion': valid_data[celltype].mean(),
                    'std_proportion': valid_data[celltype].std(),
                    'min_proportion': valid_data[celltype].min(),
                    'max_proportion': valid_data[celltype].max()
                })
            except Exception as e:
                print(f"  警告: 计算细胞类型 {celltype} 的相关性时出错: {e}")
    
    # 转换为DataFrame
    correlation_df = pd.DataFrame(correlation_results)
    
    if len(correlation_df) > 0:
        # 计算绝对值相关性
        correlation_df['abs_correlation'] = correlation_df['correlation'].abs()
        
        # 应用FDR校正
        if apply_fdr and len(correlation_df) > 1:
            p_values = correlation_df['p_value'].values
            _, fdr_corrected, _, _ = multi.multipletests(p_values, alpha=p_cutoff, method='fdr_bh')
            correlation_df['fdr_corrected'] = fdr_corrected
            correlation_df['significant_fdr'] = correlation_df['fdr_corrected'] < p_cutoff
        else:
            correlation_df['significant_fdr'] = False
        
        # 标记显著性
        correlation_df['significant_raw'] = correlation_df['p_value'] < p_cutoff
        correlation_df['meets_r_cutoff'] = correlation_df['abs_correlation'] >= r_cutoff
        correlation_df['significant_both'] = correlation_df['significant_raw'] & correlation_df['meets_r_cutoff']
        
        # 按相关性绝对值排序
        correlation_df = correlation_df.sort_values('abs_correlation', ascending=False)
        
        # 统计结果
        n_sig_raw = correlation_df['significant_raw'].sum()
        n_sig_fdr = correlation_df['significant_fdr'].sum() if 'significant_fdr' in correlation_df.columns else 0
        n_meets_r = correlation_df['meets_r_cutoff'].sum()
        n_sig_both = correlation_df['significant_both'].sum()
        
        print(f"\n相关性统计:")
        print(f"  总细胞类型数: {len(correlation_df)}")
        print(f"  通过p值阈值 (p<{p_cutoff}): {n_sig_raw}")
        if apply_fdr:
            print(f"  通过FDR校正阈值: {n_sig_fdr}")
        print(f"  通过r值阈值 (|r|>={r_cutoff}): {n_meets_r}")
        print(f"  同时通过p值和r值阈值: {n_sig_both}")
    
    print(f"计算了 {len(correlation_df)} 种细胞类型的相关性")
    
    return merged_df, correlation_df

def create_correlation_summary_plots(merged_df, correlation_df, save_dir='plot_new', 
                                    p_cutoff=0.05, r_cutoff=0.3, apply_fdr=False):
    """
    创建相关性摘要图
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 检查是否有相关性结果
    if len(correlation_df) == 0:
        print("警告: 没有计算到任何相关性，跳过摘要绘图")
        return None, None
    
    # 应用cutoff筛选显著结果
    if apply_fdr and 'significant_fdr' in correlation_df.columns:
        sig_mask = correlation_df['significant_fdr']
        sig_label = f"FDR<{p_cutoff}"
    else:
        sig_mask = correlation_df['significant_raw']
        sig_label = f"p<{p_cutoff}"
    
    # 应用r cutoff
    r_mask = correlation_df['abs_correlation'] >= r_cutoff
    
    # 同时满足两个cutoff
    both_mask = sig_mask & r_mask
    
    # 获取显著细胞类型
    sig_celltypes = correlation_df[sig_mask]['celltype'].tolist()
    r_sig_celltypes = correlation_df[r_mask]['celltype'].tolist()
    both_sig_celltypes = correlation_df[both_mask]['celltype'].tolist()
    
    print(f"\n显著细胞类型统计:")
    print(f"  通过显著性阈值 ({sig_label}): {len(sig_celltypes)}")
    print(f"  通过相关系数阈值 (|r|>={r_cutoff}): {len(r_sig_celltypes)}")
    print(f"  同时通过两个阈值: {len(both_sig_celltypes)}")
    
    # 1. 创建相关性柱状图（包含cutoff线）
    fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # 左图：相关性柱状图
    # 根据显著性上色
    colors = []
    for idx, row in correlation_df.iterrows():
        if both_mask[idx]:
            colors.append('darkred' if row['correlation'] < 0 else 'darkblue')
        elif sig_mask[idx]:
            colors.append('red' if row['correlation'] < 0 else 'blue')
        elif r_mask[idx]:
            colors.append('salmon' if row['correlation'] < 0 else 'lightblue')
        else:
            colors.append('lightgray')
    
    bars = ax1.bar(range(len(correlation_df)), correlation_df['correlation'], color=colors)
    ax1.set_xlabel('Cell Type')
    ax1.set_ylabel('Correlation Coefficient')
    ax1.set_title(f'Correlation with Clinical Parameter\nCutoffs: {sig_label}, |r|>={r_cutoff}')
    ax1.set_xticks(range(len(correlation_df)))
    ax1.set_xticklabels(correlation_df['celltype'], rotation=45, ha='right')
    
    # 添加cutoff线
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax1.axhline(y=r_cutoff, color='green', linestyle='--', linewidth=1, alpha=0.7, label=f'r={r_cutoff}')
    ax1.axhline(y=-r_cutoff, color='green', linestyle='--', linewidth=1, alpha=0.7, label=f'r=-{r_cutoff}')
    
    # 在柱状图上添加数值和显著性标记
    for i, (bar, corr, p_val, celltype) in enumerate(zip(bars, correlation_df['correlation'], 
                                                          correlation_df['p_value'], correlation_df['celltype'])):
        height = bar.get_height()
        
        # 只对显著的柱状图添加数值
        if both_mask[i] or sig_mask[i]:
            ax1.text(bar.get_x() + bar.get_width()/2., height + (0.01 if height >= 0 else -0.02),
                    f'{corr:.2f}', ha='center', va='bottom' if height >= 0 else 'top', fontsize=8)
        
        # 添加显著性标记
        sig_symbol = ""
        if p_val < 0.001:
            sig_symbol = '***'
        elif p_val < 0.01:
            sig_symbol = '**'
        elif p_val < 0.05:
            sig_symbol = '*'
        
        if sig_symbol and (both_mask[i] or sig_mask[i]):
            ax1.text(bar.get_x() + bar.get_width()/2., height + (0.04 if height >= 0 else -0.05),
                    sig_symbol, ha='center', va='bottom' if height >= 0 else 'top', 
                    fontsize=10, fontweight='bold')
    
    # 添加图例说明
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='darkblue', label=f'Significant (both cutoffs), positive'),
        Patch(facecolor='darkred', label=f'Significant (both cutoffs), negative'),
        Patch(facecolor='blue', label=f'Significant ({sig_label}) only'),
        Patch(facecolor='red', label=f'Significant ({sig_label}) only, negative'),
        Patch(facecolor='lightblue', label=f'|r|>={r_cutoff} only'),
        Patch(facecolor='lightgray', label='Not significant')
    ]
    ax1.legend(handles=legend_elements, loc='upper right', fontsize=8)
    
    # 右图：-log10(p-value)柱状图
    p_value_matrix = correlation_df[['celltype', 'p_value']].copy()
    p_value_matrix['log10_p'] = -np.log10(p_value_matrix['p_value'])
    p_value_matrix = p_value_matrix.sort_values('log10_p', ascending=False)
    
    # 根据显著性上色
    p_colors = []
    for idx, row in p_value_matrix.iterrows():
        celltype = row['celltype']
        original_idx = correlation_df[correlation_df['celltype'] == celltype].index[0]
        if both_mask[original_idx]:
            p_colors.append('darkred')
        elif sig_mask[original_idx]:
            p_colors.append('red')
        else:
            p_colors.append('gray')
    
    ax2.barh(range(len(p_value_matrix)), p_value_matrix['log10_p'], color=p_colors)
    ax2.set_xlabel('-log10(p-value)')
    ax2.set_ylabel('Cell Type')
    ax2.set_title('Statistical Significance (-log10 p-value)')
    ax2.set_yticks(range(len(p_value_matrix)))
    ax2.set_yticklabels(p_value_matrix['celltype'])
    
    # 添加cutoff线
    p_cutoff_log = -np.log10(p_cutoff)
    ax2.axvline(x=p_cutoff_log, color='red', linestyle='--', linewidth=1, alpha=0.7, 
                label=f'p={p_cutoff}')
    
    # 添加FDR线（如果应用了FDR校正）
    if apply_fdr and 'fdr_corrected' in correlation_df.columns:
        fdr_cutoff_log = -np.log10(p_cutoff)
        ax2.axvline(x=fdr_cutoff_log, color='orange', linestyle=':', linewidth=1, alpha=0.7, 
                    label=f'FDR={p_cutoff}')
    
    ax2.legend()
    
    plt.suptitle(f'Cell Type Proportion Correlation Analysis\nCutoffs applied: {sig_label}, |r|>={r_cutoff}', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # 保存柱状图
    bar_plot_png = os.path.join(save_dir, 'correlation_bar_plot_with_cutoffs.png')
    bar_plot_pdf = os.path.join(save_dir, 'correlation_bar_plot_with_cutoffs.pdf')
    fig1.savefig(bar_plot_png, dpi=300, bbox_inches='tight')
    fig1.savefig(bar_plot_pdf, bbox_inches='tight')
    print(f"保存相关性柱状图: {bar_plot_png}")
    
    plt.show()
    plt.close(fig1)
    
    # 2. 创建火山图
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    
    # 准备火山图数据
    volcano_df = correlation_df.copy()
    
    # 根据显著性分组
    colors_volcano = []
    sizes_volcano = []
    for idx, row in volcano_df.iterrows():
        if both_mask[idx]:
            colors_volcano.append('red' if row['correlation'] < 0 else 'blue')
            sizes_volcano.append(100)
        elif sig_mask[idx]:
            colors_volcano.append('salmon' if row['correlation'] < 0 else 'lightblue')
            sizes_volcano.append(60)
        elif r_mask[idx]:
            colors_volcano.append('orange' if row['correlation'] < 0 else 'cyan')
            sizes_volcano.append(40)
        else:
            colors_volcano.append('gray')
            sizes_volcano.append(20)
    
    # 绘制火山图
    scatter = ax2.scatter(volcano_df['correlation'], -np.log10(volcano_df['p_value']), 
                         c=colors_volcano, s=sizes_volcano, alpha=0.7, edgecolors='black', linewidth=0.5)
    
    # 添加cutoff线
    ax2.axvline(x=r_cutoff, color='green', linestyle='--', alpha=0.5, label=f'r={r_cutoff}')
    ax2.axvline(x=-r_cutoff, color='green', linestyle='--', alpha=0.5, label=f'r=-{r_cutoff}')
    ax2.axhline(y=-np.log10(p_cutoff), color='red', linestyle='--', alpha=0.5, label=f'p={p_cutoff}')
    
    # 标记最显著的几个点
    top_n = min(5, len(volcano_df))
    for i in range(top_n):
        row = volcano_df.iloc[i]
        ax2.annotate(row['celltype'], (row['correlation'], -np.log10(row['p_value'])), 
                    fontsize=8, alpha=0.8)
    
    ax2.set_xlabel('Correlation Coefficient (r)')
    ax2.set_ylabel('-log10(p-value)')
    ax2.set_title(f'Volcano Plot of Correlation Analysis\nCutoffs: {sig_label}, |r|>={r_cutoff}')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 添加图例
    volcano_legend = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, 
                  label=f'Significant (both cutoffs), negative'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, 
                  label=f'Significant (both cutoffs), positive'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='salmon', markersize=8, 
                  label=f'Significant ({sig_label}) only'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=6, 
                  label=f'|r|>={r_cutoff} only'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=4, 
                  label='Not significant')
    ]
    ax2.legend(handles=volcano_legend, loc='upper right')
    
    plt.tight_layout()
    
    # 保存火山图
    volcano_png = os.path.join(save_dir, 'correlation_volcano_plot.png')
    volcano_pdf = os.path.join(save_dir, 'correlation_volcano_plot.pdf')
    fig2.savefig(volcano_png, dpi=300, bbox_inches='tight')
    fig2.savefig(volcano_pdf, bbox_inches='tight')
    print(f"保存相关性火山图: {volcano_png}")
    
    plt.show()
    plt.close(fig2)
    
    # 3. 创建热图（只显示显著结果）
    if len(both_sig_celltypes) > 0:
        fig3, ax3 = plt.subplots(figsize=(12, 8))
        
        # 准备热图数据
        heatmap_data = merged_df.copy()
        
        # 标准化数据
        from sklearn.preprocessing import StandardScaler
        
        # 只显示同时通过两个cutoff的细胞类型
        display_celltypes = both_sig_celltypes[:15]  # 最多显示15个
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
        
        ax3.set_title(f'Patient-Level Cell Type Proportions and Clinical Parameter (Z-scored)\nSignificant cell types only (both cutoffs)', 
                     fontsize=14, fontweight='bold')
        ax3.set_xlabel('Variables')
        ax3.set_ylabel('Patients')
        
        plt.tight_layout()
        
        # 保存热图
        heatmap_png = os.path.join(save_dir, 'correlation_heatmap_significant.png')
        heatmap_pdf = os.path.join(save_dir, 'correlation_heatmap_significant.pdf')
        fig3.savefig(heatmap_png, dpi=300, bbox_inches='tight')
        fig3.savefig(heatmap_pdf, bbox_inches='tight')
        print(f"保存显著细胞类型热图: {heatmap_png}")
        
        plt.show()
        plt.close(fig3)
    else:
        print("警告: 没有同时通过两个cutoff的细胞类型，跳过显著热图绘制")
        fig3 = None
    
    return fig1, fig2, fig3

def create_individual_scatter_plots(merged_df, correlation_df, save_dir='plot_new', 
                                   max_per_page=9, p_cutoff=0.05, r_cutoff=0.3, apply_fdr=False):
    """
    为每个细胞类型创建单独的散点图，可应用cutoff筛选
    """
    # 创建输出文件夹
    save_dir = create_output_folder(save_dir)
    
    # 检查是否有相关性结果
    if len(correlation_df) == 0:
        print("警告: 没有计算到任何相关性，跳过散点图绘制")
        return
    
    # 应用cutoff筛选
    if apply_fdr and 'significant_fdr' in correlation_df.columns:
        sig_mask = correlation_df['significant_fdr']
    else:
        sig_mask = correlation_df['significant_raw']
    
    r_mask = correlation_df['abs_correlation'] >= r_cutoff
    both_mask = sig_mask & r_mask
    
    # 获取通过cutoff的细胞类型
    if p_cutoff is not None or r_cutoff is not None:
        # 使用筛选条件
        filtered_mask = both_mask
        if not filtered_mask.any():
            print(f"警告: 没有细胞类型通过cutoff (p<{p_cutoff}且|r|>={r_cutoff})")
            # 回退到只通过p值cutoff
            filtered_mask = sig_mask
            if not filtered_mask.any():
                print("警告: 没有细胞类型通过p值cutoff，将显示所有细胞类型")
                filtered_mask = pd.Series([True] * len(correlation_df), index=correlation_df.index)
    else:
        filtered_mask = pd.Series([True] * len(correlation_df), index=correlation_df.index)
    
    filtered_correlation_df = correlation_df[filtered_mask]
    celltypes = filtered_correlation_df['celltype'].tolist()
    n_celltypes = len(celltypes)
    
    print(f"将为 {n_celltypes} 个通过cutoff的细胞类型创建散点图")
    
    if n_celltypes == 0:
        print("警告: 没有细胞类型通过筛选条件，跳过散点图绘制")
        return
    
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
                ci_lower = celltype_stats.get('ci_lower', np.nan)
                ci_upper = celltype_stats.get('ci_upper', np.nan)
                
                # 检查是否通过cutoff
                passes_cutoff = celltype in celltypes
                
                # 提取x和y值
                x = data[celltype]
                y = data['prop']
                
                # 检查x值是否有变化
                if np.std(x) == 0 or len(np.unique(x)) == 1:
                    # x值没有变化，无法计算回归线
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
                    if not passes_cutoff:
                        title += '\n(Not significant)'
                    ax.set_title(title, fontsize=10)
                    
                    # 设置坐标轴标签
                    ax.set_xlabel(f'{celltype} Proportion (%)')
                    ax.set_ylabel('Clinical Parameter')
                    
                    # 添加网格
                    ax.grid(True, alpha=0.3, linestyle='--')
                    
                    # 添加文本标注
                    ax.text(0.05, 0.95, f'n={len(data)}\n(No variance in x)', 
                           transform=ax.transAxes, fontsize=8, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                else:
                    try:
                        # 计算回归参数
                        slope, intercept, r_value, p_val_reg, std_err = stats.linregress(x, y)
                        line_x = np.array([x.min(), x.max()])
                        line_y = intercept + slope * line_x
                        
                        # 根据是否通过cutoff选择颜色
                        if passes_cutoff:
                            point_color = 'red' if corr_coef < 0 else 'blue'
                        else:
                            point_color = 'gray'
                        
                        # 绘制散点图
                        ax.scatter(x, y, alpha=0.6, edgecolors='w', linewidth=0.5, s=50, color=point_color)
                        
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
                        if not passes_cutoff:
                            title += '\n(Not significant)'
                        ax.set_title(title, fontsize=10)
                        
                        # 设置坐标轴标签
                        ax.set_xlabel(f'{celltype} Proportion (%)')
                        ax.set_ylabel('Clinical Parameter')
                        
                        # 添加网格
                        ax.grid(True, alpha=0.3, linestyle='--')
                        
                        # 添加文本标注
                        info_text = f'n={len(data)}'
                        if not np.isnan(ci_lower) and not np.isnan(ci_upper):
                            info_text += f'\n95% CI: [{ci_lower:.2f}, {ci_upper:.2f}]'
                        ax.text(0.05, 0.95, info_text, transform=ax.transAxes,
                               fontsize=8, verticalalignment='top',
                               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                        
                    except Exception as e:
                        print(f"  警告: 计算细胞类型 {celltype} 的回归线时出错: {e}")
                        ax.scatter(x, y, alpha=0.6, edgecolors='w', linewidth=0.5, s=50, color='gray')
                        ax.set_title(f'{celltype} (Regression failed)', fontsize=10)
                        ax.set_xlabel(f'{celltype} Proportion (%)')
                        ax.set_ylabel('Clinical Parameter')
                        ax.grid(True, alpha=0.3, linestyle='--')
        
        # 隐藏未使用的子图
        for idx in range(n_subplots, len(axes_flat)):
            axes_flat[idx].axis('off')
        
        # 设置总标题
        cutoff_text = f"Cutoffs: p<{p_cutoff}, |r|>={r_cutoff}" if p_cutoff and r_cutoff else "All cell types"
        plt.suptitle(f'Cell Type Proportion vs Clinical Parameter (Page {page+1}/{n_pages})\n{cutoff_text}', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        # 保存当前页
        scatter_page_png = os.path.join(save_dir, f'scatter_plots_page_{page+1}.png')
        scatter_page_pdf = os.path.join(save_dir, f'scatter_plots_page_{page+1}.pdf')
        fig.savefig(scatter_page_png, dpi=300, bbox_inches='tight')
        fig.savefig(scatter_page_pdf, bbox_inches='tight')
        print(f"保存散点图第{page+1}页: {scatter_page_png}")
        
        plt.show()
        plt.close(fig)
    
    # 为最显著的几个细胞类型创建单独大图
    n_top = min(6, len(filtered_correlation_df))
    if n_top > 0:
        top_celltypes = filtered_correlation_df.head(n_top)['celltype'].tolist()
        
        print(f"\n为前{n_top}个最显著的细胞类型创建单独的大图...")
        
        for i, celltype in enumerate(top_celltypes):
            fig, ax = plt.subplots(figsize=(8, 6))
            
            data = merged_df[[celltype, 'prop', 'patient']].dropna()
            
            if len(data) > 0:
                # 获取细胞类型相关统计
                celltype_stats = correlation_df[correlation_df['celltype'] == celltype].iloc[0]
                corr_coef = celltype_stats['correlation']
                p_value = celltype_stats['p_value']
                ci_lower = celltype_stats.get('ci_lower', np.nan)
                ci_upper = celltype_stats.get('ci_upper', np.nan)
                
                # 提取x和y值
                x = data[celltype]
                y = data['prop']
                
                # 检查x值是否有变化
                if np.std(x) > 0 and len(np.unique(x)) > 1:
                    try:
                        # 计算回归参数
                        slope, intercept, r_value, p_val_reg, std_err = stats.linregress(x, y)
                        line_x = np.array([x.min(), x.max()])
                        line_y = intercept + slope * line_x
                        
                        # 绘制散点图
                        ax.scatter(x, y, alpha=0.7, edgecolors='k', linewidth=0.5, s=80, 
                                  color='red' if corr_coef < 0 else 'blue')
                        
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
                        if not np.isnan(ci_lower) and not np.isnan(ci_upper):
                            stats_text += f'\n95% CI: [{ci_lower:.2f}, {ci_upper:.2f}]'
                        ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                               fontsize=11, verticalalignment='top',
                               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
                        
                        # 添加图例
                        ax.legend(loc='upper right')
                        
                    except Exception as e:
                        print(f"  警告: 计算细胞类型 {celltype} 的回归线时出错: {e}")
            
            plt.tight_layout()
            
            # 保存单独的大图
            scatter_single_png = os.path.join(save_dir, f'scatter_{celltype}.png')
            scatter_single_pdf = os.path.join(save_dir, f'scatter_{celltype}.pdf')
            fig.savefig(scatter_single_png, dpi=300, bbox_inches='tight')
            fig.savefig(scatter_single_pdf, bbox_inches='tight')
            print(f"保存单独散点图: {scatter_single_png}")
            
            plt.show()
            plt.close(fig)

def save_results(correlation_df, merged_df, save_dir='plot_new', p_cutoff=0.05, r_cutoff=0.3, apply_fdr=False):
    """
    保存结果到文件
    """
    save_dir = create_output_folder(save_dir)
    
    # 应用cutoff标记
    if apply_fdr and 'significant_fdr' in correlation_df.columns:
        sig_mask = correlation_df['significant_fdr']
        sig_label = f"FDR<{p_cutoff}"
    else:
        sig_mask = correlation_df['significant_raw']
        sig_label = f"p<{p_cutoff}"
    
    r_mask = correlation_df['abs_correlation'] >= r_cutoff
    both_mask = sig_mask & r_mask
    
    # 标记cutoff状态
    correlation_df['passes_p_cutoff'] = sig_mask
    correlation_df['passes_r_cutoff'] = r_mask
    correlation_df['passes_both_cutoffs'] = both_mask
    
    # 保存所有相关性结果
    correlation_file = os.path.join(save_dir, 'correlation_results_all.csv')
    correlation_df.to_csv(correlation_file, index=False)
    print(f"保存所有相关性结果: {correlation_file}")
    
    # 保存通过cutoff的相关性结果
    if both_mask.any():
        filtered_correlation_df = correlation_df[both_mask]
        filtered_correlation_file = os.path.join(save_dir, 'correlation_results_significant.csv')
        filtered_correlation_df.to_csv(filtered_correlation_file, index=False)
        print(f"保存显著相关性结果: {filtered_correlation_file}")
    
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
        
        f.write("1. ANALYSIS PARAMETERS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Significance cutoff (p-value): {p_cutoff}\n")
        f.write(f"Correlation strength cutoff (|r|): {r_cutoff}\n")
        f.write(f"Apply FDR correction: {apply_fdr}\n")
        if apply_fdr:
            f.write(f"FDR method: Benjamini-Hochberg\n")
        f.write(f"Total patients: {merged_df['patient'].nunique()}\n")
        f.write(f"Total cell types analyzed: {len(correlation_df)}\n")
        f.write(f"Clinical parameter column: 'prop'\n\n")
        
        f.write("2. CORRELATION RESULTS SUMMARY\n")
        f.write("-" * 40 + "\n")
        
        if len(correlation_df) > 0:
            # 统计结果
            n_sig_p = sig_mask.sum()
            n_sig_r = r_mask.sum()
            n_sig_both = both_mask.sum()
            
            f.write(f"Cell types passing p-value cutoff ({sig_label}): {n_sig_p}\n")
            f.write(f"Cell types passing correlation cutoff (|r|>={r_cutoff}): {n_sig_r}\n")
            f.write(f"Cell types passing both cutoffs: {n_sig_both}\n\n")
            
            f.write("3. TOP CORRELATIONS (all cell types)\n")
            f.write("-" * 40 + "\n")
            
            # 显示所有细胞类型的top 20
            f.write("Top 20 correlations by absolute value:\n")
            f.write("-" * 30 + "\n")
            for i, row in correlation_df.head(20).iterrows():
                sig_mark = ""
                if row['p_value'] < 0.001:
                    sig_mark = "***"
                elif row['p_value'] < 0.01:
                    sig_mark = "**"
                elif row['p_value'] < 0.05:
                    sig_mark = "*"
                
                passes = "PASS" if row['passes_both_cutoffs'] else "FAIL"
                f.write(f"{row['celltype']}: r = {row['correlation']:.3f}{sig_mark} (p = {row['p_value']:.4f}) [{passes}]\n")
            
            f.write("\n4. SIGNIFICANT CORRELATIONS (passing both cutoffs)\n")
            f.write("-" * 40 + "\n")
            
            if n_sig_both > 0:
                f.write(f"Cell types passing both cutoffs (p<{p_cutoff} and |r|>={r_cutoff}):\n")
                f.write("-" * 30 + "\n")
                for i, row in correlation_df[both_mask].iterrows():
                    sig_mark = ""
                    if row['p_value'] < 0.001:
                        sig_mark = "***"
                    elif row['p_value'] < 0.01:
                        sig_mark = "**"
                    elif row['p_value'] < 0.05:
                        sig_mark = "*"
                    
                    f.write(f"{row['celltype']}: r = {row['correlation']:.3f}{sig_mark} (p = {row['p_value']:.4f})\n")
                    f.write(f"  n = {row['n_samples']}, 95% CI = [{row.get('ci_lower', 'NA'):.2f}, {row.get('ci_upper', 'NA'):.2f}]\n")
                    f.write(f"  Proportion: mean = {row.get('mean_proportion', 'NA'):.2f}%, ")
                    f.write(f"range = [{row.get('min_proportion', 'NA'):.2f}, {row.get('max_proportion', 'NA'):.2f}]%\n")
            else:
                f.write("No cell types pass both cutoffs.\n")
            
            f.write("\n5. STATISTICAL SUMMARY\n")
            f.write("-" * 40 + "\n")
            f.write(f"Mean absolute correlation: {correlation_df['abs_correlation'].mean():.3f}\n")
            f.write(f"Median absolute correlation: {correlation_df['abs_correlation'].median():.3f}\n")
            f.write(f"Maximum correlation: {correlation_df['correlation'].max():.3f}\n")
            f.write(f"Minimum correlation: {correlation_df['correlation'].min():.3f}\n")
            f.write(f"Number with p < 0.05: {(correlation_df['p_value'] < 0.05).sum()}\n")
            f.write(f"Number with p < 0.01: {(correlation_df['p_value'] < 0.01).sum()}\n")
            f.write(f"Number with p < 0.001: {(correlation_df['p_value'] < 0.001).sum()}\n")
            f.write(f"Number with |r| >= 0.3: {(correlation_df['abs_correlation'] >= 0.3).sum()}\n")
            f.write(f"Number with |r| >= 0.5: {(correlation_df['abs_correlation'] >= 0.5).sum()}\n")
        else:
            f.write("No correlation results available.\n")
        
        f.write("\n" + "=" * 60 + "\n")
        f.write("Analysis completed successfully!\n")
        f.write("=" * 60 + "\n")
    
    print(f"保存统计摘要: {summary_file}")
    
    return correlation_file, summary_file

def main():
    """
    主函数
    """
    # 文件路径
    adata_path = 'adata_with_band_by_sample_improved.h5ad'
    relate_path = 'releate.csv'
    
    # 检查文件是否存在
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"未找到h5ad文件: {adata_path}")
    if not os.path.exists(relate_path):
        raise FileNotFoundError(f"未找到相关文件: {relate_path}")
    
    # 参数设置
    p_cutoff = 0.05  # 显著性阈值
    r_cutoff = 0.3   # 相关系数绝对值阈值
    apply_fdr = False  # 是否应用FDR校正
    
    print(f"分析参数:")
    print(f"  p-value cutoff: {p_cutoff}")
    print(f"  correlation cutoff (|r|): {r_cutoff}")
    print(f"  Apply FDR correction: {apply_fdr}")
    
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
        patient_col='patient',
        p_cutoff=p_cutoff,
        r_cutoff=r_cutoff,
        apply_fdr=apply_fdr
    )
    
    # 4. 创建相关性摘要图
    print("\n创建相关性摘要图...")
    create_correlation_summary_plots(
        merged_df, 
        correlation_df, 
        save_dir='plot_new',
        p_cutoff=p_cutoff,
        r_cutoff=r_cutoff,
        apply_fdr=apply_fdr
    )
    
    # 5. 为每个细胞类型创建单独的散点图
    print("\n为每个细胞类型创建散点图...")
    create_individual_scatter_plots(
        merged_df, 
        correlation_df, 
        save_dir='plot_new', 
        max_per_page=9,
        p_cutoff=p_cutoff,
        r_cutoff=r_cutoff,
        apply_fdr=apply_fdr
    )
    
    # 6. 保存结果
    print("\n保存结果...")
    save_results(
        correlation_df, 
        merged_df, 
        save_dir='plot_new',
        p_cutoff=p_cutoff,
        r_cutoff=r_cutoff,
        apply_fdr=apply_fdr
    )
    
    # 7. 显示关键结果
    print("\n" + "=" * 60)
    print("关键结果摘要")
    print("=" * 60)
    
    if len(correlation_df) > 0:
        # 应用cutoff筛选
        if apply_fdr and 'significant_fdr' in correlation_df.columns:
            sig_mask = correlation_df['significant_fdr']
        else:
            sig_mask = correlation_df['significant_raw']
        
        r_mask = correlation_df['abs_correlation'] >= r_cutoff
        both_mask = sig_mask & r_mask
        
        n_sig_both = both_mask.sum()
        
        if n_sig_both > 0:
            print(f"\n通过双重cutoff (p<{p_cutoff}且|r|>={r_cutoff}) 的细胞类型 ({n_sig_both}个):")
            print("-" * 40)
            for i, row in correlation_df[both_mask].iterrows():
                sig_symbol = ""
                if row['p_value'] < 0.001:
                    sig_symbol = "***"
                elif row['p_value'] < 0.01:
                    sig_symbol = "**"
                elif row['p_value'] < 0.05:
                    sig_symbol = "*"
                
                print(f"{row['celltype']}: r = {row['correlation']:.3f} {sig_symbol} (p = {row['p_value']:.4f})")
        else:
            print(f"\n没有细胞类型通过双重cutoff (p<{p_cutoff}且|r|>={r_cutoff})")
            
            # 显示只通过p值cutoff的
            n_sig_p = sig_mask.sum()
            if n_sig_p > 0:
                print(f"\n只通过p值cutoff (p<{p_cutoff}) 的细胞类型 ({n_sig_p}个):")
                print("-" * 40)
                for i, row in correlation_df[sig_mask].head(5).iterrows():
                    sig_symbol = ""
                    if row['p_value'] < 0.001:
                        sig_symbol = "***"
                    elif row['p_value'] < 0.01:
                        sig_symbol = "**"
                    elif row['p_value'] < 0.05:
                        sig_symbol = "*"
                    
                    print(f"{row['celltype']}: r = {row['correlation']:.3f} {sig_symbol} (p = {row['p_value']:.4f})")
            
            # 显示只通过r值cutoff的
            n_sig_r = r_mask.sum()
            if n_sig_r > 0:
                print(f"\n只通过r值cutoff (|r|>={r_cutoff}) 的细胞类型 ({n_sig_r}个):")
                print("-" * 40)
                for i, row in correlation_df[r_mask].head(5).iterrows():
                    print(f"{row['celltype']}: r = {row['correlation']:.3f} (p = {row['p_value']:.4f})")
    else:
        print("\n没有计算到相关性结果")
    
    print("\n分析完成！所有结果已保存到 'plot' 文件夹。")
    
    return adata, relate_df, merged_df, correlation_df

# 运行主函数
if __name__ == "__main__":
    adata, relate_df, merged_df, correlation_df = main()
