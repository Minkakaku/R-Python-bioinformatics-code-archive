import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

def plot_gene_expression_violinplots(adata, condition_col='site_band', 
                                     genes=None, save_dir='gene_expression_analysis',
                                     dpi=300, figsize=(15, 10), fontsize=12,
                                     inner='box', show_means=False):
    """
    绘制指定基因在in, on, out三个条带中的表达小提琴图
    
    Args:
        adata: AnnData对象
        condition_col: 条件列名，默认'site_band'
        genes: 要分析的基因列表，如果为None则使用默认列表
        save_dir: 输出目录
        dpi: 图片分辨率
        figsize: 图片大小
        fontsize: 字体大小
        inner: 小提琴图内部显示 ('box', 'quartile', 'point', 'stick', None)
        show_means: 是否显示均值
    """
    
    # 创建输出目录
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    # 默认基因列表
    if genes is None:
        genes = [
            'MUC2', 'MUC5AC', 'MUC5B', 'MUC19', 'SPDEF', 
            'FOXA2', 'ATOH1', 'AGR2', 'EGFR', 'EGF', 
            'MAPK', 'IL4', 'IL13', 'STAT6', 'IL1',
            'IL17', 'NFKB', 'TGFB1', 'CD44', 'HA',
            'S100A8', 'S100A9', 'S100A12', 'TP63', 'GATA3'
        ]
    
    # 只保留in, on, out的数据
    valid_bands = ['in', 'on', 'out']
    mask = adata.obs[condition_col].isin(valid_bands)
    adata_filtered = adata[mask].copy()
    
    if len(adata_filtered) == 0:
        print(f"警告: 没有找到'site_band'为{valid_bands}的细胞")
        return
    
    print(f"使用 {len(adata_filtered)} 个细胞进行分析 (bands: {valid_bands})")
    print(f"条带分布: {adata_filtered.obs[condition_col].value_counts().to_dict()}")
    
    # 检查哪些基因存在
    all_genes = adata_filtered.var_names.tolist()
    available_genes = [gene for gene in genes if gene in all_genes]
    missing_genes = [gene for gene in genes if gene not in all_genes]
    
    print(f"\n可用的基因 ({len(available_genes)}个): {available_genes}")
    if missing_genes:
        print(f"缺失的基因 ({len(missing_genes)}个): {missing_genes}")
    
    if not available_genes:
        print("没有可用的基因进行分析")
        return
    
    # 获取表达矩阵
    if adata_filtered.raw is not None:
        # 使用原始矩阵
        expr_adata = adata_filtered.raw.to_adata()
        # 确保只取当前adata的细胞
        expr_adata = expr_adata[adata_filtered.obs_names]
    else:
        expr_adata = adata_filtered
    
    # 获取表达数据
    if hasattr(expr_adata.X, 'toarray'):
        expr_data = expr_adata[:, available_genes].X.toarray()
    else:
        expr_data = expr_adata[:, available_genes].X
    
    # 创建表达数据框
    expr_df = pd.DataFrame(
        expr_data,
        columns=available_genes,
        index=expr_adata.obs_names
    )
    
    # 添加条带信息
    expr_df['band'] = expr_adata.obs[condition_col].values
    
    # 转换为长格式用于绘图
    expr_long = expr_df.melt(id_vars=['band'], var_name='gene', value_name='expression')
    
    # 设置颜色
    palette = {'in': '#1f77b4', 'on': '#2ca02c', 'out': '#d62728'}
    
    # 创建统计结果数据框
    stats_results = []
    
    # 1. 创建所有基因的组合小提琴图
    print("\n创建所有基因的组合小提琴图...")
    n_genes = len(available_genes)
    n_cols = 5
    n_rows = int(np.ceil(n_genes / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(figsize[0], 2.5*n_rows))
    
    # 处理axes的形状
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    elif n_rows == 1 or n_cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    for idx, gene in enumerate(available_genes):
        if idx >= len(axes):
            break
            
        ax = axes[idx]
        gene_data = expr_long[expr_long['gene'] == gene]
        
        # 绘制小提琴图
        sns.violinplot(x='band', y='expression', data=gene_data, 
                      order=valid_bands, palette=palette, ax=ax,
                      inner=inner, linewidth=1.5, cut=0)
        
        # 如果显示均值，添加均值点
        if show_means:
            for i, band in enumerate(valid_bands):
                band_data = gene_data[gene_data['band'] == band]['expression']
                if len(band_data) > 0:
                    mean_val = band_data.mean()
                    ax.plot(i, mean_val, 'o', color='white', markersize=6, 
                           markeredgecolor='black', markeredgewidth=1)
        
        # 设置标题
        ax.set_title(f'{gene}', fontsize=fontsize-2, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('Expression', fontsize=fontsize-2)
        
        # 旋转x轴标签
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        
        # 计算统计检验
        groups = {band: gene_data[gene_data['band'] == band]['expression'].values 
                  for band in valid_bands}
        
        # 检查每组是否有数据
        groups_data = {k: v for k, v in groups.items() if len(v) > 0}
        
        if len(groups_data) >= 2:
            # Kruskal-Wallis检验（多组比较）
            try:
                kw_stat, kw_p = stats.kruskal(*[groups_data[band] for band in groups_data.keys()])
                
                # 在图上添加p值
                if kw_p < 0.001:
                    p_text = 'p < 0.001'
                elif kw_p < 0.01:
                    p_text = 'p < 0.01'
                elif kw_p < 0.05:
                    p_text = 'p < 0.05'
                else:
                    p_text = f'p = {kw_p:.3f}'
                
                # 添加统计结果
                stats_results.append({
                    'gene': gene,
                    'test': 'Kruskal-Wallis',
                    'statistic': kw_stat,
                    'p_value': kw_p,
                    'significance': '***' if kw_p < 0.001 else ('**' if kw_p < 0.01 else ('*' if kw_p < 0.05 else 'ns'))
                })
                
                # 在图上添加p值
                ax.text(0.5, 0.95, p_text, transform=ax.transAxes, 
                        ha='center', va='top', fontsize=fontsize-4,
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                
                # 如果p值小于0.05，进行两两比较
                if kw_p < 0.05 and len(groups_data) == 3:
                    # Mann-Whitney U检验（两两比较）
                    comparisons = [('in', 'on'), ('in', 'out'), ('on', 'out')]
                    
                    for i, (group1, group2) in enumerate(comparisons):
                        if group1 in groups_data and group2 in groups_data:
                            try:
                                mw_stat, mw_p = stats.mannwhitneyu(
                                    groups_data[group1], 
                                    groups_data[group2],
                                    alternative='two-sided'
                                )
                                
                                # 添加统计结果
                                stats_results.append({
                                    'gene': gene,
                                    'test': f'Mann-Whitney ({group1} vs {group2})',
                                    'statistic': mw_stat,
                                    'p_value': mw_p,
                                    'significance': '***' if mw_p < 0.001 else ('**' if mw_p < 0.01 else ('*' if mw_p < 0.05 else 'ns'))
                                })
                            except:
                                pass
            except:
                pass
    
    # 隐藏多余的子图
    for idx in range(len(available_genes), len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle(f'Gene Expression in Different Bands (n={len(adata_filtered)} cells)', 
                 fontsize=fontsize+2, fontweight='bold')
    plt.tight_layout()
    
    # 保存组合图
    fig_path_png = os.path.join(save_dir, 'all_genes_expression_violinplots.png')
    fig_path_pdf = os.path.join(save_dir, 'all_genes_expression_violinplots.pdf')
    fig.savefig(fig_path_png, dpi=dpi, bbox_inches='tight')
    fig.savefig(fig_path_pdf, dpi=dpi, bbox_inches='tight')
    print(f"保存组合小提琴图: {fig_path_png}")
    print(f"保存组合小提琴图: {fig_path_pdf}")
    
    plt.show()
    plt.close(fig)
    
    # 2. 创建每个基因的单独大图
    print("\n创建每个基因的单独大图...")
    for gene in available_genes:
        fig, ax = plt.subplots(figsize=(10, 8))
        gene_data = expr_long[expr_long['gene'] == gene]
        
        # 绘制小提琴图
        sns.violinplot(x='band', y='expression', data=gene_data, 
                      order=valid_bands, palette=palette, ax=ax, 
                      inner=inner, linewidth=2, cut=0, width=0.8)
        
        # 如果显示均值，添加均值点
        if show_means:
            for i, band in enumerate(valid_bands):
                band_data = gene_data[gene_data['band'] == band]['expression']
                if len(band_data) > 0:
                    mean_val = band_data.mean()
                    ax.plot(i, mean_val, 'o', color='white', markersize=8, 
                           markeredgecolor='black', markeredgewidth=1.5)
        
        # 设置标题和标签
        ax.set_title(f'{gene} Expression in Different Bands', 
                     fontsize=fontsize+2, fontweight='bold', pad=20)
        ax.set_xlabel('Band', fontsize=fontsize, fontweight='bold')
        ax.set_ylabel('Normalized Expression', fontsize=fontsize, fontweight='bold')
        
        # 计算并显示统计检验
        groups = {band: gene_data[gene_data['band'] == band]['expression'].values 
                  for band in valid_bands}
        
        # 检查每组是否有数据
        groups_data = {k: v for k, v in groups.items() if len(v) > 0}
        
        if len(groups_data) >= 2:
            # Kruskal-Wallis检验
            try:
                kw_stat, kw_p = stats.kruskal(*[groups_data[band] for band in groups_data.keys()])
                
                # 在图上添加p值
                p_text = f'Kruskal-Wallis test\nH = {kw_stat:.2f}, p = {kw_p:.4f}'
                if kw_p < 0.001:
                    p_text += ' ***'
                elif kw_p < 0.01:
                    p_text += ' **'
                elif kw_p < 0.05:
                    p_text += ' *'
                
                ax.text(0.5, 0.95, p_text, transform=ax.transAxes, 
                        ha='center', va='top', fontsize=fontsize-2,
                        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
                
                # 添加统计结果
                stats_results.append({
                    'gene': gene,
                    'test': 'Kruskal-Wallis',
                    'statistic': kw_stat,
                    'p_value': kw_p,
                    'significance': '***' if kw_p < 0.001 else ('**' if kw_p < 0.01 else ('*' if kw_p < 0.05 else 'ns'))
                })
                
                # 如果p值小于0.05，进行两两比较并在图上标记
                if kw_p < 0.05 and len(groups_data) == 3:
                    # 准备两两比较
                    comparisons = [('in', 'on'), ('in', 'out'), ('on', 'out')]
                    p_values = []
                    
                    for group1, group2 in comparisons:
                        if group1 in groups_data and group2 in groups_data:
                            try:
                                mw_stat, mw_p = stats.mannwhitneyu(
                                    groups_data[group1], 
                                    groups_data[group2],
                                    alternative='two-sided'
                                )
                                p_values.append((group1, group2, mw_p))
                                
                                # 添加统计结果
                                stats_results.append({
                                    'gene': gene,
                                    'test': f'Mann-Whitney ({group1} vs {group2})',
                                    'statistic': mw_stat,
                                    'p_value': mw_p,
                                    'significance': '***' if mw_p < 0.001 else ('**' if mw_p < 0.01 else ('*' if mw_p < 0.05 else 'ns'))
                                })
                            except:
                                p_values.append((group1, group2, np.nan))
                    
                    # 在图上添加两两比较结果
                    y_max = gene_data['expression'].max()
                    y_min = gene_data['expression'].min()
                    y_range = y_max - y_min
                    
                    # 定义不同比较的垂直位置，错开排列避免遮挡
                    bracket_positions = {
                        ('in', 'on'): 1,     # 中层
                        ('on', 'out'): 2,    # 上层
                        ('in', 'out'): 0     # 最上层
                    }
                    
                    # 计算每个比较的y位置
                    bracket_heights = {}
                    max_bracket_level = max(bracket_positions.values())
                    base_height = 0.05 * y_range
                    
                    for (group1, group2, p_val) in p_values:
                        if not np.isnan(p_val):
                            level = bracket_positions.get((group1, group2), 1)
                            # 不同层级的比较有不同的垂直位置
                            y = y_max + (level + 1) * base_height
                            bracket_heights[(group1, group2)] = y
                    
                    # 绘制每个比较的显著性标记
                    for (group1, group2, p_val) in p_values:
                        if not np.isnan(p_val) and (group1, group2) in bracket_heights:
                            x1 = valid_bands.index(group1)
                            x2 = valid_bands.index(group2)
                            y = bracket_heights[(group1, group2)]
                            
                            # 绘制括号
                            bracket_height = 0.02 * y_range
                            ax.plot([x1, x1, x2, x2], 
                                    [y, y+bracket_height, y+bracket_height, y], 
                                    color='black', lw=1.5)
                            
                            # 添加p值
                            if p_val < 0.001:
                                p_text = '***'
                            elif p_val < 0.01:
                                p_text = '**'
                            elif p_val < 0.05:
                                p_text = '*'
                            else:
                                p_text = f'p={p_val:.3f}'
                            
                            # 将p值放在括号上方
                            ax.text((x1+x2)/2, y+bracket_height*1.5, p_text, 
                                    ha='center', va='bottom', fontsize=fontsize-2, fontweight='bold')
            except Exception as e:
                print(f"  基因 {gene} 统计检验失败: {str(e)}")
        
        plt.tight_layout()
        
        # 保存单个基因图
        gene_safe = gene.replace('/', '_')
        fig_path_png = os.path.join(save_dir, f'{gene_safe}_expression_violinplot.png')
        fig_path_pdf = os.path.join(save_dir, f'{gene_safe}_expression_violinplot.pdf')
        fig.savefig(fig_path_png, dpi=dpi, bbox_inches='tight')
        fig.savefig(fig_path_pdf, dpi=dpi, bbox_inches='tight')
        print(f"  保存 {gene} 表达图: {fig_path_png}")
        
        plt.show()
        plt.close(fig)
    
    # 3. 创建热图展示平均表达
    print("\n创建基因表达热图...")
    
    # 计算每个基因在每个条带中的平均表达
    mean_expr = expr_df.groupby('band')[available_genes].mean()
    mean_expr = mean_expr.loc[valid_bands]  # 确保顺序
    
    # Z-score标准化
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    scaled_expr = pd.DataFrame(
        scaler.fit_transform(mean_expr.T).T,
        index=mean_expr.index,
        columns=mean_expr.columns
    )
    
    # 绘制热图
    fig, ax = plt.subplots(figsize=(max(12, len(available_genes)*0.5), 6))
    sns.heatmap(
        scaled_expr,
        cmap='coolwarm',
        center=0,
        annot=True,
        fmt='.2f',
        linewidths=0.5,
        cbar_kws={'label': 'Z-score'},
        square=True,
        ax=ax
    )
    
    ax.set_title('Average Gene Expression (Z-score) in Different Bands', 
                 fontsize=fontsize+2, fontweight='bold', pad=20)
    ax.set_xlabel('Genes', fontsize=fontsize, fontweight='bold')
    ax.set_ylabel('Bands', fontsize=fontsize, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    
    # 保存热图
    fig_path_png = os.path.join(save_dir, 'gene_expression_heatmap.png')
    fig_path_pdf = os.path.join(save_dir, 'gene_expression_heatmap.pdf')
    fig.savefig(fig_path_png, dpi=dpi, bbox_inches='tight')
    fig.savefig(fig_path_pdf, dpi=dpi, bbox_inches='tight')
    print(f"保存表达热图: {fig_path_png}")
    print(f"保存表达热图: {fig_path_pdf}")
    
    plt.show()
    plt.close(fig)
    
    # 4. 创建分组柱状图
    print("\n创建分组柱状图...")
    
    # 选择几个关键基因
    key_genes = ['MUC2', 'MUC5AC', 'MUC5B', 'EGFR', 'IL4', 'IL13', 'TP63', 'GATA3']
    key_genes = [gene for gene in key_genes if gene in available_genes]
    
    if key_genes:
        fig_bar, ax_bar = plt.subplots(figsize=(12, 8))
        
        bar_width = 0.25
        x = np.arange(len(key_genes))
        
        for i, band in enumerate(valid_bands):
            band_mean = mean_expr.loc[band, key_genes].values
            ax_bar.bar(x + i*bar_width, band_mean, width=bar_width, 
                      label=band, color=palette[band], alpha=0.8, edgecolor='black', linewidth=1.5)
        
        ax_bar.set_xlabel('Genes', fontsize=fontsize, fontweight='bold')
        ax_bar.set_ylabel('Average Expression', fontsize=fontsize, fontweight='bold')
        ax_bar.set_title('Average Gene Expression in Different Bands', 
                         fontsize=fontsize+2, fontweight='bold')
        ax_bar.set_xticks(x + bar_width)
        ax_bar.set_xticklabels(key_genes, rotation=45, ha='right')
        ax_bar.legend(title='Band', fontsize=fontsize-1)
        
        # 添加网格
        ax_bar.grid(True, alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        
        # 保存柱状图
        fig_path_png = os.path.join(save_dir, 'average_expression_barchart.png')
        fig_path_pdf = os.path.join(save_dir, 'average_expression_barchart.pdf')
        fig_bar.savefig(fig_path_png, dpi=dpi, bbox_inches='tight')
        fig_bar.savefig(fig_path_pdf, dpi=dpi, bbox_inches='tight')
        print(f"保存平均表达柱状图: {fig_path_png}")
        
        plt.show()
        plt.close(fig_bar)
    
    # 5. 保存统计结果
    if stats_results:
        stats_df = pd.DataFrame(stats_results)
        stats_file = os.path.join(save_dir, 'statistical_test_results.xlsx')
        
        # 使用Excel保存，可以保存多个sheet
        with pd.ExcelWriter(stats_file, engine='openpyxl') as writer:
            # 总体统计
            stats_df.to_excel(writer, sheet_name='All Tests', index=False)
            
            # 按基因汇总
            gene_summary = []
            for gene in available_genes:
                gene_stats = stats_df[stats_df['gene'] == gene]
                if not gene_stats.empty:
                    kw_test = gene_stats[gene_stats['test'] == 'Kruskal-Wallis']
                    if not kw_test.empty:
                        gene_summary.append({
                            'gene': gene,
                            'KW_statistic': kw_test['statistic'].values[0],
                            'KW_p_value': kw_test['p_value'].values[0],
                            'KW_significance': kw_test['significance'].values[0],
                            'num_significant_pairs': sum(gene_stats['p_value'] < 0.05) - (1 if not kw_test.empty else 0)
                        })
            
            gene_summary_df = pd.DataFrame(gene_summary)
            gene_summary_df.to_excel(writer, sheet_name='Gene Summary', index=False)
            
            # 显著差异基因
            sig_genes = gene_summary_df[gene_summary_df['KW_p_value'] < 0.05]
            sig_genes.to_excel(writer, sheet_name='Significant Genes', index=False)
        
        print(f"保存统计结果: {stats_file}")
    
    # 6. 保存表达数据
    expr_file = os.path.join(save_dir, 'gene_expression_data.csv')
    expr_df.to_csv(expr_file, index=True)
    print(f"保存表达数据: {expr_file}")
    
    # 7. 创建汇总报告
    report_file = os.path.join(save_dir, 'analysis_report.txt')
    with open(report_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write("Gene Expression Analysis Report (Violin Plots)\n")
        f.write("="*60 + "\n\n")
        
        f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Cells Analyzed: {len(adata_filtered)}\n")
        f.write(f"Bands Distribution:\n")
        for band in valid_bands:
            count = sum(expr_df['band'] == band)
            f.write(f"  {band}: {count} cells ({count/len(expr_df)*100:.1f}%)\n")
        
        f.write(f"\nGenes Analyzed: {len(available_genes)}/{len(genes)}\n")
        f.write(f"Available Genes: {', '.join(available_genes)}\n")
        if missing_genes:
            f.write(f"Missing Genes: {', '.join(missing_genes)}\n")
        
        f.write("\n" + "="*60 + "\n")
        f.write("Statistical Significance Summary\n")
        f.write("="*60 + "\n\n")
        
        if 'gene_summary_df' in locals() and not gene_summary_df.empty:
            sig_count = sum(gene_summary_df['KW_p_value'] < 0.05)
            f.write(f"Significantly Different Genes: {sig_count}/{len(gene_summary_df)} ({sig_count/len(gene_summary_df)*100:.1f}%)\n\n")
            
            f.write("Top 10 Most Significant Genes:\n")
            top_genes = gene_summary_df.nsmallest(10, 'KW_p_value')
            for _, row in top_genes.iterrows():
                f.write(f"  {row['gene']}: p = {row['KW_p_value']:.6f} {row['KW_significance']}\n")
        
        f.write("\n" + "="*60 + "\n")
        f.write("Files Generated\n")
        f.write("="*60 + "\n\n")
        
        files = [
            "all_genes_expression_violinplots.png/pdf - Combined violin plots of all genes",
            "individual gene expression violinplots - Separate violin plot for each gene",
            "gene_expression_heatmap.png/pdf - Heatmap of average expression",
            "average_expression_barchart.png/pdf - Bar chart of average expression",
            "statistical_test_results.xlsx - Results of statistical tests",
            "gene_expression_data.csv - Raw expression data",
            "analysis_report.txt - This report"
        ]
        
        for file_desc in files:
            f.write(f"  • {file_desc}\n")
        
        f.write("\n" + "="*60 + "\n")
        f.write("Analysis Complete\n")
        f.write("="*60 + "\n")
    
    print(f"保存分析报告: {report_file}")
    print(f"\n分析完成！所有结果已保存到: {os.path.abspath(save_dir)}")

# 主程序
def main():
    # 加载数据
    adata = sc.read_h5ad("adata_with_band_by_sample.h5ad")
    
    # 运行分析
    plot_gene_expression_violinplots(
        adata,
        condition_col='site_band',
        save_dir='./gene_expression_analysis_violin',
        dpi=300,
        figsize=(15, 10),
        fontsize=12,
        inner='box',  # 小提琴图内部显示箱线图
        show_means=True  # 显示均值
    )

# 使用示例
if __name__ == "__main__":
    main()
