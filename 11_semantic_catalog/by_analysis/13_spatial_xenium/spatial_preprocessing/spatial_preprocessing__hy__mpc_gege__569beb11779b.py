import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

def plot_gene_expression_by_cancertype(adata, condition_col='cancer_type', 
                                       genes=None, save_dir='gene_MPC',
                                       dpi=300, figsize=(15, 10), fontsize=12,
                                       inner='box', show_means=False,
                                       max_categories=10):
    """
    绘制指定基因在不同cancer_type中的表达小提琴图
    
    Args:
        adata: AnnData对象
        condition_col: 条件列名，默认'cancer_type'
        genes: 要分析的基因列表，如果为None则使用默认列表
        save_dir: 输出目录
        dpi: 图片分辨率
        figsize: 图片大小
        fontsize: 字体大小
        inner: 小提琴图内部显示 ('box', 'quartile', 'point', 'stick', None)
        show_means: 是否显示均值
        max_categories: 最多显示多少个cancer_type类别
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
    
    # 确保条件列存在
    if condition_col not in adata.obs.columns:
        print(f"警告: 列 '{condition_col}' 不存在")
        return
    
    # 获取所有的cancer_type类别
    cancer_types = adata.obs[condition_col].dropna().unique()
    
    if len(cancer_types) == 0:
        print(f"警告: 没有找到'{condition_col}'的数据")
        return
    
    print(f"使用 {len(adata)} 个细胞进行分析 (cancer_types: {len(cancer_types)} 种)")
    print(f"癌症类型分布:\n{adata.obs[condition_col].value_counts().to_string()}")
    
    # 如果类别太多，只取前max_categories个
    if len(cancer_types) > max_categories:
        # 按细胞数排序
        cancer_counts = adata.obs[condition_col].value_counts()
        top_cancers = cancer_counts.head(max_categories).index.tolist()
        print(f"\n注意: cancer_type 有 {len(cancer_types)} 种，只显示前 {max_categories} 种")
        print(f"显示的类型: {top_cancers}")
        cancer_types = top_cancers
    
    valid_cancers = list(cancer_types)
    
    # 检查哪些基因存在
    all_genes = adata.var_names.tolist()
    available_genes = [gene for gene in genes if gene in all_genes]
    missing_genes = [gene for gene in genes if gene not in all_genes]
    
    print(f"\n可用的基因 ({len(available_genes)}个): {available_genes}")
    if missing_genes:
        print(f"缺失的基因 ({len(missing_genes)}个): {missing_genes}")
    
    if not available_genes:
        print("没有可用的基因进行分析")
        return
    
    # 获取表达矩阵
    if adata.raw is not None:
        # 使用原始矩阵
        expr_adata = adata.raw.to_adata()
        # 确保只取当前adata的细胞
        expr_adata = expr_adata[adata.obs_names]
    else:
        expr_adata = adata
    
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
    
    # 添加cancer_type信息
    expr_df['cancer_type'] = expr_adata.obs[condition_col].values
    
    # 转换为长格式用于绘图
    expr_long = expr_df.melt(id_vars=['cancer_type'], var_name='gene', value_name='expression')
    
    # 只保留有效的cancer_type
    expr_long = expr_long[expr_long['cancer_type'].isin(valid_cancers)]
    
    # 设置颜色
    # 使用Set3颜色方案，支持多个类别
    n_colors = len(valid_cancers)
    if n_colors <= 10:
        colors = sns.color_palette("Set3", n_colors)
    else:
        colors = sns.color_palette("husl", n_colors)
    
    palette = {cancer: colors[i] for i, cancer in enumerate(valid_cancers)}
    
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
        sns.violinplot(x='cancer_type', y='expression', data=gene_data, 
                      order=valid_cancers, palette=palette, ax=ax,
                      inner=inner, linewidth=1.5, cut=0)
        
        # 如果显示均值，添加均值点
        if show_means:
            for i, cancer in enumerate(valid_cancers):
                cancer_data = gene_data[gene_data['cancer_type'] == cancer]['expression']
                if len(cancer_data) > 0:
                    mean_val = cancer_data.mean()
                    ax.plot(i, mean_val, 'o', color='white', markersize=6, 
                           markeredgecolor='black', markeredgewidth=1)
        
        # 设置标题
        ax.set_title(f'{gene}', fontsize=fontsize-2, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('Expression', fontsize=fontsize-2)
        
        # 旋转x轴标签
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        
        # 计算统计检验
        groups = {cancer: gene_data[gene_data['cancer_type'] == cancer]['expression'].values 
                  for cancer in valid_cancers}
        
        # 检查每组是否有足够数据
        groups_data = {k: v for k, v in groups.items() if len(v) > 0}
        
        if len(groups_data) >= 2:
            # Kruskal-Wallis检验（多组比较）
            try:
                kw_stat, kw_p = stats.kruskal(*[groups_data[cancer] for cancer in groups_data.keys()])
                
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
                
            except Exception as e:
                print(f"  基因 {gene} 统计检验失败: {str(e)}")
    
    # 隐藏多余的子图
    for idx in range(len(available_genes), len(axes)):
        axes[idx].set_visible(False)
    
    plt.suptitle(f'Gene Expression in Different Cancer Types (n={len(adata)} cells)', 
                 fontsize=fontsize+2, fontweight='bold')
    plt.tight_layout()
    
    # 保存组合图
    fig_path_png = os.path.join(save_dir, 'all_genes_expression_by_cancer_type.png')
    fig_path_pdf = os.path.join(save_dir, 'all_genes_expression_by_cancer_type.pdf')
    fig.savefig(fig_path_png, dpi=dpi, bbox_inches='tight')
    fig.savefig(fig_path_pdf, dpi=dpi, bbox_inches='tight')
    print(f"保存组合小提琴图: {fig_path_png}")
    print(f"保存组合小提琴图: {fig_path_pdf}")
    
    plt.show()
    plt.close(fig)
    
    # 2. 创建每个基因的单独大图
    print("\n创建每个基因的单独大图...")
    for gene in available_genes:
        fig, ax = plt.subplots(figsize=(max(10, len(valid_cancers)*0.8), 8))
        gene_data = expr_long[expr_long['gene'] == gene]
        
        # 绘制小提琴图
        sns.violinplot(x='cancer_type', y='expression', data=gene_data, 
                      order=valid_cancers, palette=palette, ax=ax, 
                      inner=inner, linewidth=2, cut=0, width=0.8)
        
        # 如果显示均值，添加均值点
        if show_means:
            for i, cancer in enumerate(valid_cancers):
                cancer_data = gene_data[gene_data['cancer_type'] == cancer]['expression']
                if len(cancer_data) > 0:
                    mean_val = cancer_data.mean()
                    ax.plot(i, mean_val, 'o', color='white', markersize=8, 
                           markeredgecolor='black', markeredgewidth=1.5)
        
        # 设置标题和标签
        ax.set_title(f'{gene} Expression in Different Cancer Types', 
                     fontsize=fontsize+2, fontweight='bold', pad=20)
        ax.set_xlabel('Cancer Type', fontsize=fontsize, fontweight='bold')
        ax.set_ylabel('Normalized Expression', fontsize=fontsize, fontweight='bold')
        
        # 旋转x轴标签
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        
        # 计算并显示统计检验
        groups = {cancer: gene_data[gene_data['cancer_type'] == cancer]['expression'].values 
                  for cancer in valid_cancers}
        
        # 检查每组是否有足够数据
        groups_data = {k: v for k, v in groups.items() if len(v) > 0}
        
        if len(groups_data) >= 2:
            # Kruskal-Wallis检验
            try:
                kw_stat, kw_p = stats.kruskal(*[groups_data[cancer] for cancer in groups_data.keys()])
                
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
                
            except Exception as e:
                print(f"  基因 {gene} 统计检验失败: {str(e)}")
        
        plt.tight_layout()
        
        # 保存单个基因图
        gene_safe = gene.replace('/', '_')
        fig_path_png = os.path.join(save_dir, f'{gene_safe}_expression_by_cancer_type.png')
        fig_path_pdf = os.path.join(save_dir, f'{gene_safe}_expression_by_cancer_type.pdf')
        fig.savefig(fig_path_png, dpi=dpi, bbox_inches='tight')
        fig.savefig(fig_path_pdf, dpi=dpi, bbox_inches='tight')
        print(f"  保存 {gene} 表达图: {fig_path_png}")
        
        plt.show()
        plt.close(fig)
    
    # 3. 创建热图展示平均表达
    print("\n创建基因表达热图...")
    
    # 计算每个基因在每个cancer_type中的平均表达
    mean_expr = expr_df.groupby('cancer_type')[available_genes].mean()
    mean_expr = mean_expr.loc[valid_cancers]  # 确保顺序
    
    # Z-score标准化
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    scaled_expr = pd.DataFrame(
        scaler.fit_transform(mean_expr.T).T,
        index=mean_expr.index,
        columns=mean_expr.columns
    )
    
    # 绘制热图
    fig, ax = plt.subplots(figsize=(max(12, len(available_genes)*0.5), max(6, len(valid_cancers)*0.5)))
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
    
    ax.set_title('Average Gene Expression (Z-score) in Different Cancer Types', 
                 fontsize=fontsize+2, fontweight='bold', pad=20)
    ax.set_xlabel('Genes', fontsize=fontsize, fontweight='bold')
    ax.set_ylabel('Cancer Types', fontsize=fontsize, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    
    # 保存热图
    fig_path_png = os.path.join(save_dir, 'gene_expression_heatmap_by_cancer_type.png')
    fig_path_pdf = os.path.join(save_dir, 'gene_expression_heatmap_by_cancer_type.pdf')
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
        fig_bar, ax_bar = plt.subplots(figsize=(max(12, len(key_genes)*0.8), 8))
        
        bar_width = 0.8 / len(valid_cancers)
        x = np.arange(len(key_genes))
        
        for i, cancer in enumerate(valid_cancers):
            cancer_mean = mean_expr.loc[cancer, key_genes].values
            ax_bar.bar(x + i*bar_width, cancer_mean, width=bar_width, 
                      label=cancer, color=palette[cancer], alpha=0.8, 
                      edgecolor='black', linewidth=1)
        
        ax_bar.set_xlabel('Genes', fontsize=fontsize, fontweight='bold')
        ax_bar.set_ylabel('Average Expression', fontsize=fontsize, fontweight='bold')
        ax_bar.set_title('Average Gene Expression in Different Cancer Types', 
                         fontsize=fontsize+2, fontweight='bold')
        ax_bar.set_xticks(x + (len(valid_cancers) - 1) * bar_width / 2)
        ax_bar.set_xticklabels(key_genes, rotation=45, ha='right')
        ax_bar.legend(title='Cancer Type', fontsize=fontsize-2, bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # 添加网格
        ax_bar.grid(True, alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        
        # 保存柱状图
        fig_path_png = os.path.join(save_dir, 'average_expression_barchart_by_cancer_type.png')
        fig_path_pdf = os.path.join(save_dir, 'average_expression_barchart_by_cancer_type.pdf')
        fig_bar.savefig(fig_path_png, dpi=dpi, bbox_inches='tight')
        fig_bar.savefig(fig_path_pdf, dpi=dpi, bbox_inches='tight')
        print(f"保存平均表达柱状图: {fig_path_png}")
        
        plt.show()
        plt.close(fig_bar)
    
    # 5. 保存统计结果
    if stats_results:
        stats_df = pd.DataFrame(stats_results)
        # 去除重复的Kruskal-Wallis检验结果
        stats_df = stats_df.drop_duplicates(subset=['gene', 'test'], keep='first')
        
        stats_file = os.path.join(save_dir, 'statistical_test_results_by_cancer_type.xlsx')
        
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
                            'KW_significance': kw_test['significance'].values[0]
                        })
            
            gene_summary_df = pd.DataFrame(gene_summary)
            gene_summary_df.to_excel(writer, sheet_name='Gene Summary', index=False)
            
            # 显著差异基因
            sig_genes = gene_summary_df[gene_summary_df['KW_p_value'] < 0.05]
            sig_genes.to_excel(writer, sheet_name='Significant Genes', index=False)
        
        print(f"保存统计结果: {stats_file}")
    
    # 6. 保存表达数据
    expr_file = os.path.join(save_dir, 'gene_expression_data_by_cancer_type.csv')
    expr_df.to_csv(expr_file, index=True)
    print(f"保存表达数据: {expr_file}")
    
    # 7. 创建汇总报告
    report_file = os.path.join(save_dir, 'analysis_report_by_cancer_type.txt')
    with open(report_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write("Gene Expression Analysis Report by Cancer Type\n")
        f.write("="*60 + "\n\n")
        
        f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Cells Analyzed: {len(adata)}\n")
        f.write(f"Cancer Types Analyzed: {len(valid_cancers)}\n")
        f.write(f"Cancer Types:\n")
        for cancer in valid_cancers:
            count = sum(expr_df['cancer_type'] == cancer)
            f.write(f"  {cancer}: {count} cells ({count/len(expr_df)*100:.1f}%)\n")
        
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
            "all_genes_expression_by_cancer_type.png/pdf - Combined violin plots of all genes",
            "individual gene expression violinplots by cancer type - Separate violin plot for each gene",
            "gene_expression_heatmap_by_cancer_type.png/pdf - Heatmap of average expression",
            "average_expression_barchart_by_cancer_type.png/pdf - Bar chart of average expression",
            "statistical_test_results_by_cancer_type.xlsx - Results of statistical tests",
            "gene_expression_data_by_cancer_type.csv - Raw expression data",
            "analysis_report_by_cancer_type.txt - This report"
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
    adata = sc.read_h5ad("filtered_with_cancer_type.h5ad")
    
    # 运行分析
    plot_gene_expression_by_cancertype(
        adata,
        condition_col='cancer_type',
        save_dir='./gene_MPC',
        dpi=300,
        figsize=(15, 10),
        fontsize=12,
        inner='box',  # 小提琴图内部显示箱线图
        show_means=True,  # 显示均值
        max_categories=10  # 最多显示10种癌症类型
    )

# 使用示例
if __name__ == "__main__":
    main()
