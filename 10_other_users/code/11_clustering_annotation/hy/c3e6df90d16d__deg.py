import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import gseapy as gp
from gseapy.plot import barplot, dotplot
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import warnings
warnings.filterwarnings('ignore')

# 设置绘图样式
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

def perform_enrichment_analysis(deg_results, organism='human', top_n=100):
    """
    对差异表达基因进行GO和KEGG富集分析[6,8](@ref)
    """
    print("\n开始进行GO和KEGG富集分析...")
    
    # 提取显著差异基因（按p值排序）
    sig_genes = deg_results[deg_results['significant']].sort_values('padj')
    
    if len(sig_genes) == 0:
        print("警告: 没有显著的差异表达基因，无法进行富集分析。")
        return None
    
    # 使用最显著的top_n个基因
    gene_list = sig_genes['gene'].head(top_n).tolist()
    print(f"使用前 {len(gene_list)} 个显著基因进行富集分析。")
    
    enrichment_results = {}
    
    try:
        # 定义要分析的基因集[8,9](@ref)
        gene_sets = {
            'GO': ['GO_Biological_Process_2023', 
                   'GO_Cellular_Component_2023', 
                   'GO_Molecular_Function_2023'],
            'KEGG': ['KEGG_2021_Human']
        }
        
        # 执行富集分析[6,9](@ref)
        enr = gp.enrichr(gene_list=gene_list,
                        gene_sets=gene_sets['GO'] + gene_sets['KEGG'],
                        organism=organism,
                        outdir=None,
                        format='png')
        
        enrichment_results['combined'] = enr
        print("富集分析完成！")
        return enrichment_results
        
    except Exception as e:
        print(f"富集分析出错: {e}")
        return None

def plot_combined_bubble_plot(enrichment_results, top_terms=10, figsize=(12, 8)):
    """
    绘制GO和KEGG结果的组合气泡图[1](@ref)
    """
    if enrichment_results is None:
        return None
    
    enr = enrichment_results['combined']
    if enr is None or enr.res2d is None or len(enr.res2d) == 0:
        print("没有富集结果可绘制")
        return None
    
    # 准备数据
    df = enr.res2d.head(top_terms * 2)  # 取前2倍项进行筛选
    
    if len(df) == 0:
        return None
    
    # 计算GeneRatio的数值
    def calculate_gene_ratio(overlap):
        try:
            parts = overlap.split('/')
            return int(parts[0]) / int(parts[1])
        except:
            return 0.0
    
    df['GeneRatio'] = df['Overlap'].apply(calculate_gene_ratio)
    df['-log10_padj'] = -np.log10(df['Adjusted P-value'])
    
    # 创建图形
    fig, ax = plt.subplots(figsize=figsize)
    
    # 为不同数据库设置颜色
    colors = {'GO_Biological_Process_2023': '#E41A1C',
              'GO_Cellular_Component_2023': '#377EB8', 
              'GO_Molecular_Function_2023': '#4DAF4A',
              'KEGG_2021_Human': '#984EA3'}
    
    # 绘制气泡图
    for idx, row in df.iterrows():
        gene_set = row['Gene_set']
        color = colors.get(gene_set, '#999999')
        
        scatter = ax.scatter(row['GeneRatio'], row['-log10_padj'],
                           s=row['GeneRatio'] * 500,  # 气泡大小与GeneRatio成正比
                           c=color, alpha=0.7, edgecolors='black', linewidth=0.5)
        
        # 添加术语标签（只添加前几个以避免过于拥挤）
        if idx < 8:
            ax.annotate(row['Term'][:30] + '...' if len(row['Term']) > 30 else row['Term'],
                       (row['GeneRatio'], row['-log10_padj']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.8)
    
    ax.set_xlabel('Gene Ratio', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log10(Adjusted P-value)', fontsize=12, fontweight='bold')
    ax.set_title('GO and KEGG Enrichment Analysis', fontsize=14, fontweight='bold')
    
    # 创建图例
    legend_patches = [mpatches.Patch(color=color, label=gene_set) 
                     for gene_set, color in colors.items()]
    ax.legend(handles=legend_patches, loc='upper right')
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return fig

def plot_treemap(enrichment_results, top_terms=15, figsize=(14, 10)):
    """
    绘制富集结果的树图[3](@ref)
    """
    if enrichment_results is None:
        return None
    
    enr = enrichment_results['combined']
    if enr is None or enr.res2d is None or len(enr.res2d) == 0:
        return None
    
    # 准备数据
    df = enr.res2d.head(top_terms).copy()
    
    # 计算富集分数和转换p值
    df['Rich Factor'] = df['Overlap'].apply(lambda x: int(x.split('/')[0]) / int(x.split('/')[1]))
    df['-log10_padj'] = -np.log10(df['Adjusted P-value'])
    
    # 创建树图
    fig, ax = plt.subplots(figsize=figsize)
    
    # 简单的树图实现
    categories = df['Gene_set'].unique()
    colors = plt.cm.Set3(np.linspace(0, 1, len(categories)))
    color_map = dict(zip(categories, colors))
    
    y_pos = 0
    for category in categories:
        cat_data = df[df['Gene_set'] == category].head(5)
        
        for i, (idx, row) in enumerate(cat_data.iterrows()):
            # 矩形大小与富集程度成正比
            size = row['Rich Factor'] * 200
            rect = mpatches.Rectangle((i * 2, y_pos), 1.5, 1,
                                    facecolor=color_map[category],
                                    alpha=0.7, edgecolor='black')
            ax.add_patch(rect)
            
            # 添加文本
            ax.text(i * 2 + 0.75, y_pos + 0.5, 
                   f"{row['Term'][:15]}...\n"
                   f"p={row['Adjusted P-value']:.2e}",
                   ha='center', va='center', fontsize=8)
        
        y_pos += 2
    
    ax.set_xlim(0, len(categories) * 2)
    ax.set_ylim(0, y_pos)
    ax.set_title('Enrichment Analysis Treemap', fontsize=14, fontweight='bold')
    ax.axis('off')  # 隐藏坐标轴
    
    # 添加图例
    legend_patches = [mpatches.Patch(color=color, label=category) 
                     for category, color in color_map.items()]
    ax.legend(handles=legend_patches, loc='upper center', 
              bbox_to_anchor=(0.5, -0.05), ncol=len(categories))
    
    plt.tight_layout()
    return fig

def plot_circular_barplot(enrichment_results, top_terms_per_category=5, figsize=(12, 10)):
    """
    绘制圆形条形图展示富集结果[4](@ref)
    """
    if enrichment_results is None:
        return None
    
    enr = enrichment_results['combined']
    if enr is None or enr.res2d is None or len(enr.res2d) == 0:
        return None
    
    # 准备数据
    df = enr.res2d.copy()
    df['-log10_padj'] = -np.log10(df['Adjusted P-value'])
    
    # 按类别分组并取每个类别的前几个术语
    top_terms = []
    for category in df['Gene_set'].unique():
        cat_data = df[df['Gene_set'] == category].nsmallest(top_terms_per_category, 'Adjusted P-value')
        top_terms.append(cat_data)
    
    if len(top_terms) == 0:
        return None
    
    plot_df = pd.concat(top_terms)
    
    # 创建圆形条形图
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, polar=True)
    
    # 计算角度
    angles = np.linspace(0, 2*np.pi, len(plot_df), endpoint=False).tolist()
    angles += angles[:1]  # 完成圆圈
    
    # 准备值
    values = plot_df['-log10_padj'].tolist()
    values += values[:1]
    
    # 绘制雷达图
    ax.plot(angles, values, 'o-', linewidth=2, label='-log10(padj)')
    ax.fill(angles, values, alpha=0.25)
    
    # 添加术语标签
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels([term[:20] + '...' for term in plot_df['Term']], fontsize=8)
    
    ax.set_title('Circular Bar Plot of Enrichment Results', size=14, fontweight='bold')
    ax.grid(True)
    
    plt.tight_layout()
    return fig

def plot_interactive_bubble_plot(enrichment_results, top_terms=20):
    """
    使用plotly创建交互式气泡图[2](@ref)
    """
    if enrichment_results is None:
        return None
    
    enr = enrichment_results['combined']
    if enr is None or enr.res2d is None or len(enr.res2d) == 0:
        return None
    
    # 准备数据
    df = enr.res2d.head(top_terms).copy()
    
    # 计算GeneRatio和富集分数
    def calculate_gene_ratio(overlap):
        try:
            parts = overlap.split('/')
            return int(parts[0]) / int(parts[1])
        except:
            return 0.0
    
    df['GeneRatio'] = df['Overlap'].apply(calculate_gene_ratio)
    df['-log10_padj'] = -np.log10(df['Adjusted P-value'])
    df['Count'] = df['Overlap'].apply(lambda x: int(x.split('/')[0]))
    
    # 创建交互式气泡图
    fig = px.scatter(df, x='GeneRatio', y='-log10_padj',
                    size='Count', color='Gene_set',
                    hover_name='Term', hover_data=['Adjusted P-value', 'Odds Ratio'],
                    size_max=30,
                    title='Interactive GO and KEGG Enrichment Analysis')
    
    fig.update_layout(
        xaxis_title='Gene Ratio',
        yaxis_title='-log10(Adjusted P-value)',
        showlegend=True,
        hovermode='closest'
    )
    
    return fig

def save_enrichment_plots(enrichment_results, output_dir='./enrichment_plots/'):
    """
    保存所有富集分析图表
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("正在生成富集分析图表...")
    
    # 1. 组合气泡图
    fig1 = plot_combined_bubble_plot(enrichment_results)
    if fig1:
        fig1.savefig(f'{output_dir}/combined_bubble_plot.png', dpi=300, bbox_inches='tight')
        print("✓ 组合气泡图已保存")
    
    # 2. 树图
    fig2 = plot_treemap(enrichment_results)
    if fig2:
        fig2.savefig(f'{output_dir}/enrichment_treemap.png', dpi=300, bbox_inches='tight')
        print("✓ 树图已保存")
    
    # 3. 圆形条形图
    fig3 = plot_circular_barplot(enrichment_results)
    if fig3:
        fig3.savefig(f'{output_dir}/circular_barplot.png', dpi=300, bbox_inches='tight')
        print("✓ 圆形条形图已保存")
    
    # 4. 交互式图表
    fig4 = plot_interactive_bubble_plot(enrichment_results)
    if fig4:
        fig4.write_html(f'{output_dir}/interactive_bubble_plot.html')
        print("✓ 交互式气泡图已保存")
    
    # 5. 使用gseapy内置绘图功能[8,9](@ref)
    try:
        enr = enrichment_results['combined']
        
        # 分别绘制GO和KEGG的结果
        for gene_set in enr.results['Gene_set'].unique():
            subset_data = enr.results[enr.results['Gene_set'] == gene_set].head(10)
            
            if len(subset_data) > 0:
                # 气泡图
                dotplot(subset_data, title=f'{gene_set} - Dot Plot',
                       cmap='viridis_r', ofname=f'{output_dir}/{gene_set}_dotplot.png')
                
                # 条形图
                barplot(subset_data, title=f'{gene_set} - Bar Plot',
                       ofname=f'{output_dir}/{gene_set}_barplot.png')
        
        print("✓ gseapy内置图表已保存")
    except Exception as e:
        print(f"gseapy内置图表保存失败: {e}")
    
    print(f"所有图表已保存至: {output_dir}")

def main():
    """
    主函数 - 集成差异表达分析和富集分析可视化
    """
    # 这里假设您已经有差异表达分析的结果
    # 如果是单独运行这个脚本，需要先加载差异表达结果
    
    print("=" * 60)
    print("GO和KEGG富集分析可视化脚本")
    print("=" * 60)
    
    # 示例：加载已有的差异表达分析结果
    try:
        deg_results = pd.read_csv('differential_expression_results.csv')
        print("差异表达结果加载成功！")
    except:
        print("请先运行差异表达分析，或确保differential_expression_results.csv文件存在")
        return
    
    # 执行富集分析
    enrichment_results = perform_enrichment_analysis(deg_results)
    
    if enrichment_results is None:
        print("富集分析失败，无法生成图表")
        return
    
    # 保存所有图表
    save_enrichment_plots(enrichment_results)
    
    # 显示图表
    plt.show()
    
    print("\n" + "=" * 60)
    print("富集分析可视化完成！")
    print("生成的文件包括:")
    print("- 组合气泡图 (combined_bubble_plot.png)")
    print("- 富集树图 (enrichment_treemap.png)") 
    print("- 圆形条形图 (circular_barplot.png)")
    print("- 交互式气泡图 (interactive_bubble_plot.html)")
    print("- 各数据库单独的气泡图和条形图")
    print("=" * 60)

if __name__ == "__main__":
    main()
