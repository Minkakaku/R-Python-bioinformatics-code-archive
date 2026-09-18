import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import os
import warnings
warnings.filterwarnings('ignore')

# 设置绘图样式
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

def load_and_validate_data(file_path):
    """
    加载h5ad文件并验证必要的列是否存在
    """
    print("正在加载数据...")
    adata = sc.read_h5ad(file_path)
    print(f"数据加载成功！细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")

    # 检查必需的列
    required_cols = ['TIME', 'celltype']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        print(f"错误: 缺少以下列: {missing_cols}")
        print(f"可用列: {list(adata.obs.columns)}")
        return None
    
    # 检查TIME分组情况
    treatment_groups = adata.obs['TIME'].unique()
    if len(treatment_groups) < 2:
        print(f"错误: 需要至少2个治疗组，但找到 {len(treatment_groups)} 个: {treatment_groups}")
        return None
    
    print(f"治疗分组: {treatment_groups}")
    print(f"细胞类型: {adata.obs['celltype'].unique()}")
    
    return adata

def create_celltype_directory(celltype, base_dir='./celltype_deg_results'):
    """
    为每个细胞类型创建独立的输出目录
    """
    # 清理细胞类型名称，移除可能无效的目录名字符
    safe_celltype = "".join(c for c in celltype if c.isalnum() or c in (' ', '-', '_')).rstrip()
    celltype_dir = os.path.join(base_dir, safe_celltype)
    
    # 创建目录
    os.makedirs(celltype_dir, exist_ok=True)
    print(f"创建目录: {celltype_dir}")
    
    return celltype_dir

def perform_de_analysis_per_celltype(adata, celltype, treatment_col='TIME', method='wilcoxon'):
    """
    对单个细胞类型进行差异表达分析
    """
    print(f"\n正在分析细胞类型: {celltype}")
    
    # 提取当前细胞类型的子集
    celltype_mask = adata.obs['celltype'] == celltype
    bdata = adata[celltype_mask].copy()
    
    print(f"细胞类型 {celltype} 的细胞数: {bdata.n_obs}")
    
    # 检查治疗组分布
    treatment_counts = bdata.obs[treatment_col].value_counts()
    print(f"治疗组分布: {treatment_counts.to_dict()}")
    
    if len(treatment_counts) < 2:
        print(f"警告: 细胞类型 {celltype} 中只有 {len(treatment_counts)} 个治疗组，跳过分析")
        return None, None, None
    
    # 获取治疗组
    groups = bdata.obs[treatment_col].unique()
    group1, group2 = groups[0], groups[1]
    
    print(f"比较组: {group1} vs {group2}")
    
    try:
        # 使用scanpy进行差异表达分析
        sc.tl.rank_genes_groups(bdata, groupby=treatment_col, method=method,
                               groups=[group1], reference=group2,
                               use_raw=True, n_genes=bdata.raw.n_vars if bdata.raw is not None else bdata.n_vars)
        
        # 提取差异表达结果
        deg_results = sc.get.rank_genes_groups_df(bdata, group=group1)
        
        # 重命名列
        deg_results = deg_results.rename(columns={
            'names': 'gene',
            'scores': 'score',
            'logfoldchanges': 'log2FC',
            'pvals': 'pvalue',
            'pvals_adj': 'padj'
        })
        
        # 添加显著性标记
        deg_results['-log10_padj'] = -np.log10(deg_results['padj'])
        deg_results['significant'] = (deg_results['padj'] < 0.05) & (np.abs(deg_results['log2FC']) > 1)
        deg_results['regulation'] = np.where(
            deg_results['significant'], 
            np.where(deg_results['log2FC'] > 0, 'up', 'down'), 
            'non-significant'
        )
        
        # 统计信息
        n_sig = deg_results['significant'].sum()
        n_up = (deg_results['regulation'] == 'up').sum()
        n_down = (deg_results['regulation'] == 'down').sum()
        
        print(f"细胞类型 {celltype} 的差异表达分析完成！")
        print(f"显著差异基因数: {n_sig} (上调: {n_up}, 下调: {n_down})")
        
        return deg_results, group1, group2
        
    except Exception as e:
        print(f"细胞类型 {celltype} 的差异表达分析失败: {e}")
        return None, None, None

def plot_volcano_per_celltype(deg_results, celltype, group1, group2, output_dir, figsize=(10, 6)):
    """
    为每个细胞类型绘制火山图
    """
    if deg_results is None or len(deg_results) == 0:
        return None
    
    print(f"正在为细胞类型 {celltype} 绘制火山图...")
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # 为不同类别的点分配颜色
    colors = {'up': '#E64B35', 'down': '#3182BD', 'non-significant': '#969696'}
    
    for reg in ['non-significant', 'down', 'up']:
        data = deg_results[deg_results['regulation'] == reg]
        if len(data) > 0:
            ax.scatter(data['log2FC'], data['-log10_padj'], 
                      c=colors[reg], label=reg, alpha=0.6, s=15, edgecolors='none')
    
    # 添加显著性阈值线
    ax.axhline(-np.log10(0.05), color='grey', linestyle='--', linewidth=1, alpha=0.8)
    ax.axvline(-1, color='grey', linestyle='--', linewidth=1, alpha=0.8)
    ax.axvline(1, color='grey', linestyle='--', linewidth=1, alpha=0.8)
    
    # 标注最显著的几个基因
    sig_genes = deg_results[deg_results['significant']]
    if len(sig_genes) > 0:
        top_genes = sig_genes.nlargest(5, '-log10_padj')
        for _, gene in top_genes.iterrows():
            ax.annotate(gene['gene'], (gene['log2FC'], gene['-log10_padj']), 
                       xytext=(5, 5), textcoords='offset points', fontsize=8, alpha=0.8)
    
    ax.set_xlabel('log2(Fold Change)', fontweight='bold')
    ax.set_ylabel('-log10(Adjusted P-value)', fontweight='bold')
    ax.set_title(f'Volcano Plot: {celltype} - {group1} vs {group2}', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # 保存图片
    volcanopath = os.path.join(output_dir, f'{celltype}_volcano_plot.png')
    plt.savefig(volcanopath, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"火山图已保存: {volcanopath}")
    
    return fig

def perform_enrichment_analysis_per_celltype(deg_results, celltype, organism='Human', top_n=100):
    """
    对单个细胞类型的差异基因进行GO和KEGG富集分析[1,4](@ref)
    """
    if deg_results is None:
        return None
    
    # 提取显著差异基因
    sig_genes = deg_results[deg_results['significant']]
    
    if len(sig_genes) == 0:
        print(f"细胞类型 {celltype} 无显著差异基因，跳过富集分析")
        return None
    
    gene_list = sig_genes['gene'].head(top_n).tolist()
    print(f"细胞类型 {celltype} 使用前 {len(gene_list)} 个显著基因进行富集分析")
    
    try:
        # 定义GO和KEGG基因集[1,4](@ref)
        gene_sets = [
            'GO_Biological_Process_2023',
            'GO_Cellular_Component_2023', 
            'GO_Molecular_Function_2023',
            'KEGG_2021_Human'
        ]
        
        enr = gp.enrichr(gene_list=gene_list,
                        gene_sets=gene_sets,
                        organism=organism,
                        outdir=None,
                        cutoff=0.05)
        
        print(f"细胞类型 {celltype} 的富集分析完成")
        return enr
        
    except Exception as e:
        print(f"细胞类型 {celltype} 的富集分析失败: {e}")
        return None

def plot_enrichment_barplot(enrichment_result, celltype, output_dir, top_terms=10, figsize=(12, 8)):
    """
    绘制富集分析条形图[6,8](@ref)
    """
    if enrichment_result is None or enrichment_result.results is None:
        return None
    
    try:
        # 提取结果数据
        df = enrichment_result.results
        
        # 按基因集类型分组并选择每个组的前top_terms个最显著项
        barplot_data = []
        for gene_set in df['Gene_set'].unique():
            subset = df[df['Gene_set'] == gene_set].nsmallest(top_terms, 'Adjusted P-value')
            barplot_data.append(subset)
        
        if len(barplot_data) == 0:
            return None
            
        plot_df = pd.concat(barplot_data)
        
        # 创建条形图
        fig, ax = plt.subplots(figsize=figsize)
        
        # 为不同基因集设置颜色
        palette = {'GO_Biological_Process_2023': '#1f77b4',
                  'GO_Cellular_Component_2023': '#ff7f0e',
                  'GO_Molecular_Function_2023': '#2ca02c',
                  'KEGG_2021_Human': '#d62728'}
        
        # 绘制条形图
        bars = ax.barh(range(len(plot_df)), plot_df['Odds Ratio'], 
                      color=[palette.get(gs, '#7f7f7f') for gs in plot_df['Gene_set']],
                      alpha=0.7)
        
        ax.set_yticks(range(len(plot_df)))
        ax.set_yticklabels(plot_df['Term'])
        ax.set_xlabel('Odds Ratio', fontweight='bold')
        ax.set_ylabel('Enriched Terms', fontweight='bold')
        ax.set_title(f'Enrichment Analysis Bar Plot: {celltype}', fontweight='bold')
        
        # 添加网格
        ax.grid(True, alpha=0.3, axis='x')
        
        # 创建图例
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=color, label=gene_set) 
                          for gene_set, color in palette.items()]
        ax.legend(handles=legend_elements, loc='lower right')
        
        plt.tight_layout()
        
        # 保存条形图
        barpath = os.path.join(output_dir, f'{celltype}_enrichment_barplot.png')
        plt.savefig(barpath, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"富集条形图已保存: {barpath}")
        
        return fig
        
    except Exception as e:
        print(f"细胞类型 {celltype} 的条形图绘制失败: {e}")
        return None

def plot_enrichment_dotplot(enrichment_result, celltype, output_dir, top_terms=10, figsize=(12, 8)):
    """
    绘制富集分析点图[6,8](@ref)
    """
    if enrichment_result is None or enrichment_result.results is None:
        return None
    
    try:
        # 使用gseapy的dotplot函数[6](@ref)
        fig = gp.dotplot(enrichment_result.results,
                        title=f'{celltype} - Enrichment Analysis Dot Plot',
                        cmap='viridis_r',
                        size=10,
                        top_term=top_terms)
        
        # 调整图形大小
        fig.set_size_inches(figsize)
        plt.tight_layout()
        
        # 保存点图
        dotpath = os.path.join(output_dir, f'{celltype}_enrichment_dotplot.png')
        plt.savefig(dotpath, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"富集点图已保存: {dotpath}")
        
        return fig
        
    except Exception as e:
        print(f"细胞类型 {celltype} 的点图绘制失败: {e}")
        # 如果gseapy的dotplot失败，使用自定义点图
        return plot_custom_dotplot(enrichment_result, celltype, output_dir, top_terms, figsize)

def plot_custom_dotplot(enrichment_result, celltype, output_dir, top_terms=10, figsize=(12, 8)):
    """
    自定义富集分析点图[8](@ref)
    """
    try:
        # 提取结果数据
        df = enrichment_result.results
        
        # 选择最显著的前top_terms个项
        plot_df = df.nsmallest(top_terms, 'Adjusted P-value')
        
        if len(plot_df) == 0:
            return None
        
        # 创建图形
        fig, ax = plt.subplots(figsize=figsize)
        
        # 计算基因比例
        def calculate_gene_ratio(overlap):
            try:
                parts = overlap.split('/')
                return int(parts[0]) / int(parts[1])
            except:
                return 0.0
        
        plot_df['GeneRatio'] = plot_df['Overlap'].apply(calculate_gene_ratio)
        plot_df['-log10_padj'] = -np.log10(plot_df['Adjusted P-value'])
        
        # 创建散点图（点图）
        scatter = ax.scatter(plot_df['GeneRatio'], plot_df['-log10_padj'],
                           s=plot_df['Odds Ratio'] * 20,  # 点大小与Odds Ratio成正比
                           c=plot_df['-log10_padj'], cmap='viridis_r', alpha=0.7)
        
        # 添加颜色条
        cbar = plt.colorbar(scatter)
        cbar.set_label('-log10(Adjusted P-value)')
        
        # 添加标签
        for i, term in enumerate(plot_df['Term']):
            ax.annotate(term[:40] + '...' if len(term) > 40 else term,
                       (plot_df['GeneRatio'].iloc[i], plot_df['-log10_padj'].iloc[i]),
                       xytext=(5, 5), textcoords='offset points', fontsize=8)
        
        ax.set_xlabel('Gene Ratio', fontweight='bold')
        ax.set_ylabel('-log10(Adjusted P-value)', fontweight='bold')
        ax.set_title(f'Enrichment Analysis Dot Plot: {celltype}', fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # 保存点图
        dotpath = os.path.join(output_dir, f'{celltype}_enrichment_dotplot.png')
        plt.savefig(dotpath, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"富集点图已保存: {dotpath}")
        
        return fig
        
    except Exception as e:
        print(f"细胞类型 {celltype} 的自定义点图绘制失败: {e}")
        return None

def save_results_per_celltype(deg_results, enrichment_result, celltype, output_dir):
    """
    保存单个细胞类型的结果文件
    """
    # 保存差异表达结果
    if deg_results is not None:
        deg_path = os.path.join(output_dir, f'{celltype}_differential_expression.csv')
        deg_results.to_csv(deg_path, index=False)
        print(f"差异表达结果已保存: {deg_path}")
    
    # 保存富集分析结果
    if enrichment_result is not None and enrichment_result.results is not None:
        enrich_path = os.path.join(output_dir, f'{celltype}_enrichment_results.csv')
        enrichment_result.results.to_csv(enrich_path, index=False)
        print(f"富集分析结果已保存: {enrich_path}")
        
        # 保存上下调基因列表
        if deg_results is not None:
            up_genes = deg_results[deg_results['regulation'] == 'up']['gene'].tolist()
            down_genes = deg_results[deg_results['regulation'] == 'down']['gene'].tolist()
            
            with open(os.path.join(output_dir, f'{celltype}_upregulated_genes.txt'), 'w') as f:
                f.write('\n'.join(up_genes))
            with open(os.path.join(output_dir, f'{celltype}_downregulated_genes.txt'), 'w') as f:
                f.write('\n'.join(down_genes))

def analyze_all_celltypes(adata, base_output_dir='./celltype_deg_results'):
    """
    分析所有细胞类型的差异表达和富集分析
    """
    print("开始按细胞类型进行差异表达分析...")
    print("=" * 60)
    
    # 获取所有细胞类型
    celltypes = adata.obs['celltype'].unique()
    print(f"发现 {len(celltypes)} 个细胞类型: {list(celltypes)}")
    
    # 存储总体结果摘要
    summary_data = []
    
    for celltype in celltypes:
        print("\n" + "=" * 50)
        
        # 为当前细胞类型创建输出目录
        celltype_dir = create_celltype_directory(celltype, base_output_dir)
        
        try:
            # 1. 差异表达分析
            deg_results, group1, group2 = perform_de_analysis_per_celltype(adata, celltype)
            
            if deg_results is None:
                print(f"跳过细胞类型 {celltype}")
                continue
            
            # 2. 绘制火山图
            volcano_fig = plot_volcano_per_celltype(deg_results, celltype, group1, group2, celltype_dir)
            
            # 3. 富集分析（GO和KEGG）
            enrichment_result = perform_enrichment_analysis_per_celltype(deg_results, celltype)
            
            # 4. 绘制富集分析图（条形图和点图）
            if enrichment_result is not None:
                # 绘制条形图
                bar_fig = plot_enrichment_barplot(enrichment_result, celltype, celltype_dir)
                
                # 绘制点图
                dot_fig = plot_enrichment_dotplot(enrichment_result, celltype, celltype_dir)
            
            # 5. 保存结果
            save_results_per_celltype(deg_results, enrichment_result, celltype, celltype_dir)
            
            # 收集摘要信息
            n_sig = deg_results['significant'].sum()
            n_up = (deg_results['regulation'] == 'up').sum()
            n_down = (deg_results['regulation'] == 'down').sum()
            
            summary_data.append({
                'CellType': celltype,
                'Total_Cells': adata[adata.obs['celltype'] == celltype].n_obs,
                'Treatment_Groups': adata[adata.obs['celltype'] == celltype].obs['TIME'].unique().tolist(),
                'Significant_Genes': n_sig,
                'Upregulated': n_up,
                'Downregulated': n_down
            })
            
            # 关闭图形以释放内存
            plt.close('all')
            
            print(f"细胞类型 {celltype} 分析完成！")
            
        except Exception as e:
            print(f"细胞类型 {celltype} 分析过程中出错: {e}")
            import traceback
            traceback.print_exc()
    
    # 保存总体摘要
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_path = os.path.join(base_output_dir, 'analysis_summary.csv')
        summary_df.to_csv(summary_path, index=False)
        print(f"\n分析摘要已保存: {summary_path}")
        
        # 打印摘要
        print("\n" + "=" * 60)
        print("分析摘要:")
        print(summary_df.to_string(index=False))
    
    return summary_data

def main():
    """
    主函数
    """
    file_path = 'extracted_specific_patients.h5ad'  # 请修改为您的文件路径
    base_output_dir = './celltype_deg_results'
    
    try:
        # 1. 加载和验证数据
        adata = load_and_validate_data(file_path)
        if adata is None:
            return
        
        # 2. 按细胞类型分析差异表达和功能富集
        summary = analyze_all_celltypes(adata, base_output_dir)
        
        print("\n" + "=" * 60)
        print("所有细胞类型的分析完成！")
        print(f"结果保存在: {base_output_dir}")
        print("每个细胞类型有独立的文件夹，包含:")
        print("- 差异表达基因表格 (CSV格式)")
        print("- 火山图 (PNG格式)") 
        print("- GO/KEGG富集分析结果 (CSV格式)")
        print("- 富集分析条形图 (Bar Plot, PNG格式)")
        print("- 富集分析点图 (Dot Plot, PNG格式)")
        print("- 上下调基因列表 (TXT格式)")
        print("=" * 60)
        
    except Exception as e:
        print(f"分析过程中出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
