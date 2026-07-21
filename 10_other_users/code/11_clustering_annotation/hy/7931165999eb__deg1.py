import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from scipy import stats
import os
import warnings
warnings.filterwarnings('ignore')

def create_directories():
    """创建保存结果的目录"""
    os.makedirs('diff_genes', exist_ok=True)
    os.makedirs('enrichment_results', exist_ok=True)
    os.makedirs('visualization', exist_ok=True)

def load_and_preprocess_data(file_path):
    """加载和预处理数据"""
    print("Loading h5ad data...")
    adata = sc.read_h5ad(file_path)
    
    # 基本质量控制
    if 'highly_variable' not in adata.var.columns:
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        adata = adata[:, adata.var.highly_variable]
        sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(adata, max_value=10)
    
    print(f"Data shape: {adata.shape}")
    print(f"Grouping variables: g1: {adata.obs['g1'].unique()}, g2: {adata.obs['g2'].unique()}")
    print(f"Cell types: {adata.obs['celltype'].unique()}")
    
    return adata

def perform_differential_expression_analysis(adata, group_var, celltype):
    """执行差异表达分析"""
    print(f"\nAnalyzing {group_var} in {celltype}...")
    
    # 提取特定细胞类型的子集
    adata_subset = adata[adata.obs['celltype'] == celltype].copy()
    
    if len(adata_subset.obs[group_var].unique()) < 2:
        print(f"Only one group found in {group_var} for {celltype}, skipping...")
        return None
    
    # 使用Wilcoxon秩和检验进行差异表达分析
    sc.tl.rank_genes_groups(adata_subset, groupby=group_var, 
                           method='wilcoxon', use_raw=False)
    
    # 提取差异表达结果
    result_dict = {}
    groups = adata_subset.uns['rank_genes_groups']['names'].dtype.names
    
    for group in groups:
        deg_df = pd.DataFrame({
            'genes': adata_subset.uns['rank_genes_groups']['names'][group],
            'scores': adata_subset.uns['rank_genes_groups']['scores'][group],
            'pvals': adata_subset.uns['rank_genes_groups']['pvals'][group],
            'pvals_adj': adata_subset.uns['rank_genes_groups']['pvals_adj'][group],
            'logfoldchanges': adata_subset.uns['rank_genes_groups']['logfoldchanges'][group]
        })
        result_dict[group] = deg_df
    
    return result_dict

def filter_significant_genes(deg_df, pval_threshold=0.05, logfc_threshold=0.5):
    """筛选显著差异表达基因"""
    significant_genes = deg_df[
        (deg_df['pvals_adj'] < pval_threshold) & 
        (abs(deg_df['logfoldchanges']) > logfc_threshold)
    ]
    return significant_genes.sort_values('logfoldchanges', ascending=False)

def perform_enrichment_analysis(gene_list, gene_sets, organism='human'):
    """执行富集分析"""
    if len(gene_list) < 5:
        print("Too few genes for enrichment analysis, skipping...")
        return None
    
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism,
            outdir=None,
            cutoff=0.05
        )
        return enr
    except Exception as e:
        print(f"Enrichment analysis error: {e}")
        return None

def plot_volcano(deg_df, group_name, celltype, group_var, pval_thresh=0.05, logfc_thresh=0.5):
    """绘制火山图"""
    plt.figure(figsize=(10, 8))
    
    # 标记显著基因
    deg_df['significance'] = 'Not significant'
    deg_df.loc[(deg_df['pvals_adj'] < pval_thresh) & 
              (abs(deg_df['logfoldchanges']) > logfc_thresh), 'significance'] = 'Significant'
    deg_df.loc[(deg_df['pvals_adj'] < pval_thresh) & 
              (deg_df['logfoldchanges'] > logfc_thresh), 'significance'] = 'Up-regulated'
    deg_df.loc[(deg_df['pvals_adj'] < pval_thresh) & 
              (deg_df['logfoldchanges'] < -logfc_thresh), 'significance'] = 'Down-regulated'
    
    # 绘制散点
    colors = {'Not significant': 'gray', 'Significant': 'blue', 
              'Up-regulated': 'red', 'Down-regulated': 'green'}
    
    for significance, color in colors.items():
        subset = deg_df[deg_df['significance'] == significance]
        plt.scatter(subset['logfoldchanges'], -np.log10(subset['pvals_adj']), 
                   c=color, alpha=0.6, s=20, label=significance)
    
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10(Adjusted P-value)')
    plt.title(f'Volcano Plot: {group_name} in {celltype} ({group_var})')
    plt.axhline(-np.log10(pval_thresh), linestyle='--', color='red', alpha=0.8)
    plt.axvline(logfc_thresh, linestyle='--', color='red', alpha=0.8)
    plt.axvline(-logfc_thresh, linestyle='--', color='red', alpha=0.8)
    plt.legend()
    plt.tight_layout()
    
    filename = f"visualization/volcano_{group_var}_{celltype}_{group_name}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def plot_enrichment_bubble(enr_results, title, top_n=10):
    """绘制富集分析气泡图"""
    if enr_results is None or len(enr_results.results) == 0:
        return
    
    # 获取前N个最显著的富集项
    top_terms = enr_results.results.head(top_n)
    
    if len(top_terms) == 0:
        return
    
    plt.figure(figsize=(12, 8))
    
    # 创建气泡图
    scatter = plt.scatter(
        -np.log10(top_terms['Adjusted P-value']),
        range(len(top_terms)),
        s=top_terms['Overlap'].str.split('/').str[0].astype(float) * 10,  # 气泡大小表示基因数量
        c=-np.log10(top_terms['Adjusted P-value']),
        cmap='viridis',
        alpha=0.7
    )
    
    plt.yticks(range(len(top_terms)), top_terms['Term'])
    plt.xlabel('-Log10(Adjusted P-value)')
    plt.title(f'Enrichment Analysis: {title}')
    plt.colorbar(scatter, label='-Log10(Adjusted P-value)')
    plt.tight_layout()
    
    filename = f"visualization/enrichment_bubble_{title.replace(' ', '_')}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def main_analysis(adata):
    """主分析函数"""
    all_results = {}
    grouping_vars = ['g1', 'g2']
    gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2023']
    
    for group_var in grouping_vars:
        print(f"\n{'='*50}")
        print(f"Analyzing grouping variable: {group_var}")
        print(f"{'='*50}")
        
        all_results[group_var] = {}
        
        for celltype in adata.obs['celltype'].unique():
            print(f"\nProcessing cell type: {celltype}")
            
            # 差异表达分析
            deg_results = perform_differential_expression_analysis(adata, group_var, celltype)
            
            if deg_results is None:
                continue
                
            all_results[group_var][celltype] = {}
            
            for group_name, deg_df in deg_results.items():
                print(f"  Group: {group_name}")
                
                # 筛选显著差异基因
                significant_genes = filter_significant_genes(deg_df)
                
                if len(significant_genes) == 0:
                    print(f"    No significant genes found in {group_name}")
                    continue
                
                # 保存差异表达结果
                filename = f"diff_genes/{group_var}_{celltype}_{group_name}_deg.csv"
                significant_genes.to_csv(filename, index=False)
                
                # 绘制火山图
                plot_volcano(deg_df, group_name, celltype, group_var)
                
                # 提取基因列表进行富集分析
                gene_list = significant_genes['genes'].head(100).tolist()
                
                # 执行富集分析
                enr_results = perform_enrichment_analysis(gene_list, gene_sets)
                
                if enr_results is not None and len(enr_results.results) > 0:
                    # 保存富集分析结果
                    enr_filename = f"enrichment_results/{group_var}_{celltype}_{group_name}_enrichment.csv"
                    enr_results.results.to_csv(enr_filename, index=False)
                    
                    # 绘制富集分析图
                    plot_title = f"{group_var}_{celltype}_{group_name}"
                    plot_enrichment_bubble(enr_results, plot_title)
                    
                    all_results[group_var][celltype][group_name] = {
                        'deg': significant_genes,
                        'enrichment': enr_results.results
                    }
                
                print(f"    Found {len(significant_genes)} significant genes")
    
    return all_results

def generate_summary_report(all_results):
    """生成总结报告"""
    print("\n" + "="*60)
    print("SUMMARY REPORT")
    print("="*60)
    
    for group_var in all_results.keys():
        print(f"\nGrouping variable: {group_var}")
        print("-" * 40)
        
        for celltype in all_results[group_var].keys():
            print(f"  Cell type: {celltype}")
            
            for group_name, results in all_results[group_var][celltype].items():
                n_genes = len(results['deg'])
                n_pathways = len(results['enrichment'])
                
                print(f"    {group_name}: {n_genes} DEGs, {n_pathways} enriched pathways")
                
                # 显示最显著的通路
                if n_pathways > 0:
                    top_pathway = results['enrichment'].iloc[0]
                    print(f"      Top pathway: {top_pathway['Term']} (p={top_pathway['Adjusted P-value']:.2e})")

# 主执行流程
if __name__ == "__main__":
    # 设置文件路径
    h5ad_file = "epi_with_g1_g2.h5ad"  # 修改为你的h5ad文件路径
    
    # 创建目录
    create_directories()
    
    try:
        # 加载数据
        adata = load_and_preprocess_data(h5ad_file)
        
        # 执行主分析
        results = main_analysis(adata)
        
        # 生成总结报告
        generate_summary_report(results)
        
        print("\nAnalysis completed successfully!")
        print("Results saved in:")
        print("  - diff_genes/ : Differential expression results")
        print("  - enrichment_results/ : Enrichment analysis results")
        print("  - visualization/ : Visualization plots")
        
    except Exception as e:
        print(f"Error occurred: {e}")
        import traceback
        traceback.print_exc()
