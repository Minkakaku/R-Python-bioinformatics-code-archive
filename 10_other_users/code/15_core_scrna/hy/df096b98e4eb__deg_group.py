import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from gseapy.plot import barplot, dotplot
import os
import warnings
warnings.filterwarnings('ignore')

# 创建输出目录
output_dir = "DEG_group"
os.makedirs(output_dir, exist_ok=True)

# 1. 差异表达分析（修复版）
def perform_dea(adata):
    # 确保分组信息正确
    adata.obs['group'] = adata.obs['group'].astype('category')
    
    # 执行差异分析
    sc.tl.rank_genes_groups(
        adata, 
        groupby='group', 
        method='wilcoxon',
        groups=['disease'],
        reference='normal',
        key_added='dea'
    )
    
    # 使用Scanpy内置方法提取结果
    deg_df = sc.get.rank_genes_groups_df(adata, group='disease', key='dea')
    return deg_df

# 2. 保存差异基因结果（修复版）
def save_deg_results(deg_df):
    # 过滤显著基因 (调整后的阈值)
    sig_genes = deg_df[
        (deg_df['pvals_adj'] < 0.05) & 
        (abs(deg_df['logfoldchanges']) > 0.5)  # 根据数据特性调整阈值
    ].sort_values('pvals_adj')
    
    # 保存结果
    sig_genes.to_csv(os.path.join(output_dir, "DEG_results.csv"), index=False)
    return sig_genes

# 3. 富集分析（保持不变）
def perform_enrichment(sig_genes):
    gene_list = sig_genes['names'].tolist()
    
    # GO富集分析
    go_results = gp.enrichr(
        gene_list=gene_list,
        gene_sets=['GO_Biological_Process_2023'],
        organism='human',
        outdir=os.path.join(output_dir, 'GO')
    )
    
    # KEGG富集分析
    kegg_results = gp.enrichr(
        gene_list=gene_list,
        gene_sets=['KEGG_2021_Human'],
        organism='human',
        outdir=os.path.join(output_dir, 'KEGG')
    )
    
    return go_results, kegg_results

# 4. 可视化（调整字段引用）
def visualize_results(adata, sig_genes, go_res, kegg_res):
    # 火山图
    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        x=sig_genes['logfoldchanges'],
        y=-np.log10(sig_genes['pvals_adj']),
        hue=np.where(
            (sig_genes['pvals_adj'] < 0.05) & 
            (abs(sig_genes['logfoldchanges']) > 0.5), 
            'Significant', 
            'Non-significant'
        ),
        palette={'Significant': 'red', 'Non-significant': 'gray'},
        s=20
    )
    plt.title('Volcano Plot')
    plt.xlabel('log2(Fold Change)')
    plt.ylabel('-log10(Adj p-value)')
    plt.savefig(os.path.join(output_dir, 'volcano_plot.png'), dpi=300)
    plt.close()

    # 热图 (Top 20 DEGs)
    top_genes = sig_genes['names'].head(20).tolist()
    sc.pl.heatmap(
        adata, 
        var_names=top_genes,
        groupby='group',
        cmap='viridis',
        show=False,
        figsize=(12, 8),
        dendrogram=False
    )
    plt.title('Top 20 DEGs Expression')
    plt.savefig(os.path.join(output_dir, 'heatmap.png'), dpi=300)
    plt.close()

    # 富集分析可视化
    for res, name in zip([go_res, kegg_res], ['GO', 'KEGG']):
        plt.figure(figsize=(12, 8))
        dotplot(
            res.res2d,
            title=name,
            cmap='viridis',
            size=10,
            figsize=(12, 8)
        )
        plt.savefig(os.path.join(output_dir, f'{name}_enrichment.png'), dpi=300)
        plt.close()

# 主流程（调整参数传递）
def main(adata):
    # 差异分析
    deg_df = perform_dea(adata)
    
    # 保存结果
    sig_genes = save_deg_results(deg_df)
    
    if not sig_genes.empty:
        # 富集分析
        go_res, kegg_res = perform_enrichment(sig_genes)
        
        # 可视化
        visualize_results(adata, sig_genes, go_res, kegg_res)
    else:
        print("未检测到显著差异表达基因！")

if __name__ == "__main__":
    # 加载数据（请替换为你的数据路径）
    adata = sc.read_h5ad("adata_celltype.h5ad")
    
    # 运行主流程
    main(adata)

