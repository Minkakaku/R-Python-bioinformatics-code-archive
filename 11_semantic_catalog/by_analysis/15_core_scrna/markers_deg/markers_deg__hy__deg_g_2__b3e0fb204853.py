import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from gseapy.plot import dotplot
import os
import warnings
warnings.filterwarnings('ignore')

# 创建分层输出目录
output_dir = "/home/hy/Project/FWSP20251095/Redo/DEG/Adipocytes"
os.makedirs(os.path.join(output_dir, "upregulated"), exist_ok=True)
os.makedirs(os.path.join(output_dir, "downregulated"), exist_ok=True)

# 1. 差异表达分析（区分上下调）
def perform_dea(adata):
    adata.obs['group'] = adata.obs['group'].astype('category')
    
    sc.tl.rank_genes_groups(
        adata, 
        groupby='group', 
        method='wilcoxon',
        groups=['disease'],
        reference='normal',
        key_added='dea'
    )
    
    deg_df = sc.get.rank_genes_groups_df(adata, group='disease', key='dea')
    
    # 添加上下调标记
    deg_df['regulation'] = np.where(
        deg_df['logfoldchanges'] > 0, 
        'upregulated', 
        'downregulated'
    )
    return deg_df

# 2. 分开展示高低表达基因
def save_split_deg(deg_df):
    # 定义显著性阈值
    sig_genes = deg_df[
        (deg_df['pvals_adj'] < 0.05) & 
        (abs(deg_df['logfoldchanges']) > 0.5)
    ]
    
    # 分别保存上下调基因
    up_genes = sig_genes[sig_genes['regulation'] == 'upregulated']
    down_genes = sig_genes[sig_genes['regulation'] == 'downregulated']
    
    up_genes.to_csv(os.path.join(output_dir, "upregulated/upregulated_genes.csv"), index=False)
    down_genes.to_csv(os.path.join(output_dir, "downregulated/downregulated_genes.csv"), index=False)
    
    return up_genes, down_genes

# 3. 分开展示富集分析
def perform_split_enrichment(up_genes, down_genes):
    # 上调基因富集
    if not up_genes.empty:
        gp.enrichr(
            gene_list=up_genes['names'].tolist(),
            gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
            organism='human',
            outdir=os.path.join(output_dir, "upregulated/enrichment")
        )
    
    # 下调基因富集
    if not down_genes.empty:
        gp.enrichr(
            gene_list=down_genes['names'].tolist(),
            gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
            organism='human',
            outdir=os.path.join(output_dir, "downregulated/enrichment")
        )

# 4. 增强可视化
def visualize_split_results(adata, deg_df, up_genes, down_genes):
    # 火山图颜色区分上下调
    plt.figure(figsize=(10, 6))
    sns.scatterplot(
        x=deg_df['logfoldchanges'],
        y=-np.log10(deg_df['pvals_adj']),
        hue=np.select(
            [
                (deg_df['pvals_adj'] < 0.05) & (deg_df['logfoldchanges'] > 0.5),
                (deg_df['pvals_adj'] < 0.05) & (deg_df['logfoldchanges'] < -0.5)
            ],
            ['Upregulated', 'Downregulated'],
            default='Non-significant'
        ),
        palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Non-significant': 'gray'},
        s=20,
        alpha=0.6
    )
    plt.title('Volcano Plot (Red: Up, Blue: Down)')
    plt.xlabel('log2(Fold Change)')
    plt.ylabel('-log10(Adj p-value)')
    plt.savefig(os.path.join(output_dir, 'split_volcano.png'), dpi=300)
    plt.close()

    # 上下调基因热图分开展示
    for direction, genes in zip(['Upregulated', 'Downregulated'], [up_genes, down_genes]):
        if len(genes) > 0:
            top_genes = genes.nlargest(10, 'logfoldchanges')['names'] if direction == 'Upregulated' \
                       else genes.nsmallest(10, 'logfoldchanges')['names']
            
            sc.pl.matrixplot(
                adata, 
                var_names=top_genes.tolist(),
                groupby='group',
                cmap='viridis',
                show=False,
                title=f'Top 10 {direction} Genes',
                dendrogram=False
            )
            plt.savefig(os.path.join(output_dir, f'{direction.lower()}_heatmap.png'), dpi=300)
            plt.close()

# 主流程
def main(adata):
    # 差异分析
    deg_df = perform_dea(adata)
    
    # 分开展示
    up_genes, down_genes = save_split_deg(deg_df)
    
    # 富集分析
    perform_split_enrichment(up_genes, down_genes)
    
    # 可视化
    visualize_split_results(adata, deg_df, up_genes, down_genes)

if __name__ == "__main__":
    adata = sc.read_h5ad("Adipocytes.h5ad") 
    
    # 建议数据预处理
    if 'log1p' not in adata.uns_keys():
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    
    main(adata)
