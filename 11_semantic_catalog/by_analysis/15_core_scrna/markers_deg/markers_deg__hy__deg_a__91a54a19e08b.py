import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import os
import warnings
from matplotlib import rcParams
warnings.filterwarnings('ignore')

# ---------------------- 全局配置 ----------------------
output_dir = "/home/hy/Project/FWSP20251095/Redo/DEG/Nerve_cells"
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['pdf.fonttype'] = 42
os.makedirs(output_dir, exist_ok=True)

# ---------------------- 核心函数 ----------------------
def perform_dea(adata):
    """执行差异表达分析"""
    sc.tl.rank_genes_groups(
        adata, 
        groupby='group', 
        method='wilcoxon',
        groups=['disease'],
        reference='normal',
        key_added='dea'
    )
    deg_df = sc.get.rank_genes_groups_df(adata, group='disease', key='dea')
    deg_df['regulation'] = np.where(
        deg_df['logfoldchanges'] > 0, 
        'Upregulated', 
        'Downregulated'
    )
    return deg_df

def enhanced_volcano(deg_df, lfc_thresh=0.5, pval_thresh=0.05):
    """绘制增强型火山图"""
    plt.figure(figsize=(10, 8))
    
    # 显著性分类
    deg_df['significance'] = np.select(
        [
            (deg_df['pvals_adj'] < pval_thresh) & (deg_df['logfoldchanges'] > lfc_thresh),
            (deg_df['pvals_adj'] < pval_thresh) & (deg_df['logfoldchanges'] < -lfc_thresh),
            True
        ],
        ['Upregulated', 'Downregulated', 'Not significant']
    )
    
    # 颜色和样式
    palette = {'Upregulated':'#E64B35', 'Downregulated':'#3182BD', 'Not significant':'#D9D9D9'}
    
    sns.scatterplot(
        data=deg_df,
        x='logfoldchanges',
        y=-np.log10(deg_df['pvals_adj']),
        hue='significance',
        palette=palette,
        s=25,
        alpha=0.7,
        edgecolor='none'
    )
    
    # 标注top基因
    for direction in ['Upregulated', 'Downregulated']:
        top_genes = deg_df.query('significance == @direction').nsmallest(3, 'pvals_adj')
        for _, row in top_genes.iterrows():
            plt.text(
                row['logfoldchanges'], 
                -np.log10(row['pvals_adj'])+0.1, 
                row['names'],
                fontsize=8,
                ha='center',
                va='bottom'
            )
    
    plt.axhline(-np.log10(pval_thresh), color='grey', ls='--', lw=1)
    plt.axvline(lfc_thresh, color='grey', ls='--', lw=1)
    plt.axvline(-lfc_thresh, color='grey', ls='--', lw=1)
    
    plt.xlabel('log2(Fold Change)', fontsize=12)
    plt.ylabel('-log10(Adj p-value)', fontsize=12)
    plt.title('Differential Expression Volcano Plot', fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig(os.path.join(output_dir, 'Volcano_Plot.pdf'), bbox_inches='tight')
    plt.close()

def optimized_heatmap(adata, deg_df, n_genes=20):
    """生成优化热图"""
    # 筛选有效基因
    valid_genes = [g for g in deg_df['names'] if g in adata.var_names]
    if not valid_genes:
        print("警告：无有效基因可供绘制热图")
        return
    
    # 筛选显著基因
    sig_genes = deg_df[
        (deg_df['pvals_adj'] < 0.05) & 
        (abs(deg_df['logfoldchanges']) > 0.5) &
        (deg_df['names'].isin(valid_genes))
    ]
    
    if sig_genes.empty:
        print("警告：无显著差异基因")
        return
    
    # 分上下调选取top基因
    up_genes = sig_genes[sig_genes['regulation'] == 'Upregulated'].nlargest(n_genes, 'logfoldchanges')
    down_genes = sig_genes[sig_genes['regulation'] == 'Downregulated'].nsmallest(n_genes, 'logfoldchanges')
    selected_genes = pd.concat([up_genes, down_genes])['names'].tolist()
    
    # 二次过滤确保基因存在
    selected_genes = [g for g in selected_genes if g in adata.var_names]
    if not selected_genes:
        print("警告：筛选后无可用基因")
        return
    
    # 创建表达矩阵
    expr_df = pd.DataFrame(
        adata[:, selected_genes].X.toarray(),
        columns=selected_genes,
        index=adata.obs['sample']
    ).groupby('sample').mean().T
    
    # 绘制热图
    plt.figure(figsize=(12, 8))
    sns.heatmap(
        expr_df,
        cmap='viridis',
        annot=False,
        fmt=".2f",
        linewidths=0.5,
        yticklabels=True
    )
    plt.title('Top Differentially Expressed Genes', fontsize=14)
    plt.xlabel('Samples', fontsize=12)
    plt.ylabel('Genes', fontsize=12)
    plt.savefig(os.path.join(output_dir, 'Heatmap.pdf'), bbox_inches='tight')
    plt.close()

def run_gsea(deg_df):
    """执行GSEA富集分析"""
    # 生成排序文件
    rank_df = deg_df[['names', 'logfoldchanges']].sort_values('logfoldchanges', ascending=False)
    rank_df.to_csv(os.path.join(output_dir, 'gsea_rank.rnk'), sep='\t', index=False, header=False)
    
    # 运行分析
    gp.prerank(
        rnk=os.path.join(output_dir, 'gsea_rank.rnk'),
        gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2023'],
        processes=4,
        permutation_num=1000,
        outdir=os.path.join(output_dir, 'GSEA_Results'),
        format='png'
    )

# ---------------------- 主流程 ----------------------
def main(adata):
    # 数据预处理
    if 'log1p' not in adata.uns_keys():
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    
    # 差异分析
    deg_df = perform_dea(adata)
    deg_df.to_csv(os.path.join(output_dir, 'DEG_Results.csv'), index=False)
    
    # 可视化
    enhanced_volcano(deg_df)
    optimized_heatmap(adata, deg_df)
    
    # GSEA分析
    run_gsea(deg_df)

if __name__ == "__main__":
    adata = sc.read_h5ad("Nerve_cells.h5ad")
    os.makedirs(os.path.join(output_dir, 'GSEA_Results'), exist_ok=True)
    main(adata)
