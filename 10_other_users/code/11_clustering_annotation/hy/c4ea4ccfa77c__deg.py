
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import issparse

# 1. 数据预处理
def preprocess_data(adata_path):
    # 加载数据和基因列表
    adata = sc.read_h5ad(adata_path)
    target_genes = ["JAK1","JAK2","JAK3","TYK2","STAT1","STAT2","STAT3","STAT4","STAT5A","STAT5B","STAT6","SOCS1","SOCS2","SOCS3","SOCS4","SOCS5","SOCS6","IL2","IL4","IL7","IL9","IL15","IL21","TSLP","IL2RA","IL2RB","IL2RG","IL4R","IL7R","IL9R","IL15RA","IL21R","TSLPR","IL3","IL5","GMCSF","IL27A","IL31","OSM","LIF","CNTF","CT1","CT2","CLCF1","IL6","IL11","IL13","IL10","IL19","IL22","IL20","IL24","IL26","IL12A","IL12B","IL23A","EPO","GH","PRL","THPO","GCSF","LEPTIN","Upd","IFNA","IFNB","IFNE","IFNK","IFNL1","IFNL2","IFNL3","IFNW","IFNG","EGF","PDGFA","PDGFB","IL6R","IL6ST","IL11RA","IL13RA1","IL13RA2","IL27RA","IL31RA","OSMR","LIFR","CNTFR","IL10RA","IL10RB","IL20RA","IL20RB","IL22RA1","IL22RA2","IL12RB1","IL12RB2","IL23R","EPOR","GHR","PRLR","TPOR","CSF3R","LEPR","IFNAR1","IFNAR2","IFNLR1","IFNGR1","IFNGR2","EGFR","PDGFRA","PDGFRB","IRF9","CISH","BCL2","MCL1","BCL2L1","PIM1","MYC","CCND1","CCND2","CCND3","CDKN1A","AOX","GFAP","STAM","PTPN2","PTPN6","EP300","PIAS1","PIAS2","PIAS3","PIAS4","FHL1","PTPN11","GRB2","SOS","HRAS","RAF1","PIK3CA","PIK3CB","PIK3CD","PIK3R1","PIK3R2","PIK3R3","AKT","MTOR"]
    
    # 处理重复的观测名[6](@ref)
    adata.obs_names_make_unique()  
    adata.var_names_make_unique()
    
    # 筛选目标基因
    gene_mask = adata.var_names.isin(target_genes)
    adata = adata[:, gene_mask].copy()
    
    return adata

# 2. 差异分析
def perform_de_analysis(adata):
    # 确保分组列存在
    if 'group' not in adata.obs.columns:
        raise KeyError("分组信息需存储在obs的'group'列")
    
    # 运行差异分析[4](@ref)
    sc.tl.rank_genes_groups(
        adata,
        groupby='group',
        groups=['disease'],
        reference='normal',
        method='wilcoxon',
        pts=True,
        key_added='diff_genes'
    )
    
    # 提取结果并合并pts值
    diff_df = sc.get.rank_genes_groups_df(adata, group='disease', key='diff_genes')
    pts_values = adata.uns['diff_genes']['pts']['disease']
    diff_df['pts'] = diff_df['names'].map(pts_values)
    
    # 筛选显著差异基因[3](@ref)
    diff_df = diff_df[
        (diff_df['pvals_adj'] < 0.05) &
        (abs(diff_df['logfoldchanges']) > 1) &
        (diff_df['pts'] > 0.1) &
        (diff_df['names'].isin(adata.var_names))  # 过滤无效基因[1](@ref)
    ]
    
    return diff_df

# 3. 结果可视化
def visualize_results(adata, diff_df):
    # 创建子数据集
    valid_genes = diff_df['names'].tolist()
    subset_adata = adata[:, valid_genes].copy()
    
    # 转换稀疏矩阵[7](@ref)
    expr_data = subset_adata.X.toarray() if issparse(subset_adata.X) else subset_adata.X
    
    # 构建数据框
    expr_df = pd.DataFrame(
        expr_data,
        columns=subset_adata.var_names,
        index=subset_adata.obs.index
    )
    expr_df['Group'] = subset_adata.obs['group']
    
    # 绘制箱线图
    plt.figure(figsize=(15, len(valid_genes)*0.5))
    sns.set(style="whitegrid", font_scale=1.2)
    
    melt_df = pd.melt(
        expr_df,
        id_vars='Group',
        var_name='Gene',
        value_name='Expression'
    )
    
    sns.boxplot(
        x='Expression', 
        y='Gene', 
        hue='Group',
        data=melt_df,
        palette={'normal':'#1f77b4', 'disease':'#ff7f0e'},
        width=0.6,
        linewidth=1.2
    )
    plt.title('Target Genes Expression (Normal vs Disease)')
    plt.xlabel('Normalized Expression')
    plt.ylabel('')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig('filtered_genes_boxplot.pdf', dpi=300)
    plt.close()

# 主流程
if __name__ == "__main__":
    # 初始化参数
    DATA_PATH = "Kerationocytes.h5ad"
    
    # 预处理数据
    adata = preprocess_data(DATA_PATH)
    
    # 执行差异分析
    diff_df = perform_de_analysis(adata)
    
    # 分割上下调基因
    upregulated = diff_df[diff_df['logfoldchanges'] > 0]
    downregulated = diff_df[diff_df['logfoldchanges'] < 0]
    
    # 保存结果
    upregulated.to_csv('filtered_upregulated.csv', index=False)
    downregulated.to_csv('filtered_downregulated.csv', index=False)
    
    # 可视化
    if not diff_df.empty:
        visualize_results(adata, diff_df)
    else:
        print("未找到显著差异基因")
