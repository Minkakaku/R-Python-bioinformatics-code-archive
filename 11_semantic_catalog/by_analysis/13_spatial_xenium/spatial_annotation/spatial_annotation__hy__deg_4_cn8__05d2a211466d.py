import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gseapy as gp
import scanpy as sc
import seaborn as sns
import os
import re
from sklearn.preprocessing import StandardScaler

def find_group_diff_genes(adata, group_col, condition_col, species=None, out_dir=os.getcwd(),
                          filter=True, logfoldchanges=0.25, pvals_adj=0.05, pct_nz_group=0.1,
                          method='wilcoxon', gene_annot=None, 
                          plot_volcano=True, plot_heatmap=True,
                          volcano_thresholds=None, heatmap_top_genes=20,
                          use_raw=True, layer=None):
    """
    查找组间差异表达基因，集成火山图和热图可视化
    
    Args:
        adata: AnnData对象
        group_col: 分组列名（如细胞类型）
        condition_col: 条件列名（如处理组/对照组）
        species: 物种（'human'或'mouse'）
        out_dir: 输出目录
        filter: 是否过滤结果
        logfoldchanges: log2倍数变化阈值
        pvals_adj: 调整后p值阈值
        pct_nz_group: 组内非零表达比例阈值
        method: 差异分析方法
        gene_annot: 基因注释
        plot_volcano: 是否绘制火山图
        plot_heatmap: 是否绘制热图
        volcano_thresholds: 火山图阈值 (logfc_thresh, pval_thresh)
        heatmap_top_genes: 热图展示的top基因数
        use_raw: 是否使用原始矩阵进行差异分析
        layer: 如果提供，则使用指定层的数据
    """
    # 确保分组列存在
    if group_col not in adata.obs.columns or condition_col not in adata.obs.columns:
        raise ValueError(f"列 '{group_col}' 或 '{condition_col}' 不存在")
    
    os.makedirs(out_dir, exist_ok=True)
    all_diff = []
    
    # 设置火山图默认阈值
    if volcano_thresholds is None:
        volcano_thresholds = (logfoldchanges, pvals_adj)
    
    groups = adata.obs[group_col].unique()
    
    for group in groups:
        group_adata = adata[adata.obs[group_col] == group].copy()
        condition_counts = group_adata.obs[condition_col].value_counts()
        
        if len(condition_counts) < 2:
            print(f"分组 '{group}' 在条件列 '{condition_col}' 中少于2个类别，跳过分析")
            continue
        
        group_dir = os.path.join(out_dir, re.sub(r'[+/\\-]', "_", group))
        os.makedirs(group_dir, exist_ok=True)
        
        # 选择用于差异分析的数据
        if use_raw and adata.raw is not None:
            # 使用原始矩阵
            print(f"分组 '{group}': 使用原始矩阵进行差异分析")
            # 从原始数据中获取对应的子集
            raw_adata = adata.raw.to_adata()
            # 获取与当前group_adata相同的细胞索引
            raw_group_adata = raw_adata[group_adata.obs_names]
        elif layer is not None and layer in adata.layers:
            # 使用指定层的数据
            print(f"分组 '{group}': 使用层 '{layer}' 进行差异分析")
            raw_group_adata = group_adata.copy()
            raw_group_adata.X = raw_group_adata.layers[layer]
        else:
            # 使用当前矩阵
            print(f"分组 '{group}': 使用当前矩阵进行差异分析")
            raw_group_adata = group_adata.copy()
        
        # 差异分析
        sc.tl.rank_genes_groups(
            raw_group_adata,  # 使用选择的数据进行差异分析
            groupby=condition_col, 
            method=method,
            reference='rest',
            pts=True
        )
        
        result = raw_group_adata.uns['rank_genes_groups']
        conditions = result['names'].dtype.names
        
        for condition in conditions:
            diff_df = sc.get.rank_genes_groups_df(raw_group_adata, group=condition)
            
            if filter:
                diff_df = diff_df[
                    (diff_df['logfoldchanges'] > logfoldchanges) &
                    (diff_df['pvals_adj'] < pvals_adj) &
                    (diff_df['pct_nz_group'] >= pct_nz_group)
                ]
            
            diff_df['group'] = group
            diff_df['condition'] = condition
            all_diff.append(diff_df)
            
            condition_file = os.path.join(group_dir, f'{condition}_diff_genes.xls')
            diff_df.to_csv(condition_file, sep='\t', index=False)
            
            # ================= 火山图绘制 =================
            if plot_volcano and not diff_df.empty:
                try:
                    # 获取当前条件所有基因结果（包括不显著基因）
                    full_df = sc.get.rank_genes_groups_df(raw_group_adata, group=condition)
                    full_df['-log10(padj)'] = -np.log10(full_df['pvals_adj'])
                    
                    # 标记基因类型
                    logfc_thresh, pval_thresh = volcano_thresholds
                    full_df['type'] = 'nosig'
                    full_df.loc[(full_df['logfoldchanges'] > logfc_thresh) & 
                               (full_df['pvals_adj'] < pval_thresh), 'type'] = 'up'
                    full_df.loc[(full_df['logfoldchanges'] < -logfc_thresh) & 
                               (full_df['pvals_adj'] < pval_thresh), 'type'] = 'down'
                    
                    # 绘制火山图
                    plt.figure(figsize=(10, 8))
                    colors = {"down": "#377EB8", "nosig": "grey", "up": "#E41A1C"}
                    
                    for gene_type, color in colors.items():
                        subset = full_df[full_df['type'] == gene_type]
                        plt.scatter(subset['logfoldchanges'], 
                                   subset['-log10(padj)'],
                                   c=color, 
                                   s=15,
                                   alpha=0.6,
                                   label=f"{gene_type} ({len(subset)})")
                    
                    # 添加阈值线
                    plt.axvline(x=logfc_thresh, color='dimgrey', linestyle='--', linewidth=1)
                    plt.axvline(x=-logfc_thresh, color='dimgrey', linestyle='--', linewidth=1)
                    plt.axhline(y=-np.log10(pval_thresh), color='dimgrey', linestyle='--', linewidth=1)
                    
                    # 标记top基因
                    top_genes = diff_df.nsmallest(5, 'pvals_adj')['names'].tolist()
                    for _, row in full_df[full_df['names'].isin(top_genes)].iterrows():
                        plt.text(row['logfoldchanges'], 
                                row['-log10(padj)'] + 0.3,
                                row['names'],
                                fontsize=9,
                                ha='center')
                    
                    plt.title(f"Volcano Plot: {group} - {condition}", fontsize=14)
                    plt.xlabel("log2(Fold Change)", fontweight='bold')
                    plt.ylabel("-log10(Adjusted P-value)", fontweight='bold')
                    plt.legend(loc='best')
                    plt.tight_layout()
                    plt.savefig(os.path.join(group_dir, f'{condition}_volcano.pdf'), dpi=300)
                    plt.close()
                    
                except Exception as e:
                    print(f"火山图绘制失败: {str(e)}")
            
            # ================= 热图绘制 =================
            if plot_heatmap and not diff_df.empty and len(diff_df) >= 5:
                try:
                    # 获取top差异基因
                    top_genes = diff_df.nsmallest(heatmap_top_genes, 'pvals_adj')['names'].tolist()
                    
                    # 从原始adata中提取表达矩阵（使用归一化后的数据展示）
                    # 这里使用原始group_adata，它包含归一化后的数据
                    if hasattr(group_adata.X, 'toarray'):
                        expr_data = group_adata[:, top_genes].X.toarray()
                    else:
                        expr_data = group_adata[:, top_genes].X
                    
                    expr_df = pd.DataFrame(
                        data=expr_data,
                        index=group_adata.obs[condition_col],
                        columns=top_genes
                    )
                    
                    # 按条件计算平均表达量
                    mean_expr = expr_df.groupby(level=0).mean()
                    
                    # Z-score标准化
                    scaler = StandardScaler()
                    scaled_expr = pd.DataFrame(
                        scaler.fit_transform(mean_expr.T).T,
                        index=mean_expr.index,
                        columns=mean_expr.columns
                    )
                    
                    # 绘制热图
                    plt.figure(figsize=(10, 0.5 * len(top_genes) + 2))
                    sns.heatmap(
                        scaled_expr,
                        cmap='coolwarm',
                        center=0,
                        annot=True,
                        fmt=".2f",
                        linewidths=0.5,
                        cbar_kws={'label': 'Z-score'}
                    )
                    
                    plt.title(f"Expression Heatmap: {group}\nTop {heatmap_top_genes} DE Genes", fontsize=12)
                    plt.xlabel("Genes")
                    plt.ylabel("Conditions")
                    plt.xticks(rotation=45)
                    plt.tight_layout()
                    plt.savefig(os.path.join(group_dir, f'{condition}_heatmap.pdf'), dpi=300)
                    plt.close()
                    
                except Exception as e:
                    print(f"热图绘制失败: {str(e)}")
            
            # 富集分析
            if not diff_df.empty and species:
                try:
                    genes = diff_df['names'].tolist()
                    if species == 'human':
                        enr = gp.enrichr(
                            gene_list=genes,
                            gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
                            organism='Human',
                            cutoff=0.05,
                            outdir=os.path.join(group_dir, f'{condition}_enrichment')
                        )
                    elif species == 'mouse':
                        enr = gp.enrichr(
                            gene_list=genes,
                            gene_sets=['GO_Biological_Process_2023', 'KEGG_2019_Mouse'],
                            organism='Mouse',
                            cutoff=0.05,
                            outdir=os.path.join(group_dir, f'{condition}_enrichment')
                        )
                except Exception as e:
                    print(f"富集分析失败: {str(e)}")
    
    # 保存所有结果
    if all_diff:
        combined_diff = pd.concat(all_diff, ignore_index=True)
        combined_diff.rename(columns={'names': 'Symbol'}, inplace=True)
        combined_diff.to_csv(os.path.join(out_dir, 'all_group_diff_genes.xls'), sep='\t', index=False)
        
        stats = combined_diff.groupby(['group', 'condition']).size().reset_index(name='count')
        stats.to_csv(os.path.join(out_dir, 'diff_gene_statistics.xls'), sep='\t', index=False)

# 使用示例

adata = sc.read_h5ad("CN8_cells.h5ad")
find_group_diff_genes(
    adata,
    group_col='annotation',  # 细胞类型分组列
    condition_col='cancer_type',  # 处理条件列
    species='human',
    out_dir='./group_diff_results_cn8',
    logfoldchanges=0.5,
    pvals_adj=0.01,
    use_raw=True  # 使用原始矩阵进行差异分析
)
