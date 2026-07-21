import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import cosg as cosg
import numpy as np
import gseapy as gp
import scanpy as sc
import os
import re

def find_group_diff_genes(adata, group_col, condition_col, species=None, out_dir=os.getcwd(),
                          filter=True, logfoldchanges=0.25, pvals_adj=0.05, pct_nz_group=0.1,
                          method='wilcoxon', gene_annot=None):
    """
    查找组间差异表达基因并进行功能富集分析。

    Args:
        adata: 输入的AnnData对象
        group_col: 分组列名（如细胞类型）
        condition_col: 条件列名（如处理组/对照组）
        species: 物种，如'human'或'mouse'
        out_dir: 输出目录
        filter: 是否过滤结果
        logfoldchanges: 最小log2倍数变化阈值
        pvals_adj: 调整后p值阈值
        pct_nz_group: 组内非零表达比例阈值
        method: 差异分析方法（默认'wilcoxon'）
        gene_annot: 基因注释信息
    """
    # 确保分组列和条件列存在
    if group_col not in adata.obs.columns or condition_col not in adata.obs.columns:
        raise ValueError(f"分组列 '{group_col}' 或条件列 '{condition_col}' 不存在于adata.obs中")
    
    # 创建主输出目录
    os.makedirs(out_dir, exist_ok=True)
    all_diff = []
    
    # 获取所有分组类型
    groups = adata.obs[group_col].unique()
    
    for group in groups:
        # 筛选当前分组的数据
        group_adata = adata[adata.obs[group_col] == group].copy()
        
        # 检查条件列是否有足够的样本
        condition_counts = group_adata.obs[condition_col].value_counts()
        if len(condition_counts) < 2:
            print(f"分组 '{group}' 在条件列 '{condition_col}' 中少于2个类别，跳过分析")
            continue
        
        # 创建分组输出目录
        group_dir = os.path.join(out_dir, re.sub(r'[+/\\-]', "_", group))
        os.makedirs(group_dir, exist_ok=True)
        
        # 设置差异分析参考组
        sc.tl.rank_genes_groups(
            group_adata, 
            groupby=condition_col, 
            method=method,
            reference='rest',  # 与所有其他组比较
            pts=True
        )
        
        # 获取差异分析结果
        result = group_adata.uns['rank_genes_groups']
        conditions = result['names'].dtype.names
        
        for condition in conditions:
            # 提取当前条件的差异基因
            diff_df = sc.get.rank_genes_groups_df(group_adata, group=condition)
            
            # 过滤结果
            if filter:
                diff_df = diff_df[
                    (diff_df['logfoldchanges'] > logfoldchanges) &
                    (diff_df['pvals_adj'] < pvals_adj) &
                    (diff_df['pct_nz_group'] >= pct_nz_group)
                ]
            
            # 添加分组和条件信息
            diff_df['group'] = group
            diff_df['condition'] = condition
            all_diff.append(diff_df)
            
            # 保存当前条件的差异基因
            condition_file = os.path.join(group_dir, f'{condition}_diff_genes.xls')
            diff_df.to_csv(condition_file, sep='\t', index=False)
            
            # 富集分析（仅当有差异基因时）
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
                    print(f"分组 '{group}' 条件 '{condition}' 富集分析失败: {str(e)}")
    
    # 保存所有结果
    if all_diff:
        combined_diff = pd.concat(all_diff, ignore_index=True)
        combined_diff.rename(columns={'names': 'Symbol'}, inplace=True)
        combined_diff.to_csv(os.path.join(out_dir, 'all_group_diff_genes.xls'), sep='\t', index=False)
        
        # 生成统计信息
        stats = combined_diff.groupby(['group', 'condition']).size().reset_index(name='count')
        stats.to_csv(os.path.join(out_dir, 'diff_gene_statistics.xls'), sep='\t', index=False)
        
        # 可视化（示例：可扩展更多可视化）
        try:
            sc.pl.rank_genes_groups_dotplot(
                adata, 
                groups=groups[:min(5, len(groups))],  # 最多显示5个组
                n_genes=5,
                groupby=group_col,
                show=False,
                save=os.path.join(out_dir, 'group_diff_dotplot.pdf')
            )
        except Exception as e:
            print(f"可视化失败: {str(e)}")

adata = sc.read_h5ad("adata_celltype.h5ad")
find_group_diff_genes(
    adata,
    group_col='celltype',  # 细胞类型分组列
    condition_col='group',  # 处理条件列
    species='human',
    out_dir='./group_diff_results',
    logfoldchanges=0.5,
    pvals_adj=0.01
)
