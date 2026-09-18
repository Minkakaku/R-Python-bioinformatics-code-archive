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

def find_band_diff_genes(adata, condition_col='mucus_zone', species=None, out_dir=os.getcwd(),
                         filter=True, logfoldchanges=0.25, pvals_adj=0.05, pct_nz_group=0.1,
                         method='wilcoxon', gene_annot=None, 
                         plot_volcano=True, plot_heatmap=True,
                         volcano_thresholds=None, heatmap_top_genes=20,
                         use_raw=True, layer=None):
    """
    查找条带间差异表达基因，集成火山图和热图可视化
    
    Args:
        adata: AnnData对象
        condition_col: 条件列名，默认'mucus_zone'
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
    # 确保条件列存在
    if condition_col not in adata.obs.columns:
        raise ValueError(f"列 '{condition_col}' 不存在")
    
    # 只保留in, on, out的数据
    valid_bands = ['inner','outer']
    mask = adata.obs[condition_col].isin(valid_bands)
    adata_subset = adata[mask].copy()
    
    if len(adata_subset) == 0:
        print(f"警告: 没有找到'mucus_zone'为{valid_bands}的细胞，跳过分析")
        return
    
    print(f"使用 {len(adata_subset)} 个细胞进行分析 (bands: {valid_bands})")
    print(f"条带分布: {adata_subset.obs[condition_col].value_counts().to_dict()}")
    
    os.makedirs(out_dir, exist_ok=True)
    
    # 设置火山图默认阈值
    if volcano_thresholds is None:
        volcano_thresholds = (logfoldchanges, pvals_adj)
    
    # 选择用于差异分析的数据
    if use_raw and adata_subset.raw is not None:
        # 使用原始矩阵
        print(f"使用原始矩阵进行差异分析")
        # 从原始数据中获取对应的子集
        raw_adata = adata_subset.raw.to_adata()
        # 获取与当前adata_subset相同的细胞索引
        raw_adata_subset = raw_adata[adata_subset.obs_names]
    elif layer is not None and layer in adata_subset.layers:
        # 使用指定层的数据
        print(f"使用层 '{layer}' 进行差异分析")
        raw_adata_subset = adata_subset.copy()
        raw_adata_subset.X = raw_adata_subset.layers[layer]
    else:
        # 使用当前矩阵
        print(f"使用当前矩阵进行差异分析")
        raw_adata_subset = adata_subset.copy()
    
    # 差异分析
    sc.tl.rank_genes_groups(
        raw_adata_subset,  # 使用选择的数据进行差异分析
        groupby=condition_col, 
        method=method,
        reference='rest',
        pts=True
    )
    
    result = raw_adata_subset.uns['rank_genes_groups']
    conditions = result['names'].dtype.names
    
    all_diff = []
    
    for condition in conditions:
        if condition not in valid_bands:
            print(f"跳过非目标条件: {condition}")
            continue
            
        diff_df = sc.get.rank_genes_groups_df(raw_adata_subset, group=condition)
        
        if filter:
            diff_df = diff_df[
                (diff_df['logfoldchanges'] > logfoldchanges) &
                (diff_df['pvals_adj'] < pvals_adj) &
                (diff_df['pct_nz_group'] >= pct_nz_group)
            ]
        
        diff_df['condition'] = condition
        all_diff.append(diff_df)
        
        condition_file = os.path.join(out_dir, f'{condition}_band_diff_genes.xls')
        diff_df.to_csv(condition_file, sep='\t', index=False)
        
        # ================= 火山图绘制 =================
        if plot_volcano and not diff_df.empty:
            try:
                # 获取当前条件所有基因结果（包括不显著基因）
                full_df = sc.get.rank_genes_groups_df(raw_adata_subset, group=condition)
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
                
                plt.title(f"Volcano Plot: {condition} band", fontsize=14)
                plt.xlabel("log2(Fold Change)", fontweight='bold')
                plt.ylabel("-log10(Adjusted P-value)", fontweight='bold')
                plt.legend(loc='best')
                plt.tight_layout()
                plt.savefig(os.path.join(out_dir, f'{condition}_band_volcano.pdf'), dpi=300)
                plt.savefig(os.path.join(out_dir, f'{condition}_band_volcano.png'), dpi=300)
                plt.close()
                
            except Exception as e:
                print(f"火山图绘制失败: {str(e)}")
        
        # ================= 热图绘制 =================
        if plot_heatmap and not diff_df.empty and len(diff_df) >= 5:
            try:
                # 获取top差异基因
                top_genes = diff_df.nsmallest(heatmap_top_genes, 'pvals_adj')['names'].tolist()
                
                # 从原始adata中提取表达矩阵（使用归一化后的数据展示）
                if hasattr(adata_subset.X, 'toarray'):
                    expr_data = adata_subset[:, top_genes].X.toarray()
                else:
                    expr_data = adata_subset[:, top_genes].X
                
                expr_df = pd.DataFrame(
                    data=expr_data,
                    index=adata_subset.obs[condition_col],
                    columns=top_genes
                )
                
                # 按条件计算平均表达量
                mean_expr = expr_df.groupby(level=0).mean()
                
                # 确保只显示in, on, out
                mean_expr = mean_expr.loc[valid_bands]
                
                # Z-score标准化
                scaler = StandardScaler()
                scaled_expr = pd.DataFrame(
                    scaler.fit_transform(mean_expr.T).T,
                    index=mean_expr.index,
                    columns=mean_expr.columns
                )
                
                # 绘制热图
                fig_height = max(6, 0.5 * len(top_genes) + 2)
                plt.figure(figsize=(12, fig_height))
                sns.heatmap(
                    scaled_expr,
                    cmap='coolwarm',
                    center=0,
                    annot=True,
                    fmt=".2f",
                    linewidths=0.5,
                    cbar_kws={'label': 'Z-score'},
                    square=True
                )
                
                plt.title(f"Expression Heatmap: {condition} band - Top {heatmap_top_genes} DE Genes", fontsize=12)
                plt.xlabel("Genes")
                plt.ylabel("Bands")
                plt.xticks(rotation=45, ha='right')
                plt.yticks(rotation=0)
                plt.tight_layout()
                plt.savefig(os.path.join(out_dir, f'{condition}_band_heatmap.pdf'), dpi=300)
                plt.savefig(os.path.join(out_dir, f'{condition}_band_heatmap.png'), dpi=300)
                plt.close()
                
            except Exception as e:
                print(f"热图绘制失败: {str(e)}")
        
        # 富集分析
        if not diff_df.empty and species:
            try:
                genes = diff_df['names'].tolist()
                enrichment_dir = os.path.join(out_dir, f'{condition}_band_enrichment')
                os.makedirs(enrichment_dir, exist_ok=True)
                
                if species == 'human':
                    enr = gp.enrichr(
                        gene_list=genes,
                        gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
                        organism='Human',
                        cutoff=0.05,
                        outdir=enrichment_dir
                    )
                elif species == 'mouse':
                    enr = gp.enrichr(
                        gene_list=genes,
                        gene_sets=['GO_Biological_Process_2023', 'KEGG_2019_Mouse'],
                        organism='Mouse',
                        cutoff=0.05,
                        outdir=enrichment_dir
                    )
                    
                # 绘制富集分析图
                if enr is not None and enr.results is not None and not enr.results.empty:
                    # 绘制top富集条目
                    for gene_set in ['GO_Biological_Process_2023', 'KEGG_2021_Human', 'KEGG_2019_Mouse']:
                        if gene_set in enr.results['Gene_set'].values:
                            subset = enr.results[enr.results['Gene_set'] == gene_set]
                            if not subset.empty:
                                # 取top 10
                                top_subset = subset.nsmallest(10, 'Adjusted P-value')
                                
                                plt.figure(figsize=(10, 6))
                                y_pos = range(len(top_subset))
                                
                                plt.barh(y_pos, -np.log10(top_subset['Adjusted P-value']))
                                plt.yticks(y_pos, top_subset['Term'])
                                plt.xlabel('-log10(Adjusted P-value)')
                                plt.title(f'{condition} band - {gene_set} Enrichment')
                                plt.tight_layout()
                                plt.savefig(os.path.join(enrichment_dir, f'{condition}_{gene_set}_enrichment.pdf'), dpi=300)
                                plt.savefig(os.path.join(enrichment_dir, f'{condition}_{gene_set}_enrichment.png'), dpi=300)
                                plt.close()
                                
            except Exception as e:
                print(f"富集分析失败: {str(e)}")
    
    # 保存所有结果
    if all_diff:
        combined_diff = pd.concat(all_diff, ignore_index=True)
        combined_diff.rename(columns={'names': 'Symbol'}, inplace=True)
        combined_diff.to_csv(os.path.join(out_dir, 'all_band_diff_genes.xls'), sep='\t', index=False)
        
        # 添加统计信息
        stats = combined_diff['condition'].value_counts().reset_index()
        stats.columns = ['band', 'num_diff_genes']
        stats.to_csv(os.path.join(out_dir, 'band_diff_gene_statistics.xls'), sep='\t', index=False)
        
        print(f"\n差异基因分析完成！")
        print(f"输出目录: {out_dir}")
        for band in valid_bands:
            band_diff = combined_diff[combined_diff['condition'] == band]
            print(f"  {band} band: {len(band_diff)} 个差异基因")
    
    # 绘制联合热图
    if plot_heatmap and len(all_diff) >= 2:
        try:
            plot_combined_heatmap(adata_subset, condition_col, out_dir, heatmap_top_genes)
        except Exception as e:
            print(f"联合热图绘制失败: {str(e)}")

def plot_combined_heatmap(adata, condition_col, out_dir, top_genes=20):
    """
    绘制所有条带差异基因的联合热图
    """
    # 收集所有条带的top差异基因
    all_top_genes = set()
    bands = ['inner','outer']
    
    for band in bands:
        diff_file = os.path.join(out_dir, f'{band}_band_diff_genes.xls')
        if os.path.exists(diff_file):
            diff_df = pd.read_csv(diff_file, sep='\t')
            if not diff_df.empty:
                top_genes_band = diff_df.nsmallest(top_genes, 'pvals_adj')['names'].tolist()
                all_top_genes.update(top_genes_band)
    
    all_top_genes = list(all_top_genes)
    if len(all_top_genes) == 0:
        print("没有找到足够的差异基因绘制联合热图")
        return
    
    # 提取表达矩阵
    if hasattr(adata.X, 'toarray'):
        expr_data = adata[:, all_top_genes].X.toarray()
    else:
        expr_data = adata[:, all_top_genes].X
    
    expr_df = pd.DataFrame(
        data=expr_data,
        index=adata.obs[condition_col],
        columns=all_top_genes
    )
    
    # 按条带计算平均表达量
    mean_expr = expr_df.groupby(level=0).mean()
    
    # 确保顺序
    mean_expr = mean_expr.loc[bands]
    
    # Z-score标准化
    scaler = StandardScaler()
    scaled_expr = pd.DataFrame(
        scaler.fit_transform(mean_expr.T).T,
        index=mean_expr.index,
        columns=mean_expr.columns
    )
    
    # 绘制联合热图
    fig_height = max(8, 0.5 * len(all_top_genes) + 2)
    plt.figure(figsize=(16, fig_height))
    
    # 创建颜色映射
    cmap = plt.cm.coolwarm
    heatmap = plt.pcolor(scaled_expr, cmap=cmap, vmin=-2, vmax=2)
    
    # 添加颜色条
    plt.colorbar(heatmap, label='Z-score')
    
    # 设置坐标轴
    plt.yticks(np.arange(0.5, len(scaled_expr.index), 1), scaled_expr.index)
    plt.xticks(np.arange(0.5, len(scaled_expr.columns), 1), scaled_expr.columns, rotation=45, ha='right')
    
    plt.title(f"Combined Heatmap: All Bands - Top {top_genes} DE Genes per Band", fontsize=14)
    plt.xlabel("Genes")
    plt.ylabel("Bands")
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'combined_band_heatmap.pdf'), dpi=300)
    plt.savefig(os.path.join(out_dir, 'combined_band_heatmap.png'), dpi=300)
    plt.close()

# 使用示例
adata = sc.read_h5ad("extracted_samples.h5ad")
find_band_diff_genes(
    adata,
    condition_col='mucus_zone',  # 条带分组列
    species='human',
    out_dir='./band_diff_results_all',
    logfoldchanges=0.5,
    pvals_adj=0.05,
    heatmap_top_genes=25,
    use_raw=True
)
