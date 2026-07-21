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

def find_pairwise_band_diff_genes(adata, condition_col='site_band', species=None, out_dir=os.getcwd(),
                                  filter=True, logfoldchanges=0.25, pvals_adj=0.05, pct_nz_group=0.1,
                                  method='wilcoxon', gene_annot=None, 
                                  plot_volcano=True, plot_heatmap=True,
                                  volcano_thresholds=None, heatmap_top_genes=20,
                                  use_raw=True, layer=None):
    """
    查找条带间差异表达基因，进行两两比较：in vs on 和 on vs out
    集成火山图和热图可视化
    
    Args:
        adata: AnnData对象
        condition_col: 条件列名，默认'site_band'
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
    
    # 定义两两比较
    comparisons = [
        ('in', 'on'),  # in vs on
        ('on', 'out')  # on vs out
    ]
    
    # 检查数据中是否包含这些条带
    valid_bands = ['in', 'on', 'out']
    band_counts = adata.obs[condition_col].value_counts()
    available_bands = [band for band in valid_bands if band in band_counts.index]
    
    if len(available_bands) < 2:
        print(f"警告: 没有足够的条带进行比较。可用条带: {available_bands}")
        return
    
    print(f"可用的条带: {available_bands}")
    print(f"条带细胞数: {band_counts[available_bands].to_dict()}")
    
    # 创建输出目录
    os.makedirs(out_dir, exist_ok=True)
    
    # 设置火山图默认阈值
    if volcano_thresholds is None:
        volcano_thresholds = (logfoldchanges, pvals_adj)
    
    all_diff = []
    stats = []  # 用于存储统计信息
    
    for case, control in comparisons:
        # 检查case和control是否都存在
        if case not in available_bands or control not in available_bands:
            print(f"跳过比较 {case} vs {control}，因为缺少数据")
            continue
        
        comparison_name = f"{case}_vs_{control}"
        comp_dir = os.path.join(out_dir, comparison_name)
        os.makedirs(comp_dir, exist_ok=True)
        
        print(f"\n正在进行比较: {case} vs {control}")
        
        # 提取两个条带的细胞
        mask = adata.obs[condition_col].isin([case, control])
        adata_pair = adata[mask].copy()
        
        # 确保有足够的细胞
        case_cells = sum(adata_pair.obs[condition_col] == case)
        control_cells = sum(adata_pair.obs[condition_col] == control)
        
        if case_cells < 3 or control_cells < 3:
            print(f"  警告: 细胞数不足 ({case}: {case_cells}, {control}: {control_cells})，跳过")
            continue
        
        print(f"  {case}: {case_cells} 个细胞")
        print(f"  {control}: {control_cells} 个细胞")
        
        # 选择用于差异分析的数据
        if use_raw and adata_pair.raw is not None:
            # 使用原始矩阵
            print(f"  使用原始矩阵进行差异分析")
            raw_adata = adata_pair.raw.to_adata()
            raw_adata_pair = raw_adata[adata_pair.obs_names].copy()  # 确保是实际对象而不是视图
        elif layer is not None and layer in adata_pair.layers:
            # 使用指定层的数据
            print(f"  使用层 '{layer}' 进行差异分析")
            raw_adata_pair = adata_pair.copy()
            raw_adata_pair.X = raw_adata_pair.layers[layer]
        else:
            # 使用当前矩阵
            print(f"  使用当前矩阵进行差异分析")
            raw_adata_pair = adata_pair.copy()
        
        # 确保条件列是分类变量
        if not pd.api.types.is_categorical_dtype(raw_adata_pair.obs[condition_col]):
            raw_adata_pair.obs[condition_col] = raw_adata_pair.obs[condition_col].astype('category')
        
        # 差异分析
        try:
            sc.tl.rank_genes_groups(
                raw_adata_pair,  # 使用选择的数据进行差异分析
                groupby=condition_col, 
                method=method,
                groups=[case],  # 只分析case组
                reference=control,  # 以control为参考
                pts=True
            )
        except Exception as e:
            print(f"  差异分析失败: {str(e)}")
            continue
        
        # 获取差异基因结果
        try:
            diff_df = sc.get.rank_genes_groups_df(raw_adata_pair, group=case)
        except Exception as e:
            print(f"  获取差异基因结果失败: {str(e)}")
            continue
        
        if filter and not diff_df.empty:
            diff_df = diff_df[
                (abs(diff_df['logfoldchanges']) > logfoldchanges) &
                (diff_df['pvals_adj'] < pvals_adj) &
                (diff_df['pct_nz_group'] >= pct_nz_group)
            ].copy()
        
        if diff_df.empty:
            print(f"  没有找到满足条件的差异基因")
            # 添加空结果的统计
            stats.append({
                'comparison': comparison_name,
                'up_genes': 0,
                'down_genes': 0,
                'total': 0
            })
            continue
        
        diff_df['comparison'] = comparison_name
        diff_df['case'] = case
        diff_df['control'] = control
        all_diff.append(diff_df)
        
        # 保存结果
        result_file = os.path.join(comp_dir, f'{comparison_name}_diff_genes.xls')
        diff_df.to_csv(result_file, sep='\t', index=False)
        
        # 计算统计信息
        up_count = sum(diff_df['logfoldchanges'] > 0)
        down_count = sum(diff_df['logfoldchanges'] < 0)
        stats.append({
            'comparison': comparison_name,
            'up_genes': up_count,
            'down_genes': down_count,
            'total': len(diff_df)
        })
        
        # ================= 火山图绘制 =================
        if plot_volcano and not diff_df.empty:
            try:
                # 获取当前比较所有基因结果（包括不显著基因）
                full_df = sc.get.rank_genes_groups_df(raw_adata_pair, group=case)
                if full_df.empty:
                    print(f"  警告: 完整差异基因结果为空，跳过火山图")
                else:
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
                        if len(subset) > 0:
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
                    up_genes = diff_df[diff_df['logfoldchanges'] > 0].nsmallest(5, 'pvals_adj')
                    down_genes = diff_df[diff_df['logfoldchanges'] < 0].nsmallest(5, 'pvals_adj')
                    top_genes = pd.concat([up_genes, down_genes])['names'].tolist()[:10]  # 最多显示10个
                    
                    for _, row in full_df[full_df['names'].isin(top_genes)].iterrows():
                        plt.text(row['logfoldchanges'], 
                                row['-log10(padj)'] + 0.3,
                                row['names'],
                                fontsize=9,
                                ha='center')
                    
                    plt.title(f"Volcano Plot: {case} vs {control}", fontsize=14)
                    plt.xlabel("log2(Fold Change)", fontweight='bold')
                    plt.ylabel("-log10(Adjusted P-value)", fontweight='bold')
                    plt.legend(loc='best')
                    plt.tight_layout()
                    plt.savefig(os.path.join(comp_dir, f'{comparison_name}_volcano.pdf'), dpi=300)
                    plt.savefig(os.path.join(comp_dir, f'{comparison_name}_volcano.png'), dpi=300)
                    plt.close()
                    
            except Exception as e:
                print(f"火山图绘制失败: {str(e)}")
        
        # ================= 热图绘制 =================
        if plot_heatmap and not diff_df.empty and len(diff_df) >= 5:
            try:
                # 获取top差异基因
                n_genes = min(heatmap_top_genes, len(diff_df))
                up_genes = diff_df[diff_df['logfoldchanges'] > 0].nsmallest(n_genes//2, 'pvals_adj')
                down_genes = diff_df[diff_df['logfoldchanges'] < 0].nsmallest(n_genes//2, 'pvals_adj')
                top_genes = pd.concat([up_genes, down_genes])['names'].tolist()
                
                if len(top_genes) == 0:
                    print(f"  没有足够的差异基因绘制热图")
                    continue
                
                # 从原始adata中提取表达矩阵
                if hasattr(adata_pair.X, 'toarray'):
                    expr_data = adata_pair[:, top_genes].X.toarray()
                else:
                    expr_data = adata_pair[:, top_genes].X
                
                expr_df = pd.DataFrame(
                    data=expr_data,
                    index=adata_pair.obs[condition_col],
                    columns=top_genes
                )
                
                # 按条件计算平均表达量
                mean_expr = expr_df.groupby(level=0).mean()
                
                # 确保顺序
                mean_expr = mean_expr.loc[[case, control]]
                
                # Z-score标准化
                scaler = StandardScaler()
                scaled_expr = pd.DataFrame(
                    scaler.fit_transform(mean_expr.T).T,
                    index=mean_expr.index,
                    columns=mean_expr.columns
                )
                
                # 绘制热图
                fig_height = max(4, 0.5 * len(top_genes) + 2)
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
                
                plt.title(f"Expression Heatmap: {case} vs {control}\nTop {len(top_genes)} DE Genes", fontsize=12)
                plt.xlabel("Genes")
                plt.ylabel("Bands")
                plt.xticks(rotation=45, ha='right')
                plt.yticks(rotation=0)
                plt.tight_layout()
                plt.savefig(os.path.join(comp_dir, f'{comparison_name}_heatmap.pdf'), dpi=300)
                plt.savefig(os.path.join(comp_dir, f'{comparison_name}_heatmap.png'), dpi=300)
                plt.close()
                
            except Exception as e:
                print(f"热图绘制失败: {str(e)}")
        
        # 富集分析
        if not diff_df.empty and species:
            try:
                # 分别对上调基因和下调基因进行富集分析
                up_genes = diff_df[diff_df['logfoldchanges'] > 0]['names'].tolist()
                down_genes = diff_df[diff_df['logfoldchanges'] < 0]['names'].tolist()
                
                enrichment_dir = os.path.join(comp_dir, 'enrichment')
                os.makedirs(enrichment_dir, exist_ok=True)
                
                # 上调基因富集分析
                if up_genes and len(up_genes) >= 5:
                    up_dir = os.path.join(enrichment_dir, 'upregulated')
                    os.makedirs(up_dir, exist_ok=True)
                    
                    if species == 'human':
                        enr_up = gp.enrichr(
                            gene_list=up_genes,
                            gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
                            organism='Human',
                            cutoff=0.05,
                            outdir=up_dir
                        )
                        plot_enrichment(enr_up, f"{case} vs {control} - Upregulated", up_dir)
                    elif species == 'mouse':
                        enr_up = gp.enrichr(
                            gene_list=up_genes,
                            gene_sets=['GO_Biological_Process_2023', 'KEGG_2019_Mouse'],
                            organism='Mouse',
                            cutoff=0.05,
                            outdir=up_dir
                        )
                        plot_enrichment(enr_up, f"{case} vs {control} - Upregulated", up_dir)
                
                # 下调基因富集分析
                if down_genes and len(down_genes) >= 5:
                    down_dir = os.path.join(enrichment_dir, 'downregulated')
                    os.makedirs(down_dir, exist_ok=True)
                    
                    if species == 'human':
                        enr_down = gp.enrichr(
                            gene_list=down_genes,
                            gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
                            organism='Human',
                            cutoff=0.05,
                            outdir=down_dir
                        )
                        plot_enrichment(enr_down, f"{case} vs {control} - Downregulated", down_dir)
                    elif species == 'mouse':
                        enr_down = gp.enrichr(
                            gene_list=down_genes,
                            gene_sets=['GO_Biological_Process_2023', 'KEGG_2019_Mouse'],
                            organism='Mouse',
                            cutoff=0.05,
                            outdir=down_dir
                        )
                        plot_enrichment(enr_down, f"{case} vs {control} - Downregulated", down_dir)
                    
            except Exception as e:
                print(f"富集分析失败: {str(e)}")
    
    # 保存所有结果
    if all_diff:
        combined_diff = pd.concat(all_diff, ignore_index=True)
        combined_diff.rename(columns={'names': 'Symbol'}, inplace=True)
        combined_diff.to_csv(os.path.join(out_dir, 'all_pairwise_diff_genes.xls'), sep='\t', index=False)
        
    # 保存统计信息
    if stats:
        stats_df = pd.DataFrame(stats)
        stats_df.to_csv(os.path.join(out_dir, 'pairwise_diff_statistics.xls'), sep='\t', index=False)
        
        print(f"\n差异基因分析完成！")
        print(f"输出目录: {out_dir}")
        for _, row in stats_df.iterrows():
            print(f"  {row['comparison']}: {row['up_genes']} 个上调基因, {row['down_genes']} 个下调基因, 共 {row['total']} 个差异基因")
    else:
        print(f"\n没有找到任何差异基因")
    
    # 绘制联合热图
    if plot_heatmap and len(all_diff) >= 1:
        try:
            plot_combined_pairwise_heatmap(adata, condition_col, out_dir, heatmap_top_genes)
        except Exception as e:
            print(f"联合热图绘制失败: {str(e)}")

def plot_enrichment(enr, title, out_dir):
    """
    绘制富集分析结果
    """
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
                    
                    # 创建条形图
                    bars = plt.barh(y_pos, -np.log10(top_subset['Adjusted P-value']))
                    plt.yticks(y_pos, [term[:60] + '...' if len(term) > 60 else term for term in top_subset['Term']])
                    plt.xlabel('-log10(Adjusted P-value)')
                    plt.title(f'{title} - {gene_set}')
                    plt.tight_layout()
                    
                    # 保存图片
                    set_name = gene_set.replace(' ', '_').replace('/', '_')
                    plt.savefig(os.path.join(out_dir, f'{set_name}_enrichment.pdf'), dpi=300, bbox_inches='tight')
                    plt.savefig(os.path.join(out_dir, f'{set_name}_enrichment.png'), dpi=300, bbox_inches='tight')
                    plt.close()

def plot_combined_pairwise_heatmap(adata, condition_col, out_dir, top_genes=20):
    """
    绘制所有两两比较差异基因的联合热图
    """
    # 收集所有比较的top差异基因
    all_top_genes = set()
    comparisons = [
        ('in', 'on'),
        ('on', 'out')
    ]
    
    for case, control in comparisons:
        comp_name = f"{case}_vs_{control}"
        diff_file = os.path.join(out_dir, comp_name, f'{comp_name}_diff_genes.xls')
        if os.path.exists(diff_file):
            diff_df = pd.read_csv(diff_file, sep='\t')
            if not diff_df.empty:
                # 取top上调基因和下调基因
                n_genes = min(top_genes, len(diff_df))
                up_genes = diff_df[diff_df['logfoldchanges'] > 0].nsmallest(n_genes//2, 'pvals_adj')
                down_genes = diff_df[diff_df['logfoldchanges'] < 0].nsmallest(n_genes//2, 'pvals_adj')
                top_genes_comp = pd.concat([up_genes, down_genes])['names'].tolist()
                all_top_genes.update(top_genes_comp)
    
    all_top_genes = list(all_top_genes)
    if len(all_top_genes) == 0:
        print("没有找到足够的差异基因绘制联合热图")
        return
    
    # 只保留in, on, out的细胞
    mask = adata.obs[condition_col].isin(['in', 'on', 'out'])
    adata_filtered = adata[mask].copy()
    
    if len(adata_filtered) == 0:
        print("没有找到符合条件的细胞")
        return
    
    # 提取表达矩阵
    if hasattr(adata_filtered.X, 'toarray'):
        expr_data = adata_filtered[:, all_top_genes].X.toarray()
    else:
        expr_data = adata_filtered[:, all_top_genes].X
    
    expr_df = pd.DataFrame(
        data=expr_data,
        index=adata_filtered.obs[condition_col],
        columns=all_top_genes
    )
    
    # 按条带计算平均表达量
    mean_expr = expr_df.groupby(level=0).mean()
    
    # 确保顺序
    available_bands = [band for band in ['in', 'on', 'out'] if band in mean_expr.index]
    if not available_bands:
        print("没有可用的条带数据")
        return
    
    mean_expr = mean_expr.loc[available_bands]
    
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
    
    # 创建热图
    heatmap = plt.pcolor(scaled_expr, cmap=plt.cm.coolwarm, vmin=-2, vmax=2)
    
    # 添加颜色条
    plt.colorbar(heatmap, label='Z-score')
    
    # 设置坐标轴
    plt.yticks(np.arange(0.5, len(scaled_expr.index), 1), scaled_expr.index)
    plt.xticks(np.arange(0.5, len(scaled_expr.columns), 1), scaled_expr.columns, rotation=45, ha='right')
    
    plt.title(f"Combined Heatmap: All Bands - Top {top_genes} DE Genes per Comparison", fontsize=14)
    plt.xlabel("Genes")
    plt.ylabel("Bands")
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'combined_pairwise_heatmap.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'combined_pairwise_heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()

# 使用示例
if __name__ == "__main__":
    adata = sc.read_h5ad("adata_with_band_by_sample.h5ad")
    find_pairwise_band_diff_genes(
        adata,
        condition_col='site_band',  # 条带分组列
        species='human',
        out_dir='./pairwise_band_diff_results',
        logfoldchanges=0.5,
        pvals_adj=0.05,
        heatmap_top_genes=20,
        use_raw=True
    )
