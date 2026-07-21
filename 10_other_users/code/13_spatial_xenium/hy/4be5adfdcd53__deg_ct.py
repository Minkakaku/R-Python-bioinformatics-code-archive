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
from itertools import combinations

def find_group_diff_genes_site_band(adata, group_col='annotation', condition_col='site_band', 
                                    comparisons=[('in', 'on'), ('in', 'out'), ('on', 'out')],
                                    species=None, out_dir=os.getcwd(),
                                    filter=True, logfoldchanges=0.25, pvals_adj=0.05, pct_nz_group=0.1,
                                    method='wilcoxon', gene_annot=None, 
                                    plot_volcano=True, plot_heatmap=True,
                                    volcano_thresholds=None, heatmap_top_genes=20,
                                    use_raw=True, layer=None):
    """
    针对site_band的三个分组（in, on, out）进行细胞类型间的差异表达分析
    
    Args:
        adata: AnnData对象
        group_col: 细胞类型分组列名（默认为'annotation'）
        condition_col: 条件列名（默认为'site_band'）
        comparisons: 要进行的比较对列表
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
    
    # 过滤只保留in, on, out的细胞
    valid_site_bands = ['in', 'on', 'out']
    mask = adata.obs[condition_col].isin(valid_site_bands)
    adata_filtered = adata[mask].copy()
    
    if adata_filtered.n_obs == 0:
        raise ValueError(f"在条件列 '{condition_col}' 中未找到有效的分组: {valid_site_bands}")
    
    print(f"过滤后细胞数: {adata_filtered.n_obs}")
    print(f"site_band分布: {adata_filtered.obs[condition_col].value_counts().to_dict()}")
    
    os.makedirs(out_dir, exist_ok=True)
    all_diff = []
    
    # 设置火山图默认阈值
    if volcano_thresholds is None:
        volcano_thresholds = (logfoldchanges, pvals_adj)
    
    # 获取所有细胞类型
    cell_types = adata_filtered.obs[group_col].unique()
    print(f"找到的细胞类型: {cell_types}")
    
    for cell_type in cell_types:
        # 提取当前细胞类型的数据
        cell_type_adata = adata_filtered[adata_filtered.obs[group_col] == cell_type].copy()
        
        if cell_type_adata.n_obs < 10:  # 如果细胞数太少，跳过
            print(f"细胞类型 '{cell_type}' 细胞数太少 ({cell_type_adata.n_obs})，跳过分析")
            continue
        
        # 检查当前细胞类型是否包含所有三个分组
        site_band_counts = cell_type_adata.obs[condition_col].value_counts()
        print(f"\n细胞类型 '{cell_type}' 的site_band分布: {site_band_counts.to_dict()}")
        
        # 对每个比较对进行分析
        for comp in comparisons:
            group1, group2 = comp
            comp_name = f"{group1}_vs_{group2}"
            
            # 检查两个分组是否存在
            if group1 not in site_band_counts.index or group2 not in site_band_counts.index:
                print(f"  跳过比较 {comp_name}: 分组不存在")
                continue
            
            # 确保每个分组至少有3个细胞
            if site_band_counts[group1] < 3 or site_band_counts[group2] < 3:
                print(f"  跳过比较 {comp_name}: 分组细胞数不足 (min=3)")
                continue
            
            # 创建输出目录
            comp_dir = os.path.join(out_dir, f'deg_ct_{comp_name}', re.sub(r'[+/\\-]', "_", cell_type))
            os.makedirs(comp_dir, exist_ok=True)
            
            # 提取两个分组的数据
            comp_mask = cell_type_adata.obs[condition_col].isin([group1, group2])
            comp_adata = cell_type_adata[comp_mask].copy()
            
            # 确保条件列是分类变量
            comp_adata.obs[condition_col] = comp_adata.obs[condition_col].astype('category')
            
            # 选择用于差异分析的数据
            if use_raw and adata.raw is not None:
                # 使用原始矩阵
                raw_adata = adata.raw.to_adata()
                raw_comp_adata = raw_adata[comp_adata.obs_names]
            elif layer is not None and layer in adata.layers:
                # 使用指定层的数据
                raw_comp_adata = comp_adata.copy()
                raw_comp_adata.X = raw_comp_adata.layers[layer]
            else:
                # 使用当前矩阵
                raw_comp_adata = comp_adata.copy()
            
            print(f"  分析 {cell_type}: {comp_name} (n={site_band_counts[group1]} vs n={site_band_counts[group2]})")
            
            try:
                # 方法1: 使用rank_genes_groups，指定groups和reference
                sc.tl.rank_genes_groups(
                    raw_comp_adata,
                    groupby=condition_col,
                    groups=[group2],  # 比较组
                    reference=group1,  # 参考组
                    method=method,
                    pts=True
                )
                
                # 获取差异表达结果
                result = raw_comp_adata.uns['rank_genes_groups']
                # 获取group2相对于group1的差异基因
                diff_df = sc.get.rank_genes_groups_df(raw_comp_adata, group=group2)
                
                if diff_df.empty:
                    print(f"    没有找到差异表达基因")
                    continue
                
                # 添加元信息
                diff_df['cell_type'] = cell_type
                diff_df['comparison'] = comp_name
                diff_df['group1'] = group1
                diff_df['group2'] = group2
                diff_df['n_group1'] = site_band_counts[group1]
                diff_df['n_group2'] = site_band_counts[group2]
                
                # 过滤差异基因
                if filter:
                    diff_df_filtered = diff_df[
                        (diff_df['logfoldchanges'].abs() > logfoldchanges) &
                        (diff_df['pvals_adj'] < pvals_adj) &
                        (diff_df['pct_nz_group'] >= pct_nz_group)
                    ].copy()
                else:
                    diff_df_filtered = diff_df.copy()
                
                # 保存结果
                condition_file = os.path.join(comp_dir, f'{cell_type}_{comp_name}_diff_genes.xls')
                diff_df_filtered.to_csv(condition_file, sep='\t', index=False)
                all_diff.append(diff_df_filtered)
                
                print(f"    发现差异基因: {len(diff_df_filtered)} 个 (过滤阈值: |logFC|>{logfoldchanges}, padj<{pvals_adj})")
                
                # ================= 火山图绘制 =================
                if plot_volcano and not diff_df_filtered.empty:
                    try:
                        # 准备火山图数据
                        volcano_df = diff_df.copy()
                        volcano_df['-log10(padj)'] = -np.log10(volcano_df['pvals_adj'])
                        
                        # 标记基因类型
                        logfc_thresh, pval_thresh = volcano_thresholds
                        volcano_df['type'] = 'nosig'
                        volcano_df.loc[(volcano_df['logfoldchanges'] > logfc_thresh) & 
                                      (volcano_df['pvals_adj'] < pval_thresh), 'type'] = 'up'
                        volcano_df.loc[(volcano_df['logfoldchanges'] < -logfc_thresh) & 
                                      (volcano_df['pvals_adj'] < pval_thresh), 'type'] = 'down'
                        
                        # 绘制火山图
                        plt.figure(figsize=(10, 8))
                        colors = {"down": "#377EB8", "nosig": "grey", "up": "#E41A1C"}
                        
                        for gene_type, color in colors.items():
                            subset = volcano_df[volcano_df['type'] == gene_type]
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
                        if not diff_df_filtered.empty:
                            # 获取上调基因中最显著的
                            up_genes = diff_df_filtered[diff_df_filtered['logfoldchanges'] > 0]
                            down_genes = diff_df_filtered[diff_df_filtered['logfoldchanges'] < 0]
                            
                            if not up_genes.empty:
                                top_up = up_genes.nsmallest(3, 'pvals_adj')['names'].tolist()
                            else:
                                top_up = []
                            
                            if not down_genes.empty:
                                top_down = down_genes.nsmallest(3, 'pvals_adj')['names'].tolist()
                            else:
                                top_down = []
                            
                            top_genes = list(set(top_up + top_down))
                            
                            for _, row in volcano_df[volcano_df['names'].isin(top_genes)].iterrows():
                                plt.text(row['logfoldchanges'], 
                                        row['-log10(padj)'] + 0.3,
                                        row['names'],
                                        fontsize=9,
                                        ha='center')
                        
                        plt.title(f"Volcano Plot: {cell_type} ({comp_name})", fontsize=14)
                        plt.xlabel(f"log2({group2}/{group1})", fontweight='bold')
                        plt.ylabel("-log10(Adjusted P-value)", fontweight='bold')
                        plt.legend(loc='best')
                        plt.tight_layout()
                        plt.savefig(os.path.join(comp_dir, f'{cell_type}_{comp_name}_volcano.pdf'), dpi=300)
                        plt.close()
                        
                    except Exception as e:
                        print(f"    火山图绘制失败: {str(e)}")
                
                # ================= 热图绘制 =================
                if plot_heatmap and not diff_df_filtered.empty and len(diff_df_filtered) >= 5:
                    try:
                        # 获取top差异基因
                        n_genes = min(heatmap_top_genes, len(diff_df_filtered))
                        # 按p值排序，取最显著的基因
                        top_genes = diff_df_filtered.nsmallest(n_genes, 'pvals_adj')['names'].tolist()
                        
                        # 从原始adata中提取表达矩阵
                        if hasattr(comp_adata.X, 'toarray'):
                            expr_data = comp_adata[:, top_genes].X.toarray()
                        else:
                            expr_data = comp_adata[:, top_genes].X
                        
                        expr_df = pd.DataFrame(
                            data=expr_data,
                            index=comp_adata.obs[condition_col],
                            columns=top_genes
                        )
                        
                        # 按条件计算平均表达量
                        mean_expr = expr_df.groupby(level=0).mean()
                        
                        # 确保两个分组都在
                        if group1 in mean_expr.index and group2 in mean_expr.index:
                            # 按group1, group2的顺序排序
                            mean_expr = mean_expr.loc[[group1, group2]]
                            
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
                            
                            plt.title(f"Expression Heatmap: {cell_type} ({comp_name})\nTop {n_genes} DE Genes", fontsize=12)
                            plt.xlabel("Genes")
                            plt.ylabel("Conditions")
                            plt.xticks(rotation=45)
                            plt.tight_layout()
                            plt.savefig(os.path.join(comp_dir, f'{cell_type}_{comp_name}_heatmap.pdf'), dpi=300)
                            plt.close()
                        else:
                            print(f"    热图绘制失败: 分组 {group1} 或 {group2} 不在数据中")
                            
                    except Exception as e:
                        print(f"    热图绘制失败: {str(e)}")
                
                # 富集分析
                if not diff_df_filtered.empty and species:
                    try:
                        genes = diff_df_filtered['names'].tolist()
                        enrichment_dir = os.path.join(comp_dir, f'enrichment_{cell_type}_{comp_name}')
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
                        print(f"    富集分析完成: {enrichment_dir}")
                    except Exception as e:
                        print(f"    富集分析失败: {str(e)}")
                
            except Exception as e:
                print(f"    差异分析失败: {str(e)}")
                import traceback
                traceback.print_exc()
    
    # 保存所有结果
    if all_diff:
        combined_diff = pd.concat(all_diff, ignore_index=True)
        combined_diff.rename(columns={'names': 'Symbol'}, inplace=True)
        combined_diff.to_csv(os.path.join(out_dir, 'all_site_band_diff_genes.xls'), sep='\t', index=False)
        
        # 统计每个比较对的差异基因数量
        stats = combined_diff.groupby(['cell_type', 'comparison', 'group1', 'group2', 'n_group1', 'n_group2']).size().reset_index(name='n_diff_genes')
        stats.to_csv(os.path.join(out_dir, 'diff_gene_statistics.xls'), sep='\t', index=False)
        
        # 创建汇总表格
        summary = pd.pivot_table(
            stats, 
            values='n_diff_genes', 
            index='cell_type', 
            columns='comparison',
            fill_value=0
        )
        summary.to_csv(os.path.join(out_dir, 'diff_gene_summary.xls'), sep='\t')
        
        print(f"\n分析完成！结果保存在: {out_dir}")
        print(f"总计发现差异基因: {len(combined_diff)}")
        print(f"统计信息:\n{summary}")
    else:
        print("没有找到任何差异表达基因")

# 使用示例
adata = sc.read_h5ad("adata_with_band_by_sample.h5ad")

# 运行site_band的差异分析
find_group_diff_genes_site_band(
    adata,
    group_col='annotation',  # 细胞类型分组列
    condition_col='site_band',  # 处理条件列
    comparisons=[('in', 'on'), ('in', 'out'), ('on', 'out')],  # 比较对
    species='human',
    out_dir='./site_band_degs',  # 输出目录
    logfoldchanges=0.5,
    pvals_adj=0.05,
    pct_nz_group=0.1,
    use_raw=True,  # 使用原始矩阵进行差异分析
    plot_volcano=True,
    plot_heatmap=True,
    heatmap_top_genes=20
)
