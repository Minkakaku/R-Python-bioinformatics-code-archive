import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import os
import re
import warnings
warnings.filterwarnings('ignore')
from sklearn.preprocessing import StandardScaler
from scipy import stats

def compare_mucus_zone_deg_within_inner_samples(adata, group_col='mucus_zone', 
                                               sample_col='sample',
                                               species=None, out_dir='./site_zone_deg',
                                               filter=True, logfoldchanges=0.5, pvals_adj=0.05, pct_nz_group=0.1,
                                               method='wilcoxon', 
                                               plot_volcano=True, plot_heatmap=True, plot_violin=True,
                                               volcano_thresholds=None, heatmap_top_genes=20, violin_top_genes=10,
                                               use_raw=True, layer=None, dpi=300, figsize=(12, 8)):
    """
    分析mucus_zone列中inner和outer组间的差异表达基因
    先提取mucus_zone为inner的样本，然后从这些样本中提取所有细胞进行分析
    
    Args:
        adata: AnnData对象
        group_col: 分组列名，默认为'mucus_zone'
        sample_col: 样本列名，默认为'sample'
        species: 物种（'human'或'mouse'）
        out_dir: 输出目录，默认为'./site_zone_deg'
        filter: 是否过滤结果
        logfoldchanges: log2倍数变化阈值
        pvals_adj: 调整后p值阈值
        pct_nz_group: 组内非零表达比例阈值
        method: 差异分析方法
        plot_volcano: 是否绘制火山图
        plot_heatmap: 是否绘制热图
        plot_violin: 是否绘制小提琴图
        volcano_thresholds: 火山图阈值 (logfc_thresh, pval_thresh)
        heatmap_top_genes: 热图展示的top基因数
        violin_top_genes: 小提琴图展示的top基因数
        use_raw: 是否使用原始矩阵进行差异分析
        layer: 如果提供，则使用指定层的数据
        dpi: 图片分辨率
        figsize: 图片大小
    """
    # 确保分组列存在
    if group_col not in adata.obs.columns:
        raise ValueError(f"列 '{group_col}' 不存在")
    
    if sample_col not in adata.obs.columns:
        raise ValueError(f"样本列 '{sample_col}' 不存在")
    
    # 检查两个分组是否存在
    groups = ('inner', 'outer')
    group1, group2 = groups
    
    # 1. 提取mucus_zone为inner的细胞对应的样本
    inner_cells_mask = adata.obs[group_col] == 'inner'
    inner_samples = adata.obs.loc[inner_cells_mask, sample_col].unique()
    
    print(f"找到 {len(inner_samples)} 个包含inner细胞的样本:")
    for sample in inner_samples:
        print(f"  - {sample}")
    
    # 2. 从这些样本中提取所有细胞
    sample_mask = adata.obs[sample_col].isin(inner_samples)
    adata_filtered = adata[sample_mask].copy()
    
    print(f"\n过滤后的数据:")
    print(f"  总细胞数: {adata_filtered.n_obs}")
    print(f"  样本数: {len(inner_samples)}")
    
    # 检查过滤后数据中两个分组是否存在
    group_counts = adata_filtered.obs[group_col].value_counts()
    for group in groups:
        if group not in adata_filtered.obs[group_col].values:
            raise ValueError(f"过滤后，分组 '{group}' 不存在于列 '{group_col}' 中")
    
    # 打印分组分布
    print(f"\n分组分布 (过滤后):")
    for group in groups:
        count = group_counts.get(group, 0)
        print(f"  {group}: {count} 个细胞")
    
    # 检查每个样本中inner和outer的分布
    print(f"\n各样本中inner/outer分布:")
    sample_group_counts = pd.crosstab(adata_filtered.obs[sample_col], adata_filtered.obs[group_col])
    print(sample_group_counts)
    
    # 创建输出目录
    os.makedirs(out_dir, exist_ok=True)
    comp_name = f"{group1}_vs_{group2}"
    comp_dir = os.path.join(out_dir, comp_name)
    os.makedirs(comp_dir, exist_ok=True)
    
    print(f"\n开始分析: {comp_name}")
    print(f"输出目录: {os.path.abspath(comp_dir)}")
    
    # 设置火山图默认阈值
    if volcano_thresholds is None:
        volcano_thresholds = (logfoldchanges, pvals_adj)
    
    # 选择用于差异分析的数据
    if use_raw and adata_filtered.raw is not None:
        # 使用原始矩阵
        print("使用原始raw矩阵进行差异分析")
        raw_adata = adata_filtered.raw.to_adata()
    elif layer is not None and layer in adata_filtered.layers:
        # 使用指定层的数据
        print(f"使用层 '{layer}' 进行差异分析")
        raw_adata = adata_filtered.copy()
        raw_adata.X = raw_adata.layers[layer]
    else:
        # 使用当前矩阵
        print("使用当前矩阵进行差异分析")
        raw_adata = adata_filtered.copy()
    
    # 差异分析
    print(f"执行差异分析: {method}")
    sc.tl.rank_genes_groups(
        raw_adata,
        groupby=group_col,
        groups=[group1],  # 比较组
        reference=group2,  # 参考组
        method=method,
        pts=True
    )
    
    # 获取差异表达结果
    result = raw_adata.uns['rank_genes_groups']
    diff_df = sc.get.rank_genes_groups_df(raw_adata, group=group1)
    
    if diff_df.empty:
        print(f"没有找到差异表达基因")
        return None
    
    # 添加元信息
    diff_df['group1'] = group1
    diff_df['group2'] = group2
    diff_df['n_group1'] = group_counts.get(group1, 0)
    diff_df['n_group2'] = group_counts.get(group2, 0)
    diff_df['analysis_samples'] = len(inner_samples)
    
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
    result_file = os.path.join(comp_dir, f'{comp_name}_diff_genes.xls')
    diff_df_filtered.to_csv(result_file, sep='\t', index=False)
    
    print(f"\n差异基因统计:")
    print(f"  总差异基因数: {len(diff_df)}")
    print(f"  过滤后差异基因数: {len(diff_df_filtered)}")
    print(f"  上调基因数: {sum(diff_df_filtered['logfoldchanges'] > 0)}")
    print(f"  下调基因数: {sum(diff_df_filtered['logfoldchanges'] < 0)}")
    print(f"  结果保存到: {result_file}")
    
    # ================= 火山图绘制 =================
    if plot_volcano and not diff_df_filtered.empty:
        try:
            print("\n绘制火山图...")
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
            
            # 统计
            n_up = sum(volcano_df['type'] == 'up')
            n_down = sum(volcano_df['type'] == 'down')
            n_nosig = sum(volcano_df['type'] == 'nosig')
            
            # 绘制火山图
            plt.figure(figsize=figsize)
            colors = {"down": "#377EB8", "nosig": "grey", "up": "#E41A1C"}
            
            for gene_type, color in colors.items():
                subset = volcano_df[volcano_df['type'] == gene_type]
                if len(subset) > 0:
                    plt.scatter(subset['logfoldchanges'], 
                              subset['-log10(padj)'],
                              c=color, 
                              s=20,
                              alpha=0.7,
                              label=f"{gene_type} ({len(subset)})",
                              edgecolors='none')
            
            # 添加阈值线
            plt.axvline(x=logfc_thresh, color='dimgrey', linestyle='--', linewidth=1.5, alpha=0.8)
            plt.axvline(x=-logfc_thresh, color='dimgrey', linestyle='--', linewidth=1.5, alpha=0.8)
            plt.axhline(y=-np.log10(pval_thresh), color='dimgrey', linestyle='--', linewidth=1.5, alpha=0.8)
            
            # 标记top基因
            if not diff_df_filtered.empty:
                # 获取上调基因中最显著的
                up_genes = diff_df_filtered[diff_df_filtered['logfoldchanges'] > 0]
                down_genes = diff_df_filtered[diff_df_filtered['logfoldchanges'] < 0]
                
                top_up = up_genes.nsmallest(5, 'pvals_adj')['names'].tolist() if not up_genes.empty else []
                top_down = down_genes.nsmallest(5, 'pvals_adj')['names'].tolist() if not down_genes.empty else []
                top_genes = list(set(top_up + top_down))
                
                for _, row in volcano_df[volcano_df['names'].isin(top_genes)].iterrows():
                    plt.text(row['logfoldchanges'], 
                            row['-log10(padj)'] + 0.2,
                            row['names'],
                            fontsize=10,
                            ha='center',
                            fontweight='bold',
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
            
            plt.title(f"Volcano Plot: {comp_name}\n({n_up} up, {n_down} down, |logFC|>{logfc_thresh}, padj<{pval_thresh})", 
                     fontsize=16, fontweight='bold', pad=20)
            plt.xlabel(f"log2({group1}/{group2})", fontsize=14, fontweight='bold')
            plt.ylabel("-log10(Adjusted P-value)", fontsize=14, fontweight='bold')
            plt.legend(loc='best', fontsize=12)
            plt.grid(True, alpha=0.3, linestyle='--')
            
            plt.tight_layout()
            
            # 保存火山图
            volcano_path_png = os.path.join(comp_dir, f'{comp_name}_volcano.png')
            volcano_path_pdf = os.path.join(comp_dir, f'{comp_name}_volcano.pdf')
            plt.savefig(volcano_path_png, dpi=dpi, bbox_inches='tight')
            plt.savefig(volcano_path_pdf, dpi=dpi, bbox_inches='tight')
            print(f"  火山图保存到: {volcano_path_png}")
            
            plt.show()
            plt.close()
            
        except Exception as e:
            print(f"火山图绘制失败: {str(e)}")
    
    # ================= 热图绘制 =================
    if plot_heatmap and not diff_df_filtered.empty and len(diff_df_filtered) >= 5:
        try:
            print("\n绘制热图...")
            # 获取top差异基因
            n_genes = min(heatmap_top_genes, len(diff_df_filtered))
            # 按p值排序，取最显著的基因
            top_genes = diff_df_filtered.nsmallest(n_genes, 'pvals_adj')['names'].tolist()
            
            # 从主数据中提取表达矩阵（使用所有基因）
            comp_adata = adata_filtered
            
            if hasattr(comp_adata.X, 'toarray'):
                expr_data = comp_adata[:, top_genes].X.toarray()
            else:
                expr_data = comp_adata[:, top_genes].X
            
            expr_df = pd.DataFrame(
                data=expr_data,
                index=comp_adata.obs[group_col],
                columns=top_genes
            )
            
            # 按组计算平均表达量
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
                plt.figure(figsize=(max(12, len(top_genes)*0.5), 6))
                sns.heatmap(
                    scaled_expr,
                    cmap='coolwarm',
                    center=0,
                    annot=True,
                    fmt=".2f",
                    linewidths=0.5,
                    linecolor='white',
                    cbar_kws={'label': 'Z-score', 'shrink': 0.8},
                    square=True
                )
                
                plt.title(f"Top {n_genes} DE Genes: {comp_name}\n(仅包含有inner细胞的样本)", fontsize=16, fontweight='bold', pad=20)
                plt.xlabel("Genes", fontsize=14, fontweight='bold')
                plt.ylabel("Mucus Zone", fontsize=14, fontweight='bold')
                plt.xticks(rotation=45, ha='right', fontsize=10)
                plt.yticks(rotation=0, fontsize=12)
                
                plt.tight_layout()
                
                # 保存热图
                heatmap_path_png = os.path.join(comp_dir, f'{comp_name}_heatmap.png')
                heatmap_path_pdf = os.path.join(comp_dir, f'{comp_name}_heatmap.pdf')
                plt.savefig(heatmap_path_png, dpi=dpi, bbox_inches='tight')
                plt.savefig(heatmap_path_pdf, dpi=dpi, bbox_inches='tight')
                print(f"  热图保存到: {heatmap_path_png}")
                
                plt.show()
                plt.close()
            else:
                print(f"  热图绘制失败: 分组 {group1} 或 {group2} 不在数据中")
                
        except Exception as e:
            print(f"热图绘制失败: {str(e)}")
    
    # ================= 小提琴图绘制 =================
    if plot_violin and not diff_df_filtered.empty:
        try:
            print("\n绘制小提琴图...")
            # 获取top差异基因
            n_genes = min(violin_top_genes, len(diff_df_filtered))
            # 按p值排序，取最显著的基因
            top_genes = diff_df_filtered.nsmallest(n_genes, 'pvals_adj')['names'].tolist()
            
            # 从主数据中提取表达矩阵
            comp_adata = adata_filtered
            
            if hasattr(comp_adata.X, 'toarray'):
                expr_data = comp_adata[:, top_genes].X.toarray()
            else:
                expr_data = comp_adata[:, top_genes].X
            
            expr_df = pd.DataFrame(
                data=expr_data,
                index=comp_adata.obs_names,
                columns=top_genes
            )
            
            # 添加分组信息
            expr_df['group'] = comp_adata.obs[group_col].values
            
            # 转换为长格式
            expr_long = expr_df.melt(id_vars=['group'], var_name='gene', value_name='expression')
            
            # 设置颜色
            palette = {group1: '#1f77b4', group2: '#d62728'}
            
            # 创建多个子图
            n_cols = 3
            n_rows = int(np.ceil(n_genes / n_cols))
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
            
            # 处理axes的形状
            if n_rows == 1 and n_cols == 1:
                axes = np.array([axes])
            elif n_rows == 1 or n_cols == 1:
                axes = axes.flatten()
            else:
                axes = axes.flatten()
            
            for idx, gene in enumerate(top_genes):
                if idx >= len(axes):
                    break
                    
                ax = axes[idx]
                gene_data = expr_long[expr_long['gene'] == gene]
                
                # 绘制小提琴图
                sns.violinplot(x='group', y='expression', data=gene_data, 
                              order=[group1, group2], palette=palette, ax=ax,
                              inner='box', linewidth=1.5, cut=0)
                
                # 计算统计检验
                groups_data = {group: gene_data[gene_data['group'] == group]['expression'].values 
                              for group in [group1, group2]}
                
                # 检查每组是否有数据
                if len(groups_data[group1]) > 0 and len(groups_data[group2]) > 0:
                    # Mann-Whitney U检验
                    try:
                        mw_stat, mw_p = stats.mannwhitneyu(
                            groups_data[group1], 
                            groups_data[group2],
                            alternative='two-sided'
                        )
                        
                        # 在图上添加p值
                        if mw_p < 0.001:
                            p_text = 'p < 0.001 ***'
                        elif mw_p < 0.01:
                            p_text = 'p < 0.01 **'
                        elif mw_p < 0.05:
                            p_text = 'p < 0.05 *'
                        else:
                            p_text = f'p = {mw_p:.3f}'
                        
                        # 添加统计结果
                        ax.text(0.5, 0.95, p_text, transform=ax.transAxes, 
                                ha='center', va='top', fontsize=10,
                                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
                    except Exception as e:
                        print(f"    基因 {gene} 统计检验失败: {e}")
                
                # 设置标题
                ax.set_title(f'{gene}', fontsize=12, fontweight='bold')
                ax.set_xlabel('')
                ax.set_ylabel('Expression', fontsize=10)
                
                # 旋转x轴标签
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            
            # 隐藏多余的子图
            for idx in range(len(top_genes), len(axes)):
                axes[idx].set_visible(False)
            
            plt.suptitle(f'Top {n_genes} DE Genes: {comp_name}\n(仅包含有inner细胞的样本)', fontsize=16, fontweight='bold')
            plt.tight_layout()
            
            # 保存小提琴图
            violin_path_png = os.path.join(comp_dir, f'{comp_name}_violinplots.png')
            violin_path_pdf = os.path.join(comp_dir, f'{comp_name}_violinplots.pdf')
            plt.savefig(violin_path_png, dpi=dpi, bbox_inches='tight')
            plt.savefig(violin_path_pdf, dpi=dpi, bbox_inches='tight')
            print(f"  小提琴图保存到: {violin_path_png}")
            
            plt.show()
            plt.close(fig)
            
        except Exception as e:
            print(f"小提琴图绘制失败: {str(e)}")
    
    # ================= 富集分析 =================
    if not diff_df_filtered.empty and species:
        try:
            print(f"\n进行富集分析...")
            genes = diff_df_filtered['names'].tolist()
            enrichment_dir = os.path.join(comp_dir, f'enrichment_analysis')
            os.makedirs(enrichment_dir, exist_ok=True)
            
            if species == 'human':
                enr = gp.enrichr(
                    gene_list=genes,
                    gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human', 'Reactome_2022'],
                    organism='Human',
                    cutoff=0.05,
                    outdir=enrichment_dir
                )
            elif species == 'mouse':
                enr = gp.enrichr(
                    gene_list=genes,
                    gene_sets=['GO_Biological_Process_2023', 'KEGG_2019_Mouse', 'Reactome_2022'],
                    organism='Mouse',
                    cutoff=0.05,
                    outdir=enrichment_dir
                )
            print(f"  富集分析完成: {enrichment_dir}")
            
            # 绘制富集分析结果
            if enr is not None and hasattr(enr, 'results'):
                # 保存富集结果
                for key, df in enr.results.items():
                    if not df.empty:
                        df.to_csv(os.path.join(enrichment_dir, f'{key}_enrichment.xls'), sep='\t', index=False)
                        
                        # 绘制top富集项
                        top_terms = df.nsmallest(10, 'Adjusted P-value')
                        if not top_terms.empty:
                            plt.figure(figsize=(12, 8))
                            plt.barh(range(len(top_terms)), -np.log10(top_terms['Adjusted P-value']))
                            plt.yticks(range(len(top_terms)), top_terms['Term'])
                            plt.xlabel('-log10(Adjusted P-value)', fontsize=12, fontweight='bold')
                            plt.title(f'Top 10 Enriched Terms: {key}', fontsize=14, fontweight='bold')
                            plt.gca().invert_yaxis()
                            plt.tight_layout()
                            
                            enrich_plot_path = os.path.join(enrichment_dir, f'{key}_enrichment.png')
                            plt.savefig(enrich_plot_path, dpi=dpi, bbox_inches='tight')
                            plt.close()
            
        except Exception as e:
            print(f"富集分析失败: {str(e)}")
    
    # ================= 创建分析报告 =================
    print(f"\n创建分析报告...")
    report_file = os.path.join(comp_dir, 'analysis_report.txt')
    with open(report_file, 'w') as f:
        f.write("="*60 + "\n")
        f.write(f"Mucus Zone Differential Expression Analysis Report\n")
        f.write("="*60 + "\n\n")
        
        f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Comparison: {group1} vs {group2}\n")
        f.write(f"Analysis Strategy: Only samples with inner cells were included\n")
        f.write(f"Total Samples Included: {len(inner_samples)}\n")
        f.write(f"Total Cells After Filtering: {adata_filtered.n_obs}\n")
        f.write(f"Group Distribution:\n")
        f.write(f"  {group1}: {group_counts.get(group1, 0)} cells\n")
        f.write(f"  {group2}: {group_counts.get(group2, 0)} cells\n")
        
        f.write(f"\nSamples Included in Analysis:\n")
        for sample in inner_samples:
            f.write(f"  - {sample}\n")
        
        f.write(f"\nDifferential Expression Analysis:\n")
        f.write(f"  Method: {method}\n")
        f.write(f"  Thresholds: |logFC| > {logfoldchanges}, padj < {pvals_adj}, pct_nz >= {pct_nz_group}\n")
        f.write(f"  Total DE genes: {len(diff_df)}\n")
        f.write(f"  Filtered DE genes: {len(diff_df_filtered)}\n")
        f.write(f"  Up-regulated in {group1}: {sum(diff_df_filtered['logfoldchanges'] > 0)}\n")
        f.write(f"  Down-regulated in {group1}: {sum(diff_df_filtered['logfoldchanges'] < 0)}\n")
        
        f.write(f"\nTop 10 Up-regulated genes in {group1}:\n")
        top_up = diff_df_filtered[diff_df_filtered['logfoldchanges'] > 0].nsmallest(10, 'pvals_adj')
        for _, row in top_up.iterrows():
            f.write(f"  {row['names']}: logFC={row['logfoldchanges']:.3f}, padj={row['pvals_adj']:.6f}\n")
        
        f.write(f"\nTop 10 Down-regulated genes in {group1}:\n")
        top_down = diff_df_filtered[diff_df_filtered['logfoldchanges'] < 0].nsmallest(10, 'pvals_adj')
        for _, row in top_down.iterrows():
            f.write(f"  {row['names']}: logFC={row['logfoldchanges']:.3f}, padj={row['pvals_adj']:.6f}\n")
        
        f.write(f"\nFiles Generated:\n")
        for file_type in ['diff_genes.xls', 'volcano.png/pdf', 'heatmap.png/pdf', 'violinplots.png/pdf']:
            f.write(f"  • {comp_name}_{file_type}\n")
        if species:
            f.write(f"  • enrichment_analysis/ directory\n")
        f.write(f"  • analysis_report.txt\n")
    
    print(f"分析报告保存到: {report_file}")
    print(f"\n分析完成！所有结果已保存到: {comp_dir}")
    
    return diff_df_filtered, adata_filtered

# 主程序
def main():
    # 加载数据
    print("加载数据...")
    adata = sc.read_h5ad("adata_with_band_by_sample.h5ad")
    
    print(f"数据信息:")
    print(f"  细胞数: {adata.n_obs}")
    print(f"  基因数: {adata.n_vars}")
    print(f"  观察值列: {adata.obs.columns.tolist()}")
    
    # 检查mucus_zone列
    if 'mucus_zone' not in adata.obs.columns:
        print("错误: 数据中未找到'mucus_zone'列")
        return
    
    # 检查sample列（根据实际情况修改列名）
    sample_column = 'sample'  # 请根据您的数据修改这个列名
    if sample_column not in adata.obs.columns:
        print(f"错误: 数据中未找到样本列 '{sample_column}'")
        print(f"可用的列: {adata.obs.columns.tolist()}")
        # 让用户输入样本列名
        sample_column = input("请输入样本列名: ")
        if sample_column not in adata.obs.columns:
            print(f"错误: 列 '{sample_column}' 不存在")
            return
    
    print(f"\nmucus_zone分布:")
    print(adata.obs['mucus_zone'].value_counts())
    
    print(f"\n样本数: {adata.obs[sample_column].nunique()}")
    
    # 运行差异表达分析
    result, filtered_adata = compare_mucus_zone_deg_within_inner_samples(
        adata,
        group_col='mucus_zone',
        sample_col=sample_column,  # 指定样本列
        species='human',  # 根据您的数据修改物种
        out_dir='./site_zone_deg',
        logfoldchanges=0.5,
        pvals_adj=0.05,
        pct_nz_group=0.1,
        method='wilcoxon',
        plot_volcano=True,
        plot_heatmap=True,
        plot_violin=True,
        heatmap_top_genes=20,
        violin_top_genes=12,
        use_raw=True,  # 使用原始raw矩阵
        dpi=300,
        figsize=(12, 8)
    )
    
    if result is not None:
        print(f"\n分析完成！")
        print(f"发现 {len(result)} 个差异表达基因")
        print(f"结果保存在: {os.path.abspath('./site_zone_deg/inner_vs_outer/')}")
        
        # 保存过滤后的数据
        filtered_data_path = './site_zone_deg/filtered_data_for_analysis.h5ad'
        filtered_adata.write(filtered_data_path)
        print(f"过滤后的分析数据保存到: {filtered_data_path}")

# 运行主程序
if __name__ == "__main__":
    main()
