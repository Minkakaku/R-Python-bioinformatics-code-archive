import scanpy as sc
import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib import pyplot as plt

# 设置中文字体（可选）
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False   # 用来正常显示负号

def perform_de_analysis(adata, groupby, method='wilcoxon', key_added=None):
    """
    执行差异表达分析
    """
    if key_added is None:
        key_added = f'rank_genes_groups_{groupby}'
    
    # 执行差异表达分析
    sc.tl.rank_genes_groups(
        adata, 
        groupby=groupby, 
        method=method,
        key_added=key_added,
        pts=True
    )
    
    return key_added

def extract_de_results(adata, key, output_dir='de_results'):
    """
    提取差异表达分析结果并保存到文件
    """
    os.makedirs(output_dir, exist_ok=True)
    
    groups = adata.uns[key]['names'].dtype.names
    all_results = []
    
    for group in groups:
        de_df = pd.DataFrame({
            'gene': adata.uns[key]['names'][group],
            'scores': adata.uns[key]['scores'][group],
            'logfoldchanges': adata.uns[key]['logfoldchanges'][group],
            'pvals': adata.uns[key]['pvals'][group],
            'pvals_adj': adata.uns[key]['pvals_adj'][group],
            'pct_expr': adata.uns[key]['pts'][group] if 'pts' in adata.uns[key] else np.nan
        })
        
        de_df['group'] = group
        significant_genes = de_df[
            (de_df['pvals_adj'] < 0.05) & 
            (abs(de_df['logfoldchanges']) > 0.5)
        ]
        
        # 保存当前组的完整结果
        de_filename = f'{output_dir}/de_{key}_{group}_all.csv'
        de_df.to_csv(de_filename, index=False)
        
        # 保存显著差异基因
        if len(significant_genes) > 0:
            sig_filename = f'{output_dir}/de_{key}_{group}_significant.csv'
            significant_genes.to_csv(sig_filename, index=False)
        
        all_results.append(de_df)
        print(f"分组 {group}: 总基因数 {len(de_df)}, 显著差异基因数 {len(significant_genes)}")
    
    # 保存所有组的综合结果
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        combined_filename = f'{output_dir}/de_{key}_combined.csv'
        combined_df.to_csv(combined_filename, index=False)
    
    return all_results

def perform_enrichment_analysis(deg_df, group_name, output_dir='enrichment_results', 
                               organism='human', gene_sets=None):
    """
    执行GO和KEGG富集分析[1,6](@ref)
    
    参数:
    deg_df: 差异表达基因DataFrame
    group_name: 分组名称
    output_dir: 结果保存目录
    organism: 物种，默认human
    gene_sets: 基因集列表，默认GO和KEGG
    """
    if gene_sets is None:
        gene_sets = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
    
    # 创建富集分析结果目录
    enrichment_dir = f'{output_dir}/{group_name}'
    os.makedirs(enrichment_dir, exist_ok=True)
    os.makedirs(f'{enrichment_dir}/plots', exist_ok=True)
    
    # 提取显著差异基因列表[6](@ref)
    significant_genes = deg_df[
        (deg_df['pvals_adj'] < 0.05) & 
        (abs(deg_df['logfoldchanges']) > 0.5)
    ]
    
    if len(significant_genes) < 5:
        print(f"分组 {group_name} 的显著差异基因数不足5个，跳过富集分析")
        return None
    
    # 取差异最显著的前100个基因[6](@ref)
    gene_list = significant_genes.sort_values('pvals_adj')['gene'].head(100).tolist()
    
    enrichment_results = {}
    
    try:
        # 执行富集分析[6](@ref)
        enr_results = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism,
            outdir=None,
            cutoff=0.05
        )
        
        # 保存富集结果
        if enr_results is not None:
            enrichment_results['enrichr'] = enr_results
            
            # 保存结果表格
            for gene_set in gene_sets:
                if gene_set in enr_results.results:
                    result_df = enr_results.results[gene_set]
                    filename = f'{enrichment_dir}/{group_name}_{gene_set}_enrichment.csv'
                    result_df.to_csv(filename, index=False)
                    
                    # 可视化前10个富集项[6](@ref)
                    if len(result_df) > 0:
                        plot_enrichment_results(result_df, gene_set, group_name, enrichment_dir)
        
        # 执行GSEA分析（基于排名）[3,4](@ref)
        print(f"执行GSEA分析 for {group_name}...")
        try:
            # 准备排名数据（基于log2 fold change）
            rank_data = deg_df.set_index('gene')['logfoldchanges'].to_dict()
            
            gsea_results = gp.gsea(
                data=rank_data,
                gene_sets=gene_sets,
                organism=organism,
                outdir=f'{enrichment_dir}/gsea',
                permutation_num=1000
            )
            
            if gsea_results is not None:
                enrichment_results['gsea'] = gsea_results
                
                # 绘制GSEA图
                plot_gsea_results(gsea_results, group_name, enrichment_dir)
                
        except Exception as e:
            print(f"GSEA分析失败: {str(e)}")
        
        return enrichment_results
        
    except Exception as e:
        print(f"富集分析出错: {str(e)}")
        return None

def plot_enrichment_results(result_df, gene_set, group_name, output_dir, top_n=10):
    """
    绘制富集分析结果图[1,6](@ref)
    """
    try:
        # 筛选显著富集的结果
        sig_results = result_df[result_df['Adjusted P-value'] < 0.05].head(top_n)
        
        if len(sig_results) == 0:
            print(f"没有显著富集的{gene_set}条目")
            return
        
        # 创建子图
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        
        # 气泡图
        plt.sca(axes[0])
        # 准备数据
        sizes = np.log10(sig_results['Odds Ratio'] + 1) * 100
        colors = -np.log10(sig_results['Adjusted P-value'] + 1e-10)
        
        scatter = plt.scatter(range(len(sig_results)), 
                             sig_results['Odds Ratio'],
                             s=sizes, c=colors, cmap='viridis', alpha=0.7)
        
        plt.colorbar(scatter, label='-log10(Adj P-value)')
        plt.xticks(range(len(sig_results)), sig_results['Term'], rotation=45, ha='right')
        plt.xlabel('Term')
        plt.ylabel('Odds Ratio')
        plt.title(f'{gene_set} - {group_name}\nBubble Plot')
        plt.tight_layout()
        
        # 条形图
        plt.sca(axes[1])
        y_pos = np.arange(len(sig_results))
        plt.barh(y_pos, -np.log10(sig_results['Adjusted P-value'] + 1e-10))
        plt.yticks(y_pos, sig_results['Term'])
        plt.xlabel('-log10(Adjusted P-value)')
        plt.title(f'{gene_set} - {group_name}\nBar Plot')
        plt.tight_layout()
        
        # 保存图片
        plot_filename = f'{output_dir}/plots/{group_name}_{gene_set}_enrichment.png'
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"富集分析图已保存: {plot_filename}")
        
    except Exception as e:
        print(f"绘制富集分析图时出错: {str(e)}")

def plot_gsea_results(gsea_results, group_name, output_dir):
    """
    绘制GSEA分析结果[3,4](@ref)
    """
    try:
        # 获取最显著的前3个通路
        top_pathways = gsea_results.results_onplot.head(3)
        
        for i, (pathway, result) in enumerate(top_pathways.iterrows()):
            try:
                # 绘制GSEA图
                plt.figure(figsize=(10, 6))
                gp.gsea_plot(gsea_results.ranking, pathway, 
                            ofile=f'{output_dir}/plots/{group_name}_gsea_{pathway}.png',
                            **result)
                plt.close()
            except Exception as e:
                print(f"绘制GSEA图 {pathway} 时出错: {str(e)}")
                
    except Exception as e:
        print(f"处理GSEA结果时出错: {str(e)}")

def create_volcano_plot(adata, key, group, output_dir='de_results'):
    """
    创建火山图可视化差异表达结果
    """
    try:
        de_df = pd.DataFrame({
            'gene': adata.uns[key]['names'][group],
            'log2fc': adata.uns[key]['logfoldchanges'][group],
            'pvals_adj': adata.uns[key]['pvals_adj'][group]
        })
        
        de_df['neg_log10_pval'] = -np.log10(de_df['pvals_adj'])
        de_df['significant'] = (de_df['pvals_adj'] < 0.05) & (abs(de_df['log2fc']) > 0.5)
        
        plt.figure(figsize=(10, 8))
        
        nonsig = de_df[~de_df['significant']]
        plt.scatter(nonsig['log2fc'], nonsig['neg_log10_pval'], 
                   alpha=0.5, s=10, color='grey', label='非显著')
        
        sig = de_df[de_df['significant']]
        plt.scatter(sig['log2fc'], sig['neg_log10_pval'], 
                   alpha=0.7, s=15, color='red', label='显著')
        
        top_genes = sig.nlargest(10, 'neg_log10_pval')
        for _, gene in top_genes.iterrows():
            plt.annotate(gene['gene'], 
                        (gene['log2fc'], gene['neg_log10_pval']),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8, alpha=0.8)
        
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.8)
        plt.axvline(x=0, color='black', linestyle='-', alpha=0.8)
        plt.xlabel('Log2 Fold Change')
        plt.ylabel('-Log10(调整后P值)')
        plt.title(f'火山图 - {key} - 分组 {group}')
        plt.legend()
        plt.tight_layout()
        
        plot_filename = f'{output_dir}/volcano_{key}_{group}.png'
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"火山图已保存: {plot_filename}")
        
    except Exception as e:
        print(f"创建火山图时出错: {str(e)}")

def main():
    # 加载h5ad数据
    print("正在加载h5ad数据...")
    adata = sc.read_h5ad('epi_with_g1_g2.h5ad')  # 请替换为你的h5ad文件路径
    
    print(f"数据形状: {adata.shape}")
    print(f"观察值(细胞)注解列: {adata.obs.columns.tolist()}")
    
    # 检查g1和g2列是否存在
    if 'g1' not in adata.obs.columns:
        print("错误: 在obs中未找到g1列")
        return
    
    if 'g2' not in adata.obs.columns:
        print("错误: 在obs中未找到g2列")
        return
    
    print(f"g1分组类别: {adata.obs['g1'].unique().tolist()}")
    print(f"g2分组类别: {adata.obs['g2'].unique().tolist()}")
    
    # 对g1分组进行差异表达分析
    print("\n正在对g1分组进行差异表达分析...")
    key_g1 = perform_de_analysis(adata, 'g1')
    results_g1 = extract_de_results(adata, key_g1, 'de_results_g1')
    
    # 对g2分组进行差异表达分析
    print("\n正在对g2分组进行差异表达分析...")
    key_g2 = perform_de_analysis(adata, 'g2')
    results_g2 = extract_de_results(adata, key_g2, 'de_results_g2')
    
    # 执行富集分析[1,6](@ref)
    print("\n正在执行富集分析...")
    
    # 为g1分组的每个比较执行富集分析
    groups_g1 = adata.uns[key_g1]['names'].dtype.names
    for group in groups_g1:
        print(f"处理g1分组 - {group}...")
        de_df = pd.DataFrame({
            'gene': adata.uns[key_g1]['names'][group],
            'logfoldchanges': adata.uns[key_g1]['logfoldchanges'][group],
            'pvals_adj': adata.uns[key_g1]['pvals_adj'][group]
        })
        perform_enrichment_analysis(de_df, f"g1_{group}", 'enrichment_results_g1')
    
    # 为g2分组的每个比较执行富集分析
    groups_g2 = adata.uns[key_g2]['names'].dtype.names
    for group in groups_g2:
        print(f"处理g2分组 - {group}...")
        de_df = pd.DataFrame({
            'gene': adata.uns[key_g2]['names'][group],
            'logfoldchanges': adata.uns[key_g2]['logfoldchanges'][group],
            'pvals_adj': adata.uns[key_g2]['pvals_adj'][group]
        })
        perform_enrichment_analysis(de_df, f"g2_{group}", 'enrichment_results_g2')
    
    # 创建火山图
    print("\n正在创建火山图...")
    for group in groups_g1:
        create_volcano_plot(adata, key_g1, group, 'de_results_g1')
    
    for group in groups_g2:
        create_volcano_plot(adata, key_g2, group, 'de_results_g2')
    
    print("\n分析完成！")
    print(f"差异表达结果保存在: de_results_g1/ 和 de_results_g2/")
    print(f"富集分析结果保存在: enrichment_results_g1/ 和 enrichment_results_g2/")
    
    # 打印摘要信息
    print("\n=== 分析摘要 ===")
    print(f"g1分组: {list(adata.obs['g1'].unique())}")
    print(f"g2分组: {list(adata.obs['g2'].unique())}")
    print(f"总基因数: {adata.n_vars}")
    print(f"总细胞数: {adata.n_obs}")

if __name__ == "__main__":
    main()
