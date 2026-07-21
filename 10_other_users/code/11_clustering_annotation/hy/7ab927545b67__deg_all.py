#!/usr/bin/env python3
"""
整体差异表达分析 + GO/KEGG 在线富集 (基于 state 分组)
输入文件：All_filter_state_samp_new_with_CNs10.h5ad
"""

import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# 导入 gseapy 用于在线富集
try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    print("gseapy 未安装，请运行: pip install gseapy")
    GSEAPY_AVAILABLE = False

# -------------------- 基础处理函数 --------------------
def save_raw_counts(adata):
    """将原始矩阵保存到 .raw"""
    adata.raw = adata.copy()
    print("原始计数已保存到 adata.raw")
    return adata

def normalize_data(adata, target_sum=1e4):
    """标准化：CPM + log1p"""
    print("数据标准化中 (target_sum=1e4 + log1p)...")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=target_sum)
    sc.pp.log1p(adata_norm)
    adata.layers['log1p_normalized'] = adata_norm.X.copy()
    print("标准化完成")
    return adata

# -------------------- 差异表达分析 --------------------
def run_dea(adata, group1, group2, group_col='state', method='wilcoxon'):
    """
    对两个组进行差异表达分析，返回结果 DataFrame
    """
    # 提取子集
    adata_sub = adata[adata.obs[group_col].isin([group1, group2])].copy()
    # 使用原始计数重新标准化
    adata_sub.X = adata_sub.raw.X
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)
    
    sc.tl.rank_genes_groups(
        adata_sub,
        groupby=group_col,
        groups=[group1],
        reference=group2,
        method=method,
        use_raw=False,
        n_genes=adata_sub.shape[1]  # 全部基因
    )
    de_results = sc.get.rank_genes_groups_df(adata_sub, group=group1)
    de_results['comparison'] = f"{group1}_vs_{group2}"
    return de_results

def plot_volcano(de_results, output_dir, comparison_name,
                 logfc_thresh=0.5, pval_thresh=0.05, top_n=10):
    """绘制火山图并保存"""
    if de_results is None or len(de_results) == 0:
        return
    fig, ax = plt.subplots(figsize=(10, 8))
    plot_data = de_results.copy()
    plot_data['neg_log10_pval'] = -np.log10(plot_data['pvals_adj'])
    # 处理无穷值
    max_val = plot_data['neg_log10_pval'].replace([np.inf, -np.inf], np.nan).max()
    if pd.isna(max_val):
        max_val = 10
    plot_data['neg_log10_pval'] = plot_data['neg_log10_pval'].replace([np.inf, -np.inf], max_val)
    
    # 标记上下调
    plot_data['regulation'] = 'Not significant'
    up = (plot_data['logfoldchanges'] > logfc_thresh) & (plot_data['pvals_adj'] < pval_thresh)
    down = (plot_data['logfoldchanges'] < -logfc_thresh) & (plot_data['pvals_adj'] < pval_thresh)
    plot_data.loc[up, 'regulation'] = f'Up (logFC > {logfc_thresh})'
    plot_data.loc[down, 'regulation'] = f'Down (logFC < -{logfc_thresh})'
    
    colors = {
        'Not significant': '#d3d3d3',
        f'Up (logFC > {logfc_thresh})': '#e74c3c',
        f'Down (logFC < -{logfc_thresh})': '#3498db'
    }
    for reg, col in colors.items():
        mask = plot_data['regulation'] == reg
        if mask.sum():
            ax.scatter(plot_data.loc[mask, 'logfoldchanges'],
                       plot_data.loc[mask, 'neg_log10_pval'],
                       c=col, s=8, alpha=0.6, label=reg, edgecolor='none')
    
    ax.axvline(x=logfc_thresh, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=-logfc_thresh, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=-np.log10(pval_thresh), color='gray', linestyle='--', alpha=0.5)
    
    # 标注 Top N 显著基因
    sig = plot_data[plot_data['regulation'] != 'Not significant'].copy()
    if len(sig) > 0:
        sig = sig.sort_values('neg_log10_pval', ascending=False)
        for _, row in sig.head(top_n * 2).iterrows():
            ax.text(row['logfoldchanges'], row['neg_log10_pval'], row['names'],
                    fontsize=8, bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.7))
    
    ax.set_xlabel('log2(Fold Change)')
    ax.set_ylabel('-log10(Adjusted P-value)')
    ax.set_title(f'{comparison_name} Volcano Plot')
    ax.legend()
    ax.grid(alpha=0.3, linestyle='--')
    
    n_up = (plot_data['regulation'] == f'Up (logFC > {logfc_thresh})').sum()
    n_down = (plot_data['regulation'] == f'Down (logFC < -{logfc_thresh})').sum()
    ax.text(0.02, 0.98, f"Significant: {n_up+n_down}\nUp: {n_up}, Down: {n_down}",
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle="round", facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    for ext in ['.png', '.pdf']:
        plt.savefig(os.path.join(output_dir, f"volcano_{comparison_name}{ext}"), dpi=300)
    plt.close()
    print(f"  火山图已保存: volcano_{comparison_name}.png/pdf")

# -------------------- 在线富集与绘图 --------------------
def run_online_enrichment(gene_list, gene_set, output_dir, organism='Human'):
    """调用 Enrichr 在线富集，结果保存到本地文件"""
    if not GSEAPY_AVAILABLE:
        print("    gseapy 不可用，跳过在线富集")
        return None
    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=[gene_set],
            organism=organism,
            outdir=output_dir,
            cutoff=0.05
        )
        expected_file = os.path.join(output_dir, f"{gene_set}.Human.enrichr.reports.txt")
        if os.path.exists(expected_file):
            print(f"    在线富集完成: {expected_file}")
            return expected_file
        else:
            print(f"    未找到预期结果文件: {expected_file}")
            return None
    except Exception as e:
        print(f"    在线富集失败: {e}")
        return None

def plot_enrichment_from_local(result_file, gene_set_name, output_dir, comparison_name, top_n=15):
    """读取本地富集结果文件，绘制条形图、气泡图、点图"""
    if not os.path.exists(result_file):
        print(f"    结果文件不存在: {result_file}")
        return
    try:
        df = pd.read_csv(result_file, sep='\t')
        # 检查并重命名列（Enrichr 输出常见格式）
        if 'Term' not in df.columns or 'Adjusted P-value' not in df.columns:
            if len(df.columns) >= 5:
                df = df.rename(columns={df.columns[0]: 'Term', df.columns[4]: 'Adjusted P-value'})
            else:
                print(f"    文件列名不匹配: {df.columns.tolist()}")
                return
        df['Adjusted P-value'] = pd.to_numeric(df['Adjusted P-value'], errors='coerce')
        df = df.dropna(subset=['Adjusted P-value'])
        if len(df) == 0:
            return
        
        top_terms = df.sort_values('Adjusted P-value').head(top_n).copy()
        if len(top_terms) < 2:
            return
        top_terms['neg_log10'] = -np.log10(top_terms['Adjusted P-value'])
        
        # ---- 条形图 ----
        fig, ax = plt.subplots(figsize=(10, 8))
        top_sorted = top_terms.sort_values('Adjusted P-value', ascending=True)
        colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(top_sorted)))
        ax.barh(range(len(top_sorted)), top_sorted['neg_log10'], color=colors)
        ax.set_yticks(range(len(top_sorted)))
        ax.set_yticklabels([t[:60]+'...' if len(t)>60 else t for t in top_sorted['Term']])
        ax.set_xlabel('-log10(Adjusted P-value)')
        ax.set_title(f"{comparison_name}\n{gene_set_name} Enrichment")
        ax.grid(axis='x', alpha=0.3)
        if 'Overlap' in top_sorted.columns:
            for i, row in enumerate(top_sorted.itertuples()):
                if isinstance(row.Overlap, str) and '/' in row.Overlap:
                    n = row.Overlap.split('/')[0]
                    ax.text(row.neg_log10+0.1, i, f'n={n}', va='center', fontsize=9)
        plt.tight_layout()
        bar_file = os.path.join(output_dir, f"enrichment_{gene_set_name}_bar.png")
        plt.savefig(bar_file, dpi=300)
        plt.savefig(bar_file.replace('.png','.pdf'), dpi=300)
        plt.close()
        print(f"      条形图已保存: {os.path.basename(bar_file)}")
        
        # ---- 气泡图 ----
        if len(top_terms) >= 3:
            fig, ax = plt.subplots(figsize=(10, 8))
            if 'Overlap' in top_terms.columns:
                gene_counts = [int(ov.split('/')[0]) if isinstance(ov,str) and '/' in ov else 1 for ov in top_terms['Overlap']]
            else:
                gene_counts = [1]*len(top_terms)
            sizes = np.array(gene_counts)*20
            scat = ax.scatter(range(len(top_terms)), top_terms['neg_log10'],
                              s=sizes, c=top_terms['neg_log10'], cmap='viridis', alpha=0.7)
            ax.set_xticks(range(len(top_terms)))
            ax.set_xticklabels([t[:40]+'...' if len(t)>40 else t for t in top_terms['Term']],
                               rotation=45, ha='right')
            ax.set_ylabel('-log10(Adjusted P-value)')
            ax.set_title(f"{gene_set_name} Bubble Plot")
            plt.colorbar(scat, label='-log10(P)')
            ax.grid(alpha=0.3)
            plt.tight_layout()
            bubble_file = os.path.join(output_dir, f"enrichment_{gene_set_name}_bubble.png")
            plt.savefig(bubble_file, dpi=300)
            plt.savefig(bubble_file.replace('.png','.pdf'), dpi=300)
            plt.close()
            print(f"      气泡图已保存: {os.path.basename(bubble_file)}")
            
            # ---- 点图 ----
            fig, ax = plt.subplots(figsize=(10, 8))
            gene_ratios = []
            for ov in top_terms['Overlap']:
                if isinstance(ov, str) and '/' in ov:
                    a, b = ov.split('/')
                    try:
                        gene_ratios.append(int(a)/int(b))
                    except:
                        gene_ratios.append(0.01)
                else:
                    gene_ratios.append(0.01)
            scat = ax.scatter(gene_ratios, top_terms['neg_log10'], s=100,
                              c=top_terms['neg_log10'], cmap='plasma', alpha=0.7)
            for i, (x, y) in enumerate(zip(gene_ratios, top_terms['neg_log10'])):
                ax.text(x, y, f'  {i+1}', fontsize=9, va='center')
            ax.set_xlabel('Gene Ratio')
            ax.set_ylabel('-log10(Adjusted P-value)')
            ax.set_title(f"{gene_set_name} Dot Plot")
            plt.colorbar(scat, label='-log10(P)')
            ax.grid(alpha=0.3)
            plt.tight_layout()
            dot_file = os.path.join(output_dir, f"enrichment_{gene_set_name}_dot.png")
            plt.savefig(dot_file, dpi=300)
            plt.savefig(dot_file.replace('.png','.pdf'), dpi=300)
            plt.close()
            print(f"      点图已保存: {os.path.basename(dot_file)}")
    except Exception as e:
        print(f"    绘图失败: {e}")

# -------------------- 主流程 --------------------
def main():
    # 输入文件（请确认路径）
    input_file = "All_filter_state_samp_new_with_CNs10.h5ad"
    output_dir = "./DEA_State_Enrichment_all"
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 60)
    print("整体差异表达 + GO/KEGG 富集分析 (按 state 分组)")
    print("=" * 60)
    
    # 1. 读取数据
    print("\n1. 读取数据...")
    adata = sc.read_h5ad(input_file)
    print(f"   数据维度: {adata.n_obs} 细胞, {adata.n_vars} 基因")
    if 'state' not in adata.obs.columns:
        raise KeyError("obs 中缺少 'state' 列")
    
    # 2. 保存原始计数并标准化
    print("\n2. 保存原始计数...")
    adata = save_raw_counts(adata)
    print("\n3. 数据标准化...")
    adata = normalize_data(adata)
    
    # 3. 获取分组信息
    groups = sorted([g for g in adata.obs['state'].unique() if str(g) != 'nan'])
    print(f"\n4. 发现的分组: {groups}")
    comparisons = list(combinations(groups, 2))
    
    # 4. 对每个比较对进行分析
    for group1, group2 in comparisons:
        comp_name = f"{group1}_vs_{group2}"
        print(f"\n>>> 处理比较组: {comp_name}")
        
        # 创建该比较的输出目录
        comp_dir = os.path.join(output_dir, comp_name)
        os.makedirs(comp_dir, exist_ok=True)
        
        # 4.1 差异表达分析
        print("  执行差异表达分析 (Wilcoxon)...")
        de_results = run_dea(adata, group1, group2)
        if de_results is None or len(de_results) == 0:
            print("  未获得差异结果，跳过")
            continue
        
        # 保存完整结果
        de_file = os.path.join(comp_dir, f"DE_results_{comp_name}.csv")
        de_results.to_csv(de_file, index=False)
        print(f"  差异结果已保存: {de_file}")
        
        # 绘制火山图
        plot_volcano(de_results, comp_dir, comp_name)
        
        # 4.2 筛选显著差异基因
        sig_genes_df = de_results[
            (de_results['pvals_adj'] < 0.05) &
            (abs(de_results['logfoldchanges']) > 0.5)
        ]
        sig_gene_list = sig_genes_df['names'].tolist()
        print(f"  显著差异基因数: {len(sig_gene_list)}")
        
        if len(sig_gene_list) < 5:
            print("  显著基因不足5个，跳过富集分析")
            continue
        
        # 保存显著基因列表
        sig_file = os.path.join(comp_dir, "significant_genes.txt")
        with open(sig_file, 'w') as f:
            f.write('\n'.join(sig_gene_list))
        
        # 4.3 GO 与 KEGG 在线富集
        gene_sets = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
        for gs in gene_sets:
            if 'GO' in gs:
                sub_dir = os.path.join(comp_dir, 'GO')
            else:
                sub_dir = os.path.join(comp_dir, 'KEGG')
            os.makedirs(sub_dir, exist_ok=True)
            
            print(f"  进行 {gs} 在线富集...")
            result_file = run_online_enrichment(sig_gene_list, gs, sub_dir, organism='Human')
            if result_file and os.path.exists(result_file):
                plot_enrichment_from_local(result_file, gs, sub_dir, comp_name)
            else:
                print(f"    {gs} 富集失败或无显著结果")
    
    # 5. 生成摘要
    with open(os.path.join(output_dir, "analysis_summary.txt"), 'w') as f:
        f.write("Analysis completed.\n")
        f.write(f"Output directory: {output_dir}\n")
        f.write(f"Comparisons performed: {comparisons}\n")
    
    print("\n" + "=" * 60)
    print("分析完成！所有结果保存在: " + output_dir)
    print("=" * 60)

if __name__ == "__main__":
    main()
