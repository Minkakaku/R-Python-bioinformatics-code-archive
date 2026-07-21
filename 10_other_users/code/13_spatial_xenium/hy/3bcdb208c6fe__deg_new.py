import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# 导入 gseapy
try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    print("gseapy未安装，请运行: pip install gseapy")
    GSEAPY_AVAILABLE = False

# ==================== 差异表达分析相关函数 ====================
def save_raw_counts(adata):
    adata.raw = adata.copy()
    print(f"已将原始计数矩阵保存到adata.raw")
    return adata

def normalize_data(adata, target_sum=1e4):
    print("开始数据标准化...")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=target_sum)
    sc.pp.log1p(adata_norm)
    adata.layers['log1p_normalized'] = adata_norm.X.copy()
    print(f"数据标准化完成")
    return adata

def get_comparison_groups(adata, group_col='group_final'):
    groups = adata.obs[group_col].unique()
    groups = [g for g in groups if str(g) != 'nan']
    comparisons = list(combinations(groups, 2))
    print(f"找到的分组: {groups}, 比较组: {comparisons}")
    return groups, comparisons

def perform_dea_for_celltype(adata, celltype, group1, group2, 
                             celltype_col='celltype', group_col='group_final',
                             method='wilcoxon'):
    celltype_data = adata[adata.obs[celltype_col] == celltype].copy()
    group1_data = celltype_data[celltype_data.obs[group_col] == group1].copy()
    group2_data = celltype_data[celltype_data.obs[group_col] == group2].copy()
    if group1_data.n_obs < 3 or group2_data.n_obs < 3:
        print(f"  警告: {celltype} 细胞数不足，跳过")
        return None
    dea_data = sc.concat([group1_data, group2_data])
    dea_data.obs['comparison_group'] = pd.Categorical([group1]*group1_data.n_obs + [group2]*group2_data.n_obs)
    dea_data.X = dea_data.raw.X
    sc.pp.normalize_total(dea_data, target_sum=1e4)
    sc.pp.log1p(dea_data)
    sc.tl.rank_genes_groups(dea_data, groupby='comparison_group', method=method, use_raw=False)
    de_results = sc.get.rank_genes_groups_df(dea_data, group=group1)
    de_results['celltype'] = celltype
    de_results['comparison'] = f"{group1}_vs_{group2}"
    return de_results

def plot_volcano_plot(de_results, output_dir, comparison_name, celltype_name, 
                      logfc_threshold=0.25, pval_threshold=0.05, top_n_genes=10):
    if de_results is None or len(de_results)==0:
        return
    fig, ax = plt.subplots(figsize=(10,8))
    plot_data = de_results.copy()
    plot_data['neg_log10_pval'] = -np.log10(plot_data['pvals_adj'])
    max_val = plot_data['neg_log10_pval'].replace([np.inf,-np.inf], np.nan).max()
    if pd.isna(max_val): max_val = 10
    plot_data['neg_log10_pval'] = plot_data['neg_log10_pval'].replace([np.inf,-np.inf], max_val)
    plot_data['regulation'] = 'Not significant'
    up = (plot_data['logfoldchanges'] > logfc_threshold) & (plot_data['pvals_adj'] < pval_threshold)
    down = (plot_data['logfoldchanges'] < -logfc_threshold) & (plot_data['pvals_adj'] < pval_threshold)
    plot_data.loc[up, 'regulation'] = f'Up (FC > {logfc_threshold})'
    plot_data.loc[down, 'regulation'] = f'Down (FC < -{logfc_threshold})'
    colors = {'Not significant':'#d3d3d3', f'Up (FC > {logfc_threshold})':'#e74c3c', f'Down (FC < -{logfc_threshold})':'#3498db'}
    for reg, col in colors.items():
        mask = plot_data['regulation'] == reg
        if mask.sum():
            ax.scatter(plot_data.loc[mask,'logfoldchanges'], plot_data.loc[mask,'neg_log10_pval'],
                       c=col, s=8, alpha=0.6, label=reg, edgecolor='none')
    ax.axvline(x=logfc_threshold, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=-logfc_threshold, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=-np.log10(pval_threshold), color='gray', linestyle='--', alpha=0.5)
    sig = plot_data[plot_data['regulation']!='Not significant'].copy()
    if len(sig)>0:
        sig = sig.sort_values('neg_log10_pval', ascending=False)
        for idx, row in sig.head(top_n_genes*2).iterrows():
            ax.text(row['logfoldchanges'], row['neg_log10_pval'], row['names'], fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.7))
    ax.set_xlabel('log2(Fold Change)')
    ax.set_ylabel('-log10(Adjusted P-value)')
    ax.set_title(f"{comparison_name} - {celltype_name} Volcano Plot")
    ax.legend()
    ax.grid(alpha=0.3, linestyle='--')
    n_up = (plot_data['regulation'] == f'Up (FC > {logfc_threshold})').sum()
    n_down = (plot_data['regulation'] == f'Down (FC < -{logfc_threshold})').sum()
    ax.text(0.02, 0.98, f"Significant genes: {n_up+n_down}\nUp: {n_up}, Down: {n_down}",
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle="round", facecolor='white', alpha=0.8))
    plt.tight_layout()
    for ext in ['.png','.pdf']:
        plt.savefig(os.path.join(output_dir, f"volcano_plot{ext}"), dpi=300)
    plt.close()
    print(f"    已保存火山图")

# ==================== 在线富集函数（保存结果到本地） ====================
def run_online_enrichment(gene_list, gene_set, output_dir, organism='Human'):
    """
    调用 Enrichr 在线富集，将结果保存为标准格式文件。
    返回：保存的文件路径，如果失败则返回 None。
    """
    if not GSEAPY_AVAILABLE:
        print("    gseapy未安装，无法进行在线富集")
        return None
    try:
        # 运行富集
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=[gene_set],
            organism=organism,
            outdir=output_dir,   # gseapy 会在此目录下生成 {gene_set}.Human.enrichr.reports.txt
            cutoff=0.05
        )
        # 检查结果文件是否生成
        expected_file = os.path.join(output_dir, f"{gene_set}.Human.enrichr.reports.txt")
        if os.path.exists(expected_file):
            print(f"      在线富集成功，结果已保存至: {expected_file}")
            return expected_file
        else:
            print(f"      在线富集完成但未找到预期文件: {expected_file}")
            return None
    except Exception as e:
        print(f"      在线富集失败: {e}")
        return None

# ==================== 从本地文件读取并绘图 ====================
def plot_enrichment_from_local(result_file, gene_set_name, output_dir, comparison_name, celltype_name, regulation_label='', top_n=15):
    """
    读取本地富集结果文件（.enrichr.reports.txt 格式），绘制条形图、气泡图、点图。
    regulation_label: 用于区分上调/下调的文字（如 'UP' 或 'DOWN'），会插入文件名和标题。
    """
    if not os.path.exists(result_file):
        print(f"      结果文件不存在: {result_file}")
        return
    try:
        df = pd.read_csv(result_file, sep='\t')
        if 'Term' not in df.columns or 'Adjusted P-value' not in df.columns:
            if len(df.columns) >= 5:
                df = df.rename(columns={df.columns[0]: 'Term', df.columns[4]: 'Adjusted P-value'})
            else:
                print(f"      文件列名不匹配: {df.columns.tolist()}")
                return
        df['Adjusted P-value'] = pd.to_numeric(df['Adjusted P-value'], errors='coerce')
        df = df.dropna(subset=['Adjusted P-value'])
        if len(df) == 0:
            print(f"      没有有效的P值")
            return
        top_terms = df.sort_values('Adjusted P-value').head(top_n).copy()
        if len(top_terms) < 2:
            print(f"      显著条目不足2个")
            return
        top_terms['neg_log10'] = -np.log10(top_terms['Adjusted P-value'])

        # 构建标题后缀和文件名前缀
        title_suffix = f" ({regulation_label})" if regulation_label else ""
        file_prefix = f"enrichment_{gene_set_name}_{regulation_label}_" if regulation_label else f"enrichment_{gene_set_name}_"
        
        # 条形图
        fig, ax = plt.subplots(figsize=(10,8))
        top_terms_sorted = top_terms.sort_values('Adjusted P-value', ascending=True)
        colors = plt.cm.viridis(np.linspace(0.3,0.9,len(top_terms_sorted)))
        ax.barh(range(len(top_terms_sorted)), top_terms_sorted['neg_log10'], color=colors)
        ax.set_yticks(range(len(top_terms_sorted)))
        ax.set_yticklabels([t[:60]+'...' if len(t)>60 else t for t in top_terms_sorted['Term']])
        ax.set_xlabel('-log10(Adjusted P-value)')
        ax.set_title(f"{comparison_name} - {celltype_name}\n{gene_set_name} Enrichment{title_suffix}")
        ax.grid(axis='x', alpha=0.3)
        if 'Overlap' in top_terms_sorted.columns:
            for i, row in enumerate(top_terms_sorted.itertuples()):
                if isinstance(row.Overlap, str) and '/' in row.Overlap:
                    n = row.Overlap.split('/')[0]
                    ax.text(row.neg_log10+0.1, i, f'n={n}', va='center', fontsize=9)
        plt.tight_layout()
        bar_file = os.path.join(output_dir, f"{file_prefix}barplot.png")
        plt.savefig(bar_file, dpi=300)
        plt.savefig(bar_file.replace('.png','.pdf'), dpi=300)
        plt.close()
        print(f"      已保存条形图: {bar_file}")
        
        # 气泡图
        if len(top_terms) >= 3:
            fig, ax = plt.subplots(figsize=(10,8))
            if 'Overlap' in top_terms.columns:
                gene_counts = [int(ov.split('/')[0]) if isinstance(ov,str) and '/' in ov else 1 for ov in top_terms['Overlap']]
            else:
                gene_counts = [1]*len(top_terms)
            sizes = np.array(gene_counts)*20
            scat = ax.scatter(range(len(top_terms)), top_terms['neg_log10'], s=sizes, c=top_terms['neg_log10'], cmap='viridis', alpha=0.7)
            ax.set_xticks(range(len(top_terms)))
            ax.set_xticklabels([t[:40]+'...' if len(t)>40 else t for t in top_terms['Term']], rotation=45, ha='right')
            ax.set_ylabel('-log10(Adjusted P-value)')
            ax.set_title(f"{gene_set_name} Bubble Plot{title_suffix}")
            plt.colorbar(scat, label='-log10(P)')
            ax.grid(alpha=0.3)
            plt.tight_layout()
            bubble_file = os.path.join(output_dir, f"{file_prefix}bubble.png")
            plt.savefig(bubble_file, dpi=300)
            plt.savefig(bubble_file.replace('.png','.pdf'), dpi=300)
            plt.close()
            print(f"      已保存气泡图: {bubble_file}")
            
            # 点图
            fig, ax = plt.subplots(figsize=(10,8))
            gene_ratios = []
            for ov in top_terms['Overlap']:
                if isinstance(ov,str) and '/' in ov:
                    a,b = ov.split('/')
                    try:
                        gene_ratios.append(int(a)/int(b))
                    except:
                        gene_ratios.append(0.01)
                else:
                    gene_ratios.append(0.01)
            scat = ax.scatter(gene_ratios, top_terms['neg_log10'], s=100, c=top_terms['neg_log10'], cmap='plasma', alpha=0.7)
            for i, (x,y) in enumerate(zip(gene_ratios, top_terms['neg_log10'])):
                ax.text(x, y, f'  {i+1}', fontsize=9, va='center')
            ax.set_xlabel('Gene Ratio')
            ax.set_ylabel('-log10(Adjusted P-value)')
            ax.set_title(f"{gene_set_name} Dot Plot{title_suffix}")
            plt.colorbar(scat, label='-log10(P)')
            ax.grid(alpha=0.3)
            plt.tight_layout()
            dot_file = os.path.join(output_dir, f"{file_prefix}dotplot.png")
            plt.savefig(dot_file, dpi=300)
            plt.savefig(dot_file.replace('.png','.pdf'), dpi=300)
            plt.close()
            print(f"      已保存点图: {dot_file}")
    except Exception as e:
        print(f"      读取或绘图失败: {e}")

# ==================== 整合：差异表达 + 在线富集（分上调/下调） + 本地绘图 ====================
def process_celltype(adata, celltype, group1, group2, comp_name, base_output_dir):
    """对单个细胞类型执行完整流程：上调基因和下调基因分别富集"""
    print(f"  细胞类型: {celltype}")
    # 1. 差异表达分析
    de_res = perform_dea_for_celltype(adata, celltype, group1, group2)
    if de_res is None or len(de_res)==0:
        return
    
    # 创建输出目录
    ct_dir = os.path.join(base_output_dir, comp_name, celltype)
    os.makedirs(ct_dir, exist_ok=True)
    # 保存差异结果
    de_res.to_csv(os.path.join(ct_dir, f"de_results_{celltype}.csv"), index=False)
    # 火山图
    plot_volcano_plot(de_res, ct_dir, comp_name, celltype)
    
    # 2. 提取上调基因和下调基因
    up_mask = (de_res['pvals_adj'] < 0.05) & (de_res['logfoldchanges'] > 0.25)
    down_mask = (de_res['pvals_adj'] < 0.05) & (de_res['logfoldchanges'] < -0.25)
    
    up_genes = de_res.loc[up_mask, 'names'].tolist()
    down_genes = de_res.loc[down_mask, 'names'].tolist()
    print(f"    上调基因数: {len(up_genes)}, 下调基因数: {len(down_genes)}")
    
    # 3. 分别对上调和下调基因进行富集分析
    gene_sets = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
    directions = [('UP', up_genes), ('DOWN', down_genes)]
    
    for direction_label, gene_list in directions:
        if len(gene_list) < 5:
            print(f"    {direction_label} 基因少于5个，跳过富集")
            continue
        
        for gs in gene_sets:
            # 为每个 direction 和基因集建立子目录
            if 'GO' in gs:
                sub_dir = os.path.join(ct_dir, direction_label, 'GO')
            else:
                sub_dir = os.path.join(ct_dir, direction_label, 'KEGG')
            os.makedirs(sub_dir, exist_ok=True)
            
            # 调用在线富集
            result_file = run_online_enrichment(gene_list, gs, sub_dir, organism='Human')
            if result_file and os.path.exists(result_file):
                # 绘图，并将 direction_label 传递给绘图函数以区分标题和文件名
                plot_enrichment_from_local(result_file, gs, sub_dir, comp_name, celltype, 
                                           regulation_label=direction_label)
            else:
                print(f"    在线富集失败，无法生成 {gs} {direction_label} 图形")

def main():
    # 输入输出路径（请根据实际情况修改）
    input_file = "/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/All_filter.h5ad"
    output_dir = "./DEA_Analysis_Publication_Quality_new_0.25"
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*60)
    print("Differential Expression + Enrichment Analysis Pipeline (Separate UP/DOWN)")
    print("="*60)
    
    # 读取数据
    print("\n1. 读取数据...")
    adata = sc.read_h5ad(input_file)
    print(f"   {adata.n_obs} cells, {adata.n_vars} genes")
    for col in ['group_final','celltype']:
        if col not in adata.obs.columns:
            raise ValueError(f"缺少列: {col}")
    
    # 保存原始计数
    print("\n2. 保存原始计数...")
    adata = save_raw_counts(adata)
    
    # 标准化
    print("\n3. 标准化...")
    adata = normalize_data(adata)
    
    # 获取比较组
    print("\n4. 获取比较组...")
    _, comparisons = get_comparison_groups(adata)
    
    # 对每个比较组和细胞类型执行分析
    print("\n5. 执行差异表达分析和富集分析（分上调/下调）...")
    for group1, group2 in comparisons:
        comp_name = f"{group1}_vs_{group2}"
        print(f"\n处理比较组: {comp_name}")
        # 获取所有细胞类型
        celltypes = adata.obs['celltype'].unique()
        celltypes = [ct for ct in celltypes if str(ct) != 'nan']
        for ct in celltypes:
            process_celltype(adata, ct, group1, group2, comp_name, output_dir)
    
    # 生成简单摘要
    with open(os.path.join(output_dir, "analysis_summary.txt"), 'w') as f:
        f.write("Analysis completed (UP and DOWN genes processed separately).\n")
        f.write(f"Output directory: {output_dir}\n")
    
    print("\n" + "="*60)
    print("分析完成！所有结果保存在: " + output_dir)
    print("="*60)

if __name__ == "__main__":
    main()
