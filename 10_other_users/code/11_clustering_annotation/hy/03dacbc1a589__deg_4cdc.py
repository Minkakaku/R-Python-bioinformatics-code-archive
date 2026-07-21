import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
import warnings
import logging
import sys

warnings.filterwarnings('ignore')
sns.set_style('whitegrid')

# 导入 gseapy
try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    print("gseapy未安装，请运行: pip install gseapy")
    GSEAPY_AVAILABLE = False

# ==================== 本地基因集文件路径映射 ====================
GENE_SET_PATHS = {
    'GO_Biological_Process_2023': '/mnt/NexData_file/home/hy/script/gsea/GO_Biological_Process_2023.txt',
    'GO_Cellular_Component_2023': '/mnt/NexData_file/home/hy/script/gsea/GO_Cellular_Component_2023.txt',
    'GO_Molecular_Function_2023': '/mnt/NexData_file/home/hy/script/gsea/GO_Molecular_Function_2023.txt',
    'KEGG_2021_Human': '/mnt/NexData_file/home/hy/script/gsea/KEGG_2021_Human.txt'
}

# ==================== 全局日志器 ====================
LOGGER_NAME = "DEA_Enrichment"
_logger_initialized = False

def _get_logger():
    logger = logging.getLogger(LOGGER_NAME)
    if not _logger_initialized:
        logger.setLevel(logging.INFO)
        logger.addHandler(logging.NullHandler())
    return logger

def setup_logger_for_celltype(log_file_path):
    logger = _get_logger()
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    file_handler = logging.FileHandler(log_file_path, mode='w', encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger

# ==================== 差异表达分析相关函数 ====================
def save_raw_counts(adata):
    adata.raw = adata.copy()
    _get_logger().info("已将原始计数矩阵保存到adata.raw")
    return adata

def normalize_data(adata, target_sum=1e4):
    logger = _get_logger()
    logger.info("开始数据标准化...")
    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=target_sum)
    sc.pp.log1p(adata_norm)
    adata.layers['log1p_normalized'] = adata_norm.X.copy()
    logger.info("数据标准化完成")
    return adata

def plot_volcano_plot(de_results, output_dir, comparison_name, celltype_name,
                      logfc_threshold=0.25, pval_threshold=0.05, top_n_genes=10):
    """绘制火山图（已处理 pvals_adj 为 0 的情况）"""
    logger = _get_logger()
    if de_results is None or len(de_results) == 0:
        return
    df = de_results.copy()
    df['pvals_adj'] = df['pvals_adj'].clip(lower=np.nextafter(0, 1))
    df['-log10_pval'] = -np.log10(df['pvals_adj'])
    df['significant'] = 'Not significant'
    up_mask = (df['pvals_adj'] < pval_threshold) & (df['logfoldchanges'] > logfc_threshold)
    down_mask = (df['pvals_adj'] < pval_threshold) & (df['logfoldchanges'] < -logfc_threshold)
    df.loc[up_mask, 'significant'] = 'Up'
    df.loc[down_mask, 'significant'] = 'Down'

    fig, ax = plt.subplots(figsize=(10, 8))
    colors = {'Not significant': 'grey', 'Up': 'red', 'Down': 'blue'}
    for label, color in colors.items():
        subset = df[df['significant'] == label]
        ax.scatter(subset['logfoldchanges'], subset['-log10_pval'], c=color, label=label, alpha=0.6, s=10)

    top_up = df[up_mask].nlargest(top_n_genes, 'logfoldchanges')
    top_down = df[down_mask].nsmallest(top_n_genes, 'logfoldchanges')
    for _, row in pd.concat([top_up, top_down]).iterrows():
        ax.text(row['logfoldchanges'], row['-log10_pval'], row['names'], fontsize=8)

    ax.axhline(-np.log10(pval_threshold), linestyle='--', color='black', alpha=0.3)
    ax.axvline(logfc_threshold, linestyle='--', color='black', alpha=0.3)
    ax.axvline(-logfc_threshold, linestyle='--', color='black', alpha=0.3)
    ax.set_xlabel('Log2 Fold Change')
    ax.set_ylabel('-Log10 Adjusted P-value')
    ax.set_title(f'{celltype_name}: {comparison_name}')
    ax.legend()
    plt.tight_layout()
    save_path = os.path.join(output_dir, f'volcano_{comparison_name}_{celltype_name}.png')
    plt.savefig(save_path, dpi=150)
    plt.close()
    logger.info(f"    火山图已保存: {save_path}")

def plot_heatmap(de_results, dea_data, output_dir, comparison_name, celltype_name,
                 top_n=10, max_cells=1000):
    """绘制 Top N 差异基因的热图"""
    logger = _get_logger()
    if de_results is None or dea_data is None:
        return

    up_mask = (de_results['pvals_adj'] < 0.05) & (de_results['logfoldchanges'] > 0.25)
    down_mask = (de_results['pvals_adj'] < 0.05) & (de_results['logfoldchanges'] < -0.25)

    up_genes = de_results[up_mask].nlargest(top_n, 'logfoldchanges')['names'].tolist()
    down_genes = de_results[down_mask].nsmallest(top_n, 'logfoldchanges')['names'].tolist()
    selected_genes = up_genes + down_genes

    if not selected_genes:
        logger.warning("没有显著差异基因，跳过热图")
        return

    expr_df = sc.get.obs_df(dea_data, keys=selected_genes).T

    n_cells = expr_df.shape[1]
    sampled_group_labels = dea_data.obs['comparison_group'].values
    if n_cells > max_cells:
        rng = np.random.RandomState(42)
        sampled_indices = rng.choice(range(n_cells), size=max_cells, replace=False)
        expr_df = expr_df.iloc[:, sampled_indices]
        sampled_group_labels = sampled_group_labels[sampled_indices]
        logger.info(f"细胞数从 {n_cells} 下采样至 {max_cells}，以绘制热图")

    gene_stds = expr_df.std(axis=1)
    constant_genes = gene_stds[gene_stds == 0].index.tolist()
    if constant_genes:
        expr_df = expr_df.drop(index=constant_genes)
        logger.debug(f"移除 {len(constant_genes)} 个恒定表达基因: {constant_genes}")
    if expr_df.shape[0] < 2:
        logger.warning("过滤后可用基因少于2个，无法聚类，跳过热图")
        return

    group_cats = dea_data.obs['comparison_group'].cat.categories
    palette = sns.color_palette("Set1", len(group_cats))
    group_to_color = dict(zip(group_cats, palette))
    col_colors = [group_to_color[g] for g in sampled_group_labels]

    width = 12
    height = max(6, expr_df.shape[0] * 0.4)
    g = sns.clustermap(expr_df, z_score=0, cmap='RdBu_r',
                       col_cluster=True, row_cluster=True,
                       col_colors=col_colors,
                       figsize=(width, height),
                       xticklabels=False, yticklabels=True,
                       dendrogram_ratio=0.1)

    import matplotlib.patches as mpatches
    handles = [mpatches.Patch(color=color, label=cat) for cat, color in zip(group_cats, palette)]
    g.ax_heatmap.legend(handles=handles, title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    g.fig.suptitle(f'{celltype_name}: {comparison_name} Top {top_n} DE genes', y=1.02)

    save_path = os.path.join(output_dir, f'heatmap_{comparison_name}_{celltype_name}.png')
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    logger.info(f"    热图已保存: {save_path}")

# ==================== 【新增】指定基因列表热图 ====================
def plot_target_gene_heatmap(dea_data, target_genes, output_dir, comparison_name, celltype_name='Global', max_cells=1000):
    """针对用户指定的候选基因绘制特异性热图"""
    logger = _get_logger()
    if not target_genes:
        return

    expr_df = sc.get.obs_df(dea_data, keys=target_genes).T
    n_cells = expr_df.shape[1]
    sampled_group_labels = dea_data.obs['comparison_group'].values

    if n_cells > max_cells:
        rng = np.random.RandomState(42)
        sampled_indices = rng.choice(range(n_cells), size=max_cells, replace=False)
        expr_df = expr_df.iloc[:, sampled_indices]
        sampled_group_labels = sampled_group_labels[sampled_indices]
        logger.info(f"细胞数从 {n_cells} 下采样至 {max_cells}，以绘制指定基因热图")

    # 移除方差为0的基因
    gene_stds = expr_df.std(axis=1)
    constant_genes = gene_stds[gene_stds == 0].index.tolist()
    if constant_genes:
        expr_df = expr_df.drop(index=constant_genes)

    if expr_df.shape[0] < 2:
        logger.warning("过滤后可用目标基因少于2个，无法绘制热图")
        return

    group_cats = dea_data.obs['comparison_group'].cat.categories
    palette = sns.color_palette("Set1", len(group_cats))
    group_to_color = dict(zip(group_cats, palette))
    col_colors = [group_to_color[g] for g in sampled_group_labels]

    width = 12
    height = max(6, expr_df.shape[0] * 0.4)
    g = sns.clustermap(expr_df, z_score=0, cmap='RdBu_r',
                       col_cluster=True, row_cluster=True,
                       col_colors=col_colors,
                       figsize=(width, height),
                       xticklabels=False, yticklabels=True,
                       dendrogram_ratio=0.1)

    import matplotlib.patches as mpatches
    handles = [mpatches.Patch(color=color, label=cat) for cat, color in zip(group_cats, palette)]
    g.ax_heatmap.legend(handles=handles, title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    g.fig.suptitle(f'Target Genes of Interest: {comparison_name}', y=1.02)

    save_path = os.path.join(output_dir, f'heatmap_target_genes_{comparison_name}.png')
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    logger.info(f"    指定基因热图已保存: {save_path}")

# ==================== 本地 GMT 解析与富集函数 ====================
def parse_gmt_file(filepath):
    gene_sets = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            if len(parts) == 2:
                term = parts[0]
                genes = parts[1].split('\t')
            else:
                term = parts[0]
                genes = parts[2:]
            genes = [g.strip() for g in genes if g.strip()]
            if genes:
                gene_sets[term] = genes
    return gene_sets

def run_local_enrichment(gene_list, gene_set, output_dir, organism='Human'):
    logger = _get_logger()
    if not GSEAPY_AVAILABLE:
        logger.error("gseapy未安装，无法进行富集分析")
        return None

    local_path = GENE_SET_PATHS.get(gene_set)
    if local_path is None:
        logger.error(f"未找到基因集 {gene_set} 的本地文件路径")
        return None
    if not os.path.exists(local_path):
        logger.error(f"本地基因集文件不存在: {local_path}")
        return None

    try:
        gene_set_dict = parse_gmt_file(local_path)
        logger.debug(f"成功加载 {len(gene_set_dict)} 条通路")

        enr = gp.enrich(
            gene_list=gene_list,
            gene_sets=gene_set_dict,
            outdir=None,
            background=None,
            no_plot=True,
            verbose=False
        )

        result_file = os.path.join(output_dir, f"{gene_set}.Human.enrichr.reports.txt")
        enr.results.to_csv(result_file, sep='\t', index=False)
        logger.info(f"本地富集成功，结果已保存至: {result_file}")
        return result_file
    except Exception as e:
        logger.error(f"本地富集失败: {e}")
        return None

def plot_enrichment_from_local(result_file, gene_set_name, output_dir, comparison_name, celltype_name, regulation_label='', top_n=15):
    logger = _get_logger()
    if not os.path.exists(result_file):
        logger.warning(f"结果文件不存在: {result_file}")
        return
    try:
        df = pd.read_csv(result_file, sep='\t')
        if df.empty:
            logger.warning("富集结果为空，跳过绘图")
            return

        df_sig = df[df['Adjusted P-value'] < 0.05].head(top_n)
        if df_sig.empty:
            logger.info("没有显著富集的通路 (Adjusted P-value < 0.05)")
            return

        fig, ax = plt.subplots(figsize=(10, max(6, len(df_sig)*0.4)))
        df_sig = df_sig.sort_values('Adjusted P-value', ascending=False)
        ax.barh(df_sig['Term'], -np.log10(df_sig['Adjusted P-value']), color='steelblue')
        ax.set_xlabel('-Log10 Adjusted P-value')
        ax.set_title(f'{celltype_name}: {regulation_label} {gene_set_name}')
        plt.tight_layout()
        save_path = os.path.join(output_dir, f'enrichment_{gene_set_name}_{regulation_label}.png')
        plt.savefig(save_path, dpi=150)
        plt.close()
        logger.info(f"    富集图已保存: {save_path}")
    except Exception as e:
        logger.error(f"读取或绘图失败: {e}")

# ==================== 【修改】全局差异与富集分析函数 ====================
def process_global_comparison(adata, group1, group2, comp_name, base_output_dir):
    comp_dir = os.path.join(base_output_dir, comp_name)
    os.makedirs(comp_dir, exist_ok=True)
    log_file = os.path.join(comp_dir, "global_analysis.log")
    logger = setup_logger_for_celltype(log_file)

    logger.info(f"开始全局差异分析: {comp_name}")
    
    # 过滤出这两个组的所有细胞（不分细胞类型）
    adata_sub = adata[adata.obs['state'].isin([group1, group2])].copy()
    n1 = sum(adata_sub.obs['state'] == group1)
    n2 = sum(adata_sub.obs['state'] == group2)
    if n1 < 3 or n2 < 3:
        logger.warning(f"组 {group1} 或 {group2} 细胞数不足 (组1:{n1}, 组2:{n2})，跳过差异分析")
        return

    # 准备差异分析标签
    adata_sub.obs['comparison_group'] = pd.Categorical(
        [group1]*n1 + [group2]*n2
    )

    # 归一化与对数转换
    if adata_sub.raw is None:
        adata_sub.raw = adata_sub.copy()
    adata_sub.X = adata_sub.raw.X
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)
    
    # Wilcoxon 差异分析
    sc.tl.rank_genes_groups(adata_sub, groupby='comparison_group', method='wilcoxon', use_raw=False)
    de_res = sc.get.rank_genes_groups_df(adata_sub, group=group1)
    de_res['comparison'] = comp_name
    de_res['pvals_adj'] = de_res['pvals_adj'].clip(lower=np.nextafter(0, 1)) # 修复无限大问题

    # 保存结果与基础绘图
    de_res.to_csv(os.path.join(comp_dir, f"de_results_global.csv"), index=False)
    plot_volcano_plot(de_res, comp_dir, comp_name, "Global")
    plot_heatmap(de_res, adata_sub, comp_dir, comp_name, "Global", top_n=10)

    # ================== 【新增】指定基因热图绘制 ==================
    target_genes_candidates = [
        'IGFBP4', 'WFDC2', 'TFF1', 'LTBP2', 'COL18A1', 'CCN3', 'SELL', 'PON3', 
        'ADGRG2', 'IGF2R', 'PZP', 'CDON', 'TNXB', 'APLP1', 'CR2', 'GDF15'
    ]
    valid_genes = [g for g in target_genes_candidates if g in adata_sub.var_names]
    if valid_genes:
        logger.info(f"检测到 {len(valid_genes)} 个目标基因，正在绘制专用热图...")
        plot_target_gene_heatmap(adata_sub, valid_genes, comp_dir, comp_name, "Global")
    else:
        logger.warning("目标基因列表均未在当前数据集中找到。")

    # ================== 功能富集分析 ==================
    up_mask = (de_res['pvals_adj'] < 0.05) & (de_res['logfoldchanges'] > 0.25)
    down_mask = (de_res['pvals_adj'] < 0.05) & (de_res['logfoldchanges'] < -0.25)
    up_genes = de_res.loc[up_mask, 'names'].tolist()
    down_genes = de_res.loc[down_mask, 'names'].tolist()
    logger.info(f"全局上调基因: {len(up_genes)}, 下调基因: {len(down_genes)}")

    gene_sets = ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023',
                 'GO_Molecular_Function_2023', 'KEGG_2021_Human']
    directions = [('UP', up_genes), ('DOWN', down_genes)]

    for direction_label, gene_list in directions:
        if len(gene_list) < 5:
            logger.warning(f"{direction_label} 基因少于5个，跳过富集分析")
            continue
        for gs in gene_sets:
            suffix = gs.replace('_2023','').replace('_2021_Human','')
            sub_dir = os.path.join(comp_dir, direction_label, suffix)
            os.makedirs(sub_dir, exist_ok=True)

            result_file = run_local_enrichment(gene_list, gs, sub_dir)
            if result_file and os.path.exists(result_file):
                plot_enrichment_from_local(result_file, gs, sub_dir, comp_name, "Global",
                                           regulation_label=direction_label)
            else:
                logger.error(f"  {gs} {direction_label} 本地富集未生成有效文件")

    logger.info(f"全局比较 {comp_name} 分析完成\n")

# ==================== 主程序 ====================
def main():
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    input_file = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/All_filter_state_samp_new_with_CNs10-15.h5ad"
    
    # === 【修改】结果保存路径 ===
    output_dir = "./DEG_CRC_4_CDC"
    os.makedirs(output_dir, exist_ok=True)

    print("="*60)
    print("Differential Expression + Enrichment Analysis Pipeline (Local GMT)")
    print("Modified: Global comparison (bulk-like), A vs N, C vs A, C vs N")
    print("Output Path:", output_dir)
    print("="*60)

    adata = sc.read_h5ad(input_file)
    adata = save_raw_counts(adata)
    adata = normalize_data(adata)

    desired_comparisons = [('A', 'N'), ('C', 'A'), ('C', 'N')]
    available_states = set(adata.obs['state'].unique())
    comparisons = []
    for g1, g2 in desired_comparisons:
        if g1 in available_states and g2 in available_states:
            comparisons.append((g1, g2))
        else:
            print(f"警告: 状态 {g1} 或 {g2} 不在数据中，跳过该比较")

    if not comparisons:
        print("错误: 没有有效的比较组，请检查 state 列中的标签")
        return

    for group1, group2 in comparisons:
        comp_name = f"{group1}_vs_{group2}"
        print(f"\n处理比较组: {comp_name}")
        
        # === 【修改】直接调用全局差异分析函数，不再按细胞类型循环 ===
        process_global_comparison(adata, group1, group2, comp_name, output_dir)

    summary_log = os.path.join(output_dir, "analysis_summary.log")
    with open(summary_log, 'w') as f:
        f.write("Analysis completed. Global logs per comparison are in their respective directories.\n")

    print("\n" + "="*60)
    print("分析完成！所有结果保存在: " + output_dir)
    print("="*60)

if __name__ == "__main__":
    main()
