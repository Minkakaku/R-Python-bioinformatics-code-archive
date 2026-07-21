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

def get_comparison_groups(adata, group_col='state'):
    groups = adata.obs[group_col].unique()
    groups = [g for g in groups if str(g) != 'nan']
    comparisons = list(combinations(groups, 2))
    _get_logger().info(f"找到的分组: {groups}, 比较组: {comparisons}")
    return groups, comparisons

def perform_dea_for_celltype(adata, celltype, group1, group2, 
                             celltype_col='celltype', group_col='state',
                             method='wilcoxon'):
    logger = _get_logger()
    celltype_data = adata[adata.obs[celltype_col] == celltype].copy()
    group1_data = celltype_data[celltype_data.obs[group_col] == group1].copy()
    group2_data = celltype_data[celltype_data.obs[group_col] == group2].copy()
    
    if group1_data.n_obs < 3 or group2_data.n_obs < 3:
        logger.warning(f"{celltype} 细胞数不足 (组1:{group1_data.n_obs}, 组2:{group2_data.n_obs})，跳过差异分析")
        return None
    
    dea_data = sc.concat([group1_data, group2_data])
    if dea_data.raw is None:
        logger.debug("concat 导致 raw 丢失，手动补充 raw")
        dea_data.raw = dea_data.copy()
    dea_data.X = dea_data.raw.X
    dea_data.obs['comparison_group'] = pd.Categorical(
        [group1]*group1_data.n_obs + [group2]*group2_data.n_obs
    )
    
    sc.pp.normalize_total(dea_data, target_sum=1e4)
    sc.pp.log1p(dea_data)
    sc.tl.rank_genes_groups(dea_data, groupby='comparison_group', method=method, use_raw=False)
    
    de_results = sc.get.rank_genes_groups_df(dea_data, group=group1)
    de_results['celltype'] = celltype
    de_results['comparison'] = f"{group1}_vs_{group2}"
    
    logger.info(f"差异分析完成: {celltype} ({group1} vs {group2})，共 {len(de_results)} 个基因")
    return de_results

def plot_volcano_plot(de_results, output_dir, comparison_name, celltype_name, 
                      logfc_threshold=0.25, pval_threshold=0.05, top_n_genes=10):
    logger = _get_logger()
    if de_results is None or len(de_results) == 0:
        return
    # 准备数据
    df = de_results.copy()
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
    
    # 标注 top N 基因
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
        logger.info(f"本地富集成功: {result_file}")
        return result_file
    except Exception as e:
        logger.error(f"本地富集失败: {e}")
        return None

# ==================== 新增：组合气泡图与条形图 ====================
def plot_combined_enrichment(up_file, down_file, gene_set_name, output_dir, 
                             comparison_name, celltype_name, up_counts, down_counts, top_n=15):
    """为上/下调富集结果绘制气泡图（红/蓝）和竖直条形图（右/左）。"""
    logger = _get_logger()
    frames = []
    # 读取上调结果
    if up_file and os.path.exists(up_file):
        up_df = pd.read_csv(up_file, sep='\t')
        up_df['Direction'] = 'UP'
        up_df['Input_size'] = up_counts
        frames.append(up_df)
    else:
        logger.warning(f"上调富集文件缺失: {up_file}")
    # 读取下调结果
    if down_file and os.path.exists(down_file):
        down_df = pd.read_csv(down_file, sep='\t')
        down_df['Direction'] = 'DOWN'
        down_df['Input_size'] = down_counts
        frames.append(down_df)
    else:
        logger.warning(f"下调富集文件缺失: {down_file}")
    
    if not frames:
        logger.error("无可用富集结果，无法绘图")
        return
    
    df = pd.concat(frames, ignore_index=True)
    # 解析基因重叠数
    def parse_overlap(ov):
        try:
            return int(ov.split('/')[0])
        except:
            return 0
    df['Hit'] = df['Overlap'].apply(parse_overlap)
    df['GeneRatio'] = df['Hit'] / df['Input_size']
    df['-log10_padj'] = -np.log10(df['Adjusted P-value'])
    df['Signed_log10_padj'] = df.apply(lambda row: row['-log10_padj'] if row['Direction']=='UP' else -row['-log10_padj'], axis=1)
    
    # 筛选显著通路
    sig = df[df['Adjusted P-value'] < 0.05].copy()
    if sig.empty:
        logger.info("无显著富集通路 (p.adj < 0.05)，跳过组合图绘制")
        return
    
    # 分别取上/下调 top_n
    up_sig = sig[sig['Direction']=='UP'].nsmallest(top_n, 'Adjusted P-value')
    down_sig = sig[sig['Direction']=='DOWN'].nsmallest(top_n, 'Adjusted P-value')
    plot_data = pd.concat([up_sig, down_sig], ignore_index=True)
    if plot_data.empty:
        logger.info("无满足条件的上下调通路")
        return
    
    # 颜色映射
    colors = {'UP': '#d62728', 'DOWN': '#1f77b4'}  # 红/蓝
    
    # ---------- 气泡图 ----------
    fig, ax = plt.subplots(figsize=(10, max(6, len(plot_data)*0.4)))
    # 按方向分组绘制，便于控制颜色
    for direction, grp in plot_data.groupby('Direction'):
        ax.scatter(grp['GeneRatio'], grp['Term'], s=grp['Hit']*20, 
                   c=colors[direction], alpha=0.7, edgecolors='black', linewidth=0.5,
                   label=direction)
    ax.set_xlabel('Gene Ratio')
    ax.set_ylabel('')
    ax.set_title(f'{celltype_name}: {gene_set_name} ({comparison_name})')
    ax.legend(title='Regulation')
    # 添加图例说明气泡大小 = 基因数
    handles, labels = ax.get_legend_handles_labels()
    # 简单气泡大小图例
    bubble_sizes = [5, 10, 20]
    dummy = [plt.scatter([], [], s=s*20, c='gray', alpha=0.6, edgecolors='black') for s in bubble_sizes]
    legends = [f'{s} genes' for s in bubble_sizes]
    # 创建一个图例用于气泡大小，放在右侧
    legend2 = ax.legend(dummy, legends, title='Gene Count', loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0)
    ax.add_artist(legend2)
    plt.tight_layout()
    bubble_path = os.path.join(output_dir, f'bubble_{gene_set_name}.png')
    plt.savefig(bubble_path, dpi=150, bbox_inches='tight')
    plt.close()
    logger.info(f"    气泡图已保存: {bubble_path}")
    
    # ---------- 竖直条形图（左下调，右上调） ----------
    fig, ax = plt.subplots(figsize=(10, max(6, len(plot_data)*0.4)))
    plot_data_sorted = plot_data.sort_values('Signed_log10_padj', ascending=True)
    bar_colors = [colors[d] for d in plot_data_sorted['Direction']]
    ax.barh(plot_data_sorted['Term'], plot_data_sorted['Signed_log10_padj'], color=bar_colors)
    ax.axvline(0, color='black', linewidth=0.8)
    ax.set_xlabel('Signed -log10(Adjusted P-value)')
    ax.set_ylabel('')
    ax.set_title(f'{celltype_name}: {gene_set_name} ({comparison_name})')
    plt.tight_layout()
    bar_path = os.path.join(output_dir, f'bar_{gene_set_name}.png')
    plt.savefig(bar_path, dpi=150, bbox_inches='tight')
    plt.close()
    logger.info(f"    条形图已保存: {bar_path}")

# ==================== 整合流程 ====================
def process_celltype(adata, celltype, group1, group2, comp_name, base_output_dir):
    ct_log_dir = os.path.join(base_output_dir, comp_name, celltype)
    os.makedirs(ct_log_dir, exist_ok=True)
    log_file = os.path.join(ct_log_dir, "enrichment_analysis.log")
    logger = setup_logger_for_celltype(log_file)
    
    logger.info(f"开始处理细胞类型: {celltype}, 比较组: {comp_name}")
    de_res = perform_dea_for_celltype(adata, celltype, group1, group2)
    if de_res is None or len(de_res) == 0:
        logger.warning(f"无差异表达结果，终止 {celltype} 分析")
        return
    
    de_res.to_csv(os.path.join(ct_log_dir, f"de_results_{celltype}.csv"), index=False)
    plot_volcano_plot(de_res, ct_log_dir, comp_name, celltype)
    
    up_mask = (de_res['pvals_adj'] < 0.05) & (de_res['logfoldchanges'] > 0.25)
    down_mask = (de_res['pvals_adj'] < 0.05) & (de_res['logfoldchanges'] < -0.25)
    up_genes = de_res.loc[up_mask, 'names'].tolist()
    down_genes = de_res.loc[down_mask, 'names'].tolist()
    up_counts = len(up_genes)
    down_counts = len(down_genes)
    logger.info(f"上调基因: {up_counts}, 下调基因: {down_counts}")
    
    gene_sets = ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023',
                 'GO_Molecular_Function_2023', 'KEGG_2021_Human']
    directions = [('UP', up_genes), ('DOWN', down_genes)]
    
    # 存储富集结果文件路径： {gene_set: {'UP': path, 'DOWN': path}}
    enrich_files = {gs: {} for gs in gene_sets}
    
    # 1) 运行所有富集分析
    for direction_label, gene_list in directions:
        if len(gene_list) < 5:
            logger.warning(f"{direction_label} 基因少于5个，跳过此方向的所有富集分析")
            continue
        
        for gs in gene_sets:
            suffix = gs.replace('_2023','').replace('_2021_Human','')
            sub_dir = os.path.join(ct_log_dir, direction_label, suffix)
            os.makedirs(sub_dir, exist_ok=True)
            
            result_file = run_local_enrichment(gene_list, gs, sub_dir)
            if result_file:
                enrich_files[gs][direction_label] = result_file
            else:
                logger.error(f"  {gs} {direction_label} 本地富集未成功")
    
    # 2) 为每个基因集绘制组合图（气泡图、条形图）
    for gs in gene_sets:
        up_file = enrich_files[gs].get('UP')
        down_file = enrich_files[gs].get('DOWN')
        if up_file or down_file:
            plot_combined_enrichment(up_file=up_file, down_file=down_file,
                                     gene_set_name=gs, output_dir=ct_log_dir,
                                     comparison_name=comp_name, celltype_name=celltype,
                                     up_counts=up_counts, down_counts=down_counts)
    
    logger.info(f"细胞类型 {celltype} 分析完成\n")

# ==================== 主程序 ====================
def main():
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    input_file = "/mnt/NexData_file/home/hy/Project/Project/FWSP20255000-lanxi/sub/celltype_subsets/combine/All_filter_state_samp_new_with_CNs10-15.h5ad"
    output_dir = "./DEA_Analysis_Publication_Quality_new_0.25"
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*60)
    print("Differential Expression + Enrichment Analysis Pipeline (Local GMT)")
    print("="*60)
    
    adata = sc.read_h5ad(input_file)
    adata = save_raw_counts(adata)
    adata = normalize_data(adata)
    _, comparisons = get_comparison_groups(adata)
    
    for group1, group2 in comparisons:
        comp_name = f"{group1}_vs_{group2}"
        print(f"\n处理比较组: {comp_name}")
        celltypes = adata.obs['celltype'].unique()
        celltypes = [ct for ct in celltypes if str(ct) != 'nan']
        for ct in celltypes:
            process_celltype(adata, ct, group1, group2, comp_name, output_dir)
    
    summary_log = os.path.join(output_dir, "analysis_summary.log")
    with open(summary_log, 'w') as f:
        f.write("Analysis completed. Detailed logs per celltype are in their respective directories.\n")
    
    print("\n" + "="*60)
    print("分析完成！所有结果保存在: " + output_dir)
    print("="*60)

if __name__ == "__main__":
    main()
