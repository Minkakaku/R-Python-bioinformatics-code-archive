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
    'KEGG_2021_Human': '/mnt/NexData_file/home/hy/script/gsea/KEGG_2021_Human.txt'
}

# ==================== 全局日志器 ====================
LOGGER_NAME = "DEA_Enrichment"
_logger_initialized = False

def _get_logger():
    """获取或初始化全局 logger"""
    logger = logging.getLogger(LOGGER_NAME)
    if not _logger_initialized:
        logger.setLevel(logging.INFO)
        # 先添加一个空的 NullHandler 占位
        logger.addHandler(logging.NullHandler())
    return logger

def setup_logger_for_celltype(log_file_path):
    """为当前细胞类型配置日志文件和控制台输出"""
    logger = _get_logger()
    # 移除所有已有的 handler（上次细胞类型留下的）
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # 创建文件 handler
    file_handler = logging.FileHandler(log_file_path, mode='w', encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)  # 文件记录更详细级别
    
    # 创建控制台 handler
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.INFO)
    
    # 统一格式
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
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

def get_comparison_groups(adata, group_col='group_final'):
    groups = adata.obs[group_col].unique()
    groups = [g for g in groups if str(g) != 'nan']
    comparisons = list(combinations(groups, 2))
    _get_logger().info(f"找到的分组: {groups}, 比较组: {comparisons}")
    return groups, comparisons

def perform_dea_for_celltype(adata, celltype, group1, group2, 
                             celltype_col='celltype', group_col='group_final',
                             method='wilcoxon'):
    logger = _get_logger()
    celltype_data = adata[adata.obs[celltype_col] == celltype].copy()
    group1_data = celltype_data[celltype_data.obs[group_col] == group1].copy()
    group2_data = celltype_data[celltype_data.obs[group_col] == group2].copy()
    
    if group1_data.n_obs < 3 or group2_data.n_obs < 3:
        logger.warning(f"{celltype} 细胞数不足 (组1:{group1_data.n_obs}, 组2:{group2_data.n_obs})，跳过差异分析")
        return None
    
    dea_data = sc.concat([group1_data, group2_data])
    # 修复 raw 丢失的问题
    if dea_data.raw is None:
        # 尝试从原始 adata 还原 raw
        logger.debug("concat 导致 raw 丢失，手动补充 raw")
        # 直接使用原始计数（未归一化）
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
    if de_results is None or len(de_results)==0:
        return
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
    """
    手动解析 GMT 文件，返回 {通路ID: [基因列表]}
    支持标准三列（Term\\tDescription\\tGene1\\tGene2...）及两列（无Description）格式。
    """
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
                # 只有 Term 和 Genes（无描述列）
                term = parts[0]
                genes = parts[1].split('\t')
            else:
                term = parts[0]
                genes = parts[2:]  # 跳过 description 列
            genes = [g.strip() for g in genes if g.strip()]
            if genes:
                gene_sets[term] = genes
    return gene_sets

def run_local_enrichment(gene_list, gene_set, output_dir, organism='Human'):
    """使用本地 GMT 文件进行富集分析（手动解析 GMT 传入字典）"""
    logger = _get_logger()
    if not GSEAPY_AVAILABLE:
        logger.error("gseapy未安装，无法进行富集分析")
        return None
    
    local_path = GENE_SET_PATHS.get(gene_set)
    if local_path is None:
        logger.error(f"未找到基因集 {gene_set} 的本地文件路径，请在 GENE_SET_PATHS 中添加")
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
    """从本地结果文件绘制富集分析条形图"""
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

# ==================== 整合：差异表达 + 本地富集 + 绘图 ====================
def process_celltype(adata, celltype, group1, group2, comp_name, base_output_dir):
    """对单个细胞类型执行完整流程"""
    # 为当前细胞类型配置独立日志文件
    ct_log_dir = os.path.join(base_output_dir, comp_name, celltype)
    os.makedirs(ct_log_dir, exist_ok=True)
    log_file = os.path.join(ct_log_dir, "enrichment_analysis.log")
    logger = setup_logger_for_celltype(log_file)
    
    logger.info(f"开始处理细胞类型: {celltype}, 比较组: {comp_name}")
    # 1. 差异表达分析
    de_res = perform_dea_for_celltype(adata, celltype, group1, group2)
    if de_res is None or len(de_res)==0:
        logger.warning(f"无差异表达结果，终止 {celltype} 分析")
        return
    
    # 保存差异结果
    de_res.to_csv(os.path.join(ct_log_dir, f"de_results_{celltype}.csv"), index=False)
    # 火山图
    plot_volcano_plot(de_res, ct_log_dir, comp_name, celltype)
    
    # 2. 提取上调/下调基因
    up_mask = (de_res['pvals_adj'] < 0.05) & (de_res['logfoldchanges'] > 0.25)
    down_mask = (de_res['pvals_adj'] < 0.05) & (de_res['logfoldchanges'] < -0.25)
    up_genes = de_res.loc[up_mask, 'names'].tolist()
    down_genes = de_res.loc[down_mask, 'names'].tolist()
    logger.info(f"上调基因: {len(up_genes)}, 下调基因: {len(down_genes)}")
    
    # 3. 富集分析 —— 使用本地文件
    gene_sets = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
    directions = [('UP', up_genes), ('DOWN', down_genes)]
    
    for direction_label, gene_list in directions:
        if len(gene_list) < 5:
            logger.warning(f"{direction_label} 基因少于5个，跳过富集分析")
            continue
        
        for gs in gene_sets:
            suffix = 'GO' if 'GO' in gs else 'KEGG'
            sub_dir = os.path.join(ct_log_dir, direction_label, suffix)
            os.makedirs(sub_dir, exist_ok=True)
            
            result_file = run_local_enrichment(gene_list, gs, sub_dir, organism='Human')
            if result_file and os.path.exists(result_file):
                plot_enrichment_from_local(result_file, gs, sub_dir, comp_name, celltype, 
                                           regulation_label=direction_label)
            else:
                logger.error(f"  {gs} {direction_label} 本地富集未生成有效文件")
    
    logger.info(f"细胞类型 {celltype} 分析完成\n")

# ==================== 主程序 ====================
def main():
    # 全局基本日志设置（在配置细胞类型前只输出到控制台）
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    input_file = "/mnt/NexData_file/home/hy/Project/Project/quzhouCRC/xenium/sub/combine/All_filter.h5ad"
    output_dir = "./DEA_Analysis_Publication_Quality_new_0.25"
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*60)
    print("Differential Expression + Local GMT Enrichment Analysis Pipeline")
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
    
    # 生成简单摘要（同样写入总输出目录的 log）
    summary_log = os.path.join(output_dir, "analysis_summary.log")
    with open(summary_log, 'w') as f:
        f.write("Analysis completed. Detailed logs per celltype are in their respective directories.\n")
    
    print("\n" + "="*60)
    print("分析完成！所有结果保存在: " + output_dir)
    print("="*60)

if __name__ == "__main__":
    main()
