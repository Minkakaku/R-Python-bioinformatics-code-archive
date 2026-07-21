import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from scipy import stats
import scipy.sparse
import os

# 1. 读取h5ad数据并加载聚类信息
adata = sc.read_h5ad("adata_celltype.h5ad")  # 替换为你的h5ad文件路径


# --------------------------
# 核心分析函数
# --------------------------

def differential_expression_analysis(adata, group_key='group', target_genes=None, 
                                    fdr_threshold=0.05, figsize=(10,6), output_dir="results"):
    """
    完整的差异表达分析流程
    参数：
        adata: AnnData对象
        group_key: 分组列名
        target_genes: 目标基因列表
        fdr_threshold: 显著性阈值
        figsize: 图表尺寸
        output_dir: 结果输出目录
    """
    # 初始化结果目录
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # ========== 数据预处理 ==========
        adata = load_and_preprocess(adata)
        
        # ========== 数据验证 ==========
        validate_inputs(adata, group_key, target_genes)
        
        # ========== 差异分析 ==========
        results = perform_diff_analysis(adata, group_key)
        
        if results.empty:
            print("警告：未发现任何差异表达基因")
            return pd.DataFrame(), None
        
        # ========== 结果处理 ==========
        processed_results = process_results(results, fdr_threshold)
        
        # ========== 可视化 ==========
        generate_visualizations(adata, processed_results, group_key, output_dir, figsize)
        
        # ========== 保存结果 ==========
        save_path = save_all_results(adata, processed_results, output_dir, fdr_threshold)
        
        return processed_results, save_path
    
    except Exception as e:
        print(f"分析流程中断，错误信息：{str(e)}")
        return pd.DataFrame(), None

# --------------------------
# 子功能模块
# --------------------------

def load_and_preprocess(adata):
    """数据加载与预处理"""
    # 修复细胞名重复问题
    if not adata.obs_names.is_unique:
        adata.obs_names_make_unique(join="-cell_")
        print("已修复重复细胞名，使用-cell_作为分隔符")
    
    # 预处理检查
    if 'raw_counts' not in adata.layers:
        print("检测到原始数据，开始预处理...")
        adata.layers['raw_counts'] = adata.X.copy()
        
        # 基础过滤
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=10)
        
        # 处理负值
        if scipy.sparse.issparse(adata.X):
            adata.X.data[adata.X.data < 0] = 0
            adata.X.eliminate_zeros()
        else:
            adata.X[adata.X < 0] = 0
        
        # 标准化
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.layers['normalized'] = adata.X.copy()
    else:
        print("使用已预处理数据")
        adata.X = adata.layers['normalized']
    
    return adata

def validate_inputs(adata, group_key, target_genes):
    """输入验证"""
    if group_key not in adata.obs:
        raise ValueError(f"分组列 {group_key} 不存在")
    
    groups = adata.obs[group_key].unique()
    if len(groups) != 2:
        raise ValueError(f"需要两组比较，当前分组：{groups}")
    
    if target_genes:
        missing_genes = [g for g in target_genes if g not in adata.var_names]
        if missing_genes:
            print(f"警告：以下基因不存在，已忽略：{missing_genes[:5]}{'...' if len(missing_genes)>5 else ''}")
        adata = adata[:, [g for g in target_genes if g in adata.var_names]]
    
    if adata.n_vars == 0:
        raise ValueError("无有效基因可供分析")

def perform_diff_analysis(adata, group_key):
    """执行差异分析"""
    groups = adata.obs[group_key].cat.categories
    results = pd.DataFrame(index=adata.var_names,
                          columns=['logFC', 'pval', 'FDR', 'test_method', 'zero_ratio_diff'])
    
    for gene in adata.var_names:
        expr = adata[:, gene].X.flatten()
        group1 = expr[adata.obs[group_key] == groups[0]]
        group2 = expr[adata.obs[group_key] == groups[1]]
        
        # 跳过全零基因
        if np.all(expr == 0):
            continue
            
        try:
            # 计算零值比例差异
            zero_ratio_diff = np.mean(group2 == 0) - np.mean(group1 == 0)
            
            # 计算logFC
            mean1 = np.mean(group1[group1 > 0]) + 1e-9
            mean2 = np.mean(group2[group2 > 0]) + 1e-9
            logFC = np.log2(mean2) - np.log2(mean1)
            
            # 选择统计方法
            if len(group1) > 5000 or len(group2) > 5000:
                _, pval = stats.mannwhitneyu(group1, group2)
                method = 'mannwhitneyu (large)'
            else:
                if check_normality(group1) and check_normality(group2):
                    _, pval = stats.ttest_ind(group1, group2)
                    method = 't-test'
                else:
                    _, pval = stats.mannwhitneyu(group1, group2)
                    method = 'mannwhitneyu'
            
            results.loc[gene] = [logFC, pval, np.nan, method, zero_ratio_diff]
        
        except Exception as e:
            print(f"基因 {gene} 分析失败: {str(e)}")
            continue
    
    # FDR校正
    valid = results.dropna(subset=['pval'])
    if not valid.empty:
        valid['FDR'] = sm.stats.multipletests(valid['pval'], method='fdr_bh')[1]
        results.update(valid)
    
    return results

def process_results(results, fdr_threshold):
    """结果后处理"""
    processed = results.replace([np.inf, -np.inf], np.nan).dropna()
    processed['significant'] = (processed['FDR'] < fdr_threshold) & (abs(processed['logFC']) > 0.5)
    return processed

def generate_visualizations(adata, results, group_key, output_dir, figsize):
    """生成可视化图表"""
    if results.empty:
        return
    
    # 火山图
    plt.figure(figsize=figsize)
    significant = results['significant']
    plt.scatter(results['logFC'], -np.log10(results['pval']), 
                c=significant, cmap='viridis', alpha=0.6, s=30)
    plt.xlabel('log2(Fold Change)')
    plt.ylabel('-log10(p-value)')
    plt.title('Volcano Plot')
    plt.savefig(f"{output_dir}/volcano.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 表达分布图
    top_genes = results[results['significant']].nsmallest(5, 'FDR').index
    for gene in top_genes:
        plt.figure(figsize=(8,5))
        df = pd.DataFrame({
            'Expression': adata[:, gene].X.flatten(),
            'Group': adata.obs[group_key]
        })
        sns.violinplot(x='Group', y='Expression', data=df, cut=0)
        plt.title(f"{gene} Expression Distribution")
        plt.savefig(f"{output_dir}/{gene}_distribution.png", dpi=150)
        plt.close()

def save_all_results(adata, results, output_dir, fdr_threshold):
    """保存所有结果"""
    # 保存差异分析结果
    results.to_csv(f"{output_dir}/differential_expression_results.csv")
    
    # 保存元数据
    with open(f"{output_dir}/analysis_report.txt", 'w') as f:
        f.write(f"分析报告\n{'='*20}\n")
        f.write(f"总基因数: {adata.n_vars}\n")
        f.write(f"显著基因数 (FDR<{fdr_threshold}): {sum(results['significant'])}\n")
        f.write(f"最大logFC: {results['logFC'].abs().max():.2f}\n")
        f.write(f"最小FDR: {results['FDR'].min():.2e}\n")
    
    # 保存可视化数据
    adata.write_h5ad(f"{output_dir}/processed_data.h5ad")
    
    return os.path.abspath(output_dir)

# --------------------------
# 辅助函数
# --------------------------

def check_normality(data, alpha=0.05):
    """改进的正态性检验"""
    if len(data) < 3 or len(data) > 5000:
        return False
    try:
        _, p = stats.shapiro(data)
        return p > alpha
    except:
        return False

def run_diagnostics(adata, group_key):
    """运行数据诊断"""
    print("\n=== 数据诊断报告 ===")
    print(f"总细胞数: {adata.n_obs}")
    print(f"疾病组细胞数: {sum(adata.obs[group_key] == adata.obs[group_key].cat.categories[0])}")
    print(f"正常组细胞数: {sum(adata.obs[group_key] == adata.obs[group_key].cat.categories[1])}")
    print(f"检测基因数: {adata.n_vars}")
    print(f"数据矩阵类型: {'稀疏矩阵' if scipy.sparse.issparse(adata.X) else '密集矩阵'}")
    print(f"标准化状态: {'已标准化' if 'normalized' in adata.layers else '未标准化'}")
    print("===================\n")

# --------------------------
# 示例使用
# --------------------------

if __name__ == "__main__":
    # 示例数据加载（需替换为实际数据）
    # adata = sc.read_h5ad("your_data.h5ad")
    
    # 目标基因列表
    target_genes = [
        "APOL2","GBP1","JAK1","JAK2","JAK3",
        "STAT1","STAT2","STAT3","STAT4","STAT5A",
        "STAT5B","STAT6","WARS","IFNG","IFNGR1","IFNGR2"
    ]
    
    # 运行分析
    results, save_path = differential_expression_analysis(
        adata,
        group_key='group',
        target_genes=target_genes,
        output_dir="DEG_Results",
        fdr_threshold=0.05
    )
    
    # 运行诊断
    run_diagnostics(adata, 'group')
    
    if not results.empty:
        print(f"分析完成，结果保存至：{save_path}")
        print("显著基因示例：")
        print(results[results['significant']].head())
    else:
        print("分析完成，但未发现显著差异基因")
