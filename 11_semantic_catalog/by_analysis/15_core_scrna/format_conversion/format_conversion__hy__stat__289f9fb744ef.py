import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse

def calculate_sample_metrics(adata, sample_col='sample'):
    """
    统计每个样本的基因中位数和转录本中位数
    
    参数:
    adata: AnnData对象
    sample_col: obs中存储样本信息的列名
    
    返回:
    sample_stats: 包含统计结果的DataFrame
    """
    
    # 检查样本列是否存在
    if sample_col not in adata.obs.columns:
        raise ValueError(f"样本列 '{sample_col}' 不存在于 adata.obs 中")
    
    # 获取所有样本
    samples = adata.obs[sample_col].unique()
    print(f"找到 {len(samples)} 个样本: {list(samples)}")
    
    # 存储结果的列表
    results = []
    
    # 对每个样本进行计算
    for sample in samples:
        # 获取当前样本的细胞索引
        sample_mask = adata.obs[sample_col] == sample
        sample_cells = np.sum(sample_mask)
        
        if sample_cells == 0:
            print(f"警告: 样本 {sample} 没有细胞")
            continue
        
        # 提取当前样本的数据
        sample_data = adata[sample_mask, :]
        
        # 计算每个细胞的基因数（非零值数量）
        if sparse.issparse(sample_data.X):
            # 稀疏矩阵的处理方式
            genes_per_cell = np.array(sample_data.X.getnnz(axis=1)).flatten()
        else:
            # 密集矩阵的处理方式
            genes_per_cell = np.count_nonzero(sample_data.X, axis=1)
        
        # 计算每个细胞的转录本数（表达量总和）
        if sparse.issparse(sample_data.X):
            transcripts_per_cell = np.array(sample_data.X.sum(axis=1)).flatten()
        else:
            transcripts_per_cell = sample_data.X.sum(axis=1)
        
        # 计算中位数
        median_genes = np.median(genes_per_cell)
        median_transcripts = np.median(transcripts_per_cell)
        
        # 计算其他统计量用于参考
        mean_genes = np.mean(genes_per_cell)
        mean_transcripts = np.mean(transcripts_per_cell)
        
        results.append({
            'sample': sample,
            'n_cells': sample_cells,
            'median_genes': round(median_genes, 2),
            'median_transcripts': round(median_transcripts, 2),
            'mean_genes': round(mean_genes, 2),
            'mean_transcripts': round(mean_transcripts, 2),
            'min_genes': np.min(genes_per_cell),
            'max_genes': np.max(genes_per_cell),
            'min_transcripts': np.min(transcripts_per_cell),
            'max_transcripts': np.max(transcripts_per_cell)
        })
        
        print(f"样本 {sample}: {sample_cells} 个细胞, "
              f"基因中位数: {median_genes:.0f}, "
              f"转录本中位数: {median_transcripts:.0f}")
    
    # 创建结果DataFrame
    sample_stats = pd.DataFrame(results)
    
    return sample_stats

# 主程序
if __name__ == "__main__":
    # 读取h5ad文件（请将路径替换为您的实际文件路径）
    h5ad_file_path = "adata_raw.h5ad"  # 替换为您的文件路径
    adata = sc.read_h5ad(h5ad_file_path)
    
    print("数据基本信息:")
    print(f"细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
    print(f"样本列信息: {adata.obs.columns.tolist()}")
    
    # 检查数据矩阵类型
    print(f"数据矩阵类型: {type(adata.X)}")
    
    # 计算样本统计量
    sample_stats = calculate_sample_metrics(adata, sample_col='sample')
    
    # 显示详细结果
    print("\n" + "="*80)
    print("各样本统计摘要:")
    print("="*80)
    print(sample_stats.to_string(index=False))
    
    # 保存结果到CSV文件
    output_file = "sample_quality_metrics.csv"
    sample_stats.to_csv(output_file, index=False)
    print(f"\n结果已保存到: {output_file}")
    
    # 可选：生成简要统计摘要
    print("\n" + "="*80)
    print("整体统计摘要:")
    print("="*80)
    print(f"总样本数: {len(sample_stats)}")
    print(f"总细胞数: {sample_stats['n_cells'].sum()}")
    print(f"所有样本的基因中位数范围: {sample_stats['median_genes'].min():.0f} - {sample_stats['median_genes'].max():.0f}")
    print(f"所有样本的转录本中位数范围: {sample_stats['median_transcripts'].min():.0f} - {sample_stats['median_transcripts'].max():.0f}")
    
    # 按基因中位数排序显示前5个样本
    print("\n基因中位数最高的5个样本:")
    top_genes = sample_stats.nlargest(5, 'median_genes')[['sample', 'median_genes', 'n_cells']]
    print(top_genes.to_string(index=False))
