import scanpy as sc
import pandas as pd
import numpy as np

def count_unknown_celltypes(h5ad_file):
    """
    统计annotation列为Unknown的细胞在celltype列上的分布
    
    参数:
    h5ad_file: h5ad文件路径
    """
    # 读取h5ad文件
    print("读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_file)
    print(f"数据形状: {adata.shape}")
    print(f"obs列名: {list(adata.obs.columns)}")
    
    # 检查必要的列是否存在
    if 'annotation' not in adata.obs.columns:
        raise ValueError("h5ad文件的obs中不存在'annotation'列")
    
    if 'celltype' not in adata.obs.columns:
        raise ValueError("h5ad文件的obs中不存在'celltype'列")
    
    # 统计Unknown细胞的数量
    unknown_mask = adata.obs['annotation'] == 'Unknown'
    unknown_count = unknown_mask.sum()
    total_cells = adata.n_obs
    
    print(f"\n总体统计:")
    print(f"总细胞数: {total_cells}")
    print(f"Unknown注释的细胞数: {unknown_count}")
    print(f"Unknown细胞占比: {unknown_count/total_cells*100:.2f}%")
    
    if unknown_count == 0:
        print("没有找到annotation为'Unknown'的细胞")
        return
    
    # 获取Unknown细胞的celltype分布
    unknown_celltypes = adata.obs.loc[unknown_mask, 'celltype']
    
    print(f"\nUnknown细胞的celltype分布:")
    print("=" * 50)
    
    # 统计各celltype的数量和比例
    celltype_counts = unknown_celltypes.value_counts(dropna=False)
    celltype_percentages = (celltype_counts / unknown_count * 100).round(2)
    
    # 创建统计结果DataFrame
    stats_df = pd.DataFrame({
        '细胞数量': celltype_counts,
        '在Unknown中的占比(%)': celltype_percentages,
        '在总细胞中的占比(%)': (celltype_counts / total_cells * 100).round(2)
    })
    
    # 打印详细统计
    print(stats_df)
    
    # 打印更直观的统计信息
    print(f"\n详细统计:")
    print("-" * 30)
    for i, (celltype, count) in enumerate(celltype_counts.items(), 1):
        pct_in_unknown = celltype_percentages[celltype]
        pct_total = (count / total_cells * 100).round(2)
        
        if pd.isna(celltype):
            celltype_str = "NaN/空值"
        else:
            celltype_str = str(celltype)
            
        print(f"{i}. {celltype_str}: {count}个细胞 "
              f"({pct_in_unknown}%的Unknown细胞, {pct_total}%的总细胞)")
    
    return stats_df

if __name__ == "__main__":
    # 文件路径
    h5ad_file = "adata_celltype_RCTD.h5ad"  # 替换为您的文件路径
    
    try:
        result = count_unknown_celltypes(h5ad_file)
        print("\n统计完成!")
        
        # 可选：保存统计结果到CSV文件
        # if result is not None:
        #     result.to_csv('unknown_celltype_stats.csv')
        #     print("统计结果已保存至: unknown_celltype_stats.csv")
            
    except Exception as e:
        print(f"处理过程中出现错误: {str(e)}")
