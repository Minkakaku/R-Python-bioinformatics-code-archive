import scanpy as sc
import pandas as pd

def extract_samples_and_save(h5ad_path, target_samples, output_filename="TMA.h5ad"):
    """
    从h5ad文件中提取指定样本并保存
    
    参数:
    h5ad_path: h5ad文件路径
    target_samples: 要提取的样本列表
    output_filename: 输出文件名
    """
    # 读取h5ad文件
    adata = sc.read_h5ad(h5ad_path)
    
    print(f"原始数据: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    print(f"原始样本数量: {adata.obs['sample'].nunique()}")
    
    # 检查sample列是否存在
    if 'sample' not in adata.obs.columns:
        raise ValueError("在adata.obs中未找到'sample'列。可用的列有: " + ", ".join(adata.obs.columns))
    
    # 检查哪些样本存在于数据中
    existing_samples = [sample for sample in target_samples if sample in adata.obs['sample'].values]
    missing_samples = [sample for sample in target_samples if sample not in adata.obs['sample'].values]
    
    print(f"\n=== 样本匹配情况 ===")
    print(f"目标样本数: {len(target_samples)}")
    print(f"找到的样本数: {len(existing_samples)}")
    print(f"缺失的样本数: {len(missing_samples)}")
    
    if missing_samples:
        print(f"未找到的样本: {missing_samples}")
    
    if len(existing_samples) == 0:
        print("错误: 没有找到任何目标样本！")
        print("数据集中存在的样本:", list(adata.obs['sample'].unique())[:20])
        return None
    
    # 提取目标样本
    adata_subset = adata[adata.obs['sample'].isin(existing_samples), :].copy()
    
    print(f"\n=== 提取结果 ===")
    print(f"提取后的细胞数: {adata_subset.n_obs}")
    print(f"提取后的基因数: {adata_subset.n_vars}")
    print("\n样本分布统计:")
    sample_stats = adata_subset.obs['sample'].value_counts()
    print(sample_stats)
    
    # 保存为新的h5ad文件
    adata_subset.write(output_filename)
    print(f"\n数据已保存为: {output_filename}")
    
    return adata_subset

# 要提取的样本列表
target_samples = ["F1803468-2","Z1635989-12","Z1802235-10","Z1814675-2","Z1826868-6","Z1952876-9","Z2057500-35","Z2164107-11","Z2171766-3","Z2323818-8","Z2336500-6","Z2345459-7","Z2353385-10","ZF2304103-8"]
# 执行提取
if __name__ == "__main__":
    h5ad_path = "adata_celltype_RCTD_removed.h5ad"  # 替换为您的文件路径
    extract_samples_and_save(h5ad_path, target_samples, "TMA_5.h5ad")
