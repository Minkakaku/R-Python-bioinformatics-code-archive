import scanpy as sc
import pandas as pd
import numpy as np

def extract_cells_by_res4(h5ad_path, target_categories, output_csv='Fibro_res4_removed.csv'):
    """
    从h5ad文件中提取指定res_4类别的细胞obs信息
    
    参数:
    h5ad_path: h5ad文件路径
    target_categories: 要保留的res_4类别列表
    output_csv: 输出CSV文件名
    """
    try:
        # 读取h5ad文件
        print(f"正在读取文件: {h5ad_path}")
        adata = sc.read_h5ad(h5ad_path)
        
        print("原始数据信息:")
        print(f"细胞总数: {adata.n_obs}")
        print(f"基因总数: {adata.n_vars}")
        
        # 检查res_4列是否存在
        if 'res_4' not in adata.obs.columns:
            available_cols = list(adata.obs.columns)
            print(f"错误: obs中未找到'res_4'列")
            print(f"可用的列: {available_cols}")
            return None
        
        # 将res_4列转换为字符串类型
        adata.obs['res_4'] = adata.obs['res_4'].astype(str)
        
        # 筛选细胞
        mask = adata.obs['res_4'].isin(target_categories)
        filtered_adata = adata[mask].copy()
        
        print(f"\n筛选结果:")
        print(f"保留的细胞数: {filtered_adata.n_obs}")
        print(f"过滤掉的细胞数: {adata.n_obs - filtered_adata.n_obs}")
        
        # 详细的类别统计
        original_counts = adata.obs['res_4'].value_counts()
        filtered_counts = filtered_adata.obs['res_4'].value_counts()
        
        print("\n各类别统计:")
        print("类别\t原始数量\t保留数量")
        for category in target_categories:
            original_count = original_counts.get(category, 0)
            filtered_count = filtered_counts.get(category, 0)
            print(f"{category}\t{original_count}\t\t{filtered_count}")
        
        # 保存obs信息
        filtered_adata.obs.to_csv(output_csv)
        print(f"\n筛选后的obs信息已保存为: {output_csv}")
        
        # 返回筛选后的数据，便于后续处理
        return filtered_adata
        
    except Exception as e:
        print(f"处理过程中发生错误: {str(e)}")
        return None

# 使用示例
if __name__ == "__main__":
    # 定义参数
    h5ad_file = "adata_cluster.h5ad"  # 替换为您的文件路径
    target_res4_categories = ['0','1','3','5','7','12','13','16','19','20',
                             '23','25','28','30','35','36','38','39','41','46','47',
                             '49','50','51','53','54','56','57','58','59','61']
    output_filename = "Fibro_res4_removed.csv"
    
    # 执行提取
    result = extract_cells_by_res4(h5ad_file, target_res4_categories, output_filename)
    
    if result is not None:
        print("\n处理完成！")
        # 可以在这里添加后续处理代码
    else:
        print("\n处理失败，请检查错误信息。")
