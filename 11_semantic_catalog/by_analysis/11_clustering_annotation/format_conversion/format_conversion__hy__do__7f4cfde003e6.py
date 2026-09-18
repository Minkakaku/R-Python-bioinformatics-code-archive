import scanpy as sc
import pandas as pd
import numpy as np
import re

def add_group_new_column(h5ad_file, output_file=None):
    """
    为h5ad文件添加group_new列（修复版）
    """
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_file)
    print(f"成功读取文件: {h5ad_file}")
    
    # 检查必要的列是否存在
    required_cols = ['group', 'sample']
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        print(f"错误: 缺少必要的列: {missing_cols}")
        return None
    
    print(f"\n=== 数据质量检查 ===")
    print(f"group列数据类型: {adata.obs['group'].dtype}")
    print(f"sample列数据类型: {adata.obs['sample'].dtype}")
    
    # 处理数据类型问题
    adata.obs['group'] = adata.obs['group'].astype(str)
    adata.obs['sample'] = adata.obs['sample'].astype(str)
    
    # 检查空值
    null_group = adata.obs['group'].isnull().sum()
    null_sample = adata.obs['sample'].isnull().sum()
    
    if null_group > 0 or null_sample > 0:
        print(f"发现空值 - group: {null_group}, sample: {null_sample}")
        adata = adata[~adata.obs['group'].isnull() & ~adata.obs['sample'].isnull()].copy()
    
    print(f"\n=== 创建group_new列 ===")
    
    # 使用向量化操作替代循环，更高效且避免类型错误
    def create_group_new(row):
        sample_val = str(row['sample']).lower()
        group_val = str(row['group'])
        
        if 'her2' in sample_val:
            return f"{group_val}_Her2_0"
        elif 'low' in sample_val:
            return f"{group_val}_Her2_low"
        else:
            return group_val
    
    # 应用函数
    adata.obs['group_new'] = adata.obs.apply(create_group_new, axis=1)
    
    # 统计结果
    her2_count = (adata.obs['sample'].str.lower().str.contains('her2', na=False)).sum()
    low_count = (adata.obs['sample'].str.lower().str.contains('low', na=False) & 
                ~adata.obs['sample'].str.lower().str.contains('her2', na=False)).sum()
    other_count = len(adata) - her2_count - low_count
    
    print(f"匹配到 'Her2' 模式的细胞数量: {her2_count}")
    print(f"匹配到 'low' 模式的细胞数量: {low_count}")
    print(f"未匹配到特定模式的细胞数量: {other_count}")
    print(f"group_new列的唯一值: {adata.obs['group_new'].unique()}")
    
    # 验证结果
    print("\n=== 数据验证（前5行）===")
    print(adata.obs[['group', 'sample', 'group_new']].head())
    
    if output_file:
        print(f"\n正在保存文件到: {output_file}")
        adata.write(output_file, compression='gzip')
        print("文件保存成功！")
    
    return adata

# 替代方案：使用字符串方法（更简洁）
def add_group_new_column_alternative(h5ad_file, output_file=None):
    """
    替代方案：使用pandas字符串方法
    """
    adata = sc.read_h5ad(h5ad_file)
    
    # 确保数据类型正确
    adata.obs['group'] = adata.obs['group'].astype(str)
    adata.obs['sample'] = adata.obs['sample'].astype(str)
    
    # 使用str.contains方法（需要处理NaN值）
    sample_lower = adata.obs['sample'].str.lower()
    
    # 创建条件掩码
    her2_mask = sample_lower.str.contains('her2', na=False)
    low_mask = sample_lower.str.contains('low', na=False) & ~her2_mask
    
    # 应用条件
    adata.obs['group_new'] = adata.obs['group'].copy()
    adata.obs.loc[her2_mask, 'group_new'] = adata.obs.loc[her2_mask, 'group'] + '_Her2_0'
    adata.obs.loc[low_mask, 'group_new'] = adata.obs.loc[low_mask, 'group'] + '_Her2_low'
    
    return adata

def main():
    """
    主执行函数
    """
    input_file = "merged_single_cell_data.h5ad"
    output_file = "merged_single_cell_data_with_group_new.h5ad"
    
    try:
        # 使用修复后的函数
        adata = add_group_new_column(input_file, output_file)
        
        if adata is not None:
            print("\n=== 处理完成 ===")
            print(f"总细胞数: {adata.n_obs}")
            print(f"group_new列统计:")
            print(adata.obs['group_new'].value_counts())
            
    except Exception as e:
        print(f"处理过程中出错: {e}")
        print("尝试使用替代方案...")
        
        # 尝试替代方案
        try:
            adata = add_group_new_column_alternative(input_file, output_file)
            if adata is not None:
                print("替代方案执行成功！")
        except Exception as e2:
            print(f"替代方案也失败: {e2}")

if __name__ == "__main__":
    main()
