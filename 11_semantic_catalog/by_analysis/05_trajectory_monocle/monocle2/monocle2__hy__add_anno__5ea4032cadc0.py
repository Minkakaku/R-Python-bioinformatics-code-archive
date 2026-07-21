import scanpy as sc
import pandas as pd
import argparse

def add_annotation_to_h5ad(h5ad_file, meta_file, output_file=None):
    """
    将meta.file中的annotation信息添加到h5ad文件的obs中
    
    参数:
    h5ad_file: h5ad文件路径
    meta_file: meta.file文件路径
    output_file: 输出文件路径(可选，如不指定则覆盖原文件)
    """
    # 读取h5ad文件 [4,6](@ref)
    print("读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_file)
    print(f"原始数据形状: {adata.shape}")
    print(f"obs列名: {list(adata.obs.columns)}")
    
    # 读取meta.file文件(无表头) [11](@ref)
    print("\n读取meta.file文件...")
    meta_df = pd.read_csv(meta_file, header=None, 
                         names=['cell_id', 'sample', 'annotation'])
    print(f"meta.file行数: {len(meta_df)}")
    print(f"meta.file前几行:\n{meta_df.head()}")
    
    # 检查h5ad的obs中是否已存在cell_id和sample列
    if 'cell_id' not in adata.obs.columns or 'sample' not in adata.obs.columns:
        raise ValueError("h5ad文件的obs中必须包含'cell_id'和'sample'列")
    
    # 备份原始的obs索引
    original_index = adata.obs.index
    
    # 将cell_id和sample作为临时索引用于合并
    adata.obs['temp_index'] = adata.obs['cell_id'].astype(str) + '_' + adata.obs['sample'].astype(str)
    meta_df['temp_index'] = meta_df['cell_id'].astype(str) + '_' + meta_df['sample'].astype(str)
    
    # 在meta_df中去除可能的重复项(保留第一个) [11](@ref)
    meta_df = meta_df.drop_duplicates(subset=['temp_index'], keep='first')
    
    # 创建用于合并的临时DataFrame [9,10](@ref)
    merge_df = pd.DataFrame(index=adata.obs.index)
    merge_df['temp_index'] = adata.obs['temp_index']
    
    # 使用merge进行数据合并 [10,11](@ref)
    merged_obs = pd.merge(merge_df, meta_df[['temp_index', 'annotation']], 
                         on='temp_index', how='left')
    
    # 将合并后的annotation列添加到原始obs中
    adata.obs['annotation'] = merged_obs['annotation'].values
    
    # 清理临时列
    adata.obs.drop('temp_index', axis=1, inplace=True)
    
    # 统计匹配情况
    matched_count = adata.obs['annotation'].notna().sum()
    print(f"\n匹配统计:")
    print(f"总细胞数: {adata.n_obs}")
    print(f"成功匹配annotation的细胞数: {matched_count}")
    print(f"未匹配的细胞数: {adata.n_obs - matched_count}")
    
    # 显示annotation的分布
    if matched_count > 0:
        print(f"\nAnnotation分布:")
        print(adata.obs['annotation'].value_counts(dropna=False))
    
    # 保存结果 [3](@ref)
    if output_file is None:
        output_file = h5ad_file
    
    adata.write(output_file, compression="gzip")
    print(f"\n结果已保存至: {output_file}")
    
    return adata

if __name__ == "__main__":
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='将meta.file中的annotation添加到h5ad文件')
    parser.add_argument('--h5ad', required=True, help='输入h5ad文件路径')
    parser.add_argument('--meta', required=True, help='meta.file文件路径')
    parser.add_argument('--output', help='输出文件路径(可选)')
    
    args = parser.parse_args()
    
    try:
        result = add_annotation_to_h5ad(args.h5ad, args.meta, args.output)
        print("处理完成!")
    except Exception as e:
        print(f"处理过程中出现错误: {str(e)}")
