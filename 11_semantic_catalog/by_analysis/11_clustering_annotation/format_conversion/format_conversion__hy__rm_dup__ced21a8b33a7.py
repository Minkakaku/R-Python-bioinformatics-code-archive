import scanpy as sc
import pandas as pd
import anndata as ad

def remove_duplicate_cells_by_cell_id(h5ad_file, output_file=None):
    """
    删除h5ad文件中cell_id列有重复的细胞，保留每个cell_id的第一次出现
    
    参数:
    h5ad_file: 输入的h5ad文件路径
    output_file: 输出的h5ad文件路径，如果为None则不保存
    """
    # 1. 读取h5ad文件
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_file)
    print(f"原始数据形状: {adata.shape} (细胞数 x 基因数)")
    
    # 2. 检查cell_id列是否存在
    if 'cell_id' not in adata.obs.columns:
        raise KeyError("在adata.obs中未找到'cell_id'列，请检查数据。")
    
    # 3. 检查重复的cell_id
    duplicate_mask = adata.obs.duplicated(subset=['cell_id'], keep='first')
    duplicate_count = duplicate_mask.sum()
    
    print(f"发现 {duplicate_count} 个重复的cell_id（将删除这些细胞）")
    print(f"唯一cell_id数量: {adata.obs['cell_id'].nunique()}")
    
    if duplicate_count == 0:
        print("没有发现重复的cell_id，无需处理。")
        if output_file:
            adata.write_h5ad(output_file)
            print(f"数据已保存至: {output_file}")
        return adata
    
    # 4. 显示重复cell_id的详细信息
    duplicate_cell_ids = adata.obs.loc[duplicate_mask, 'cell_id'].unique()
    print(f"\n重复的cell_id示例（前10个）:")
    for i, cell_id in enumerate(duplicate_cell_ids[:10]):
        duplicate_indices = adata.obs[adata.obs['cell_id'] == cell_id].index.tolist()
        print(f"  {cell_id}: 出现在 {len(duplicate_indices)} 个细胞中")
        if i >= 4:  # 只显示前5个的详细信息
            break
    
    # 5. 删除重复的细胞（保留每个cell_id的第一次出现）
    print("\n正在删除重复的细胞...")
    
    # 方法1: 使用布尔索引（推荐）
    unique_mask = ~duplicate_mask
    adata_unique = adata[unique_mask, :]
    
    # 方法2: 也可以使用drop_duplicates（需要重建AnnData对象）
    # unique_obs = adata.obs.drop_duplicates(subset=['cell_id'], keep='first')
    # adata_unique = adata[unique_obs.index, :]
    
    print(f"去重后数据形状: {adata_unique.shape} (细胞数 x 基因数)")
    print(f"删除了 {adata.n_obs - adata_unique.n_obs} 个重复细胞")
    
    # 6. 验证去重结果
    remaining_duplicates = adata_unique.obs.duplicated(subset=['cell_id']).sum()
    if remaining_duplicates == 0:
        print("✅ 验证通过: 所有cell_id重复项已成功删除")
    else:
        print(f"⚠️ 警告: 仍有 {remaining_duplicates} 个重复cell_id存在")
    
    # 7. 保存结果
    if output_file:
        print(f"\n正在保存去重后的数据到: {output_file}")
        adata_unique.write_h5ad(output_file)
        print("保存完成！")
    else:
        print("\n注意: 未指定输出文件，去重结果未保存")
    
    return adata_unique

def remove_duplicates_advanced(h5ad_file, output_file=None, keep='first'):
    """
    高级去重功能：可以选择保留第一个或最后一个重复项
    
    参数:
    h5ad_file: 输入的h5ad文件路径
    output_file: 输出的h5ad文件路径
    keep: 保留策略 'first'（保留第一个）或 'last'（保留最后一个）
    """
    adata = sc.read_h5ad(h5ad_file)
    print(f"原始数据形状: {adata.shape}")
    
    if 'cell_id' not in adata.obs.columns:
        raise KeyError("在adata.obs中未找到'cell_id'列")
    
    # 使用drop_duplicates方法，可以灵活选择保留策略
    unique_obs = adata.obs.drop_duplicates(subset=['cell_id'], keep=keep)
    adata_unique = adata[unique_obs.index, :]
    
    print(f"去重后数据形状: {adata_unique.shape}")
    print(f"删除策略: 保留每个cell_id的{keep}次出现")
    
    if output_file:
        adata_unique.write_h5ad(output_file)
        print(f"已保存到: {output_file}")
    
    return adata_unique

if __name__ == "__main__":
    # 设置文件路径
    input_file = "merged_celltypes.h5ad"
    output_file = "merged_celltypes_unique.h5ad"
    
    try:
        # 基本去重（保留第一个出现的细胞）
        result_adata = remove_duplicate_cells_by_cell_id(input_file, output_file)
        
        # 如果需要高级功能，可以取消注释下面的代码
        # result_adata = remove_duplicates_advanced(input_file, output_file, keep='first')
        # result_adata = remove_duplicates_advanced(input_file, output_file, keep='last')
        
        print("\n" + "="*50)
        print("去重处理完成!")
        print("="*50)
        
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
