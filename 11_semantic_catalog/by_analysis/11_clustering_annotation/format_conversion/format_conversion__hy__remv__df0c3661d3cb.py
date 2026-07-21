import scanpy as sc
import pandas as pd

def modify_celltype_names(h5ad_file, output_file=None):
    """
    修改h5ad文件中celltype列的分类名称
    
    Parameters:
    -----------
    h5ad_file : str
        输入的h5ad文件路径
    output_file : str, optional
        输出的h5ad文件路径，如果为None则不保存
    """
    
    # 读取h5ad文件
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_file)
    print(f"成功读取文件: {h5ad_file}")
    print(f"数据形状: {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    
    # 检查celltype列是否存在
    if 'celltype' not in adata.obs.columns:
        print("错误: 数据中不存在 'celltype' 列")
        print(f"可用的列: {list(adata.obs.columns)}")
        return None
    
    # 显示修改前的细胞类型分布
    print("\n=== 修改前的细胞类型分布 ===")
    original_counts = adata.obs['celltype'].value_counts()
    print(original_counts)
    
    # 执行重命名操作
    print("\n=== 执行细胞类型重命名 ===")
    
    # 创建重命名映射字典
    rename_mapping = {
        'CCL28+ Stem Tumor': 'CCL28+ Tumor',
        'Tumor-Immu': 'CXCL9+ Tumor', 
        'Mixed': 'Therapy-sensitive Tumor'
    }
    
    # 应用重命名
    adata.obs['celltype'] = adata.obs['celltype'].replace(rename_mapping)
    
    # 统计修改情况
    modified_cells = 0
    for old_name, new_name in rename_mapping.items():
        count = (adata.obs['celltype'] == new_name).sum()
        print(f"'{old_name}' -> '{new_name}': {count} 个细胞")
        modified_cells += count
    
    # 显示修改后的细胞类型分布
    print("\n=== 修改后的细胞类型分布 ===")
    new_counts = adata.obs['celltype'].value_counts()
    print(new_counts)
    
    # 验证修改结果
    print("\n=== 修改验证 ===")
    for old_name in rename_mapping.keys():
        if old_name not in adata.obs['celltype'].values:
            print(f"✓ '{old_name}' 已成功替换")
        else:
            print(f"✗ '{old_name}' 仍然存在，替换可能失败")
    
    # 显示前几行数据作为验证
    print("\n=== 数据验证（前5行）===")
    print(adata.obs[['celltype']].head())
    
    # 如果需要保存文件
    if output_file:
        # 确保数据类型正确（避免保存时的TypeError问题）
        if adata.obs_names.name is None:
            adata.obs_names.name = 'cell_barcode'
        if adata.var_names.name is None:
            adata.var_names.name = 'gene_name'
        
        # 将字符串列转换为category类型以节省空间
        if adata.obs['celltype'].dtype == 'object':
            adata.obs['celltype'] = adata.obs['celltype'].astype('category')
        
        print(f"\n正在保存文件到: {output_file}")
        adata.write(output_file, compression='gzip')
        print("文件保存成功！")
    
    return adata

def check_data_quality(adata):
    """检查数据质量"""
    print("\n=== 数据质量检查 ===")
    print(f"总细胞数: {adata.n_obs}")
    print(f"总基因数: {adata.n_vars}")
    print(f"celltype列数据类型: {adata.obs['celltype'].dtype}")
    print(f"celltype列唯一值数量: {adata.obs['celltype'].nunique()}")
    print(f"celltype列唯一值: {list(adata.obs['celltype'].unique())}")
    
    # 检查是否有空值
    null_count = adata.obs['celltype'].isnull().sum()
    if null_count > 0:
        print(f"警告: celltype列有 {null_count} 个空值")

def main():
    """主执行函数"""
    # 设置文件路径
    input_file = "Epi.h5ad"
    output_file = "Epi_renamed.h5ad"  # 输出文件路径
    
    try:
        # 执行修改操作
        adata = modify_celltype_names(input_file, output_file)
        
        if adata is not None:
            # 检查数据质量
            check_data_quality(adata)
            
            # 总结报告
            print("\n" + "="*50)
            print("修改完成!")
            print("="*50)
            print(f"修改了 {adata.n_obs} 个细胞中的细胞类型标签")
            print(f"新的细胞类型分布:")
            for celltype, count in adata.obs['celltype'].value_counts().items():
                percentage = (count / adata.n_obs) * 100
                print(f"  {celltype}: {count} 个细胞 ({percentage:.1f}%)")
            
            if output_file:
                print(f"结果已保存到: {output_file}")
    
    except Exception as e:
        print(f"处理过程中出错: {e}")
        import traceback
        traceback.print_exc()

# 如果直接运行此脚本
if __name__ == "__main__":
    main()
