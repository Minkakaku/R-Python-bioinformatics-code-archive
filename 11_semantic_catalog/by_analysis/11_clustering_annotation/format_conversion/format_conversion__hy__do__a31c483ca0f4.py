import scanpy as sc
import pandas as pd
import numpy as np

def extract_celltypes_with_raw(adata, target_celltypes, output_filename="CNV.h5ad"):
    """
    从AnnData对象中提取指定细胞类型，并正确保留raw矩阵
    
    参数:
    adata: AnnData对象
    target_celltypes: 目标细胞类型列表
    output_filename: 输出文件名
    """
    
    # 1. 验证输入数据
    if 'cell_type' not in adata.obs.columns:
        raise ValueError("错误: adata.obs中不存在'cell_type'列")
    
    # 检查目标细胞类型是否存在
    available_celltypes = adata.obs['cell_type'].unique()
    target_celltypes = [ct for ct in target_celltypes if ct in available_celltypes]
    
    if not target_celltypes:
        raise ValueError("错误: 没有找到任何目标细胞类型")
    
    print(f"目标细胞类型: {target_celltypes}")
    
    # 2. 创建细胞选择掩码
    cell_mask = adata.obs['cell_type'].isin(target_celltypes)
    print(f"符合条件的细胞数: {cell_mask.sum()}/{adata.n_obs}")
    
    # 3. 提取细胞子集
    adata_subset = adata[cell_mask].copy()
    
    # 4. 正确处理raw矩阵（关键修正）
    if adata.raw is not None:
        print("检测到raw矩阵，正在处理...")
        
        # 方法1: 直接使用相同的细胞索引（最安全）
        try:
            # 确保raw对象的细胞顺序与主对象一致
            if len(adata.raw.obs_names) == len(adata.obs_names) and \
               all(adata.raw.obs_names == adata.obs_names):
                
                # 使用相同的cell_mask索引raw对象
                raw_subset = adata.raw[cell_mask].copy()
                adata_subset.raw = raw_subset
                print("成功使用直接索引法保留raw矩阵")
                
            else:
                # 方法2: 使用细胞名称匹配（当索引不完全一致时）
                subset_cell_names = adata_subset.obs_names
                raw_cell_mask = adata.raw.obs_names.isin(subset_cell_names)
                
                if raw_cell_mask.sum() == len(subset_cell_names):
                    raw_subset = adata.raw[raw_cell_mask].copy()
                    # 确保细胞顺序一致
                    raw_subset = raw_subset[subset_cell_names].copy()
                    adata_subset.raw = raw_subset
                    print("成功使用名称匹配法保留raw矩阵")
                else:
                    print("警告: raw矩阵细胞不匹配，将创建新的raw引用")
                    # 方法3: 创建新的raw引用
                    adata_subset.raw = adata_subset
                    
        except Exception as e:
            print(f"处理raw矩阵时出现错误: {e}")
            print("将使用备用方案: 创建新的raw引用")
            adata_subset.raw = adata_subset
    else:
        print("未检测到raw矩阵，将当前数据设置为raw")
        adata_subset.raw = adata_subset
    
    # 5. 验证结果
    print(f"提取后数据: {adata_subset.n_obs}个细胞, {adata_subset.n_vars}个基因")
    if adata_subset.raw is not None:
        print(f"Raw矩阵: {adata_subset.raw.n_obs}个细胞, {adata_subset.raw.n_vars}个基因")
    
    # 显示细胞类型分布
    celltype_counts = adata_subset.obs['cell_type'].value_counts()
    print("细胞类型分布:")
    for celltype, count in celltype_counts.items():
        percentage = count / adata_subset.n_obs * 100
        print(f"  {celltype}: {count}个细胞 ({percentage:.1f}%)")
    
    # 6. 保存结果
    adata_subset.write(output_filename)
    print(f"结果已保存为: {output_filename}")
    
    return adata_subset

# 使用示例
if __name__ == "__main__":
    # 请先加载您的数据
    adata1 = sc.read_h5ad('result.h5ad')
    adata = adata1.raw.to_adata()
    
    # 定义目标细胞类型
    target_celltypes = ['Epithelial cells', 'TNK cells','Endothelial cells']
    
    # 执行提取
    cnv_adata = extract_celltypes_with_raw(
        adata=adata,
        target_celltypes=target_celltypes,
        output_filename="CNV.h5ad"
    )
