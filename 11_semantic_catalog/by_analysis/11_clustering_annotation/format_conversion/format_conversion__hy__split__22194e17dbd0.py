import scanpy as sc
import os
import pandas as pd
import argparse

def split_anndata_by_celltype(input_file, output_dir="./cell_type_subsets"):
    """
    按cell_type列分割单细胞数据，保存为多个h5ad文件，并保留raw数据
    
    参数:
    input_file (str): 输入h5ad文件路径
    output_dir (str): 输出目录路径，默认为'./cell_type_subsets'
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 加载数据
    print(f"正在加载数据: {input_file}")
    adata = sc.read_h5ad(input_file)
    print(f"原始数据维度: {adata.n_obs} 细胞, {adata.n_vars} 基因")
    
    # 检查cell_type列是否存在
    if 'cell_type' not in adata.obs.columns:
        raise ValueError("错误: adata.obs中不存在'cell_type'列。请检查数据。")
    
    # 检查是否有raw数据
    has_raw = adata.raw is not None
    if has_raw:
        print("检测到raw数据，将自动保留在子集中")
    else:
        print("未检测到raw数据，将仅保存处理后的数据")
    
    # 获取唯一的cell_type值
    unique_cell_types = adata.obs['cell_type'].unique()
    print(f"发现 {len(unique_cell_types)} 个细胞类型: {list(unique_cell_types)}")
    
    # 遍历每个cell_type并保存子集
    for cell_type in unique_cell_types:
        # 创建子集（选择特定cell_type的细胞）
        cell_mask = adata.obs['cell_type'] == cell_type
        adata_subset = adata[cell_mask].copy()
        
        # 生成安全的文件名（替换特殊字符）
        safe_name = str(cell_type).replace(" ", "_").replace("/", "_").replace("\\", "_")
        filename = f"cell_type_{safe_name}.h5ad"
        filepath = os.path.join(output_dir, filename)
        
        # 保存子集
        adata_subset.write(filepath)
        print(f"保存 {cell_type}:\t{adata_subset.n_obs} 个细胞 -> {filepath}")
    
    print(f"\n分割完成！所有文件已保存到: {output_dir}")
    print(f"每个文件包含对应细胞类型的子集，并保留了原始数据结构和raw数据。")

if __name__ == "__main__":
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="按cell_type分割单细胞数据")
    parser.add_argument("-i", "--input", required=True, help="输入h5ad文件路径")
    parser.add_argument("-o", "--output", default="./cell_type_subsets", help="输出目录路径")
    
    args = parser.parse_args()
    
    # 运行主函数
    split_anndata_by_celltype(args.input, args.output)
