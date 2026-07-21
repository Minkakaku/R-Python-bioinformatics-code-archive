import scanpy as sc
import os
import pandas as pd
import argparse

def split_anndata_by_celltype(input_file, output_dir="./celltype_subsets"):
    """
    按celltype列分割单细胞数据，保存为多个h5ad文件，并保留raw数据
    
    参数:
    input_file (str): 输入h5ad文件路径
    output_dir (str): 输出目录路径，默认为'./celltype_subsets'
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 加载数据
    print(f"正在加载数据: {input_file}")
    adata = sc.read_h5ad(input_file)
    print(f"原始数据维度: {adata.n_obs} 细胞, {adata.n_vars} 基因")
    
    # 检查celltype列是否存在
    if 'celltype' not in adata.obs.columns:
        raise ValueError("错误: adata.obs中不存在'celltype'列。请检查数据。")
    
    # 检查是否有raw数据
    has_raw = adata.raw is not None
    if has_raw:
        print("检测到raw数据，将自动保留在子集中")
    else:
        print("未检测到raw数据，将仅保存处理后的数据")
    
    # 获取唯一的celltype值
    unique_celltypes = adata.obs['celltype'].unique()
    print(f"发现 {len(unique_celltypes)} 个细胞类型: {list(unique_celltypes)}")
    
    # 遍历每个celltype并保存子集
    for celltype in unique_celltypes:
        # 创建子集（选择特定celltype的细胞）
        cell_mask = adata.obs['celltype'] == celltype
        adata_subset = adata[cell_mask].copy()
        
        # 生成安全的文件名（替换特殊字符）
        safe_name = str(celltype).replace(" ", "_").replace("/", "_").replace("\\", "_")
        filename = f"celltype_{safe_name}.h5ad"
        filepath = os.path.join(output_dir, filename)
        
        # 保存子集
        adata_subset.write(filepath)
        print(f"保存 {celltype}:\t{adata_subset.n_obs} 个细胞 -> {filepath}")
    
    print(f"\n分割完成！所有文件已保存到: {output_dir}")
    print(f"每个文件包含对应细胞类型的子集，并保留了原始数据结构和raw数据。")

if __name__ == "__main__":
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="按celltype分割单细胞数据")
    parser.add_argument("-i", "--input", required=True, help="输入h5ad文件路径")
    parser.add_argument("-o", "--output", default="./celltype_subsets", help="输出目录路径")
    
    args = parser.parse_args()
    
    # 运行主函数
    split_anndata_by_celltype(args.input, args.output)
