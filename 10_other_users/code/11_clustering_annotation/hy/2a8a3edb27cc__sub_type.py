import scanpy as sc
import os

def split_h5ad_by_celltype(input_path, output_dir, celltype_column='cell_type'):
    """
    功能：按细胞类型拆分h5ad文件
    参数：
        input_path: 输入h5ad文件路径
        output_dir: 输出目录路径
        celltype_column: obs中存储细胞类型的列名（需替换为实际列名）
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取h5ad文件
    adata = sc.read_h5ad(input_path)
    
    try:
        # 检查细胞类型列是否存在
        if celltype_column not in adata.obs.columns:
            raise ValueError(f"列名'{celltype_column}'不存在于obs中，请检查可用列：{adata.obs.columns.tolist()}")
        
        # 获取唯一细胞类型列表
        unique_celltypes = adata.obs[celltype_column].unique()
        print(f"发现{len(unique_celltypes)}种细胞类型")
        
        # 遍历每个细胞类型
        for ct in unique_celltypes:
            # 过滤当前细胞类型
            mask = adata.obs[celltype_column] == ct
            subset = adata[mask].copy()
            
            # 生成安全文件名（替换特殊字符）
            safe_ct_name = ct.replace(" ", "_").replace("/", "-")
            output_path = os.path.join(output_dir, f"{safe_ct_name}.h5ad")
            
            # 保存子集
            subset.write(output_path)
            print(f"保存 {ct}: {subset.n_obs} 个细胞 -> {output_path}")
            
    except Exception as e:
        print(f"错误发生：{str(e)}")

# 使用示例
if __name__ == "__main__":
    split_h5ad_by_celltype(
        input_path="adata_celltype.h5ad",  # 替换为实际路径
        output_dir=".",         # 输出目录
        celltype_column="celltype"             # 根据实际列名修改
    )
