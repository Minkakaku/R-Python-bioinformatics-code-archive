import pandas as pd
import scanpy as sc

def extract_cell_info(adata, output_file="cell_info_res0.4.csv"):
    """
    从AnnData对象中提取细胞名称和res_0.6列信息
    
    参数:
    adata: AnnData对象
    output_file: 输出文件名
    """
    
    # 1. 检查res_0.6列是否存在
    if 'res_0.4' not in adata.obs.columns:
        print("错误: adata.obs中不存在'res_0.6'列")
        print("可用的列有:", list(adata.obs.columns))
        return None
    
    # 2. 提取细胞名称和res_0.6列
    cell_info = pd.DataFrame({
        'cell_id': adata.obs_names,
        'res_0.4': adata.obs['res_0.4']
    })
    
    # 3. 显示基本信息
    print(f"提取了 {len(cell_info)} 个细胞的信息")
    print(f"res_0.6的聚类情况:")
    print(cell_info['res_0.4'].value_counts().sort_index())
    
    # 4. 保存为CSV文件
    cell_info.to_csv(output_file, index=False)
    print(f"细胞信息已保存到: {output_file}")
    
    return cell_info

# 使用示例
if __name__ == "__main__":
    # 加载您的数据（请替换为实际路径）
    adata = sc.read_h5ad('adata_cluster.h5ad')
    
    # 提取细胞信息
    cell_data = extract_cell_info(adata)
    
    # 可选：查看前几行数据
    if cell_data is not None:
        print("\n前5行数据预览:")
        print(cell_data.head())
