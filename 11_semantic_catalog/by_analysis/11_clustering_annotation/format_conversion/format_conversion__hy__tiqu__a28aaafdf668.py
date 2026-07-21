import pandas as pd
import scanpy as sc
import logging

# 设置日志记录
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def extract_cell_data(h5ad_path, output_csv=None):
    """
    从H5AD文件中安全提取cell_id和res_0.4列
    
    参数:
        h5ad_path: H5AD文件路径
        output_csv: 输出CSV文件路径(可选)
    
    返回:
        pandas.DataFrame: 提取的数据框
    """
    try:
        # 读取H5AD文件
        logger.info(f"读取文件: {h5ad_path}")
        adata = sc.read_h5ad(h5ad_path)
        
        # 检查必要列是否存在
        required_cols = ['cell_id', 'res_0.4']
        missing_cols = [col for col in required_cols if col not in adata.obs.columns]
        
        if missing_cols:
            available_cols = list(adata.obs.columns)
            raise ValueError(f"缺少必要列 {missing_cols}。可用列: {available_cols}")
        
        # 提取数据
        extracted_data = adata.obs[required_cols].copy()
        logger.info(f"成功提取 {len(extracted_data)} 行数据")
        
        # 处理分类数据问题（避免.min()/.max()错误）
        if hasattr(extracted_data['res_0.4'], 'cat'):
            logger.info("res_0.4列为分类数据类型")
            # 显示分类信息而不计算极值
            unique_values = extracted_data['res_0.4'].unique()
            logger.info(f"res_0.4唯一值: {list(unique_values)}")
        
        # 设置输出路径
        if output_csv is None:
            import os
            base_name = os.path.splitext(h5ad_path)[0]
            output_csv = f"{base_name}_extracted.csv"
        
        # 保存为CSV
        extracted_data.to_csv(output_csv, index=False)
        logger.info(f"数据已保存至: {output_csv}")
        
        return extracted_data
        
    except Exception as e:
        logger.error(f"处理失败: {str(e)}")
        raise

# 使用示例
if __name__ == "__main__":
    # 示例用法
    input_file = "adata_cluster.h5ad"  # 替换为您的文件路径
    output_file = "extracted_data.csv"  # 可选输出路径
    
    try:
        data = extract_cell_data(input_file, output_file)
        print("提取成功！前5行数据:")
        print(data.head())
    except Exception as e:
        print(f"错误: {e}")
