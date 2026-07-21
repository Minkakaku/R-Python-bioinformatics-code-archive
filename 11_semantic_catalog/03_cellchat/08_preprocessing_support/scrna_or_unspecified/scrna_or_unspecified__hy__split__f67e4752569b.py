import anndata as ad
import os

def extract_cn_subgroups(input_path, output_dir):
    """
    从h5ad文件中提取指定CNs_10类别的数据子集
    
    参数：
    input_path: 输入h5ad文件路径
    output_dir: 输出目录
    """
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取h5ad数据
    adata = ad.read_h5ad(input_path)
    
    # 检查必需的列是否存在
    if 'sub_CNs_10' not in adata.obs.columns:
        raise ValueError("obs中不存在'CNs_10'列！")
    
    # 定义要提取的类别和对应文件名
    cn_groups = {
        "sub_CN1": "sub_CN1",
        "sub_CN2": "sub_CN2",
        "sub_CN3": "sub_CN3",
        "sub_CN5": "sub_CN5",
        "sub_CN7": "sub_CN7",
        "sub_CN9": "sub_CN9",
    }
    
    # 提取并保存每个子集
    for cn_value, file_prefix in cn_groups.items():
        # 创建逻辑索引 - 兼容整数和字符串类型
        mask = adata.obs['sub_CNs_10'].astype(str) == str(cn_value)
        
        # 检查是否有足够的数据点
        if mask.sum() == 0:
            print(f"警告: 未找到CNs_10值为{cn_value}的细胞")
            continue
            
        # 创建子集
        subset = adata[mask].copy()
        
        # 设置输出文件名
        filename = os.path.join(output_dir, f"{file_prefix}.h5ad")
        
        # 保存数据
        subset.write(filename)
        print(f"已保存: {filename} (包含{mask.sum()}个细胞)")

if __name__ == "__main__":
    # 配置路径 - 替换为您的实际路径
    input_h5ad = "new.h5ad"
    output_directory = "."
    
    # 执行提取
    extract_cn_subgroups(input_h5ad, output_directory)
    print("处理完成！")
