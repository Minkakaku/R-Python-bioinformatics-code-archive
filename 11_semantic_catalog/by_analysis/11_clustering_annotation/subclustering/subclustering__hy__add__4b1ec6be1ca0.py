import scanpy as sc
import pandas as pd
import os

def integrate_subcluster_annotations(main_h5ad_path, subcluster_dir, output_path):
    """
    将亚群h5ad的celltype_sub信息整合到总h5ad文件
    
    参数:
        main_h5ad_path: 总h5ad文件路径
        subcluster_dir: 存放亚群h5ad的目录
        output_path: 输出文件路径
    """
    # 1. 加载总h5ad文件
    adata_main = sc.read_h5ad(main_h5ad_path)
    print(f"总数据集细胞数: {adata_main.n_obs}")
    
    # 2. 初始化存储字典
    subcluster_dict = {}
    
    # 3. 遍历亚群目录
    for filename in os.listdir(subcluster_dir):
        if filename.endswith(".h5ad"):
            filepath = os.path.join(subcluster_dir, filename)
            print(f"处理亚群文件: {filename}")
            
            # 加载亚群数据
            adata_sub = sc.read_h5ad(filepath)
            
            # 提取celltype_sub信息
            if 'celltype_sub' in adata_sub.obs.columns:
                # 添加到字典 {细胞ID: 亚群类型}
                subcluster_dict.update(adata_sub.obs['celltype_sub'].to_dict())
            else:
                print(f"警告: {filename} 中未找到celltype_sub列")
    
    # 4. 创建新的注释列
    adata_main.obs['celltype_final'] = "unassigned"  # 初始化
    
    # 5. 匹配细胞ID并赋值
    matched = 0
    for cell_id in adata_main.obs_names:
        if cell_id in subcluster_dict:
            adata_main.obs.at[cell_id, 'celltype_final'] = subcluster_dict[cell_id]
            matched += 1
    
    print(f"成功匹配细胞数: {matched}/{adata_main.n_obs} ({matched/adata_main.n_obs:.1%})")
    
    # 6. 保存结果
    adata_main.write(output_path)
    print(f"结果已保存至: {output_path}")
    return adata_main

# 使用示例
if __name__ == "__main__":
    # 配置路径
    main_file = "/home/hy/Project/FWSP20251095/Redo/output/data/adata_celltype.h5ad"        # 总h5ad文件路径
    sub_dir = "/home/hy/Project/FWSP20251095/Redo/subcluster/sub/combine/"            # 亚群文件目录
    output_file = "integrated_data.h5ad"  # 输出文件路径
    
    # 执行整合
    integrated_adata = integrate_subcluster_annotations(main_file, sub_dir, output_file)
    
    # 验证结果
    print("\n细胞类型分布:")
    print(integrated_adata.obs['celltype_final'].value_counts())
