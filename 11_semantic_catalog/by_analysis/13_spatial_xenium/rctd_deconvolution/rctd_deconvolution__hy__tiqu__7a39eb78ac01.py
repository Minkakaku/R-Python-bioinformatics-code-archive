import scanpy as sc
import pandas as pd
import numpy as np

def process_h5ad_with_releate(h5ad_path, releate_path, output_path=None):
    """
    根据releate文件处理h5ad数据
    
    Parameters:
    -----------
    h5ad_path : str
        h5ad文件路径
    releate_path : str 
        releate文件路径（第一列patient，第二列T分组）
    output_path : str, optional
        输出文件路径，如不指定则不保存
    """
    
    # 1. 读取h5ad文件[4,6,8](@ref)
    print("正在读取h5ad文件...")
    adata = sc.read_h5ad(h5ad_path)
    print(f"原始数据形状: {adata.shape}")
    print(f"观察变量(obs)列: {list(adata.obs.columns)}")
    
    # 2. 读取releate文件[1](@ref)
    print("\n正在读取releate文件...")
    # 尝试不同的分隔符读取
    try:
        # 先尝试制表符分隔
        releate_df = pd.read_csv(releate_path, sep='\t', header=None, names=['patient', 'T_group'])
    except:
        try:
            # 尝试逗号分隔
            releate_df = pd.read_csv(releate_path, sep=',', header=None, names=['patient', 'T_group'])
        except:
            # 尝试空格分隔
            releate_df = pd.read_csv(releate_path, sep='\s+', header=None, names=['patient', 'T_group'])
    
    print(f"releate文件内容:\n{releate_df}")
    print(f"releate文件中唯一的patient数量: {len(releate_df['patient'].unique())}")
    
    # 3. 检查h5ad数据中是否有patient信息[5,8](@ref)
    if 'patient' not in adata.obs.columns:
        raise ValueError("h5ad文件的obs中未找到'patient'列，请检查数据格式")
    
    # 4. 筛选patient在releate文件中的细胞[8](@ref)
    print("\n正在筛选细胞...")
    patient_in_releate = releate_df['patient'].unique()
    
    # 确保patient类型一致（字符串比较）
    adata.obs['patient'] = adata.obs['patient'].astype(str)
    patient_in_releate = [str(p) for p in patient_in_releate]
    
    # 创建筛选掩码
    mask = adata.obs['patient'].isin(patient_in_releate)
    adata_filtered = adata[mask, :].copy()
    
    print(f"筛选后数据形状: {adata_filtered.shape}")
    print(f"筛选掉的细胞数量: {adata.shape[0] - adata_filtered.shape[0]}")
    
    # 5. 添加MPCpercent分组信息[5](@ref)
    print("\n正在添加MPCpercent分组信息...")
    
    # 创建patient到T分组的映射字典
    patient_to_group = dict(zip(
        releate_df['patient'].astype(str), 
        releate_df['T_group'].astype(str)
    ))
    
    # 将分组信息添加到obs中
    adata_filtered.obs['MPCpercent'] = adata_filtered.obs['patient'].map(patient_to_group)
    
    # 检查是否有未映射的细胞
    unmapped_cells = adata_filtered.obs['MPCpercent'].isna().sum()
    if unmapped_cells > 0:
        print(f"警告: 有{unmapped_cells}个细胞未能成功映射分组信息")
        # 删除未映射的细胞
        adata_filtered = adata_filtered[~adata_filtered.obs['MPCpercent'].isna(), :].copy()
        print(f"清理后数据形状: {adata_filtered.shape}")
    
    # 6. 统计信息
    print("\n=== 处理结果统计 ===")
    print(f"最终细胞数量: {adata_filtered.shape[0]}")
    print(f"最终基因数量: {adata_filtered.shape[1]}")
    print(f"包含的patient数量: {adata_filtered.obs['patient'].nunique()}")
    print(f"MPCpercent分组分布:")
    print(adata_filtered.obs['MPCpercent'].value_counts())
    
    # 7. 保存结果[5](@ref)
    if output_path:
        print(f"\n正在保存结果到: {output_path}")
        adata_filtered.write_h5ad(output_path)
        print("保存完成!")
    
    return adata_filtered

def validate_data_integrity(adata_filtered):
    """
    验证数据完整性
    """
    print("\n=== 数据完整性验证 ===")
    
    # 检查基本属性
    print(f"数据形状: {adata_filtered.shape}")
    print(f"观察变量列: {list(adata_filtered.obs.columns)}")
    print(f"变量特征列: {list(adata_filtered.var.columns)}")
    
    # 检查X矩阵数据类型
    print(f"X矩阵数据类型: {type(adata_filtered.X)}")
    print(f"X矩阵形状: {adata_filtered.X.shape}")
    
    # 检查稀疏矩阵情况
    from scipy.sparse import issparse
    if issparse(adata_filtered.X):
        print("X矩阵为稀疏矩阵格式")
        print(f"非零元素数量: {adata_filtered.X.nnz}")
    else:
        print("X矩阵为稠密矩阵格式")

# 主程序
if __name__ == "__main__":
    # 文件路径配置
    h5ad_file = "adata_celltype_RCTD_removed_with_patient.h5ad"  # 替换为您的h5ad文件路径
    releate_file = "releate"  # 替换为您的releate文件路径
    output_file = "adata_filtered_with_MPCpercent.h5ad"  # 输出文件路径
    
    try:
        # 执行处理
        adata_result = process_h5ad_with_releate(
            h5ad_path=h5ad_file,
            releate_path=releate_file, 
            output_path=output_file
        )
        
        # 数据验证
        validate_data_integrity(adata_result)
        
        print("\n=== 处理完成 ===")
        print("下一步建议:")
        print("1. 使用 sc.read_h5ad('adata_filtered_with_MPCpercent.h5ad') 读取处理后的数据")
        print("2. 使用 adata.obs['MPCpercent'] 访问分组信息进行后续分析")
        print("3. 可使用 sc.pl.umap(adata, color=['MPCpercent', 'cell_type']) 可视化结果")
        
    except FileNotFoundError as e:
        print(f"文件未找到错误: {e}")
        print("请检查文件路径是否正确")
    except Exception as e:
        print(f"处理过程中出现错误: {e}")
        print("请检查数据格式是否符合要求")
