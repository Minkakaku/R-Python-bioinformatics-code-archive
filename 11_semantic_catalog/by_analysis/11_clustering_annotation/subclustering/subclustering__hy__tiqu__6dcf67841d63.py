import scanpy as sc
import pandas as pd
import numpy as np
import os

def extract_patient_subset(adata_path, patient_list, output_path=None):
    """
    从h5ad文件中提取指定Patient的样本，并保留raw矩阵
    """
    
    print("开始提取指定Patient的样本数据...")
    
    # 1. 读取h5ad文件
    print(f"读取文件: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    print(f"原始数据形状: {adata.shape}")
    
    # 检查patient列名
    patient_col = None
    for col in ['patient', 'Patient', 'sample', 'Sample']:
        if col in adata.obs.columns:
            patient_col = col
            break
    
    if patient_col is None:
        raise ValueError("在obs中未找到patient、Patient、sample或Sample列")
    
    print(f"使用列: {patient_col} 作为Patient标识")
    
    # 检查raw矩阵
    if adata.raw is None:
        print("警告: 原始数据中没有raw矩阵")
        has_raw = False
    else:
        print("✅ 检测到raw矩阵")
        print(f"Raw矩阵形状: {adata.raw.X.shape}")
        has_raw = True
    
    # 2. 检查指定的Patient是否存在
    available_patients = adata.obs[patient_col].unique()
    target_patients = [str(patient).strip() for patient in patient_list]
    existing_patients = [pat for pat in target_patients if pat in adata.obs[patient_col].values]
    
    if len(existing_patients) == 0:
        raise ValueError(f"指定的Patient在数据中均未找到: {target_patients}")
    
    print(f"提取Patient: {existing_patients}")
    
    # 3. 创建细胞选择掩码
    patient_mask = adata.obs[patient_col].isin(existing_patients)
    subset_cell_indices = np.where(patient_mask)[0]
    
    # 4. 提取主数据子集
    adata_subset = adata[patient_mask, :].copy()
    print(f"提取后数据形状: {adata_subset.shape}")
    
    # 5. 处理raw矩阵（关键修正部分）
    if has_raw:
        print("处理raw矩阵的子集化...")
        
        try:
            # 方法1: 直接使用Raw对象的X矩阵创建子集
            raw_x_subset = adata.raw.X[subset_cell_indices, :]
            
            # 创建新的AnnData对象作为raw子集
            # 注意: Raw对象通常没有obs信息，只包含X和var
            raw_subset = sc.AnnData(
                X=raw_x_subset,
                var=adata.raw.var.copy()  # 保留基因信息
            )
            
            # 设置raw属性
            adata_subset.raw = raw_subset
            print("✅ Raw矩阵子集化完成")
            
        except Exception as e:
            print(f"Raw矩阵处理警告: {e}")
            print("尝试备用方法处理raw矩阵...")
            
            # 备用方法: 使用to_adata()方法
            try:
                raw_adata = adata.raw.to_adata()
                raw_subset = raw_adata[subset_cell_indices, :].copy()
                adata_subset.raw = raw_subset
                print("✅ 备用方法Raw矩阵处理完成")
            except Exception as e2:
                print(f"备用方法也失败: {e2}")
                print("将只保留主数据，不包含raw矩阵")
                adata_subset.raw = None
    
    # 6. 验证结果
    print(f"\n提取结果验证:")
    print(f"总细胞数: {adata_subset.n_obs}")
    patient_counts = adata_subset.obs[patient_col].value_counts()
    print(f"各Patient细胞数量分布:")
    for patient, count in patient_counts.items():
        print(f"  {patient}: {count} 个细胞")
    
    if adata_subset.raw is not None:
        print(f"✅ 子集raw矩阵形状: {adata_subset.raw.X.shape}")
    else:
        print("⚠️ 子集数据中没有raw矩阵")
    
    # 7. 保存结果
    if output_path is None:
        base_name = os.path.splitext(adata_path)[0]
        output_path = f"{base_name}_subset_patients.h5ad"
    
    print(f"保存结果到: {output_path}")
    adata_subset.write_h5ad(output_path)
    
    return adata_subset

# 使用示例
if __name__ == "__main__":
    # 设置参数
    adata_file = "../merged_data_with_annotations_v2.h5ad"  # 您的h5ad文件
    target_patients = ["PRE_5", "PRE_6", "PRE_8", "PRE_14", "PRE_11", "PRE_13"]
    output_file = "extracted_specific_patients.h5ad"
    
    try:
        extracted_adata = extract_patient_subset(
            adata_path=adata_file,
            patient_list=target_patients,
            output_path=output_file
        )
        
        print("\n 数据提取完成！")
        print(f"输出文件: {output_file}")
        
    except Exception as e:
        print(f"\n❌ 处理过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
