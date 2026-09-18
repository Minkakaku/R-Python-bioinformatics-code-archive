import scanpy as sc
import pandas as pd

# 1. 读取h5ad文件，请将路径替换为您的实际文件路径
adata = sc.read_h5ad('merged_single_cell_data_with_group_new.h5ad')

# 2. 定义要提取的样本列表
samples_to_extract = ["Her2_0_11_1","Her2_0_11_2","Her2_0_11_3","Her2_0_3","Her2_0_4_1","Her2_0_4_2","Her2_0_6_1","Her2_0_6_2","low_10_1","low_10_2","low_4_1","low_4_2","low_6","low_7_1","low_7_2","low_9_1","low_9_2","Her2_017","Her2_019","Her2_020","Her2_022","low_16","low_17","low_21","low_22","low_23"]
# 3. 检查'sample'列是否存在
if 'sample' not in adata.obs.columns:
    print("错误：在adata.obs中未找到'sample'列")
    print("可用的列有：", adata.obs.columns.tolist())
else:
    # 4. 检查请求的样本是否存在于数据中
    existing_samples = set(adata.obs['sample'])
    samples_found = [s for s in samples_to_extract if s in existing_samples]
    samples_not_found = [s for s in samples_to_extract if s not in existing_samples]
    
    if samples_not_found:
        print("警告：以下样本在数据中未找到，将被忽略：")
        print(samples_not_found)
    
    if samples_found:
        # 5. 使用布尔索引提取样本子集
        # 这将保留所有与这些样本相关的细胞
        adata_subset = adata[adata.obs['sample'].isin(samples_found), :]
        
        # 6. 显示提取结果
        print(f"\n=== 提取成功 ===")
        print(f"原始数据细胞数: {adata.n_obs}")
        print(f"提取后细胞数: {adata_subset.n_obs}")
        print(f"基因数: {adata_subset.n_vars}")
        
        # 显示提取的样本及其细胞数量
        print("\n=== 提取的样本统计 ===")
        sample_counts = adata_subset.obs['sample'].value_counts()
        for sample, count in sample_counts.items():
            print(f"{sample}: {count} 个细胞")
        
        # 7. 保存提取后的数据（可选）
        output_file = "extracted_samples.h5ad"
        adata_subset.write(output_file)
        print(f"\n提取的数据已保存为: {output_file}")
        
        # 8. 验证保存的文件
        verify_adata = sc.read_h5ad(output_file)
        print(f"验证: 保存的文件包含 {verify_adata.n_obs} 个细胞")
    else:
        print("错误：未找到任何指定的样本")
        print("数据中存在的样本示例：", list(existing_samples)[:10])
