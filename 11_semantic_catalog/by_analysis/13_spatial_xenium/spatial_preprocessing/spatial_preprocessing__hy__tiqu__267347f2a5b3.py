import scanpy as sc
import pandas as pd

# 1. 读取h5ad文件，请将路径替换为您的实际文件路径
adata = sc.read_h5ad('merged_TMA_data.h5ad')

# 2. 定义要提取的样本列表
samples_to_extract = ["Z2038176-15","Q2324341-1","Z1530765-11","F1706967-11","F1805049-2","Z1915172-7","Z1954693-9","ZF2207170-8","ZF2202916-8","ZF2208790-1","Z2262137-5","Z2271577-12","ZF2304103-8","Z2323818-8","Z2336500-6","Z2345459-7","Z2332034-5","Z2417954-8","Q2236500-1","Q2406275-2","Z2101519-1","Z2123354-1","Z1642635-14","Z2224197-10","Z1736838-15","Z1944219-14","Z1430967-9","Z2057500-35","Z2171766-3","Z2353385-10","Z1635989-12","Z1814675-2","Z1826868-6","Z1952876-9","Z1802235-10"]
# 3. 检查'patient'列是否存在
if 'patient' not in adata.obs.columns:
    print("错误：在adata.obs中未找到'patient'列")
    print("可用的列有：", adata.obs.columns.tolist())
else:
    # 4. 检查请求的样本是否存在于数据中
    existing_samples = set(adata.obs['patient'])
    samples_found = [s for s in samples_to_extract if s in existing_samples]
    samples_not_found = [s for s in samples_to_extract if s not in existing_samples]
    
    if samples_not_found:
        print("警告：以下样本在数据中未找到，将被忽略：")
        print(samples_not_found)
    
    if samples_found:
        # 5. 使用布尔索引提取样本子集
        # 这将保留所有与这些样本相关的细胞
        adata_subset = adata[adata.obs['patient'].isin(samples_found), :]
        
        # 6. 显示提取结果
        print(f"\n=== 提取成功 ===")
        print(f"原始数据细胞数: {adata.n_obs}")
        print(f"提取后细胞数: {adata_subset.n_obs}")
        print(f"基因数: {adata_subset.n_vars}")
        
        # 显示提取的样本及其细胞数量
        print("\n=== 提取的样本统计 ===")
        sample_counts = adata_subset.obs['patient'].value_counts()
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
