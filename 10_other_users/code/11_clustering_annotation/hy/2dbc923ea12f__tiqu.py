import scanpy as sc
import pandas as pd
import numpy as np

# 读取h5ad文件
adata = sc.read_h5ad('adata_celltype.h5ad')  # 替换为您的文件路径

# 定义表格数据（从图片中提取）
table_data = [
    ("Her2_0_10","Her2_0","Her2_0_K"),
    ("Her2_0_11_1","Her2_0","Her2_0_C"),
    ("Her2_0_11_2","Her2_0","Her2_0_C"),
    ("Her2_0_11_3","Her2_0","Her2_0_C"),
    ("Her2_0_1_1","Her2_0","Her2_0_K"),
    ("Her2_0_12_1","Her2_0","Her2_0_K"),
    ("Her2_0_12_2","Her2_0","Her2_0_K"),
    ("Her2_0_12_3","Her2_0","Her2_0_K"),
    ("Her2_0_12_4","Her2_0","Her2_0_K"),
    ("Her2_0_1_2","Her2_0","Her2_0_K"),
    ("Her2_0_13_1","Her2_0","Her2_0_K"),
    ("Her2_0_13_2","Her2_0","Her2_0_K"),
    ("Her2_0_1_3","Her2_0","Her2_0_K"),
    ("Her2_0_14","Her2_0","Her2_0_K"),
    ("Her2_0_15_1","Her2_0","Her2_0_K"),
    ("Her2_0_15_2","Her2_0","Her2_0_K"),
    ("Her2_0_16","Her2_0","Her2_0_K"),
    ("Her2_0_2_1","Her2_0","Her2_0_C"),
    ("Her2_0_2_2","Her2_0","Her2_0_C"),
    ("Her2_032_1","Her2_0","Her2_0_K"),
    ("Her2_032_2","Her2_0","Her2_0_K"),
    ("Her2_034","Her2_0","Her2_0_K"),
    ("Her2_035","Her2_0","Her2_0_K"),
    ("Her2_036","Her2_0","Her2_0_K"),
    ("Her2_0_3","Her2_0","Her2_0_C"),
    ("Her2_0_4_1","Her2_0","Her2_0_C"),
    ("Her2_0_4_2","Her2_0","Her2_0_C"),
    ("Her2_0_5_1","Her2_0","Her2_0_C"),
    ("Her2_0_5_2","Her2_0","Her2_0_C"),
    ("Her2_0_6_1","Her2_0","Her2_0_C"),
    ("Her2_0_6_2","Her2_0","Her2_0_C"),
    ("Her2_0_7_1","Her2_0","Her2_0_C"),
    ("Her2_0_7_2","Her2_0","Her2_0_C"),
    ("Her2_0_8_1","Her2_0","Her2_0_C"),
    ("Her2_0_8_2","Her2_0","Her2_0_C"),
    ("Her2_0_8_3","Her2_0","Her2_0_C"),
    ("Her2_0_8_4","Her2_0","Her2_0_C"),
    ("Her2_0_9_1","Her2_0","Her2_0_C"),
    ("Her2_0_9_2","Her2_0","Her2_0_C"),
    ("Her2_0_9_3","Her2_0","Her2_0_C"),
    ("low_10_1","Her2_LOW","Her2_LOW_C"),
    ("low_10_2","Her2_LOW","Her2_LOW_C"),
    ("low_1_1","Her2_LOW","Her2_LOW_K"),
    ("low_1_2","Her2_LOW","Her2_LOW_K"),
    ("low_12","Her2_LOW","Her2_LOW_K"),
    ("low_13_1","Her2_LOW","Her2_LOW_C"),
    ("low_13_2","Her2_LOW","Her2_LOW_C"),
    ("low_13_3","Her2_LOW","Her2_LOW_C"),
    ("low_1_3","Her2_LOW","Her2_LOW_K"),
    ("low_14_1","Her2_LOW","Her2_LOW_C"),
    ("low_14_2","Her2_LOW","Her2_LOW_C"),
    ("low_14_3","Her2_LOW","Her2_LOW_C"),
    ("low_14_4","Her2_LOW","Her2_LOW_C"),
    ("low_14_5","Her2_LOW","Her2_LOW_C"),
    ("low_2_1","Her2_LOW","Her2_LOW_K"),
    ("low_2_2","Her2_LOW","Her2_LOW_K"),
    ("low_2_3","Her2_LOW","Her2_LOW_K"),
    ("low_27_1","Her2_LOW","Her2_LOW_C"),
    ("low_27_2","Her2_LOW","Her2_LOW_C"),
    ("low_28_1","Her2_LOW","Her2_LOW_K"),
    ("low_28_2","Her2_LOW","Her2_LOW_K"),
    ("low_29_1","Her2_LOW","Her2_LOW_K"),
    ("low_29_2","Her2_LOW","Her2_LOW_K"),
    ("low_3_1","Her2_LOW","Her2_LOW_C"),
    ("low_3_2","Her2_LOW","Her2_LOW_C"),
    ("low_4_1","Her2_LOW","Her2_LOW_C"),
    ("low_4_2","Her2_LOW","Her2_LOW_C"),
    ("low_6","Her2_LOW","Her2_LOW_C"),
    ("low_7_1","Her2_LOW","Her2_LOW_C"),
    ("low_7_2","Her2_LOW","Her2_LOW_C"),
    ("low_8","Her2_LOW","Her2_LOW_C"),
    ("low_9_1","Her2_LOW","Her2_LOW_C"),
    ("low_9_2","Her2_LOW","Her2_LOW_C")
]

# 创建映射字典
g1_map = {}
g2_map = {}
for sample, g1, g2 in table_data:
    g1_map[sample] = g1
    g2_map[sample] = g2

# 检查adata.obs中是否有sample列
if 'sample' not in adata.obs.columns:
    print("错误：adata.obs中没有找到'sample'列")
    print(f"可用的列有：{list(adata.obs.columns)}")
else:
    # 1. 筛选出表格中存在的样本
    table_samples = set(g1_map.keys())
    adata_samples = set(adata.obs['sample'].unique())
    
    # 找到在表格中存在的样本
    samples_to_keep = list(table_samples.intersection(adata_samples))
    
    if len(samples_to_keep) == 0:
        print("警告：在adata.obs中没有找到表格中的任何样本")
        print(f"表格中的样本：{list(table_samples)[:10]}...")
        print(f"adata中的样本：{list(adata_samples)[:10]}...")
    else:
        print(f"找到 {len(samples_to_keep)} 个匹配的样本")
        
        # 筛选adata，只保留表格中存在的样本
        adata_filtered = adata[adata.obs['sample'].isin(samples_to_keep), :].copy()
        
        # 2. 添加g1和g2列
        adata_filtered.obs['g1'] = adata_filtered.obs['sample'].map(g1_map)
        adata_filtered.obs['g2'] = adata_filtered.obs['sample'].map(g2_map)
        
        # 检查是否有未映射的样本
        unmapped_samples = adata_filtered.obs[adata_filtered.obs['g1'].isna() | adata_filtered.obs['g2'].isna()]['sample'].unique()
        if len(unmapped_samples) > 0:
            print(f"警告：以下样本没有映射到g1或g2值：{unmapped_samples}")
        
        # 3. 保留raw矩阵
        # 检查是否有raw矩阵
        if adata_filtered.raw is not None:
            print("已保留raw矩阵")
        else:
            print("警告：原始数据没有raw矩阵")
        
        # 查看结果
        print(f"\n处理后的数据形状：{adata_filtered.shape}")
        print(f"\n样本数量统计：")
        print(adata_filtered.obs['sample'].value_counts())
        print(f"\ng1值分布：")
        print(adata_filtered.obs['g1'].value_counts())
        print(f"\ng2值分布：")
        print(adata_filtered.obs['g2'].value_counts())
        
        # 保存处理后的数据
        output_file = 'epi_with_g1_g2.h5ad'
        adata_filtered.write(output_file)
        print(f"\n处理完成，已保存为：{output_file}")
        
        # 可选：导出CSV查看映射结果
        adata_filtered.obs[['sample', 'g1', 'g2']].to_csv('sample_g1_g2_mapping.csv')
        print("映射结果已保存为：sample_g1_g2_mapping.csv")
