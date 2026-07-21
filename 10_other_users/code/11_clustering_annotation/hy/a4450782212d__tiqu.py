import scanpy as sc
import pandas as pd
import numpy as np

# 1. 读取h5ad文件
print("读取h5ad文件...")
adata = sc.read_h5ad("All_filter_state_samp_new_with_CNs10-15.h5ad")
print(f"原始数据形状: {adata.shape}")
print(f"obs列: {adata.obs.columns.tolist()}")

# 检查必需的列是否存在
if 'sample_new' not in adata.obs.columns or 'state' not in adata.obs.columns:
    print("错误: 数据中缺少 'sample_new' 或 'state' 列")
    exit(1)

# 2. 定义筛选条件
# 要筛选的sample_new值列表
target_samples = [
    "C-P1-1", "C-P2-1", "C-P2-2", "C-P3-1", 
    "C-P3-2", "C-P3-3", "C-Ind-P5", "C-P6-1", "C-P6-2"
]

# 定义T列的映射关系
t_mapping = {
    "C-P1-1": "T2",
    "C-P3-1": "T2",
    "C-P3-2": "T2",
    "C-P3-3": "T2",
    "C-P2-1": "T3",
    "C-P2-2": "T3",
    "C-Ind-P5": "T3",
    "C-P6-1": "T3",
    "C-P6-2": "T3"
}

# 3. 筛选细胞
print(f"\n筛选条件:")
print(f"- sample_new 在: {target_samples}")
print(f"- state 为: 'C'")

# 创建筛选掩码
sample_mask = adata.obs['sample_new'].isin(target_samples)
state_mask = adata.obs['state'] == 'C'
combined_mask = sample_mask & state_mask

# 应用筛选
adata_filtered = adata[combined_mask, :].copy()

print(f"筛选后的数据形状: {adata_filtered.shape}")

# 4. 检查筛选结果
if adata_filtered.shape[0] == 0:
    print("警告: 筛选后没有细胞!")
    exit(1)

# 统计每个样本的细胞数
print("\n筛选后各样本细胞数:")
sample_counts = adata_filtered.obs['sample_new'].value_counts()
for sample, count in sample_counts.items():
    print(f"  {sample}: {count}")

# 5. 添加T列
print("\n添加T列映射...")

# 初始化T列为空字符串
adata_filtered.obs['T'] = ''

# 根据映射关系填充T列
for sample, t_value in t_mapping.items():
    mask = adata_filtered.obs['sample_new'] == sample
    adata_filtered.obs.loc[mask, 'T'] = t_value

# 检查是否有未映射的样本
unmapped_samples = adata_filtered.obs[adata_filtered.obs['T'] == '']['sample_new'].unique()
if len(unmapped_samples) > 0:
    print(f"警告: 以下样本没有映射到T值: {unmapped_samples}")
else:
    print("所有样本都已成功映射到T值")

# 统计T值的分布
print("\nT值分布:")
t_counts = adata_filtered.obs['T'].value_counts()
for t_value, count in t_counts.items():
    print(f"  {t_value}: {count} 个细胞")

# 6. 保存结果
output_file = "All_filter_state_samp_new_with_CNs10-15_T.h5ad"
adata_filtered.write_h5ad(output_file)
print(f"\n结果已保存到: {output_file}")

# 7. 额外信息输出
print("\n最终数据信息:")
print(f"细胞数: {adata_filtered.shape[0]}")
print(f"基因数: {adata_filtered.shape[1]}")
print(f"样本数: {adata_filtered.obs['sample_new'].nunique()}")
print(f"T值类别: {adata_filtered.obs['T'].nunique()}")

# 保存详细的样本-T映射关系
sample_t_mapping_df = pd.DataFrame({
    'sample_new': adata_filtered.obs['sample_new'],
    'state': adata_filtered.obs['state'],
    'T': adata_filtered.obs['T']
}).drop_duplicates().sort_values('sample_new')

sample_t_mapping_df.to_csv("sample_T_mapping.csv", index=False, encoding='utf-8-sig')
print(f"样本-T映射关系已保存到: sample_T_mapping.csv")

# 打印样本-T映射表
print("\n样本-T映射表:")
print("="*30)
print(f"{'样本':<15} {'状态':<5} {'T':<5}")
print("-"*30)
for _, row in sample_t_mapping_df.iterrows():
    print(f"{row['sample_new']:<15} {row['state']:<5} {row['T']:<5}")
