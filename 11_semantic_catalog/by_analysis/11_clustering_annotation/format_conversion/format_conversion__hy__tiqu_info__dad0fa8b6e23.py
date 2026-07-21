import scanpy as sc
import pandas as pd
import os

# 1. 设置文件路径
input_file = "All_filter_state_samp_new_with_CNs10-15.h5ad"
output_file = "meta.aaa"

# 2. 读取 h5ad 文件
print(f"正在读取 {input_file} ...")
adata = sc.read_h5ad(input_file)

# 3. 指定需要提取的列
columns_to_extract = ['cell_id', 'sample', 'celltype', 'state', 'sample_new']

# 4. 检查这些列是否存在于 obs 中
missing_cols = [col for col in columns_to_extract if col not in adata.obs.columns]
if missing_cols:
    print(f"警告: 以下列在 obs 中不存在: {missing_cols}")
    # 仅提取存在的列
    columns_to_extract = [col for col in columns_to_extract if col in adata.obs.columns]

if not columns_to_extract:
    raise ValueError("没有找到任何有效的列，无法保存。")

# 5. 提取数据并保存为 TSV 格式 (后缀为 .aaa)
# sep='\t' 表示制表符分隔，如果你想用逗号分隔，改为 sep=','
print(f"正在提取列 {columns_to_extract} 并保存到 {output_file} ...")
meta_df = adata.obs[columns_to_extract]
meta_df.to_csv(output_file, sep='\t', index=False)

print("提取完成！")
