import scanpy as sc
import pandas as pd
import numpy as np

# 读取h5ad文件
print("正在读取h5ad文件...")
adata = sc.read_h5ad("adata_celltype.h5ad")

# 检查所需的列是否存在
required_cols = ['cell_id', 'sample', 'celltype']
missing_cols = [col for col in required_cols if col not in adata.obs.columns]

if missing_cols:
    print(f"警告: 以下列不存在: {missing_cols}")
    print(f"可用的列: {adata.obs.columns.tolist()}")
    
    # 尝试寻找类似的列
    for missing_col in missing_cols:
        # 寻找相似的列名
        similar_cols = [col for col in adata.obs.columns if missing_col in col.lower()]
        if similar_cols:
            print(f"  '{missing_col}' 可能的替代列: {similar_cols}")
    
    # 询问用户是否继续
    response = input("是否继续提取其他列？(y/n): ")
    if response.lower() != 'y':
        exit()
else:
    print("所有必需的列都存在")

# 提取指定的列
print("正在提取指定列...")
# 使用reindex确保列的顺序
meta_df = adata.obs[required_cols].copy()

# 查看提取的数据
print(f"提取的数据形状: {meta_df.shape}")
print(f"\n前5行数据:")
print(meta_df.head())
print(f"\n列信息:")
print(meta_df.dtypes)
print(f"\n唯一值数量:")
for col in meta_df.columns:
    print(f"  {col}: {meta_df[col].nunique()}")

# 保存为meta.file
output_file = "meta.file"
meta_df.to_csv(output_file, sep=',', index=True)  # 保存为制表符分隔文件
print(f"\n元数据已保存到: {output_file}")

# 另外保存为excel格式（可选）
try:
    excel_file = "meta.xlsx"
    meta_df.to_excel(excel_file, index=True)
    print(f"元数据（Excel格式）已保存到: {excel_file}")
except ImportError:
    print("未安装openpyxl，跳过Excel格式保存")

# 保存为CSV格式（可选）
csv_file = "meta.csv"
meta_df.to_csv(csv_file, index=True)
print(f"元数据（CSV格式）已保存到: {csv_file}")

# 统计信息
print("\n=== 统计信息 ===")
for col in ['sample', 'celltype']:
    if col in meta_df.columns:
        print(f"\n{col}分布:")
        counts = meta_df[col].value_counts()
        for i, (value, count) in enumerate(counts.items()):
            if i < 10:  # 只显示前10个
                print(f"  {value}: {count} ({count/len(meta_df):.2%})")
        if len(counts) > 10:
            print(f"  ... 共 {len(counts)} 个类别")
