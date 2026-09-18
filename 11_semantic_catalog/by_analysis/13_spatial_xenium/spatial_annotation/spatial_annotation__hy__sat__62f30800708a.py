import pandas as pd

# 读取Excel文件
df = pd.read_csv('CNs_10_all_cell_features.csv',sep = "\t")

# 检查必要的列是否存在
if 'CNs_10' not in df.columns or 'annotation' not in df.columns:
    raise ValueError("数据文件中必须包含'CNs_10'和'annotation'列")

# 统计每个CNs_10中细胞类型的数量
annotation_counts = df.groupby(['CNs_10', 'annotation']).size().reset_index(name='cell_count')

# 计算每个CNs_10的总细胞数
total_cells_per_cns = df.groupby('CNs_10').size().reset_index(name='total_cells')

# 合并数据并计算占比
result = pd.merge(annotation_counts, total_cells_per_cns, on='CNs_10')
result['percentage'] = (result['cell_count'] / result['total_cells']) * 100

# 按CNs_10和占比排序
result = result.sort_values(['CNs_10', 'percentage'], ascending=[True, False])

# 添加排名列
result['rank'] = result.groupby('CNs_10')['percentage'].rank(method='dense', ascending=False).astype(int)

# 选择并重排列
result = result[['CNs_10', 'annotation', 'cell_count', 'total_cells', 'percentage', 'rank']]

# 保存结果到CSV文件
result.to_csv('stat.csv', index=False, float_format='%.2f')

print("统计完成！结果已保存到stat.csv")
print(f"共处理了 {len(result['CNs_10'].unique())} 个CNs_10类别")
print(f"共包含 {len(result['annotation'].unique())} 种细胞类型")
print(f"总细胞数: {len(df)}")

# 显示每个CNs_10类别的前3种细胞类型
print("\n各CNs_10类别中占比前三的细胞类型:")
for cns in result['CNs_10'].unique():
    top3 = result[result['CNs_10'] == cns].head(3)
    print(f"\nCNs_10 {cns}:")
    for _, row in top3.iterrows():
        print(f"  {row['annotation']}: {row['percentage']:.2f}% (排名: {row['rank']})")
