import scanpy as sc
import pandas as pd
import numpy as np

# 1. 读取h5ad数据
print("正在读取h5ad文件...")
adata = sc.read('adata_counts.h5ad')
print(f"h5ad数据形状: {adata.shape}")
print(f"h5ad的obs列: {list(adata.obs.columns)}")

# 2. 读取all.meta文件
print("\n正在读取all.meta文件...")
try:
    # 尝试不同的分隔符
    meta_df = pd.read_csv('all.meta', sep='\t')
except:
    try:
        meta_df = pd.read_csv('all.meta', sep=',')
    except:
        print("请检查all.meta文件的分隔符，默认使用制表符")

print(f"all.meta数据形状: {meta_df.shape}")
print(f"all.meta列名: {list(meta_df.columns)}")
print("\nall.meta前5行:")
print(meta_df.head())

# 检查all.meta文件的列名是否符合预期
# 可能列名有空格或其他字符，这里做标准化处理
meta_df.columns = meta_df.columns.str.strip()  # 去除列名两端的空格

# 确保列名是预期的格式
required_cols = ['cell_id', 'sample', 'cell_type']
if not all(col in meta_df.columns for col in required_cols):
    print(f"警告: all.meta文件列名不匹配，当前列名: {list(meta_df.columns)}")
    print("尝试自动重命名列...")
    
    # 尝试根据列数自动匹配
    if len(meta_df.columns) >= 3:
        meta_df.columns = ['cell_id', 'sample', 'cell_type'] + list(meta_df.columns[3:])
        print(f"重命名后的列名: {list(meta_df.columns)}")
    else:
        raise ValueError("all.meta文件需要至少3列: cell_id, sample, cell_type")

# 检查h5ad中是否有sample列
if 'sample' not in adata.obs.columns:
    print("错误: h5ad数据的obs中没有'sample'列!")
    print(f"可用的列: {list(adata.obs.columns)}")
    # 尝试查找可能的样本列
    possible_cols = [col for col in adata.obs.columns if 'sample' in col.lower() or 'group' in col.lower()]
    if possible_cols:
        print(f"可能的样本列: {possible_cols}")
        sample_col = possible_cols[0]
        # 重命名列以便后续处理
        adata.obs['sample'] = adata.obs[sample_col]
        print(f"使用列 '{sample_col}' 作为样本列")
    else:
        raise ValueError("未找到样本列，无法进行匹配")
else:
    sample_col = 'sample'

# 3. 准备匹配
print("\n准备进行匹配...")
print(f"h5ad中的细胞数: {adata.shape[0]}")
print(f"all.meta中的记录数: {meta_df.shape[0]}")

# 创建h5ad数据的匹配键
h5ad_match_keys = pd.DataFrame({
    'cell_id': adata.obs.index,
    'sample': adata.obs[sample_col]
})
h5ad_match_keys['key'] = h5ad_match_keys['cell_id'].astype(str) + '_' + h5ad_match_keys['sample'].astype(str)

# 创建meta数据的匹配键
meta_match_keys = meta_df[['cell_id', 'sample', 'cell_type']].copy()
meta_match_keys['key'] = meta_match_keys['cell_id'].astype(str) + '_' + meta_match_keys['sample'].astype(str)

# 创建从匹配键到cell_type的映射字典
key_to_celltype = dict(zip(meta_match_keys['key'], meta_match_keys['cell_type']))

# 4. 进行匹配
print("\n开始匹配...")
# 初始化celltype列为'NA'
adata.obs['celltype'] = 'NA'

# 统计匹配结果
matched_count = 0
unmatched_count = 0
matched_cell_ids = []  # 存储匹配成功的细胞ID

# 遍历h5ad中的细胞
for idx, row in h5ad_match_keys.iterrows():
    cell_key = row['key']
    
    if cell_key in key_to_celltype:
        # 找到匹配，赋值celltype
        adata.obs.loc[row['cell_id'], 'celltype'] = key_to_celltype[cell_key]
        matched_count += 1
        matched_cell_ids.append(row['cell_id'])
    else:
        # 未找到匹配，保持为'NA'
        unmatched_count += 1

# 5. 输出匹配结果统计
print("\n匹配结果统计:")
print(f"匹配成功的细胞数: {matched_count}")
print(f"未匹配的细胞数: {unmatched_count}")
print(f"匹配成功率: {matched_count/adata.shape[0]*100:.2f}%")

# 6. 创建只包含匹配成功细胞的新AnnData对象
print("\n创建只包含匹配成功细胞的新数据集...")
adata_matched = adata[matched_cell_ids, :].copy()
print(f"匹配成功数据集形状: {adata_matched.shape}")

# 7. 显示匹配后的细胞类型分布
print("\n匹配成功的细胞类型分布:")
celltype_counts = adata_matched.obs['celltype'].value_counts()
print(celltype_counts)

# 8. 保存更新后的h5ad文件
# 8.1 保存包含所有细胞但标记了celltype的原始数据
output_h5ad_all = 'adata_counts_all_with_celltype.h5ad'
adata.write(output_h5ad_all)
print(f"\n包含所有细胞的h5ad文件（已添加celltype列）已保存为: {output_h5ad_all}")

# 8.2 保存只包含匹配成功细胞的数据
output_h5ad_matched = 'adata_counts_matched_only.h5ad'
adata_matched.write(output_h5ad_matched)
print(f"只包含匹配成功细胞的h5ad文件已保存为: {output_h5ad_matched}")

# 9. 可选：将匹配结果保存为CSV文件
meta_output = pd.DataFrame({
    'cell_id': adata_matched.obs.index,
    'sample': adata_matched.obs[sample_col],
    'celltype': adata_matched.obs['celltype']
})
meta_output.to_csv('matched_celltypes_only.csv', index=False)
print(f"匹配成功的细胞信息已保存为: matched_celltypes_only.csv")

# 10. 显示前几行匹配结果
print("\n匹配成功的细胞示例 (前10行):")
print(meta_output.head(10))

# 11. 检查未匹配的细胞
if unmatched_count > 0:
    print(f"\n有 {unmatched_count} 个细胞未匹配成功，不会被保存到匹配成功的数据集中。")
    print("未匹配的细胞示例如下:")
    unmatched_cells = adata.obs[adata.obs['celltype'] == 'NA']
    if len(unmatched_cells) > 0:
        print(unmatched_cells[['sample']].head(10))
        
        # 分析未匹配的可能原因
        print("\n未匹配的可能原因分析:")
        # 1. 检查样本名称是否一致
        h5ad_samples = set(adata.obs[sample_col].unique())
        meta_samples = set(meta_df['sample'].unique())
        sample_diff = h5ad_samples - meta_samples
        if sample_diff:
            print(f"h5ad中有但all.meta中没有的样本: {sample_diff}")
        
        # 2. 检查细胞ID格式
        print(f"\nh5ad细胞ID示例: {list(adata.obs.index[:3])}")
        print(f"all.meta细胞ID示例: {list(meta_df['cell_id'].head(3))}")
        
        # 3. 检查是否有空格或其他字符问题
        print(f"\n注意: 请确保两个文件中的细胞ID格式完全一致，包括大小写、空格等")
        
        # 4. 保存未匹配细胞的信息用于调试
        unmatched_output = pd.DataFrame({
            'cell_id': adata.obs[adata.obs['celltype'] == 'NA'].index,
            'sample': adata.obs[adata.obs['celltype'] == 'NA'][sample_col]
        })
        unmatched_output.to_csv('unmatched_cells.csv', index=False)
        print(f"未匹配细胞的信息已保存为: unmatched_cells.csv")
else:
    print("\n所有细胞都匹配成功!")

print("\n脚本执行完成!")
print(f"主要输出文件:")
print(f"1. {output_h5ad_all} - 包含所有细胞（已标记celltype）")
print(f"2. {output_h5ad_matched} - 只包含匹配成功的细胞")
print(f"3. matched_celltypes_only.csv - 匹配成功细胞的元数据")
