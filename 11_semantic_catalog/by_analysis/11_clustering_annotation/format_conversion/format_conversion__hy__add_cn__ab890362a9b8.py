import scanpy as sc
import pandas as pd

# 1. 读取relation文件
relation_df = pd.read_csv('relation.tsv', sep='\t')  # 假设文件名为relation.tsv

# 2. 创建复合键映射字典
# 确保cell_id和sample列都是字符串类型
relation_df['cell_id'] = relation_df['cell_id'].astype(str)
relation_df['sample'] = relation_df['sample'].astype(str)
relation_df['composite_key'] = relation_df['cell_id'] + '_' + relation_df['sample']

# 构建映射字典
cn_dict = dict(zip(relation_df['composite_key'], relation_df['CNs_10']))

# 3. 读取h5ad数据
adata = sc.read_h5ad('All_filter_state_samp_new.h5ad')

# 4. 确保adata的cell_id和sample列为字符串
adata.obs['cell_id'] = adata.obs['cell_id'].astype(str)
adata.obs['sample'] = adata.obs['sample'].astype(str)

# 5. 创建复合键
adata.obs['composite_key'] = adata.obs['cell_id'] + '_' + adata.obs['sample']

# 6. 使用映射添加CNs_10列
adata.obs['CNs_10'] = adata.obs['composite_key'].map(cn_dict)

# 7. 删除临时复合键列
adata.obs = adata.obs.drop(columns=['composite_key'])

# 8. 检查匹配情况
matched_count = adata.obs['CNs_10'].notnull().sum()
total_count = adata.shape[0]
print(f"总细胞数: {total_count}")
print(f"成功匹配的细胞数: {matched_count}")
print(f"未匹配的细胞数: {total_count - matched_count}")

# 9. 可选：保存数据
adata.write('All_filter_state_samp_new_with_CNs10.h5ad')
