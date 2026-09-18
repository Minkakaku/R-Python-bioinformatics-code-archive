import pandas as pd
import scanpy as sc

# ------------------------------------------------------------
# 1. 读取 CSV 文件（假设无表头，制表符分隔，第一列 cell_id，第二列 sample）
# ------------------------------------------------------------
csv_path = "xenium_T.csv"
target = pd.read_csv(csv_path, header=None, sep='\t', names=['cell_id', 'sample'])

# 如果你的 CSV 是用空格分隔的，可以把 sep='\t' 改成 sep=r'\s+'

# ------------------------------------------------------------
# 2. 读取 h5ad 文件
# ------------------------------------------------------------
adata = sc.read_h5ad("Xenium_with_XY.h5ad")

# ------------------------------------------------------------
# 3. 确认 obs 中确实包含 'cell_id' 和 'sample' 列
# ------------------------------------------------------------
assert 'cell_id' in adata.obs.columns, "obs 中缺少 'cell_id' 列"
assert 'sample' in adata.obs.columns, "obs 中缺少 'sample' 列"

# ------------------------------------------------------------
# 4. 构建匹配键 → 只保留同时匹配两列的细胞
# ------------------------------------------------------------
# 为目标表创建唯一的复合键
target['match_key'] = target['cell_id'].astype(str) + '|||' + target['sample'].astype(str)
target_keys = set(target['match_key'])

# 为 AnnData 的 obs 创建同样的复合键
adata.obs['match_key'] = adata.obs['cell_id'].astype(str) + '|||' + adata.obs['sample'].astype(str)

# 布尔筛选
mask = adata.obs['match_key'].isin(target_keys)
adata_filtered = adata[mask].copy()

# 删除临时辅助列，保持 obs 整洁
del adata_filtered.obs['match_key']

# 打印匹配情况
print(f"目标细胞数：{len(target)}")
print(f"实际匹配到细胞数：{adata_filtered.n_obs}")

# ------------------------------------------------------------
# 5. 保存筛选后的 AnnData
# ------------------------------------------------------------
output_name = "XeniumT_with_XY.h5ad"
adata_filtered.write_h5ad(output_name)
print(f"已保存为：{output_name}")
