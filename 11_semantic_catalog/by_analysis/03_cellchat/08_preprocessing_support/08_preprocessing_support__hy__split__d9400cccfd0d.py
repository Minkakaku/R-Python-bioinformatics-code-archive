import scanpy as sc
import numpy as np
from scipy.spatial import cKDTree

# ============================
# 1. 读取数据
# ============================
input_file = "merged_annotated_final_group_reamin_new_with_CNs.h5ad"
adata = sc.read_h5ad(input_file)
print(f"原始数据细胞数: {adata.n_obs}")

# ---- 检查必要的列是否存在 ----
required_cols = ['singler_best', 'group_new']
missing = [col for col in required_cols if col not in adata.obs.columns]
if missing:
    print("错误：以下必要列在 obs 中缺失：", missing)
    print("当前 obs 列名：", adata.obs.columns.tolist())
    raise ValueError(f"缺失列: {missing}")

if 'spatial' not in adata.obsm:
    raise ValueError("数据中缺少 'spatial' 列，请检查 obsm")

coords = adata.obsm['spatial']

# ============================
# 2. 挑选中心细胞（仅 myCAF）
# ============================
center_mask = adata.obs['singler_best'] == 'myCAF'
center_indices = np.where(center_mask)[0]
print(f"中心细胞（myCAF）数量: {len(center_indices)}")

if len(center_indices) == 0:
    raise ValueError("未找到任何 singler_best 为 'myCAF' 的细胞")

# ============================
# 3. 半径 100 μm 搜索，并构建 base_adata（关键修改）
# ============================
tree = cKDTree(coords)
radius = 50.0
all_neighbor_indices = set()
for idx in center_indices:
    all_neighbor_indices.update(tree.query_ball_point(coords[idx], radius))

print(f"半径 50 μm 内去重后细胞数: {len(all_neighbor_indices)}")

neighbor_adata = adata[list(all_neighbor_indices)].copy()
her2_low_mask = neighbor_adata.obs['group_new'].str.contains('Her2_low', na=False)
base_adata = neighbor_adata[her2_low_mask].copy()
print(f"邻域内且含 Her2_low 的细胞数: {base_adata.n_obs}")

# =============== 关键修改开始 ===============
# 提取所有中心细胞（myCAF），不添加任何特殊标签
center_adata = adata[center_indices].copy()
# 确保 cellchat_group 列存在（初始为空），避免后续函数报错
if 'cellchat_group' not in center_adata.obs.columns:
    center_adata.obs['cellchat_group'] = ''

# 将中心细胞强制加入 base_adata（不管它们是否含 Her2_low）
# 这样它们就能参与后续所有的模式匹配筛选
base_adata = sc.concat([base_adata, center_adata])
base_adata = base_adata[~base_adata.obs_names.duplicated()].copy()
print(f"强制加入所有 myCAF 后，base_adata 总细胞数: {base_adata.n_obs}")
# =============== 关键修改结束 ===============

if base_adata.n_obs == 0:
    raise ValueError("没有找到任何满足条件的细胞，请检查数据")

# ============================
# 4. 定义通用函数（用于操作1-3）
# ============================
def create_cellchat_group(adata_sub, prefix_patterns, label_map):
    """
    通用函数：在 adata_sub 中筛选含 "_P_" 的细胞，
    然后根据 prefix_patterns 匹配标记 cellchat_group，
    最后保留标记后的细胞。
    """
    has_p = adata_sub.obs['group_new'].str.contains('_P_', na=False)
    adata_p = adata_sub[has_p].copy()
    adata_p.obs['cellchat_group'] = ''
    for pattern, label in label_map.items():
        mask = adata_p.obs['group_new'].str.contains(pattern, na=False)
        adata_p.obs.loc[mask, 'cellchat_group'] = label
    valid_mask = adata_p.obs['cellchat_group'] != ''
    return adata_p[valid_mask].copy()

# ============================
# 5. 操作1-3（基于 _P_ 和 NpCR/pCR）
# ============================
# 操作1: 全局 NpCR/pCR
label_map1 = {'_NpCR_': 'NpCR', '_pCR_': 'pCR'}
res1 = create_cellchat_group(base_adata, label_map1.keys(), label_map1)
res1.write('pre_all_cellchat.h5ad')
print(f"操作1 输出细胞数: {res1.n_obs}, 已保存至 pre_all_cellchat.h5ad")

# 操作2: K 前缀
label_map2 = {'K_NpCR': 'K_NpCR', 'K_pCR': 'K_pCR'}
res2 = create_cellchat_group(base_adata, label_map2.keys(), label_map2)
res2.write('pre_all_K_cellchat.h5ad')
print(f"操作2 输出细胞数: {res2.n_obs}, 已保存至 pre_all_K_cellchat.h5ad")

# 操作3: C 前缀
label_map3 = {'C_NpCR': 'C_NpCR', 'C_pCR': 'C_pCR'}
res3 = create_cellchat_group(base_adata, label_map3.keys(), label_map3)
res3.write('pre_all_C_cellchat.h5ad')
print(f"操作3 输出细胞数: {res3.n_obs}, 已保存至 pre_all_C_cellchat.h5ad")

# ============================
# 6. 操作4-6（基于 _NpCR_ 和 P/A）
# ============================
# 操作4: 全局 NpCR → pre/post
mask_npcr = base_adata.obs['group_new'].str.contains('_NpCR_', na=False)
adata_npcr = base_adata[mask_npcr].copy()
adata_npcr.obs['cellchat_group'] = ''
# 标记 pre（含 _P_）和 post（含 _A_）
mask_p = adata_npcr.obs['group_new'].str.contains('_P_', na=False)
mask_a = adata_npcr.obs['group_new'].str.contains('_A_', na=False)
adata_npcr.obs.loc[mask_p, 'cellchat_group'] = 'pre'
adata_npcr.obs.loc[mask_a, 'cellchat_group'] = 'post'
# 只保留 pre 或 post
adata_npcr = adata_npcr[adata_npcr.obs['cellchat_group'].isin(['pre', 'post'])].copy()
adata_npcr.write('all_npcr_cellchat.h5ad')
print(f"操作4 输出细胞数: {adata_npcr.n_obs}, 已保存至 all_npcr_cellchat.h5ad")

# 操作5: K_NpCR → K_pre/K_post
mask_k_npcr = (base_adata.obs['group_new'].str.contains('K_', na=False) &
               base_adata.obs['group_new'].str.contains('_NpCR_', na=False))
adata_k_npcr = base_adata[mask_k_npcr].copy()
adata_k_npcr.obs['cellchat_group'] = ''
mask_p = adata_k_npcr.obs['group_new'].str.contains('_P_', na=False)
mask_a = adata_k_npcr.obs['group_new'].str.contains('_A_', na=False)
adata_k_npcr.obs.loc[mask_p, 'cellchat_group'] = 'K_pre'
adata_k_npcr.obs.loc[mask_a, 'cellchat_group'] = 'K_post'
adata_k_npcr = adata_k_npcr[adata_k_npcr.obs['cellchat_group'].isin(['K_pre', 'K_post'])].copy()
adata_k_npcr.write('all_npcr_K_cellchat.h5ad')
print(f"操作5 输出细胞数: {adata_k_npcr.n_obs}, 已保存至 all_npcr_K_cellchat.h5ad")

# 操作6: C_NpCR → C_pre/C_post
mask_c_npcr = (base_adata.obs['group_new'].str.contains('C_', na=False) &
               base_adata.obs['group_new'].str.contains('_NpCR_', na=False))
adata_c_npcr = base_adata[mask_c_npcr].copy()
adata_c_npcr.obs['cellchat_group'] = ''
mask_p = adata_c_npcr.obs['group_new'].str.contains('_P_', na=False)
mask_a = adata_c_npcr.obs['group_new'].str.contains('_A_', na=False)
adata_c_npcr.obs.loc[mask_p, 'cellchat_group'] = 'C_pre'
adata_c_npcr.obs.loc[mask_a, 'cellchat_group'] = 'C_post'
adata_c_npcr = adata_c_npcr[adata_c_npcr.obs['cellchat_group'].isin(['C_pre', 'C_post'])].copy()
adata_c_npcr.write('all_npcr_C_cellchat.h5ad')
print(f"操作6 输出细胞数: {adata_c_npcr.n_obs}, 已保存至 all_npcr_C_cellchat.h5ad")

print("全部6个操作完成！")
