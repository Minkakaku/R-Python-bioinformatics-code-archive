#!/usr/bin/env python3
"""
增强版：从 h5ad 文件中先按 obs['group_final'] 筛选 'PRE-NR' 和 'PRE-R'，
再按 obs['site_band'] 提取 'in', 'on', 'out'，分别保存为三个 h5ad，
并将三者合并为一个包含 site_band 的 h5ad。
"""

import anndata as ad
import sys

# ========== 用户参数 ==========
INPUT_FILE = "adata_with_band_by_sample.h5ad"
OUT_IN     = "adata_in.h5ad"
OUT_ON     = "adata_on.h5ad"
OUT_OUT    = "adata_out.h5ad"
OUT_MERGED = "adata_in_on_out.h5ad"
# =============================

def main():
    print(f"读取数据: {INPUT_FILE}")
    adata = ad.read_h5ad(INPUT_FILE)

    # ---------- 新增步骤：按 group_final 筛选 ----------
    REQUIRED_GROUP_COL = 'group_final'
    if REQUIRED_GROUP_COL not in adata.obs.columns:
        print(f"错误: obs 中没有 '{REQUIRED_GROUP_COL}' 列")
        sys.exit(1)

    target_groups = ['PRE-NR', 'PRE-R']
    mask_group = adata.obs[REQUIRED_GROUP_COL].isin(target_groups)
    print(f"根据 '{REQUIRED_GROUP_COL}' 筛选出 {mask_group.sum()} 个细胞 "
          f"(目标值: {target_groups})")
    if mask_group.sum() == 0:
        print("错误: 没有细胞属于 'PRE-NR' 或 'PRE-R'，终止执行。")
        sys.exit(1)

    # 复制一份过滤后的数据，避免后续操作影响原变量
    adata_filtered = adata[mask_group].copy()
    # --------------------------------------------------------

    # 检查 site_band 列
    if 'site_band' not in adata_filtered.obs.columns:
        print("错误: obs 中没有 'site_band' 列")
        sys.exit(1)

    # 提取三个 site_band 子集
    adata_in  = adata_filtered[adata_filtered.obs['site_band'] == 'in'].copy()
    adata_on  = adata_filtered[adata_filtered.obs['site_band'] == 'on'].copy()
    adata_out = adata_filtered[adata_filtered.obs['site_band'] == 'out'].copy()

    # 分别保存
    print(f"保存 {OUT_IN} (细胞数: {adata_in.n_obs})")
    adata_in.write_h5ad(OUT_IN)
    print(f"保存 {OUT_ON} (细胞数: {adata_on.n_obs})")
    adata_on.write_h5ad(OUT_ON)
    print(f"保存 {OUT_OUT} (细胞数: {adata_out.n_obs})")
    adata_out.write_h5ad(OUT_OUT)

    # 合并 in, on, out
    adata_merged = ad.concat(
        [adata_in, adata_on, adata_out],
        axis=0,
        join='outer',      # 保留所有基因（如基因集一致也可用 'inner'）
        index_unique=None  # 若细胞条形码可能重复，可改为 '-' 以防冲突
    )

    print(f"合并数据: {adata_merged.n_obs} 个细胞, {adata_merged.n_vars} 个基因")
    print(f"保存 {OUT_MERGED}")
    adata_merged.write_h5ad(OUT_MERGED)

    print("完成。")

if __name__ == "__main__":
    main()
