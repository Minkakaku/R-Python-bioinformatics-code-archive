import scanpy as sc
import os

# ========== 配置 ==========
input_file = "All_filter_state_samp_new_with_CNs10-15.h5ad"  # 输入文件路径
output_dir = "./split_groups"                                # 输出文件夹，可自定义
# ==========================

# 创建输出目录
os.makedirs(output_dir, exist_ok=True)

# 读取数据
print(f"正在读取 {input_file} ...")
adata = sc.read_h5ad(input_file)

# 检查 state 列是否存在
if 'state' not in adata.obs.columns:
    raise KeyError("obs 中缺少 'state' 列，请检查数据。")

# 定义三组分组
groups = {
    "A_C": ["A", "C"],
    "A_N": ["A", "N"],
    "C_N": ["C", "N"]
}

# 分割并保存
for group_name, states in groups.items():
    # 筛选出对应 state 的细胞
    adata_sub = adata[adata.obs['state'].isin(states)].copy()
    print(f"分组 {group_name}: 保留 {adata_sub.n_obs} 个细胞 (原始 {adata.n_obs} 个)")

    # 构建输出路径，保留原始文件名前缀
    base = os.path.splitext(os.path.basename(input_file))[0]
    output_file = os.path.join(output_dir, f"{base}_state_{group_name}.h5ad")
    adata_sub.write(output_file)
    print(f"已保存: {output_file}")

print("所有分组保存完成！")
