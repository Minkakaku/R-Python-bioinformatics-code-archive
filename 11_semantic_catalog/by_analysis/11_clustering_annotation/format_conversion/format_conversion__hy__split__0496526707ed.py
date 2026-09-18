import scanpy as sc

# 输入文件名
input_file = "All_filter_state_samp_new_with_CNs10-15.h5ad"

# 读取数据
print(f"正在读取: {input_file}")
adata = sc.read_h5ad(input_file)

# 检查 state 列是否存在
if 'state' not in adata.obs.columns:
    raise KeyError("obs 中不存在 'state' 列，请检查数据。")

# 需要拆分的类别
categories = ['N', 'A', 'C']

for cat in categories:
    # 提取子集
    subset = adata[adata.obs['state'] == cat].copy()
    n_cells = subset.n_obs
    print(f"state = {cat}: {n_cells} 个细胞")

    if n_cells == 0:
        print(f"警告: 没有 state = {cat} 的细胞，跳过保存。")
        continue

    # 生成输出文件名
    base_name = input_file.replace(".h5ad", "")
    output_file = f"{base_name}_{cat}.h5ad"

    # 保存
    subset.write(output_file)
    print(f"已保存至: {output_file}")

print("拆分完成！")
