import os
import anndata as ad

# ==================== 配置参数 ====================
H5AD_PATH = "All_filter_state_samp_new_with_CNs10-15.h5ad"   # 输入文件路径
OUTPUT_BASE_DIR = "."                                         # 输出根目录，可修改为指定路径

# 样本映射字典：文件夹名 -> 样本名列表
SAMPLE_MAPPING = {
    "output-XETG00305__0086299__TMA2__20251220__034037": [
        "SR202500005C", "SR202500005D", "SR202500009D-1", "SR202500009D-2",
        "SR2025000012F", "SR2025000012L-1", "SR2025000012L-2"
    ],
    "output-XETG00546__0066018__009D-1__20250919__053836": [
        "009D-1"
    ],
    "output-XETG00546__0066018__SR202500002C__20250919__053836": [
        "SR202500002C-1", "SR202500002C-2"
    ],
    "output-XETG00546__0066018__SR202500003B_009D-2__20250919__053835": [
        "SR202500003B", "009D-2"
    ],
    "output-XETG00546__0066018__SR202500005F-1_SR202500005F-2__20250919__053836": [
        "SR202500005F-1", "SR202500005F-2"
    ],
    "output-XETG00546__0066411__SR2025000014B__20260129__075125": [
        "SR2025000014B"
    ],
    "output-XETG00546__0066411__SR2025000014D_SR2025000014F_SR2025000016J__20260129__075125": [
        "SR2025000014D", "SR2025000014F", "SR2025000016J"
    ],
    "output-XETG00546__0066411__SR2025000016K__20260129__075125": [
        "SR2025000016K"
    ],
}

# ==================== 主处理逻辑 ====================
def main():
    print(f"正在读取 AnnData 文件: {H5AD_PATH}")
    adata = ad.read_h5ad(H5AD_PATH)
    print(f"原始数据维度: {adata.shape}")
    print(f"obs 中包含的样本: {adata.obs['sample'].unique().tolist()}")

    # 检查所有待提取样本是否存在于 obs['sample'] 中
    all_requested_samples = set()
    for samples in SAMPLE_MAPPING.values():
        all_requested_samples.update(samples)

    existing_samples = set(adata.obs['sample'].unique())
    missing_samples = all_requested_samples - existing_samples
    if missing_samples:
        print(f"警告: 以下样本在数据中不存在，将被跳过: {missing_samples}")

    # 遍历每个输出文件夹
    for folder_name, target_samples in SAMPLE_MAPPING.items():
        # 过滤出实际存在于数据中的样本
        valid_samples = [s for s in target_samples if s in existing_samples]
        if not valid_samples:
            print(f"跳过文件夹 '{folder_name}': 所有目标样本均不存在")
            continue

        # 提取子集
        subset_adata = adata[adata.obs['sample'].isin(valid_samples), :].copy()
        print(f"文件夹 '{folder_name}': 提取了 {subset_adata.n_obs} 个细胞, 来自样本 {valid_samples}")

        # 创建输出目录
        output_dir = os.path.join(OUTPUT_BASE_DIR, folder_name)
        os.makedirs(output_dir, exist_ok=True)

        # 保存为 h5ad 文件
        output_path = os.path.join(output_dir, "subset_data.h5ad")
        subset_adata.write_h5ad(output_path)
        print(f"已保存: {output_path}\n")

    print("全部处理完成！")

if __name__ == "__main__":
    main()
