import anndata as ad
import pandas as pd
import sys

# ===================== 用户配置 =====================
H5AD_PATH = "All_filter_state_samp_new_with_CNs10.h5ad"   # 修改为实际文件路径
CELLTYPE_COL = "celltype"          # 细胞类型列名
TARGET_CELLTYPE = "matri-prod CAF" # 目标细胞类型名称（注意大小写和空格）
STATE_COL = "state"                # state列名
SAMPLE_COL = "sample"              # sample列名
OUTPUT_CSV = "matri-prod_CAF_counts.csv"  # 输出CSV文件路径，若不需要可设为None
# ===================================================

def main():
    print(f"正在读取 {H5AD_PATH} ...")
    try:
        adata = ad.read_h5ad(H5AD_PATH)
    except Exception as e:
        print(f"读取文件失败: {e}")
        sys.exit(1)

    # 检查必需的列是否存在
    for col in [CELLTYPE_COL, STATE_COL, SAMPLE_COL]:
        if col not in adata.obs.columns:
            print(f"错误: obs中缺少列 '{col}'")
            sys.exit(1)

    # 筛选目标细胞类型
    mask = adata.obs[CELLTYPE_COL] == TARGET_CELLTYPE
    n_total = mask.sum()
    print(f"总细胞数: {len(adata.obs)}")
    print(f"目标细胞类型 '{TARGET_CELLTYPE}' 的数量: {n_total}")

    if n_total == 0:
        print("警告: 未找到任何目标细胞，请检查细胞类型名称是否准确。")
        # 可打印所有celltype取值供参考
        print("可用的细胞类型:", adata.obs[CELLTYPE_COL].unique())
        sys.exit(0)

    # 筛选后的obs
    obs_filtered = adata.obs.loc[mask]

    # ---- 按 state 分组统计 ----
    print("\n=== 按 state 分组统计 ===")
    state_counts = obs_filtered[STATE_COL].value_counts()
    print(state_counts)

    # ---- 按 sample 分组统计 ----
    print("\n=== 按 sample 分组统计 ===")
    sample_counts = obs_filtered[SAMPLE_COL].value_counts()
    print(sample_counts)

    # ---- 可选：保存为CSV ----
    if OUTPUT_CSV:
        # 构造一个DataFrame，包含两种分组结果（方便保存）
        # 先转为DataFrame并重置索引
        df_state = state_counts.reset_index()
        df_state.columns = [STATE_COL, "count"]
        df_state["group_type"] = "state"

        df_sample = sample_counts.reset_index()
        df_sample.columns = [SAMPLE_COL, "count"]
        df_sample["group_type"] = "sample"

        # 合并（由于列名不同，分开保存更清晰）
        with pd.ExcelWriter(OUTPUT_CSV.replace(".csv", ".xlsx")) as writer:  # 若需Excel
            df_state.to_excel(writer, sheet_name="by_state", index=False)
            df_sample.to_excel(writer, sheet_name="by_sample", index=False)
        print(f"\n结果已保存至 {OUTPUT_CSV.replace('.csv','.xlsx')} (Excel格式)")

        # 也可保存为两个CSV
        df_state.to_csv("state_counts.csv", index=False)
        df_sample.to_csv("sample_counts.csv", index=False)
        print("同时保存了 state_counts.csv 和 sample_counts.csv")

if __name__ == "__main__":
    main()
