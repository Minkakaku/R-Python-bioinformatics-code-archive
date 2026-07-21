import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

# ============================
# 配置区域
# ============================
INPUT_FILE = "all_with_group.csv"
OUTPUT_DIR = "score_heatmap"
GROUP_ORDER = ["N", "A", "C"]
CMAP = "coolwarm"

TARGET_CELLTYPES = [
    "CXCL5+ Inf stem tumor Epi",
    "IGF2+ prolif stem tumor Epi",
    "Hypo Stem tumor Epi",
    "WNT+ ISC",
    "ESC-like Epi",
    "Prol Metab Epi",
    "EMT-like Inv Ent",
    "Inf Hyp Epi",
    "Early stress Epi",
    "Metab stress Epi",
    "Mat Ent",
    "MUC5B+ CGC",
    "MUC5AC+ MGC",
    "PCSK1+ AGC",
    "Hypo Metab Epi"
]

PATHWAY_COLS = [
    "DNA_repair_score",
    "EMT_score",
    "angiogenesis_score",
    "cellular_metabolism_score",
    "cycling_score",
    "dedifferentiation_score",
    "hypoxia_score",
    "inflammation_score",
    "oncogenic_TFs_score",
    "stemness_score"
]

# ============================
# 主程序
# ============================
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ---------- 读取数据 ----------
    df = pd.read_csv(INPUT_FILE, low_memory=False)

    for col in PATHWAY_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    df = df[df["celltype"].isin(TARGET_CELLTYPES)].copy()

    # 用于存储所有细胞类型和统计量的长格式数据
    all_rows = []

    for celltype, sub_df in df.groupby("celltype"):
        # ---------- 计算中位数 ----------
        median_df = sub_df.groupby("group")[PATHWAY_COLS].median()
        median_df = median_df.reindex(GROUP_ORDER)
        median_df = median_df.T

        # ---------- 计算平均数 ----------
        mean_df = sub_df.groupby("group")[PATHWAY_COLS].mean()
        mean_df = mean_df.reindex(GROUP_ORDER)
        mean_df = mean_df.T

        for stat_name, matrix in [("median", median_df), ("mean", mean_df)]:
            # 跳过全空矩阵
            if matrix.isnull().all().all():
                print(f"跳过 {celltype} {stat_name}：所有分组均无有效数值")
                continue

            # ----- 原有热图绘制 -----
            fig, ax = plt.subplots(figsize=(6, 8))
            sns.heatmap(
                matrix,
                annot=True,
                fmt=".3f",
                cmap=CMAP,
                center=0,
                linewidths=0.5,
                ax=ax,
                cbar_kws={"shrink": 0.8}
            )
            ax.set_title(f"{celltype}\n{stat_name.capitalize()} Score", fontsize=12)
            ax.set_ylabel("Pathway")
            ax.set_xlabel("Group")
            plt.tight_layout()

            safe_celltype = celltype.replace("/", "_").replace(" ", "_").replace("+", "p")
            base_name = f"{safe_celltype}_{stat_name}_heatmap"
            png_path = os.path.join(OUTPUT_DIR, f"{base_name}.png")
            pdf_path = os.path.join(OUTPUT_DIR, f"{base_name}.pdf")
            plt.savefig(png_path, dpi=300, bbox_inches="tight")
            plt.savefig(pdf_path, bbox_inches="tight")
            plt.close(fig)

            # ----- 收集 CSV 数据 -----
            stacked = matrix.stack().reset_index()
            stacked.columns = ["Pathway", "Group", "Score"]
            stacked["Celltype"] = celltype
            stacked["Statistic"] = stat_name
            all_rows.append(stacked)

    # ---------- 输出汇总 CSV ----------
    if all_rows:
        final_df = pd.concat(all_rows, ignore_index=True)
        # 调整列顺序，便于查看
        final_df = final_df[["Celltype", "Statistic", "Pathway", "Group", "Score"]]
        csv_path = os.path.join(OUTPUT_DIR, "all_scores.csv")
        final_df.to_csv(csv_path, index=False)
        print(f"汇总评分表已保存至 {csv_path}")
    else:
        print("警告：没有有效数据可输出为 CSV。")

    print(f"热图已全部生成，保存至 {OUTPUT_DIR} 文件夹。")

if __name__ == "__main__":
    main()
