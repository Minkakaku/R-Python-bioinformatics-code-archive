import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")  # 可选，忽略无关警告

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
    # low_memory=False 避免列类型混合警告
    df = pd.read_csv(INPUT_FILE, low_memory=False)

    # 将通路评分列强制转为数值，无法解析的变为 NaN
    for col in PATHWAY_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # 筛选目标细胞类型
    df = df[df["celltype"].isin(TARGET_CELLTYPES)].copy()

    # 对每种细胞类型绘图
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
            # 如果整个矩阵全是 NaN（该细胞类型的该分组无数据），跳过
            if matrix.isnull().all().all():
                print(f"跳过 {celltype} {stat_name}：所有分组均无有效数值")
                continue

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

    print(f"热图已全部生成，保存至 {OUTPUT_DIR} 文件夹。")

if __name__ == "__main__":
    main()
