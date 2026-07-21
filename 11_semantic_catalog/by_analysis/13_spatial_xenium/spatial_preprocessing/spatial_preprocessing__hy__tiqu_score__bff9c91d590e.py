import scanpy as sc
import pandas as pd

# 读取数据
adata = sc.read_h5ad("All_filter_score.h5ad")

# 需要提取的列
columns_needed = [
    "cell_id", "sample", "celltype", "group_final",
    "Effector_CD8T_5k_score",
    "Proliferating_CD8T_5k_score",
    "TerminalExhausted_CD8T_5k_score",
    "IFN_Response_5k_score",
    "Hypoxia_Stress_5k_score"
]

# 检查列是否完整
missing_cols = [col for col in columns_needed if col not in adata.obs.columns]
if missing_cols:
    print(f"警告：以下列在 .obs 中不存在，将被跳过：{missing_cols}")
    columns_needed = [col for col in columns_needed if col in adata.obs.columns]

# 提取并保存
df = adata.obs[columns_needed]
df.to_csv("extracted_scores.csv", index=False)
print(f"已保存 {df.shape[0]} 行 × {df.shape[1]} 列 到 extracted_scores.csv")
