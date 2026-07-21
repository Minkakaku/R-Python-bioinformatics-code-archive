import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 构建数据
data = [
    ("B_follicular", "MPC", 2),
    ("B_germinal_center", "MPC", 2),
    ("B_germinal_center", "TNBC", 9),
    ("B_memory", "MPC", 2),
    ("B_memory", "TNBC", 7),
    ("CD4+_Tcm", "MPC", 2),
    ("CD4+_Tcm", "TNBC", 9),
    ("CD4+_Tex", "MPC", 1),
    ("CD4+_Tex", "TNBC", 24),
    ("CD8+_Tcm", "MPC", 1),
    ("CD8+_Tcm", "TNBC", 7),
    ("CD8+_Tem", "TNBC", 22),        # MPC 不存在，计数为 0
    ("CD8+_Tex", "MPC", 4),
    ("CD8+_Tex", "TNBC", 33),
    ("EMT_M2", "MPC", 34),
    ("EMT_M2", "TNBC", 264),
    ("M2", "MPC", 10),
    ("M2", "TNBC", 88),
    ("MKI67+_Plasma", "MPC", 5),
    ("MKI67+_Plasma", "TNBC", 136),
    ("MKI67+_T", "MPC", 4),
    ("MKI67+_T", "TNBC", 62),
    ("NK", "MPC", 2),
    ("NK", "TNBC", 6),
    ("Neutrophil", "MPC", 1),
    ("Neutrophil", "TNBC", 15),
    ("Plasma", "MPC", 6),
    ("Plasma", "TNBC", 203),
    ("TAM", "MPC", 21),
    ("TAM", "TNBC", 227),
    ("Treg", "MPC", 1),
    ("Treg", "TNBC", 18),
    ("cDC", "MPC", 4),
    ("cDC", "TNBC", 11),
    ("pDC", "MPC", 1),
    ("pDC", "TNBC", 3),
]

df = pd.DataFrame(data, columns=["group", "condition", "count"])

# 补全缺失的组合：如 CD8+_Tem 在 MPC 中为 0
# 通过 pivot 补全，再用 melt 还原
full = df.pivot(index="group", columns="condition", values="count").fillna(0).reset_index()
df_full = full.melt(id_vars="group", var_name="condition", value_name="count")

# 按 TNBC 数量降序排列细胞类型
order = df_full[df_full["condition"] == "TNBC"].sort_values("count", ascending=False)["group"]

# 绘图风格
sns.set_style("whitegrid")
plt.figure(figsize=(14, 6))

ax = sns.barplot(data=df_full, x="group", y="count", hue="condition",
                 order=order, palette=["#1f77b4", "#d62728"], edgecolor="black")

ax.set_xlabel("Cell Type", fontsize=13)
ax.set_ylabel("Number of High-Expression Genes", fontsize=13)
ax.set_title("Differential Expression Gene Counts: TNBC vs MPC", fontsize=15)
plt.xticks(rotation=45, ha="right", fontsize=10)
plt.legend(title="Condition")
plt.tight_layout()
plt.savefig("gene_counts_per_celltype.png", dpi=300)
plt.show()
