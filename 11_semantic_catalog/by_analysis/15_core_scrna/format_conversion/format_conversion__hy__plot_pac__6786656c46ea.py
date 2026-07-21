import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

adata=sc.read_h5ad("adata_celltype.h5ad")
adata.obs_names_make_unique("-")
# 提取PCA坐标并转换为DataFrame
pca_df = pd.DataFrame(
    adata.obsm["X_pca"][:, :2], 
    columns=["PC1", "PC2"]
)
# 合并样本分组信息
pca_df["sample"] = adata.obs["sample"].astype("category")

plt.figure(figsize=(10, 6))
sns.scatterplot(
    data=pca_df,
    x="PC1",
    y="PC2",
    hue="sample",
    palette="tab20",  # 使用区分度高的调色板
    s=20,             # 点大小
    alpha=0.7,        # 透明度
    edgecolor="none"  # 去除边缘线
)

# 添加标题和坐标轴标签（含方差解释率）
pca_var = adata.uns["pca"]["variance_ratio"]
plt.xlabel(f"PC1 ({pca_var[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca_var[1]*100:.1f}%)")
plt.title("样本在PC1-PC2空间的分布", fontsize=14)

# 调整图例位置
plt.legend(
    bbox_to_anchor=(1.05, 1), 
    loc="upper left",
    title="Sample"
)
plt.tight_layout()
plt.savefig("pca_1_2.png", facecolor='white', bbox_inches="tight", dpi=300)
