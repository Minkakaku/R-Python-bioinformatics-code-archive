# 1. 环境配置
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
adata = sc.read_h5ad("adata_celltype.h5ad")
# 中文字体设置（需提前安装）
plt.rcParams["font.sans-serif"] = ["SimHei"]
plt.rcParams["axes.unicode_minus"] = False

# 2. 数据预处理
# 假设原始数据已加载为adata
if not adata.obs_names.is_unique:
    adata.obs_names_make_unique(join="_sample_")  # 自定义分隔符

# 3. 可视化绘图
pca_df = pd.DataFrame(
    adata.obsm["X_pca_harmony"][:, :2],
    columns=["PC1", "PC2"],
    index=adata.obs.index
)
pca_df["sample"] = adata.obs["sample"]
print(pca_df)

plt.figure(figsize=(10, 6))
ax = sns.scatterplot(
    data=pca_df,
    x="PC1",
    y="PC2",
    hue="sample",
    palette="tab20",
    s=20,
    alpha=0.7
)

# 添加带中文的标题和标签
ax.set_title("样本在PC1-PC2空间的分布", fontsize=14)
ax.set_xlabel(f"PC1方差贡献率({adata.uns['pca']['variance_ratio'][0]*100:.1f}%)")
ax.set_ylabel(f"PC2方差贡献率({adata.uns['pca']['variance_ratio'][1]*100:.1f}%)")

# 优化图例显示
plt.legend(
    title="样本来源",
    labels=adata.obs["sample"].unique(),
    bbox_to_anchor=(1.05, 1),
    loc="upper left"
)

# 4. 保存图像（避免中文乱码）
plt.savefig(
    "pca_plot.png",
    dpi=300,
    bbox_inches="tight",
    facecolor="white"
)
