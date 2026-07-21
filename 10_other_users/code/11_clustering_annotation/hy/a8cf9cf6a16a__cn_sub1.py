import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import umap  # 新增UMAP库

# 加载h5ad数据
adata = sc.read_h5ad('adata_celltype_modified.h5ad')  # 请替换为您的文件路径

# 检查数据基本结构
print("数据形状:", adata.shape)
print("细胞类型列名:", adata.obs.columns)
print("细胞类型类别:", adata.obs['celltype'].unique())

# 基于基因表达数据计算每个细胞的10个最近邻
# 使用adata.X作为特征矩阵（行是细胞，列是基因）
X = adata.X
if hasattr(X, 'toarray'):  # 如果是稀疏矩阵，转换为密集矩阵（注意内存使用）
    X = X.toarray()

# 计算最近邻，n_neighbors=11（包括自身）
nbrs = NearestNeighbors(n_neighbors=11, algorithm='auto', metric='euclidean')
nbrs.fit(X)
distances, indices = nbrs.kneighbors(X)  # indices包含每个细胞及其邻居的索引

# 获取所有细胞的celltype标签
celltypes = adata.obs['celltype']
unique_celltypes = celltypes.unique()
n_cells = adata.n_obs
n_types = len(unique_celltypes)

# 创建细胞类型到索引的映射
type_to_idx = {ctype: idx for idx, ctype in enumerate(unique_celltypes)}

# 初始化比例矩阵：每行代表一个细胞，每列代表一种细胞类型的比例
proportion_matrix = np.zeros((n_cells, n_types))

# 计算每个细胞的10个最近邻的细胞类型比例
for i in range(n_cells):
    # 取第2到第11个邻居（排除自身，索引0是自身）
    neighbor_indices = indices[i, 1:11]
    neighbor_types = celltypes.iloc[neighbor_indices]
    
    # 统计邻居中每种细胞类型的数量
    type_counts = neighbor_types.value_counts()
    
    # 计算比例
    for ctype, count in type_counts.items():
        j = type_to_idx[ctype]
        proportion_matrix[i, j] = count / 10.0

# 将比例矩阵转换为DataFrame便于查看
proportion_df = pd.DataFrame(proportion_matrix, 
                             index=adata.obs.index,
                             columns=unique_celltypes)
print("比例矩阵前5行:")
print(proportion_df.head())

# 1. 输出proportion_df为CSV文件
proportion_df.to_csv('proportion_matrix.csv', index=True, encoding='utf-8-sig')
print("比例矩阵已保存为 proportion_matrix.csv")

# 提取celltype为"Epithelials"的细胞亚群
epithelial_mask = adata.obs['celltype'] == 'Epithelials'
epithelial_indices = np.where(epithelial_mask)[0]

if len(epithelial_indices) == 0:
    # 尝试常见变体
    possible_names = ['Epithelials', '上皮', 'epithelial', 'Epithelial']
    for name in possible_names:
        epithelial_mask = adata.obs['celltype'] == name
        epithelial_indices = np.where(epithelial_mask)[0]
        if len(epithelial_indices) > 0:
            print(f"使用细胞类型标签: {name}")
            break

if len(epithelial_indices) == 0:
    raise ValueError("未找到上皮细胞亚群，请检查数据中celltype列的实际标签名称")

epithelial_proportions = proportion_matrix[epithelial_indices, :]
epithelial_cell_barcodes = adata.obs.index[epithelial_indices]

print(f"找到上皮细胞数量: {len(epithelial_indices)}")

# 使用比例矩阵进行降维
# 首先使用PCA降维
pca = PCA(n_components=min(50, epithelial_proportions.shape[1]))
X_pca = pca.fit_transform(epithelial_proportions)
print(f"PCA降维后形状: {X_pca.shape}")

# 3. 生成PCA碎石图并保存[6,7](@ref)
plt.figure(figsize=(10, 6))
# 计算累计方差贡献率
explained_variance_ratio = pca.explained_variance_ratio_
cumulative_variance = np.cumsum(explained_variance_ratio)

# 绘制碎石图
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))

# 左侧：方差贡献率
ax1.plot(range(1, len(explained_variance_ratio) + 1), explained_variance_ratio, 'bo-')
ax1.set_xlabel('主成分')
ax1.set_ylabel('方差贡献率')
ax1.set_title('PCA方差贡献率')
ax1.grid(True)

# 右侧：累计方差贡献率
ax2.plot(range(1, len(cumulative_variance) + 1), cumulative_variance, 'ro-')
ax2.set_xlabel('主成分数量')
ax2.set_ylabel('累计方差贡献率')
ax2.set_title('PCA累计方差贡献率')
ax2.grid(True)

# 标记80%和90%方差线
ax2.axhline(y=0.8, color='g', linestyle='--', alpha=0.8, label='80%方差')
ax2.axhline(y=0.9, color='purple', linestyle='--', alpha=0.8, label='90%方差')
ax2.legend()

plt.tight_layout()
plt.savefig('pca_scree_plot.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"前10个主成分的累计方差贡献率: {cumulative_variance[:10]}")

# 使用t-SNE进一步降维到2D用于可视化
tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(epithelial_indices)-1))
X_tsne = tsne.fit_transform(X_pca)
print(f"t-SNE降维后形状: {X_tsne.shape}")

# 4. 使用UMAP进行降维可视化[9](@ref)
umap_reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
X_umap = umap_reducer.fit_transform(epithelial_proportions)
print(f"UMAP降维后形状: {X_umap.shape}")

# 聚类
# 使用肘部法则或先验知识确定聚类数量，这里假设3类
n_clusters = min(5, len(epithelial_indices))  # 增加最大聚类数到5
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
clusters = kmeans.fit_predict(X_pca)  # 基于PCA结果聚类

# 5. 将聚类ID命名为cluster_cn
cluster_cn = clusters

# 可视化结果
plt.figure(figsize=(18, 6))

# 子图1: t-SNE可视化着色按聚类
plt.subplot(1, 3, 1)
scatter1 = plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=cluster_cn, cmap='tab10', s=20, alpha=0.7)
plt.colorbar(scatter1, label='Cluster')
plt.xlabel('t-SNE 1')
plt.ylabel('t-SNE 2')
plt.title('t-SNE可视化 (基于邻居比例矩阵)')

# 子图2: UMAP可视化着色按聚类[9](@ref)
plt.subplot(1, 3, 2)
scatter2 = plt.scatter(X_umap[:, 0], X_umap[:, 1], c=cluster_cn, cmap='tab10', s=20, alpha=0.7)
plt.colorbar(scatter2, label='Cluster')
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.title('UMAP可视化 (基于邻居比例矩阵)')

# 子图3: 显示每个聚类中主要细胞类型比例
plt.subplot(1, 3, 3)
# 计算每个聚类中原始细胞类型的平均比例（在邻居中）
cluster_profiles = []
for cluster_id in range(n_clusters):
    cluster_mask = cluster_cn == cluster_id
    cluster_avg_proportions = epithelial_proportions[cluster_mask].mean(axis=0)
    cluster_profiles.append(cluster_avg_proportions)

cluster_profiles = np.array(cluster_profiles)
im = plt.imshow(cluster_profiles.T, aspect='auto', cmap='YlOrRd')
plt.colorbar(im, label='平均比例')
plt.xlabel('聚类ID')
plt.ylabel('细胞类型')
plt.xticks(range(n_clusters))
plt.yticks(range(len(unique_celltypes)), unique_celltypes, fontsize=8)
plt.title('各聚类中邻居细胞类型比例')

plt.tight_layout()
plt.savefig('epithelial_clustering_results.png', dpi=300, bbox_inches='tight')
plt.show()

# 输出聚类结果统计
print("\n聚类结果统计:")
for i in range(n_clusters):
    cluster_size = np.sum(cluster_cn == i)
    print(f"聚类 {i}: {cluster_size} 个细胞 ({cluster_size/len(epithelial_indices)*100:.1f}%)")

# 6. 保存数据为新的h5ad
# 创建上皮亚群的AnnData对象
epithelial_adata = sc.AnnData(X=epithelial_proportions,
                              obs=pd.DataFrame(index=epithelial_cell_barcodes))
epithelial_adata.obs['cluster_cn'] = cluster_cn  # 使用新的列名
epithelial_adata.obsm['X_pca'] = X_pca
epithelial_adata.obsm['X_tsne'] = X_tsne
epithelial_adata.obsm['X_umap'] = X_umap  # 添加UMAP坐标

# 保存PCA相关结果
epithelial_adata.uns['pca_variance_ratio'] = explained_variance_ratio
epithelial_adata.uns['pca_cumulative_variance'] = cumulative_variance
epithelial_adata.varm['PCs'] = pca.components_.T  # 主成分载荷

# 保存到新的h5ad文件
epithelial_adata.write('epithelial_subcluster_analysis.h5ad')
print("上皮细胞亚群分析结果已保存为 epithelial_subcluster_analysis.h5ad")

# 同时保存聚类结果的CSV文件
results_df = pd.DataFrame({
    'cell_barcode': epithelial_cell_barcodes,
    'cluster_cn': cluster_cn,
    'tsne_1': X_tsne[:, 0],
    'tsne_2': X_tsne[:, 1],
    'umap_1': X_umap[:, 0],
    'umap_2': X_umap[:, 1]
})
results_df.to_csv('epithelial_clustering_results.csv', index=False, encoding='utf-8-sig')
print("聚类详细结果已保存为 epithelial_clustering_results.csv")

print("分析完成！所有结果文件已保存。")
