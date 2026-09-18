import os
os.environ['OPENBLAS_NUM_THREADS'] = '0'
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import seaborn as sns
from natsort import natsorted
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.colors import ListedColormap
from matplotlib.pyplot import rc_context
from sklearn.cluster import MiniBatchKMeans
from sklearn.neighbors import radius_neighbors_graph, kneighbors_graph

plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
sc.settings.set_figure_params(dpi=100, dpi_save=200, figsize=(5, 5), facecolor='white')


def calculate_cn(
        adata, 
        mode='knn', 
        ident="sample", 
        column_ct="celltype", 
        radius=50, 
        k_neighbor=20, 
        n_clusters=10,
        output_dir=None):
    """
    多样本 Cellular Neighborhood (CN) 分析
    """
    os.makedirs(output_dir, exist_ok=True)

    samples = adata.obs[ident].unique().tolist()
    neighbors_dict = {}
    all_celltypes = sorted(adata.obs[column_ct].unique().tolist())  # 确保所有样本对齐

    for sample in samples:
        tmp = adata[adata.obs[ident] == sample]
        x_y_coordinates = tmp.obsm["spatial"]
        cts = tmp.obs[column_ct]

        # 构建邻居图
        if mode == 'radius':
            connectivity_matrix = radius_neighbors_graph(
                x_y_coordinates, radius, mode='connectivity', include_self=False
            )
        elif mode == 'knn':
            connectivity_matrix = kneighbors_graph(
                x_y_coordinates, k_neighbor, mode='connectivity', include_self=False
            )
        else:
            raise ValueError("Mode should be either 'radius' or 'knn'.")

        # one-hot，确保与全局 celltypes 对齐
        one_hot_matrix = pd.get_dummies(cts).reindex(columns=all_celltypes, fill_value=0).values
        one_hot_matrix = pd.get_dummies(cts).reindex(columns=all_celltypes, fill_value=0).values.astype(np.float32)
        cell_by_neighborhood_matrix = connectivity_matrix.dot(one_hot_matrix)

        # 归一化
        row_sums = cell_by_neighborhood_matrix.sum(axis=1, keepdims=True)
        norm_matrix = np.divide(
            cell_by_neighborhood_matrix, row_sums, where=row_sums != 0
        )

        # 如果某个 cell 没有邻居 → 用自己的 cell type one-hot 代替
        zero_idx = np.where(row_sums.flatten() == 0)[0]
        if len(zero_idx) > 0:
            norm_matrix[zero_idx] = one_hot_matrix[zero_idx]

        # 构建 DataFrame
        df = pd.DataFrame(norm_matrix, columns=all_celltypes, index=tmp.obs_names)
        df["celltype"] = cts.values
        df[["x", "y"]] = x_y_coordinates
        df[ident] = sample
        neighbors_dict[sample] = df

    # 合并所有样本
    all_neighbors_df = pd.concat(neighbors_dict.values())

    # 提取特征列（所有 cell type 列）
    feature_cols = all_celltypes
    values = all_neighbors_df[feature_cols].values

    # KMeans 聚类
    km = MiniBatchKMeans(n_clusters=n_clusters, random_state=0)
    labels = km.fit_predict(values)
    adata.obs[f'CNs_{n_clusters}'] = [f'CN{i+1}' for i in labels]
    cn_labels = [f"CN{i+1}" for i in labels]
    all_neighbors_df["CN"] = cn_labels
    k_centroids = km.cluster_centers_

    # fold-change 分析 (log2 ratio)
    tissue_avgs = values.mean(axis = 0)
    fc = np.log2(((k_centroids+tissue_avgs)/(k_centroids+tissue_avgs).sum(axis = 1, keepdims = True))/tissue_avgs)
    fc = pd.DataFrame(fc, columns = all_celltypes, index = ["CN"+ str(i+1) for i in range(n_clusters)])

    # CN热图
    plt.figure(figsize=(6, 4))
    sns.set_style("white")
    sns.clustermap(fc, vmin=-2, vmax=2, cmap="vlag", row_cluster=False, linewidths=0.5)
    plt.savefig(os.path.join(output_dir, f"CNs_{n_clusters}.png"), bbox_inches='tight')

    # CN 空间原位映射
    if len(adata.obs[f'CNs_{n_clusters}'].unique()) <= 10:
        palette = dict(zip(natsorted(adata.obs[f'CNs_{n_clusters}'].drop_duplicates().tolist()), list(map(colors.to_hex, cm.tab10.colors))))
    else:
        palette = dict(zip(natsorted(adata.obs[f'CNs_{n_clusters}'].drop_duplicates().tolist()), list(map(colors.to_hex, cm.tab20.colors))))
        
    library_names = adata.obs[ident].drop_duplicates().to_list()
    for library in library_names:
        adata_sub = adata[adata.obs[ident] == library, :]
        cn_sort_list = natsorted(adata_sub.obs[f'CNs_{n_clusters}'].drop_duplicates().tolist())
        current_palette  = ListedColormap([palette[i] for i in cn_sort_list])
        with rc_context({'figure.figsize': (8, 8)}):
            sq.pl.spatial_scatter(
                adata_sub,
                library_id="spatial",
                color=[f'CNs_{n_clusters}'],
                title=f'{library} CNs',
                img=False,
                shape=None,
                palette=current_palette
                # size=1,
                # use_raw=True
            )                       
            plt.savefig(os.path.join(output_dir, f'{library}_spatial_CNs_{n_clusters}.png'), facecolor='white', bbox_inches="tight", dpi=300)   

    # 分组柱状图
    # 计算每个样本中每个聚类的细胞数量
    cluster_counts = adata.obs.groupby(['sample', f'CNs_{n_clusters}']).size().unstack(fill_value=0)
    # 计算每个样本中每个聚类的细胞占比
    cluster_proportions = cluster_counts.div(cluster_counts.sum(axis=1), axis=0)
    # 转换为长格式以便绘图
    cluster_proportions_long = cluster_proportions.reset_index().melt(id_vars='sample', var_name='Cluster', value_name='Proportion')
    # 添加组信息
    cluster_proportions_long = cluster_proportions_long.merge(adata.obs[['sample', 'group']].drop_duplicates(), on='sample')
    cluster_proportions_long.to_csv(os.path.join(output_dir, f"CNs_{n_clusters}_cell_proportion.txt"), sep='\t', index=None)

    # 提取 Cluster 列中的数字部分并排序
    cluster_order = sorted(cluster_proportions_long['Cluster'].unique(), key=lambda x: int(x.replace('CN', '')))
    # 绘制箱线图
    plt.figure(figsize=(12, 6))
    sns.boxplot(
        x='Cluster',
        y='Proportion',
        hue='group',
        data=cluster_proportions_long,
        order=cluster_order
    )
    sns.stripplot(
        x='Cluster',
        y='Proportion',
        hue='group',
        data=cluster_proportions_long,
        dodge=True,
        jitter=False,
        color='gray',
        alpha=0.3,
        legend=False,
        order=cluster_order
    )
    plt.savefig(os.path.join(output_dir, f"CNs_{n_clusters}_cell_proportion.png"), dpi=300, bbox_inches='tight')

    return adata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cellular Neighborhood (CN) analysis.")
    parser.add_argument("--adata_h5ad", type=str, required=True, help="Input AnnData file (.h5ad)")
    parser.add_argument("--column_ct", type=str, default="celltype", help="Cell type column name")
    parser.add_argument("--ident", type=str, default="sample", help="Sample identity column")
    parser.add_argument("--mode", type=str, default="knn", choices=["knn", "radius"], help="Neighbor mode")
    parser.add_argument("--radius", type=float, default=50, help="Radius for radius-neighbors")
    parser.add_argument("--n_neighbors", type=int, default=25, help="Neighbor numbers for kNN")
    parser.add_argument("--n_clusters", type=int, nargs="+", default=[10], help="Number(s) of CN clusters")
    parser.add_argument("--output_dir", type=str, default="CN_results", help="Output directory")

    args = parser.parse_args()

    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)

    # 读取数据
    adata = sc.read_h5ad(args.adata_h5ad)

    # 多参数循环
    for n_clusters in args.n_clusters:
        print(f"Running CN analysis with k={args.n_neighbors}, n_clusters={n_clusters} ...")
        adata_cn = calculate_cn(
            adata,
            mode=args.mode,
            ident=args.ident,
            column_ct=args.column_ct,
            radius=args.radius,
            k_neighbor=args.n_neighbors,
            n_clusters=n_clusters,
            output_dir=args.output_dir
        )

