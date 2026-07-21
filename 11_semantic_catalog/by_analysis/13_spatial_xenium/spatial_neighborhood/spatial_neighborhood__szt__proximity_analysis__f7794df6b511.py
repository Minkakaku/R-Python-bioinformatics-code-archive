import os
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from scipy.spatial import KDTree

plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
sc.settings.set_figure_params(dpi=100, dpi_save=200, figsize=(5, 5), facecolor='white')


def calculate_proximity(adata, distance_threshold=50, celltype_column='celltype', target_celltype=None, out_dir=None):
    group_list = adata.obs['group'].unique()
    palette = dict(zip(group_list, sns.color_palette("Set2", n_colors=len(group_list))))
    
    all_types = adata.obs[celltype_column].unique()
    proximity_data = []
    for sample in adata.obs['sample'].unique():
        ad_s = adata[adata.obs['sample'] == sample].copy()
        sq.gr.spatial_neighbors(ad_s, coord_type='generic', radius=distance_threshold)
        cells_of_type = ad_s.obs_names[ad_s.obs[celltype_column] == target_celltype]
        for cell_label in cells_of_type:
            cell_idx = ad_s.obs_names.get_loc(cell_label)
            neighbors = ad_s.obsp['spatial_connectivities'][cell_idx].nonzero()[1]
            neighbor_labels = ad_s.obs_names[neighbors]
            neighbor_types = ad_s.obs.loc[neighbor_labels, celltype_column]
            # 去掉自己
            neighbor_types = neighbor_types[neighbor_types.index != cell_label]
            counts = neighbor_types.value_counts(normalize=True).reindex(all_types, fill_value=0.0)
            for ct, prop in counts.items():
                proximity_data.append({
                    'celltype': ct,
                    'proportion': prop,
                    'sample': sample,
                    'group': ad_s.obs.loc[cell_label, 'group']
                })

    results_df = pd.DataFrame(proximity_data)
    results_df.to_csv(os.path.join(out_dir, f'nearly_{distance_threshold}um_distance_celltype_proportion.xls'), sep='\t', index=None)

    plt.figure(figsize=(10, 6))
    ax = sns.boxplot(
        x='celltype', 
        y='proportion', 
        hue='group', 
        palette=palette,
        data=results_df, 
        showfliers=False
    )
    ax.grid(visible=False)

    max_proportion = results_df['proportion'].max()
    # max_proportion = 1
    # y_offset = max_proportion * 0.1

    # 进行统计检验并添加 p 值
    celltypes = results_df['celltype'].unique()
    for celltype in celltypes:
        data_subset = results_df[results_df['celltype'] == celltype]
        groups = data_subset['group'].unique()       
        if len(groups) == 2:  # 确保只有两个组
            group1 = data_subset[data_subset['group'] == groups[0]]['proportion']
            group2 = data_subset[data_subset['group'] == groups[1]]['proportion']

            # 使用 Mann-Whitney U 检验计算 p 值
            if len(group1) > 0 and len(group2) > 0:
                stat, p_value = mannwhitneyu(group1, group2)

                # 判断显著性并标注
                if p_value < 0.01:
                    annotation = '**'
                elif p_value < 0.05:
                    annotation = '*'
                else:
                    annotation = 'ns'

                # 添加注释，标注显著性水平，仅显示符号，放在图形顶部
                x_position = list(celltypes).index(celltype)
                ax.annotate(f'{annotation}', 
                            xy=(x_position, max_proportion-max_proportion * 0.05),
                            ha='center', color='black', fontsize=12)

        # ax.set_ylim(0, max_proportion * 1.2)
    
    ax.legend(title="", bbox_to_anchor=(1.02, 0.5), loc='upper left', borderaxespad=0.)
    
    plt.xlabel("")
    plt.title(f'Celltype proportion of {ct} nearly {distance_threshold}um distance')
    plt.xticks(rotation=45)
    plt.tight_layout()

    plt.savefig(os.path.join(out_dir, 'proximity_50um.png'))
    plt.savefig(os.path.join(out_dir, 'proximity_50um.pdf'))
    plt.clf()


def calculate_proximity_range(adata, target_celltype, celltype_column='celltype', distance_ranges=[0, 50, 100, 150, 200], out_dir=None):
    proximity_data = []
    for sample in np.unique(adata.obs['sample'].values):
        ad = adata[adata.obs['sample']==sample].copy()
        coords = ad.obsm['spatial']
        cell_types = np.asarray(ad.obs[celltype_column])
        groups = np.asarray(ad.obs['group'])

        # 如果样本里没有目标细胞，跳过
        cells_indices = np.where(cell_types == target_celltype)[0]
        if len(cells_indices) == 0:
                continue

        # 按样本建 KDTree
        tree = KDTree(coords)

        # 对样本中每个目标细胞统计不同距离区间内的邻居
        for idx in cells_indices:
            idx_coord = coords[idx] 
            previous_indices = set()
            # 遍历区间
            for i in range(len(distance_ranges) - 1):
                # 获取在当前最大距离范围内的邻居
                indices_within_max_range = set(tree.query_ball_point(idx_coord, distance_ranges[i + 1]))
                # 去掉自己
                indices_within_max_range.discard(idx)           
                # 当前范围内的邻居是最大范围内的邻居减去之前的邻居
                current_indices = indices_within_max_range - previous_indices            
                # 更新之前的邻居集合
                previous_indices.update(current_indices)
                # 统计当前范围内的邻居类型
                neighbor_types = cell_types[list(current_indices)]
                counts = pd.Series(neighbor_types).value_counts()
                # 将计数结果存储到 proximity_data 中，并计算均值
                for neighbor_type in counts.index:
                        proximity_data.append({
                        'distance_range': f"{distance_ranges[i]}-{distance_ranges[i + 1]}",
                        'celltype': neighbor_type,
                        'count': counts[neighbor_type],
                        'group': groups[idx]
                        })

    # 转换为 DataFrame
    proximity_df = pd.DataFrame(proximity_data)

    # 按照 distance_range 和 celltype 分组，并计算均值
    mean_proximity_df = proximity_df.groupby(['distance_range', 'celltype', 'group'], as_index=False).mean()

    # Ensure correct order of distance ranges
    mean_proximity_df['distance_range'] = pd.Categorical(
        mean_proximity_df['distance_range'], 
        categories=["0-50", "50-100", "100-150", "150-200"], 
        ordered=True
    )
    mean_proximity_df.sort_values('distance_range', inplace=True)
    mean_proximity_df.to_csv(os.path.join(out_dir, 'proximity_distance_range.txt'), sep='\t', index=None)

    # 绘制折线图，按细胞类型分面展示，分成4列
    celltypes = sorted(mean_proximity_df['celltype'].unique().tolist())
    n_cols = 4  # 每行显示4个子图
    n_rows = (len(celltypes) + n_cols - 1) // n_cols  # 计算行数

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(20, n_rows * 4), sharex=True)

    for ax, celltype in zip(axs.flatten(), celltypes):
        subset = mean_proximity_df[mean_proximity_df['celltype'] == celltype]
        for group in sorted(subset['group'].unique().tolist()):
            group_data = subset[subset['group'] == group]
            ax.plot(group_data['distance_range'], group_data['count'], marker='o', label=group)
        
        ax.set_title(f'{celltype}')
        ax.set_ylabel('Cell count')
        # ax.set_xlabel('Distance (um)')
        # ax.grid()
        ax.legend()

    # 隐藏多余的子图（如果有）
    for i in range(len(celltypes), n_rows * n_cols):
        fig.delaxes(axs.flatten()[i])

    plt.xticks(rotation=45)
    plt.tight_layout()

    plt.savefig(os.path.join(out_dir, 'proximity_distance_range.png'))
    plt.savefig(os.path.join(out_dir, 'proximity_distance_range.pdf'))
    plt.clf()








