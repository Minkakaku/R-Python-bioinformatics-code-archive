import warnings
warnings.filterwarnings("ignore", message="No data for colormapping provided via 'c'")
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import os
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from scipy.spatial import KDTree
from anndata import AnnData
from typing import Literal, Optional, Union, List
from tqdm.auto import tqdm
import logging as logg
from itertools import combinations
import argparse
from statannotations.Annotator import Annotator

format="%(asctime)s-%(levelname)s-%(message)s"
logg.basicConfig(format=format, level=logg.INFO)

plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
sc.settings.set_figure_params(dpi=100, dpi_save=300, figsize=(7, 7), facecolor='white')


def _auto_figsize_by_celltypes(
    celltypes,
    min_width: float = 8,
    width_per_celltype: float = 0.7,
    padding_width: float = 4,
    max_width: float = 80,
    min_height: float = 6,
    height_per_celltype: float = 0.12,
    padding_height: float = 5,
    max_height: float = 60
) -> tuple:
    """根据细胞类型数量自动计算图片宽度和高度。"""
    n_celltypes = max(len(pd.unique(celltypes)), 1)
    width = padding_width + n_celltypes * width_per_celltype
    height = padding_height + n_celltypes * height_per_celltype
    width = min(max(width, min_width), max_width)
    height = min(max(height, min_height), max_height)
    return (width, height)


def calculate_proximity(
    adata: AnnData,
    target_celltype: str,
    ident: str = "sample",
    group_key: str = "group",
    celltype_column: str = 'celltype',
    distance_threshold: int = 50,
    output_dir: str = "./",
    fig_size: Optional[tuple] = None
) -> pd.DataFrame:
    """
    计算目标细胞类型周围一定距离内各细胞类型的比例（样本级别聚合）

    Parameters
    ----------
    adata : AnnData
        AnnData object
    target_celltype : str
        目标细胞类型
    ident : str
        adata.obs中区分不同样本的列名
    group_key : str
        adata.obs中区分实验组的列名
    celltype_column : str
        adata.obs中细胞类型注释的列名
    distance_threshold : int
        邻居检测半径(微米)
    output_dir : str
        输出目录路径
    fig_size : tuple, optional
        图形尺寸；默认根据细胞类型数量自动调整图片宽度和高度

    Returns
    -------
    sample_df : pd.DataFrame
        包含样本级别proximity分析结果的DataFrame
    """
    logg.info(f"Calculating proximity for {target_celltype} with distance threshold {distance_threshold}um...")
    
    group_list = adata.obs[group_key].unique()
    palette = dict(zip(group_list, sns.color_palette("Set2", n_colors=len(group_list))))
    
    all_types = adata.obs[celltype_column].unique()
    proximity_data = []
    
    samples = adata.obs[ident].unique()
    for sample in tqdm(samples, desc="Processing samples", unit="sample"):
        ad_s = adata[adata.obs[ident] == sample].copy()
        sq.gr.spatial_neighbors(ad_s, coord_type='generic', radius=distance_threshold)
        cells_of_type = ad_s.obs_names[ad_s.obs[celltype_column] == target_celltype]

        for cell_label in cells_of_type:
            cell_idx = ad_s.obs_names.get_loc(cell_label)
            neighbors = ad_s.obsp['spatial_connectivities'][cell_idx].nonzero()[1]
            neighbor_labels = ad_s.obs_names[neighbors]
            neighbor_types = ad_s.obs.loc[neighbor_labels, celltype_column]
            # 去掉自己
            # neighbor_types = neighbor_types[neighbor_types.index != cell_label]
            counts = neighbor_types.value_counts(normalize=True).reindex(all_types, fill_value=0.0)
            for ct, prop in counts.items():
                proximity_data.append({
                    'celltype': ct,
                    'proportion': prop,
                    ident: sample,
                    group_key: ad_s.obs.loc[cell_label, group_key]
                })

    # 细胞级别结果
    cell_df = pd.DataFrame(proximity_data)
    cell_df.to_csv(os.path.join(output_dir, f'{target_celltype}_nearly_{distance_threshold}um_cell_level.xls'), sep='\t', index=None)

    # 样本级别聚合：计算每个样本中每种细胞类型的平均比例
    sample_df = cell_df.groupby([ident, 'celltype', group_key], as_index=False)['proportion'].mean()
    sample_df.to_csv(os.path.join(output_dir, f'{target_celltype}_nearly_{distance_threshold}um_sample_level.xls'), sep='\t', index=None)

    plot_fig_size = fig_size or _auto_figsize_by_celltypes(
        cell_df['celltype'],
        min_height=8,
        padding_height=7,
        height_per_celltype=0.15
    )

    # ========== 细胞级别绘图 ==========
    logg.info(f"Plotting proximity boxplot (cell-level)...")
    plt.figure(figsize=plot_fig_size)
    ax_cell = sns.boxplot(
        x='celltype', 
        y='proportion', 
        hue=group_key, 
        palette=palette,
        data=cell_df, 
        showfliers=False
    )
    ax_cell.grid(visible=False)
    ax_cell.legend(title="", bbox_to_anchor=(1.02, 0.5), loc='upper left', borderaxespad=0.)
    plt.xlabel("")
    plt.title(f'Celltype proportion of {target_celltype} nearly {distance_threshold}um distance (cell-level)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{target_celltype}_proximity_{distance_threshold}um_cell_level.png'))
    plt.savefig(os.path.join(output_dir, f'{target_celltype}_proximity_{distance_threshold}um_cell_level.pdf'))
    plt.close()

    # ========== 样本级别绘图 ==========
    logg.info(f"Plotting proximity boxplot (sample-level)...")
    plt.figure(figsize=plot_fig_size)
    ax = sns.boxplot(
        x='celltype', 
        y='proportion', 
        hue=group_key, 
        palette=palette,
        data=sample_df, 
        showfliers=False
    )
    sns.stripplot(
        x='celltype', 
        y='proportion', 
        hue=group_key, 
        palette=palette,
        data=sample_df, 
        dodge=True,
        alpha=0.6,
        size=5,
        ax=ax,
        legend=False
    )
    ax.grid(visible=False)

    # 进行统计检验并添加 p 值（样本级别，避免伪重复）
    if group_key and group_key in adata.obs.columns:
        unique_groups = sample_df[group_key].unique()
        if len(unique_groups) >= 2:
            box_pairs = []
            celltypes = sample_df['celltype'].unique()
            # 遍历每一个 CellType
            for ct in celltypes:
                # 在当前 CellType 内，生成组的所有两两组合
                for g1, g2 in combinations(unique_groups, 2):
                    box_pairs.append(((ct, g1), (ct, g2)))
            try:
                annotator = Annotator(
                    ax, data=sample_df, x='celltype', y='proportion', hue=group_key,
                    pairs=box_pairs,
                    verbose=0,
                    comparisons_correction='BH' # 多组比较建议加上校正，防止假阳性
                )
                annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', hide_non_significant=True, line_height=0.005, line_width=1.2)
                annotator.apply_and_annotate()
            except Exception as e:
                logg.warning(f"Stat annotation skipped: {e}")

    ax.legend(title="", bbox_to_anchor=(1.02, 0.5), loc='upper left', borderaxespad=0.)
    
    plt.xlabel("")
    plt.title(f'Celltype proportion of {target_celltype} nearly {distance_threshold}um distance (sample-level)')
    plt.xticks(rotation=45)
    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, f'{target_celltype}_proximity_{distance_threshold}um_sample_level.png'))
    plt.savefig(os.path.join(output_dir, f'{target_celltype}_proximity_{distance_threshold}um_sample_level.pdf'))
    plt.close()

    return cell_df, sample_df


def calculate_proximity_range(
    adata: AnnData,
    target_celltype: str,
    ident: str = "sample",
    group_key: str = "group",
    celltype_column: str = 'celltype',
    distance_ranges: List[int] = [0, 50, 100, 150, 200],
    output_dir: str = "./",
    figsize: Optional[tuple] = None
) -> pd.DataFrame:
    """
    计算目标细胞类型在不同距离区间内邻居细胞类型的数量。

    Parameters
    ----------
    adata : AnnData
        AnnData object
    target_celltype : str
        目标细胞类型
    ident : str
        adata.obs中区分不同样本的列名
    group_key : str
        adata.obs中区分实验组的列名
    celltype_column : str
        adata.obs中细胞类型注释的列名
    distance_ranges : List[int]
        距离区间列表，例如 [0, 50, 100, 150, 200]
    output_dir : str
        输出目录路径
    figsize : tuple, optional
        图形尺寸；默认根据细胞类型数量自动调整图片宽度和高度

    Returns
    -------
    mean_proximity_df : pd.DataFrame
        按距离区间、细胞类型和分组汇总的平均邻居数量
    """
    logg.info(f"Calculating proximity range for {target_celltype}...")
    
    proximity_data = []
    samples = np.unique(adata.obs[ident].values)
    
    for sample in tqdm(samples, desc="Processing samples", unit="sample"):
        ad = adata[adata.obs[ident] == sample].copy()
        coords = ad.obsm['spatial']
        cell_types = np.asarray(ad.obs[celltype_column])
        groups = np.asarray(ad.obs[group_key])

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
                        group_key: groups[idx]
                    })

    # 转换为 DataFrame
    proximity_df = pd.DataFrame(proximity_data)

    # 按照 distance_range 和 celltype 分组，并计算均值
    mean_proximity_df = proximity_df.groupby(['distance_range', 'celltype', group_key], as_index=False).mean()

    # 构建距离区间类别
    distance_categories = [f"{distance_ranges[i]}-{distance_ranges[i+1]}" for i in range(len(distance_ranges)-1)]
    mean_proximity_df['distance_range'] = pd.Categorical(
        mean_proximity_df['distance_range'], 
        categories=distance_categories, 
        ordered=True
    )
    mean_proximity_df.sort_values('distance_range', inplace=True)
    mean_proximity_df.to_csv(os.path.join(output_dir, f'{target_celltype}_proximity_distance_range.txt'), sep='\t', index=None)

    # 绘制折线图，按细胞类型分面展示，分成4列
    logg.info(f"Plotting proximity range lineplot...")
    celltypes = sorted(mean_proximity_df['celltype'].unique().tolist())
    n_cols = 4  # 每行显示4个子图
    n_rows = (len(celltypes) + n_cols - 1) // n_cols  # 计算行数
    # distance_range 是分面图：宽度按每行列数控制，高度按行数增加，避免细胞类型多时横向过长。
    plot_figsize = figsize or (max(4.2 * n_cols, 12), max(5.2 * n_rows, 5.2))

    fig, axs = plt.subplots(n_rows, n_cols, figsize=plot_figsize, sharex=True)

    for ax, celltype in zip(axs.flatten(), celltypes):
        subset = mean_proximity_df[mean_proximity_df['celltype'] == celltype]
        for group in sorted(subset[group_key].unique().tolist()):
            group_data = subset[subset[group_key] == group]
            ax.plot(group_data['distance_range'], group_data['count'], marker='o', label=group)
        
        ax.set_title(f'{celltype}')
        ax.set_ylabel('Cell count')
        ax.legend()

    # 隐藏多余的子图（如果有）
    for i in range(len(celltypes), n_rows * n_cols):
        fig.delaxes(axs.flatten()[i])

    plt.xticks(rotation=45)
    plt.tight_layout()

    plt.savefig(os.path.join(output_dir, f'{target_celltype}_proximity_distance_range.png'))
    plt.savefig(os.path.join(output_dir, f'{target_celltype}_proximity_distance_range.pdf'))
    plt.close()

    return mean_proximity_df


# --- Snakemake 统一入口 ---
def run_proximity_wrapper(
    input_h5ad: AnnData,
    output_dir: str,
    target_celltype: str,
    ident: str = 'sample',
    group_key: str = 'group',
    celltype_column: str = 'celltype',
    distance_threshold: int = 50,
    distance_ranges: List[int] = [0, 50, 100, 150, 200]
):
    """
    Snakemake统一入口函数，执行proximity分析

    Parameters
    ----------
    input_h5ad : AnnData
        已读取的AnnData对象
    output_dir : str
        输出目录
    target_celltype : str
        目标细胞类型
    ident : str
        样本ID列名
    group_key : str
        分组列名
    celltype_column : str
        细胞类型列名
    distance_threshold : int
        proximity分析的距离阈值
    distance_ranges : List[int]
        proximity range分析的距离区间
    """
    calculate_proximity(
        adata=input_h5ad,
        target_celltype=target_celltype,
        ident=ident,
        group_key=group_key,
        celltype_column=celltype_column,
        distance_threshold=distance_threshold,
        output_dir=output_dir
    )
    
    calculate_proximity_range(
        adata=input_h5ad,
        target_celltype=target_celltype,
        ident=ident,
        group_key=group_key,
        celltype_column=celltype_column,
        distance_ranges=distance_ranges,
        output_dir=output_dir
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Proximity Analysis")
    parser.add_argument("--input_h5ad", type=str, required=True, help="Input AnnData file (.h5ad)")
    parser.add_argument("--output_dir", type=str, default="./", help="Output directory")
    parser.add_argument("--target_celltype", type=str, required=True, nargs="+", help="Target cell type(s) for analysis")
    parser.add_argument("--ident", type=str, default="sample", help="Column name for sample ID")
    parser.add_argument("--group_key", type=str, default="group", help="Column name for experimental groups")
    parser.add_argument("--celltype_column", type=str, default="celltype", help="Column name for cell type annotations")
    parser.add_argument("--distance_threshold", type=int, default=50, help="Distance threshold for proximity analysis (um)")
    parser.add_argument("--distance_ranges", type=int, nargs="+", default=[0, 50, 100, 150, 200], help="Distance ranges for proximity range analysis")


    args = parser.parse_args()

    # 读取数据
    logg.info(f"Reading data from {args.input_h5ad}...")
    adata = sc.read_h5ad(args.input_h5ad)

    os.makedirs(args.output_dir, exist_ok=True)

    # 处理 target_celltype 为 'all' 的情况
    if "all" in args.target_celltype:
        target_celltypes = sorted(adata.obs[args.celltype_column].unique().tolist())
        logg.info(f"Detected 'all' in target_celltype. Processing all {len(target_celltypes)} cell types.")
    else:
        target_celltypes = args.target_celltype

    # 多目标细胞类型循环
    for target_ct in target_celltypes:
        logg.info(f"Running proximity analysis for target celltype: {target_ct}...")
        run_proximity_wrapper(
            input_h5ad=adata,
            output_dir=args.output_dir,
            target_celltype=target_ct,
            ident=args.ident,
            group_key=args.group_key,
            celltype_column=args.celltype_column,
            distance_threshold=args.distance_threshold,
            distance_ranges=args.distance_ranges
        )

    done_file = os.path.join(args.output_dir, "done")
    open(done_file, "w").close()
    logg.info("Done.")