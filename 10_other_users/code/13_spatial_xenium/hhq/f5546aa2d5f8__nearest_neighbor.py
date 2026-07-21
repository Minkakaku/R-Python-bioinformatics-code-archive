import warnings
warnings.filterwarnings("ignore", message="No data for colormapping provided via 'c'")
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

import os
import itertools
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from anndata import AnnData
from typing import Optional
from tqdm.auto import tqdm
from statannotations.Annotator import Annotator
import logging as logg
import argparse

format="%(asctime)s-%(levelname)s-%(message)s"
logg.basicConfig(format=format, level=logg.INFO)

plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'
sc.settings.set_figure_params(dpi=100, dpi_save=300, figsize=(7, 7), facecolor='white')


def _normalize_group_col(s: pd.Series) -> pd.Series:
    """Convert group values to clean strings, handling tuples/lists/MultiIndex-like values."""
    def to_str(x):
        if isinstance(x, (tuple, list)):
            if len(x) == 0:
                return ""
            return str(x[0])
        return str(x)
    return s.map(to_str)


def nearest_neighbor_analysis(
    adata: AnnData,
    ident: str = "sample",
    group_key: str = "group",
    celltype_column: str = "celltype",
    output_dir: str = "./"
) -> None:
    """
    计算每个细胞类型到最近邻的距离分布

    Parameters
    ----------
    adata : AnnData
        AnnData object
    ident : str
        adata.obs中区分不同样本的列名
    group_key : str
        adata.obs中区分实验组的列名
    celltype_column : str
        adata.obs中细胞类型注释的列名
    output_dir : str
        输出目录路径
    """
    logg.info("Starting nearest neighbor analysis...")
    
    group_list = adata.obs[group_key].unique()
    palette = dict(zip(group_list, sns.color_palette("Set2", n_colors=len(group_list))))

    all_cell_types = np.unique(adata.obs[celltype_column].values)

    for cell_type in tqdm(all_cell_types, desc="Processing cell types", unit="type"):
        out = os.path.join(output_dir, str(cell_type))
        os.makedirs(out, exist_ok=True)
        proportions = []

        # 按 sample 分组，避免跨样本计算距离
        for sample in adata.obs[ident].unique():
            ad = adata[adata.obs[ident] == sample].copy()

            # 提前缓存数组，提高效率
            groups = ad.obs[group_key].values
            samples = ad.obs[ident].values
            celltypes = ad.obs[celltype_column].values
            coords = ad.obsm['spatial']

            # 构建 KDTree
            tree = KDTree(coords)
            distances, indices = tree.query(coords, k=2)

            # 找到当前 cell_type 的细胞
            mask = (celltypes == cell_type)
            idxs = np.where(mask)[0]

            if len(idxs) == 0:
                continue

            # 生成结果 dataframe
            for i in idxs:
                nn_idx = indices[i, 1]        # 最近邻（排除自己）
                nn_type = celltypes[nn_idx]
                nn_group = groups[nn_idx]
                nn_sample = samples[nn_idx]

                proportions.append([
                    cell_type,
                    nn_type,
                    distances[i, 1],
                    nn_sample,
                    nn_group
                ])

        # 汇总 dataframe
        df = pd.DataFrame(
            proportions,
            columns=['source_type', 'neighbor_type', 'distance', ident, group_key]
        )

        if df.empty:
            logg.warning(f"No data for cell type: {cell_type}")
            continue

        # 输出结果文件
        df.to_csv(os.path.join(out, "nearest_neighbor_proportion.tsv"), sep='\t', index=False)

        # 绘制密度图
        g = sns.FacetGrid(df, col="neighbor_type", hue=group_key, col_wrap=4, palette=palette)
        g.map(sns.kdeplot, "distance", multiple="stack", alpha=0.5, linewidth=0)

        g.set_titles("{col_name}", size=14)
        g.add_legend(title="Group", bbox_to_anchor=(1, 1), loc='upper left')

        for ax in g.axes.flatten():
            ax.yaxis.set_tick_params(labelleft=True)
            ax.grid(False)
            ax.set_xlim(0, 30)

        g.fig.subplots_adjust(wspace=0.3)

        g.savefig(os.path.join(out, "nearest_neighbor_proportion.png"), dpi=300, bbox_inches='tight')
        plt.close(g.fig)

    logg.info("Nearest neighbor analysis completed.")


def plot_distance_group_facet(
    df: pd.DataFrame,
    out: str,
    x_col: str = "group",
    y_col: str = "distance",
    src_col: str = "source_type",
    neighbor_col: str = "neighbor_type",
    test: str = "Mann-Whitney",
    text_format: str = "star",
    hide_non_significant: bool = True,
    col_wrap: int = 4,
    height: float = 4.0,
    aspect: float = 1.0,
    palette: str = "Set2",
    showfliers: bool = False,
    dpi: int = 300,
):
    """
    为每个source_type创建FacetGrid，绘制箱线图并标注p值

    Parameters
    ----------
    df : pd.DataFrame
        包含距离数据的DataFrame
    out : str
        输出目录
    x_col : str
        分组列名
    y_col : str
        距离列名
    src_col : str
        源细胞类型列名
    neighbor_col : str
        邻居细胞类型列名
    test : str
        统计检验方法
    text_format : str
        p值显示格式
    hide_non_significant : bool
        是否隐藏非显著结果
    col_wrap : int
        每行显示的facet数
    height : float
        图高度
    aspect : float
        宽高比
    palette : str
        调色板
    showfliers : bool
        是否显示离群点
    dpi : int
        图片分辨率
    """
    os.makedirs(out, exist_ok=True)

    # Basic validation
    required = {src_col, neighbor_col, x_col, y_col}
    missing = required - set(df.columns)
    if missing:
        logg.warning(f"Missing required columns in df: {sorted(missing)}")
        return

    # Normalize group column
    df = df.copy()
    df[x_col] = _normalize_group_col(df[x_col])

    # Ensure numeric y
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    df = df.dropna(subset=[y_col, x_col, src_col, neighbor_col])

    # Iterate per source type
    for src in df[src_col].unique():
        tmp_src = df[df[src_col] == src].copy()
        if tmp_src.empty:
            continue

        # Create FacetGrid per neighbor_type
        g = sns.FacetGrid(
            tmp_src,
            col=neighbor_col,
            sharey=False,
            col_wrap=col_wrap,
            height=height,
            aspect=aspect,
        )

        # Draw boxplot in each facet
        g.map_dataframe(
            sns.boxplot,
            x=x_col,
            y=y_col,
            palette=palette,
            showfliers=showfliers,
        )

        # Annotate p-values per facet
        for ax, neighbor in zip(g.axes.flat, g.col_names):
            tmp_neighbor = tmp_src[tmp_src[neighbor_col] == neighbor].copy()
            if tmp_neighbor.empty:
                ax.set_title(f"{neighbor} (empty)")
                continue

            tmp_neighbor[x_col] = _normalize_group_col(tmp_neighbor[x_col])

            groups = pd.unique(tmp_neighbor[x_col]).tolist()
            groups = [g_ for g_ in groups if g_ != "" and pd.notna(g_)]

            if len(groups) < 2:
                ax.set_title(f"{neighbor} (n_groups<2)")
                continue

            group_pairs = list(itertools.combinations(groups, 2))
            pairs = [(g1, g2) for g1, g2 in group_pairs]

            try:
                annot = Annotator(
                    ax,
                    pairs,
                    data=tmp_neighbor,
                    x=x_col,
                    y=y_col,
                    order=groups,
                )
                annot.configure(
                    test=test,
                    text_format=text_format,
                    hide_non_significant=hide_non_significant,
                )
                annot.apply_and_annotate()
            except Exception as e:
                logg.warning(f"Annotation failed for {neighbor}: {e}")

        g.set_axis_labels("", y_col.capitalize())
        g.set_titles("{col_name}")
        g.tight_layout()

        out_png = os.path.join(out, "nearest_neighbor_proportion_boxplot.png")
        g.savefig(out_png, dpi=dpi, bbox_inches="tight")
        plt.close(g.fig)


def nearest_neighbor_boxplot(
    adata: AnnData,
    output_dir: str,
    celltype_column: str = "celltype",
    group_key: str = "group",
):
    """
    为每个细胞类型绘制最近邻距离箱线图

    Parameters
    ----------
    adata : AnnData
        AnnData object
    output_dir : str
        输出目录
    celltype_column : str
        细胞类型列名
    group_key : str
        分组列名
    """
    logg.info("Generating nearest neighbor boxplots...")
    sns.set_theme(style="white")

    cell_types = pd.unique(adata.obs[celltype_column]).tolist()

    for cell_type in tqdm(cell_types, desc="Plotting boxplots", unit="type"):
        cell_type_str = str(cell_type)
        out = os.path.join(output_dir, cell_type_str)
        
        tsv_path = os.path.join(out, "nearest_neighbor_proportion.tsv")
        if not os.path.exists(tsv_path):
            logg.warning(f"File not found: {tsv_path}")
            continue

        df = pd.read_csv(tsv_path, sep="\t")
        plot_distance_group_facet(df, out=out, x_col=group_key)

    logg.info("Boxplot generation completed.")


# --- Snakemake 统一入口 ---
def run_nearest_neighbor_wrapper(
    input_h5ad: AnnData,
    output_dir: str,
    ident: str = "sample",
    group_key: str = "group",
    celltype_column: str = "celltype"
):
    """
    Snakemake统一入口函数，执行最近邻分析

    Parameters
    ----------
    input_h5ad : AnnData
        已读取的AnnData对象
    output_dir : str
        输出目录
    ident : str
        样本ID列名
    group_key : str
        分组列名
    celltype_column : str
        细胞类型列名
    """
    # 最近邻分析
    nearest_neighbor_analysis(
        adata=input_h5ad,
        ident=ident,
        group_key=group_key,
        celltype_column=celltype_column,
        output_dir=output_dir
    )
    
    # 生成箱线图
    nearest_neighbor_boxplot(
        adata=input_h5ad,
        output_dir=output_dir,
        celltype_column=celltype_column,
        group_key=group_key
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Nearest Neighbor Analysis")
    parser.add_argument("--input_h5ad", type=str, required=True, help="Input AnnData file (.h5ad)")
    parser.add_argument("--output_dir", type=str, default="./", help="Output directory")
    parser.add_argument("--ident", type=str, default="sample", help="Column name for sample ID")
    parser.add_argument("--group_key", type=str, default="group", help="Column name for experimental groups")
    parser.add_argument("--celltype_column", type=str, default="celltype", help="Column name for cell type annotation")

    args = parser.parse_args()

    # 读取数据
    logg.info(f"Reading data from {args.input_h5ad}...")
    adata = sc.read_h5ad(args.input_h5ad)

    os.makedirs(args.output_dir, exist_ok=True)

    # 执行分析
    run_nearest_neighbor_wrapper(
        input_h5ad=adata,
        output_dir=args.output_dir,
        ident=args.ident,
        group_key=args.group_key,
        celltype_column=args.celltype_column
    )

    done_file = os.path.join(args.output_dir, "done")
    open(done_file, "w").close()
    logg.info("Done.")