import os
os.environ['OPENBLAS_NUM_THREADS'] = "64"
os.environ['MKL_NUM_THREADS'] = "64"
os.environ['OMP_NUM_THREADS'] = "64"
import shutil
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scanpy.external as sce
from matplotlib.pyplot import rc_context
from matplotlib import font_manager as fm
sc.settings.verbosity = 3
sc.settings.set_figure_params(figsize=(5,5), facecolor='white')
plt.tight_layout()

def save_fig(plt, out_dir, file_name, fig_size={'width':5, 'height':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['width'], fig_size['height'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def set_scanpy_color_palette(adata, color_key, color_palette=None):
    """
    为 Scanpy 对象设置统一配色方案。
    
    参数:
    - adata: AnnData 对象
    - color_key: str，obs 中的列名，用于分组着色（如 'cell_type'）
    - color_palette: list，可选，自定义的颜色列表（最多 50 个）

    默认色板支持最多 50 种颜色，超出会报错。
    """
    default_palette = [
        "#ea7070", "#fdc4b6", "#e59572", "#2694ab", "#96ceb4", "#ffeead", 
        "#d9534f", "#ffad60", "#57D1C9", "#ED5485", "#7F95D1", "#FFE869", 
        "#de4307", "#f29c2b", "#f6d04d", "#8bc24c", "#4695d6", "#fed95c", 
        "#fa6e57", "#f69e53", "#248888", "#999999", "#E7475E", "#F0D879", 
        "#29A2C6", "#FFCB18", "#73B66B", "#FF6D31", "#f1ac9d", "#f06966", 
        "#dee2d1", "#6abe83", "#FCF4D9", "#FFB85F", "#FF7A5A", "#8ED2C9", 
        "#a5dff9", "#ef5285", "#60c5ba", "#feee7d", "#e8a0b8", "#ffc300", 
        "#bccf3d", "#02c9c9", "#fbf579", "#005995", "#fa625f", "#600473", 
        "#f17d80", "#737495"
    ]
    palette = color_palette if color_palette else default_palette

    if color_key not in adata.obs:
        raise KeyError(f"'{color_key}' 不存在于 adata.obs 中")

    categories = sorted(adata.obs[color_key].astype(str).unique().tolist())
    n_categories = len(categories)

    if n_categories > len(palette):
        raise ValueError(f"类别数为 {n_categories}，超过了可用颜色数 {len(palette)}，请扩展色板。")

    # 设置为 Categorical 类型并指定顺序
    adata.obs[color_key] = adata.obs[color_key].astype("category")
    adata.obs[color_key] = adata.obs[color_key].cat.set_categories(categories)

    # 设置颜色
    adata.uns[f"{color_key}_colors"] = palette[:n_categories]

import matplotlib as mpl

mpl.rcParams.update({
    'legend.handleheight': 1.2,    # 全局色块高度
    'legend.handlelength': 1.8,    # 全局色块长度
    'legend.labelspacing': 0.4,     # 全局标签间距
    'legend.handletextpad': 0.4,    # 全局色块文本间距
    'legend.borderpad': 0.6,        # 全局内边距
})

# 如果使用seaborn，需要重置
import seaborn as sns
sns.reset_orig()  # 清除seaborn样式影响
def plot_celltype_proportion(adata, out_dir, type='sample', fig_size=None):
    # 获取细胞类型颜色映射
    keys = adata.obs['cell_type'].cat.categories.to_list()
    values = list(adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))

    # 计算比例表
    ctab = pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index')

    # 设置图像尺寸
    if fig_size is None:
        fig_width = max(len(ctab.index) * 1.5, 6)
        fig_size = (fig_width, 6)

    # 创建图像和子图
    fig, ax = plt.subplots(figsize=fig_size)

    # 绘制堆积条形图
    ctab.plot(kind='bar', stacked=True, ax=ax,
              color=[cell_type_colors[cell] for cell in ctab.columns],
              edgecolor='black', width=0.8)
    ax.grid(False)
    # 设置坐标轴
    ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')
    ax.tick_params(axis='x', labelsize=12, rotation=90)
    ax.tick_params(axis='y', labelsize=12)

    # 加粗边框
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)

    # 设置图例
    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(
        bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0.,
        fontsize=10, title="Cell Types", title_fontsize=12,
        prop=legend_font, frameon=False,
    handleheight=1.5,      # 增加色块高度（默认0.7）
    handlelength=1.5,      # 增加色块长度
    labelspacing=0.5,      # 增加标签间距
    handletextpad=0.5,     # 增加色块与文本间距
    borderpad=0.8,         # 增加图例内边距
    )

    # 调整子图布局，确保 legend 不遮挡主图
    plt.subplots_adjust(right=0.75)
    file_name=f"celltype_proportion_{type}"
    # 保存图像
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)
def dim_reduction_cluster(adata, method, n_pcs, n_neighbors, resolution, batch_effect, out_dir):
    adata.layers["counts"] = adata.X.copy()
    # data normalization and log10 transform
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata
    # select highly variable gene and data scale 
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
    sc.pl.highly_variable_genes(adata, show=False)
    save_fig(plt, out_dir, file_name='highly_variable_genes', fig_size={'width':16, 'height':8})
    adata = adata[:, adata.var['highly_variable']]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca(adata, color='sample', show=False)
    save_fig(plt, out_dir, file_name='pca_sample')
    sc.pl.pca(adata, color='group', show=False)
    save_fig(plt, out_dir, file_name='pca_group')
    sc.pl.pca_variance_ratio(adata, log=True, show=False)
    save_fig(plt, out_dir, file_name='pca_variance_ratio')

    if batch_effect == 'harmony':
        sce.pp.harmony_integrate(adata, 'sample', adjusted_basis='X_pca')
    # cluster
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)    

    # Remove batch effect
    if batch_effect == 'bbknn':
        sce.pp.bbknn(adata, batch_key='sample', n_pcs=n_pcs)  # running bbknn 1.3.6


    sc.tl.umap(adata)
    if method == 'leiden':
        sc.tl.leiden(adata, resolution=resolution, key_added='cluster')
    if method == 'louvain':
        sc.tl.louvain(adata, resolution=resolution, key_added='cluster')



    # UMAP
    set_scanpy_color_palette(adata, color_key="cluster")
    sc.pl.umap(adata, 
               color='cluster',
               title='Cluster',
               legend_loc='on data',
               legend_fontsize=10,
               legend_fontoutline=2, 
               show=False)
    save_fig(plt, out_dir, file_name='umap_dimplot')
    sc.pl.umap(adata, color='sample', title='', show=False)
    save_fig(plt, out_dir, file_name='umap_sample_dimplot')
    sc.pl.umap(adata, color='group', title='', show=False)
    save_fig(plt, out_dir, file_name='umap_group_dimplot')
      
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(adata, 
                   color=['cluster', "total_counts", "n_genes_by_counts", "total_counts_mt"],
                   title=['Cluster', 'total_counts', 'n_genes_by_counts', "total_counts_mt"],
                   legend_loc='on data',
                   frameon=True,
                   legend_fontsize=9, 
                   legend_fontoutline=3, 
                   wspace=0.2, 
                   ncols=2, 
                   cmap='coolwarm', 
                   use_raw=True, 
                   show=False)
               
        plt.savefig(os.path.join(out_dir, 'cluster_total_counts_n_genes_by_counts.png'), facecolor='white', bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(out_dir, 'cluster_total_counts_n_genes_by_counts.pdf'))

    
    # write all cells features
    adata.obs.to_csv(os.path.join(out_dir,'all_cells_features.xls'), sep='\t')
    
    # write cluster cell statistics result
    result = adata.obs.groupby(['cluster', 'sample']).size().reset_index(name='number').pivot(index='cluster', columns='sample', values='number').reset_index()
    result.columns.name = None    
    result.to_csv(os.path.join(out_dir,'cluster_cells_statistics.xls'), sep='\t', index=None)
    
    # cluster correlation
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.correlation_matrix(adata, 'cluster', show=False)
        plt.savefig(os.path.join(out_dir, 'cluster_correlation.png'), facecolor='white', bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(out_dir, 'cluster_correlation.pdf'))


    groups = sorted(adata.obs['group'].drop_duplicates().tolist())
    if len(groups) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
        for i in range(len(groups)):
            if i == len(groups) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cluster", title=groups[i], ax=axs[i], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        save_fig(plt, out_dir, file_name='umap_split_group_cluster', fig_size={'width':len(groups)*5, 'height':5})

    samples = sorted(adata.obs['sample'].drop_duplicates().tolist())
    if len(samples) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
        fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
        for i in range(len(samples)):
            if i == len(samples) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(adata[adata.obs["sample"]==samples[i]], color="cluster", title=samples[i], ax=axs[i], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        save_fig(plt, out_dir, file_name='umap_split_sample_cluster', fig_size={'width':len(samples)*5, 'height':5})

    return adata


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Dimension reduction and clustering")
    parser.add_argument("--input_file", required=True, help="输入文件 h5ad")
    parser.add_argument("--out_dir", required=True, help="输出文件夹路径")
    parser.add_argument("--method", default="leiden", help="聚类方法：leiden 或 louvain")
    parser.add_argument("--n_pcs", type=int, default=30, help="PCA主成分数")
    parser.add_argument("--n_neighbors", type=int, default=20, help="邻居数量")
    parser.add_argument("--resolution", type=float, default=0.5, help="聚类分辨率")
    parser.add_argument("--batch_effect", default="none", help="批次校正方法：none/harmony/bbknn")
    args = parser.parse_args()
    Dimension_reduction_cluster_out_dir = os.path.join(args.out_dir, '04-Dimension_reduction_cluster')
    if not os.path.exists(Dimension_reduction_cluster_out_dir):
        os.makedirs(Dimension_reduction_cluster_out_dir)
    else:
        os.system(f'rm -rf {Dimension_reduction_cluster_out_dir}')
        os.makedirs(Dimension_reduction_cluster_out_dir)
    adata = sc.read_h5ad(args.input_file)
    adata_clustered = dim_reduction_cluster(adata, args.method, args.n_pcs, args.n_neighbors, args.resolution, args.batch_effect, Dimension_reduction_cluster_out_dir)
    adata_clustered.write_h5ad(os.path.join(Dimension_reduction_cluster_out_dir,'adata_clustered.h5ad'))

    csv_path = os.path.join(args.out_dir, "Parameter_dict.csv")
    if os.path.exists(csv_path):
        param_df = pd.read_csv(csv_path)
        Parameter_dict = dict(zip(param_df["Key"], param_df["Value"]))
    else:
        Parameter_dict = {}

    VALUE_KEYS = [
    "method",
    "n_pcs",
    "n_neighbors",
    "resolution",
    "batch_effect"
    ]

    run_info = {k: getattr(args, k) for k in VALUE_KEYS}
    Parameter_dict.update(run_info)

    pd.DataFrame(
        Parameter_dict.items(),
        columns=["Key", "Value"]
        ).to_csv(csv_path, index=False)
