import os
import pandas as pd
import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import matplotlib.font_manager as fm
import seaborn as sns

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
        "#ea7070", "#fdc4b6",  "#2694ab", "#96ceb4", "#ffeead", 
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

def plot_celltype_proportion(adata, out_dir,type='sample', fig_size=None):
    # 获取细胞类型颜色映射
    keys = adata.obs['cell_type'].cat.categories.to_list()
    values = list(adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))

    # 创建堆积条形图
    ax = pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').plot(
        kind='bar', xlabel='', stacked=True,
        color=[cell_type_colors[cell] for cell in pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').columns],
        edgecolor='black', width=0.8
    )

    # 设置X轴和Y轴标签
    ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')  # 设置X轴标签
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')  # 设置Y轴标签

    # 设置图例
    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=10, title="Cell Types", 
        title_fontsize=12, prop=legend_font, frameon=False, handletextpad=1.5, labelspacing=1.0, 
        columnspacing=1.0, markerscale=1,  handlelength=2, handleheight=2
    )

    # 设置坐标轴标签的字体大小和颜色
    ax.tick_params(axis='x', labelsize=12, labelrotation=45, labelcolor='black')
    ax.tick_params(axis='y', labelsize=12, labelcolor='black')

    # 设置背景颜色
    ax.set_facecolor('white')

    # 设置图表边框加粗
    for spine in ax.spines.values():
        spine.set_edgecolor('black')  # 设置边框颜色
        spine.set_linewidth(1.5)

    # 设置图像大小，如果没有传入fig_size，则使用默认大小
    if fig_size is None:
        fig_size = {'width': len(set(adata.obs[type])) * 2, 'height': 6}

    # 调整布局，确保图例不遮挡图形
    plt.tight_layout()
    file_name=f"celltype_proportion_{type}"
    # 保存图像
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)


def celltype_annotation(adata, out_dir, annotation_dict):
    set_scanpy_color_palette(adata, color_key="cell_type")
    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.umap(adata, color='cell_type', show=False, title='', legend_loc='on data', legend_fontsize=9, legend_fontweight='normal', legend_fontoutline=2)
        plt.savefig(os.path.join(out_dir, 'umap_dimplot_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)

    groups = sorted(adata.obs['group'].drop_duplicates().tolist())
    if len(groups) > 1:
        axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
        fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
        for i in range(len(groups)):
            if i == len(groups) - 1:
                legend = 'right margin'
            else:
                legend = None
            sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", title=groups[i], ax=axs[i], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        save_fig(plt, out_dir, file_name='umap_split_group_celltype', fig_size={'width':len(groups)*5, 'height':5})

    result = adata.obs.groupby(['cell_type', 'sample']).size().reset_index(name='number').pivot(index='cell_type', columns='sample', values='number').reset_index()
    result.columns.name = None
    # result['Total'] = result.sum(axis=1)
    # result['Percent'] = round((result['Total'] / sum(result['Total'])) * 100, 1)

    result.to_csv(os.path.join(out_dir, 'celltype_statistics.xls'), sep='\t', index=None)

    adata.obs.to_csv(os.path.join(out_dir, 'all_cells_features.xls'), sep='\t')
    plot_celltype_proportion(adata, out_dir=out_dir,type='group')
    plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
    pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['UMAP_1', 'UMAP_2']).to_csv(os.path.join(out_dir, 'data_umap.xls'), sep='\t')

    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.umap(adata, color=["cluster","cell_type"], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
        plt.savefig(os.path.join(out_dir, 'umap_cluster_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
        plt.savefig(os.path.join(out_dir, 'umap_cluster_celltype.pdf'))

    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.dotplot(adata, genes,groupby='cell_type',cmap='coolwarm', use_raw=True, show=False,standard_scale='var')
        plt.savefig(os.path.join(out_dir, 'anno_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=200)
        plt.savefig(os.path.join(out_dir, 'anno_dotplot.pdf'))
    if organism == 'mouse':
        adata.var.index = adata.var['gene']

    adata.write(os.path.join(out_dir, 'result.h5ad'))

    return adata, annotation_dict

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="annotation")
    parser.add_argument("--input_file", required=True, help="输入文件 h5ad")
    parser.add_argument("--out_dir", required=True, help="输出文件夹路径")
    parser.add_argument("--anno_file", required=True, help="输入注释文件")
    args = parser.parse_args()
    anno_out_dir = os.path.join(args.out_dir, '06-Cell_type_annotation')
    if not os.path.exists(anno_out_dir):
        os.makedirs(anno_out_dir)
    adata = sc.read_h5ad(args.input_file)
    df = pd.read_csv("args.anno_file", sep="\t")
    anno_dict = dict(zip(df["cluster"], param_df["cell_type"]))
    adata_final,annotation_dict=celltype_annotation(adata, anno_out_dir, anno_dict)
    anno_df = pd.DataFrame(list(annotation_dict.items()), columns=['cluster', 'cell_type'])
    anno_df.to_csv(os.path.join(anno_out_dir, 'annotation_cluster.xls'), sep='\t', index=False)