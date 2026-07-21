from matplotlib import font_manager as fm
import os
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()


def plot_celltype_proportion(_adata, out_dir,type='sample', fig_size=None):
    # 获取细胞类型颜色映射
    keys = _adata.obs['cell_type'].cat.categories.to_list()
    values = list(_adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))

    # 创建堆积条形图
    ax = pd.crosstab(_adata.obs[type], _adata.obs['cell_type'], normalize='index').plot(
        kind='bar', xlabel='', stacked=True,
        color=[cell_type_colors[cell] for cell in pd.crosstab(_adata.obs[type], _adata.obs['cell_type'], normalize='index').columns],
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
        columnspacing=1.0, markerscale=1, alignment='left', handlelength=2, handleheight=2
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
        fig_size = {'width': len(set(_adata.obs[type])) * 2, 'height': 6}

    # 调整布局，确保图例不遮挡图形
    plt.tight_layout()
    file_name=f"celltype_proportion_{type}"
    # 保存图像
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)
