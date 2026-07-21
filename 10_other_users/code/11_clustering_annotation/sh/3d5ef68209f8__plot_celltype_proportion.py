# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def plot_celltype_proportion(_adata, out_dir, group_key='sample', fig_size=None):
    # 获取细胞类型颜色映射
    keys = _adata.obs['cell_type'].cat.categories.to_list()
    values = list(_adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))

    # 计算比例表
    ctab = pd.crosstab(_adata.obs[group_key], _adata.obs['cell_type'], normalize='index')

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

    # 设置坐标轴
    ax.set_xlabel(group_key.capitalize(), fontsize=14, family='sans-serif')
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
        prop=legend_font, frameon=False
    )

    # 调整子图布局，确保 legend 不遮挡主图
    plt.subplots_adjust(right=0.75)

    # 创建输出目录
    os.makedirs(out_dir, exist_ok=True)
    outfile1 = os.path.join(out_dir, f'celltype_proportion_{group_key}.png')
    outfile2 = os.path.join(out_dir, f'celltype_proportion_{group_key}.pdf')
    plt.savefig(outfile1, dpi=300, bbox_inches='tight')
    plt.savefig(outfile2)

    plt.close()
    print(f"✅ 图像已保存至: {outfile1}")

def main():
    parser = argparse.ArgumentParser(description="绘制细胞类型比例堆积图")
    parser.add_argument('--h5ad', required=True, help='输入 h5ad 文件路径')
    parser.add_argument('--outdir', required=True, help='输出目录')
    parser.add_argument('--group_key', default='sample', help='用于分组的 obs 字段名（默认: sample）')
    args = parser.parse_args()

    # 读取 h5ad 文件
    adata = sc.read_h5ad(args.h5ad)

    # 检查字段存在
    if 'cell_type' not in adata.obs.columns:
        raise ValueError("❌ obs 中未找到 'cell_type' 字段")
    if args.group_key not in adata.obs.columns:
        raise ValueError(f"❌ obs 中未找到 '{args.group_key}' 字段")

    plot_celltype_proportion(adata, args.outdir, group_key=args.group_key)

if __name__ == '__main__':
    main()
