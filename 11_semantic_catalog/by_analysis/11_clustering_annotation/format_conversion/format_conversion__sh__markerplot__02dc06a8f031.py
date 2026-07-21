import os, sys
import numpy as np
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import pandas as pd
from matplotlib import font_manager as fm

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    '''
    out_dir : 输出文件夹地址
    file_name : 文件名称
    fig_size : 输出图片的大小
    '''
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

# adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2508294648-01/02.Mission/Fibroblasts/Main_mapping.h5ad")
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2508294648-01/03.Mission/Pericytes/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2508294648-01/03.Mission/Zbp1_plot_new/"
sc.pl.umap(adata, color=['Zbp1'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_Zbp1.png'))
plt.savefig(os.path.join(out_dir, 'umap_Zbp1.pdf'))
sc.pl.dotplot(adata, var_names=['Zbp1'], groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Zbp1_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Zbp1_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Zbp1'], groupby='cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_Zbp1_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Zbp1_scale.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Zbp1'], groupby=['group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Zbp1_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Zbp1_group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Zbp1'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Zbp1_cell_type_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Zbp1_cell_type_group.pdf'),bbox_inches='tight')