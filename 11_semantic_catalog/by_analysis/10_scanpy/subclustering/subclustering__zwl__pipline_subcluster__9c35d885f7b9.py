import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
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

adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2512027077-01/01.Mission/Main_mapping/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512027077-01/05.Mission/gene_plot/"
sc.pl.umap(adata, color=["App","Cd74"], cmap='coolwarm', ncols=3)
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby='sub_cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sub_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sub_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby='group', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby='sample', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sample.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_sample.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=["App","Cd74"], groupby=['cell_type','sample'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-sample.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type-sample.pdf'),bbox_inches='tight')