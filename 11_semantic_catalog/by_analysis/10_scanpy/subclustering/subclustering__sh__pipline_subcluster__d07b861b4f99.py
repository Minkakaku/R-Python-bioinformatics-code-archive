import os, sys
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["OMP_NUM_THREADS"] = "16"
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
import pandas as pd
from matplotlib import font_manager as fm
import seaborn as sns

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

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/tmp/"

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/B-Plasma/B_Plasma.h5ad")
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))

bdata1 = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/02.Mission/Plasma cells/Plasma.h5ad")
bdata2 = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/02.Mission/B cells/Bcells.h5ad")
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(bdata1, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_oldp.png'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(bdata2, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_oldB.png'))
bdata1.obs['cell_type'] = 'old_' + bdata1.obs['cell_type'].astype(str)
bdata2.obs['cell_type'] = 'old_' + bdata2.obs['cell_type'].astype(str)

adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
cell_type_mapping = bdata1.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
cell_type_mapping = bdata2.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')