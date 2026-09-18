import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["BLIS_NUM_THREADS"] = "1"
import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib import font_manager as fm

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

# load functions
import sys
sys.path.append('/mnt/NexData_file/home/xzx/scripts/scanpy_fun')
from subcluster import recluster
from find_marker import find_marker_gene
from draw_group import draw_group
from draw_bar import plot_celltype_proportion

import pandas as pd
# adata=sc.read_h5ad('/mnt/NexData_file/scRNA/scVision_output/AE2508064183-01/06-Cell_type_annotation/result.h5ad')
# outdir = '/mnt/NexData_file/home/xzx/project/AE2508064183/shouhou6/recluster'
# os.makedirs(outdir,exist_ok=True)

# cellnames = ['Neutrophil']
# for cells in cellnames:
#     celltype = cells
#     outdir1 = os.path.join(outdir,celltype)
#     os.makedirs(outdir1,exist_ok=True)
#     adata_tmp = adata[adata.obs['cell_type'].isin([celltype])].copy()
#     #sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
#     adata_tmp = recluster(adata_tmp,species='human',resolution = 0.5,n_pcs = 10,neighbors = 15,batch_effect=2,output=outdir1)
#     adata_tmp = adata_tmp.sub_cluster()

###过滤重分群
adata=sc.read_h5ad('/mnt/NexData_file/home/xzx/project/AE2508064183/shouhou6/recluster/Neutrophil/result.h5ad')
outdir = '/mnt/NexData_file/home/xzx/project/AE2508064183/shouhou6/recluster'
#os.makedirs(outdir,exist_ok=True)
cellnames = ['Neutrophil_filter']
for cells in cellnames:
    celltype = cells
    outdir1 = os.path.join(outdir,celltype)
    os.makedirs(outdir1,exist_ok=True)
    adata_tmp = adata[~adata.obs['leiden'].isin(['6','7'])].copy()
    #sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
    adata_tmp = recluster(adata_tmp,species='human',resolution = 0.5,n_pcs = 8,neighbors = 15,batch_effect=2,output=outdir1)
    adata_tmp = adata_tmp.sub_cluster()


### 髓系
# adata=sc.read_h5ad('/mnt/NexData_file/home/xzx/project/AE2508064183/shouhou5/Myeloids/result.h5ad')
# outdir = '/mnt/NexData_file/home/xzx/project/AE2508064183/shouhou5'

# cellnames = ['Myeloids_filter']
# for cells in cellnames:
#     print(cells)
#     celltype = cells
#     outdir1 = os.path.join(outdir,celltype)
#     os.makedirs(outdir1,exist_ok=True)
#     #adata_tmp = adata[~adata.obs['leiden'].isin(['9'])].copy()

#     adata_tmp = adata[~adata.obs['leiden'].isin(['10'])].copy()
#     #sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
#     adata_tmp = recluster(adata_tmp,species='human',resolution = 0.5,n_pcs = 8,batch_effect=2,output=outdir1,neighbors = 15)
#     adata_tmp = adata_tmp.sub_cluster()