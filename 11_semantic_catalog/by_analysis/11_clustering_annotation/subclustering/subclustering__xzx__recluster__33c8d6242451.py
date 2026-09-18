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

adata=sc.read_h5ad('/mnt/NexData_file/scRNA/scVision_output/SC2512117402/06-Cell_type_annotation/result.h5ad')
outdir = '/mnt/NexData_file/home/xzx/project/SC2512117402/shouhou2'
os.makedirs(outdir,exist_ok=True)



outdir1 = os.path.join(outdir,'Epi')
os.makedirs(outdir1,exist_ok=True)

adata_tmp = adata[adata.obs['cell_type'].isin(['Epithelial cells','Proximal tubule cells','Podocytes'])].copy()
#sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
adata_tmp = recluster(adata_tmp,species='human',resolution = 0.5,n_pcs = 10,neighbors = 15,batch_effect=2,output=outdir1)
adata_tmp = adata_tmp.sub_cluster()