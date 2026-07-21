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


# 加载亚群细分的函数
import sys
sys.path.append('/mnt/NexData_file/share/scRNA-tools/scanpy_fun')
from subcluster import recluster
from find_marker import find_marker_gene
from draw_group import draw_group
from draw_bar import plot_celltype_proportion
os.chdir('/mnt/NexData_file/home/wmx/Project/SC2411133007')
import pandas as pd

# outdir = '/mnt/NexData_file/home/wmx/Project/SC2411133007/subcluster/'
# adata=sc.read_h5ad('/mnt/NexData_file/home/wmx/Project/SC2411133007/subcluster/T cells/result.h5ad')
# celltype = 'T_delete'
# outdir1 = os.path.join(outdir,celltype)
# os.makedirs(outdir1,exist_ok=True)
# adata_tmp = adata[~adata.obs['leiden'].isin(['3','7'])].copy()
# sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
# adata_tmp = recluster(adata_tmp,species='human',resolution = 0.5,n_pcs = 10,batch_effect=2,neighbors = 15,output=outdir1)
# adata_tmp = adata_tmp.sub_cluster()

# outdir = '/mnt/NexData_file/home/wmx/Project/SC2411133007/subcluster/'
# adata=sc.read_h5ad('/mnt/NexData_file/home/wmx/Project/SC2411133007/subcluster/T cells/result.h5ad')
# celltype = 'T_delete_only_c7'
# outdir1 = os.path.join(outdir,celltype)
# os.makedirs(outdir1,exist_ok=True)
# adata_tmp = adata[~adata.obs['leiden'].isin(['7'])].copy()
# sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
# adata_tmp = recluster(adata_tmp,species='human',resolution = 0.5,n_pcs = 10,batch_effect=2,neighbors = 15,output=outdir1)
# adata_tmp = adata_tmp.sub_cluster()

outdir = '/mnt/NexData_file/home/wmx/Project/SC2411133007/subcluster/'
adata=sc.read_h5ad('/mnt/NexData_file/home/wmx/Project/SC2411133007/subcluster/T_delete/result.h5ad')
celltype = 'T_delete_res_1.2'
outdir1 = os.path.join(outdir,celltype)
os.makedirs(outdir1,exist_ok=True)
sc.pl.umap(adata,color = ['cell_type','sample'])
adata = recluster(adata,species='human',resolution = 1.2,n_pcs = 8,batch_effect=2,neighbors = 15,output=outdir1)
adata = adata.sub_cluster()



