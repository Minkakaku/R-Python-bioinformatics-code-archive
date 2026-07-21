# import diopy
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
# os.chdir('/lvdata/wzb/scRNA/FW2023-656/2024.01.02/')
import pandas as pd
adata=sc.read_h5ad('/mnt/NexData_file/home/wmx/Project/AE2507033447-08/20260413/subcluster/B_cells_and_Plasma/result.h5ad')
outdir = '/mnt/NexData_file/home/wmx/Project/AE2507033447-08/20260413/subcluster/'

celltype = 'B_and_Plasma_delete'
outdir1 = os.path.join(outdir,celltype)
os.makedirs(outdir1,exist_ok=True)
adata_tmp = adata[~adata.obs['leiden'].isin(['6','8'])].copy()
sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
adata_tmp = recluster(adata_tmp,species='human',resolution = 0.8,n_pcs = 15,batch_effect=2,output=outdir1) # 调整一下分辨率
adata_tmp = adata_tmp.sub_cluster()
