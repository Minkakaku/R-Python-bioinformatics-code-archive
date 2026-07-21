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
sys.path.append('/mnt/NexData_file/home/wmx/SC2601280737/SC2601280737_shouhou1/scanpy_fun')
from subcluster import recluster
from find_marker import find_marker_gene
from draw_group import draw_group
from draw_bar import plot_celltype_proportion
# os.chdir('/lvdata/wzb/scRNA/FW2023-656/2024.01.02/') # 改变工作路径
import pandas as pd
# adata=sc.read_h5ad('/mnt/NexData_file/home/wmx/SC2601280737/06-Cell_type_annotation/result.h5ad')
outdir = '/mnt/NexData_file/home/wmx/SC2601280737/SC2601280737_shouhou1/subcluster/'

target_celltypes = ['T cells', 'Epithelial'] # 售后目前只需要做上皮 和 T细胞
# # 确保 outdir 存在
# os.makedirs(outdir, exist_ok=True)
# for celltype in target_celltypes:
#     outdir1 = os.path.join(outdir,celltype)
#     os.makedirs(outdir1,exist_ok=True)
#     adata_tmp = adata[adata.obs['cell_type'] == celltype].copy()
#     sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
#     adata_tmp = recluster(adata_tmp,species='mouse',resolution = 0.5,n_pcs = 20,batch_effect=0,output=outdir1)
#     adata_tmp = adata_tmp.sub_cluster()
    
adata=sc.read_h5ad('/mnt/NexData_file/home/wmx/SC2601280737/SC2601280737_shouhou1/subcluster/Epithelial/result.h5ad')
target_celltypes = ['Epithelial_delete'] # 去掉上皮中的6亚群，是质量差的亚群
# 确保 outdir 存在
os.makedirs(outdir, exist_ok=True)
for celltype in target_celltypes:
    outdir1 = os.path.join(outdir,celltype)
    os.makedirs(outdir1,exist_ok=True)
    adata_tmp = adata[adata.obs['leiden'] != '6'].copy()
    sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
    adata_tmp = recluster(adata_tmp,species='mouse',resolution = 0.5,n_pcs = 20,batch_effect=0,output=outdir1)
    adata_tmp = adata_tmp.sub_cluster()

