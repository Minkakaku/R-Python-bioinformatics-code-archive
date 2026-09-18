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
# adata=sc.read_h5ad('/mnt/NexData_file/home/xzx/project/SC2509164991/shouhou1/T cells/result.h5ad')
outdir = '/mnt/NexData_file/home/wmx/SC2509164991/subcluster/'
#for celltype in list(set(adata.obs['cell_type'])): 

# celltype = 'T_delete'
# outdir1 = os.path.join(outdir,celltype)
# os.makedirs(outdir1,exist_ok=True)
# adata_tmp = adata[~adata.obs['leiden'].isin(['0', '3','6','7', '8'])].copy()
# sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
# adata_tmp = recluster(adata_tmp,species='human',resolution = 0.5,n_pcs = 20,batch_effect=0,output=outdir1)
# adata_tmp = adata_tmp.sub_cluster()

# adata=sc.read_h5ad('/mnt/NexData_file/scRNA/scVision_output/SC2509164991/06-Cell_type_annotation/result.h5ad')
# celltype = 'B_and_plasma'
# outdir1 = os.path.join(outdir,celltype)
# os.makedirs(outdir1,exist_ok=True)
# adata_tmp = adata[adata.obs['cell_type'].isin(['B cells','Plasma cells'])].copy()
# sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
# adata_tmp = recluster(adata_tmp,species='human',resolution = 0.5,n_pcs = 15,batch_effect=2,output=outdir1)
# adata_tmp = adata_tmp.sub_cluster()

# adata=sc.read_h5ad('/mnt/NexData_file/home/wmx/SC2509164991/subcluster/B_and_plasma/result.h5ad')
# celltype = 'B_and_Plasma_delete'
# outdir1 = os.path.join(outdir,celltype)
# os.makedirs(outdir1,exist_ok=True)
# adata_tmp = adata[~adata.obs['leiden'].isin(['6','7','9','10'])].copy()
# sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
# adata_tmp = recluster(adata_tmp,species='human',resolution = 0.5,n_pcs = 10,batch_effect=2,neighbors=10,output=outdir1) # 调分辨率，B全是一群，看能不能分开
# adata_tmp = adata_tmp.sub_cluster()

# Pericyte
adata=sc.read_h5ad('/mnt/NexData_file/scRNA/scVision_output/SC2509164991/06-Cell_type_annotation/result.h5ad')
celltype = 'Pericyte'
outdir1 = os.path.join(outdir,celltype)
os.makedirs(outdir1,exist_ok=True)
adata_tmp = adata[~adata.obs['cell_type'].isin(['Pericyte'])].copy()
sc.pl.umap(adata_tmp,color = ['cell_type','sample'])
adata_tmp = recluster(adata_tmp,species='human',resolution = 0.5,n_pcs = 20,batch_effect=2,neighbors=20,output=outdir1) 
adata_tmp = adata_tmp.sub_cluster()