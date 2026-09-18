#export OPENBLAS_NUM_THREADS=1
#export MKL_NUM_THREADS=1
#export OMP_NUM_THREADS=1
#export NUMEXPR_NUM_THREADS=1

import os
# 强制所有计算库使用单线程
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

#print(adata.obs['cell_type'].unique())

# outdir='/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou1/'

# labels = adata.obs['cell_type'].unique()
# for i in labels:
#     outdir1= os.path.join(outdir,i)
#     os.makedirs(outdir1, exist_ok=True)
#     adata_tmp=adata[adata.obs['cell_type']==i]
#     #sc.pl.umap(adata_tmp,color = ['cell_type','group'])
#     adata_tmp = recluster(adata_tmp,species='human',neighbors=10,resolution=0.3,n_pcs = 10,batch_effect=2,output=outdir1)
#     adata_tmp = adata_tmp.sub_cluster()

# adata = sc.read_h5ad("/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou1/adata_pericyte.h5ad")
# outdir = "/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou1/stromal_all"
# os.makedirs(outdir,exist_ok=True)
# adata = recluster(adata,species='human',resolution=0.7,n_pcs = 8,batch_effect=1,neighbors=10,output=outdir)
# adata = adata.sub_cluster()

# adata = sc.read_h5ad("/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou1/adata_mye.h5ad")
# print(adata.obs['cell_type'].unique())
# outdir = "/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou1/Myeloids"
# os.makedirs(outdir,exist_ok=True)
# adata = recluster(adata,species='human',resolution=0.9,n_pcs = 8,batch_effect=1,output=outdir,neighbors=10)
# adata = adata.sub_cluster()

# adata = sc.read_h5ad("/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou1/adata_stromal.h5ad")
# outdir = "/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou1/stromal"
# os.makedirs(outdir,exist_ok=True)
# adata = recluster(adata,species='human',resolution=0.5,n_pcs = 8,batch_effect=1,neighbors=10,output=outdir)
# adata = adata.sub_cluster()

# adata = sc.read_h5ad("/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou1/adata_pericytes_only.h5ad")
# print(adata.obs['cell_type'].unique())
# outdir = "/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou1/pericytes"
# os.makedirs(outdir,exist_ok=True)
# adata = recluster(adata,species='human',resolution=0.5,n_pcs = 8,batch_effect=1,output=outdir,neighbors=10)
# adata = adata.sub_cluster()



# adata = sc.read_h5ad("/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou2/Myeloids/adata_mye_filter.h5ad")
# print(adata.obs['cell_type'].unique())
# outdir = "/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou2/Myeloids/recluster"
# os.makedirs(outdir,exist_ok=True)
# adata = recluster(adata,species='human',resolution=0.7,n_pcs = 8,batch_effect=1,output=outdir,neighbors=10)
# adata = adata.sub_cluster()

# adata = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2509195142/06-Cell_type_annotation/result.h5ad")
# #print(adata.obs['cell_type'].unique())
# outdir = "/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou2/Fibroblasts/recluster"
# os.makedirs(outdir,exist_ok=True)
# adata_tmp = adata[adata.obs['cell_type'].isin(['Fibroblasts'])]
# adata_tmp = recluster(adata_tmp,species='human',resolution=0.5,n_pcs = 8,batch_effect=1,output=outdir,neighbors=10)
# adata_tmp = adata_tmp.sub_cluster()

adata = sc.read_h5ad("/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou3/adata_stromal.h5ad")
print(adata.obs['cell_type'].unique())
outdir = "/mnt/NexData_file/home/xzx/project/AE2509195142/shouhou3/stromal/recluster"
os.makedirs(outdir,exist_ok=True)
adata = recluster(adata,species='human',resolution=0.7,n_pcs = 10,batch_effect=1,output=outdir,neighbors=10)
adata = adata.sub_cluster()
