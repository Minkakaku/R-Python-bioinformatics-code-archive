import os
import cupy as cp
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import rapids_singlecell as rsc
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

import warnings
warnings.filterwarnings("ignore")

plt.style.use('default')  # 使用默认样式
plt.rcParams['figure.facecolor'] = 'white'  # 确保图像背景为白色
sc.settings.set_figure_params(dpi=100, dpi_save=200, figsize=(5, 5), facecolor='white')

import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(
    managed_memory=True,  # Allows oversubscription
    pool_allocator=False, # default is False
    devices=0,            # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)

def cluster_gpu(adata, n_hvg=2000, n_neighbors=20, n_pcs=15, res=[0.5]):
    rsc.get.anndata_to_GPU(adata, convert_all=True)

    # Normalize
    rsc.pp.normalize_total(adata, inplace=True)
    rsc.pp.log1p(adata)
    adata.raw = adata.copy()

    # HVG
    rsc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_hvg,
        flavor="seurat",
        batch_key="sample",
    )

    adata = adata[:, adata.var['highly_variable']]
    rsc.pp.regress_out(adata, ['total_counts'])
    rsc.pp.scale(adata, max_value=10)
    rsc.pp.pca(adata)
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=30)


    rsc.pp.harmony_integrate(adata, key="sample")

    # Cluster
    rsc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, metric='cosine', use_rep='X_pca_harmony', random_state=10)
    rsc.tl.umap(adata, min_dist=0.3, spread=1, random_state=10)
    for r in res:
        rsc.tl.leiden(adata, resolution=r, key_added=f'res_{r}', random_state=10)

    rsc.get.anndata_to_CPU(adata, convert_all=True)

    if 120000/adata.obs.shape[0] < 1:
        size = (120000/adata.obs.shape[0]) * 5
    else:
        size = 120000/adata.obs.shape[0]
        
    with rc_context({'figure.figsize': (6, 6)}):
        sc.pl.umap(
            adata, 
            color=['sample'] + [f'res_{r}' for r in res],
            legend_loc='on data', 
            legend_fontsize=9, 
            legend_fontoutline=3,
            frameon=True,
            size=size,
            show=False
        )

    return adata

#os.chdir('/mnt/NexData_file/home/hy/Project/Project/liuxiaoxiao/xenium/output')

res = [0.2,0.4,0.5,0.6,0.8,1]
adata = sc.read_h5ad('/mnt/NexData_file/home/hy/Project/Project/liuxiaoxiao/xenium/output/data/adata_filtered.h5ad')
adata = cluster_gpu(adata, n_hvg=2000, n_neighbors=20, n_pcs=10, res=res)

adata.write(os.path.join('/mnt/NexData_file/home/hy/Project/Project/liuxiaoxiao/xenium/output/data/adata_cluster_0925_new.h5ad'))
