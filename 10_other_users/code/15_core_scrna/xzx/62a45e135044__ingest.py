import os
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
import pandas as pd
import scanpy as sc
import anndata as ad
import scanpy.external as sce
import numpy as np
import os
import shutil
import random
import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
import streamlit as st
from PIL import Image
import celltypist
from celltypist import models
import scanpy.external as sce
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from matplotlib import font_manager as fm
import warnings
# adata_obj = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/AE2509104884/03-Cell_filter/adata_counts.h5ad")
# adata_ref = sc.read_h5ad("/mnt/NexData_file/home/xzx/project/result.h5ad")

# var_names = adata_ref.var_names.intersection(adata_obj.var_names)
# adata_ref = adata_ref[:, var_names]
# adata_obj = adata_obj[:, var_names]
# sc.tl.ingest(adata_obj, adata_ref, obs="cell_type")
# adata_obj.uns["cell_type_colors"] = adata_ref.uns["cell_type_colors"] 

# sc.pp.normalize_total(adata_obj)
# sc.pp.log1p(adata_obj)
# sc.pp.scale(adata_obj)

# sc.tl.pca(adata_obj, svd_solver='arpack')
# #sce.pp.bbknn(adata_obj, batch_key='sample', n_pcs=15)
# #sc.external.pp.bbknn(adata_obj, batch_key="sample") 
# sce.pp.harmony_integrate(adata_obj, key ='sample',adjusted_basis="rep")
# sc.pp.neighbors(adata_obj, n_neighbors=15, n_pcs=20)
# sc.tl.umap(adata_obj)

# output = '/mnt/NexData_file/home/xzx/project/AE2509104884'
# sc.pl.umap(adata_obj, color=['cell_type','sample'], show=False)
# save_fig(plt, output, file_name='cluster_umap')
# sc.pl.umap(adata_obj, color=['cell_type','sample'],legend_loc ='on data', show=False)
# save_fig(plt, output, file_name='cluster_umap_withlabel')

# adata_obj.write_h5ad('/mnt/NexData_file/home/xzx/project/AE2509104884/adata_anno.h5ad')


#sce.pp.harmony_integrate(adata_concat, key ='sample',adjusted_basis="X_pca")
# sc.pp.neighbors(adata_concat, n_neighbors=20, n_pcs=20)
# sce.pp.bbknn(adata_concat, batch_key='sample', n_pcs=20)
# sc.tl.umap(adata_concat)
# sc.pl.umap(adata_concat, color=['cell_type','sample'],legend_loc ='on data', show=False)
# output = '/mnt/NexData_file/home/xzx/project/AE2509104884'
# save_fig(plt, output, file_name='cluster_umap_withlabel')
# sc.pl.umap(adata_concat, color=['cell_type','sample'], show=False)
# save_fig(plt, output, file_name='cluster_umap')

adata_count_filtered = sc.read_h5ad('/mnt/NexData_file/home/xzx/project/AE2509104884/adata_count_filtered.h5ad')

adata_count_filtered.layers["counts"] = adata_count_filtered.X.copy()
# data normalization and log10 transform
sc.pp.normalize_total(adata_count_filtered)
sc.pp.log1p(adata_count_filtered)
adata_count_filtered.raw = adata_count_filtered
# select highly variable gene and data scale 
sc.pp.highly_variable_genes(adata_count_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=2000)
sc.pl.highly_variable_genes(adata_count_filtered, show=False)
#save_fig(plt, out_dir, file_name='highly_variable_genes', fig_size={'width':16, 'height':8})
adata_count_filtered = adata_count_filtered[:, adata_count_filtered.var['highly_variable']]
sc.pp.regress_out(adata_count_filtered, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata_count_filtered, max_value=10)
# PCA
sc.tl.pca(adata_count_filtered, svd_solver='arpack')
#sc.pl.pca(adata_count_filtered, color='sample', show=False)
#save_fig(plt, out_dir, file_name='pca_sample')
#sc.pl.pca(_adata, color='group', show=False)
#save_fig(plt, out_dir, file_name='pca_group')
#sc.pl.pca_variance_ratio(_adata, log=True, show=False)
#save_fig(plt, out_dir, file_name='pca_variance_ratio')

#if batch_effect == '**harmony**':
#    sce.pp.harmony_integrate(_adata, 'sample', adjusted_basis='X_pca')
# cluster

sce.pp.harmony_integrate(adata_count_filtered, 'sample', adjusted_basis='X_pca')
sc.pp.neighbors(adata_count_filtered, n_neighbors=20, n_pcs=20)    

# Remove batch effect
#if batch_effect == '**bbknn**':
#sce.pp.bbknn(adata_count_filtered, batch_key='sample', n_pcs=20)  # running bbknn 1.3.6
sc.tl.umap(adata_count_filtered)


out_dir = '/mnt/NexData_file/home/xzx/project/AE2509104884'
adata_count_filtered.write_h5ad(os.path.join(out_dir, 'result_filtered.h5ad'))
with rc_context({'figure.figsize': (10, 10)}):
    sc.pl.umap(adata_count_filtered, 
                color=['cell_type', "sample"],
                title=['cell_type', "sample"],
                legend_loc='on data',
                frameon=True,
                legend_fontsize=9, 
                legend_fontoutline=3, 
                wspace=0.2, 
                ncols=2, 
                #cmap='coolwarm', 
                use_raw=True, 
                show=False)
            
    plt.savefig(os.path.join(out_dir, 'umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
