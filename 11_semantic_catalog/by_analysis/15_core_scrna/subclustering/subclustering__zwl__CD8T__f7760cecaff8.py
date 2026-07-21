import pandas as pd
import anndata as ad
import numpy as np
import gc
import sys
gc.collect()
# 保存图片
import os
import matplotlib.pyplot as plt
import scanpy as sc
sc.settings.set_figure_params(dpi=80, facecolor="white")
# sys.path.append('/mnt/NexData_file/share/scRNA-tools/scanpy_fun')
sys.path.append('/mnt/NexData_file/home/zwl/Project/personal')
from subcluster import *
from draw_group import draw_group
# sys.path.append('/mnt/NexData_file/home/zwl/Project/personal')
from commen_def import *
# from find_marker import find_marker_gene



os.getcwd()
os.chdir('/mnt/NexData_file/home/zwl/Project/FWSC20240999')
adata_CD8 = sc.read_h5ad('/mnt/NexData_file/share/lst/FWSC20240999/CD8T.h5ad')

#'T cells', 'Epithelial cells', 'Macrophages', 'Endothelial cells','Fibroblasts','B cells'
i='CD8T'
outdir = '/mnt/NexData_file/home/zwl/Project/FWSC20240999/no_bbknn'
outdir = os.path.join(outdir,i+'_bbknn')
os.makedirs(outdir,exist_ok=True)


# 过滤细胞数过少的批次
batch_counts = adata_CD8.obs['sample'].value_counts()
keep_batches = batch_counts[batch_counts >= 2].index.tolist()
# 过滤批次
adata_i_filtered = adata_CD8[adata_CD8.obs['sample'].isin(keep_batches), :]
sc.pl.umap(adata_i_filtered,color = ['cell_type','sample'])

adata_i_filtered = recluster(adata_i_filtered,species='human',resolution = 0.8,batch_effect = 1,output=outdir)
adata_i_filtered = adata_i_filtered.sub_cluster()


