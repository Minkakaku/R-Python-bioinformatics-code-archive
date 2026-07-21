# 载入包
import scanpy as sc
import numpy as np
import pandas as pd
import os
from matplotlib import pyplot as plt
import scipy
from scanpy_fun_subcluster import recluster

# recluster
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20250297/result.h5ad")
#adata_raw = adata_raw[adata_raw.obs.cell_type=="T cells", :]
adata_sub = recluster(adata=adata_raw,species='mouse',neighbors=20,resolution=1,n_pcs = 15,output="/data/NAS03/yrc/FWSC20250297/")
adata = adata_sub.sub_cluster()
adata = adata.raw.to_adata()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=5,
    show=True
)
plt.savefig("/data/NAS03/yrc/FWSC20250297/top5_diff_genes_dotplot.png", bbox_inches="tight", dpi=300)
plt.savefig("/data/NAS03/yrc/FWSC20250297/top5_diff_genes_dotplot.pdf", bbox_inches="tight", dpi=300)

# marker鉴定
markergenes = ['Cd3d','Cd8a','Cd4','Foxp3','Nkg7','Klrd1']
#vilnplot
sc.pl.violin(adata, markergenes, groupby='leiden')
plt.savefig('/data/NAS03/yrc/FWSC20250297/marker_vlnplot.png',bbox_inches = 'tight')
plt.savefig('/data/NAS03/yrc/FWSC20250297/marker_vlnplot.pdf',bbox_inches = 'tight')
#feature plot
sc.pl.umap(adata, color=markergenes, ncols=3,cmap='coolwarm')
plt.savefig(f'/data/NAS03/yrc/FWSC20250297/marker_featureplot.png',bbox_inches = 'tight')
plt.savefig(f'/data/NAS03/yrc/FWSC20250297/marker_featureplot.pdf',bbox_inches = 'tight')
#dotplot
sc.pl.dotplot(adata, markergenes,groupby='leiden',show=False)
plt.savefig(f'/data/NAS03/yrc/FWSC20250297/marker_dotplot.png',bbox_inches = 'tight')
plt.savefig(f'/data/NAS03/yrc/FWSC20250297/marker_dotplot.pdf',bbox_inches = 'tight')

# 假设我们想要找到cluster 0相对于其他所有clusters中的差异表达基因
sc.tl.rank_genes_groups(adata, groupby='leiden', groups=['0'], reference='rest', method='wilcoxon')


# 注释
annotation = {
    '0': 'CD4 T',
    '1': 'CD14+ Mono',
    '2': 'B',
    '3': 'CD8 T',
    '4': 'FCGR3A+ Mono',
    '5': 'NK',
    '6': 'DC'
}
adata.obs['celltype'] = adata.obs['leiden'].map(annotation).astype('category')

sc.pl.umap(adata, color=['leiden'], legend_loc='on data',title='celltype', add_outline=False,
           frameon=False, legend_fontsize=10, legend_fontoutline=2)
sc.pl.umap(adata, color='celltype', legend_loc='on data',title='celltype', add_outline=True,
           frameon=False, legend_fontsize=10, legend_fontoutline=2)
plt.savefig(f'/data/NAS03/yrc/FWSC20250297/umap_celltype.png',bbox_inches = 'tight', dpi=300)
plt.savefig(f'/data/NAS03/yrc/FWSC20250297/umap_celltype.pdf',bbox_inches = 'tight', dpi=300)