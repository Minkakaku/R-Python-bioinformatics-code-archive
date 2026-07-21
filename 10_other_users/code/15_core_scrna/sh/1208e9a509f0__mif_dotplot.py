import os, sys
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["OMP_NUM_THREADS"] = "16"
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
import pandas as pd
from matplotlib import font_manager as fm

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    '''
    out_dir : 输出文件夹地址
    file_name : 文件名称
    fig_size : 输出图片的大小
    '''
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()
    
adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/mapping_raw/Main_mapping.h5ad")
adata.obs['sample'].cat.categories
adata_sub = adata[adata.obs['sample'].isin(['iMCD-PB_3','iMCD-PB_7','iMCD-PB_17','iMCD-PB_20','iMCD-PB_21','iMCD-PB_5']), :]
adata_sub.obs['sample'].cat.categories

group_dict={
    'iMCD-PB_3':'3-7-20-21',
    'iMCD-PB_5':'5-17',
    'iMCD-PB_7':'3-7-20-21',
    'iMCD-PB_17':'5-17',
    'iMCD-PB_20':'3-7-20-21',
    'iMCD-PB_21':'3-7-20-21'
}
adata_sub.obs['group'] = [group_dict[samples] for samples in adata_sub.obs['sample']]
adata_sub.obs['group'] = adata_sub.obs['group'].astype("category")
adata_sub.obs['sample'] = adata_sub.obs['sample'].cat.reorder_categories(['iMCD-PB_3','iMCD-PB_7','iMCD-PB_17','iMCD-PB_20','iMCD-PB_21','iMCD-PB_5'], ordered=True)

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/iMCD_plot/"
sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['sample'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_sample.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_sample.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.violin(adata_sub, ['MIF'], jitter=0, groupby='sample', rotation=45, multi_panel=True, show=False)
plt.savefig(os.path.join(out_dir, 'violin_sample.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'violin_sample.pdf'), facecolor='white', bbox_inches="tight")

sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_group.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_group.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.violin(adata_sub, ['MIF'], jitter=0, groupby='group', rotation=45, multi_panel=True, show=False)
plt.savefig(os.path.join(out_dir, 'violin_group.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'violin_group.pdf'), facecolor='white', bbox_inches="tight")

adata=sc.read("/mnt/NexData_file/scRNA/scVision_output/iMCD_after/result.h5ad")
sc.pl.dotplot(adata, var_names=['MIF'], groupby=['sample'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_after_sample.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_after_sample.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.violin(adata, ['MIF'], jitter=0, groupby='sample', rotation=45, multi_panel=True, show=False)
plt.savefig(os.path.join(out_dir, 'violin_after_sample.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'violin_after_sample.pdf'), facecolor='white', bbox_inches="tight")