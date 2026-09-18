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

adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12_mod.Mission/mapping_raw/Main_mapping.h5ad")
adata.obs['sample'].cat.categories

adata_sub = adata[adata.obs['sample'].isin(['iMCD_3','iMCD_7','iMCD_17','iMCD_20','iMCD_21','iMCD_5']), :]
adata_sub.obs['sample'].cat.categories
adata_sub.obs['sample'] = pd.Categorical(
    adata_sub.obs['sample'],                       # 原值保留
    categories=['iMCD_3','iMCD_7','iMCD_17','iMCD_20','iMCD_21','iMCD_5'],                     # 新顺序
    ordered=False                                 # 无需数学大小比较就 False
)

group_dict={
    'iMCD_3':'severe',
    'iMCD_5':'non-severe',
    'iMCD_7':'severe',
    'iMCD_17':'non-severe',
    'iMCD_20':'severe',
    'iMCD_21':'severe'
}
adata_sub.obs['group'] = [group_dict[samples] for samples in adata_sub.obs['sample']]
adata_sub.obs['group'] = adata_sub.obs['group'].astype("category")

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/17.Mission/fig3_F_severe_plot/"
sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['group'], dendrogram=False, show=False, dot_max=0.8)
plt.savefig(os.path.join(out_dir, 'group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)



bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/proportion/Main_mapping.h5ad")
bdata.obs['sample'].cat.categories
bdata_sub = bdata[bdata.obs['sample'].isin(['iMCD-PB_1','iMCD-PB_2','iMCD-PB_8','iMCD-PB_11','iMCD-PB_15','iMCD-PB_16']), :]

adata =  adata_sub.concatenate(bdata_sub)
adata.obs['sample'].cat.categories
group_dict={
    'iMCD_3':'severe',
    'iMCD_5':'non-severe',
    'iMCD_7':'severe',
    'iMCD_17':'non-severe',
    'iMCD_20':'severe',
    'iMCD_21':'severe',
    'iMCD-PB_1':'non-severe',
    'iMCD-PB_2':'non-severe',
    'iMCD-PB_8':'severe',
    'iMCD-PB_11':'non-severe',
    'iMCD-PB_15':'non-severe',
    'iMCD-PB_16':'non-severe'
}
adata.obs['group'] = [group_dict[samples] for samples in adata.obs['sample']]
adata.obs['group'] = adata.obs['group'].astype("category")

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/17.Mission/fig3_F_severe_plot/"
sc.pl.dotplot(adata, var_names=['MIF'], groupby=['group'], dendrogram=False, show=False, dot_max=0.8)
plt.savefig(os.path.join(out_dir, 'group_bubble_new.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'group_bubble_new.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_sub=adata[adata.obs['cell_type'].isin(['Plasma cells']), :]
sc.pl.dotplot(adata_sub, var_names=['MIF'], groupby=['group'], dendrogram=False, show=False, dot_max=0.8,expression_cutoff=0.5)
plt.savefig(os.path.join(out_dir, 'group_bubble_new_Plasma.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'group_bubble_new_Plasma.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.dotplot(adata, var_names=['MIF'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble_new.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble_new.pdf'), facecolor='white', bbox_inches="tight", dpi=300)