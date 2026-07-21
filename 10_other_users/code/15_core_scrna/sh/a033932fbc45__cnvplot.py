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

adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251249/04.Mission/mapping_raw/Main_mapping.h5ad")
csv_df = pd.read_csv('/mnt/NexData_file/home/sh/Projects/FWAE20251249/05.Mission/CNV/metadata.csv')
csv_df = csv_df.drop(columns=set(csv_df.columns) & set(adata.obs.columns))
csv_df = csv_df.rename(columns={'Unnamed: 0': 'barcode'})
csv_df = csv_df.set_index('barcode')
adata.obs = adata.obs.join(csv_df, how='left')

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251249/05.Mission/CNV/"
adata_sub=adata[adata.obs.cell_type=="Tumor-like cells", :]
sc.pl.dotplot(adata_sub, var_names=['cnv_score_cell'], groupby=['group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_cnv_score_cell-group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_cnv_score_cell-group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata_sub, var_names=['cnv_score_cell'], groupby=['sample'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_cnv_score_cell-sample.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_cnv_score_cell-sample.pdf'),bbox_inches='tight')