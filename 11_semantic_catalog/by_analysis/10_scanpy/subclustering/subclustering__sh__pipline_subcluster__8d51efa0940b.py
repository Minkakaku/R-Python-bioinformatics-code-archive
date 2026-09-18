import os, sys
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

# 映射细分，重新聚类
adata_sub = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2507163761-01/02.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2507163761-01/03.Mission/Main_markerplot/"

from matplotlib.colors import LinearSegmentedColormap
# top_marker=["Slc12a3","Gramd1b","Pecam1","Krt19","Lrp1","Atp6v1g3","Ctss","Akap12","Tns3","Ptpro","Phactr1","Slc17a3","Slc12a1"]
cmap = LinearSegmentedColormap.from_list('cut_gray_red', ['#cccccc', '#cccccc', '#ff0000'], N=256)
# for i in top_marker:
#         sc.pl.umap(adata_sub, color=i, title=i, show=False, cmap=cmap, vmin=None)
#         plt.savefig(os.path.join(out_dir, f'umap/{i}_umap.png'),bbox_inches='tight')
#         plt.savefig(os.path.join(out_dir, f'umap/{i}_umap.pdf'),bbox_inches='tight')

select_marker=["Slc17a3","Cubn","Cda","Slc12a3","Trpm6","Abca13","Pecam1","Eng","Kdr","Krt19","Krt18","Krt7","Lrp1","Pdgfrb","Col1a2","Atp6v1g3","Atp6v0d2","Clcnkb","Ctss","Ms4a6c","Adgre1","Ptpro","Podxl","Nphs2","Phactr1","Scnn1b","Hsd11b2","Slc12a1",'Fxyd2','Kcnj1']
sc.pl.dotplot(adata_sub, var_names=select_marker, groupby='cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type_scale.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata_sub, var_names=select_marker, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type.pdf'),bbox_inches='tight')


sc.pl.umap(adata_sub, color='Slc17a3', title='Slc17a3', show=False, cmap='Reds', vmin=None)
plt.savefig(os.path.join(out_dir, 'umap/Slc17a3_umap.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap/Slc17a3_umap.pdf'),bbox_inches='tight')
sc.pl.umap(adata_sub, color='Slc12a3', title='Slc12a3', show=False, cmap='Reds', vmin=None)
plt.savefig(os.path.join(out_dir, 'umap/Slc12a3_umap.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap/Slc12a3_umap.pdf'),bbox_inches='tight')
sc.pl.umap(adata_sub, color='Pecam1', title='Pecam1', show=False, cmap='Reds', vmin=None)
plt.savefig(os.path.join(out_dir, 'umap/Pecam1_umap.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap/Pecam1_umap.pdf'),bbox_inches='tight')
sc.pl.umap(adata_sub, color='Krt19', title='Krt19', show=False, cmap='Reds', vmin=None)
plt.savefig(os.path.join(out_dir, 'umap/Krt19_umap.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap/Krt19_umap.pdf'),bbox_inches='tight')
sc.pl.umap(adata_sub, color='Lrp1', title='Lrp1', show=False, cmap='Reds', vmin=None)
plt.savefig(os.path.join(out_dir, 'umap/Lrp1_umap.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap/Lrp1_umap.pdf'),bbox_inches='tight')
sc.pl.umap(adata_sub, color='Atp6v1g3', title='Atp6v1g3', show=False, cmap='Reds', vmin=None)
plt.savefig(os.path.join(out_dir, 'umap/Atp6v1g3_umap.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap/Atp6v1g3_umap.pdf'),bbox_inches='tight')
sc.pl.umap(adata_sub, color='Ctss', title='Ctss', show=False, cmap='Reds', vmin=None)
plt.savefig(os.path.join(out_dir, 'umap/Ctss_umap.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap/Ctss_umap.pdf'),bbox_inches='tight')
sc.pl.umap(adata_sub, color='Ptpro', title='Ptpro', show=False, cmap='Reds', vmin=None)
plt.savefig(os.path.join(out_dir, 'umap/Ptpro_umap.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap/Ptpro_umap.pdf'),bbox_inches='tight')
sc.pl.umap(adata_sub, color='Phactr1', title='Phactr1', show=False, cmap='Reds', vmin=None)
plt.savefig(os.path.join(out_dir, 'umap/Phactr1_umap.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap/Phactr1_umap.pdf'),bbox_inches='tight')
sc.pl.umap(adata_sub, color='Slc12a1', title='Slc12a1', show=False, cmap='Reds', vmin=None)
plt.savefig(os.path.join(out_dir, 'umap/Slc12a1_umap.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap/Slc12a1_umap.pdf'),bbox_inches='tight')