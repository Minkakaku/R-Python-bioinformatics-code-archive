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
import seaborn as sns

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
    plt.savefig(os.path.join(out_dir, f'{file_name}.tiff'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def plot_celltype_proportion(adata, out_dir,type='sample', fig_size=None,prefix='F'):
    keys = adata.obs['cell_type'].cat.categories.to_list()
    values = list(adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))
    ax = pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').plot(
        kind='bar', xlabel='', stacked=True,
        color=[cell_type_colors[cell] for cell in pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').columns],
        edgecolor='black', width=0.8
    )
    ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')
    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=10, title="Cell Types", 
        title_fontsize=12, prop=legend_font, frameon=False, handletextpad=1.5, labelspacing=1.0, 
        columnspacing=1.0, markerscale=1, alignment='left', handlelength=2, handleheight=2
    )
    ax.tick_params(axis='x', labelsize=12, labelrotation=90, labelcolor='black')
    ax.tick_params(axis='y', labelsize=12, labelcolor='black')
    ax.grid(False)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    if fig_size is None:
        fig_size = {'width': len(set(adata.obs[type])) * 2, 'height': 6}
    plt.tight_layout()
    file_name=f"{prefix}_celltype_proportion_{type}"
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/fig1_A_umap/"

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/proportion/Main_mapping.h5ad")
adata_raw = adata_raw[~adata_raw.obs.cell_type.isin(['Erythroid cells','Neutrophils']), :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,wspace=0.1, legend_fontsize=12)
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_celltype.pdf'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_celltype.tiff'))
# sc.pl.umap(adata_raw, color="cell_type", title='Cell Type',size=8, wspace=0.1, legend_fontsize=9, legend_fontweight='normal', legend_fontoutline=2, show=False)
# save_fig(plt, out_dir, file_name='Main_dimplot_umap_celltype', fig_size={'width':5, 'height':5})
# 'B cells', 'Erythroid cells', 'Mast cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'

adata_raw.uns['group_colors']=["#F8766D", "#00BF7D", "#00B0F6"]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='group', title='group', show=False)
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_group.pdf'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_group.tiff'))

desired_order=["iMCD-PB_1","iMCD-PB_2","iMCD-PB_8","iMCD-PB_11","iMCD-PB_15","iMCD-PB_16","iMCD-PB_3","iMCD-PB_5","iMCD-PB_7","iMCD-PB_17","iMCD-PB_20","iMCD-PB_21","iMCD-PB_6","iMCD-PB_18","iMCD-PB_25"]
adata_raw.obs['sample'] = (
    pd.Categorical(adata_raw.obs['sample'],
                   categories=desired_order,
                   ordered=True)
)
name_map = {
    'iMCD-PB_1': 'N_1',
    'iMCD-PB_2': 'N_2',
    'iMCD-PB_8': 'N_3',
    'iMCD-PB_11': 'N_4',
    'iMCD-PB_15': 'N_5',
    'iMCD-PB_16': 'N_6',
    'iMCD-PB_3': 'P_1',
    'iMCD-PB_5': 'P_2',
    'iMCD-PB_7': 'P_3',
    'iMCD-PB_17': 'P_4',
    'iMCD-PB_20': 'P_5',
    'iMCD-PB_21': 'P_6',
    'iMCD-PB_6': 'AID_1',
    'iMCD-PB_18': 'AID_2',
    'iMCD-PB_25': 'AID_3'
}
adata_raw.obs['sample'] = adata_raw.obs['sample'].map(name_map).astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sample', title='sample', show=False)
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_sample.pdf'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_sample.tiff'))

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir, "Main_all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',prefix='Main')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',prefix='Main')

method = "group"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'_color_dict.xls'),sep="\t",index=True,index_label='group')

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
df.to_csv(os.path.join(out_dir,method+'_color_dict.xls'),sep="\t",index=True,index_label='cell_type')

method = "sub_cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
df.to_csv(os.path.join(out_dir,method+'_color_dict.xls'),sep="\t",index=True,index_label='sub_cell_type')

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/15.Mission/fig1_B_markerplot/"
marker_genes=['CD79A','CD79B','MS4A1','ITGA2B','PF4','PPBP','LYZ','VCAN','FCN1','NKG7','GNLY','GZMA','FGFBP2','CD7',
'IL32','JCHAIN','IGKC','MZB1','CD3D','CD3E','CD3G','IL3RA','LILRA4','TCF4']
my_colors = sns.color_palette('Set3', n_colors=10).as_hex()
# 或者自己写列表：['#1f77b4', '#ff7f0e', '#279e68', '#d62728']
adata_raw.uns['cell_type_colors'] = my_colors
sc.pl.stacked_violin(adata_raw, var_names=marker_genes, groupby='cell_type', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, 'Main_stacked_violin.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Main_stacked_violin.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Main_stacked_violin.tiff'),bbox_inches='tight')

marker_genes=['CD79A','ITGA2B','LYZ','NKG7','FGFBP2','JCHAIN','CD3D','IL3RA']
my_colors = sns.color_palette('Set3', n_colors=10).as_hex()
# 或者自己写列表：['#1f77b4', '#ff7f0e', '#279e68', '#d62728']
adata_raw.uns['cell_type_colors'] = my_colors
sc.pl.stacked_violin(adata_raw, var_names=marker_genes, groupby='cell_type', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, 'Main_stacked_violin_top1.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Main_stacked_violin_top1.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Main_stacked_violin_top1.tiff'),bbox_inches='tight')