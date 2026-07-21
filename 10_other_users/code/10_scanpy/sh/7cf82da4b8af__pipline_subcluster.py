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

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/proportion/"

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/05.Mission/mapping_raw/Main_mapping.h5ad")
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_celltype.png'))
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Negative', 'N')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Positive', 'P')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Rheumatic Diseases', 'AID')
# 'B cells', 'Erythroid cells', 'Mast cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
adata_raw.obs.to_csv(os.path.join(out_dir, "Main_all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',prefix='Main')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',prefix='Main')
sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups_dotplot( adata_raw, groupby='cell_type', key='rank_genes_groups', n_genes=3, show=False)
plt.savefig(os.path.join(out_dir, 'Main_top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'Main_top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))

# "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"
adata_raw.uns['group_colors']=["#F8766D", "#00BF7D", "#00B0F6"]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='group', title='group', show=False)
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_group.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sample', title='sample', show=False)
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_sample.pdf'))
marker_genes=['CD79A','HBB','ITGA2B','LYZ','GNLY','NCAM1','CSF3R','JCHAIN','CD3D','IL3RA']
my_colors = sns.color_palette('Set3', n_colors=10).as_hex()
# 或者自己写列表：['#1f77b4', '#ff7f0e', '#279e68', '#d62728']
adata_raw.uns['cell_type_colors'] = my_colors
sc.pl.stacked_violin(adata_raw, var_names=marker_genes, groupby='cell_type', show=False, swap_axes=True, dendrogram=False)
sc.pl.stacked_violin(
    adata_raw,
    var_names=marker_genes,
    groupby='cell_type',
    show=False,
    swap_axes=True,
    dendrogram=False
)
plt.savefig(os.path.join(out_dir, 'Main_stacked_violin.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'Main_stacked_violin.png'),bbox_inches='tight')
for i in range(len(marker_genes)):
    i=marker_genes[i]
    with rc_context({'figure.figsize': (12, 2), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
        sc.pl.violin(adata_raw,keys=i,groupby='cell_type',jitter=False,ylabel=i)
        plt.savefig(os.path.join(out_dir, f'Main_{i}_violin.pdf'),bbox_inches='tight')
        plt.savefig(os.path.join(out_dir, f'Main_{i}_violin.png'),bbox_inches='tight')

from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('gray_red', ['#D3D3D3', '#FF0000'])
for i in range(len(marker_genes)):
    i=marker_genes[i]
    with rc_context({'figure.figsize': (7, 7), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
        sc.pl.umap(adata_raw, color=i, title=i, show=False,cmap=cmap,size=5)
        plt.savefig(os.path.join(out_dir, f'Main_{i}_umap.pdf'),bbox_inches='tight')
        plt.savefig(os.path.join(out_dir, f'Main_{i}_umap.png'),bbox_inches='tight')

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,"Main_"+method+'_color_dict.xls'),sep="\t",index=True,index_label='cell_type')

method = "group"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'_color_dict.xls'),sep="\t",index=True,index_label='group')



adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/05.Mission/Myeloid cells/Myeloid.h5ad")
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'Myeloid_dimplot_umap_celltype.png'))
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Negative', 'N')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Positive', 'P')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Rheumatic Diseases', 'AID')
# 'B cells', 'Erythroid cells', 'Mast cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
adata_raw.obs.to_csv(os.path.join(out_dir, "Myeloid_all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',prefix='Myeloid')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',prefix='Myeloid')
sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups_dotplot( adata_raw, groupby='cell_type', key='rank_genes_groups', n_genes=3, show=False)
plt.savefig(os.path.join(out_dir, 'Myeloid_top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'Myeloid_top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()
adata_raw.write(os.path.join(out_dir, 'Myeloid_mapping.h5ad'))

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,"Myeloid_"+method+'_color_dict.xls'),sep="\t",index=True,index_label='cell_type')


adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/02.Mission/Plasma cells/Plasma.h5ad")
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'Plasma_dimplot_umap_celltype.png'))
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Negative', 'N')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Positive', 'P')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Rheumatic Diseases', 'AID')
# 'B cells', 'Erythroid cells', 'Mast cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
adata_raw.obs.to_csv(os.path.join(out_dir, "Plasma_all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',prefix='Plasma')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',prefix='Plasma')
sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups_dotplot( adata_raw, groupby='cell_type', key='rank_genes_groups', n_genes=3, show=False)
plt.savefig(os.path.join(out_dir, 'Plasma_top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'Plasma_top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()
adata_raw.write(os.path.join(out_dir, 'Plasma_mapping.h5ad'))

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,"Plasma_"+method+'_color_dict.xls'),sep="\t",index=True,index_label='cell_type')


adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/02.Mission/B cells/Bcells.h5ad")
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'Bcells_dimplot_umap_celltype.png'))
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Negative', 'N')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Positive', 'P')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('Rheumatic Diseases', 'AID')
# 'B cells', 'Erythroid cells', 'Mast cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
adata_raw.obs.to_csv(os.path.join(out_dir, "Bcells_all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',prefix='Bcells')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',prefix='Bcells')
sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups_dotplot( adata_raw, groupby='cell_type', key='rank_genes_groups', n_genes=3, show=False)
plt.savefig(os.path.join(out_dir, 'Bcells_top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'Bcells_top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()
adata_raw.write(os.path.join(out_dir, 'Bcells_mapping.h5ad'))

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,"Bcells_"+method+'_color_dict.xls'),sep="\t",index=True,index_label='cell_type')