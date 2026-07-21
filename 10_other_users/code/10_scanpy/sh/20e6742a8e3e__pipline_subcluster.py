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

def plot_celltype_proportion(adata, out_dir,type='sample', fig_size=None):
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
    file_name=f"celltype_proportion_{type}"
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)

# 调整大类注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2507163761_fin/result.h5ad")
# 'Distal tubule cells', 'Endothelial cells', 'Epithelial cells','Fibroblasts', 'Intercalated cells', 'Macrophages', 'Mesangial cells','Podocytes', 'Principal cells', 'Proximal tubule cells','loop of Henle'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2507163761-01/01.Mission/Main/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cluster', title='Cluster', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Cluster.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Cluster.pdf'))


marker=["Slc5a2","Slc5a12","Slc22a6","Slco1a6","Slc17a3","Eci3","Slc22a7","Slc22a13","Slc13a1","Slc7a13","Cubn","Lrp2","Aldob","Cda","Kcnj15","Slc6a18","Them7","Cyp2j5", ####Proximal tubule cells
                "Slc12a3","Slc8a1","Calb1",'Cdh1','Slc16a7','Tmem52b',"Wfdc15b","Cdh16","Gpx6","Clcnka","Kng2","Pvalb","Hoxb5os","Emx1","Scnn1a","Lhx1os","Tfrc","Napsa","Acsm2",###Distal tubule cells
                "Akap12","Ccn2","Myl9","Mylk","Fhl2","Tns3","Pdgfb",###Mesangial cells
                "Slc12a1","Umod","Cldn16",'Kcnj1','Fxyd2',"Egf",##loop of Henle
                "Pecam1","Eng","Kdr","Cd34","Ackr1","Vwf","Cldn5",####Endothelial cells
                "Dcn","Col1a1","Col4a1","Col4a2","Col1a2",####Fibroblasts
                "Aqp2","Aqp3",######Principal cells
                "Ctss","Ms4a6c","Adgre1","Plac8","Cd68","C1qa","C1qb","Mrc1",#Macrophages
                "Avpr1a","Atp6v1g3","Atp6v0d2",#Intercalated cells
                "Nphs2","Nphs1","Podxl","Ptpro"]###Podocytes
sc.pl.dotplot(adata_raw, var_names=marker, groupby='cluster', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata_raw, var_names=marker, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type.pdf'),bbox_inches='tight')
sc.pl.umap(adata_raw, color=marker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_marker.png'))
plt.savefig(os.path.join(out_dir, 'umap_marker.pdf'))

adata_raw.uns['log1p'] = {'base': None}
sc.tl.rank_genes_groups(adata_raw,groupby='cell_type', method='wilcoxon',pts=True)
find_group_marker(adata=adata_raw,species='mouse',output=out_dir)

sc.pl.rank_genes_groups_dotplot(
    adata_raw,
    groupby='cell_type',
    key='rank_genes_groups',
    n_genes=3,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

select_marker=["Slc12a3","Slc8a1","Calb1","Pecam1","Eng","Kdr","Krt19","Krt18","Krt7","Lrp1","Pdgfrb","Col1a2","Atp6v1g3","Atp6v0d2","Clcnkb","Ctss","Ms4a6c","Adgre1","Akap12","Tns3","Pdgfb","Ptpro","Podxl","Nphs2","Phactr1","Scnn1b","Tmem45b","Slc17a3","Cubn","Cda","Slc12a1",'Fxyd2','Kcnj1']
sc.pl.dotplot(adata_raw, var_names=select_marker, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata_raw, var_names=select_marker, groupby='cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type_scale.pdf'),bbox_inches='tight')

groups = ['CMC', 'UDA']
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=None, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', legend_loc=None, title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_nolegend.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_nolegend.pdf'))

top_marker=["Slc12a3","Pecam1","Krt19","Lrp1","Atp6v1g3","Ctss","Akap12","Tns3","Ptpro","Phactr1","Slc17a3","Slc12a1"]
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('cut_gray_red', ['#cccccc', '#cccccc', '#ff0000'], N=256)
for i in top_marker:
        sc.pl.umap(adata_raw, color=i, title=i, show=False, cmap=cmap, vmin=None)
        plt.savefig(os.path.join(out_dir, f'umap/{i}_umap.png'),bbox_inches='tight')
        plt.savefig(os.path.join(out_dir, f'umap/{i}_umap.pdf'),bbox_inches='tight')

top_marker=["Pecam1","Krt19","Lrp1","Ctss"]
cmap = LinearSegmentedColormap.from_list('gray_red', ['#cccccc', '#ff0000'], N=256)
for i in top_marker:
        sc.pl.umap(adata_raw, color=i, title=i, show=False, cmap=cmap, vmin=None)
        plt.savefig(os.path.join(out_dir, f'umap/{i}_umap.png'),bbox_inches='tight')
        plt.savefig(os.path.join(out_dir, f'umap/{i}_umap.pdf'),bbox_inches='tight')

adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
method = "group"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')