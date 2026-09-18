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

# DTC注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2507163761_fin/result.h5ad")
# 'Distal tubule cells', 'Endothelial cells', 'Epithelial cells','Fibroblasts', 'Intercalated cells', 'Macrophages', 'Mesangial cells','Podocytes', 'Principal cells', 'Proximal tubule cells','loop of Henle'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2507163761-01/02.Mission/DTC/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub = adata_raw[adata_raw.obs.cell_type=="Distal tubule cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=1,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['total_counts','n_genes_by_counts'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_marker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_marker.pdf'),bbox_inches='tight')
marker=["Slc12a3","Slc8a1","Calb1",'Cdh1','Slc16a7','Tmem52b',"Wfdc15b","Cdh16","Gpx6","Clcnka","Kng2","Pvalb","Hoxb5os","Emx1","Scnn1a","Lhx1os","Tfrc","Napsa","Acsm2"]
sc.pl.dotplot(adata, var_names=marker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=marker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_marker.png'))
plt.savefig(os.path.join(out_dir, 'umap_marker.pdf'))

adata_sub = adata[adata.obs['leiden'].isin(['3','7','11','14','15']), :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
    
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.2,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

# sc.tl.dendrogram(adata, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

marker=["Slc12a3","Slc8a1","Calb1",'Cdh1','Slc16a7','Tmem52b',"Wfdc15b","Cdh16","Gpx6","Clcnka","Kng2","Pvalb","Hoxb5os","Emx1","Scnn1a","Lhx1os","Tfrc","Napsa","Acsm2",'Slc13a1']
sc.pl.dotplot(adata, var_names=marker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=marker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_marker.png'))
plt.savefig(os.path.join(out_dir, 'umap_marker.pdf'))

anno_dict={
    '0':'Slc12a3 DTC',
    '1':'Slc13a1 DTC',
    '2':'Slc12a3 DTC',
    '3':'Slc13a1 DTC',
    '4':'Slc12a3 DTC',
    '5':'Slc13a1 DTC'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

group = adata.obs['group'].drop_duplicates().tolist()
groups = ['CMC', 'UDA']
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", legend_loc=None, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'DTC.h5ad'))

# 映射细分，重新聚类
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2507163761_fin/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2507163761-01/02.Mission/mapping_raw/"

adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_sub = adata_raw[~adata_raw.obs.sub_cell_type.isin(['Distal tubule cells']), :]
adata_sub.obs['cell_type'] = adata_sub.obs['cell_type'].replace('Mesangial cells', 'Fibroblasts')
adata_sub.obs['sub_cell_type'] = adata_sub.obs['sub_cell_type'].replace('Mesangial cells', 'Fibroblasts')
desired_order=['Proximal tubule cells', 'Distal tubule cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Intercalated cells', 'Macrophages', 'Podocytes', 'Principal cells', 'loop of Henle']
adata_sub.obs['cell_type'] = (
    pd.Categorical(adata_sub.obs['cell_type'],
                   categories=desired_order,
                   ordered=True)
)

sc.pl.umap(adata_sub, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        sc.pl.umap(adata_sub[adata_sub.obs["group"]==groups[i]], color="cell_type", legend_loc=None, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_sub.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
plot_celltype_proportion(adata_sub, out_dir=out_dir,type='group')
plot_celltype_proportion(adata_sub, out_dir=out_dir,type='sample')
adata_sub.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')

method = "cell_type"
keys = adata_sub.obs[method].cat.categories.to_list()
values = list(adata_sub.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

from matplotlib.colors import LinearSegmentedColormap
top_marker=["Slc12a3","Gramd1b","Pecam1","Krt19","Lrp1","Atp6v1g3","Ctss","Akap12","Tns3","Ptpro","Phactr1","Slc17a3","Slc12a1"]
cmap = LinearSegmentedColormap.from_list('cut_gray_red', ['#cccccc', '#cccccc', '#ff0000'], N=256)
for i in top_marker:
        sc.pl.umap(adata_sub, color=i, title=i, show=False, cmap=cmap, vmin=None)
        plt.savefig(os.path.join(out_dir, f'umap/{i}_umap.png'),bbox_inches='tight')
        plt.savefig(os.path.join(out_dir, f'umap/{i}_umap.pdf'),bbox_inches='tight')

adata_sub.uns['log1p'] = {'base': None}
sc.tl.rank_genes_groups(adata_sub,groupby='cell_type', method='wilcoxon',pts=True)
find_group_marker(adata=adata_sub,species='mouse',output=out_dir)

sc.pl.rank_genes_groups_dotplot(
    adata_sub,
    groupby='cell_type',
    key='rank_genes_groups',
    n_genes=5,
    show=False)
plt.savefig(os.path.join(out_dir, 'top5_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top5_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

select_marker=["Slc17a3","Cubn","Cda","Slc12a3","Trpm6","Abca13","Pecam1","Eng","Kdr","Krt19","Krt18","Krt7","Lrp1","Pdgfrb","Col1a2","Atp6v1g3","Atp6v0d2","Clcnkb","Ctss","Ms4a6c","Adgre1","Ptpro","Podxl","Nphs2","Phactr1","Scnn1b","Tmem45b","Slc12a1",'Fxyd2','Kcnj1']
sc.pl.dotplot(adata_sub, var_names=select_marker, groupby='cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type_scale.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata_sub, var_names=select_marker, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker_cell_type.pdf'),bbox_inches='tight')
