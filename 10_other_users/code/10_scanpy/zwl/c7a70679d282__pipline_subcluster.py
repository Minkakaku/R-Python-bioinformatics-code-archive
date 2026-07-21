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

# 对Macrophages重新注释细分
adata_sub = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2505152218-01/01.Mission/Macrophages/Macrophages.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2505152218-01/02.Mission/Macrophages/"
# adata_sub = adata_raw[adata_raw.obs.cell_type=="Macrophages", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=30,resolution=0.6,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['Itgam','Adgre1','Cd86','Mrc1','Ly6g','Ly6c1','Cd14','Fcgr4'], cmap='coolwarm',show=False)
plt.savefig(os.path.join(out_dir, 'umap_marker.png'))
plt.savefig(os.path.join(out_dir, 'umap_marker.pdf'))
sc.pl.dotplot(adata, var_names=['Itgam','Adgre1','Cd86','Mrc1'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata, var_names=['Itgam','Adgre1','Cd86','Mrc1'], groupby=['leiden','group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_group-marker.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_group-marker.pdf'), facecolor='white', bbox_inches="tight")

anno_dict={
    '0':'Classical Monocytes',
    '1':'Non-classical Monocytes',
    '2':'Non-classical Monocytes',
    '3':'M1',
    '4':'M2',
    '5':'M2',
    '6':'M2',
    '7':'Classical Monocytes',
    '8':'lowQ',
    '9':'M1'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

adata = adata[adata.obs.cell_type!="lowQ", :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
adata.obs['group'] = adata.obs['group'].cat.reorder_categories(groups, ordered=True)
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'Macrophages.h5ad'))

# 映射所有细分
out_dir="/mnt/NexData_file/home/sh/Projects/SC2505152218-01/02.Mission/mapping_raw/"
adata_raw = sc.read("/mnt/NexData_file/home/sh/Projects/SC2505152218-01/01.Mission/mapping_raw/Main_mapping.h5ad")
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories
adata_raw = adata_raw[adata_raw.obs.sub_cell_type!="Monocyte", :]
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin', show=False)
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))

sc.pl.umap(adata_raw, color=['Gls'], cmap='coolwarm',show=False)
plt.savefig(os.path.join(out_dir, 'umap_Gls.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_Gls.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.umap(adata_raw, color=['Gls'], cmap='coolwarm',show=False,vmin=1)
plt.savefig(os.path.join(out_dir, 'umap_Gls_fiter.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_Gls_fiter.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata_raw, var_names=['Gls'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Gls-cell_type_group.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_Gls-cell_type_group.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata_raw, var_names=['Gls'], groupby='sub_cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Gls-sub_cell_type.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_Gls-sub_cell_type.pdf'), facecolor='white', bbox_inches="tight")