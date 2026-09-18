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

# Enteroendocrine cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509235064-01/02.Mission/mapping_raw/Main_mapping.h5ad")
# 'B cells', 'Endothelial cells', 'Enteroendocrine cells', 'Epithelial cells', 'Fibroblasts', 'Macrophages', 'Neutrophils', 'Pericyte', 'Plasma cells', 'SMC', 'Schwann cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509235064-01/03.Mission/Enteroendocrine/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

ECmarker=['Gcg','Pyy','Pcsk1','Pcsk2','Cck','Sct','Sst','Tph1','Slc18a1','Fev','Nts','Hdc']
sc.pl.dotplot(adata_raw, var_names=ECmarker, groupby='cluster', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_EC.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_EC.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata_raw, var_names=ECmarker, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_EC_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_EC_cell_type.pdf'),bbox_inches='tight')
sc.pl.umap(adata_raw, color=ECmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_EC_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_EC_diffgenes.pdf'))


adata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509235064-01/01.Mission/Enteroendocrine/Enteroendocrine.h5ad")
sc.pl.dotplot(adata, var_names=ECmarker, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_ECsub_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_ECsub_cell_type.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=ECmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_ECsub_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_ECsub_diffgenes.pdf'))

adata_sub=adata.copy()
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=1.2,n_pcs = 20,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()
sc.pl.dotplot(adata, var_names=ECmarker, groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_ECsub.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_ECsub.pdf'),bbox_inches='tight')

sc.tl.dendrogram(adata,groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'EC cells',
    '1':'Kit Enteroendocrine',
    '2':'Enterochromaffin-like cells',
    '3':'Kit Enteroendocrine',
    '4':'I cells',
    '5':'L cells',
    '6':'Kit Enteroendocrine',
    '7':'L cells',
    '8':'L cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['CON','DSS','RI']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':7, 'h':5})
adata.write(os.path.join(out_dir, 'Enteroendocrine.h5ad'))

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)


#### 基因表达展示
geneset_plot=['Tph1','Ido1','Ahr','Htr4','Slc6a4','Tdo2','Kat2a','Kmo','Haao','Qprt']
pathways = {'geneset_plot': geneset_plot}

base_dir = "/mnt/NexData_file/home/sh/Projects/SC2509235064-01/03.Mission/"
for name, genes in pathways.items():
    out_dir = os.path.join(base_dir, name)
    os.makedirs(out_dir, exist_ok=True)
    sc.pl.dotplot(adata, var_names=genes, groupby='cell_type', dendrogram=False, show=False)
    plt.savefig(os.path.join(out_dir, 'genelist_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'genelist_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
    sc.pl.umap(adata, color=genes, cmap='coolwarm', ncols=3)
    plt.savefig(os.path.join(out_dir, 'genelist_umap.png'))
    plt.savefig(os.path.join(out_dir, 'genelist_umap.pdf'))
    for g in genes:
        sc.pl.dotplot(adata, var_names=g, groupby=['cell_type','group'], dendrogram=False, show=False)
        plt.savefig(os.path.join(out_dir, f'{g}_dotplot_group_celltype.png'), facecolor='white', bbox_inches='tight', dpi=300)
        plt.savefig(os.path.join(out_dir, f'{g}_dotplot_group_celltype.pdf'), facecolor='white', bbox_inches='tight', dpi=300)
        plt.close()


# 映射所有细分 下面都没做
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2509235064/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509235064-01/01.Mission/mapping_raw/"

adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509235064-01/01.Mission/Enteroendocrine/Enteroendocrine.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')