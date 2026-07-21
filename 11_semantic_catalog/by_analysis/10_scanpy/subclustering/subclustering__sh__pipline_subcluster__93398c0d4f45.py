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

# 对Principal cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/wzb/backup/scRNA/SC2509094836_2/2026.03.19/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'B cells', 'Basal cells', 'Distal tubule cells', 'Ductal cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Myeloid cells', 'NK cells', 'Neutrophils', 'Principal cells', 'Proximal tubule cells', 'Smooth muscle cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094836-05/03.Mission/Principal/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["Principal cells"]), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

cell_markers = {
    "Cortical Collecting Duct Principal Cell":["SCNN1G", "SCNN1B", "FXYD4", "SOX5", "PDE10A", "SLC25A29", "ST6GAL1", "PAPPA"],
    "Outer Medullary Collecting Duct Principal Cell":["SCNN1G", "SCNN1B", "FXYD4", "SOX5", "SYK", "FAM81A", "PROM1", "KCNK13"],
    "Inner Medullary Collecting Duct Cell":["FXYD4", "SOX5", "PHACTR1", "PCDH7", "SLC14A2", "HS3ST5"],
    "Papillary Epithelial Cells":["TACSTD2", "TP63", "GPX2", "FXYD3", "KRT5"]
}
markers=["SCNN1G", "SCNN1B", "FXYD4", "SOX5", "PDE10A", "SLC25A29", "ST6GAL1", "PAPPA","SYK", "FAM81A", "PROM1", "KCNK13","PHACTR1", "PCDH7", "SLC14A2", "HS3ST5","TACSTD2", "TP63", "GPX2", "FXYD3", "KRT5"]
sc.pl.umap(adata, color=markers, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

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
    '0':'IMCD',
    '1':'IMCD',
    '2':'OMCD-PC',
    '3':'IC',
    '4':'CCD-PC',
    '5':'OMCD-PC'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'PC_IC.h5ad'))

# sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
# save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
# sc.pl.umap(adata_raw, color='sample', title='sample', show=False)
# save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

# sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


# 对Epithelial cells,Basal cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/wzb/backup/scRNA/SC2509094836_2/2026.03.19/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'B cells', 'Basal cells', 'Distal tubule cells', 'Ductal cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Myeloid cells', 'NK cells', 'Neutrophils', 'Principal cells', 'Proximal tubule cells', 'Smooth muscle cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094836-05/03.Mission/Epithelial_Basal/"
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['Epithelial cells','Basal cells']), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.2,n_pcs = 10,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

cell_markers = {
    "Cortical Collecting Duct Principal Cell":["SCNN1G", "SCNN1B", "FXYD4", "SOX5", "PDE10A", "SLC25A29", "ST6GAL1", "PAPPA"],
    "Outer Medullary Collecting Duct Principal Cell":["SCNN1G", "SCNN1B", "FXYD4", "SOX5", "SYK", "FAM81A", "PROM1", "KCNK13"],
    "Inner Medullary Collecting Duct Cell":["FXYD4", "SOX5", "PHACTR1", "PCDH7", "SLC14A2", "HS3ST5"],
    "Papillary Epithelial Cells":["TACSTD2", "TP63", "GPX2", "FXYD3", "KRT5"]
}
markers=["SCNN1G", "SCNN1B", "FXYD4", "SOX5", "PDE10A", "SLC25A29", "ST6GAL1", "PAPPA","SYK", "FAM81A", "PROM1", "KCNK13","PHACTR1", "PCDH7", "SLC14A2", "HS3ST5","TACSTD2", "TP63", "GPX2", "FXYD3", "KRT5"]
sc.pl.umap(adata, color=markers, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

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
    '0':'PapE',
    '1':'OMCD-PC',
    '2':'OMCD-PC',
    '3':'OMCD-PC',
    '4':'PapE',
    '5':'IMCD'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Epithelial_Basal.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)


# 对T cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/wzb/backup/scRNA/SC2509094836_2/2026.03.19/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'B cells', 'Basal cells', 'Distal tubule cells', 'Ductal cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Myeloid cells', 'NK cells', 'Neutrophils', 'Principal cells', 'Proximal tubule cells', 'Smooth muscle cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094836-05/03.Mission/Tcells/"
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['T cells']), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
adata_raw.obs['cell_type'].cat.categories
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub = adata[~adata.obs['cell_type'].isin(['6']), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

markers=['CD4','CD8A','CD8B','FOXP3','GZMA','NKG7','CD3D','CD3E','NCAM1','SELL','LEF1','CD69','RGS1','IL17A','TOX','PDCD1','TRDC',"TRGC1", "TRGC2",'MKI67']
sc.pl.umap(adata, color=markers, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=markers, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

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
    '0':'CD4T naive',
    '1':'CD4T memory',
    '2':'CD8T cytotoxic',
    '3':'CD8T cytotoxic',
    '4':'Treg',
    '5':'CD8T exhaust',
    '6':'CD4T naive'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Tcells.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/wzb/backup/scRNA/SC2509094836_2/2026.03.19/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094836-05/03.Mission/mapping_raw/"
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509094836-05/03.Mission/Epithelial_Basal/Epithelial_Basal.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509094836-05/03.Mission/Principal/PC_IC.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['T cells']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')