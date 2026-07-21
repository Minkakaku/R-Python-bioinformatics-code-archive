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

# 对Macrophages注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20253251/02.Mission/mapping_raw/Main_mapping.h5ad")
# 'Proximal tubule cells', 'Macrophages', 'Distal tubule cells','Proliferation', 'Endothelial cells', 'T cells', 'Neutrophils','Dendritic cells', 'Fibroblasts', 'Erythrocyte', 'Principal cells','B cells', 'Smooth muscle cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/Macrophages/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Macrophages", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.8,n_pcs = 10,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['Lyz2','Vcan','Ly6c2','Cd14','Fcgr3','Apoe','C1qa','Cd86','Cd274','Mrc1','Cd163','Clec9a','Cd209a','Siglech','Tcf4'], cmap='coolwarm',show=False)
plt.savefig(os.path.join(out_dir, 'umap_gene.png'))
plt.savefig(os.path.join(out_dir, 'umap_gene.pdf'))
sc.pl.dotplot(adata, var_names=['Lyz2','Vcan','Ly6c2','Cd14','Fcgr3','Apoe','C1qa','Cd86','Cd274','Mrc1','Cd163','Clec9a','Cd209a','Siglech','Tcf4'], groupby='leiden', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_gene.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_gene.pdf'), facecolor='white', bbox_inches="tight")

anno_dict={
    '0':'Jun Macrophage',
    '1':'Gpx3 Macrophage',
    '2':'Gda Macrophage',
    '3':'Jun Macrophage',
    '4':'Inpp5d Macrophage',
    '5':'Aoah Macrophage',
    '6':'Gpx3 Macrophage',
    '7':'Inpp5d Macrophage',
    '8':'Gda Macrophage',
    '9':'Gpx3 Macrophage',
    '10':'Gpx3 Macrophage'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

# sc.tl.rank_genes_groups(adata,groupby='cell_type', method='wilcoxon', pts=True)
# sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False,show=False)
# save_fig(plt, out_dir, file_name='rank_genes_groups_tmp', fig_size={'width':10, 'height':12})

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['H0', 'H2', 'H8', 'H24', 'H72']
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

sc.tl.rank_genes_groups(adata,groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)


# 对T cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20253251/02.Mission/mapping_raw/Main_mapping.h5ad")
# 'Proximal tubule cells', 'Macrophages', 'Distal tubule cells','Proliferation', 'Endothelial cells', 'T cells', 'Neutrophils','Dendritic cells', 'Fibroblasts', 'Erythrocyte', 'Principal cells','B cells', 'Smooth muscle cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/Tcells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="T cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

# 分出来的是Basophil
adata_sub = adata[~adata.obs.leiden.isin(["14"]), :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.8,n_pcs = 10,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

# 分出来是双胞
adata_sub = adata[~adata.obs.leiden.isin(["9","5","13"]), :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

sc.tl.dendrogram(adata,groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

marker=["Cd3d", "Cd4", "Cd8a","Cd8b1",'Cxcr5','Bcl6','Batf','Foxp3',"Sell", "Ccr7", "Tcf7", "Lef1",'Cd69','Rgs1',"Pdcd1", "Tox", "Havcr2", "Lag3", "Prf1", "Gzmb", "Ifng","Gzmk",'Gzma',"Ncam1",'Mki67','Top2a','Klrd1','Xcl1']
marker=['Cd79a','Ms4a1','Jchain','Ly6c2','C1qa','Lyz2','Vcan','Apoe']
sc.pl.dotplot(adata, var_names=marker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=marker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_marker.png'))
plt.savefig(os.path.join(out_dir, 'umap_marker.pdf'))
sc.pl.umap(adata, color=['total_counts','n_genes_by_counts'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_counts.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_counts.pdf'),bbox_inches='tight')

anno_dict={
    '0':'NK',
    '1':'Treg',
    '2':'Cytotoxic T',
    '3':'Th1',
    '4':'Th1',
    '5':'Tfh',
    '6':'Activated T',
    '7':'NK',
    '8':'Th2',
    '9':'Naive T',
    '10':'Proliferative T',
    '11':'lowQ',
    '12':'NK'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata = adata[~adata.obs.cell_type.isin(["lowQ"]), :]

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['H0', 'H2', 'H8', 'H24', 'H72']
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
adata.write(os.path.join(out_dir, 'Tcells.h5ad'))

sc.tl.rank_genes_groups(adata,groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

# 对Neutrophils注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20253251/02.Mission/mapping_raw/Main_mapping.h5ad")
# 'Proximal tubule cells', 'Macrophages', 'Distal tubule cells','Proliferation', 'Endothelial cells', 'T cells', 'Neutrophils','Dendritic cells', 'Fibroblasts', 'Erythrocyte', 'Principal cells','B cells', 'Smooth muscle cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/Neutrophils/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Neutrophils", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'Pglyrp1 Neutrophils',
    '1':'Il1rn Neutrophils',
    '2':'Cmss1 Neutrophils',
    '3':'Il1rn Neutrophils',
    '4':'Il1rn Neutrophils',
    '5':'Pglyrp1 Neutrophils',
    '6':'Cmss1 Neutrophils',
    '7':'Pglyrp1 Neutrophils',
    '8':'Cmss1 Neutrophils'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.tl.rank_genes_groups(adata,groupby='cell_type', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False,show=False)
save_fig(plt, out_dir, file_name='rank_genes_groups_tmp', fig_size={'width':10, 'height':12})

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['H0', 'H2', 'H8', 'H24', 'H72']
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
adata.write(os.path.join(out_dir, 'Neutrophils.h5ad'))

sc.tl.rank_genes_groups(adata,groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20253251/02.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/mapping_raw/"

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/Macrophages/Macrophages.h5ad")
bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/Neutrophils/Neutrophils.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs.sub_cell_type.isin(["T cells"]), :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin', show=False)
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))