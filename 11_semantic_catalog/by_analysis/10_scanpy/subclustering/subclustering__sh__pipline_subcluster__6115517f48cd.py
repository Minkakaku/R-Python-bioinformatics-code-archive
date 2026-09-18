import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
import pandas as pd
from harmonypy import run_harmony
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

# 对Astrocyte注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20250765/result.h5ad")
# 'Astrocyte', 'Choroid plexus cells', 'Endothelial cells','Ependymal cells', 'Microglia', 'Oligodendrocytes', 'Smooth muscle cells','T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250765/01.Mission/Astrocyte/"
# adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['B cells', 'Plasma cells'])]
adata_sub = adata_raw[adata_raw.obs.cell_type=="Astrocyte", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=30,resolution=0.1,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['Timp1','Gap43','Hspb1','Lcn2','Gfap','Ifitm3','H2-K1','H2-T23','H2-D1','Aqp4','Gfap','Vim','Id1','Id3','Id4'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Astrocyte_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Astrocyte_diffgenes.pdf'))

anno_dict={
    '0':'Astrocyte C1',
    '1':'Astrocyte C2',
    '2':'Astrocyte C3',
    '3':'Astrocyte C4',
		'4':'Astrocyte C5',
    '5':'Astrocyte C6'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

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
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
adata.write(os.path.join(out_dir, 'Astrocyte.h5ad'))



# 对Microglia注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20250765/result.h5ad")
# 'Astrocyte', 'Choroid plexus cells', 'Endothelial cells','Ependymal cells', 'Microglia', 'Oligodendrocytes', 'Smooth muscle cells','T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250765/01.Mission/Microglia/"
# adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['B cells', 'Plasma cells'])]
adata_sub = adata_raw[adata_raw.obs.cell_type=="Microglia", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=30,resolution=0.5,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['Itgax','Spp1'], cmap='coolwarm')
sc.pl.umap(adata, color=['Acr','Itgam','Trem2','Cx3cr1','Aif1','P2ry12'], cmap='coolwarm')
sc.pl.umap(adata, color=['Lck','Cx3cr1','P2ry12','Ptprc'], cmap='coolwarm')
sc.pl.umap(adata, color=['Cd163','Apoe','Lpl','Ifi27','Ifitm3','Cx3cr1','Selplg'], cmap='coolwarm')
sc.pl.umap(adata, color=['Cd163','Ifi27','Ifitm3'], cmap='coolwarm')
sc.pl.umap(adata, color=['Itgax','Spp1','Itgam','Trem2','Cx3cr1','Aif1','P2ry12','Lck','P2ry12','Ptprc','Cd163','Apoe','Ifitm3','Selplg'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Microglia_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Microglia_diffgenes.pdf'))

anno_dict={
    '0':'repair-associated microglia',
    '1':'repair-associated microglia',
    '2':'repair-associated microglia',
    '3':'repair-associated microglia',
    '4':'CD11c Microglia',
    '5':'glioma-associated microglia',
    '6':'CD11c Microglia',
    '7':'repair-associated microglia',
    '8':'glioma-associated microglia',
    '9':'glioma-associated microglia'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

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
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
adata.write(os.path.join(out_dir, 'Microglia.h5ad'))


# 对Microglia注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20250765/result.h5ad")
# 'Astrocyte', 'Choroid plexus cells', 'Endothelial cells','Ependymal cells', 'Microglia', 'Oligodendrocytes', 'Smooth muscle cells','T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250765/01.Mission/T cells/"
# adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['B cells', 'Plasma cells'])]
adata_sub = adata_raw[adata_raw.obs.cell_type=="T cells", :]
sc.tl.leiden(adata_sub, resolution=0.1)
# adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.1,n_pcs = 15,output=out_dir)
# adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='leiden', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'sub_cluster_umap.png'))
    plt.savefig(os.path.join(out_dir, 'sub_cluster_umap.pdf'))

sc.pl.umap(adata_sub, color=['Cd4','Cd3d','Cd3e','Cd8a','Ncam1','Ncr1','Nkg7','Klrd1','Gzma'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Tcells_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Tcells_diffgenes.pdf'))

# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20250765/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250765/01.Mission/mapping_raw/"
adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250765/01.Mission/Astrocyte/Astrocyte.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250765/01.Mission/Microglia/Microglia.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.write(os.path.join(out_dir, 'mapping_raw.h5ad'))
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group')