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

# 对T cells注释细分
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250829/00.data/recluster/anno_result.h5ad")
# adata_raw = adata_raw.raw.to_adata()
# 'B cells', 'Mast cells', 'Monocytes', 'NK cells', 'Neutrophils', 'Plasmacytoid dendritic cells', 'Platelets', 'T cells'
out_dir="/data/NAS02/sh/Projects/FWSC20250829/01.Mission/T cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="T cells", :]
sc.pp.log1p(adata_sub)
adata_sub.obs['platform'] = adata_sub.obs['platform'].astype(str).astype(object).astype("category")
del adata_sub.obs['RNA_snn_res.0.3']
del adata_sub.obs['seurat_clusters']
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=30,resolution=0.8,n_pcs = 20,batch_effect = 3,output=out_dir)
adata = adata_recluster.sub_cluster()

# 运行 Harmony
# adata_harmony = run_harmony(adata, adata.obs, ['platform'])
# sc.pp.neighbors(adata_harmony)
# sc.tl.umap(adata_harmony)
# sc.pl.umap(adata, color='leiden')
# plt.savefig(os.path.join(out_dir, 'harmony_umap_celltype.png'))
# plt.savefig(os.path.join(out_dir, 'harmony_umap_celltype.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['CD4','CD8A','TRGC1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_T_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_T_diffgenes.pdf'))

anno_dict={
    '0':'CD4 Tcells',
    '1':'CD8 Tcells',
    '2':'CD4 Tcells',
    '3':'CD8 Tcells',
		'4':'CD4 Tcells',
    '5':'CD4 Tcells',
    '6':'CD8 Tcells',
    '7':'Gamma-delta Tcells',
		'8':'CD8 Tcells',
    '9':'CD4 Tcells',
    '10':'CD4 Tcells',
		'11':'CD4 Tcells'
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

samples = adata.obs['sample'].drop_duplicates().tolist()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["sample"]==samples[i]], color="cell_type", legend_loc=legend, title=samples[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'T cells.h5ad'))

# 对CD4 Tcells注释细分
out_dir="/data/NAS02/sh/Projects/FWSC20250829/01.Mission/T cells/CD4 Tcells/"
adata_sub = adata[adata.obs.cell_type=="CD4 Tcells", :]
sc.pl.umap(adata_sub, color='leiden', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

groups = adata_sub.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_sub[adata_sub.obs["group"]==groups[i]], color="leiden", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.umap(adata_sub, color=['CCR4', 'FOXP3', 'IL2RA'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'CD4 Treg_umap.png'))
plt.savefig(os.path.join(out_dir, 'CD4 Treg_umap.pdf'))
sc.pl.stacked_violin(adata, var_names=['CCR4', 'FOXP3', 'IL2RA'], groupby='leiden', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, f'CD4 Treg_violin.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, f'CD4 Treg_violin.png'),bbox_inches='tight')
sc.pl.umap(adata_sub, color=['CCR7','SELL','CD5'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'CD4 naive Tcells_umap.png'))
plt.savefig(os.path.join(out_dir, 'CD4 naive Tcells_umap.pdf'))
sc.pl.stacked_violin(adata_sub, var_names=['CCR7','SELL','CD5'], groupby='leiden', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, f'CD4 naive Tcells_violin.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, f'CD4 naive Tcells_violin.png'),bbox_inches='tight')
sc.pl.umap(adata_sub, color=['IL15','MBD2','PRDM1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'CD4 Tem_umap.png'))
plt.savefig(os.path.join(out_dir, 'CD4 Tem_umap.pdf'))
sc.pl.stacked_violin(adata_sub, var_names=['IL15','MBD2','PRDM1'], groupby='leiden', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, f'CD4 Tem_violin.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, f'CD4 Tem_violin.png'),bbox_inches='tight')
anno_dict={
    '0':'CD4 Tem',
    '2':'CD4 naive Tcells',
		'4':'CD4 naive Tcells',
    '5':'CD4 Temra',
    '9':'CD4 Treg',
    '10':'CD4 naive Tcells',
		'11':'CD4 naive Tcells'
}
adata_sub.obs['cell_type'] = [anno_dict[clust] for clust in adata_sub.obs['leiden']]
adata_sub.obs['cell_type'] = adata_sub.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))


# 对CD8 Tcells注释细分
out_dir="/data/NAS02/sh/Projects/FWSC20250829/01.Mission/T cells/CD8 Tcells/"
adata_sub = adata[adata.obs.cell_type=="CD8 Tcells", :]
sc.pl.umap(adata_sub, color='leiden', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

groups = adata_sub.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_sub[adata_sub.obs["group"]==groups[i]], color="leiden", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.umap(adata_sub, color=['CCR7', 'SELL', 'IL7R','CD27','PRF1','GZMA','CCL5','GPR183','S1PR1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'CD8 TCM_umap.png'))
plt.savefig(os.path.join(out_dir, 'CD8 TCM_umap.pdf'))
sc.pl.stacked_violin(adata, var_names=['CCR7', 'SELL', 'IL7R','CD27','PRF1','GZMA','CCL5','GPR183','S1PR1'], groupby='leiden', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, f'CD8 TCM_violin.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, f'CD8 TCM_violin.png'),bbox_inches='tight')
sc.pl.umap(adata_sub, color=['CCR7','SELL','CD5'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'CD8 naive Tcells_umap.png'))
plt.savefig(os.path.join(out_dir, 'CD8 naive Tcells_umap.pdf'))
sc.pl.stacked_violin(adata_sub, var_names=['CCR7','SELL','CD5'], groupby='leiden', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, f'CD8 naive Tcells_violin.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, f'CD8 naive Tcells_violin.png'),bbox_inches='tight')
sc.pl.umap(adata_sub, color=['GZMK','CXCR4','CXCR3'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'CD8 Tem_umap.png'))
plt.savefig(os.path.join(out_dir, 'CD8 Tem_umap.pdf'))
sc.pl.stacked_violin(adata_sub, var_names=['GZMK','CXCR4','CXCR3'], groupby='leiden', show=False, swap_axes=True, dendrogram=False)
plt.savefig(os.path.join(out_dir, f'CD8 Tem_violin.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, f'CD8 Tem_violin.png'),bbox_inches='tight')

anno_dict={
    '1':'CD8 Cytotoxic Tcells',
    '3':'CD8 naive Tcells',
		'6':'CD8 Tcm',
    '8':'MAIT'
}
adata_sub.obs['cell_type'] = [anno_dict[clust] for clust in adata_sub.obs['leiden']]
adata_sub.obs['cell_type'] = adata_sub.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# 映射CD4/CD8 Tcells注释细分
out_dir="/data/NAS02/sh/Projects/FWSC20250829/01.Mission/T cells/"
anno_dict={
    '0':'CD4 Tem',
    '1':'CD8 Cytotoxic Tcells',
    '2':'CD4 naive Tcells',
    '3':'CD8 naive Tcells',
		'4':'CD4 naive Tcells',
    '5':'CD4 Cytotoxic Tcells',
    '6':'CD8 Tcm',
    '7':'Gamma-delta Tcells',
    '8':'MAIT',
    '9':'CD4 Treg',
    '10':'CD4 naive Tcells',
		'11':'CD4 naive Tcells'
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

samples = adata.obs['sample'].drop_duplicates().tolist()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["sample"]==samples[i]], color="cell_type", legend_loc=legend, title=samples[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'Tcells_anno.h5ad'))

sc.pl.umap(adata_sub, color=['TNFRSF8'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'CD8_umap.png'))



# 对NK cells注释细分
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250829/00.data/recluster/anno_result.h5ad")
# adata_raw = adata_raw.raw.to_adata()
# 'B cells', 'Mast cells', 'Monocytes', 'NK cells', 'Neutrophils', 'Plasmacytoid dendritic cells', 'Platelets', 'T cells'
out_dir="/data/NAS02/sh/Projects/FWSC20250829/01.Mission/NK cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="NK cells", :]
sc.pl.umap(adata_sub, color=['NCAM1','FCGR3A'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_NK_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_NK_diffgenes.pdf'))



# 对B cells注释细分
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250829/00.data/recluster/anno_result.h5ad")
# adata_raw = adata_raw.raw.to_adata()
# 'B cells', 'Mast cells', 'Monocytes', 'NK cells', 'Neutrophils', 'Plasmacytoid dendritic cells', 'Platelets', 'T cells'
out_dir="/data/NAS02/sh/Projects/FWSC20250829/01.Mission/B cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="B cells", :]
sc.pp.log1p(adata_sub)
adata_sub.obs['platform'] = adata_sub.obs['platform'].astype(str).astype(object).astype("category")
del adata_sub.obs['RNA_snn_res.0.3']
del adata_sub.obs['seurat_clusters']
adata_sub.write(os.path.join(out_dir, 'result.h5ad'))
# adata_recluster = recluster(adata=adata_sub,species='human',neighbors=30,resolution=0.3,n_pcs = 20,batch_effect = 0,output=out_dir)
# adata = adata_recluster.sub_cluster()



# 对Monocytes注释细分
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250829/00.data/recluster/anno_result.h5ad")
# adata_raw = adata_raw.raw.to_adata()
# 'B cells', 'Mast cells', 'Monocytes', 'NK cells', 'Neutrophils', 'Plasmacytoid dendritic cells', 'Platelets', 'T cells'
out_dir="/data/NAS02/sh/Projects/FWSC20250829/01.Mission/Monocytes/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Monocytes", :]
sc.pp.log1p(adata_sub)
adata_sub.obs['platform'] = adata_sub.obs['platform'].astype(str).astype(object).astype("category")
del adata_sub.obs['RNA_snn_res.0.3']
del adata_sub.obs['seurat_clusters']
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=30,resolution=0.5,n_pcs = 20,batch_effect = 3,output=out_dir)
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

sc.pl.umap(adata, color=['FCGR3A','CD14'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Monocytes_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Monocytes_diffgenes.pdf'))

cell_type_mapping = adata.obs[['leiden']].to_dict()['leiden']
adata_raw.obs['test'] = adata_raw.obs['cell_type']
adata_raw.obs['test'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['test'] = adata_raw.obs['test'].astype('category')
sc.pl.umap(adata_raw, color='test', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_Monocytes_diffgenes.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_Monocytes_diffgenes.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

anno_dict={
    '0':'Classical Monocytes',
    '1':'Classical Monocytes',
    '2':'Classical Monocytes',
    '3':'Non-classical Monocytes',
		'4':'Macrophage',
    '5':'Classical Monocytes',
    '6':'lowQ',
    '7':'DC'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
adata = adata[~adata.obs['leiden'].isin(['6']), :]
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

samples = adata.obs['sample'].drop_duplicates().tolist()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["sample"]==samples[i]], color="cell_type", legend_loc=legend, title=samples[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'Myeloid.h5ad'))



# 映射所有细分
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250829/00.data/recluster/anno_result.h5ad")
out_dir="/data/NAS02/sh/Projects/FWSC20250829/01.Mission/mapping_raw/"
adata_raw.obs['test'] = adata_raw.obs['cell_type']

bdata = sc.read("/data/NAS02/sh/Projects/FWSC20250829/01.Mission/Monocytes/Myeloid.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['test'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['test']), axis=1)

bdata = pd.read_csv("/data/NAS02/sh/Projects/FWSC20250829/01.Mission/B cells/all_cells_features.xls", sep='\t')
cell_type_mapping = bdata[['cell_type']].to_dict()['cell_type']
adata_raw.obs['test'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['test']), axis=1)

bdata = sc.read("/data/NAS02/sh/Projects/FWSC20250829/01.Mission/T cells/Tcells_anno.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['test'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['test']), axis=1)

adata_raw.obs['test'] = adata_raw.obs['test'].astype('category')
adata_raw.obs['test'].cat.categories

sc.pl.umap(adata_raw, color='test', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping_diffgenes.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping_diffgenes.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw = adata_raw[~adata_raw.obs['test'].isin(['Monocytes']), :]
sc.pl.umap(adata_raw, color='test', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping_raw_diffgenes.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping_raw_diffgenes.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.write(os.path.join(out_dir, 'result_mapping_raw.h5ad'))
adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')

sc.pl.umap(adata_raw, color='seurat_clusters', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping_diffgenes.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping_diffgenes.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw = adata_raw[~adata_raw.obs['seurat_clusters'].isin(['11','16']), :]
adata_raw.obs['sub_cell_type'] = adata_raw.obs['test']
del adata_raw.obs['test']
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping_diffgenes.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping_diffgenes.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping_sub_diffgenes.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping_sub_diffgenes.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.write(os.path.join(out_dir, 'result_mapping_sub.h5ad'))
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample')



adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250829/01.Mission/mapping_raw/result_mapping_sub.h5ad")
bdata = pd.read_csv("/mnt/NexData_file/home/sh/Projects/FWSC20250829/01.Mission/mapping_raw/all_cells_features_clean.xls", sep='\t')
common_rows = bdata.index.intersection(adata.obs.index)
adata_filtered = adata[adata.obs.index.isin(common_rows)]
out_dir = "/mnt/NexData_file/home/sh/Projects/FWSC20250829/01.Mission/mapping_raw/result_mapping_clean/"
sc.pl.umap(adata_filtered, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_filtered, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_sub_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_sub_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

groups = adata_filtered.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_filtered[adata_filtered.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

samples = adata_filtered.obs['sample'].drop_duplicates().tolist()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_filtered[adata_filtered.obs["sample"]==samples[i]], color="cell_type", legend_loc=legend, title=samples[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_filtered.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_filtered, out_dir=out_dir,type='group')
plot_celltype_proportion(adata_filtered, out_dir=out_dir,type='sample')
adata_filtered.write(os.path.join(out_dir, 'result_mapping_clean.h5ad'))