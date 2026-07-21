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

# 对Monocytes,Plasmacytoid dendritic cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20250829/01.Mission/mapping_raw/result_mapping_clean/result_mapping_clean.h5ad")
# adata_raw = adata_raw.raw.to_adata()
# 'B cells', 'Mast cells', 'Monocytes', 'NK cells', 'Plasmacytoid dendritic cells', 'Platelets', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/Myeloid/"
os.makedirs(out_dir, exist_ok=True)
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['Monocytes','Plasmacytoid dendritic cells']), :]
sc.pp.log1p(adata_sub)
adata_sub.obs['platform'] = adata_sub.obs['platform'].astype(str).astype(object).astype("category")
del adata_sub.obs['RNA_snn_res.0.3']
del adata_sub.obs['seurat_clusters']
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 10,batch_effect = 3,output=out_dir)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250829/01.Mission/Monocytes/Myeloid.h5ad")
cell_type_mapping = bdata.obs[['leiden']].to_dict()['leiden']
adata.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.umap(adata, color=['NFKBIA','RETN','CD302','IFI44'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Monocytes_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Monocytes_diffgenes.pdf'))

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
    '0':'Classical Monocytes NFKBIA',
    '1':'Classical Monocytes CD302',
    '2':'Classical Monocytes RETN',
    '3':'Non-classical Monocytes',
		'4':'Macrophage',
    '5':'Classical Monocytes IFI44',
    '6':'Classical Monocytes RETN',
    '7':'Classical Monocytes NFKBIA',
		'8':'cDCs',
    '9':'pDCs'
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

adata.obs['sample2'] = adata.obs['sample']
adata.obs['sample'] = adata.obs['sample'].astype(str)
adata.obs['sample'] = adata.obs['sample'].replace('TB102-case_1', 'case_3')
adata.obs['sample'] = adata.obs['sample'].replace('TB113-case_2', 'control_3')
adata.obs['sample'] = adata.obs['sample'].replace('TB125-control_1', 'case_4')
adata.obs['sample'] = adata.obs['sample'].replace('TB133-control_2', 'control_4')
adata.obs['sample'] = adata.obs['sample'].replace('case-1', 'case-1')
adata.obs['sample'] = adata.obs['sample'].replace('case-2', 'case-2')
adata.obs['sample'] = adata.obs['sample'].replace('control-1', 'control-1')
adata.obs['sample'] = adata.obs['sample'].replace('control-2', 'control-2')
adata.obs['sample'] = adata.obs['sample'].astype('category')
del adata.obs['sample2']

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

del adata.obs['sub_cell_type']
adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'Myeloid.h5ad'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='platform', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata.obs['cell_type']=adata.obs['cell_type'].cat.set_categories([
    'Classical Monocytes CD302', 'Classical Monocytes IFI44', 'Classical Monocytes NFKBIA', 'Classical Monocytes RETN',
    'Macrophage', 
    'Non-classical Monocytes', 
    'cDCs', 'pDCs'])
adata.uns['cell_type_colors'] = ["#1F77B4","#4E9BD4","#7FB6E1","#A5C8ED",
                                 "#FF7F0E", 
                                 "#9467BD", 
                                 "#2CA02C", 
                                 "#D62728"
                                 ]
adata.obs['group'] = adata.obs['group'].cat.set_categories([
         'control'  , 'case'    ])
adata.uns['group_colors'] = ['#2F6BFF','#E74C3C']

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

umap_coords = adata.obsm['X_umap']
umap_df = pd.DataFrame(umap_coords, columns=['UMAP_1', 'UMAP_2'], index=adata.obs_names)
umap_df.to_csv(os.path.join(out_dir, "all_cells_umap.xls"), sep='\t')


# 对CD4 Tcells注释细分
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/CD4 Tcells/"
os.makedirs(out_dir, exist_ok=True)
adata_sub = adata_raw[adata_raw.obs['sub_cell_type'].isin(['CD4 Cytotoxic Tcells','CD4 Tem','CD4 Treg','CD4 naive Tcells']), :]
sc.pl.umap(adata_sub, color='sub_cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pp.log1p(adata_sub)
adata_sub.obs['platform'] = adata_sub.obs['platform'].astype(str).astype(object).astype("category")
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.3,n_pcs = 10,batch_effect = 3,output=out_dir)
adata = adata_recluster.sub_cluster()

anno_dict={
    '0':'CD4 Tem',
    '1':'CD4 naive Tcells',
		'2':'CD4 naive Tcells',
    '3':'CD4 naive Tcells',
    '4':'CD4 Cytotoxic Tcells',
    '5':'CD4 Tem',
		'6':'CD4 Treg'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata.obs['sample2'] = adata.obs['sample']
adata.obs['sample'] = adata.obs['sample'].astype(str)
adata.obs['sample'] = adata.obs['sample'].replace('TB102-case_1', 'case_3')
adata.obs['sample'] = adata.obs['sample'].replace('TB113-case_2', 'control_3')
adata.obs['sample'] = adata.obs['sample'].replace('TB125-control_1', 'case_4')
adata.obs['sample'] = adata.obs['sample'].replace('TB133-control_2', 'control_4')
adata.obs['sample'] = adata.obs['sample'].replace('case-1', 'case-1')
adata.obs['sample'] = adata.obs['sample'].replace('case-2', 'case-2')
adata.obs['sample'] = adata.obs['sample'].replace('control-1', 'control-1')
adata.obs['sample'] = adata.obs['sample'].replace('control-2', 'control-2')
adata.obs['sample'] = adata.obs['sample'].astype('category')
del adata.obs['sample2']
del adata.obs['sub_cell_type']
del adata.obs['RNA_snn_res.0.3']
del adata.obs['seurat_clusters']

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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'CD4T.h5ad'))

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250829/01.Mission/T cells/Tcells_anno.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/T cells/"
os.makedirs(out_dir, exist_ok=True)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(bdata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
bdata.obs['cell_type'] = bdata.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(bdata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
bdata.write(os.path.join(out_dir, 'Tcells_anno.h5ad'))

sc.pl.umap(adata, color=["TRGC1","TRAC","TRBC1","TRBC2"], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'DNT_umap.png'))
plt.savefig(os.path.join(out_dir, 'DNT_umap.pdf'))


# 对NK cells注释细分
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/NK cells/"
os.makedirs(out_dir, exist_ok=True)
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['NK cells']), :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
sc.pp.log1p(adata_sub)
adata_sub.obs['platform'] = adata_sub.obs['platform'].astype(str).astype(object).astype("category")
del adata_sub.obs['RNA_snn_res.0.3']
del adata_sub.obs['seurat_clusters']
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.2,n_pcs = 10,batch_effect = 3,output=out_dir)
adata = adata_recluster.sub_cluster()


# 04.24注释位置
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['FGFBP2','CD4','CD8A','CD3D','CD3G','CD3E'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'NKT_umap.png'))
plt.savefig(os.path.join(out_dir, 'NKT_umap.pdf'))
sc.pl.umap(adata, color=['NCAM1','CCR7','CXCR3','KLRC1','SELL'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'Immature_umap.png'))
plt.savefig(os.path.join(out_dir, 'Immature_umap.pdf'))
sc.pl.umap(adata, color=['FCGR3A','CX3CR1','CXCR1','KLRG1','EOMES'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'Mature_umap.png'))
plt.savefig(os.path.join(out_dir, 'Mature_umap.pdf'))

anno_dict={
    '0':'Mature NK c1',
    '1':'Mature NK c2',
		'2':'NKT',
    '3':'Immature NK'
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

adata.obs['sample2'] = adata.obs['sample']
adata.obs['sample'] = adata.obs['sample'].astype(str)
adata.obs['sample'] = adata.obs['sample'].replace('TB102-case_1', 'case_3')
adata.obs['sample'] = adata.obs['sample'].replace('TB113-case_2', 'control_3')
adata.obs['sample'] = adata.obs['sample'].replace('TB125-control_1', 'case_4')
adata.obs['sample'] = adata.obs['sample'].replace('TB133-control_2', 'control_4')
adata.obs['sample'] = adata.obs['sample'].replace('case-1', 'case-1')
adata.obs['sample'] = adata.obs['sample'].replace('case-2', 'case-2')
adata.obs['sample'] = adata.obs['sample'].replace('control-1', 'control-1')
adata.obs['sample'] = adata.obs['sample'].replace('control-2', 'control-2')
adata.obs['sample'] = adata.obs['sample'].astype('category')
del adata.obs['sample2']
del adata.obs['sub_cell_type']

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
adata.write(os.path.join(out_dir, 'NKcells_anno.h5ad'))


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20250829/01.Mission/mapping_raw/result_mapping_clean/result_mapping_clean.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/mapping_raw/"
os.makedirs(out_dir, exist_ok=True)

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/Myeloid/Myeloid.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping_diffgenes.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping_diffgenes.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.obs['sample2'] = adata_raw.obs['sample']
adata_raw.obs['sample'] = adata_raw.obs['sample'].astype(str)
adata_raw.obs['sample'] = adata_raw.obs['sample'].replace('TB102-case_1', 'case_3')
adata_raw.obs['sample'] = adata_raw.obs['sample'].replace('TB113-case_2', 'control_3')
adata_raw.obs['sample'] = adata_raw.obs['sample'].replace('TB125-control_1', 'case_4')
adata_raw.obs['sample'] = adata_raw.obs['sample'].replace('TB133-control_2', 'control_4')
adata_raw.obs['sample'] = adata_raw.obs['sample'].replace('case-1', 'case-1')
adata_raw.obs['sample'] = adata_raw.obs['sample'].replace('case-2', 'case-2')
adata_raw.obs['sample'] = adata_raw.obs['sample'].replace('control-1', 'control-1')
adata_raw.obs['sample'] = adata_raw.obs['sample'].replace('control-2', 'control-2')
adata_raw.obs['sample'] = adata_raw.obs['sample'].astype('category')
del adata_raw.obs['sample2']
del adata_raw.obs['RNA_snn_res.0.3']
del adata_raw.obs['seurat_clusters']

adata_raw.write(os.path.join(out_dir, 'result_mapping_clean.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')

sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/mapping_raw/result_mapping_clean.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250829/03.Mission/mapping_raw/"
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample')