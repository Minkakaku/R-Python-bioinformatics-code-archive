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

# 对Myeloid cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2508294648/result.h5ad")
# 'B cells', 'Cardiomyocytes', 'Cyclin cells', 'Doublets', 'Endothelial cells', 'Erythrocyte', 'Fibroblasts', 'Myeloid cells', 'Neutrophils', 'Pericytes', 'SMC', 'Schwann cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2508294648-01/01.Mission/Myeloid/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Myeloid cells", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
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

sc.pl.umap(adata, color=['Clec9a','Cd209a','Mrc1','Cd68','Cd86','Ccl5','Cd14','Fcgr3'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['Clec9a','Cd209a','Mrc1','Cd68','Cd86','Ccl5','Cd14','Fcgr3'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Clec9a','Xcr1','Ckb','Naga','Ccl22','Ccr7','Fscn1','Cacnb3','Cd209a','Clec10a','Mgl2','Clec12a','Klrk1','Esam','Tbx21','Dtx1'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_DC.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_DC.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Cd80','Cd86','Fcgr1','Fcgr2b','Cd163','Mrc1','Cd68'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_M1_2.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_M1_2.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Macrophage',
    '1':'Monocyte',
    '2':'cDC C1',
    '3':'Monocyte',
		'4':'Gpnmb Macrophage',
    '5':'Macrophage',
    '6':'cDC C2',
    '7':'Macrophage'
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':14, 'h':5})
adata.write(os.path.join(out_dir, 'Myeloid.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)


# 对Fibroblasts注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2508294648/result.h5ad")
# 'B cells', 'Cardiomyocytes', 'Cyclin cells', 'Doublets', 'Endothelial cells', 'Erythrocyte', 'Fibroblasts', 'Myeloid cells', 'Neutrophils', 'Pericytes', 'SMC', 'Schwann cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2508294648-01/01.Mission/Fibroblasts/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Fibroblasts", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['Thbs4','Postn','Clip1','Clip2','Clip3','Clip4'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_gene.png'))
plt.savefig(os.path.join(out_dir, 'umap_gene.pdf'))
sc.pl.dotplot(adata, var_names=['Thbs4','Postn','Clip1','Clip2','Clip3','Clip4'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_gene.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene.pdf'),bbox_inches='tight')

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
    '0':'Casp4 Fib',
    '1':'Pi16 Fib',
    '2':'Casp4 Fib',
    '3':'Thbs4 Fib',
    '4':'Tex14 Fib',
    '5':'Flt1 Fib',
    '6':'Pi16 Fib'
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':14, 'h':5})
adata.write(os.path.join(out_dir, 'Fibroblast.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)
sc.pl.umap(adata, color=['Thbs4'], cmap='coolwarm',vmax='p95')
plt.savefig(os.path.join(out_dir, 'umap_Thbs4.png'))
plt.savefig(os.path.join(out_dir, 'umap_Thbs4.pdf'))

# 对Endothelial cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2508294648/result.h5ad")
# 'B cells', 'Cardiomyocytes', 'Cyclin cells', 'Doublets', 'Endothelial cells', 'Erythrocyte', 'Fibroblasts', 'Myeloid cells', 'Neutrophils', 'Pericytes', 'SMC', 'Schwann cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2508294648-01/01.Mission/Endothelial/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Endothelial cells", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 1)
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

sc.pl.umap(adata, color=['Sema3c','Kit','Gpihbp1','Ccl21a', 'Prox1','Hey1', 'Igfbp3','Cd36', 'Car4','Ackr1','Mki67','Top2a'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['Sema3c','Kit','Gpihbp1','Ccl21a', 'Prox1','Hey1', 'Igfbp3','Cd36', 'Car4','Ackr1','Mki67','Top2a'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'capillary EC3',
    '1':'capillary EC1',
    '2':'capillary EC2',
    '3':'capillary EC3',
    '4':'artery EC',
    '5':'capillary EC4',
    '6':'capillary EC1',
    '7':'vein EC',
    '8':'lymphatic EC'
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':14, 'h':5})
adata.write(os.path.join(out_dir, 'Endothelial.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)


# 对Cardiomyocytes注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2508294648/result.h5ad")
# 'B cells', 'Cardiomyocytes', 'Cyclin cells', 'Doublets', 'Endothelial cells', 'Erythrocyte', 'Fibroblasts', 'Myeloid cells', 'Neutrophils', 'Pericytes', 'SMC', 'Schwann cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2508294648-01/01.Mission/Cardiomyocytes/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Cardiomyocytes", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.2,n_pcs = 15,output=out_dir,batch_effect = 1)
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

anno_dict={
    '0':'Cryab Cardiomyocytes',
    '1':'Sptbn1 Cardiomyocytes'
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':14, 'h':5})
adata.write(os.path.join(out_dir, 'Cardiomyocytes.h5ad'))


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2508294648/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2508294648-01/01.Mission/mapping_raw/"
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2508294648-01/01.Mission/Endothelial/Endothelial.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2508294648-01/01.Mission/Fibroblasts/Fibroblast.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2508294648-01/01.Mission/Myeloid/Myeloid.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['Doublets']), :]
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
sc.tl.umap(adata_raw)
sc.tl.leiden(adata_raw, resolution=0.5)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='on data',legend_fontsize=5)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')

sc.pl.umap(adata_raw, color=['Thbs4','Postn','Clip1','Clip2','Clip3','Clip4','Clu'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_gene.png'))
plt.savefig(os.path.join(out_dir, 'umap_gene.pdf'))
sc.pl.dotplot(adata_raw, var_names=['Thbs4','Postn','Clip1','Clip2','Clip3','Clip4'], groupby='cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata_raw, var_names=['Thbs4','Postn','Clip1','Clip2','Clip3','Clip4'], groupby=['cell_type','group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type-group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene_cell_type-group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata_raw, var_names=['Thbs4','Postn','Clip1','Clip2','Clip3','Clip4'], groupby='sub_cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_gene_sub_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene_sub_cell_type.pdf'),bbox_inches='tight')