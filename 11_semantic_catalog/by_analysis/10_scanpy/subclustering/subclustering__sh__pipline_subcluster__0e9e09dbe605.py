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

adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/AE2507033447-07/result.h5ad")
adata_raw =  adata_raw[~adata_raw.obs['sample'].isin(['ADC3', 'SCC5']),].copy()
adata=sc.read("/mnt/NexData_file/home/zwl/Project/AE2507033447-07/02.Mission/Epithelial/Epithelial.h5ad")
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/zwl/Project/AE2507033447-07/02.Mission/Fibroblasts/Fibroblasts.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/04.Mission/"
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')



# 对Epithelial cells注释细分 #### 后面没搞
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/AE2507033447-07/result.h5ad")
# adata_raw.obs.cell_type.cat.categories
adata_sub = adata_raw[adata_raw.obs.cell_type=="Epithelial cells", :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/01.Mission/Epithelial/"
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=1,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['KRT6A','POSTN','MKI67','EGLN3','SPRR3','C15orf48','SPRR2A','LINC01980','LINC02432','IGHM','AKR1C3','SPINT2','MGST1','CXCL10','GBP1','PARP14'], size=25, show=False,cmap='RdBu_r',sort_order=True, vmax='p99')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_marker.png'))
plt.savefig(os.path.join(out_dir, 'dimplot_umap_marker.pdf'))

genes=['APOH', 'APP', 'CCND2', 'COL3A1', 'COL5A2', 'CXCL6', 'FGFR1', 'FSTL1', 'ITGAV', 'JAG1', 'JAG2', 'KCNJ8', 'LPL', 'LRPAP1', 'LUM', 'MSX1', 'NRP1', 'OLR1', 'PDGFA', 'PF4', 'PGLYRP1', 'POSTN', 'PRG2', 'PTK2', 'S100A4', 'SERPINA5', 'SLCO2A1', 'SPP1', 'STC1', 'THBD', 'TIMP1', 'TNFRSF21', 'VAV2' ,'VCAN', 'VEGFA', 'VTN']
sc.tl.score_genes(adata, gene_list=genes, score_name='ANGIOGENESIS')

anno_dict={
    '0':'KRT5 Epi',
    '1':'PHGR1 Epi',
    '2':'PHGR1 Epi',
    '3':'KRT5 Epi',
    '4':'SPINK1 Epi',
    '5':'WFDC2 Epi',
    '6':'KRT6A Epi',
    '7':'MUC16 Epi',
    '8':'FABP5 Epi',
    '9':'SPRR3 Epi',
    '10':'FABP5 Epi',
    '11':'MUC16 Epi',
    '12':'FABP5 Epi',
    '13':'MUC16 Epi'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.umap(adata, color='ANGIOGENESIS', size=25, show=False,cmap='RdBu_r')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_ANGIOGENESIS', fig_size={'w': 7, 'h': 6})
sc.pl.violin(adata, keys='ANGIOGENESIS', groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_ANGIOGENESIS', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['ANGIOGENESIS'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_ANGIOGENESIS', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['ANGIOGENESIS'], groupby=['cell_type','group'],  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_ANGIOGENESIS_group', fig_size={'w': 6, 'h': 12})

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
adata_raw.obs['cell_type'].cat.categories
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':30, 'h':5})
adata.write(os.path.join(out_dir, 'Epithelial.h5ad'))

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

sc.pl.matrixplot(adata, var_names=genes, groupby='cell_type', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_cell_type.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_cell_type.pdf'),bbox_inches='tight')
plt.close()
sc.pl.matrixplot(adata, var_names=genes, groupby='group', standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_group.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_group.pdf'),bbox_inches='tight')
plt.close()
sc.pl.matrixplot(adata, var_names=genes, groupby=['cell_type','group'], standard_scale='var', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_cell_type-group.png'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'ANGIOGENESIS_heatmap_cell_type-group.pdf'),bbox_inches='tight')
plt.close()

# out_dir="/mnt/NexData_file/home/sh/Projects/SC2507033447-36/01.Mission/AAE_gene/"
# adata = adata[adata.obs.group=='ADC', :]
# sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
# find_group_marker(adata=adata,species='human',output=out_dir)

# marker_genes=[]
# sc.pl.stacked_violin(adata, var_names=marker_genes, groupby='cell_type', show=False, swap_axes=True, dendrogram=False)
# plt.savefig(os.path.join(out_dir, 'AAEs_violin.pdf'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'AAEs_violin.png'),bbox_inches='tight')


# 对Fibroblasts注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/AE2507033447-07/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Macrophages', 'Mast cells', 'Neutrophil', 'Plasma cells', 'Smooth muscle cells', 'T cells'
adata_sub = adata_raw[adata_raw.obs.cell_type=="Fibroblasts", :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-07/01.Mission/Fibroblasts/"
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=1,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

sc.pl.umap(adata, color=['ITGA8','POSTN','FAP','MMP2','S100A4','PDGFRA','COL12A1','CD74','HLA-DPA1','HLA-DRA','HLA-DQA1','SLPI','CCL2','TNF','CD34','PLA2G2A','CFD','CD248','MMP11','CTHRC1','FN1','LUM','ACTA2', 'MCAM', 'RGS5','PDPN', 'MME', 'ENO1', 'NDRG1', 'CA9','HSPH1','IDO1', 'COL1A1','CCL21','CCL19','MKI67','COL1A2','PI16','CXCL8','LRRC15','CD37','CD24','ICAM1','PLAT','HGF','LIF','DLG2'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['ITGA8','POSTN','FAP','MMP2','S100A4','PDGFRA','COL12A1','CD74','HLA-DPA1','HLA-DRA','HLA-DQA1','SLPI','CCL2','TNF','CD34','PLA2G2A','CFD','CD248','MMP11','CTHRC1','FN1','LUM','ACTA2', 'MCAM', 'RGS5','PDPN', 'MME', 'ENO1', 'NDRG1', 'CA9','HSPH1','IDO1', 'COL1A1','CCL21','CCL19','MKI67','COL1A2','PI16','CXCL8','LRRC15','CD37','CD24','ICAM1','PLAT','HGF','LIF','DLG2'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

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

anno_dict={
    '0':'Tumour-like CAF',
    '1':'Tumour-like CAF',
    '2':'Matrix CAF',
    '3':'Inflammatory CAF',
    '4':'Neuronal CAF',
    '5':'Tumour-like CAF',
    '6':'Inflammatory CAF',
    '7':'Antigen-presenting CAF',
    '8':'Inflammatory CAF'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':30, 'h':5})
adata.write(os.path.join(out_dir, 'Fibroblasts.h5ad'))

sc.pl.umap(adata, color=['EPHA3'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_EPHA3.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_EPHA3.pdf'), bbox_inches="tight")
sc.pl.dotplot(adata, var_names=['EPHA3'], groupby='cell_type', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_EPHA3.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_EPHA3.pdf'), bbox_inches="tight")


# 对Endothelial cells注释细分 #### 后面都没细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2507033447-01/8vs8/01.Mission/Main/Main.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'Plasma cells', 'Macrophages', 'Epithelial cells', 'T/NK cells','Mast cells', 'Fibroblasts', 'Endothelial cells', 'Neutrophils', 'Smooth muscle cells', 'B cells'
adata_sub = adata_raw[adata_raw.obs.cell_type=="Endothelial cells", :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2507033447-01/8vs8/01.Mission/Endothelial/"
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.8,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

sc.pl.umap(adata, color=['SEMA3C','KIT','GPIHBP1','CCL21', 'PROX1','HEY1', 'IGFBP3','CD36', 'CA4','ACKR1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['SEMA3C','KIT','GPIHBP1','CCL21', 'PROX1','HEY1', 'IGFBP3','CD36', 'CA4','ACKR1'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

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

anno_dict={
    '0':'CLU vein EC',
    '1':'COL4A1 vein EC',
    '2':'artery EC',
    '3':'CLU vein EC',
    '4':'capillary EC',
    '5':'CLU vein EC',
    '6':'CLU vein EC',
    '7':'capillary EC',
    '8':'COL4A1 vein EC',
    '9':'lymphatic EC',
    '10':'COL4A1 vein EC'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':30, 'h':5})
adata.write(os.path.join(out_dir, 'Endothelial.h5ad'))


# 映射所有细分 #### 后面还没映射
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/mapping_raw/"
# adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Bcells/Bcells.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/CD4T/CD4T.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/CD8T/CD8T.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Endothelial/EC.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Epithelial cells/Epithelial cells.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Fibroblasts/Fibroblasts.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Plasma/Plasma.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['lowQ']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')