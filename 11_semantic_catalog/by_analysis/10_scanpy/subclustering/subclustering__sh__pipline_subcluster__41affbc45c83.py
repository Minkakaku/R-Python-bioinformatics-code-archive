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

# 对肾注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/zwl/Project/SC2512227685-01/Kidney_res.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512227685-01/03.Mission/Kidney/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_recluster = recluster(adata=adata_raw,species='mouse',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cell_markers = {
    "POD":['Ptprq','Wt1','Ntng1','Nphs1','Nphs2','Clic5','Podxl'],
    "PEC":['Cldn1','Vcam1','Cfh','Rbfox1','Aldh1a2'],
    "PT":['Lrp2','Cubn','Slc13a1'],
    "TL":['Cryab','Tacstd2','Slc44a5','Klrg2','Col26a1','Boc'],
    "TAL":['Casr','Slc12a1','Umod'],
    "DCT":['Slc12a3','Cnnm2','Fgf13','Klhl3','Trpm6'],
    "CNT":['Slc8a1','Scn2a'],
    "PC":['Gata3','Aqp2','Aqp3'],
    "IC":['Atp6v0d2','Atp6v1c2','Tmem213','Clnk'],
    "EC":['Cd34','Pecam1','Ptprb','Meis2','Flt1','Emcn'],
    "VSM/P":['Notch3','Pdgfrb','Itga8'],
    "FIB":['Col1a1','Col1a2','C7','Negr1','Fbln5','Dcn','Cdh11'],
    "IMM":['Ptprc'],
}
markers=['Ptprq','Wt1','Ntng1','Nphs1','Nphs2','Clic5','Podxl','Cldn1','Vcam1','Cfh','Rbfox1','Aldh1a2','Lrp2','Cubn','Slc13a1','Cryab','Tacstd2','Slc44a5','Klrg2','Col26a1','Boc','Casr','Slc12a1','Umod','Slc12a3','Cnnm2','Fgf13','Klhl3','Trpm6','Slc8a1','Scn2a','Gata3','Aqp2','Aqp3','Atp6v0d2','Atp6v1c2','Tmem213','Clnk','Cd34','Pecam1','Ptprb','Meis2','Flt1','Emcn','Notch3','Pdgfrb','Itga8','Col1a1','Col1a2','C7','Negr1','Fbln5','Dcn','Cdh11','Ptprc']
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
    '0':'Epithelial cells',
    '1':'Epithelial cells',
    '2':'Epithelial cells',
    '3':'Endothelial cells',
    '4':'Epithelial cells',
    '5':'Epithelial cells',
    '6':'Epithelial cells',
    '7':'Epithelial cells',
    '8':'Smooth muscle cells',
    '9':'Epithelial cells',
    '10':'Endothelial cells',
    '11':'Immune cells',
    '12':'Epithelial cells',
    '13':'Epithelial cells',
    '14':'Epithelial cells',
    '15':'Epithelial cells',
    '16':'Epithelial cells',
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Main.h5ad'))

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2512227685-01/06.Mission/Kidney/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512227685-01/07.Mission/"

result = adata_raw.uns['rank_genes_groups']
groups = result['names'].dtype.names
all_diff = []
for c in groups:
    tmp = sc.get.rank_genes_groups_df(adata_raw, group=c)
    tmp = tmp[(tmp['logfoldchanges']>1) & (tmp['pvals_adj']<0.05) & (tmp['pct_nz_group'] >= 0.25)]
    tmp['cluster'] = c
    all_diff.append(tmp)

all_diff = pd.concat(all_diff, axis=0)
all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
all_diff.to_excel(os.path.join(out_dir, '02_cluster_marker_0_16.xlsx'), index=None)
adata_raw.obs[['leiden','sample','group']].value_counts().sort_index().to_excel(os.path.join(out_dir, '03_cluster_proportion_sample.xlsx'))
pd.concat([adata_raw.obs, pd.DataFrame(adata_raw.obsm['X_umap'], index=adata_raw.obs_names, columns=['UMAP1', 'UMAP2'])], axis=1).to_csv(os.path.join(out_dir, "04_cell_metadata.xls"), sep='\t')



# 对肺注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/zwl/Project/SC2512227685-01/03.Mission/Lung_res.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512227685-01/03.Mission/Lung/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_recluster = recluster(adata=adata_raw,species='mouse',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cell_markers = {
    "Epithelial": ["Epcam", "Krt8", "Krt18"],
    "Endothelial": ["Pecam1", "Cldn5", "Vwf"],
    "Fibroblast": ["Dcn", "Col1a1"],
    "T cells": ["Cd3d", "Cd3e", "Cd3g", "Nkg7"],
    "NK cells": ["Tyrobp", "Ncam1"],
    "B cells": ["Ms4a1", "Cd79a", "Cd19"],
    "Plasma cells": ["Tnfrsf17", "Jchain"],
    "Myeloid": ["C1qa", "Itgam"],
    "DC": ["Itgax", "Clec9a", "Lamp3"],
    "Mast cells": ["Ms4a2", "Cpa3"],
    "Neutrophil": ["S100a9", "Csf3r"],
    "APC": ["H2-Aa","H2-Ab1"],  # 小鼠对应 H2-Aa / H2-Ab1 等，保留或替换
    "Pericyte_SMC": ["Rgs5", "Acta2", "Myh11"]
}
markers=["Epcam", "Krt8", "Krt18", "Pecam1", "Cldn5", "Vwf", "Dcn", "Col1a1", "Cd3d", "Cd3e", "Cd3g", "Nkg7", "Tyrobp", "Ncam1", "Ms4a1", "Cd79a", "Cd19", "Tnfrsf17", "Jchain", "C1qa", "Itgam","Itgax", "Clec9a", "Lamp3","Ms4a2", "Cpa3", "S100a9", "Csf3r", "Rgs5", "Acta2", "Myh11"]
sc.pl.umap(adata, color=markers, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

epimarkers=['Itgb4',#Lnep
                 'Calca',#Neuroendocrine 16？ √
                 'Scgb3a1','Muc5b',#Goblet ×
                 'Ccdc153','Foxj1','Cdhr3',#Ciliated √
                 'Scgb1a1',#Club
                 'Krt8',#Krt8+ 17？ √
                 'Sftpd','Sftpc','Sftpa1',#At2 √
                 'Pdpn','Ager','Aqp5']#AT1 √
sc.pl.dotplot(adata, var_names=epimarkers, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_epimarkers.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_epimarkers.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=epimarkers, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_epimarkers.png'))
plt.savefig(os.path.join(out_dir, 'umap_epimarkers.pdf'))

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
    '0':'Epithelial cells',
    '1':'Epithelial cells',
    '2':'Epithelial cells',
    '3':'Epithelial cells',
    '4':'Epithelial cells',
    '5':'Macrophage',
    '6':'Endothelial cells',
    '7':'Smooth muscle cells',
    '8':'T cells'
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Main.h5ad'))

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')