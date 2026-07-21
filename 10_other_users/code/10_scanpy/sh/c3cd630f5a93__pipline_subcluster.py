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

# 对Epithelial cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Plasma cells', 'Schwann cells', 'Smooth muscle cells', 'T/NK', 'lowQ'
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Epithelial cells/"
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/Epithelial cells/Epithelial cells.h5ad")

anno_dict={
    '0':'Stem cells',
    '1':'Stem cells',
    '2':'Tumor-like cells',
    '3':'Tumor-like cells',
    '4':'Tumor-like cells',
    '5':'Tumor-like cells',
    '6':'Tumor-like cells',
    '7':'Goblet cells',
    '8':'Colonocyte',
    '9':'DCS cells',
    '10':'Stem cells',
    '11':'Colonocyte',
    '12':'DCS cells',
    '13':'Stem cells',
    '14':'Tumor-like cells',
    '15':'DCS cells',
    '16':'Tumor-like cells',
    '17':'Tuft cells',
    '18':'TA cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_celltype', fig_size={'w':6, 'h':5})
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':15, 'h':5})

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

sc.pl.dotplot(adata, var_names=['ALPI','SULT1B1','MEP1A','REG4','TFF3','SMIM14','SPDEF','B3GNT6','HMGCS2','MKI67','TOP2A','CENPF','COL1A1','COL1A2','COL3A1','RNF43','MLXIPL','AC004687.2'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'),bbox_inches='tight')

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Epithelial cells.h5ad'))

groups = ['LARC', 'N']
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

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':22, 'h':5})

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 对T/NK注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Plasma cells', 'Schwann cells', 'Smooth muscle cells', 'T/NK', 'lowQ'
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/T_NK/"
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/T_NK/T.h5ad")

sc.pl.umap(adata, color=['CD4','TRDC','FOXP3','CTLA4','S1PR4','SELL','LEF1','TCF7','CCR7','CD69','NR4A1','KLF2','TOB1','RGS1','ANXA1','TBX21','TNF','IFNG','RORC','IL17A','IL23R','CXCR5','BCL6','BATF','CD8A','CD8B','CD27','CD28','GZMA','GZMH','GZMK','FGFBP2','FCGR3A','KLRG1','NKG7','PRF1','HAVCR2','PDCD1','TOX','CXCR3','ENTPD1','MKI67','KLRD1','XCL1','GNLY','NCAM1','NCR1'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['CD4','TRDC','FOXP3','CTLA4','S1PR4','SELL','LEF1','TCF7','CCR7','CD69','NR4A1','KLF2','TOB1','RGS1','ANXA1','TBX21','TNF','IFNG','RORC','IL17A','IL23R','CXCR5','BCL6','BATF','CD8A','CD8B','CD27','CD28','GZMA','GZMH','GZMK','FGFBP2','FCGR3A','KLRG1','NKG7','PRF1','HAVCR2','PDCD1','TOX','CXCR3','ENTPD1','MKI67','KLRD1','XCL1','GNLY','NCAM1','NCR1'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'), bbox_inches="tight")


# 对CD4T注释细分
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/CD4T/"
adata_sub= adata[adata.obs.leiden.isin(['0','1','2','7','8','9']), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['CD4','TRDC','FOXP3','CTLA4','S1PR4','SELL','LEF1','TCF7','CCR7','CD69','NR4A1','KLF2','TOB1','RGS1','ANXA1','TBX21','TNF','IFNG','RORC','IL17A','IL23R','CXCR5','BCL6','BATF'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['CD4','TRDC','FOXP3','CTLA4','S1PR4','SELL','LEF1','TCF7','CCR7','CD69','NR4A1','KLF2','TOB1','RGS1','ANXA1','TBX21','TNF','IFNG','RORC','IL17A','IL23R','CXCR5','BCL6','BATF'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'), bbox_inches="tight")

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
    '0':'ANXA1 CD4T memory',
    '1':'LAIR2 CD4T memory',
    '2':'CD4T naive',
    '3':'Treg',
    '4':'CD4T naive',
    '5':'Th17',
    '6':'CD4T naive'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':2, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'CD4T.h5ad'))

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':15, 'h':5})

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

sc.pl.dotplot(adata, var_names=['ANXA1','CD69','NR4A1','LEF1','TCF7','CCR7','RGS1','TOB1','LAIR2','RORC','IL17A','IL23R','FOXP3','CTLA4','S1PR4'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'),bbox_inches='tight')

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


# 对CD8T注释细分
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/CD8T/"
adata_sub= adata[adata.obs.leiden.isin(['3','4','5','6']), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['CD8A','CD8B','SELL','LEF1','TCF7','CD27','CD28','GZMA','GZMH','GZMK','FGFBP2','FCGR3A','KLRG1','NKG7','PRF1','HAVCR2','PDCD1','TOX','CXCR3','ENTPD1','MKI67','KLRD1','XCL1','GNLY','NCAM1','NCR1'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['CD8A','CD8B','SELL','LEF1','TCF7','CD27','CD28','GZMA','GZMH','GZMK','FGFBP2','FCGR3A','KLRG1','NKG7','PRF1','HAVCR2','PDCD1','TOX','CXCR3','ENTPD1','MKI67','KLRD1','XCL1','GNLY','NCAM1','NCR1'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'), bbox_inches="tight")

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
    '0':'terminally exhausted CD8T cells',
    '1':'CD8T exhaust',
    '2':'CD8T cytotoxic',
    '3':'CD8T naive',
    '4':'precursor exhausted CD8T cells',
    '5':'precursor exhausted CD8T cells',
    '6':'CD8T memory',
    '7':'CD8T effector',
    '8':'precursor exhausted CD8T cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':2, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'CD8T.h5ad'))

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data',legend_fontsize=8)
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.pl.dotplot(adata, var_names=['NKG7','PRF1','IFNG','GZMH','FGFBP2','FCGR3A','HAVCR2','PDCD1','TOX','CD27','CD28','TCF7','SELL','LEF1','GZMK','CXCR3','ENTPD1'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


# 对B cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Plasma cells', 'Schwann cells', 'Smooth muscle cells', 'T/NK', 'lowQ'
adata_sub = adata_raw[adata_raw.obs.cell_type=="B cells", :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Bcells/"
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['CD19','IL2RA','TNFRSF8','CD19','MS4A1','CD19','TCL1A','IL4R','CD27','AIM2','TNFRSF13B','CD80','PDCD1LG2','CXCR3','CD22','IGHD','CR2','FCER2','CD5','CD24','TLR4','TGFB1','CD19','CD27','NOTCH2','CR2','MME','AICDA','SUGCT','MKI67'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['CD19','IL2RA','TNFRSF8','CD19','MS4A1','CD19','TCL1A','IL4R','CD27','AIM2','TNFRSF13B','CD80','PDCD1LG2','CXCR3','CD22','IGHD','CR2','FCER2','CD5','CD24','TLR4','TGFB1','CD19','CD27','NOTCH2','CR2','MME','AICDA','SUGCT','MKI67'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'), bbox_inches="tight")

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
    '0':'Follicular B',
    '1':'naive B',
    '2':'memory B',
    '3':'memory B',
    '4':'B reg',
    '5':'Germinal center B'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':2, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'Bcells.h5ad'))

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

sc.pl.dotplot(adata, var_names=['TGFB1','CD5','IGHD','FCER2','CR2','MME','AICDA','SUGCT','AIM2','TNFRSF13B','CD80','CD19','TCL1A','IL4R'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'),bbox_inches='tight')

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


# 对Plasma cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Plasma cells', 'Schwann cells', 'Smooth muscle cells', 'T/NK', 'lowQ'
adata_sub = adata_raw[adata_raw.obs.cell_type=="Plasma cells", :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Plasma/"
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['IGHA1','IGHG1','IGHG3','IGHG4','IRF4','PRDM1','CD38'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['IGHA1','IGHG1','IGHG3','IGHG4','IRF4','PRDM1','CD38'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'), bbox_inches="tight")

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
    '0':'Plasmablast',
    '1':'IgA Plasma',
    '2':'IgG Plasma',
    '3':'IgA Plasma',
    '4':'IgA Plasma',
    '5':'IgA Plasma',
    '6':'IgA Plasma',
    '7':'IgA Plasma'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':2, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'Plasma.h5ad'))

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

sc.pl.dotplot(adata, var_names=['IGHA1','IGHG1','IGHG3','IGHG4','IRF4','PRDM1','CD38'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'),bbox_inches='tight')

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


# 对Endothelial cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Plasma cells', 'Schwann cells', 'Smooth muscle cells', 'T/NK', 'lowQ'
adata_sub = adata_raw[adata_raw.obs.cell_type=="Endothelial cells", :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Endothelial/"
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['SEMA3C','KIT','GPIHBP1','CCL21', 'PROX1','HEY1', 'IGFBP3','CD36', 'CA4','ACKR1','SOCS3','HLA-DRA','GADD45B','PLVAP','IGFBP3','RGCC','COL4A1','COL4A2','PLPP3'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['SEMA3C','KIT','GPIHBP1','CCL21', 'PROX1','HEY1', 'IGFBP3','CD36', 'CA4','ACKR1','SOCS3','HLA-DRA','GADD45B','PLVAP','IGFBP3','RGCC','COL4A1','COL4A2','PLPP3'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'), bbox_inches="tight")

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
    '0':'PLVAP EC',
    '1':'HLA-DRA EC',
    '2':'PLVAP EC',
    '3':'PLVAP EC',
    '4':'HLA-DRA EC',
    '5':'PLVAP EC',
    '6':'lymphatic EC',
    '7':'PLVAP EC',
    '8':'PLVAP EC',
    '9':'PLVAP EC'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':2, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'EC.h5ad'))

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

sc.pl.dotplot(adata, var_names=['PECAM1','ACKR1','SOCS3','HLA-DRA','GADD45B','PLVAP','IGFBP3','RGCC','COL4A1','COL4A2','PLPP3','CCL21', 'PROX1'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'),bbox_inches='tight')

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


# 对Fibroblasts注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Monocytes', 'Plasma cells', 'Schwann cells', 'Smooth muscle cells', 'T/NK', 'lowQ'
adata_sub = adata_raw[adata_raw.obs.cell_type=="Fibroblasts", :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Fibroblasts/"
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['ITGA8','POSTN','FAP','MMP2','S100A4','PDGFRA','COL12A1','CD74','HLA-DPA1','HLA-DRA','HLA-DQA1','SLPI','CCL2','TNF','CD34','PLA2G2A','CFD','CD248','POSTN','MMP11','CTHRC1','FN1','LUM','ACTA2', 'MCAM', 'RGS5','PDPN', 'MME', 'ENO1', 'NDRG1', 'CA9','HSPH1','IDO1', 'COL1A1','CCL21','CCL19','MKI67'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['ITGA8','POSTN','FAP','MMP2','S100A4','PDGFRA','COL12A1','CD74','HLA-DPA1','HLA-DRA','HLA-DQA1','SLPI','CCL2','TNF','CD34','PLA2G2A','CFD','CD248','POSTN','MMP11','CTHRC1','FN1','LUM','ACTA2', 'MCAM', 'RGS5','PDPN', 'MME', 'ENO1', 'NDRG1', 'CA9','HSPH1','IDO1', 'COL1A1','CCL21','CCL19','MKI67'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'), bbox_inches="tight")

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
    '0':'Interferon response CAFs',
    '1':'Inflammatory CAF',
    '2':'matrix CAF',
    '3':'matrix CAF',
    '4':'Tumour-like CAFs',
    '5':'Vascular CAFs',
    '6':'Vascular CAFs',
    '7':'Antigen presenting CAFs',
    '8':'Inflammatory CAF'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':27, 'h':5})
adata.write(os.path.join(out_dir, 'Fibroblasts.h5ad'))

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

sc.pl.dotplot(adata, var_names=['CD74','HLA-DPA1','HLA-DQA1','CD34','PLA2G2A','CFD','COL1A1','IDO1','PDPN','ENO1','NDRG1','MCAM','RGS5','ACTA2','POSTN','FN1','LUM'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'),bbox_inches='tight')

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


# 对Monocytes注释细分
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/01.Mission/Monocytes/Macrophage.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/Monocytes/"

sc.pl.umap(adata, color=['CD14','FCGR3A','NOS2','CD86','CCL5','IL1B','CD274','MRC1','CD163','CCND1','CKB','GSN','CXCL9','SPP1','VEGFA','APOE','FN1','S100A9','CD1C','CLEC4C','CD209','CLEC9A','LILRA4','TCF4','BATF3','XCR1','ITGAM','FCER1A','CLEC10A','CCL22','LAMP3','FSCN1','CCL19','CCR7'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['CD14','FCGR3A','NOS2','CD86','CCL5','IL1B','CD274','MRC1','CD163','CCND1','CKB','GSN','CXCL9','SPP1','VEGFA','APOE','FN1','S100A9','CD1C','CLEC4C','CD209','CLEC9A','LILRA4','TCF4','BATF3','XCR1','ITGAM','FCER1A','CLEC10A','CCL22','LAMP3','FSCN1','CCL19','CCR7'], groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'), bbox_inches="tight")

anno_dict={
    '0':'M2 Macro',
    '1':'M1 Macro',
    '2':'Infiltrating macrophages',
    '3':'DC',
    '4':'mregDC',
    '5':'DC'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':27, 'h':5})
adata.write(os.path.join(out_dir, 'Macrophage.h5ad'))

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='sample', title='sample', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_sample', fig_size={'w':8, 'h':5})

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

sc.pl.dotplot(adata, var_names=['CD1C','FCER1A','CLEC10A','S100A9','FN1','CCL5','MRC1','CD163','CCND1','CCL22','LAMP3','FSCN1'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'),bbox_inches='tight')

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')


# 映射所有细分
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