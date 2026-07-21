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
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2512227685-01/06.Mission/Kidney/Main.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512227685-01/07.Mission/05_cluster_fine_annotation/"

# cell_markers = {
#     "POD":['Ptprq','Wt1','Ntng1','Nphs1','Nphs2','Clic5','Podxl'],
#     "PEC":['Cldn1','Vcam1','Cfh','Rbfox1','Aldh1a2'],
#     "PT":['Lrp2','Cubn','Slc13a1'],
#     "TL":['Cryab','Tacstd2','Slc44a5','Klrg2','Col26a1','Boc'],
#     "TAL":['Casr','Slc12a1','Umod'],
#     "DCT":['Slc12a3','Cnnm2','Fgf13','Klhl3','Trpm6'],
#     "CNT":['Slc8a1','Scn2a'],
#     "PC":['Gata3','Aqp2','Aqp3'],
#     "IC":['Atp6v0d2','Atp6v1c2','Tmem213','Clnk'],
#     "EC":['Cd34','Pecam1','Ptprb','Meis2','Flt1','Emcn'],
#     "VSM/P":['Notch3','Pdgfrb','Itga8'],
#     "FIB":['Col1a1','Col1a2','C7','Negr1','Fbln5','Dcn','Cdh11'],
#     "IMM":['Ptprc'],
# }
# markers=['Ptprq','Wt1','Ntng1','Nphs1','Nphs2','Clic5','Podxl','Cldn1','Vcam1','Cfh','Rbfox1','Aldh1a2','Lrp2','Cubn','Slc13a1','Cryab','Tacstd2','Slc44a5','Klrg2','Col26a1','Boc','Casr','Slc12a1','Umod','Slc12a3','Cnnm2','Fgf13','Klhl3','Trpm6','Slc8a1','Scn2a','Gata3','Aqp2','Aqp3','Atp6v0d2','Atp6v1c2','Tmem213','Clnk','Cd34','Pecam1','Ptprb','Meis2','Flt1','Emcn','Notch3','Pdgfrb','Itga8','Col1a1','Col1a2','C7','Negr1','Fbln5','Dcn','Cdh11','Ptprc']
# sc.pl.umap(adata, color=markers, cmap='coolwarm',vmax='p99')
# plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
# plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
# sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
# plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
# plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

# sc.pl.rank_genes_groups_dotplot(
#     adata,
#     groupby='leiden',
#     key='rank_genes_groups',
#     n_genes=10,
#     show=False)
# plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
# plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
# plt.close()

anno_dict={
    '0':'PT',
    '1':'PT',
    '2':'PT',
    '3':'EC',
    '4':'PT',
    '5':'PT',
    '6':'TAL',
    '7':'DCT',
    '8':'FIB',
    '9':'PC',
    '10':'EC',
    '11':'IMM',
    '12':'IC',
    '13':'TAL',
    '14':'POD',
    '15':'TAL',
    '16':'PT',
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['ConK', 'Tre40K', 'Tre400K']
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
adata.obs['group'] = (
    pd.Categorical(adata.obs['group'],
                   categories=groups,
                   ordered=True)
)
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Main.h5ad'))

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False, legend_loc='on data')
save_fig(plt, out_dir, file_name='dimplot_umap_celltype_label', fig_size={'w':5, 'h':5})
sc.pl.umap(adata, color='group', title='group', show=False)
save_fig(plt, out_dir, file_name='dimplot_umap_group', fig_size={'w':5, 'h':5})

cell_markers = {
    "DCT":['Slc12a3','Cnnm2','Fgf13','Klhl3','Trpm6'],
    "EC":['Cd34','Pecam1','Ptprb','Meis2','Flt1','Emcn'],
    "FIB":['Col1a1','Col1a2','C7','Negr1','Fbln5','Dcn','Cdh11'],
    "IC":['Atp6v0d2','Atp6v1c2','Tmem213','Clnk'],
    "IMM":['Ptprc'],
    "PC":['Gata3','Aqp2','Aqp3'],
    "POD":['Ptprq','Wt1','Ntng1','Nphs1','Nphs2','Clic5','Podxl'],
    "PT":['Lrp2','Cubn','Slc13a1'],
    "TAL":['Casr','Slc12a1','Umod'],
}
sc.pl.dotplot(adata, var_names=cell_markers, groupby='cell_type', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), bbox_inches="tight")

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'_color_dict.xls'),sep="\t",index=True,index_label='cell_type')

method = "group"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'_color_dict.xls'),sep="\t",index=True,index_label='group')

# sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
# find_group_marker(adata=adata,species='mouse',output=out_dir)

# sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')