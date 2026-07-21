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

# 对Fibroblast注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2604012626-01/02.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'Epithelial', 'Macrophage', 'Fibroblast', 'T cell', 'Neutrophil', 'B cell', 'Plasma cell'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2604012626-01/08.Mission/Fibroblast/"
# with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["Fibroblast"]), :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

cell_markers = {
    "Endothelial": ["Pecam1", "Kdr", "Cdh5"],
    "Lymphatic endothelial": ["Lyve1", "Prox1", "Pdpn"],
    "Fibroblast": ["Col1a1", "Col1a2", "Dcn", "Lum"],
    "Myofibroblast": ["Acta2", "Tagln", "Myl9"],
    "Smooth muscle": ["Myh11", "Acta2"],
    "Neuronal": ["Plp1", "Sox10", "Elavl3"]
}
marker=[gene for markers in cell_markers.values() for gene in markers]
sc.pl.umap(adata, color=marker, cmap='coolwarm',vmax='p99')
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
    '0':'Bgn Fibroblast',
    '1':'Bgn Fibroblast',
    '2':'Endothelial',
    '3':'Bgn Fibroblast',
    '4':'Slpi Fibroblast',
    '5':'Myofibroblast',
    '6':'Bgn Fibroblast',
    '7':'Lymphatic endothelial',
    '8':'Neuronal'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['NC', 'HLH', 'iMCD-TAFRO']
# adata.obs['group'] = (
#     pd.Categorical(adata.obs['group'],
#                    categories=groups,
#                    ordered=True)
# )
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
adata.write(os.path.join(out_dir, 'Fibroblast.h5ad'))


#### 评分
out_dir="/mnt/NexData_file/home/sh/Projects/SC2604012626-01/08.Mission/genesets_plot/"
genesets = {
    "ISC_identity": ["Lgr5", "Olfm4", "Ascl2", "Smoc2", "Axin2", "Sox9"],
    "TFRC_iron_metabolism": ["Tfrc", "Slc11a2", "Slc40a1", "Fth1", "Ftl1", "Hmox1", "Ncoa4", "Steap3", "Ireb2", "Aco1"],
    "ferroptosis": ["Gpx4", "Slc7a11", "Slc3a2", "Acsl4", "Lpcat3", "Alox15", "Alox5", "Sat1", "Chac1", "Nfe2l2", "Hmox1"],
    "proliferation": ["Mki67", "Top2a", "Pcna", "Mcm2", "Mcm3", "Mcm4", "Mcm5", "Mcm6", "Mcm7", "Ccna2", "Ccnb1"],
    "DNA_damage": ["Trp53", "Cdkn1a", "Gadd45a", "Bax", "Bbc3", "Mdm2", "Ddit4"],
    "oxidative_stress": ["Nqo1", "Gclm", "Gclc", "Hmox1", "Txnrd1", "Srxn1", "Prdx1", "Sod1", "Sod2"],
    "apoptosis": ["Bax", "Bbc3", "Pmaip1", "Casp3", "Casp7"],
    "senescence": ["Cdkn1a", "Cdkn2a", "Cdkn2b", "Serpine1", "Igfbp7", "Glb1"],
    "regeneration": ["Clu", "Anxa1", "Tacstd2", "Ly6a", "Krt19", "Yap1", "Tead1", "Areg"]
}
for set_name, genes in genesets.items():
    sc.tl.score_genes(adata, gene_list=genes, score_name=set_name)

for set_name, genes in genesets.items():
    sc.pl.umap(adata, color=set_name, size=25, show=False,cmap='RdBu_r')
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_{set_name}_score', fig_size={'w': 7, 'h': 6})
    sc.pl.violin(adata, keys=set_name, groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_{set_name}_score_cell_type', fig_size={'w': 10, 'h': 6})
    sc.pl.violin(adata, keys=set_name, groupby='group', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_{set_name}_score_group', fig_size={'w': 10, 'h': 6})
    sc.pl.dotplot(adata, var_names=[set_name], groupby='group', standard_scale='var',  cmap='RdBu_r',  show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_score', fig_size={'w': 6, 'h': 6})


# 映射所有细分 后面没映射
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/mapping_raw/"
# adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/Monocytes/Macrophage.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD4T/CD4T.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD8T/CD8T.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

# adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['lowQ']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')