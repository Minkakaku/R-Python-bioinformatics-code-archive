import os, sys
import numpy as np
import scanpy as sc
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import pandas as pd
import seaborn as sns

adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_recluster_GABAergic-Other/subcluster.h5ad")
adata_raw = adata_raw[adata_raw.obs.cell_type=="GABAergic neuron", :]
adata_sub = recluster(adata=adata_raw,species='mouse',neighbors=15,resolution=0.2,n_pcs = 15,output="/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_subcluster/")
adata = adata_sub.sub_cluster()

adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_subcluster/all_cells_features.xls"), sep='\t')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=True
)
plt.savefig("/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_subcluster/top10_diff_genes_dotplot.png", dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_subcluster/top10_diff_genes_dotplot.pdf", dpi=300)
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=5,
    show=True
)
plt.savefig("/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_subcluster/top5_diff_genes_dotplot.png", dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_subcluster/top5_diff_genes_dotplot.pdf", dpi=300)

###后面没做
#2.10位置，后续进行细分
adata = sc.read("/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_subcluster/result.h5ad")
all_diff_genes_df = sc.get.rank_genes_groups_df(adata, group=None)
print(all_diff_genes_df.head())

# cell_types = adata.obs['leiden'].cat.categories.tolist()
# top10_genes_per_cell_type = {}
# for cell_type in cell_types:
#     diff_genes_df = sc.get.rank_genes_groups_df(adata, group=cell_type)
#     top10_genes = diff_genes_df.head(10)['names'].tolist()
#     top10_genes_per_cell_type[cell_type] = top10_genes
#     sc.pl.umap(adata, color=top10_genes, ncols=5,cmap='coolwarm')
#     plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_subcluster_anno/cluster{cell_type}/umap_top10_genes.png')
#     plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_subcluster_anno/cluster{cell_type}/umap_top10_genes.pdf')

anno_dict={
    '0':'PN:Meis2/Bcl11b',
    '1':'IN:Nnat/Nr2f2',
    '2':'IN:Cck/Reln',
    '3':'IN:Calb2/Nxph1',
    '4':'PN:Bcl11b/Foxp1',
    '5':'PN:Tshz1/Pbx3',
    '6':'IN:Nfib/Tcf4',
    '7':'PN:Six3/Foxp1'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
out_dir="/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_subcluster_anno/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
    sc.pl.umap(adata, color='cell_type',legend_loc ='on data', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_withlabel.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_withlabel.pdf'))

# pl = pd.crosstab(adata.obs['group'],adata.obs['cell_type'], normalize='index')
# pl.plot.bar(stacked=True).legend(loc='upper right')
# plt.savefig(os.path.join(out_dir, 'barplot_celltype_group.png'))
# plt.savefig(os.path.join(out_dir, 'barplot_celltype_group.pdf'))

def plot_stacked_bar(adata, out_dir, group_col='group', cell_type_col='cell_type'):
    pl = pd.crosstab(adata.obs[group_col], adata.obs[cell_type_col], normalize='index')
    sns.set(context='talk',style="ticks")
    fig, ax = plt.subplots(figsize=(10, 8))
    pl.plot(kind='bar', stacked=True, ax=ax, width=0.7, legend=False, edgecolor='black', linewidth=0.8)
    ax.margins(x=0.01)
    if len(pl) > 10:
        for bar in ax.patches:
            bar.set_width(0.7 / len(pl))
    ax.set_ylabel('Proportion', fontsize=18, fontweight='bold')
    ax.set_xlabel('')
    ax.tick_params(axis='x', which='major', labelsize=16, width=2, length=6)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{int(y * 100)}%'))
    ax.set_ylim(bottom=-0.02)
    plt.xticks(rotation=45, ha='center')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, title='CellType', loc='center left', bbox_to_anchor=(1, 0.5), 
              frameon=False, handleheight=2, handlelength=2, prop={'size': 14}, markerscale=1.5)
    plt.tight_layout()
    png_path = os.path.join(out_dir, f'barplot_celltype_{group_col}.png')
    pdf_path = os.path.join(out_dir, f'barplot_celltype_{group_col}.pdf')
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)
    proportion_table = pl.applymap(lambda x: f'{x:.2%}')
    table_filename = f'proportion_table_{group_col}.csv'
    table_path = os.path.join(out_dir, table_filename)
    proportion_table.to_csv(table_path)

plot_stacked_bar(adata, out_dir)
plot_stacked_bar(adata, out_dir, group_col='sample')

cafmarkers={'PN:Meis2/Bcl11b':['Meis2','Bcl11b'],
            'IN:Nnat/Nr2f2':['Nnat','Nr2f2'],
            'IN:Cck/Reln':['Cck','Reln'],
            'IN:Calb2/Nxph1':['Calb2','Nxph1'],
            'PN:Bcl11b/Foxp1':['Bcl11b','Foxp1'],
            'PN:Tshz1/Pbx3':['Tshz1','Pbx3'],
            'IN:Nfib/Tcf4':['Nfib','Tcf4'],
            'PN:Six3/Foxp1':['Six3','Foxp1']
            }
sc.pl.dotplot(adata,cafmarkers,groupby='cell_type',show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker_celltype.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_marker_celltype.pdf'),bbox_inches='tight')

samples = adata.obs['group'].drop_duplicates().tolist()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==samples[i]], color="cell_type", legend_loc=legend, title=samples[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
    #plt.savefig(os.path.join(out_dir, 'umap_group_celltype.pdf'))

adata.write(os.path.join(out_dir, 'subcluster.h5ad'))