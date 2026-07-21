import os, sys
import numpy as np
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import pandas as pd

adata = sc.read_h5ad("/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/Neurons/subcluster.h5ad")

anno_dict={
    '0':'Excitatory neuron',
    '1':'Excitatory neuron',
    '2':'Inhibitory neuron',
    '3':'Excitatory neuron',
    '4':'Other',
    '5':'Other',
    '6':'Other',
    '7':'Other',
    '8':'Other',
    '9':'Inhibitory neuron',
    '10':'Excitatory neuron'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
out_dir="/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_recluster_Excitatory-Inhibitory/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata,
                        groupby='cell_type',
                        method='wilcoxon',
                        pts=True)
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False,show=False)
plt.savefig(os.path.join(out_dir, 'rank_genes_groups.pdf'))
plt.savefig(os.path.join(out_dir, 'rank_genes_groups.png'))

rank_results = adata.uns['rank_genes_groups']
all_groups_results = pd.DataFrame()
 
group_labels = adata.obs[rank_results['params']['groupby']].unique()
for group_label in group_labels:
    group_df = sc.get.rank_genes_groups_df(adata, group=group_label)
    group_df = group_df.sort_values(by="scores", ascending=False)
    group_df = group_df[abs(group_df['logfoldchanges']) >= 0.25]
    group_df = group_df[group_df['pvals'] <= 0.05]
    pts = adata.uns['rank_genes_groups']['pts'][group_label]
    selected_genes = pts[pts >= 0.1].index.tolist()
    group_df = group_df[group_df['names'].isin(selected_genes)]
    group_df['group'] = group_label
    all_groups_results = pd.concat([all_groups_results, group_df], ignore_index=True)
 
all_groups_results.to_csv(os.path.join(out_dir, 'differential_expression_all_cell_type.csv'), index=False)

# cafmarkers={'Bergmann glial cell':['Grid2','Mcc','Pdzrn3'],
#             'Excitatory neuron':['Camk2a','Slc17a7','Camk2a'],
#             'GABAergic neuron':['Gad1','Gad2','Pvalb','Slc32a1','Vip'],
#             'Inhibitory neuron':['Zfpm2'],
#             'Interneuron':['Gad1','Gad2','Bcl11b'],
#             'Neuroblast':['Meis2','Sox5','Tubb3'],
#             }
# sc.pl.dotplot(adata,cafmarkers,groupby='leiden',show=False)
# fig = matplotlib.pyplot.gcf()
# fig.set_size_inches(8, 8)
# plt.savefig(os.path.join(out_dir, 'dotplot_marker_celltype.png'),dpi=300)
# plt.savefig(os.path.join(out_dir, 'dotplot_marker_celltype.pdf'),dpi=300)

# ax = sc.pl.stacked_violin(adata, var_names=cafmarkers, groupby='cell_type', show=False)
# fig = matplotlib.pyplot.gcf()
# fig.set_size_inches(10, 6)
# plt.savefig(os.path.join(out_dir, f'violin_marker_celltype.pdf'),dpi=300)
# plt.savefig(os.path.join(out_dir, f'violin_marker_celltype.png'),dpi=300)
# plt.close()


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



adata = sc.read_h5ad("/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/Neurons/subcluster.h5ad")

anno_dict={
    '0':'Other',
    '1':'Other',
    '2':'Other',
    '3':'Other',
    '4':'Other',
    '5':'GABAergic neuron',
    '6':'GABAergic neuron',
    '7':'Other',
    '8':'Other',
    '9':'Other',
    '10':'Other'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
out_dir="/data/NAS02/sh/Projects/FWAE20242338/05.Mission/result_recluster_GABAergic-Other/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))


sc.tl.rank_genes_groups(adata,
                        groupby='cell_type',
                        method='wilcoxon',
                        pts=True)
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False,show=False)
plt.savefig(os.path.join(out_dir, 'rank_genes_groups.pdf'))
plt.savefig(os.path.join(out_dir, 'rank_genes_groups.png'))

rank_results = adata.uns['rank_genes_groups']
all_groups_results = pd.DataFrame()
 
group_labels = adata.obs[rank_results['params']['groupby']].unique()
for group_label in group_labels:
    group_df = sc.get.rank_genes_groups_df(adata, group=group_label)
    group_df = group_df.sort_values(by="scores", ascending=False)
    group_df = group_df[abs(group_df['logfoldchanges']) >= 0.25]
    group_df = group_df[group_df['pvals'] <= 0.05]
    pts = adata.uns['rank_genes_groups']['pts'][group_label]
    selected_genes = pts[pts >= 0.1].index.tolist()
    group_df = group_df[group_df['names'].isin(selected_genes)]
    group_df['group'] = group_label
    all_groups_results = pd.concat([all_groups_results, group_df], ignore_index=True)
 
all_groups_results.to_csv(os.path.join(out_dir, 'differential_expression_all_cell_type.csv'), index=False)


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