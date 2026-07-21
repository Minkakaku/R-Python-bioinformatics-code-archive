import os, sys
import numpy as np
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

adata = sc.read_h5ad("/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/Neurons/subcluster.h5ad")

anno_dict={
    '0':'Excitatory neuron',
    '1':'Excitatory neuron',
    '2':'Inhibitory neuron',
    '3':'Excitatory neuron',
    '4':'Interneuron',
    '5':'Neuroblast',
    '6':'GABAergic neuron',
    '7':'Interneuron',
    '8':'Bergmann glial cell',
    '9':'Inhibitory neuron',
    '10':'Excitatory neuron'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
out_dir="/data/NAS02/sh/Projects/FWAE20242338/04.Mission/result_recluster/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cafmarkers={'Bergmann glial cell':['Grid2','Mcc','Pdzrn3'],
            'Excitatory neuron':['Camk2a','Slc17a7','Camk2a'],
            'GABAergic neuron':['Gad1','Gad2','Pvalb','Slc32a1','Vip'],
            'Inhibitory neuron':['Zfpm2'],
            'Interneuron':['Gad1','Gad2','Bcl11b'],
            'Neuroblast':['Meis2','Sox5','Tubb3'],
            }
sc.pl.dotplot(adata,cafmarkers,groupby='leiden',show=False)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(8, 8)
plt.savefig(os.path.join(out_dir, 'dotplot_marker_celltype.png'),dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_marker_celltype.pdf'),dpi=300)

ax = sc.pl.stacked_violin(adata, var_names=cafmarkers, groupby='cell_type', show=False)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(10, 6)
plt.savefig(os.path.join(out_dir, f'violin_marker_celltype.pdf'),dpi=300)
plt.savefig(os.path.join(out_dir, f'violin_marker_celltype.png'),dpi=300)
plt.close()


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