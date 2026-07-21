import os, sys
import numpy as np
import scanpy as sc
from scanpy_fun_draw_group import draw_group
from scanpy_fun_draw_group import getcolor
from scanpy_fun_find_marker import find_marker_gene
from scanpy_fun_subcluster import recluster
from scanpy_fun_subcluster import find_group_marker
from scanpy_fun_subcluster import draw_cell_type
from scanpy_fun_subcluster import save_fig
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

adata_raw = sc.read_h5ad("/home/yrc/FWSC20242077/result.h5ad")
adata_raw = adata_raw[adata_raw.obs.cell_type=="T cells", :]
adata_sub = recluster(adata=adata_raw,species='mouse',neighbors=20,resolution=0.5,n_pcs = 20,output="/home/yrc/FWSC20242077/T cells subcluster/")
adata = adata_sub.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=5,
    show=True
)
plt.savefig("/home/yrc/FWSC20242077/T cells subcluster/top5_diff_genes_dotplot.png", bbox_inches="tight", dpi=300)
plt.savefig("/home/yrc/FWSC20242077/T cells subcluster/top5_diff_genes_dotplot.pdf", bbox_inches="tight", dpi=300)


markergenes = ['Cd8a', 'Cd4', 'Foxp3','Ncam1','Nkg7','Gzma']

sc.pl.umap(adata, color=markergenes, ncols=3,cmap='coolwarm')
plt.savefig(f'/home/yrc/FWSC20242077/T cells subcluster/marker_feature.png')
plt.savefig(f'/home/yrc/FWSC20242077/T cells subcluster/marker_feature.pdf')

sc.pl.dotplot(adata,markergenes,groupby='leiden',show=False)
plt.savefig(f'/home/yrc/FWSC20242077/T cells subcluster/marker_dot.png')
plt.savefig(f'/home/yrc/FWSC20242077/T cells subcluster/marker_dot.pdf')

anno_dict={
    '0':'CD8 T',
    '1':'CD4 T',
    '2':'CD8 T',
    '3':'Treg',
    '4':'NK'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
out_dir="/home/yrc/FWSC20242077/T cells subcluster/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cafmarkers={'Excitatory neuron':['Camk2a','Slc17a7','Camk2a'],
            'Inhibitory neuron':['Zfpm2'],
            'Interneuron':['Gad1','Gad2','Bcl11b'],
            'Neuroblast':['Meis2','Sox5','Tubb3'],
            'Bergmann glial cell':['Grid2','Mcc','Pdzrn3']
            }
sc.pl.dotplot(adata,cafmarkers,groupby='leiden',show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_marker_celltype.png'))
plt.savefig(os.path.join(out_dir, 'dotplot_marker_celltype.pdf'))

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

adata.write(os.path.join(out_dir, 'T.h5ad'))

