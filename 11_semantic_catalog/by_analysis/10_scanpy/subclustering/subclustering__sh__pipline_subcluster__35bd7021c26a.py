import os, sys
import numpy as np
import scanpy as sc
# from scanpy_fun_draw_group import draw_group
# from scanpy_fun_draw_group import getcolor
# from scanpy_fun_find_marker import find_marker_gene
from scanpy_fun_subcluster import recluster
# from scanpy_fun_subcluster import find_group_marker
# from scanpy_fun_subcluster import draw_cell_type
# from scanpy_fun_subcluster import save_fig
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

adata_raw = sc.read_h5ad("/data/outputdata/outputdata-sc/scVision_output/FWAE20242338/result.h5ad")
adata_raw = adata_raw[adata_raw.obs.cell_type=="Neurons", :]
adata_sub = recluster(adata=adata_raw,species='mouse',neighbors=20,resolution=0.2,n_pcs = 15,output="/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/")
adata = adata_sub.sub_cluster()

adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/all_cells_features.xls"), sep='\t')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=True
)
plt.savefig("/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/top10_diff_genes_dotplot.png", dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/top10_diff_genes_dotplot.pdf", dpi=300)
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=5,
    show=True
)
plt.savefig("/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/top5_diff_genes_dotplot.png", dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/top5_diff_genes_dotplot.pdf", dpi=300)

#12.25位置，后续进行细分
all_diff_genes_df = sc.get.rank_genes_groups_df(adata, group=None)
print(all_diff_genes_df.head())

cell_types = adata.obs['leiden'].cat.categories.tolist()
top10_genes_per_cell_type = {}
for cell_type in cell_types:
    diff_genes_df = sc.get.rank_genes_groups_df(adata, group=cell_type)
    top10_genes = diff_genes_df.head(10)['names'].tolist()
    top10_genes_per_cell_type[cell_type] = top10_genes
    sc.pl.umap(adata, color=top10_genes, ncols=5,cmap='coolwarm')
    plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/cluster{cell_type}/umap_top10_genes.png')
    plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/cluster{cell_type}/umap_top10_genes.pdf')

markergenes={'0': ['Rgs6', 'Lingo2', 'B230216N24Rik', 'Ntrk3', 'Camk2a'], '1': ['Il1rapl2', 'Slc24a2', 'St6galnac5', 'Gm49678', 'Pigk'], '2': ['Hs3st4', 'Cdh18', 'Zfpm2', 'Tle4', 'Pdzrn4'], '3': ['Kcnh5', 'Rorb', 'Brinp3', 'R3hdm1', 'Hs6st3'], '4': ['Kcnc2', 'Sox6', 'Kcnmb2', 'Rbms3', 'Ptprm'], '5': ['Meis2', 'Pbx3', 'Pcbp3', 'Ppm1e', 'Dpp6'], '6': ['Adarb2', 'Zfp536', 'Tcf4', 'Dlx6os1', 'Kcnip1'], '7': ['Tafa1', 'Pex5l', 'Ryr3', 'Gm2164', 'Bcl11b'], '8': ['Tafa2', 'Htr2c', 'Brinp3', 'Khdrbs2', 'Pdzrn3'], '9': ['Vwc2l', 'Dcc', 'Adcy2', 'Grik1', 'Olfm3'], '10': ['Pcdh15', 'Ntng1', 'Rora', 'Hs6st3', 'Gm15398']}
for keys, values in markergenes.items():
    print(f"marker genes for {keys}: {values}")

sc.pl.umap(adata, color=['Meg3', 'Ralyl', 'Kcnd2', 'Chrm2', 'Pcdh9','Tcf4','Zfpm2','Cntn4','Lingo2','Cntnap2'], ncols=5,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/cluster8/umap_top20_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/cluster8/umap_top20_genes.pdf')
sc.pl.umap(adata, color=['Bcan','Cdh19','Fabp7','Kcnj10'], ncols=4,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/cluster8/umap_Glial cell_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/cluster8/umap_Glial cell_genes.pdf')
sc.pl.umap(adata, color=['Slc17a7','Snap25','Gad1','Gad2','Cd24a','Dcx','Stmn2'], ncols=4,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/umap_othercell_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/umap_othercell_genes.pdf')

sc.pl.umap(adata, color=['Camk2a','Slc17a7','Camk2a','Gad1','Gad2','Bcl11b','Zfpm2','Meis2','Sox5','Tubb3','Grid2','Mcc','Pdzrn3'], ncols=4,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/umap_othercell_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/umap_othercell_genes.pdf')
sc.pl.umap(adata, color=['Sp9','Stmn1','Stmn2','Tubb3','Eomes','Neurod1','Sox11'], ncols=4,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/umap_Neuroblast_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/umap_Neuroblast_genes.pdf')
sc.pl.umap(adata, color=['Shroom3','Slc1a3','Slc1a4','Sned1','Sparcl1','Stk32a','Tcf7l1','Tnc','Usp6nl','Vcan'], ncols=4,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/umap_Bergmann glial cell_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/umap_Bergmann glial cell_genes.pdf')

markergenes = ['Camk2a','Slc17a7','Camk2a','Gad1','Gad2','Bcl11b','Zfpm2','Meis2','Sox5','Tubb3','Grid2','Mcc','Pdzrn3']
sc.pl.dotplot(adata,markergenes,groupby='leiden',show=False)
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/umap_othercell_genes_dot.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/umap_othercell_genes_dot.pdf')

anno_dict={
    '0':'Excitatory neuron',
    '1':'Excitatory neuron',
    '2':'Inhibitory neuron',
    '3':'Excitatory neuron',
    '4':'Interneuron',
    '5':'Neuroblast',
    '6':'Interneuron',
    '7':'Interneuron',
    '8':'Bergmann glial cell',
    '9':'Inhibitory neuron',
    '10':'Excitatory neuron'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
out_dir="/data/NAS02/sh/Projects/FWAE20242338/01.Mission/result_subcluster/Neurons/"
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

adata.write(os.path.join(out_dir, 'subcluster.h5ad'))