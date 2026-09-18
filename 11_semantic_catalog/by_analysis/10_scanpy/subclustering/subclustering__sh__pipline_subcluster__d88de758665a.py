import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

# 对巨噬进行细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20250579/result.h5ad")
adata_raw = adata_raw[adata_raw.obs.cell_type=="Macrophages", :]
adata_sub = recluster(adata=adata_raw,species='human',neighbors=20,resolution=0.3,n_pcs = 15,output="/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/")
adata = adata_sub.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/top10_diff_genes_dotplot.png",bbox_inches='tight', dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/top10_diff_genes_dotplot.pdf",bbox_inches='tight')
plt.close()

anno_dict={
    '0':'Macrophages_C1qa',
    '1':'Macrophages_Plac8',
    '2':'Macrophages_H2-Aa',
    '3':'Macrophages_Mki67',
		'4':'Macrophages_Cbfa2t3'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/all_cells_features.xls"), sep='\t')
adata.write('/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/Macrophages.h5ad')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig('/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/dimplot_umap_celltype.png',bbox_inches='tight', dpi=300)
plt.savefig('/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/dimplot_umap_celltype.pdf',bbox_inches='tight')

# 对T进行细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20250579/result.h5ad")
adata_raw = adata_raw[adata_raw.obs.cell_type=="T cells", :]
adata_sub = recluster(adata=adata_raw,species='human',neighbors=15,resolution=0.3,n_pcs = 15,output="/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/")
adata = adata_sub.sub_cluster()

adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/all_cells_features.xls"), sep='\t')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/top10_diff_genes_dotplot.png",bbox_inches='tight', dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/top10_diff_genes_dotplot.pdf",bbox_inches='tight')
plt.close()

# 对肿瘤进行细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20250579/result.h5ad")
adata_raw = adata_raw[adata_raw.obs.cell_type=="tumor-like cell", :]
adata_sub = recluster(adata=adata_raw,species='human',neighbors=15,resolution=0.3,n_pcs = 15,output="/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_tumor-like cell/")
adata = adata_sub.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_tumor-like cell/top10_diff_genes_dotplot.png",bbox_inches='tight', dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_tumor-like cell/top10_diff_genes_dotplot.pdf",bbox_inches='tight')
plt.close()

anno_dict={
    '0':'tumor-like cell_Mki67',
    '1':'tumor-like cell_Col6a2',
    '2':'tumor-like cell_Rps26',
    '3':'tumor-like cell_Fcer1g'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_tumor-like cell/all_cells_features.xls"), sep='\t')
adata.write('/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_tumor-like cell/tumor-like cells.h5ad')

out_dir="/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_tumor-like cell/"
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

sc.tl.rank_genes_groups(adata,groupby='cell_type',method='wilcoxon',pts=True)
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='cell_type',
    key='rank_genes_groups',
    n_genes=5,
    show=False)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_tumor-like cell/top5_diff_genes_dotplot.png",bbox_inches='tight', dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_tumor-like cell/top5_diff_genes_dotplot.pdf",bbox_inches='tight')
plt.close()

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
    plt.savefig(os.path.join(out_dir, 'umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


bdata = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_tumor-like cell/tumor-like cells.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')

out_dir="/data/NAS02/sh/Projects/FWSC20250579/01.Mission/mapping_recluster/"
sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')
adata_raw.write(os.path.join(out_dir, 'mapping_recluster.h5ad'))


# 重新对巨噬和T进行细分
adata = sc.read("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/mapping_recluster/mapping_recluster.h5ad")
bdata = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/Macrophages.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata.obs['cell_type'] = adata.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
bdata = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/result.h5ad")
bdata.obs['cell_type'] = bdata.obs['leiden']
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata.obs['cell_type'] = adata.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

out_dir="/data/NAS02/sh/Projects/FWSC20250579/01.Mission/mapping_recluster/"
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

# 对T进行细分
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/result.h5ad")
adata_raw = adata_raw[~adata_raw.obs['leiden'].isin(['4']), :]
adata_sub = recluster(adata=adata_raw,species='human',neighbors=15,resolution=0.3,n_pcs = 15,output="/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/")
adata = adata_sub.sub_cluster()

anno_dict={
    '0':'CD8T cell',
    '1':'CD4T cell',
    '2':'CD8T cell',
    '3':'Treg'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
adata.write('/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/Tcells.h5ad')

adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/all_cells_features.xls"), sep='\t')

out_dir="/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/"
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

sc.tl.rank_genes_groups(adata,groupby='cell_type',method='wilcoxon',pts=True)
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='cell_type',
    key='rank_genes_groups',
    n_genes=5,
    show=False)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/top5_diff_genes_dotplot.png",bbox_inches='tight', dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/top5_diff_genes_dotplot.pdf",bbox_inches='tight')
plt.close()

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
    plt.savefig(os.path.join(out_dir, 'umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

# 对巨噬进行细分
adata = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/Macrophages.h5ad")

anno_dict={
    '0':'Macrophage',
    '1':'DC',
    '2':'Monocyte',
    '3':'Macrophage',
		'4':'Monocyte'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/all_cells_features.xls"), sep='\t')
adata.write('/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/Macrophages.h5ad')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig('/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/dimplot_umap_celltype.png',bbox_inches='tight', dpi=300)
plt.savefig('/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/dimplot_umap_celltype.pdf',bbox_inches='tight')

sc.tl.rank_genes_groups(adata,groupby='cell_type',method='wilcoxon',pts=True)
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='cell_type',
    key='rank_genes_groups',
    n_genes=5,
    show=False)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/top5_diff_genes_dotplot.png",bbox_inches='tight', dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/top5_diff_genes_dotplot.pdf",bbox_inches='tight')
plt.close()

out_dir="/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/"
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
    plt.savefig(os.path.join(out_dir, 'umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw = sc.read("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/mapping_recluster/mapping_recluster.h5ad")

bdata = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Macrophages/Macrophages.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
bdata = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/Tcells.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
adata_raw = adata_raw[~adata_raw.obs['cell_type'].isin(['T cells']), :]

out_dir="/data/NAS02/sh/Projects/FWSC20250579/01.Mission/mapping_recluster/"
sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')
adata_raw.write(os.path.join(out_dir, 'mapping_recluster.h5ad'))