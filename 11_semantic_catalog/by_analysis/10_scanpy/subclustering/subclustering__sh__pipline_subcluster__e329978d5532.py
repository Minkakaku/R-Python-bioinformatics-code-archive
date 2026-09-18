import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context


# 在大类增加系膜细胞
adata_raw = sc.read_h5ad("/home/lst/platform/scVision4/tmp_out/FWSC20241287/06-Cell_type_annotation/result.h5ad")
anno_dict={
    '12':'Mesangial cell'
}
mask = adata_raw.obs['cluster'].isin(anno_dict.keys())
# 获取新的类别列表
new_categories = list(adata_raw.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata_raw.obs['cell_type'].cat.categories]
# 更新 categories
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.set_categories(new_categories)
adata_raw.obs.loc[mask, 'cell_type'] = adata_raw.obs.loc[mask, 'cluster'].map(anno_dict)
out_dir="/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_recluster/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

samples = adata_raw.obs['group'].drop_duplicates().tolist()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==samples[i]], color="cell_type", legend_loc=legend, title=samples[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
    #plt.savefig(os.path.join(out_dir, 'umap_group_celltype.pdf'))

adata_raw.write(os.path.join(out_dir, 'recluster.h5ad'))


# 对成纤维进行细分
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_recluster/recluster.h5ad")
adata_raw = adata_raw[adata_raw.obs.cell_type=="Fibroblasts", :]
adata_sub = recluster(adata=adata_raw,species='mouse',neighbors=20,resolution=0.2,n_pcs = 15,output="/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Fibroblasts/")
adata = adata_sub.sub_cluster()

adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Fibroblasts/all_cells_features.xls"), sep='\t')
with plt.rc_context({'figure.figsize': (16, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.3)
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby='leiden',
        key='rank_genes_groups',
        n_genes=10,
        show=False,
        ax=ax
    )

plt.savefig("/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Fibroblasts/top10_diff_genes_dotplot.png", dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Fibroblasts/top10_diff_genes_dotplot.pdf")
plt.close()
with plt.rc_context({'figure.figsize': (12, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.3)
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby='leiden',
        key='rank_genes_groups',
        n_genes=5,
        show=False,
        ax=ax
    )

plt.savefig("/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Fibroblasts/top5_diff_genes_dotplot.png", dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Fibroblasts/top5_diff_genes_dotplot.pdf")
plt.close()
#1.20位置，后续进行细分
all_diff_genes_df = sc.get.rank_genes_groups_df(adata, group=None)
print(all_diff_genes_df.head())

cell_types = adata.obs['leiden'].cat.categories.tolist()
top10_genes_per_cell_type = {}
for cell_type in cell_types:
    diff_genes_df = sc.get.rank_genes_groups_df(adata, group=cell_type)
    top10_genes = diff_genes_df.head(10)['names'].tolist()
    top10_genes_per_cell_type[cell_type] = top10_genes
    sc.pl.umap(adata, color=top10_genes, ncols=5,cmap='coolwarm')
    plt.savefig(f'/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Fibroblasts/cluster{cell_type}/umap_top10_genes.png')
    plt.savefig(f'/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Fibroblasts/cluster{cell_type}/umap_top10_genes.pdf')


anno_dict={
    '0':'Fibroblasts_Cfh',
    '1':'Fibroblasts_Slc6a6',
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
out_dir="/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Fibroblasts/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

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

adata.write(os.path.join(out_dir, 'subcluster_Fibroblasts.h5ad'))


# 对巨噬进行细分
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_recluster/recluster.h5ad")
adata_raw = adata_raw[adata_raw.obs.cell_type=="Macrophages", :]
adata_sub = recluster(adata=adata_raw,species='mouse',neighbors=15,resolution=0.5,n_pcs = 30,output="/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Macrophages/")
adata = adata_sub.sub_cluster()

adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Macrophages/all_cells_features.xls"), sep='\t')
with plt.rc_context({'figure.figsize': (18, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.3)
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby='leiden',
        key='rank_genes_groups',
        n_genes=5,
        show=False,
        ax=ax,
        figsize=(18, 4),
        dendrogram=True
    )

plt.savefig("/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Macrophages/top5_diff_genes_dotplot.png", dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Macrophages/top5_diff_genes_dotplot.pdf")
plt.close()
# 进行细分
all_diff_genes_df = sc.get.rank_genes_groups_df(adata, group=None)
print(all_diff_genes_df.head())

cell_types = adata.obs['leiden'].cat.categories.tolist()
top10_genes_per_cell_type = {}
for cell_type in cell_types:
    diff_genes_df = sc.get.rank_genes_groups_df(adata, group=cell_type)
    top10_genes = diff_genes_df.head(10)['names'].tolist()
    top10_genes_per_cell_type[cell_type] = top10_genes
    sc.pl.umap(adata, color=top10_genes, ncols=5,cmap='coolwarm')
    plt.savefig(f'/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Macrophages/cluster{cell_type}/umap_top10_genes.png')
    plt.savefig(f'/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Macrophages/cluster{cell_type}/umap_top10_genes.pdf')


anno_dict={
    '0':'Pro-inflammatory Mac',
    '1':'Immature Mac',
    '2':'Immature Mac',
    '3':'ECM remodeling Mac',
    '4':'Pro-repair Mac',
    '5':'Immature Mac',
    '6':'Immature Mac',
    '7':'DC',
    '8':'Monocyte',
    '9':'Proliferative Mac',
    '10':'Pro-inflammatory Mac'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
out_dir="/data/NAS02/sh/Projects/FWSC20241287/01.Mission/result_subcluster_Macrophages/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

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

adata.write(os.path.join(out_dir, 'subcluster_Macrophages.h5ad'))