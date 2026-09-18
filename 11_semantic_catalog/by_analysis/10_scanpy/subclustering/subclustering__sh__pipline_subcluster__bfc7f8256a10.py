import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

# 对T行细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20250070/result.h5ad")
# adata_raw = adata_raw.raw.to_adata()
# 'B cells', 'Distal tubule cells', 'Endothelial cells','Epithelial cells', 'Macrophages', 'Neutrophils', 'Proximal tubule cells','Smooth muscle cells', 'T cells'
adata_sub = adata_raw[adata_raw.obs.cell_type=="T cells", :]
adata_sub = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.2,n_pcs = 15,output="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/T cells")
adata = adata_sub.sub_cluster()

sc.pl.umap(adata, color=['Cd4', 'Cd8a', 'Mki67', 'Gzma','Gzmb'], ncols=5,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/T cells/umap_marker_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/T cells/umap_marker_genes.pdf')

adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWSC20250070/01.Mission/T cells/all_cells_features.xls"), sep='\t')
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

plt.savefig("/data/NAS02/sh/Projects/FWSC20250070/01.Mission/T cells/top5_diff_genes_dotplot.png", dpi=300)
plt.savefig("/data/NAS02/sh/Projects/FWSC20250070/01.Mission/T cells/top5_diff_genes_dotplot.pdf")
plt.close()

anno_dict={
    '0':'T cells',
    '1':'T cells',
    '2':'NK cells',
    '3':'T cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
out_dir="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/T cells/"
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

adata.write(os.path.join(out_dir, 'subcluster_T cells.h5ad'))


# 对B cells行细分：分不了
# 'B cells', 'Distal tubule cells', 'Endothelial cells','Epithelial cells', 'Macrophages', 'Neutrophils', 'Proximal tubule cells','Smooth muscle cells', 'T cells'
adata_sub = adata_raw[adata_raw.obs.cell_type=="B cells", :]
adata_sub = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.2,n_pcs = 15,output="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/B cells")
adata = adata_sub.sub_cluster()


# 对Macrophages行细分：分不了
# 'B cells', 'Distal tubule cells', 'Endothelial cells','Epithelial cells', 'Macrophages', 'Neutrophils', 'Proximal tubule cells','Smooth muscle cells', 'T cells'
adata_sub = adata_raw[adata_raw.obs.cell_type=="Macrophages", :]
adata_sub = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.2,n_pcs = 15,output="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Macrophages")
adata = adata_sub.sub_cluster()

sc.pl.umap(adata, color=['Areg','Cd14','Cd300e','Adgre1', 'Cd86', 'Arg1', 'Mki67','Cd209a','Itgax'], ncols=5,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Macrophages/umap_marker_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Macrophages/umap_marker_genes.pdf')
sc.pl.umap(adata, color=["Il1a","Il1b","Il6","Cd80","Cd86","Fcgr2b","Cd163","Mrc1","Cd68"], ncols=5,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Macrophages/umap_M1M2_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Macrophages/umap_M1M2_genes.pdf')


# 对肾实质行细分：分不了
# 'B cells', 'Distal tubule cells', 'Endothelial cells','Epithelial cells', 'Macrophages', 'Neutrophils', 'Proximal tubule cells','Smooth muscle cells', 'T cells'
subcluster=['Distal tubule cells','Endothelial cells','Epithelial cells','Proximal tubule cells','Smooth muscle cells']
adata_sub = adata_raw[adata_raw.obs.cell_type.isin(subcluster), :]

adata_sub = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.3,n_pcs = 15,output="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/kidney cells")
adata = adata_sub.sub_cluster()

sc.pl.umap(adata, color=['Slc34a1', 'Lrp2', 'Nphs1','Nphs2', 'Acta2','Umod','Slc12a1','Slc12a3','Pvalb','Atp6v0d2','Atp6v1g3','Aqp2','Flt1','Ehd3','Eln'], ncols=5,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/kidney cells/umap_marker_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/kidney cells/umap_marker_genes.pdf')


# 在大类增加细胞
anno_dict={
    '22':'Podocyte',
    '11':'Myofibroblast',
    '8':'Loop of henle',
    '14':'Distal tuble cells',
    '12':'Collecting duct_intercalated cells',
    '7':'Collecting duct_principal cells',
}
mask = adata_raw.obs['cluster'].isin(anno_dict.keys())
# 获取新的类别列表
new_categories = list(adata_raw.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata_raw.obs['cell_type'].cat.categories]
# 更新 categories
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.set_categories(new_categories)
adata_raw.obs.loc[mask, 'cell_type'] = adata_raw.obs.loc[mask, 'cluster'].map(anno_dict)
out_dir="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/kidney cells_sub/"
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

adata_raw.write(os.path.join(out_dir, 'kidney cells_sub.h5ad'))


#02.23位置，明天 PTC 细分
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250070/01.Mission/kidney cells_sub/kidney cells_sub.h5ad")
subcluster=['Proximal tubule cells','Epithelial cells']
adata_sub = adata_raw[adata_raw.obs.cell_type.isin(subcluster), :]
adata_sub = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 30,output="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/PTC cells/")
adata = adata_sub.sub_cluster()

adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWSC20250070/01.Mission/PTC cells/all_cells_features.xls"), sep='\t')

sc.pl.umap(adata, color=['Slc34a1','Slc5a2','Slc22a6','Havcr1','Fxyd5','Krt20','Mki67'], ncols=4,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/PTC cells/umap_marker_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/PTC cells/umap_marker_genes.pdf')

anno_dict={
    '0':'Healthy',
    '1':'PTC_Btg2',
    '2':'PTC_Mt1',
    '3':'Healthy',
    '4':'Healthy',
    '5':'Healthy',
    '6':'Healthy',
    '7':'PTC_Miox',
    '8':'PTC_Cdkn1a',
    '9':'Healthy',
    '10':'PTC_Bcam',
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
out_dir="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/PTC cells/"
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

adata.write(os.path.join(out_dir, 'PTC cells.h5ad'))


# 对Endothelial cells进行细分
# 'B cells', 'Distal tubule cells', 'Endothelial cells','Epithelial cells', 'Macrophages', 'Neutrophils', 'Proximal tubule cells','Smooth muscle cells', 'T cells'
subcluster=['Endothelial cells']
adata_sub = adata_raw[adata_raw.obs.cell_type.isin(subcluster), :]
adata_sub = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.3,n_pcs = 15,output="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Endothelial cells")
adata = adata_sub.sub_cluster()

adata.obs.to_csv(os.path.join("/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Endothelial cells/all_cells_features.xls"), sep='\t')

sc.pl.umap(adata, color=['Flt1','Ehd3','Eln'], ncols=5,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Endothelial cells/umap_marker_genes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Endothelial cells/umap_marker_genes.pdf')

anno_dict={
    '0':'EC_1',
    '1':'EC_1',
    '2':'EC_2',
    '3':'EC_1',
    '4':'EC_1',
    '5':'EC_1'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
out_dir="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Endothelial cells/"
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

adata.write(os.path.join(out_dir, 'Endothelial cells.h5ad'))



# 合并，recluster
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250070/01.Mission/kidney cells_sub/kidney cells_sub.h5ad")
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.remove_unused_categories()
out_dir="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/kidney cells_sub/"
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

adata_raw.write(os.path.join(out_dir, 'kidney cells_sub.h5ad'))

adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250070/01.Mission/kidney cells_sub/kidney cells_sub.h5ad")
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].replace({'Epithelial cells': 'Proximal tubule cells'})

bdata = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250070/01.Mission/T cells/subcluster_T cells.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
bdata = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250070/01.Mission/Endothelial cells/result.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')

out_dir="/data/NAS02/sh/Projects/FWSC20250070/01.Mission/mapping_recluster/"
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

adata_raw.write(os.path.join(out_dir, 'mapping_recluster.h5ad'))

out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250070/01.Mission/mapping_recluster/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='n_genes_by_counts', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_n_genes_by_counts.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_n_genes_by_counts.pdf'))