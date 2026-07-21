import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
import pandas as pd

# 对大类重新注释harmony
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20250829/result.h5ad")
adata_raw.obs['platform'] = adata_raw.obs['sample']
adata_raw.obs['platform'] = adata_raw.obs['platform'].astype(str)
adata_raw.obs['platform'] = adata_raw.obs['platform'].replace('TB102-case_1', 'a')
adata_raw.obs['platform'] = adata_raw.obs['platform'].replace('TB113-case_2', 'a')
adata_raw.obs['platform'] = adata_raw.obs['platform'].replace('TB125-control_1', 'a')
adata_raw.obs['platform'] = adata_raw.obs['platform'].replace('TB133-control_2', 'a')
adata_raw.obs['platform'] = adata_raw.obs['platform'].replace('case-1', 'b')
adata_raw.obs['platform'] = adata_raw.obs['platform'].replace('case-2', 'b')
adata_raw.obs['platform'] = adata_raw.obs['platform'].replace('control-1', 'b')
adata_raw.obs['platform'] = adata_raw.obs['platform'].replace('control-2', 'b')
adata_raw.obs['platform'] = adata_raw.obs['platform'].astype('category')
# adata_raw = adata_raw.raw.to_adata()
# 'B cells', 'Distal tubule cells', 'Endothelial cells','Epithelial cells', 'Macrophages', 'Neutrophils', 'Proximal tubule cells','Smooth muscle cells', 'T cells'
# adata_sub = adata_raw[adata_raw.obs.cell_type=="T cells", :]
# 改了recluster 的代码，加了一个去批次的选择
# adata_sub = recluster(adata=adata_raw,species='human',neighbors=30,resolution=0.3,n_pcs = 20,batch_effect = 1,output="/data/NAS02/sh/Projects/FWSC20250829/00.data/recluster/")
# adata = adata_sub.sub_cluster()

out_dir="/data/NAS02/sh/Projects/FWSC20250829/00.data/recluster/"

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

method='seurat_clusters'
adata = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250829/00.data/recluster/anno.h5ad")
sc.pl.umap(adata, color=method,legend_loc ='on data', show=False)
save_fig(plt, out_dir, file_name='sub_cluster_umap')
method='cell_type'
sc.pl.umap(adata, color=method,legend_loc ='on data', show=False)
save_fig(plt, out_dir, file_name='sub_cluster_umap_withlabel')
# sc.tl.tsne(adata)
# sc.pl.tsne(adata, color=method, size=25, show=False)
# save_fig(plt, out_dir, file_name='sub_cluster_tsne')
# sc.pl.tsne(adata, color=method,size=25, legend_loc ='on data', show=False)
# save_fig(plt, out_dir, file_name='sub_cluster_tsne_withlabel')
method='seurat_clusters'
print(adata.obs.columns.tolist())
# adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype(str).astype(object).astype("category")
adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype("category")
sc.tl.rank_genes_groups(adata,
                        groupby=f"{method}",
                        method='wilcoxon',
                        pts=True)
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False,show=False)
save_fig(plt, out_dir, file_name='rank_genes_groups', fig_size={'width':20, 'height':12})

adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype(str).astype("category")
# 生成 dendrogram 数据
sc.tl.dendrogram(adata, groupby='seurat_clusters', key_added='dendrogram_seurat_clusters')

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='seurat_clusters',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, "top10_diff_genes_dotplot.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, "top10_diff_genes_dotplot.pdf"),bbox_inches='tight')
plt.close()

adata.write(os.path.join(out_dir, 'anno_result.h5ad'))

sc.pl.umap(adata, color=['FGFBP2'], ncols=5,cmap='coolwarm')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250829/00.data/recluster/umap_T_diffgenes.png')
plt.savefig(f'/data/NAS02/sh/Projects/FWSC20250829/00.data/recluster/umap_T_diffgenes.pdf')

find_group_marker(adata=adata,species='human',output=out_dir)
adata.obs['cell_type'].cat.categories
anno_dict={
    '3':'T cells',
    '4':'T cells',
    '5':'T cells',
    '2':'T cells',
    '18':'T cells',
    '8':'T cells',
    '14':'Mast cells',
    '17':'Plasmacytoid dendritic cells'
}
mask = adata.obs['seurat_clusters'].isin(anno_dict.keys())
# 获取新的类别列表
new_categories = list(adata.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata.obs['cell_type'].cat.categories]
# 更新 categories
adata.obs['cell_type'] = adata.obs['cell_type'].cat.set_categories(new_categories)
adata.obs.loc[mask, 'cell_type'] = adata.obs.loc[mask, 'seurat_clusters'].map(anno_dict)

adata.obs['cell_type'] = adata.obs['cell_type'].str.replace('T memory cells', 'T cells')

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

samples = adata.obs['sample'].drop_duplicates().tolist()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["sample"]==samples[i]], color="cell_type", legend_loc=legend, title=samples[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'umap_sample_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'umap_sample_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')

adata.write(os.path.join(out_dir, 'anno_result.h5ad'))