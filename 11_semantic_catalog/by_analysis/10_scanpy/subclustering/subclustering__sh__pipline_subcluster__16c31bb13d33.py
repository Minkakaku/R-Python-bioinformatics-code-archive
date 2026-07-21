import os, sys
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
import pandas as pd

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

adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20250326/06-Cell_type_annotation/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Ciliated cells', 'Cycling cells', 'Endothelial cells','Fibroblasts', 'Macrophages', 'Monocytes', 'NK cells', 'Neutrophils','Pericytes', 'Plasma cells', 'Pulmonary alveolar cells','Smooth muscle cells', 'T cells'
group_dict={
    'control1':'Control',
    'control2':'Control',
    'BLM1':'Case_control',
    'BLM2':'Case_control',
    'BLM3':'Case_control',
    'BLM4':'Case_control',
    'Treatment1':'Treatment1',
    'Treatment2':'Treatment1',
    'Treatment3':'Treatment2',
    'Treatment4':'Treatment3',
    'Treatment5':'Treatment3',
    'Treatment6':'Treatment3'
}
adata_raw.obs['group'] = [group_dict[samples] for samples in adata_raw.obs['sample']]
adata_raw.obs['group'] = adata_raw.obs['group'].astype("category")

# Fibroblasts细分
out_dir="/data/NAS02/sh/Projects/FWSC20250326/01.Mission/subcluster_Fibroblasts"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Fibroblasts", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=30,resolution=0.3,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()


# 肺泡上皮细胞细分
out_dir="/data/NAS02/sh/Projects/FWSC20250326/01.Mission/subcluster_Pulmonary alveolar cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Pulmonary alveolar cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=30,resolution=0.4,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['Lrp5','Muc5b','Tm4sf1','Ager','Cav1','Pdpn','Lamp3','Sftpa1','Sftpd','Sftpc'], ncols=5,cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'marker_plot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'marker_plot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata, color=['Dcn','Acta2','Col1a1','Col1a2'], ncols=4,cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'marker_plot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'marker_plot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

anno_dict={
    '0':'AT2',
    '1':'AT1',
    '2':'AT1',
    '3':'AT1',
		'4':'AT2',
    '5':'AT1',
    '6':'AT1',
    '7':'AT1',
    '8':'AT1',
		'9':'AT1',
		'10':'AT2'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'AT.h5ad'))

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
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


# 巨噬细胞细分
out_dir="/data/NAS02/sh/Projects/FWSC20250326/01.Mission/subcluster_Macrophages/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Macrophages", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=30,resolution=0.5,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

adata.write(os.path.join(out_dir, 'Macrophages_ori.h5ad'))
adata_sub = adata[adata.obs['leiden'] != '8', :]

anno_dict={
    '0':'M1',
    '1':'M2',
    '2':'M2',
    '3':'M1',
		'4':'M1',
    '5':'DC',
    '6':'M2',
    '7':'M1',
		'9':'DC',
		'10':'M2'
}
adata_sub.obs['cell_type'] = [anno_dict[clust] for clust in adata_sub.obs['leiden']]
adata_sub.obs['cell_type'] = adata_sub.obs['cell_type'].astype('category')

sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

adata_sub.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata_sub.write(os.path.join(out_dir, 'Macrophages.h5ad'))

groups = adata_sub.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_sub[adata_sub.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250326/00.data/result.h5ad")
bdata = sc.read_h5ad("/data/NAS02/sh/Projects/FWSC20250326/01.Mission/subcluster_Macrophages/Macrophages_ori.h5ad")
adata_raw.obs['leiden'] = adata_raw.obs['cell_type']
cell_type_mapping = bdata.obs[['leiden']].to_dict()['leiden']
adata_raw.obs['leiden'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['leiden']), axis=1)
adata_raw.obs['leiden'] = adata_raw.obs['leiden'].astype('category')
adata_raw.obs['leiden'].cat.categories
adata_raw = adata_raw[adata_raw.obs['leiden'] != '8', :]

anno_dict={
    '0':'Macrophages',
    '1':'Macrophages',
    '2':'Macrophages',
    '3':'Macrophages',
		'4':'Macrophages',
    '5':'Macrophages',
    '6':'Macrophages',
    '7':'Macrophages',
		'9':'Macrophages',
		'10':'Macrophages'
}
mask = adata_raw.obs['leiden'].isin(anno_dict.keys())
# 获取新的类别列表
new_categories = list(adata_raw.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata_raw.obs['cell_type'].cat.categories]
# 更新 categories
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.set_categories(new_categories)
adata_raw.obs.loc[mask, 'cell_type'] = adata_raw.obs.loc[mask, 'leiden'].map(anno_dict)
del adata_raw.obs['leiden']
adata_raw.write("/data/NAS02/sh/Projects/FWSC20250326/00.data/result.h5ad")