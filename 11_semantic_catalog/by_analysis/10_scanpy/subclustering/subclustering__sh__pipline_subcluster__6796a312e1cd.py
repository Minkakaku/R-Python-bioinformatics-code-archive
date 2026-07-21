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

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/00.data/result.h5ad")
adata_raw.obs['sample'].cat.categories
# 'B cells', 'Ciliated cells', 'Cycling cells', 'Endothelial cells','Fibroblasts', 'Macrophages', 'Monocytes', 'NK cells', 'Neutrophils','Pericytes', 'Plasma cells', 'Pulmonary alveolar cells','Smooth muscle cells', 'T cells'
samplelist = ['iMCD-PB_1', 'iMCD-PB_2', 'iMCD-PB_3', 'iMCD-PB_5', 'iMCD-PB_6', 'iMCD-PB_7', 'iMCD-PB_8', 'iMCD-PB_11', 'iMCD-PB_15', 'iMCD-PB_16', 'iMCD-PB_17', 'iMCD-PB_18', 'iMCD-PB_20', 'iMCD-PB_21', 'iMCD-PB_24', 'iMCD-PB_25']
adata_sub = adata_raw[adata_raw.obs['sample'].isin(samplelist), :]
adata_sub.obs['sample'].cat.categories
adata_sub.obs['cell_type'].cat.categories
group_dict={
    'iMCD-PB_1':'Negative',
    'iMCD-PB_2':'Negative',
    'iMCD-PB_3':'Positive',
    'iMCD-PB_5':'Positive',
    'iMCD-PB_6':'Rheumatic Diseases',
    'iMCD-PB_7':'Positive',
    'iMCD-PB_8':'Negative',
    'iMCD-PB_11':'Negative',
    'iMCD-PB_15':'Negative',
    'iMCD-PB_16':'Negative',
    'iMCD-PB_17':'Positive',
    'iMCD-PB_18':'Rheumatic Diseases',
    'iMCD-PB_20':'Positive',
    'iMCD-PB_21':'Positive',
    'iMCD-PB_24':'Rheumatic Diseases',
    'iMCD-PB_25':'Rheumatic Diseases'
}
adata_sub.obs['group'] = [group_dict[samples] for samples in adata_sub.obs['sample']]
adata_sub.obs['group'] = adata_sub.obs['group'].astype("category")

out_dir="/data/NAS02/sh/Projects/FWAE20251216/00.data/Main/"
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.3,n_pcs = 20,output=out_dir)
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

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

# sc.pl.umap(adata, color=['Lrp5','Muc5b','Tm4sf1','Ager','Cav1','Pdpn','Lamp3','Sftpa1','Sftpd','Sftpc'], ncols=5,cmap='coolwarm')
# plt.savefig(os.path.join(out_dir, 'marker_plot.png'), facecolor='white', bbox_inches="tight", dpi=300)
# plt.savefig(os.path.join(out_dir, 'marker_plot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
anno_dict={
    '0':'T cells',
    '1':'Myeloid cells',
    '3':'T cells',
    '4':'Neutrophils',
    '8':'T cells',
    '2':'NK cells',
		'9':'NK cells',
    '5':'B cells'
}
mask = adata.obs['leiden'].isin(anno_dict.keys())
# 获取新的类别列表
new_categories = list(adata.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata.obs['cell_type'].cat.categories]
# 更新 categories
adata.obs['cell_type'] = adata.obs['cell_type'].cat.set_categories(new_categories)
adata.obs.loc[mask, 'cell_type'] = adata.obs.loc[mask, 'leiden'].map(anno_dict)

adata.obs['test'] = adata.obs['cell_type']
adata.obs['test'] = adata.obs['leiden'].astype(str) + '_' + adata.obs['test'].astype(str)

anno_dict={
    '6_HSPCs':'T cells',
    '6_Neutrophils':'T cells',
    '6_NK cells':'T cells',
    '6_Erythroid cells':'T cells'
}
mask = adata.obs['test'].isin(anno_dict.keys())
# 获取新的类别列表
new_categories = list(adata.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata.obs['cell_type'].cat.categories]
# 更新 categories
adata.obs['cell_type'] = adata.obs['cell_type'].cat.set_categories(new_categories)
adata.obs.loc[mask, 'cell_type'] = adata.obs.loc[mask, 'test'].map(anno_dict)
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_test.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_test.pdf"),bbox_inches='tight')
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Main.h5ad'))

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
    
adata_raw = sc.read_h5ad("/data/NAS02/sh/Projects/FWAE20251216/00.data/Main/Main.h5ad")
del adata_raw.obs['test']
out_dir="/data/NAS02/sh/Projects/FWAE20251216/00.data/Main/"
adata_raw.write(os.path.join(out_dir, 'Main.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')

out_dir = "/data/NAS02/sh/Projects/FWAE20251216/00.data/script/"
method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')