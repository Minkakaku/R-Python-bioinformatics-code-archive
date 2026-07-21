import os, sys
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["OMP_NUM_THREADS"] = "16"
import numpy as np
import scanpy as sc
import matplotlib
from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from scanpy_fun_subcluster import find_group_marker
import pandas as pd
from matplotlib import font_manager as fm

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

def plot_celltype_proportion(adata, out_dir,type='sample', fig_size=None):
    keys = adata.obs['cell_type'].cat.categories.to_list()
    values = list(adata.uns['cell_type_colors'])
    cell_type_colors = dict(zip(keys, values))
    ax = pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').plot(
        kind='bar', xlabel='', stacked=True,
        color=[cell_type_colors[cell] for cell in pd.crosstab(adata.obs[type], adata.obs['cell_type'], normalize='index').columns],
        edgecolor='black', width=0.8
    )
    ax.set_xlabel(type.capitalize(), fontsize=14, family='sans-serif')
    ax.set_ylabel("Proportion", fontsize=14, family='sans-serif')
    legend_font = fm.FontProperties(family='sans-serif', size=10)
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., fontsize=10, title="Cell Types", 
        title_fontsize=12, prop=legend_font, frameon=False, handletextpad=1.5, labelspacing=1.0, 
        columnspacing=1.0, markerscale=1, alignment='left', handlelength=2, handleheight=2
    )
    ax.tick_params(axis='x', labelsize=12, labelrotation=90, labelcolor='black')
    ax.tick_params(axis='y', labelsize=12, labelcolor='black')
    ax.grid(False)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    if fig_size is None:
        fig_size = {'width': len(set(adata.obs[type])) * 2, 'height': 6}
    plt.tight_layout()
    file_name=f"celltype_proportion_{type}"
    save_fig(plt, out_dir, file_name=file_name, fig_size=fig_size)

adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2511266943-05/03.Mission/Main_mapping/Main_mapping_filter.h5ad")
adata.obs['cell_type'].cat.categories
adata.obs['sub_cell_type'].cat.categories
bdata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2511266943-05/04.Mission/Main_mapping/Main_mapping_filter.h5ad")
bdata.obs['cell_type'].cat.categories
bdata.obs['sub_cell_type'].cat.categories
# 'B cells', 'Basophils', 'Erythroid cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2511266943-05/05.Mission/"
# with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.umap(bdata, color='cell_type', title='Cell Type', show=False)
#     plt.savefig(os.path.join(out_dir, 'bdata_dimplot_umap_celltype.png'))
# bdata = bdata[bdata.obs['sample'].isin(['iMCD_1_after','iMCD_15_after']), :]

# cdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/proportion/Main_mapping.h5ad")
# cdata.obs['sample'].cat.categories
# cdata = cdata[cdata.obs['group'].isin(['N']), :]

cdata =  adata.concatenate(bdata)
cdata.write(os.path.join(out_dir, 'Main_mapping_combine.h5ad'))



ddata.obs['group'].cat.categories
fdata =  ddata.concatenate(cdata)
fdata.obs['group'].cat.categories

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(fdata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'fdata_dimplot_umap_celltype.png'))

fdata.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
fdata.uns['log1p']={'base':None}
adata_recluster = recluster(adata=fdata,species='human',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
adata.obs['raw_cell_type'] = adata.obs['cell_type']

anno_dict = {
    'Activated B':'B cells',
    'B cells':'B cells',
    'Basophils':'Basophils',
    'CCL classical monocyte':'Myeloid cells',
    'CD4 Tcm':'T cells',
    'CD4 Teff/Temra':'T cells',
    'CD4 Tn':'T cells',
    'CD4 Treg':'T cells',
    'CD8 Teff/Temra':'T cells',
    'CD8 Tem':'T cells',
    'CD8 Tex':'T cells',
    'CD8 Tn':'T cells',
    'CX3CR1 NK':'NK cells',
    'Dividing T':'T cells',
    'Erythroid cells':'Erythroid cells',
    'GZMK NK':'NK cells',
    'IFI classical monocyte':'Myeloid cells',
    'IgA Plasma C1':'Plasma cells',
    'IgA Plasma C2':'Plasma cells',
    'IgG Plasma C1':'Plasma cells',
    'IgG Plasma C2':'Plasma cells',
    'KLRC2 NK':'NK cells',
    'LRMDA classical monocyte':'Myeloid cells',
    'LRMDA non-classical monocyte':'Myeloid cells',
    'Megakaryocyte':'Megakaryocyte',
    'Megakaryocytes':'Megakaryocyte',
    'Memory B':'B cells',
    'Myeloid cells':'Myeloid cells',
    'NK cells':'NK cells',
    'NKT':'NKT cells',
    'NKT cells':'NKT cells',
    'Naive B':'B cells',
    'Neutrophils':'Neutrophils',
    'PPBP classical monocyte':'Myeloid cells',
    'Plasma cells':'Plasma cells',
    'RETN classical monocyte':'Myeloid cells',
    'T cells':'T cells',
    'T/NK':'T/NK',
    'cDC':'Myeloid cells',
    'classical monocyte':'Myeloid cells',
    'non-classical monocyte':'Myeloid cells',
    'pDC':'pDCs',
    'pDCs':'pDCs',
    'γδ T':'T cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['raw_cell_type']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color=['CD3D','CD4','CD8A','NCAM1','GZMA'], title='Cell Type', show=False,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_gene.png'))

sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color='leiden', show=False)
save_fig(plt, out_dir, file_name='sub_cluster_umap')
sc.pl.umap(adata, color='leiden',legend_loc ='on data', show=False)
save_fig(plt, out_dir, file_name='sub_cluster_umap_withlabel')

adata = adata[~adata.obs['leiden'].isin(['6','12','13']), :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
adata.write(os.path.join(out_dir, 'Main_mapping_clean.h5ad'))


# 对B/Plasma cells重新注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/18.Mission/Main_mapping/Main_mapping_clean.h5ad")
# 'B cells', 'Basophils', 'Erythroid cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/18.Mission/B-Plasma/"
adata_sub = adata_raw[adata_raw.obs['leiden'].isin(['7']), :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.8,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'B cells',
    '1':'Plasma cells',
    '2':'B cells',
    '3':'Plasma cells',
    '4':'B cells',
    '5':'B cells',
    '6':'B cells',
    '7':'B cells',
    '8':'B cells',
    '9':'pDCs',
    '10':'B cells',
    '11':'B cells',
    '12':'B cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_map = adata[adata.obs['cell_type'].isin(['Plasma cells']), :]
sc.tl.leiden(adata_map, resolution=0.5)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_map, color='raw_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'Plasma cells',
    '1':'Plasma cells',
    '2':'Plasma cells',
    '3':'B cells',
    '4':'Plasma cells'
}
adata_map.obs['cell_type'] = [anno_dict[clust] for clust in adata_map.obs['leiden']]
adata_map.obs['cell_type'] = adata_map.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_map, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cell_type_mapping = adata_map.obs[['cell_type']].to_dict()['cell_type']
adata.obs['cell_type'] = adata.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].cat.categories.to_list()
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

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'B_Plasma.h5ad'))

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 对T/NK重新注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/18.Mission/Main_mapping/Main_mapping_clean.h5ad")
# 'B cells', 'Basophils', 'Erythroid cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/18.Mission/T-NK/"
adata_sub = adata_raw[adata_raw.obs['leiden'].isin(['0','1','2','5','8']), :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.umap(adata, color=['CD4','CD8A','FOXP3','GZMA','NKG7','CD3D','CD3E','NCAM1'], cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))

anno_dict={
    '0':'T cells',
    '1':'NKT cells',
    '2':'T cells',
    '3':'NK cells',
    '4':'T cells',
    '5':'T cells',
    '6':'T cells',
    '7':'T cells',
    '8':'T cells',
    '9':'T cells',
    '10':'T cells',
    '11':'T cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].cat.categories.to_list()
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

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'T_NK.h5ad'))


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/18.Mission/Main_mapping/Main_mapping_clean.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/18.Mission/mapping_raw/"

adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/18.Mission/T-NK/T_NK.h5ad")
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)

bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/18.Mission/B-Plasma/B_Plasma.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)

adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
adata_raw.obs['cell_type'].cat.categories

# sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
# plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
# plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

mask = adata_raw.obs['leiden'].isin(['3','4','9','10'])          # 找到 Monocytes 行
adata_raw.obs.loc[mask, 'cell_type'] = 'Myeloid cells'
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.remove_unused_categories()
adata_raw = adata_raw[~adata_raw.obs['cell_type'].isin(['Basophils','T/NK']), :]

sc.pl.umap(adata_raw, color='cell_type', show=False)
plt.savefig(os.path.join(out_dir, 'cell_type_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', show=False,legend_loc ='on data')
plt.savefig(os.path.join(out_dir, 'cell_type_umap_withlabel.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type_umap_withlabel.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')

adata_raw.obs['group'] = adata_raw.obs['group'].cat.add_categories('N_after')
mask = adata_raw.obs['group'] == 'iMCD_after'
adata_raw.obs.loc[mask, 'group'] = 'N_after'
adata_raw.obs['group'] = adata_raw.obs['group'].cat.remove_unused_categories()
adata_raw.obs['group'] = (
    pd.Categorical(adata_raw.obs['group'],
                   categories=['HD', 'N', 'N_after', 'P', 'P-after'],
                   ordered=True)
)

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/18.Mission/gene_plot/"
sc.pl.dotplot(adata_raw, var_names=['PDK1','AKT1','MTOR'], groupby=['group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'group_bubble_new.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'group_bubble_new.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata_raw, var_names=['PDK1','AKT1','MTOR'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble_new.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble_new.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_sub=adata_raw[adata_raw.obs['cell_type'].isin(['Plasma cells']), :]
sc.pl.dotplot(adata_sub, var_names=['PDK1','AKT1','MTOR'], groupby=['group'], dendrogram=False, show=False)#, dot_max=0.8,expression_cutoff=0.5
plt.savefig(os.path.join(out_dir, 'group_bubble_new_Plasma.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'group_bubble_new_Plasma.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

# method = "cell_type"
# keys = adata_raw.obs[method].cat.categories.to_list()
# values = list(adata_raw.uns[method+'_colors'])
# color_dict = dict(zip(keys,values))
# df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# # df.columns.names = ["cell_type"]
# df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

# groups = ['HD', 'P', 'P-after']
# if len(groups) > 1:
#     axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
#     fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
#     for i in range(len(groups)):
#         if i == len(groups) - 1:
#             legend = 'right margin'
#         else:
#             legend = None
#         sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)