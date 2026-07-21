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

adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12_mod.Mission/mapping_raw/Main_mapping.h5ad")
# ['B cells', 'Basophils', 'Erythroid cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_A_umap/"
adata = adata[~adata.obs.cell_type.isin(['Erythroid cells','Neutrophils','Basophils']), :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cluster', title='Cell Type', show=False,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color=['CD3D','CD4','CD8A','NCAM1','GZMA'], title='Cell Type', show=False,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))

anno_dict = {
    '0':'T cells',
    '1':'Myeloid cells',
    '2':'T cells',
    '3':'NK cells',
    '4':'T cells',
    '5':'Myeloid cells',
    '6':'NKT cells',
    '8':'B cells',
    '9':'Myeloid cells',
    '10':'T cells',
    '11':'Megakaryocyte',
    '12':'T cells',
    '14':'Myeloid cells',
    '15':'T cells',
    '19':'Plasma cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['cluster']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

adata.obs['cell_type'] = adata.obs['cell_type'].cat.add_categories('pDCs')
adata.obs.loc[adata.obs['sub_cell_type'] == 'pDCS', 'cell_type'] = 'pDCs'

adata.uns['cell_type_colors']=['#1f77b4', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#b5bd61', '#17becf', '#aec7e8']
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.tiff'))

name_map = {
    'HD': 'HD',
    'P': 'P',
    'P-after': 'P_after',
}
adata.obs['group'] = adata.obs['group'].map(name_map).astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='group', title='group', show=False)
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_group.pdf'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_group.tiff'))

desired_order=['NC_1', 'NC_2', 'NC_3', 'iMCD_3', 'iMCD_5', 'iMCD_7', 'iMCD_17', 'iMCD_20', 'iMCD_21', 'iMCD_5_after', 'iMCD_20_after','iMCD_21_after']
adata.obs['sample'] = (
    pd.Categorical(adata.obs['sample'],
                   categories=desired_order,
                   ordered=True)
)
name_map = {
    'NC_1': 'HD_1',
    'NC_2': 'HD_2',
    'NC_3': 'HD_3',
    'iMCD_3': 'P_1',
    'iMCD_5': 'P_2',
    'iMCD_7': 'P_3',
    'iMCD_17': 'P_4',
    'iMCD_20': 'P_5',
    'iMCD_21': 'P_6',
    'iMCD_5_after': 'P_2_after',
    'iMCD_20_after': 'P_5_after',
    'iMCD_21_after': 'P_6_after'
}
adata.obs['sample'] = adata.obs['sample'].map(name_map).astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='sample', show=False)
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_sample.pdf'))
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_sample.tiff'))

adata.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata.obs.to_csv(os.path.join(out_dir, "Main_all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

adata_raw=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_A_umap/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_B_ROE/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='sub_cell_type', show=False)
    plt.savefig(os.path.join(out_dir, 'Main_dimplot_umap_sub_cell_type.png'))

method = "group"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')

adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['Myeloid cells','pDCs']), :]
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_B_ROE/Myeloid/"
adata_sub.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
method = "sub_cell_type"
keys = adata_sub.obs[method].cat.categories.to_list()
values = list(adata_sub.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_C_barplot/Myeloid/"
sc.tl.rank_genes_groups(adata_sub, groupby="sub_cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata_sub,species='human',output=out_dir)

adata_sub = adata_raw[adata_raw.obs['sub_cell_type'].isin(['Activated B','Memory B C1', 'Memory B C2', 'Naive B']), :]
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_B_ROE/Bcells/"
adata_sub.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
method = "sub_cell_type"
keys = adata_sub.obs[method].cat.categories.to_list()
values = list(adata_sub.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_C_barplot/Bcells/"
sc.tl.rank_genes_groups(adata_sub, groupby="sub_cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata_sub,species='human',output=out_dir)

adata_sub = adata_raw[adata_raw.obs['sub_cell_type'].isin(['IgA Plasma', 'IgG Plasma C1','IgG Plasma C2', 'IgM Plasma']), :]
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_B_ROE/Plasma/"
adata_sub.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
method = "sub_cell_type"
keys = adata_sub.obs[method].cat.categories.to_list()
values = list(adata_sub.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')




############################
anno_dict={
    '0':'VCAN Classical monocyte',
    '1':'VCAN Classical monocyte',
    '2':'DUSP1 Classical monocyte',
    '3':'LRMDA Classical monocyte',
		'4':'VCAN Classical monocyte',
    '5':'Non-classical monocyte',
    '6':'LRMDA Classical monocyte',
    '7':'LRMDA Classical monocyte',
    '8':'DUSP1 Classical monocyte',
    '9':'VCAN Classical monocyte',
    '10':'CCL2 Classical monocyte',
    '11':'VCAN Classical monocyte',
    '12':'cDCs',
    '13':'VCAN Classical monocyte',
    '14':'LRMDA Classical monocyte',
    '15':'LRMDA Classical monocyte',
    '16':'pDCS',
    '17':'LRMDA Classical monocyte',
    '18':'Non-classical monocyte',
    '19':'CCL2 Classical monocyte',
    '20':'DUSP1 Classical monocyte'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['HD', 'P', 'P-after']
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
adata.write(os.path.join(out_dir, 'Myeloid.h5ad'))

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 对B/Plasma cells重新注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/mapping_raw/Main_mapping.h5ad")
# ['B cells', 'Basophils', 'Erythroid cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/12_mod.Mission/B-Plasma/"
adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/B-Plasma/B_Plasma.h5ad")

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_map = adata[adata.obs['cell_type'].isin(['IgG Plasma C1','IgG Plasma C2']), :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_map, color='leiden', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '11':'IgG Plasma C2',
    '12':'IgG Plasma C1',
    '15':'IgG Plasma C1',
    '16':'IgG Plasma C2'
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

groups = ['HD', 'P', 'P-after']
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


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/12_mod.Mission/mapping_raw/"

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12_mod.Mission/Myeloid cells/Myeloid.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.pl.umap(adata_raw, color='cell_type', show=False)
plt.savefig(os.path.join(out_dir, 'cell_type_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', show=False,legend_loc ='on data')
plt.savefig(os.path.join(out_dir, 'cell_type_umap_withlabel.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type_umap_withlabel.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

groups = ['HD', 'P', 'P-after']
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)



#### 画il6评分
adata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_A_umap/Main_mapping.h5ad")

adata = adata.raw.to_adata()
adata.raw = adata
adata.obs_names_make_unique()

group = adata.obs['group'].cat.categories.to_list()

def read_gmt(gmt_file_path):
    data = []
    with open(gmt_file_path, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue  # 跳过格式不正确的行
            pathway_name = parts[0]
            pathway_id = parts[1]
            genes = parts[2:]
            data.append([pathway_name, pathway_id] + genes)
    return pd.DataFrame(data)
gene_sets = read_gmt("/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/script/IL6score.gmt")
gene_sets.columns = ['Pathway Name', 'Pathway ID'] + [f'Gene {i+1}' for i in range(len(gene_sets.columns) - 2)]
gene_sets.set_index('Pathway Name', inplace=True)

out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/16.Mission/fig4_E_scoreplot/"
output_data = pd.DataFrame()
# For each gene set term, perform scoring and generate UMAP plots
num_rows = len(gene_sets)
for i in range(num_rows):
    term = gene_sets.index[i]
    genes = gene_sets.iloc[i, 2:].dropna().to_list()  # Get gene list for the gene set
    # Score genes using Scanpy
    sc.tl.score_genes(adata, gene_list=genes, score_name=term)
output_data = adata.obs[['cell_type', 'group']]  # 获取cell_type和group列
for i in range(len(gene_sets)):
    term = gene_sets.index[i]
    output_data[term] = adata.obs[term]     # 将每个基因集得分列加入到结果中        
output_data.to_csv(os.path.join(out_dir, 'sc_score.xls'), sep='\t', index=True)

import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu

plt.figure(figsize=(6, 4))
ax = sns.violinplot(data=output_data,
                    x='cell_type',
                    y='IL6_Signaling_Score',
                    hue='group',
                    palette='Set2',
                    inner='box',
                    scale='width')

plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'IL6_Signaling_Score_vln.png'), facecolor='white', bbox_inches="tight")

plt.figure(figsize=(6, 4))
ax = sns.violinplot(data=output_data,
                    x='cell_type',
                    y='BIOCARTA_IL6_PATHWAY',
                    hue='group',
                    palette='Set2',
                    inner='box',
                    scale='width')

plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'BIOCARTA_IL6_PATHWAY.png'), facecolor='white', bbox_inches="tight")