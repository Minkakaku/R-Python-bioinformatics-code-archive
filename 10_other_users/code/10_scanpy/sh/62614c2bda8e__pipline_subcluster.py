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

adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/05.Mission/mapping_raw/Main_mapping.h5ad")
# adata_raw.raw[:, 'ISG15'].X.toarray()
# adata_raw[:, 'ISG15'].X.toarray()
adata_raw.obs['sample'].cat.categories
bdata = sc.read("/data/outputdata/outputdata-sc/scVision_output/FW2024-225_05/06-Cell_type_annotation/result.h5ad")
# bdata.raw[:, 'ISG15'].X.toarray()
# bdata[:, 'ISG15'].X.toarray()
bdata.obs['sample'].cat.categories
bdata_map = bdata[bdata.obs['sample'].isin(['NC_1','NC_2','NC_3','NC_4']), :]
bdata_map.obs['sample'].cat.categories
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/Main/"

adata =  adata_raw.concatenate(bdata_map)
# sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
# plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
# plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')
adata_sub = adata.copy()
columns_to_drop = [ 'n_cells_by_counts-0', 'mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'n_cells-0', 'means-0', 'dispersions-0', 'dispersions_norm-0', 'n_cells_by_counts-1', 'mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'n_cells-1', 'means-1', 'dispersions-1', 'dispersions_norm-1', 'mean-1', 'std-1']
adata_sub.var.drop(columns=columns_to_drop, inplace=True)
adata_sub.uns['log1p'] = {}
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=30,resolution=0.3,n_pcs = 10,output=out_dir,batch_effect = 3)
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

anno_dict={
    '0':'T/NK',
    '1':'Myeloid cells',
    '2':'T/NK',
    '3':'T/NK',
    '4':'Neutrophils',
    '5':'B cells',
		'6':'Myeloid cells',
    '7':'T/NK',
    '8':'Myeloid cells',
    '9':'Plasma cells',
    '10':'Erythroid cells',
    '11':'Megakaryocytes'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.pl.dotplot(adata, var_names=['MIF'], groupby=['cell_type','sample'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), facecolor='white', bbox_inches="tight")
dp=sc.pl.dotplot(adata, var_names=['MIF'], groupby=['cell_type','group'], dendrogram=False, show=False,return_fig=True)
mean_df = dp.dot_color_df
frac_df = dp.dot_size_df

sc.pl.umap(adata, color=['MIF'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_marker.png'))
plt.savefig(os.path.join(out_dir, 'umap_marker.pdf'))

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

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')


# 对Myeloid cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/Main/Main.h5ad")
# 'B cells', 'Erythroid cells', 'Megakaryocytes', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/Myeloid cells/"
# adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['B cells', 'Plasma cells'])]
adata_sub = adata_raw[adata_raw.obs.cell_type=="Myeloid cells", :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_recluster = recluster(adata=adata_sub,species='human',neighbors=30,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata,groupby="leiden",method='wilcoxon',pts=True)
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'VCAN Classical monocyte',
    '1':'VCAN Classical monocyte',
    '2':'VCAN Classical monocyte',
    '3':'DUSP1 Classical monocyte',
		'4':'lowQ',
    '5':'Non-classical monocyte',
    '6':'ISG15 Classical monocyte',
    '7':'Neutrophils',
    '8':'cDCs',
    '9':'VCAN Classical monocyte',
    '10':'pDCs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub = adata[~adata.obs['leiden'].isin(['4','7']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=30,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

sc.tl.dendrogram(adata, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'VCAN Classical monocyte',
    '1':'VCAN Classical monocyte',
    '2':'DUSP1 Classical monocyte',
    '3':'Non-classical monocyte',
		'4':'ISG15 Classical monocyte',
    '5':'VCAN Classical monocyte',
    '6':'cDCs',
    '7':'VCAN Classical monocyte',
    '8':'pDCs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['NC', 'Negative', 'Positive', 'Rheumatic Diseases']
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

adata_sub.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_sub, out_dir=out_dir,type='group')
plot_celltype_proportion(adata_sub, out_dir=out_dir,type='sample')
adata_sub.write(os.path.join(out_dir, 'Myeloid.h5ad'))

method = "cell_type"
keys = adata_sub.obs[method].cat.categories.to_list()
values = list(adata_sub.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

gene_use = ['MIF']
outdir_tmp = "/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/Myeloid_plot/"
df = sc.get.obs_df(adata, keys=gene_use, use_raw=True)  # cells × 1 (Irf5)
meta = adata.obs
ax1 = sc.pl.dotplot(
 adata,
    var_names=gene_use,  # 显示所有分组
    groupby=['cell_type','group'],  # 按细胞类型分组
    show=True,
    return_fig=True
)
expr_matrix = pd.DataFrame(
    np.nan,  # 初始化为NaN
    index=df.index,
    columns=meta['group'].cat.categories.to_list()
)
for group in expr_matrix.columns:
    group_cells = meta[meta['group'] == group].index
    expr_matrix.loc[group_cells, group] = df.loc[group_cells, gene_use[0]]
new_adata = sc.AnnData(
    X=expr_matrix.values,  # 使用numpy数组
    obs=meta.drop(columns=['group']),  # 保留其他metadata
    var=pd.DataFrame(index=expr_matrix.columns)  # var是分组名
)
a = sc.pl.dotplot(
    adata,
    var_names=gene_use,  # 显示所有分组
    groupby=['cell_type','group'],  # 按细胞类型分组
    show=True,
    return_fig=True
)
b = a.dot_size_df
b['cell_type_group'] = b.index
b[['cell_type', 'group']] = b['cell_type_group'].str.split('_', n=1,expand=True)
b = b.drop(columns='cell_type_group')
b.columns = ['PCT', 'cell_type', 'group']
b = b.pivot(index='cell_type', columns='group', values='PCT')
b.columns.name  = ""
b = b[meta['group'].cat.categories.to_list()]
b = b.loc[meta['cell_type'].cat.categories.to_list()]
sc.pl.dotplot(
    new_adata,
    var_names=new_adata.var_names,  # 显示所有分组
    groupby='cell_type',  # 按细胞类型分组
    dot_size_df = b,
    title=f'{gene_use[0]} expression',
    show=False
)
save_fig(plt, out_dir = outdir_tmp, file_name=f'{gene_use[0]}_dotplot')

sc.tl.rank_genes_groups(adata,groupby="cell_type",method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

# 映射所有细分
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/mapping_raw/"
adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['Myeloid cells']), :]
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))


# 对Plasma cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/Plasma cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Plasma cells", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=50,resolution=0.3,n_pcs = 10,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

adata_test= adata[~adata.obs['sample'].isin(['NC_4']), :]
sc.pl.dotplot(adata, var_names=['MIF'], groupby=['leiden','sample'], dendrogram=False, show=False) #,standard_scale="var"
plt.savefig(os.path.join(out_dir, 'dotplot_sample-Plasma.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_sample-Plasma.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata, var_names=['MIF'], groupby=['group'], dendrogram=False, show=False) #,standard_scale="var"
plt.savefig(os.path.join(out_dir, 'dotplot_group-Plasma.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_group-Plasma.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata_test, var_names=['MIF'], groupby=['group'], dendrogram=False, show=False) #,standard_scale="var"
plt.savefig(os.path.join(out_dir, 'dotplot_group-Plasma_filter.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_group-Plasma_filter.pdf'), facecolor='white', bbox_inches="tight")

sc.pl.umap(adata, color=['IGHA1','IGHG1','IGHM','IGHD','IGHE'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.pdf'))

sc.tl.dendrogram(adata, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'IgA Plasma',
    '1':'IgA Plasma',
    '2':'IgA Plasma',
    '3':'IgA Plasma',
    '4':'IgG Plasma C1',
    '5':'IgA Plasma',
    '6':'IgG Plasma C2',
    '7':'IgA Plasma',
    '8':'IgG Plasma C2',
    '9':'IgG Plasma C1'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['NC', 'Negative', 'Positive', 'Rheumatic Diseases']
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
adata.write(os.path.join(out_dir, 'Plasma.h5ad'))


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/mapping_raw/"

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
method = "group"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')



# 重新提子集做注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/05.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs['sample'].cat.categories
bdata = sc.read("/data/outputdata/outputdata-sc/scVision_output/FW2024-225_05/06-Cell_type_annotation/result.h5ad")
bdata.obs['sample'].cat.categories
bdata_map = bdata[bdata.obs['sample'].isin(['NC_1','NC_2','NC_3']), :]
bdata_map.obs['sample'].cat.categories
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/Main_recluster/"

adata =  adata_raw.concatenate(bdata_map)
adata_sub = adata.copy()
columns_to_drop = [ 'n_cells_by_counts-0', 'mean_counts-0', 'pct_dropout_by_counts-0', 'total_counts-0', 'n_cells-0', 'means-0', 'dispersions-0', 'dispersions_norm-0', 'n_cells_by_counts-1', 'mean_counts-1', 'pct_dropout_by_counts-1', 'total_counts-1', 'n_cells-1', 'means-1', 'dispersions-1', 'dispersions_norm-1', 'mean-1', 'std-1']
adata_sub.var.drop(columns=columns_to_drop, inplace=True)
adata_sub.uns['log1p'] = {}
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=30,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 3)
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

anno_dict={
    '0':'T/NK',
    '1':'Myeloid cells',
    '2':'T/NK',
    '3':'T/NK',
    '4':'T/NK',
    '5':'T/NK',
    '6':'Neutrophils',
    '7':'Myeloid cells',
    '8':'Myeloid cells',
    '9':'B cells',
    '10':'Myeloid cells',
    '11':'Plasma cells',
    '12':'Megakaryocytes',
    '13':'Erythroid cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.pl.dotplot(adata, var_names=['MIF'], groupby=['cell_type','group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), facecolor='white', bbox_inches="tight")
dp=sc.pl.dotplot(adata, var_names=['MIF'], groupby=['cell_type','group'], dendrogram=False, show=False,return_fig=True)
mean_df = dp.dot_color_df
frac_df = dp.dot_size_df

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Main.h5ad'))

groups = ['NC', 'Negative', 'Positive', 'Rheumatic Diseases']
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

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')


# 对Myeloid cells重新注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/Main_recluster/Main.h5ad")
# 'B cells', 'Erythroid cells', 'Megakaryocytes', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/Myeloid_recluster/"
# adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['B cells', 'Plasma cells'])]
adata_sub = adata_raw[adata_raw.obs.cell_type=="Myeloid cells", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=30,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'VCAN Classical monocyte',
    '1':'VCAN Classical monocyte',
    '2':'DUSP1 Classical monocyte',
    '3':'lowQ',
		'4':'ISG15 Classical monocyte',
    '5':'Non-classical monocyte',
    '6':'Neutrophils',
    '7':'cDCs',
    '8':'VCAN Classical monocyte',
    '9':'lowQ',
    '10':'pDCs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub = adata_raw[adata_raw.obs.cell_type.isin(['VCAN Classical monocyte','DUSP1 Classical monocyte','ISG15 Classical monocyte','Non-classical monocyte','cDCs','pDCs']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=30,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

sc.tl.dendrogram(adata, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'VCAN Classical monocyte',
    '1':'VCAN Classical monocyte',
    '2':'DUSP1 Classical monocyte',
    '3':'VCAN Classical monocyte',
		'4':'Non-classical monocyte',
    '5':'ISG15 Classical monocyte',
    '6':'DUSP1 Classical monocyte',
    '7':'cDCs',
    '8':'VCAN Classical monocyte',
    '9':'pDCs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['NC', 'Negative', 'Positive', 'Rheumatic Diseases']
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

gene_use = ['MIF']
outdir_tmp = "/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/Myeloid_recluster_plot/"
df = sc.get.obs_df(adata, keys=gene_use, use_raw=True)  # cells × 1 (Irf5)
meta = adata.obs
ax1 = sc.pl.dotplot(
 adata,
    var_names=gene_use,  # 显示所有分组
    groupby=['cell_type','group'],  # 按细胞类型分组
    show=True,
    return_fig=True
)
expr_matrix = pd.DataFrame(
    np.nan,  # 初始化为NaN
    index=df.index,
    columns=meta['group'].cat.categories.to_list()
)
for group in expr_matrix.columns:
    group_cells = meta[meta['group'] == group].index
    expr_matrix.loc[group_cells, group] = df.loc[group_cells, gene_use[0]]
new_adata = sc.AnnData(
    X=expr_matrix.values,  # 使用numpy数组
    obs=meta.drop(columns=['group']),  # 保留其他metadata
    var=pd.DataFrame(index=expr_matrix.columns)  # var是分组名
)
a = sc.pl.dotplot(
    adata,
    var_names=gene_use,  # 显示所有分组
    groupby=['cell_type','group'],  # 按细胞类型分组
    show=True,
    return_fig=True
)
b = a.dot_size_df
b['cell_type_group'] = b.index
b[['cell_type', 'group']] = b['cell_type_group'].str.split('_', n=1,expand=True)
b = b.drop(columns='cell_type_group')
b.columns = ['PCT', 'cell_type', 'group']
b = b.pivot(index='cell_type', columns='group', values='PCT')
b.columns.name  = ""
b = b[meta['group'].cat.categories.to_list()]
b = b.loc[meta['cell_type'].cat.categories.to_list()]
sc.pl.dotplot(
    new_adata,
    var_names=new_adata.var_names,  # 显示所有分组
    groupby='cell_type',  # 按细胞类型分组
    dot_size_df = b,
    title=f'{gene_use[0]} expression',
    show=False
)
save_fig(plt, out_dir = outdir_tmp, file_name=f'{gene_use[0]}_dotplot')


# 映射所有细分
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/mapping_raw_recluster/"
adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['Myeloid cells']), :]
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))


# 对B/Plasma cells重新注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/mapping_raw_recluster/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/B-Plasma_recluster/"
adata_sub = adata_raw[adata_raw.obs.cell_type.isin(['B cells','Plasma cells']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.8,n_pcs = 10,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

sc.pl.dotplot(adata, var_names=['MIF'], groupby=['group'], dendrogram=False, show=False) #,standard_scale="var"
plt.savefig(os.path.join(out_dir, 'dotplot_group-Plasma.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_group-Plasma.pdf'), facecolor='white', bbox_inches="tight")

sc.pl.umap(adata, color=['IGHA1','IGHG1','IGHM','IGHD','IGHE','CD3D','CD3E','CD3G'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.pdf'))
sc.pl.umap(adata, color=['CD27','CD38','IGHD','FCER2','TCL1A','IL4R'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_NaiveB_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_NaiveB_diffgenes.pdf'))
sc.pl.umap(adata, color=['CD27','CD38','AIM2','TNFRSF13B'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_MemoryB_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_MemoryB_diffgenes.pdf'))
sc.pl.umap(adata, color=['PCNA','CD19','MS4A1','CD27','PTPRC','CD38','CD40','CD83','TNFRSF13B','TNFRSF13C','CD69'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_ActivatedB_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_ActivatedB_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['CD69'], groupby=['leiden'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_ActivatedB.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_ActivatedB.pdf'), facecolor='white', bbox_inches="tight")

sc.tl.dendrogram(adata, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
mask = adata_raw.obs['sub_cell_type'] != '1'
adata_raw.obs['test']=adata_raw.obs['sub_cell_type']
adata_raw.obs.loc[mask, 'test'] = np.nan
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='test', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# 1是低质量
anno_dict={
    '0':'Naive B C1',
    '1':'lowQ',
    '2':'Memory B C1',
    '3':'Naive B C2',
    '4':'Memory B C1',
    '5':'IgA Plasma C1',
    '6':'Naive B C2',
    '7':'Memory B C2',
    '8':'Naive B C2',
    '9':'Naive B C2',
    '10':'Naive B C1',
    '11':'IgA Plasma C2',
    '12':'Activated B',
    '13':'IgG Plasma C2',
    '14':'IgG Plasma C1',
    '15':'IgM Plasma C2',
    '16':'Memory B C2',
    '17':'IgA Plasma C1',
    '18':'Activated B',
    '19':'IgM Plasma C1'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub = adata[~adata.obs.cell_type.isin(['lowQ']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.8,n_pcs = 10,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['IGHA1','IGHG1','IGHM','IGHD','IGHE','CD3D','CD3E','CD3G'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Plasma_diffgenes.pdf'))
sc.pl.umap(adata, color=['CD27','CD38','IGHD','FCER2','TCL1A','IL4R'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_NaiveB_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_NaiveB_diffgenes.pdf'))
sc.pl.umap(adata, color=['CD27','CD38','AIM2','TNFRSF13B'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_MemoryB_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_MemoryB_diffgenes.pdf'))
sc.pl.umap(adata, color=['PCNA','CD19','MS4A1','CD27','PTPRC','CD38','CD40','CD83','TNFRSF13B','TNFRSF13C','CD69'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_ActivatedB_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_ActivatedB_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['CD69'], groupby=['leiden'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_ActivatedB.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_ActivatedB.pdf'), facecolor='white', bbox_inches="tight")

sc.tl.dendrogram(adata, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'Naive B C1',
    '1':'Memory B C1',
    '2':'Naive B C2',
    '3':'Naive B C1',
    '4':'IgA Plasma',
    '5':'Memory B C1',
    '6':'Memory B C2',
    '7':'Naive B C1',
    '8':'Naive B C2',
    '9':'Activated B',
    '10':'IgA Plasma',
    '11':'Activated B',
    '12':'IgG Plasma C2',
    '13':'Naive B C2',
    '14':'IgG Plasma C1',
    '15':'IgM Plasma C2',
    '16':'Memory B C2',
    '17':'IgA Plasma',
    '18':'IgM Plasma C1',
    '19':'Activated B',
    '20':'Naive B C2'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['NC', 'Negative', 'Positive', 'Rheumatic Diseases']
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
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/mapping_raw_recluster/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/09.Mission/mapping_raw_recluster/"

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_sub = adata_raw[adata_raw.obs.sub_cell_type!="lowQ", :]
sc.pl.umap(adata_sub, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_sub.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
method = "group"
keys = adata_sub.obs[method].cat.categories.to_list()
values = list(adata_sub.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')