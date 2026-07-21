import os, sys
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

# 调整大类注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/lung_cancer/anno.h5ad")
# 'Distal tubule cells', 'Endothelial cells', 'Epithelial cells','Fibroblasts', 'Intercalated cells', 'Macrophages', 'Mesangial cells','Podocytes', 'Principal cells', 'Proximal tubule cells','loop of Henle'
out_dir="/mnt/NexData_file/home/sh/Projects/lung_cancer/Main/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='group', title='group', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='patient', title='sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

groups = ['Normal', 'Tumor']
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

adata_raw.uns['log1p'] = {'base': None}
sc.tl.rank_genes_groups(adata_raw,groupby='cell_type', method='wilcoxon',pts=True)
find_group_marker(adata=adata_raw,species='human',output=out_dir)

sc.pl.rank_genes_groups_dotplot(
    adata_raw,
    groupby='cell_type',
    key='rank_genes_groups',
    n_genes=3,
    show=False)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

adata_raw.obs['sample']=adata_raw.obs['patient']
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample')
adata_raw.write(os.path.join(out_dir, 'Main.h5ad'))

method = "cell_type"
method = "group"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')


# 对B注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/lung_cancer/Main/Main.h5ad")
# 'B', 'Endothelial', 'Epithelial', 'Fibroblast', 'Mast', 'Myeloid', 'Plasma', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/lung_cancer/Bcells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="B", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.3,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()

Bmarker=['MS4A1','CXCR4',]#PB
sc.pl.dotplot(adata, var_names=Bmarker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Bmarker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Bmarker.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=Bmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Bmarker.png'))
plt.savefig(os.path.join(out_dir, 'umap_Bmarker.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=3,
    show=False)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'MemoryB',
    '1':'ISGB',
    '2':'MemoryB'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups =  adata.obs['group'].drop_duplicates().tolist()
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
adata.write(os.path.join(out_dir, 'B.h5ad'))

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

# 对Endothelial注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/lung_cancer/Main/Main.h5ad")
# 'B', 'Endothelial', 'Epithelial', 'Fibroblast', 'Mast', 'Myeloid', 'Plasma', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/lung_cancer/Endothelial/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Endothelial", :]
adata_sub.raw=adata_sub
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

Endothelialmarker=['LYVE1','PROX1','CCL21','FLT4',   'CLDN5','ACKR1',   'FBLN5','GJA5','DKK2',   'PTGS1','RGS5','HDAC9',   'CA4','EPAS1','NTRK2']
sc.pl.dotplot(adata, var_names=Endothelialmarker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Endothelialmarker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Endothelialmarker.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=Endothelialmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Endothelialmarker.png'))
plt.savefig(os.path.join(out_dir, 'umap_Endothelialmarker.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=3,
    show=False)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'vein ECs',
    '1':'artery ECs',
    '2':'vein ECs',
    '3':'capillary ECs',
    '4':'artery ECs',
    '5':'lymphatic ECs',
    '6':'capillary ECs',
    '7':'artery ECs',
    '8':'vein ECs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups =  adata.obs['group'].drop_duplicates().tolist()
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
adata.write(os.path.join(out_dir, 'ECs.h5ad'))

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

# 对Epithelial注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/lung_cancer/Main/Main.h5ad")
# 'B', 'Endothelial', 'Epithelial', 'Fibroblast', 'Mast', 'Myeloid', 'Plasma', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/lung_cancer/Epithelial/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Epithelial", :]
adata_sub.raw=adata_sub
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.2,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

Epithelialmarker=['AGER','CAV1',   'SFTPC','ABCA3',   'SCGB1A1','SCGB3A1','CP'   ,'KRT15','TP63','KRT5',   'MUC5B','MUC2',   'FOXJ1','CAPS','PIFO']
sc.pl.dotplot(adata, var_names=Epithelialmarker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Epithelialmarker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Epithelialmarker.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=Epithelialmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Epithelialmarker.png'))
plt.savefig(os.path.join(out_dir, 'umap_Epithelialmarker.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=3,
    show=False)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'Goblet cells',
    '1':'Proliferating',
    '2':'Club cells',
    '3':'AT2',
    '4':'AT1',
    '5':'Basal cells',
    '6':'Proliferating',
    '7':'Ciliated cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups =  adata.obs['group'].drop_duplicates().tolist()
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
adata.write(os.path.join(out_dir, 'Epithelial.h5ad'))

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

# 对Fibroblast注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/lung_cancer/Main/Main.h5ad")
# 'B', 'Endothelial', 'Epithelial', 'Fibroblast', 'Mast', 'Myeloid', 'Plasma', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/lung_cancer/Fibroblast/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Fibroblast", :]
adata_sub.raw=adata_sub
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir)
adata = adata_recluster.sub_cluster()

Fibroblastmarker=['COL10A1','COL4A1','PLA2G2A','MMP3','VEGFD','CCL2']
sc.pl.dotplot(adata, var_names=Fibroblastmarker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Fibroblastmarker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Fibroblastmarker.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=Fibroblastmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Fibroblastmarker.png'))
plt.savefig(os.path.join(out_dir, 'umap_Fibroblastmarker.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=3,
    show=False)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'Fibroblast C1',
    '1':'Fibroblast C2',
    '2':'Fibroblast C3',
    '3':'Fibroblast C4',
    '4':'Fibroblast C5',
    '5':'Fibroblast C2',
    '6':'Fibroblast C3',
    '7':'Fibroblast C4'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups =  adata.obs['group'].drop_duplicates().tolist()
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
adata.write(os.path.join(out_dir, 'Fibroblast.h5ad'))

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 对Myeloid注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/lung_cancer/Main/Main.h5ad")
# 'B', 'Endothelial', 'Epithelial', 'Fibroblast', 'Mast', 'Myeloid', 'Plasma', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/lung_cancer/Myeloid/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Myeloid", :]
adata_sub.raw=adata_sub
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=1,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

Myeloidmarker=['S100A12','CD207','FCGR3A','CYTIP','CLEC9A','XCR1','FCER1A','CD163','IFITM3','RGCC']
sc.pl.dotplot(adata, var_names=Myeloidmarker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Myeloidmarker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Myeloidmarker.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=Myeloidmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Myeloidmarker.png'))
plt.savefig(os.path.join(out_dir, 'umap_Myeloidmarker.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=3,
    show=False)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'Lung-associated Macrophages',
    '1':'Tumour-associated Macrophages',
    '2':'Tumour-associated Macrophages',
    '3':'Lung-associated Macrophages',
    '4':'Lung-associated Macrophages',
    '5':'Langerhans cells',
    '6':'Cross-presenting DCs',
    '7':'Tumour-associated Macrophages',
    '8':'Lung-associated Macrophages',
    '9':'Monocyte-derived DCs',
    '10':'Lung-associated Macrophages',
    '11':'Tumour-associated Macrophages',
    '12':'Langerhans cells',
    '13':'Tumour-associated Macrophages',
    '14':'Tumour-associated Macrophages',
    '15':'Granulocytes',
    '16':'Tumour-associated Macrophages',
    '17':'Tumour-associated Macrophages'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups =  adata.obs['group'].drop_duplicates().tolist()
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
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

# 对Plasma注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/lung_cancer/Main/Main.h5ad")
# 'B', 'Endothelial', 'Epithelial', 'Fibroblast', 'Mast', 'Myeloid', 'Plasma', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/lung_cancer/Plasma/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Plasma", :]
adata_sub.raw=adata_sub
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

Plasmamarker=['IGHA1','IGHG1','IGHM','IGHD','IGHE']
sc.pl.dotplot(adata, var_names=Plasmamarker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Plasmamarker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Plasmamarker.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=Plasmamarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Plasmamarker.png'))
plt.savefig(os.path.join(out_dir, 'umap_Plasmamarker.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=3,
    show=False)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'IgG Plasma',
    '1':'IgA Plasma',
    '2':'IgA Plasma',
    '3':'IgA Plasma',
    '4':'IgG Plasma',
    '5':'IgG Plasma',
    '6':'IgA Plasma',
    '7':'IgG Plasma',
    '8':'IgG Plasma',
    '9':'IgG Plasma',
    '10':'IgG Plasma'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups =  adata.obs['group'].drop_duplicates().tolist()
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

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

# 对T/NK注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/lung_cancer/Main/Main.h5ad")
# 'B', 'Endothelial', 'Epithelial', 'Fibroblast', 'Mast', 'Myeloid', 'Plasma', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/lung_cancer/T_NK/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="T/NK", :]
adata_sub.raw=adata_sub
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=1,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

T_NKmarker=['FOXP3','FGFBP2','CD8A','CD4','LAG3','MKI67']
sc.pl.dotplot(adata, var_names=T_NKmarker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_T_NKmarker.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_T_NKmarker.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=T_NKmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_T_NKmarker.png'))
plt.savefig(os.path.join(out_dir, 'umap_T_NKmarker.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=3,
    show=False)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

anno_dict={
    '0':'Exhausted Tcells',
    '1':'CD4T',
    '2':'CD8T',
    '3':'CD8T',
    '4':'CD4T',
    '5':'NK cells',
    '6':'CD4T',
    '7':'Tregs',
    '8':'CD8T',
    '9':'CD4T',
    '10':'CD8T',
    '11':'CD8T',
    '12':'Tregs',
    '13':'Tregs',
    '14':'CD8T'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups =  adata.obs['group'].drop_duplicates().tolist()
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

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/lung_cancer/Main/Main.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/lung_cancer/mapping_raw/"

adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/lung_cancer/Bcells/B.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/lung_cancer/Endothelial/ECs.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/lung_cancer/Epithelial/Epithelial.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/lung_cancer/Fibroblast/Fibroblast.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/lung_cancer/Myeloid/Myeloid.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/lung_cancer/Plasma/Plasma.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')