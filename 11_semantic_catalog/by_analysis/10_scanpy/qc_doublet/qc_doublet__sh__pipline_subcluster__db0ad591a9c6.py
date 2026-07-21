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

# 对Epithelial cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252753/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Epithelial cells', 'T', 'DC', 'B', 'Neutrophil', 'Fibroblast', 'Macrophage', 'Endothelial cells', 'Plasma', 'Schwann cells', 'Myoepithelial cells', 'Acinar cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Epithelial cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Epithelial cells", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 1)
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

sc.pl.umap(adata, color=['Lgr5','Hmgcs2','Cdx2','Mki67','Top2a','Cenpf','Reg4','Mptx1','Kit','Dclk1','Col1a1','Col1a2','Alpi','Wfdc2','Fgfbp1','Chga','Chgb','Neurod1','Spdef','Smim14','Cgref1'], cmap='coolwarm', ncols=3)
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['Lgr5','Hmgcs2','Cdx2','Mki67','Top2a','Cenpf','Reg4','Mptx1','Kit','Dclk1','Col1a1','Col1a2','Alpi','Wfdc2','Fgfbp1','Chga','Chgb','Neurod1','Spdef','Smim14','Cgref1'], groupby='leiden', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Colonocytes',
    '1':'Colonocytes',
    '2':'DCS cells',
    '3':'TA cells',
    '4':'Stem cells',
    '5':'Tuft cells',
    '6':'TA cells',
    '7':'Tuft cells',
    '8':'Goblet cells',
    '9':'Erythrocyte'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Epithelial cells.h5ad'))

groups = ['Baseline', 'Day1', 'Day3', 'Day5', 'Day7', 'Day14', 'Day21']
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

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':22, 'h':5})

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 对T注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252753/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Tcells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="T", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.8,n_pcs = 10,output=out_dir,batch_effect = 1)
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

sc.pl.umap(adata, color=['Cd4','Cd8a','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['Cd4','Cd8a','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'))

adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'Treg',
    '1':'Th17',
    '2':'lowQ',
    '3':'NK',
    '4':'CD4 naive cells',
    '5':'CD4 naive cells',
    '6':'Treg',
    '7':'lymphoid progenitor',
    '8':'lowQ'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub = adata[~adata.obs.cell_type.isin(['lowQ']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.8,n_pcs = 10,output=out_dir,batch_effect = 1)
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

sc.pl.umap(adata, color=['Cd4','Cd8a','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra','Il17a'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['Cd4','Cd8a','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra','Il17a'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'))

anno_dict={
    '0':'Treg',
    '1':'CD8T',
    '2':'CD4Tn',
    '3':'DNT',
    '4':'CD4Tn',
    '5':'Treg',
    '6':'Th17',
    '7':'Th17',
    '8':'NK',
    '9':'Proliferative Tcell'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['Baseline', 'Day1', 'Day3', 'Day5', 'Day7', 'Day14', 'Day21']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'T.h5ad'))

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 对B注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252753/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Bcells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="B", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'Naive B',
    '1':'Naive B',
    '2':'Memory B',
    '3':'lowQ',
    '4':'GCB'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub = adata[~adata.obs.cell_type.isin(['lowQ']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

sc.tl.dendrogram(adata,groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['Cd27','Cd38','Fcer2a','Tcl1','Il4ra'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_NaiveB_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_NaiveB_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['Cd27','Cd38','Fcer2a','Tcl1','Il4ra'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_NaiveB.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_NaiveB.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=['Cd27','Cd38','Aim2','Tnfrsf13b'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_MemoryB_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_MemoryB_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['Cd27','Cd38','Aim2','Tnfrsf13b'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_MemoryB.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_MemoryB.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=['Pcna','Cd19','Ms4a1','Cd27','Ptprc','Cd38','Cd40','Cd83','Tnfrsf13b','Tnfrsf13c','Cd69'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_ActivatedB_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_ActivatedB_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['Pcna','Cd19','Ms4a1','Cd27','Ptprc','Cd38','Cd40','Cd83','Tnfrsf13b','Tnfrsf13c','Cd69'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_ActivatedB.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_ActivatedB.pdf'),bbox_inches='tight')

anno_dict={
    '0':'MemoryB',
    '1':'NaiveB',
    '2':'NaiveB',
    '3':'ActivatedB'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['Baseline', 'Day1', 'Day3', 'Day5', 'Day7', 'Day14', 'Day21']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'B.h5ad'))


# 对Macrophage注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252753/result.h5ad")
# 'Epithelial cells', 'T', 'DC', 'B', 'Neutrophil', 'Fibroblast', 'Macrophage', 'Endothelial cells', 'Plasma', 'Schwann cells', 'Myoepithelial cells', 'Acinar cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Macrophage/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Macrophage", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 8,output=out_dir,batch_effect = 1)
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

sc.pl.umap(adata, color=['Clec4a3','Il3ra','Itgam','Itgax','Lamp3','Lyz1','Ccl22','Clec9a','Fscn1','Itgae','Tnip3','Wdfy4','Bst2','Ccr9','Cd14','Fcgr4'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['Clec4a3','Il3ra','Itgam','Itgax','Lamp3','Lyz1','Ccl22','Clec9a','Fscn1','Itgae','Tnip3','Wdfy4','Bst2','Ccr9','Cd14','Fcgr4'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Clec9a','Xcr1','Ckb','Naga','Ccl22','Ccr7','Fscn1','Cacnb3','Cd209a','Clec10a','Mgl2','Clec12a','Klrk1','Esam','Tbx21','Dtx1'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_DC.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_DC.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Cd80','Cd86','Fcgr1','Fcgr2b','Cd163','Mrc1','Cd68'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_M1_2.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_M1_2.pdf'),bbox_inches='tight')

anno_dict={
    '0':'M2',
    '1':'M1',
    '2':'Monocyte',
    '3':'M1',
		'4':'Doublet',
    '5':'M2',
    '6':'M1',
    '7':'Monocyte'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

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

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'Macrophage.h5ad'))


# 对Fibroblast注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252753/result.h5ad")
# 'Epithelial cells', 'T', 'DC', 'B', 'Neutrophil', 'Fibroblast', 'Macrophage', 'Endothelial cells', 'Plasma', 'Schwann cells', 'Myoepithelial cells', 'Acinar cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Fibroblast/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Fibroblast", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.6,n_pcs = 20,output=out_dir,batch_effect = 2)
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

sc.pl.umap(adata, color=['Mmp11', 'Cdh11','Postn','Col1a1','Cfd','C3', 'Cd34','Cxcl12', 'Cd248','Acta2', 'Mcam','Rgs5','Pdpn', 'Mme', 'Eno1', 'Ndrg1', 'Car9','Hsph1','Ido1','Cd74', 'Ccl19','Mki67'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_CAF.png'))
plt.savefig(os.path.join(out_dir, 'umap_CAF.pdf'))
sc.pl.dotplot(adata, var_names=['Mmp11', 'Cdh11','Postn','Col1a1','Cfd','C3', 'Cd34','Cxcl12', 'Cd248','Acta2', 'Mcam','Rgs5','Pdpn', 'Mme', 'Eno1', 'Ndrg1', 'Car9','Hsph1','Ido1','Cd74', 'Ccl19','Mki67'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_CAF.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_CAF.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Inflammatory CAFs',
    '1':'Inflammatory CAFs',
    '2':'Antigen-presenting CAFs',
    '3':'Tumour-like CAFs',
    '4':'Matrix CAFs',
    '5':'Inflammatory CAFs',
    '6':'Matrix CAFs',
    '7':'Vascular CAFs',
    '8':'Vascular CAFs',
    '9':'Antigen-presenting CAFs',
    '10':'Antigen-presenting CAFs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

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

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'Fibroblast.h5ad'))


# 对Endothelial cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252753/result.h5ad")
# 'Epithelial cells', 'T', 'DC', 'B', 'Neutrophil', 'Fibroblast', 'Macrophage', 'Endothelial cells', 'Plasma', 'Schwann cells', 'Myoepithelial cells', 'Acinar cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Endothelial cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Endothelial cells", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
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

sc.pl.umap(adata, color=['Sema3c','Kit','Gpihbp1','Ccl21a', 'Prox1','Hey1', 'Igfbp3','Cd36', 'Car4','Ackr1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['Sema3c','Kit','Gpihbp1','Ccl21a', 'Prox1','Hey1', 'Igfbp3','Cd36', 'Car4','Ackr1'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'lymphatic ECs',
    '1':'artery ECs',
    '2':'artery ECs',
    '3':'vein ECs',
    '4':'artery ECs',
    '5':'capillary ECs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['Baseline', 'Day1', 'Day3', 'Day5', 'Day7', 'Day14', 'Day21']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'Endothelial.h5ad'))


# 对DC注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252753/result.h5ad")
# 'Epithelial cells', 'T', 'DC', 'B', 'Neutrophil', 'Fibroblast', 'Macrophage', 'Endothelial cells', 'Plasma', 'Schwann cells', 'Myoepithelial cells', 'Acinar cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/DC/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="DC", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.2,n_pcs = 10,output=out_dir)
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
    '0':'cDC2-Cd209a',
    '1':'cDC1-Clec9a',
    '2':'cDC1-Ccr7',
    '3':'pDC',
    '4':'cDC1-Ccr7',
    '5':'cDC2-Cd209a'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['Baseline', 'Day1', 'Day3', 'Day5', 'Day7', 'Day14', 'Day21']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'DC.h5ad'))

# 对Plasma注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252753/result.h5ad")
# 'Epithelial cells', 'T', 'DC', 'B', 'Neutrophil', 'Fibroblast', 'Macrophage', 'Endothelial cells', 'Plasma', 'Schwann cells', 'Myoepithelial cells', 'Acinar cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Plasma/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Plasma", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.1,n_pcs = 10,output=out_dir)
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
    '0':'Plasma C1',
    '1':'Plasma C2',
    '2':'Plasma C3'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['Baseline', 'Day1', 'Day3', 'Day5', 'Day7', 'Day14', 'Day21']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})
adata.write(os.path.join(out_dir, 'Plasma.h5ad'))


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20252753/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/mapping_raw/"
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Plasma/Plasma.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Tcells/T.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Macrophage/Macrophage.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Epithelial cells/Epithelial cells.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Endothelial cells/Endothelial.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/DC/DC.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Bcells/B.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['B','T']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')