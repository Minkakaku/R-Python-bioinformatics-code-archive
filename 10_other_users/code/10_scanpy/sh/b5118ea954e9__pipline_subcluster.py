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

# 调整Macrophage注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/mapping_raw/Main_mapping.h5ad")
# 'Epithelial cells', 'T', 'DC', 'B', 'Neutrophil', 'Fibroblast', 'Macrophage', 'Endothelial cells', 'Plasma', 'Schwann cells', 'Myoepithelial cells', 'Acinar cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/02.Mission/Macrophage/"
adata_sub = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20252753/01.Mission/Macrophage/Macrophage.h5ad")
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
adata_sub =  adata_sub[~adata_sub.obs.leiden.isin(['4','6']), :]
adata_sub.write(os.path.join(out_dir, 'Macrophage.h5ad'))

# groups = adata_sub.obs['group'].drop_duplicates().tolist()
groups = ['Baseline', 'Day1', 'Day3', 'Day5', 'Day7', 'Day14', 'Day21']
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
plot_celltype_proportion(adata_sub, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata_sub, out_dir=out_dir,type='sample',fig_size={'w':20, 'h':5})

# 重新降维聚类
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/02.Mission/mapping_raw/"

cell_type_mapping = adata_sub.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['M1','Doublet']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
sc.tl.umap(adata_raw)
sc.tl.leiden(adata_raw, resolution=0.5)
# plot cluster
sc.pl.umap(adata_raw, color='leiden', show=False)
plt.savefig(os.path.join(out_dir, 'sub_cluster_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'sub_cluster_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.write(os.path.join(out_dir, 'Main_recluster.h5ad'))
sc.pl.umap(adata_raw, color='cell_type', show=False)
plt.savefig(os.path.join(out_dir, 'cell_type_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='sub_cell_type', show=False)
plt.savefig(os.path.join(out_dir, 'sub_cell_type_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'sub_cell_type_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', show=False,legend_loc ='on data')
plt.savefig(os.path.join(out_dir, 'cell_type_umap_withlabel.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type_umap_withlabel.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

cell_type_mapping = adata_sub.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.write(os.path.join(out_dir, 'Main_recluster.h5ad'))


# 对B注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20252753/02.Mission/mapping_raw/Main_recluster.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/02.Mission/Bcells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="B", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=1,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

Bmarker=["Tcl1","Fcer2a","Il4ra" , ####Naiveb
                "Ifit3","Ifi44l","Stat1","Isg15",###Ifn
                "Hspa1a","Dnajb1",###Stress
                "Slc30a1",##Mt1x
                "Egr1","Dusp2",####Acb1
                "Nr4a2","Cd69","Cd83",####Acb2
                "Ccr7","Pim3","Nfkbid",######Acb3
                "S100a10","Crip1","S100a4","Itgb1","Cd27","Cr2","Aim2","Gpr183",#Bm
                "Dusp4","Fcrl5","Zeb2","Itgax","Fgr","Cd274",#Atm
                "Nme1","Psme2","Eno1","Fabp5",###Pregc
                "Cxcr4","Rack1","Atp5l","Sugct",###Dzgc And Gc
                "Lmo2","Gmds","Prpsap2","Marcksl1",###Lzgc
                "Stmn1","Tubb5","Hmgb2","Tuba1b","Mki67"]#PB
sc.pl.dotplot(adata, var_names=Bmarker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_B.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_B.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=Bmarker, groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_B_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_B_scale.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=Bmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_B_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_B_diffgenes.pdf'))

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

anno_dict={
    '0':'Naiveb',
    '1':'Pregc',
    '2':'Acb2',
    '3':'Acb3',
    '4':'Pregc',
    '5':'Bm',
    '6':'Pregc',
    '7':'Stress',
    '8':'Acb3',
    '9':'Acb1',
    '10':'Pregc'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_map = adata[adata.obs['leiden'].isin(['5']), :]
sc.tl.leiden(adata_map, resolution=0.3)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_map, color='leiden', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'Bm',
    '1':'PB'
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


# 对Neutrophil注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20252753/02.Mission/mapping_raw/Main_recluster.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/02.Mission/Neutrophil/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Neutrophil", :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

Neutrophilmarker=['Fut4','Itgam','Itga4','Cd34','Cd38','Cd101','Fcgr4','Tfrc','Kit','Cxcr4','Cxcr2']
sc.pl.dotplot(adata, var_names=Neutrophilmarker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Neutrophil.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Neutrophil.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=Neutrophilmarker, groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_Neutrophil_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Neutrophil_scale.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=Neutrophilmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Neutrophil_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Neutrophil_diffgenes.pdf'))

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
    '0':'Mature neutrophil',
    '1':'Immature neutrophil',
    '2':'Por-neutrophil',
    '3':'Pre-neutrophil',
    '4':'Mature neutrophil'
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
adata.write(os.path.join(out_dir, 'Neutrophil.h5ad'))


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20252753/02.Mission/mapping_raw/Main_recluster.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20252753/02.Mission/mapping_raw/"

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20252753/02.Mission/Bcells/B.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')