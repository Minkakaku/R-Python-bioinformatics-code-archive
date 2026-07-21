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


adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/Figure1D/"
sc.tl.rank_genes_groups(adata_raw,groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata_raw,species='mouse',output=out_dir)

sc.pl.dotplot(adata_raw, var_names=['Acsm2', 'Fbp1','Miox','C1qb', 'Cd74', 'Fcer1g','Defb1', 'Klk1','Umod','Stmn1','Top2a','Mki67','Adgrf5','Emcn','Adgrf5','Cd3d','Cd3e','Cd3g','Clec4d', 'S100a8', 'S100a9','Cd74', 'Cd83', 'Lsp1','Col4a1','Col4a2','Dcn','Hba-a1','Hba-a2','Hbb-bs','Aqp3','Hsd11b2','Scnn1a','Cd79a','Cd79b','Ms4a1','Rgs5','Myh11','Mylk'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), facecolor='white', bbox_inches="tight")

sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=3, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'heatmap_top3genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top3genes.pdf'),bbox_inches='tight')

sc.tl.rank_genes_groups(adata_raw,groupby='group', method='wilcoxon', pts=True)
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/Figure1E/"
find_group_marker(adata=adata_raw,species='mouse',output=out_dir)


# 对B注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/Bcells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="B cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

Bmarker=["Fcer2a","Il4ra" , ####Naiveb
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
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=Bmarker, groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_genes_scale.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_scale.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=Bmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))

# sc.tl.dendrogram(adata,groupby='leiden')
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
    '0':'Activate B',
    '1':'Naive B',
    '2':'Memory B',
    '3':'Memory B',
    '4':'Activate B',
    '5':'Naive B',
    '6':'Activate B',
    '7':'pDC',
    '8':'Memory B',
    '9':'Activate B'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['H0', 'H2', 'H8', 'H24', 'H72']
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

adata.obs['group'] = (
    pd.Categorical(adata.obs['group'],
                   categories=groups,
                   ordered=True)
)
adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Bcells.h5ad'))

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

# 对DC注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20253251/04.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/07.Mission/DC/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Dendritic cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 1)
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

sc.pl.dotplot(adata, var_names=['Cd209a','Clec9a','Ccr7'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=['Cd209a','Clec9a','Ccr7'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))

anno_dict={
    '0':'cDC2-Cd209a',
    '1':'cDC2-Cd209a',
    '2':'cDC1-Ccr7',
    '3':'cDC1-Ccr7',
    '4':'cDC1-Ccr7',
    '5':'cDC1-Ccr7',
    '6':'cDC1-Clec9a'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['H0', 'H2', 'H8', 'H24', 'H72']
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

adata.obs['group'] = (
    pd.Categorical(adata.obs['group'],
                   categories=groups,
                   ordered=True)
)
adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'DC.h5ad'))

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)