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

# 重新换颜色出图
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2509225162/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Main/"
ct_colors = pd.Series(adata_raw.uns['cell_type_colors'], index=adata_raw.obs['cell_type'].cat.categories)  # 对应类别顺序
ct_colors['Endothelial cells']     = '#57D1C9' 
ct_colors['Basophils'] = '#7F95D1' 
adata_raw.uns['cell_type_colors'] = ct_colors.values
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,legend_loc='on data', wspace=0.1, legend_fontsize=8, legend_fontweight='normal', legend_fontoutline=2)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_ondata.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_ondata.pdf'))

groups = adata_raw.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        sc.pl.umap(adata_raw[adata_raw.obs["group"]==groups[i]], color="cell_type",legend_loc='on data', wspace=0.1, legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2,title=groups[i], ax=axs[i], show=False)
    save_fig(plt, out_dir, file_name='dimplot_umap_group_celltype', fig_size={'width':len(groups)*5, 'height':5})

adata_raw.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',fig_size={'w':7, 'h':5})
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':14, 'h':5})
adata_raw.write(os.path.join(out_dir, 'Main.h5ad'))


# 对Tumor-like cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Main/Main.h5ad")
# 'Tumor-like cells', 'T/NK', 'Myeloid cells', 'Fibroblasts', 'Neutrophils', 'Endothelial cells', 'Basophils'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Tumor-like cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Tumor-like cells", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir)
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
    '0':'Dcn Tumor-like cells',
    '1':'Hist1h2ap Tumor-like cells',
    '2':'Lyz2 Tumor-like cells',
    '3':'Dcn Tumor-like cells',
    '4':'Cdc20 Tumor-like cells',
    '5':'Dcn Tumor-like cells',
    '6':'Zfpm2 Tumor-like cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':14, 'h':5})
adata.write(os.path.join(out_dir, 'Tumor.h5ad'))


# 对T/NK注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Main/Main.h5ad")
# 'Tumor-like cells', 'T/NK', 'Myeloid cells', 'Fibroblasts', 'Neutrophils', 'Endothelial cells', 'Basophils'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/T_NK/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="T/NK", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.8,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['Cd4','Cd8a','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra','Il17a','Trdc','Sell','Ccr7','Cd69','Rgs1','Rorc','Ifng','Tox','Ncr1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_gene.png'))
plt.savefig(os.path.join(out_dir, 'umap_gene.pdf'))
sc.pl.dotplot(adata, var_names=['Cd4','Cd8a','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra','Il17a','Trdc','Sell','Ccr7','Cd69','Rgs1','Rorc','Ifng','Tox','Ncr1'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_gene.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_gene.pdf'),bbox_inches='tight')

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
    '0':'CD8T',
    '1':'CD4T',
    '2':'CD8T',
    '3':'CD4T',
    '4':'NK',
    '5':'Treg',
    '6':'ProliferativeT',
    '7':'CD4T',
    '8':'Th17',
    '9':'CD4T',
    '10':'CD4T',
    '11':'CD8T',
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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':14, 'h':5})
adata.write(os.path.join(out_dir, 'T_NK.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)


# 对Myeloid cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Main/Main.h5ad")
# 'Tumor-like cells', 'T/NK', 'Myeloid cells', 'Fibroblasts', 'Neutrophils', 'Endothelial cells', 'Basophils'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Myeloid/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Myeloid cells", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.3,n_pcs = 25,output=out_dir)
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

sc.pl.umap(adata, color=['Cd80','Cd86','Fcgr1','Fcgr2b','Cd163','Mrc1','Cd68','Lyz2','Clec9a','Cd209a','Tnip3','Cd14','Fcgr4'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata, var_names=['Cd80','Cd86','Fcgr1','Fcgr2b','Cd163','Mrc1','Cd68','Lyz2','Clec9a','Cd209a','Tnip3','Cd14','Fcgr4'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Cd74 Macrophage',
    '1':'Pf4 Macrophage',
    '2':'Spp1 Macrophage',
    '3':'Cd74 Macrophage',
    '4':'Proliferative Macrophage',
    '5':'cDC',
    '6':'Pf4 Macrophage',
    '7':'Cd74 Macrophage',
    '8':'pDC'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':14, 'h':5})
adata.write(os.path.join(out_dir, 'Myeloid.h5ad'))


# 对Fibroblasts注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Main/Main.h5ad")
# 'Tumor-like cells', 'T/NK', 'Myeloid cells', 'Fibroblasts', 'Neutrophils', 'Endothelial cells', 'Basophils'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Fibroblasts/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Fibroblasts", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.4,n_pcs = 20,output=out_dir)
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
    '0':'iCAF',
    '1':'mCAF',
    '2':'myCAF',
    '3':'mCAF',
    '4':'prolCAF',
    '5':'vCAF',
    '6':'devCAF'
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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':14, 'h':5})
adata.write(os.path.join(out_dir, 'Fibroblasts.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)


# 对Neutrophils注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Main/Main.h5ad")
# 'Tumor-like cells', 'T/NK', 'Myeloid cells', 'Fibroblasts', 'Neutrophils', 'Endothelial cells', 'Basophils'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Neutrophils/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Neutrophils", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.3,n_pcs = 10,output=out_dir)
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

Neutrophilmarker=['Fut4','Itgam','Itga4','Cd34','Cd38','Cd101','Fcgr4','Tfrc','Kit','Cxcr4','Cxcr2']
sc.pl.dotplot(adata, var_names=Neutrophilmarker, groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Neutrophil.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_Neutrophil.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=Neutrophilmarker, cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Neutrophil_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Neutrophil_diffgenes.pdf'))

anno_dict={
    '0':'Gbp2 Neutrophils',
    '1':'Cxcl3 Neutrophils',
    '2':'S100a9 Neutrophils',
    '3':'lowQ'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

adata=adata[~adata.obs.cell_type.isin(["lowQ"]), :]
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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':14, 'h':5})
adata.write(os.path.join(out_dir, 'Neutrophils.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2508294648/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/mapping_raw/"
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Fibroblasts/Fibroblasts.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Myeloid/Myeloid.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/T_NK/T_NK.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509225162-02/01.Mission/Tumor-like cells/Tumor.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['Neutrophils']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='on data',legend_fontsize=5)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'), facecolor='white', bbox_inches="tight", dpi=300)