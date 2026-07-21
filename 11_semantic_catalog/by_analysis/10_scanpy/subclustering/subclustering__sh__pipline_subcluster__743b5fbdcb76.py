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

# 调整Neutrophils注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2511266943-04/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Neutrophils', 'GMP', 'Myeloid cells', 'T cells', 'Erythrocytes', 'HSC', 'Plasma cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2511266943-05/04.Mission/Neutrophils/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Neutrophils", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['Cd34','Sox4','Elane','Mpo','Prtn3','Rpl12','Camp','Chil3','Fcnb','Ltf','Ngp','Tuba1b','Camp','Cxcl2','Mmp8','Ngp','Retnlg','Ccl6','Cxcl2','Mmp8','Retnlg','Ccl6','Cstdc4','Stfa2l1','Ifit3','Isg15','Rsad2','Fgl2','Gm2a','Gngt2'], cmap='coolwarm', ncols=3)
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['Cd34','Sox4','Elane','Mpo','Prtn3','Rpl12','Camp','Chil3','Fcnb','Ltf','Ngp','Tuba1b','Camp','Cxcl2','Mmp8','Ngp','Retnlg','Ccl6','Cxcl2','Mmp8','Retnlg','Ccl6','Cstdc4','Stfa2l1','Ifit3','Isg15','Rsad2','Fgl2','Gm2a','Gngt2'], groupby='leiden', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False,legend_loc='on data',legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

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
    '0':'Neutrophil_G3',
    '1':'Neutrophil_G2',
    '2':'Neutrophil_G5',
    '3':'Neutrophil_G4',
    '4':'Neutrophil_G5',
    '5':'Neutrophil_G1',
    '6':'Neutrophil_G4',
    '7':'Neutrophil_G2',
    '8':'Neutrophil_G1',
    '9':'lowQ'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

adata = adata[~adata.obs['cell_type'].isin(["lowQ"]), :]
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Neutrophils.h5ad'))

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

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 调整Myeloid cells注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2511266943-04/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Neutrophils', 'GMP', 'Myeloid cells', 'T cells', 'Erythrocytes', 'HSC', 'Plasma cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2511266943-05/04.Mission/Myeloid/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Myeloid cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['Lyz2','Apoe','Cd14','Fcgr3','Vcan','Mki67','Top2a','Cd86','Mrc1','Clec9a','Cd209a'], cmap='coolwarm', ncols=3)
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['Lyz2','Apoe','Cd14','Fcgr3','Vcan','Mki67','Top2a','Cd86','Mrc1','Clec9a','Cd209a'], groupby='leiden', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Ifitm3 Monocytes',
    '1':'Ifitm3 Monocytes',
    '2':'Myeloid progenitor cells',
    '3':'Cycling Monocytes',
    '4':'cDC',
    '5':'pDC',
    '6':'Ifitm3 Monocytes',
    '7':'cDC',
    '8':'Non-classical Monocyte',
    '9':'cDC',
    '10':'Cycling DC',
    '11':'Lrmda Monocyte',
    '12':'Lrmda Monocyte'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Myeloid.h5ad'))

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

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 调整T cells注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2511266943-04/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'Neutrophils', 'GMP', 'Myeloid cells', 'T cells', 'Erythrocytes', 'HSC', 'Plasma cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2511266943-05/04.Mission/Tcells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="T cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.8,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()

adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['Cd4','Cd8a','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Sell','Ccr7','Lef1','Cd27','Cd28','Prf1','Rgs1','Cd69','Tbx21','Tnf','Rorc','Il23r','Bcl6','Batf','Klrg1','Ifng','Tox','Pdcd1'], cmap='coolwarm', ncols=3)
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['Cd4','Cd8a','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Sell','Ccr7','Lef1','Cd27','Cd28','Prf1','Rgs1','Cd69','Tbx21','Tnf','Rorc','Il23r','Bcl6','Batf','Klrg1','Ifng','Tox','Pdcd1'], groupby='leiden', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Treg',
    '1':'CD8T naive',
    '2':'NKT',
    '3':'CD8T memory',
    '4':'CD8T cytotoxic',
    '5':'NK',
    '6':'CD8T naive',
    '7':'lowQ',
    '8':'Th17'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

adata = adata[~adata.obs['cell_type'].isin(["lowQ"]), :]
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Tcells.h5ad'))

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

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 调整B cells注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2511266943-04/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Epithelial cells', 'B cells', 'T cells', 'Macrophages', 'Enterocytes', 'Mesothelial cells', 'Endothelial cells', 'Fibroblasts', 'Paneth cells'  'Tuft cells', 'SMCs', 'Plasma cells', 'Schwann cells']
out_dir="/mnt/NexData_file/home/sh/Projects/SC2511266943-05/03.Mission/Bcells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="B cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.2,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['Cd24a','Dntt','Vpreb1','Cd19','Cd24a','Dntt','Igll1','Vpreb1','Rag1'], cmap='coolwarm', ncols=3)
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['Cd24a','Dntt','Vpreb1','Cd19','Cd24a','Dntt','Igll1','Vpreb1','Rag1'], groupby='leiden', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Mature B',
    '1':'Plasmablast',
    '2':'Pre-B',
    '3':'Pro-Pre-B',
    '4':'Immature B',
    '5':'Pro-B',
    '6':'lowQ'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

adata = adata[~adata.obs['cell_type'].isin(["lowQ"]), :]
sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Bcells.h5ad'))

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

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2511266943-04/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2511266943-05/04.Mission/Main_mapping/"

adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2511266943-05/04.Mission/Bcells/Bcells.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2511266943-05/04.Mission/Neutrophils/Neutrophils.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2511266943-05/04.Mission/Tcells/Tcells.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['B cells','T cells','Neutrophils']), :]
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
# mask = adata_raw.obs['sub_cell_type'].isin(['B cells', 'Pre-B'])          # 找到 Monocytes 行
# adata_raw.obs.loc[mask, 'cell_type'] = 'B cells'                    # 把对应 cell_type 改为 A
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping_filter.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')

method = "group"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')

# out="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Macrophages/cluster2/"
# species='human'
# os.system(f'/mnt/NexData_file/share/scRNA-tools/envs/pip_diffgene_enrich/bin/Rscript /mnt/NexData_file/home/lst/scrna-code/diff/do_enrichment.R -f {os.path.join(out, "differentially_expressed_genes.xls")} -s {species} -o /mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Macrophages/cluster2/')