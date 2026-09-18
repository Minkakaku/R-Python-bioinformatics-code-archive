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

# 对Proximal tubule cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWAE20254246/adata.h5ad")
# AnnData object with n_obs × n_vars = 25540923 × 36601
# obs: 'sample', 'original_barcode', 'sample_id_x', 'group_id_x', 'patient_id_x', 'level1_x', 'level2_x', 'level3_x', 'sample_id_y', 'group_id_y', 'patient_id_y', 'l_y', 'level2_y', 'level3_y'
# 删除：sample_id_x', 'group_id_x', 'level1_x', 'level2_x', 'level3_x'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20254246/00.data/Main/"
sc.pl.violin(adata_raw, ['sample', 'l_y'], jitter=0, groupby='sample', rotation=45, multi_panel=True, show=False)
axes = plt.gcf().axes
for ax in axes:
    ax.set_xticklabels(ax.get_xticklabels(), horizontalalignment='right')
save_fig(plt, out_dir, file_name='orig_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})

sc.pl.umap(adata_raw, color=['sample'], cmap='coolwarm',show=False)
sc.pl.violin(adata_raw, ['sample'], jitter=0, groupby='sample', rotation=45, multi_panel=True, show=False)
plt.savefig("/mnt/NexData_file/home/sh/Projects/FWAE20254246/00.data/Main/umnp.png")
plt.savefig("/mnt/NexData_file/home/sh/Projects/FWAE20254246/00.data/Main/umap.pdf")

adata_sub = adata_raw[adata_raw.obs.cell_type=="Proximal tubule cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=30,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 3)
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

sc.pl.umap(adata, color=['Miox','Ass1','Fxyd2','Fbp1','Gsta2','Slc5a12','Slc7a13','Havcr1','Mki67','Slc22a30'], cmap='coolwarm',show=False)
plt.savefig(os.path.join(out_dir, 'umap_Brush border cell.png'))
plt.savefig(os.path.join(out_dir, 'umap_Brush border cell.pdf'))
sc.pl.dotplot(adata, var_names=['Lrp2','Hnf4a','Slc5a12','Slc5a2','Slc13a1','Slco1a1','Ugt2b37','Slc7a13','Cyp7b1','Slc22a7','Slc7a12','Spag5','Mki67','Havcr1','Cdh13','Cdh6','Cald1','Pdgfd','Kcnip4','Vcam1','Ccl2'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_PTC.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_PTC.pdf'), facecolor='white', bbox_inches="tight")

anno_dict={
    '0':'PTCs1',
    '1':'PTCs4',
    '2':'PTCs2',
    '3':'PTCs3',
    '4':'PTCs3',
    '5':'PTCs4',
    '6':'PTCs4',
    '7':'PTCs4',
    '8':'PTCs2'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type',method='wilcoxon',pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

anno_dict={
    '0':'Pde4d PTC',
    '1':'PTC',
    '2':'Slc7a13 PTC',
    '3':'Gbp2 PTC',
    '4':'Gbp2 PTC',
    '5':'PTC',
    '6':'PTC',
    '7':'PTC',
    '8':'Slc7a13 PTC'
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

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
adata.obs['group'] = adata.obs['group'].cat.reorder_categories(groups, ordered=True)
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'PTC.h5ad'))


# 映射所有细分
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/01.Mission/mapping_raw/"
adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin', show=False)
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))


# 对Distal tubule cells注释细分
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/01.Mission/Distal tubule cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Distal tubule cells", :]
adata_sub.obs['group'] = adata_sub.obs['group'].cat.reorder_categories(groups, ordered=True)
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=30,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 3)
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
    '0':'DTCs1',
    '1':'DTCs2',
    '2':'DTCs4',
    '3':'DTCs2',
    '4':'DTCs3',
    '5':'DTCs1',
    '6':'DTCs1',
    '7':'DTCs4'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type',method='wilcoxon',pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

anno_dict={
    '0':'DTC',
    '1':'Calb1 DTC',
    '2':'Erbb4 DTC',
    '3':'Calb1 DTC',
    '4':'Trpm6 DTC',
    '5':'DTC',
    '6':'DTC',
    '7':'Erbb4 DTC'
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

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
adata.obs['group'] = adata.obs['group'].cat.reorder_categories(groups, ordered=True)
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'DTC.h5ad'))


# 映射所有细分
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/01.Mission/mapping_raw/"
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin', show=False)
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))


# 提取ARG2表达
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/01.Mission/ARG2/"
adata_sub = adata_raw[adata_raw.obs.cell_type.isin(["Proximal tubule cells","Distal tubule cells"]), :]
adata_sub.obs['group'] = adata_sub.obs['group'].cat.reorder_categories(groups, ordered=True)
adata = adata_sub
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_subcelltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_subcelltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.umap(adata, color=['Arg2'], cmap='coolwarm',show=False)
plt.savefig(os.path.join(out_dir, 'umap_Arg2.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_Arg2.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.umap(adata, color=['Arg2'], cmap='coolwarm',show=False,vmax=1)
plt.savefig(os.path.join(out_dir, 'umap_Arg2_fiter.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_Arg2_fiter.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata, var_names=['Arg2'], groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Arg2-cell_type.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_Arg2-cell_type.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata, var_names=['Arg2'], groupby='sub_cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Arg2-sub_cell_type.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_Arg2-sub_cell_type.pdf'), facecolor='white', bbox_inches="tight")

adata = adata_raw[adata_raw.obs.cell_type.isin(["Distal tubule cells"]), :]
adata=adata.raw.to_adata()
expr = adata[:, 'Arg2'].X.toarray().ravel()
adata.obs['Arg2_gene'] = pd.Categorical(
    np.where(expr >= 1, 'high', 'low'),
    categories=['high', 'low']          # 固定顺序，防止后续错位
)
adata.uns['Arg2_gene_colors'] = ['#1f77b4', 'lightgray']
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='Arg2_gene', title='Arg2_gene', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Arg2_gene.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Arg2_gene.pdf'))

groups = ['H0', 'H2', 'H8', 'H24', 'H72']
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="Arg2_gene", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_Arg2_gene.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_Arg2_gene.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.tl.rank_genes_groups(adata, groupby='Arg2_gene',method='wilcoxon',pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
adata.obs['group'] = adata.obs['group'].cat.reorder_categories(groups, ordered=True)
adata.write(os.path.join(out_dir, 'Arg2_gene.h5ad'))
adata.obs['cell_type'] = adata.obs['Arg2_gene'] 
adata.uns['cell_type_colors'] = adata.uns['Arg2_gene_colors'] 
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
