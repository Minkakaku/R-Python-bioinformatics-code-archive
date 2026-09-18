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

# 对Myeloid cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWAE20251216/result.h5ad")
# ['B cells', 'Basophils', 'Erythroid cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/Myeloid cells/"
# adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['B cells', 'Plasma cells'])]
adata_sub = adata_raw[adata_raw.obs.cell_type=="Myeloid cells", :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.8,n_pcs = 10,output=out_dir,batch_effect = 2)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['CD14','FCGR3A','VCAN','DUSP1','ISG15','LYZ','APOE','C1QA','CLEC9A','DPYD','LRMDA'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Myeloid_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Myeloid_diffgenes.pdf'))
sc.pl.umap(adata, color=['CD3D','CD3G','CD3E'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_T_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_T_diffgenes.pdf'))

adata_sub = adata[~adata.obs.leiden.isin(['7']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=1,n_pcs = 10,output=out_dir,batch_effect = 3)
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

sc.pl.umap(adata, color=['CD14','FCGR3A','VCAN','DUSP1','ISG15','LYZ','APOE','C1QA','CLEC9A','DPYD','LRMDA'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Myeloid_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Myeloid_diffgenes.pdf'))

sc.pl.umap(adata, color=["LRMDA",'sample'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_sample.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_sample.pdf'),bbox_inches='tight')

genes = ['MTOR','RPTOR','RICTOR','RHEB','TSC1','TSC2']
term='mTOR pathway'
sc.tl.score_genes(adata, gene_list=genes, score_name=term)
groups = adata.obs['group'].cat.categories.to_list()
fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups) * 6, 6)) if len(groups) > 1 else plt.subplots(figsize=(6, 6))
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups)))) #加上这两句没问题了
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*6, 6))
    # Generate and save UMAP for each group
    for j in range(len(groups)):
        tmp = adata[adata.obs["group"] == groups[j]]
        sc.pl.umap(tmp, size=20, color=term, title=groups[j], ax=axs[j], show=False) #
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_{term}_split_groups', fig_size={'w': len(groups) * 6, 'h': 6})

sc.pl.umap(adata, color=['mTOR pathway'], cmap='Reds',vmin=0)
plt.savefig(os.path.join(out_dir, 'mTOR_umap.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'mTOR_umap.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata, var_names=['mTOR pathway'], groupby=['leiden','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'mTOR_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'mTOR_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
# sc.pl.violin(adata, ['mTOR pathway'], jitter=0, groupby='group', rotation=45, multi_panel=True, show=False)
# plt.savefig(os.path.join(out_dir, 'mTOR_violin.png'), facecolor='white', bbox_inches="tight", dpi=300)
# plt.savefig(os.path.join(out_dir, 'mTOR_violin.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.obs['group'] = adata_raw.obs['group'].replace('NC', 'HD')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('iMCD', 'P')
adata_raw.obs['group'] = adata_raw.obs['group'].replace('iMCD_after', 'P-after')
adata.obs['group'] = adata.obs['group'].replace('NC', 'HD')
adata.obs['group'] = adata.obs['group'].replace('iMCD', 'P')
adata.obs['group'] = adata.obs['group'].replace('iMCD_after', 'P-after')

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
    '10':'ISG15 Classical monocyte',
    '11':'VCAN Classical monocyte',
    '12':'cDCs',
    '13':'VCAN Classical monocyte',
    '14':'LRMDA Classical monocyte',
    '15':'LRMDA Classical monocyte',
    '16':'pDCS',
    '17':'LRMDA Classical monocyte',
    '18':'Non-classical monocyte',
    '19':'ISG15 Classical monocyte',
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


# 对B/Plasma cells重新注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWAE20251216/result.h5ad")
# ['B cells', 'Basophils', 'Erythroid cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/B-Plasma/"
adata_sub = adata_raw[adata_raw.obs.cell_type.isin(['B cells','Plasma cells']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=1.2,n_pcs = 10,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

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
sc.pl.dotplot(adata, var_names=['PCNA','CD19','MS4A1','CD27','PTPRC','CD38','CD40','CD83','TNFRSF13B','TNFRSF13C','CD69'], groupby=['leiden'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'ActivatedB_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'ActivatedB_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_sub = adata[~adata.obs.leiden.isin(['4','5']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=1,n_pcs = 10,output=out_dir,batch_effect = 3)
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

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'Naive B',
    '1':'Naive B',
    '2':'Memory B C1',
    '3':'Memory B C1',
    '4':'Memory B C2',
    '5':'Naive B',
    '6':'Memory B C2',
    '7':'Naive B',
    '8':'Activated B',
    '9':'Activated B',
    '10':'Memory B C2',
    '11':'Plasma',
    '12':'Plasma',
    '13':'Activated B',
    '14':'Memory B C2',
    '15':'Plasma',
    '16':'Plasma',
    '17':'Memory B C1',
    '18':'Plasma',
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_map = adata[adata.obs['cell_type'].isin(['Plasma']), :]
sc.tl.leiden(adata_map, resolution=0.8)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_map, color='leiden', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'IgA Plasma',
    '1':'IgA Plasma',
    '2':'IgA Plasma',
    '3':'IgG Plasma C1',
    '4':'IgG Plasma C2',
    '5':'IgM Plasma',
    '6':'IgA Plasma',
    '7':'IgG Plasma C2',
    '8':'IgA Plasma',
    '9':'IgA Plasma',
    '10':'IgA Plasma',
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

adata.obs['group'] = adata.obs['group'].replace('NC', 'HD')
adata.obs['group'] = adata.obs['group'].replace('iMCD', 'P')
adata.obs['group'] = adata.obs['group'].replace('iMCD_after', 'P-after')
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


# 对T/NK注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWAE20251216/result.h5ad")
# 'B cells', 'Erythroid cells', 'Mast cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/T_NK/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="T/NK", :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.3,n_pcs = 10,output=out_dir,batch_effect = 3)
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

sc.pl.umap(adata, color=['CD4','CD8A','FOXP3','GZMA','NKG7','CD3D','CD3E','NCAM1','FCGR3A','CTLA4','IL2RA'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata, var_names=['CD4','CD8A','FOXP3','GZMA','NKG7','CD3D','CD3E','NCAM1','FCGR3A','CTLA4','IL2RA'], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

anno_dict={
    '0':'CD4T',
    '1':'NK',
    '2':'CD4T',
    '3':'CD8T',
		'4':'CD8T',
    '5':'CD4T',
    '6':'CD4T',
    '7':'CD8T',
    '8':'CD4T',
    '9':'NK',
    '10':'CD4T'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata.obs['group'] = adata.obs['group'].replace('NC', 'HD')
adata.obs['group'] = adata.obs['group'].replace('iMCD', 'P')
adata.obs['group'] = adata.obs['group'].replace('iMCD_after', 'P-after')
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
adata.write(os.path.join(out_dir, 'T_NK.h5ad'))

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWAE20251216/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/mapping_raw/"

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/Myeloid cells/Myeloid.h5ad")
bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/B-Plasma/B_Plasma.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_sub = adata_raw[~adata_raw.obs.sub_cell_type.isin(['B cells','Myeloid cells']), :]
sc.pl.umap(adata_sub, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

# sc.tl.umap(adata_sub)
# sc.tl.leiden(adata_sub, resolution=0.5)
# plot cluster
# sc.pl.umap(adata_sub, color='leiden', show=False)
# plt.savefig(os.path.join(out_dir, 'sub_cluster_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
# plt.savefig(os.path.join(out_dir, 'sub_cluster_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_sub, color='cell_type', show=False)
plt.savefig(os.path.join(out_dir, 'cell_type_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
# sc.pl.umap(adata_sub, color='sub_cell_type', show=False)
# plt.savefig(os.path.join(out_dir, 'sub_cell_type_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
# plt.savefig(os.path.join(out_dir, 'sub_cell_type_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_sub, color='cell_type', show=False,legend_loc ='on data')
plt.savefig(os.path.join(out_dir, 'cell_type_umap_withlabel.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type_umap_withlabel.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_sub.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_sub.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
method = "cell_type"
keys = adata_sub.obs[method].cat.categories.to_list()
values = list(adata_sub.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

adata_sub.obs['group'] = adata_sub.obs['group'].replace('NC', 'HD')
adata_sub.obs['group'] = adata_sub.obs['group'].replace('iMCD', 'P')
adata_sub.obs['group'] = adata_sub.obs['group'].replace('iMCD_after', 'P-after')
groups = ['HD', 'P', 'P-after']
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

sc.pl.umap(adata_sub, color=['MIF','CD74','CXCR2','CXCR4','ACKR3'], cmap='Reds',)
plt.savefig(os.path.join(out_dir, 'umap_MIFgenes.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_MIFgenes.pdf'), facecolor='white', bbox_inches="tight")


adata = sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/12.Mission/B-Plasma/B_Plasma.h5ad")
adata = adata[~adata.obs.cell_type.isin(['Activated B', 'Memory B C1', 'Memory B C2', 'Naive B']), :]
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/geneplot/"
sc.settings.dotplot_size = 150
sc.pl.dotplot(adata, var_names=['IL6','MIF'], groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type-group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.settings.dotplot_size =300
sc.pl.dotplot(adata, var_names=['IL6','MIF'], groupby=['group'], dendrogram=False, show=False)
# save_fig(plt, out_dir=os.path.join(out_dir), file_name='group_bubble', fig_size={'w': 6, 'h': 6})
plt.savefig(os.path.join(out_dir, 'group_bubble.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'group_bubble.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


gene_use = ['ACKR3'] #['MIF','CD74','CXCR2','CXCR4','ACKR3']
outdir_tmp ="/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/geneplot/MIF/"
df = sc.get.obs_df(adata_sub, keys=gene_use, use_raw=True)  # cells × 1 (Irf5)
meta = adata_sub.obs
ax1 = sc.pl.dotplot(
 adata_sub,
    var_names=gene_use,  # 显示所有分组
    groupby=['cell_type','group'],  # 按细胞类型分组
    show=True,
    return_fig=True
)

expr_matrix = pd.DataFrame(
    np.nan,  # 初始化为np.nan
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
    adata_sub,
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


####之前的注释结果再做一下
adata_sub=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/05.Mission/mapping_raw/Main_mapping.h5ad")
outdir_tmp="/mnt/NexData_file/home/sh/Projects/FWAE20251216/13.Mission/geneplot/MIF_NPAID/"
sc.pl.umap(adata_sub, color=['cell_type'], cmap='Reds',)
plt.savefig(os.path.join(out_dir, 'umap_MIFgenes.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_MIFgenes.pdf'), facecolor='white', bbox_inches="tight")

adata_sub.obs['group'] = adata_sub.obs['group'].replace('NC', 'HD')
adata_sub.obs['group'] = adata_sub.obs['group'].replace('iMCD', 'P')
adata_sub.obs['group'] = adata_sub.obs['group'].replace('iMCD_after', 'P-after')
adata_sub.obs['raw_cell_type']=adata_sub.obs['cell_type']
adata_sub.obs['cell_type']=adata_sub.obs['group']
adata_sub.obs['group']=adata_sub.obs['raw_cell_type']
gene_use = ['ACKR3'] #['MIF','CD74','CXCR2','CXCR4','ACKR3']
df = sc.get.obs_df(adata_sub, keys=gene_use, use_raw=True)  # cells × 1 (Irf5)
meta = adata_sub.obs
ax1 = sc.pl.dotplot(
 adata_sub,
    var_names=gene_use,  # 显示所有分组
    groupby=['cell_type','group'],  # 按细胞类型分组
    show=True,
    return_fig=True
)

expr_matrix = pd.DataFrame(
    np.nan,  # 初始化为np.nan
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
    adata_sub,
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
