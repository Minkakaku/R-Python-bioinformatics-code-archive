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
# 'B cells', 'Erythroid cells', 'Mast cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/Myeloid cells/"
# adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['B cells', 'Plasma cells'])]
adata_sub = adata_raw[adata_raw.obs.cell_type=="Myeloid cells", :]
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.3,n_pcs = 15,output=out_dir,batch_effect = 3)
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

sc.pl.umap(adata, color=['CD14','FCGR3A','VCAN','DUSP1','ISG15','LYZ','APOE','C1QA','CLEC9A','DPYD','LRMDA'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Myeloid_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_Myeloid_diffgenes.pdf'))

adata_sub = adata[~adata.obs['leiden'].isin(['6']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 3)
adata = adata_recluster.sub_cluster()

sc.tl.dendrogram(adata, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight', dpi=300)
plt.close()

anno_dict={
    '0':'VCAN Classical monocyte',
    '1':'VCAN Classical monocyte',
    '2':'Non-classical monocyte',
    '3':'VCAN Classical monocyte',
		'4':'DUSP1 Classical monocyte',
    '5':'LRMDA Classical monocyte',
    '6':'ISG15 Classical monocyte',
    '7':'VCAN Classical monocyte',
    '8':'cDCs',
    '9':'VCAN Classical monocyte',
    '10':'pDCs',
    '11':'Non-classical monocyte',
    '12':'DUSP1 Classical monocyte',
    '13':'cDCs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['NC', 'iMCD', 'iMCD_after']
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
# 'B cells', 'Erythroid cells', 'Mast cells', 'Megakaryocyte', 'Myeloid cells', 'Neutrophils', 'Plasma cells', 'T/NK'
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/B-Plasma/"
adata_sub = adata_raw[adata_raw.obs.cell_type.isin(['B cells','Plasma cells']), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=15,resolution=0.3,n_pcs = 10,output=out_dir,batch_effect = 3)
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

adata_sub = adata[~adata.obs.leiden.isin(['2']), :]
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

anno_dict={
    '0':'Naive B',
    '1':'Memory B C1',
    '2':'Memory B C1',
    '3':'Naive B',
    '4':'Naive B',
    '5':'Memory B C2',
    '6':'Naive B',
    '7':'Naive B',
    '8':'Activated B',
    '9':'Naive B',
    '10':'Memory B C1',
    '11':'IgA Plasma',
    '12':'IgA Plasma',
    '13':'Activated B',
    '14':'Memory B C2',
    '15':'Naive B',
    '16':'IgA Plasma',
    '17':'Memory B C1',
    '18':'IgA Plasma',
    '19':'IgA Plasma',
    '20':'Memory B C1'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_map = adata[adata.obs['cell_type'].isin(['IgA Plasma']), :]
sc.tl.leiden(adata_map, resolution=0.5)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_map, color='leiden', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

anno_dict={
    '0':'IgA Plasma',
    '1':'IgA Plasma',
    '2':'IgA Plasma',
    '3':'IgG Plasma C1',
    '4':'IgA Plasma',
    '5':'IgG Plasma C2',
    '6':'IgG Plasma C2',
    '7':'IgM Plasma',
    '8':'IgA Plasma',
    '9':'IgA Plasma'
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

groups = ['NC', 'iMCD', 'iMCD_after']
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
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/T_NK/"
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
    '2':'CD8T',
    '3':'CD8T',
		'4':'CD4T',
    '5':'CD4T',
    '6':'NK',
    '7':'NK'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

groups = ['NC', 'iMCD', 'iMCD_after']
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
out_dir="/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/mapping_raw/"

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/Myeloid cells/Myeloid.h5ad")
bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/T_NK/T_NK.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_sub = adata_raw[~adata_raw.obs.sub_cell_type.isin(['2','Myeloid cells']), :]
sc.pl.umap(adata_sub, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_map = adata_sub[adata_sub.obs['cell_type'].isin(['Mast cells']), :]
sc.tl.leiden(adata_map, resolution=0.3)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_map, color='leiden', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.umap(adata_sub, color=['CD34'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'), facecolor='white', bbox_inches="tight")

anno_dict={
    '0':'HSPCs',
    '1':'Mast cells',
    '2':'Mast cells',
    '3':'Mast cells',
    '4':'HSPCs'
}
adata_map.obs['cell_type'] = [anno_dict[clust] for clust in adata_map.obs['leiden']]
adata_map.obs['cell_type'] = adata_map.obs['cell_type'].astype('category')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_map, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cell_type_mapping = adata_map.obs[['cell_type']].to_dict()['cell_type']
adata_sub.obs['cell_type'] = adata_sub.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_sub.obs['sub_cell_type'] = adata_sub.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/Myeloid cells/Myeloid.h5ad")
bdata=sc.read("/mnt/NexData_file/home/sh/Projects/FWAE20251216/11.Mission/T_NK/T_NK.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_sub.obs['sub_cell_type'] = adata_sub.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_sub.obs['sub_cell_type'] = adata_sub.obs['sub_cell_type'].astype('category')
adata_sub.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_sub, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.tl.umap(adata_sub)
sc.tl.leiden(adata_sub, resolution=0.5)
# plot cluster
sc.pl.umap(adata_sub, color='leiden', show=False)
plt.savefig(os.path.join(out_dir, 'sub_cluster_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'sub_cluster_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_sub, color='cell_type', show=False)
plt.savefig(os.path.join(out_dir, 'cell_type_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'cell_type_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_sub, color='sub_cell_type', show=False)
plt.savefig(os.path.join(out_dir, 'sub_cell_type_umap.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'sub_cell_type_umap.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
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