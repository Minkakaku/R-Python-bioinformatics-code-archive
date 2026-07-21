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

# 重新大类注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094836-05/03.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'B cells', 'Basal cells', 'Distal tubule cells', 'Ductal cells', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Mast cells', 'Myeloid cells', 'NK cells', 'Neutrophils', 'Principal cells', 'Proximal tubule cells', 'Smooth muscle cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094836-05/04.Mission/Main/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,legend_loc='on data',legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cluster', title='Cluster', show=False,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_cluster.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_cluster.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

cell_markers = {
    "POD":['PTPRQ','WT1','NTNG1','NPHS1','NPHS2','CLIC5','PODXL'],
    "PEC":['CLDN1','VCAM1','CFH','RBFOX1','ALDH1A2'],
    "PT":['LRP2','CUBN','SLC13A1'],
    "TL":['CRYAB','TACSTD2','SLC44A5','KLRG2','COL26A1','BOC'],
    "TAL":['CASR','SLC12A1','UMOD'],
    "DCT":['SLC12A3','CNNM2','FGF13','KLHL3','LHX1','TRPM6'],
    "CNT":['SLC8A1','SCN2A','HSD11B2','CALB1'],
    "PC":['GATA3','AQP2','AQP3'],
    "IC":['ATP6V0D2','ATP6V1C2','TMEM213','CLNK'],
    "EC":['CD34','PECAM1','PTPRB','MEIS2','FLT1','EMCN'],
    "VSM/P":['NOTCH3','PDGFRB','ITGA8'],
    "FIB":['COL1A1','COL1A2','C7','NEGR1','FBLN5','DCN','CDH11'],
    "IMM":['PTPRC'],
}
markers=['PTPRQ','WT1','NTNG1','NPHS1','NPHS2','CLIC5','PODXL','CLDN1','VCAM1','CFH','RBFOX1','ALDH1A2','LRP2','CUBN','SLC13A1','CRYAB','TACSTD2','SLC44A5','KLRG2','COL26A1','BOC','CASR','SLC12A1','UMOD','SLC12A3','CNNM2','FGF13','KLHL3','LHX1','TRPM6','SLC8A1','SCN2A','HSD11B2','CALB1','GATA3','AQP2','AQP3','ATP6V0D2','ATP6V1C2','TMEM213','CLNK','CD34','PECAM1','PTPRB','MEIS2','FLT1','EMCN','NOTCH3','PDGFRB','ITGA8','COL1A1','COL1A2','C7','NEGR1','FBLN5','DCN','CDH11','PTPRC']
sc.pl.umap(adata_raw, color=markers, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata_raw, var_names=cell_markers, groupby='cluster', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

anno_dict={
    '0':'EC',
    '1':'EC',
    '2':'VSM/P',
    '3':'IMM',
    '4':'TL',
    '5':'TAL',
    '6':'FIB',
    '7':'TL',
    '8':'IMM',
    '9':'TAL',
    '10':'PC',
    '11':'IMM',
    '12':'PC',
    '13':'EC',
    '14':'FIB',
    '15':'VSM/P',
    '16':'PC',
    '17':'IMM',
    '18':'IMM',
    '19':'IMM',
    '20':'IMM',
    '21':'PT'
}
adata_raw.obs['cell_type'] = [anno_dict[clust] for clust in adata_raw.obs['cluster']]
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')

adata=sc.read("/mnt/NexData_file/home/sh/Projects/SC2509094836-05/03.Mission/mapping_raw/Main_mapping.h5ad")
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sub_cell_type', title='Cell Type', show=False,legend_loc='on data',legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.add_categories('IC')
mask = adata.obs['sub_cell_type'] == 'IC'
adata_raw.obs.loc[mask, 'cell_type'] = 'IC'
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.remove_unused_categories()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,legend_loc='on data',legend_fontsize=8)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

groups = adata_raw.obs['group'].drop_duplicates().tolist()
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

adata_raw.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata_raw.write(os.path.join(out_dir, 'Main.h5ad'))

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

method = "group"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')

cross_tab = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv(os.path.join(out_dir, "celltype_proportions_by_group.csv"))

cell_type_markers = {
    "EC":['CD34','PECAM1','PTPRB','MEIS2','FLT1','EMCN'],
    "FIB":['COL1A1','COL1A2','C7','NEGR1','FBLN5','DCN','CDH11'],
    "IMM":['PTPRC'],
    "PC":['GATA3','AQP2','AQP3'],
    "PT":['LRP2','CUBN','SLC13A1'],
    "TAL":['CASR','SLC12A1','UMOD'],
    "TL":['CRYAB','TACSTD2','SLC44A5','KLRG2','COL26A1','BOC'],
    "VSM/P":['NOTCH3','PDGFRB','ITGA8'],
    "IC":['ATP6V0D2','ATP6V1C2','TMEM213','CLNK']
}
sc.pl.dotplot(adata_raw, var_names=cell_type_markers, groupby='cell_type', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), bbox_inches="tight")


# 对PC注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094836-05/04.Mission/Main/Main.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'EC', 'FIB', 'IMM', 'PC', 'PT', 'TAL', 'TL', 'VSM/P', 'IC'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094836-05/04.Mission/PC/"
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["PC"]), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

cell_markers = {
    "CCD-PC":["SCNN1G", "SCNN1B", "FXYD4", "SOX5", "PDE10A", "SLC25A29", "ST6GAL1", "PAPPA"],
    "IMCD":["FXYD4", "SOX5", "PHACTR1", "PCDH7", "SLC14A2", "HS3ST5"],
    "OMCD-PC":["SCNN1G", "SCNN1B", "FXYD4", "SOX5", "SYK", "FAM81A", "PROM1", "KCNK13"],
    "PapE":["TACSTD2", "TP63", "GPX2", "FXYD3", "KRT5"]
}
markers=["SCNN1G", "SCNN1B", "FXYD4", "SOX5", "PDE10A", "SLC25A29", "ST6GAL1", "PAPPA","SYK", "FAM81A", "PROM1", "KCNK13","PHACTR1", "PCDH7", "SLC14A2", "HS3ST5","TACSTD2", "TP63", "GPX2", "FXYD3", "KRT5"]
sc.pl.umap(adata, color=markers, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

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
    '0':'OMCD-PC',
    '1':'PapE',
    '2':'CCD-PC',
    '3':'IMCD',
    '4':'OMCD-PC',
    '5':'IMCD'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_loc='on data',legend_fontsize=8)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'PC.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

cross_tab = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv(os.path.join(out_dir, "celltype_proportions_by_group.csv"))

sc.pl.dotplot(adata, var_names=cell_markers, groupby='cell_type', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), bbox_inches="tight")


# 对EC注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094836-05/04.Mission/Main/Main.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'EC', 'FIB', 'IMM', 'PC', 'PT', 'TAL', 'TL', 'VSM/P', 'IC'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094836-05/04.Mission/EC/"
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["EC"]), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

cell_markers = {
    "EC-GC":["EMCN","HECW2","PLAT","ITGA8","EHD3","KDR","SOST"],
    "EC-AEA":["BTNL9","ADAMTS6","PALMD","AQP1","TM4SF1","VEGFC","CCDC3","CDH5","SERPINE2","FBLN5","CXCL12","SOX17"],
    "EC-DVR":["BTNL9","ADAMTS6","PALMD","AQP1","TM4SF1","MCTP1","SLC14A1","ENPP2","LYPD6B"],
    "EC-PTC":["CEACAM1","DNASE1L3","PLVAP","PITPNC1","GRB10","SLCO2A1","RAPGEF4"],
    "EC-AVR":["CEACAM1","DNASE1L3","PLVAP","GPM6A","EDIL3","TLL1","ZNF385D","NR2F2"],
    "EC-LYM":["MMRN1","CD36","TBX1","PKHD1L1","PROX1"]
}
markers=["EMCN","HECW2","PLAT","ITGA8","EHD3","KDR","SOST","BTNL9","ADAMTS6","PALMD","AQP1","TM4SF1","VEGFC","CCDC3","CDH5","SERPINE2","FBLN5","CXCL12","SOX17","BTNL9","ADAMTS6","PALMD","AQP1","TM4SF1","MCTP1","SLC14A1","ENPP2","LYPD6B","CEACAM1","DNASE1L3","PLVAP","PITPNC1","GRB10","SLCO2A1","RAPGEF4","CEACAM1","DNASE1L3","PLVAP","GPM6A","EDIL3","TLL1","ZNF385D","NR2F2","MMRN1","CD36","TBX1","PKHD1L1","PROX1"]
sc.pl.umap(adata, color=markers, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

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
    '0':'EC-DVR',
    '1':'EC-GC',
    '2':'EC-AEA',
    '3':'EC-GC',
    '4':'EC-GC',
    '5':'EC-AVR',
    '6':'EC-PTC'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_loc='on data',legend_fontsize=8)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'EC.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

cross_tab = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv(os.path.join(out_dir, "celltype_proportions_by_group.csv"))

cell_type_markers = {
    "EC-AEA":["BTNL9","ADAMTS6","PALMD","AQP1","TM4SF1","VEGFC","CCDC3","CDH5","SERPINE2","FBLN5","CXCL12","SOX17"],
    "EC-AVR":["CEACAM1","DNASE1L3","PLVAP","GPM6A","EDIL3","TLL1","ZNF385D","NR2F2"],
    "EC-DVR":["BTNL9","ADAMTS6","PALMD","AQP1","TM4SF1","MCTP1","SLC14A1","ENPP2","LYPD6B"],
    "EC-GC":["EMCN","HECW2","PLAT","ITGA8","EHD3","KDR","SOST"],
    "EC-PTC":["CEACAM1","DNASE1L3","PLVAP","PITPNC1","GRB10","SLCO2A1","RAPGEF4"]
}
sc.pl.dotplot(adata, var_names=cell_type_markers, groupby='cell_type', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), bbox_inches="tight")


# 对VSM/P注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094836-05/04.Mission/Main/Main.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'EC', 'FIB', 'IMM', 'PC', 'PT', 'TAL', 'TL', 'VSM/P', 'IC'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094836-05/04.Mission/VSM_P/"
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["VSM/P"]), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

cell_markers = {
    "MC":["PIP5K1B","ROBO1","PIEZO2","DAAM2","PHTF2","GATA3","POSTN"],
    "REN":["PIP5K1B","ROBO1","REN","PDE10A","ABCC8","COL13A1","GRID2"],
    "VSMC":["NTRK3","MYH11","RGS6","ADRA1A","LDB3","MCAM"],
    "VSMC/P":["NTRK3","CCDC102B","RGS5","ABCC9","ADCY3","ADGRB3"]
}
markers=["PIP5K1B","ROBO1","PIEZO2","DAAM2","PHTF2","GATA3","POSTN","PIP5K1B","ROBO1","REN","PDE10A","ABCC8","COL13A1","GRID2","NTRK3","MYH11","RGS6","ADRA1A","LDB3","MCAM","NTRK3","CCDC102B","RGS5","ABCC9","ADCY3","ADGRB3"]
sc.pl.umap(adata, color=markers, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

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
    '0':'VSMC/P',
    '1':'VSMC/P',
    '2':'REN',
    '3':'VSMC/P',
    '4':'VSMC',
    '5':'VSMC/P',
    '6':'MC'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_loc='on data',legend_fontsize=8)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'VSM_P.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

cross_tab = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv(os.path.join(out_dir, "celltype_proportions_by_group.csv"))

sc.pl.dotplot(adata, var_names=cell_markers, groupby='cell_type', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), bbox_inches="tight")


# 对FIB注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094836-05/04.Mission/Main/Main.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'EC', 'FIB', 'IMM', 'PC', 'PT', 'TAL', 'TL', 'VSM/P', 'IC'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094836-05/04.Mission/FIB/"
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["FIB"]), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

cell_markers = {
    "FIB":["LAMA2","GGT5","LUM","AEBP1","C1S","SFRP1","MEG3","CXCL12"],
    "M-FIB":["SYT1","TNC","PLCXD3","GABRG3","GREB1L","KCNK2","ADAMTSL1","CLMP","NCAM1","FREM1","IGFBP5","MGP","BGN","IGFBP2"],
    "MYOF":["SYNPO2","PCDH7","KCNMA1","LMOD1","TTLL7","DTNA","COL14A1","COL16A1","SULF1","COL6A3","NTM","GLI2","COL5A1","FN1","SLC24A3","CDH13","SMOC2","ANO3","RXFP1"]
}
markers=["LAMA2","GGT5","LUM","AEBP1","C1S","SFRP1","MEG3","CXCL12","SYT1","TNC","PLCXD3","GABRG3","GREB1L","KCNK2","ADAMTSL1","CLMP","NCAM1","FREM1","IGFBP5","MGP","BGN","IGFBP2","SYNPO2","PCDH7","KCNMA1","LMOD1","TTLL7","DTNA","COL14A1","COL16A1","SULF1","COL6A3","NTM","GLI2","COL5A1","FN1","SLC24A3","CDH13","SMOC2","ANO3","RXFP1"]
sc.pl.umap(adata, color=markers, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

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
    '0':'M-FIB',
    '1':'M-FIB',
    '2':'M-FIB',
    '3':'MYOF',
    '4':'FIB'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_loc='on data',legend_fontsize=8)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'FIB.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

cross_tab = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv(os.path.join(out_dir, "celltype_proportions_by_group.csv"))

sc.pl.dotplot(adata, var_names=cell_markers, groupby='cell_type', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), bbox_inches="tight")


# 对Myeloid cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094836-05/04.Mission/Main/Main.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'EC', 'FIB', 'IMM', 'PC', 'PT', 'TAL', 'TL', 'VSM/P', 'IC'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094836-05/04.Mission/Myeloid/"
adata_sub = adata_raw[adata_raw.obs['sub_cell_type'].isin(["Myeloid cells"]), :]
adata_recluster = recluster(adata=adata_sub,species='human',neighbors=20,resolution=0.5,n_pcs = 20,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

cell_markers = {
    "MAC-M1":["NOS2","CD86","CCL5","IL1B","CD274"],
    "MAC-M2":["F13A1","MRC1","CD163","STAB1","SLC1A3","CD14","FOLR2"],
    "MDC":["MSR1","ITGAX","HLA-DQA1","HLA-DRB1","CSF2RA","CD14","TRPM2"],
    "cDC":["ITGAX","HLA-DQA1","HLA-DRA","CSF2RA","CIITA","WDFY4","FLT3","ZNF366","CADM1","ZBTB46","CLEC9A"],
    "pDC":["IRF8","CUX2","P2RY14","IL3RA","CLEC4C"],
    "ncMON":["CTSS","IRAK3","TCF7L2","TNFRSF1B","FCN1","HLA-DRA","FCGR3A"]
}
markers=["NOS2","CD86","CCL5","IL1B","CD274","F13A1","MRC1","CD163","STAB1","SLC1A3","CD14","FOLR2","MSR1","ITGAX","HLA-DQA1","HLA-DRB1","CSF2RA","CD14","TRPM2","ITGAX","HLA-DQA1","HLA-DRA","CSF2RA","CIITA","WDFY4","FLT3","ZNF366","CADM1","ZBTB46","CLEC9A","IRF8","CUX2","P2RY14","IL3RA","CLEC4C","CTSS","IRAK3","TCF7L2","TNFRSF1B","FCN1","HLA-DRA","FCGR3A"]
sc.pl.umap(adata, color=markers, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'), bbox_inches="tight")

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
    '0':'ncMON',
    '1':'MAC-M1',
    '2':'MAC-M2',
    '3':'cDC',
    '4':'MAC-M2',
    '5':'MAC-M1'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_loc='on data',legend_fontsize=8)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))

groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['LARC', 'N']
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
plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})
adata.write(os.path.join(out_dir, 'Myeloid.h5ad'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

cross_tab = pd.crosstab(adata_raw.obs['group'], adata_raw.obs['cell_type'])
# 计算每个分组内的占比（行归一化）
proportion_df = cross_tab.div(cross_tab.sum(axis=1), axis=0)
proportion_df.to_csv(os.path.join(out_dir, "celltype_proportions_by_group.csv"))

cell_type_markers = {
    "MAC-M1":["NOS2","CD86","CCL5","IL1B","CD274"],
    "MAC-M2":["F13A1","MRC1","CD163","STAB1","SLC1A3","CD14","FOLR2"],
    "cDC":["ITGAX","HLA-DQA1","HLA-DRA","CSF2RA","CIITA","WDFY4","FLT3","ZNF366","CADM1","ZBTB46","CLEC9A"],
    "ncMON":["CTSS","IRAK3","TCF7L2","TNFRSF1B","FCN1","HLA-DRA","FCGR3A"]
}
sc.pl.dotplot(adata, var_names=cell_type_markers, groupby='cell_type', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), bbox_inches="tight")


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/wzb/backup/scRNA/SC2509094836_2/2026.03.19/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094836-05/03.Mission/mapping_raw/"
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509094836-05/03.Mission/Epithelial_Basal/Epithelial_Basal.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509094836-05/03.Mission/Principal/PC_IC.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['T cells']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')