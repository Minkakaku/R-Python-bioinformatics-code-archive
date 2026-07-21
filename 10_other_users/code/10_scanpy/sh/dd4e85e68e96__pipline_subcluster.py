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
import seaborn as sns

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

# 调整Main注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/12.Mission/Main/Main.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Distal tubule cells', 'Endothelial cells', 'Fibroblasts', 'IMM', 'Intercalated cells', 'LOH', 'Principal cells', 'Proximal tubule cells', 'Urothelial cells'
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/13.Mission/"
sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

# 调整Fibroblasts注释 
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/12.Mission/Main/Main.h5ad")
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["Fibroblasts"]), :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/13.Mission/Fibroblasts/"
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
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

# 如果 FIB1 和 FIB2 要合并为同一个细胞类型
cell_markers = {
    "VSM/P": ["Notch3", "Pdgfrb", "Itga8"],
    "MC": ["Pip5k1b", "Robo1", "Piezo2", "Daam2", "Phtf2", "Gata3", "Postn"],
    "REN": ["Pip5k1b", "Robo1", "Pde10a", "Abcc8", "Col13a1", "Grid2"],
    "VSMC": ["Ntrk3", "Myh11", "Rgs6", "Adra1a", "Ldb3", "Mcam"],
    "VSMC/P": ["Ntrk3", "Rgs5", "Abcc9", "Adcy3", "Adgrb3"],
    "FIB": ["Col1a1", "Col1a2", "C7", "Negr1", "Fbln5", "Dcn", "Cdh11", "Lama2", "Ggt5", "Lum", "Aebp1", "Sfrp1", "Meg3", "Cxcl12"],
    "M-FIB": ["Syt1", "Tnc", "Plcxd3", "Gabrg3", "Kcnk2", "Adamtsl1", "Clmp", "Ncam1", "Frem1", "Igfbp5", "Mgp", "Bgn", "Igfbp2"],
    "MYOF": ["Synpo2", "Pcdh7", "Kcnma1", "Lmod1", "Ttll7", "Dtna", "Col14a1", "Lama2", "Col16a1", "Sulf1", "Col6a3", "Ntm", "Gli2", "Col5a1", "Fn1", "Slc24a3", "Cdh13", "Smoc2", "Ano3", "Rxfp1"]
}
all_genes = []
for markers in cell_markers.values():
    all_genes.extend(markers)
sc.pl.umap(adata, color=all_genes, cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Myofibroblast',
    '1':'Fibroblast',
    '2':'Pericyte',
    '3':'Pericyte',
    '4':'Myofibroblast',
    '5':'VSMC',
    '6':'Fibroblast',
    '7':'Fibroblast',
    '8':'Medullary Fibroblast',
    '9':'Fibroblast',
    '10':'Fibroblast',
    '11':'Myofibroblast',
    '12':'Mesangial cell',
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Fibroblasts.h5ad'))

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

plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='sample', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='group', title='group', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.pdf'))


# 调整IMM注释 
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/12.Mission/Main/Main.h5ad")
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["IMM"]), :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/13.Mission/IMM/"
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
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

cell_markers = {
    "B cell": ["Bank1", "Blk", "Ms4a1", "Bach2"],
    "Plasma cell": ["Tent5c", "Mzb1", "Fcrl5", "Cd38", "Jchain"],
    "T cell": ["Cd3d","Cd3e","Cd4","Cd8a","Foxp3","Cd96", "Cd247", "Themis", "Bcl11b", "Camk4", "Il7r","Nkg7","Ncam1"],
    "MAST": ["Ms4a2", "Cpa3", "Kit"],
    "MAC-M2": ["Mrc1", "Cd163"],
    "MDC": ["Msr1", "Itgax", "H2-Aa", "H2-Ab1", "Csf2ra", "Cd14", "Trpm2"],
    "cDC": ["Itgax", "H2-Aa", "H2-Ab1", "Csf2ra", "Ciita", "Wdfy4", "Flt3", "Zbtb46", "Clec9a"],
    "pDC": ["Irf8", "Cux2", "P2ry14", "Il3ra"],
    "ncMON": ["Ctss", "Irak3", "Tcf7l2", "Tnfrsf1b", "H2-Aa"],
    "Neutrophil": ["S100a9", "S100a8", "Ifitm2"]
}
all_genes = []
for markers in cell_markers.values():
    all_genes.extend(markers)
sc.pl.umap(adata, color=all_genes, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=cell_markers, groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

sc.pl.umap(adata, color='sample', title='sample', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_sample.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_sample.pdf"),bbox_inches='tight')

anno_dict={
    '0':'Macrophage',
    '1':'Macrophage',
    '2':'T cell',
    '3':'T cell',
    '4':'cDC',
    '5':'B cell',
    '6':'Plasma cell',
    '7':'cDC',
    '8':'ncMON',
    '9':'T cell',
    '10':'Neutrophil',
    '11':'pDC',
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'IMM.h5ad'))

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

plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='sample', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='group', title='group', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.pdf'))


# 调整Principal cells注释 
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/12.Mission/Main/Main.h5ad")
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["Principal cells"]), :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/13.Mission/Principal/"
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()


# 调整DTC注释 
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/12.Mission/Main/Main.h5ad")
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["Distal tubule cells"]), :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/13.Mission/DTC/"
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()


# 调整Intercalated注释 
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/12.Mission/Main/Main.h5ad")
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["Intercalated cells"]), :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/13.Mission/Intercalated/"
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/12.Mission/Main/Main.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/13.Mission/Main_mapping/"
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/13.Mission/IMM/IMM.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/12.Mission/PTC/PTC.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/12.Mission/Endothelial/Endothelial.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Sub Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_subcelltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_subcelltype.pdf'))