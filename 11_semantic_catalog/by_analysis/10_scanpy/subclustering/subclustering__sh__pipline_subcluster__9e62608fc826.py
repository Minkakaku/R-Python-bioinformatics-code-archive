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

# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2510235895/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/07.Mission/"
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Fibroblasts/Fibroblasts.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Macrophages/Macrophages.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Neutrophils/Neutrophils.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2510235895/03.Mission/mapping_T/Tcells.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

# adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['lowQ']), :]

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')


#### 基因集评分
adata=sc.read("/mnt/NexData_file/home/sh/Projects/SC2510235895/07.Mission/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/07.Mission/LACTATE_geneset/"
GOBP_LACTATE_BIOSYNTHETIC_PROCESS=['Htt','Ldha','Ldhb','Ldhc','Park7','Per2']
GOBP_LACTATE_METABOLIC_PROCESS=['Actn3','Gatd1','Glo1','Hagh','Haghl','Hif1a','Htt','Ldha','Ldhal6b','Ldhb','Ldhc','Ldhd','Mrs2','Mtch2','Park7','Per2','Pfkfb2','Pnkd','Slc25a12','Slc37a4','Tigar','Trp53']
GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT=['Acacb','Emb','Myc','Slc16a1','Slc16a3','Slc16a7','Slc5a8']
sc.tl.score_genes(adata, gene_list=GOBP_LACTATE_BIOSYNTHETIC_PROCESS, score_name='GOBP_LACTATE_BIOSYNTHETIC_PROCESS')
sc.pl.umap(adata, color='GOBP_LACTATE_BIOSYNTHETIC_PROCESS', size=25, show=False,cmap='RdBu_r')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_GOBP_LACTATE_BIOSYNTHETIC_PROCESS', fig_size={'w': 7, 'h': 6})
# sc.pl.violin(adata, keys='GOBP_LACTATE_BIOSYNTHETIC_PROCESS', groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
# save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_GOBP_LACTATE_BIOSYNTHETIC_PROCESS', fig_size={'w': 8, 'h': 6})
sc.pl.dotplot(adata, var_names=['GOBP_LACTATE_BIOSYNTHETIC_PROCESS'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_GOBP_LACTATE_BIOSYNTHETIC_PROCESS_cell_type', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['GOBP_LACTATE_BIOSYNTHETIC_PROCESS'], groupby='group', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_GOBP_LACTATE_BIOSYNTHETIC_PROCESS_group', fig_size={'w': 6, 'h': 6})
# sc.pl.dotplot(adata, var_names=['GOBP_LACTATE_BIOSYNTHETIC_PROCESS'], groupby=['cell_type','group'], standard_scale='var',  cmap='RdBu_r',  show=False)
# save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_GOBP_LACTATE_BIOSYNTHETIC_PROCESS_cell_type-group', fig_size={'w': 6, 'h': 6})
sc.tl.score_genes(adata, gene_list=GOBP_LACTATE_METABOLIC_PROCESS, score_name='GOBP_LACTATE_METABOLIC_PROCESS')
sc.pl.umap(adata, color='GOBP_LACTATE_METABOLIC_PROCESS', size=25, show=False,cmap='RdBu_r')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_GOBP_LACTATE_METABOLIC_PROCESS', fig_size={'w': 7, 'h': 6})

sc.pl.dotplot(adata, var_names=['GOBP_LACTATE_METABOLIC_PROCESS'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_GOBP_LACTATE_METABOLIC_PROCESS_cell_type', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['GOBP_LACTATE_METABOLIC_PROCESS'], groupby='group', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_GOBP_LACTATE_METABOLIC_PROCESS_group', fig_size={'w': 6, 'h': 6})

sc.tl.score_genes(adata, gene_list=GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT, score_name='GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT')
sc.pl.umap(adata, color='GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT', size=25, show=False,cmap='RdBu_r')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT', fig_size={'w': 7, 'h': 6})

sc.pl.dotplot(adata, var_names=['GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT_cell_type', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata, var_names=['GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT'], groupby='group', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT_group', fig_size={'w': 6, 'h': 6})

adata_sub=sc.read("/mnt/NexData_file/home/sh/Projects/SC2510235895/03.Mission/mapping_T/Tcells.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/07.Mission/LACTATE_geneset/Tsub/"
sc.tl.score_genes(adata_sub, gene_list=GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT, score_name='GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT')
sc.pl.umap(adata_sub, color='GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT', size=25, show=False,cmap='RdBu_r')
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT', fig_size={'w': 7, 'h': 6})
sc.pl.dotplot(adata_sub, var_names=['GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT'], groupby='cell_type', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT_cell_type', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata_sub, var_names=['GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT'], groupby='group', standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT_group', fig_size={'w': 6, 'h': 6})
sc.pl.dotplot(adata_sub, var_names=['GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT'], groupby=['cell_type','group'], standard_scale='var',  cmap='RdBu_r',  show=False)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT_cell_type-group', fig_size={'w': 6, 'h': 6})