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

# 对T Cell注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2604012626-02/01.Mission/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'B Cell', 'DC', 'Endothelial', 'Epithelial', 'Fibroblast', 'Macrophage', 'Mast Cell', 'Neutrophil', 'Plasma Cell', 'T Cell'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2604012626-02/02.Mission/TCell/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["T Cell"]), :]
adata_sub.obs.rename(columns={'Sample': 'sample'}, inplace=True)
adata_sub.obs.rename(columns={'Group': 'group'}, inplace=True)
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.8,n_pcs = 15,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['Cd4','Cd8a','Cd8b1','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra','Trdc',"Trgc1", "Trgc2",'Mki67','Il17a','Rorc','Il23r',"Sell","Lef1","Cd69","Cxcr5","Bcl6","Batf","Tox","Pdcd1"], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'),bbox_inches='tight', dpi=300)
sc.pl.dotplot(adata, var_names=['Cd4','Cd8a','Cd8b1','Foxp3','Gzma','Nkg7','Cd3d','Cd3e','Ncam1','Fcgr3','Ctla4','Il2ra','Trdc',"Trgc1", "Trgc2",'Mki67','Il17a','Rorc','Il23r',"Sell","Lef1","Cd69","Cxcr5","Bcl6","Batf","Tox","Pdcd1"], groupby='leiden', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight', dpi=300)

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

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
    '0':'Treg',
    '1':'CD4T',
    '2':'GDT',
    '3':'Treg',
    '4':'NKT',
    '5':'GDT',
    '6':'Treg',
    '7':'CD8T',
    '8':'CD4T',
    '9':'CD4T',
    '10':'Treg',
    '11':'Proliferative cells',
    '12':'CD4T'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

groups = adata.obs['group'].drop_duplicates().tolist()
# groups = ['NC', 'HLH', 'iMCD-TAFRO']
# adata.obs['group'] = (
#     pd.Categorical(adata.obs['group'],
#                    categories=groups,
#                    ordered=True)
# )
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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':5, 'h':5})
adata.write(os.path.join(out_dir, 'TCell.h5ad'))


#### 评分/基因表达
out_dir="/mnt/NexData_file/home/sh/Projects/SC2604012626-02/02.Mission/TCell_plot/"
sc.pl.dotplot(adata, var_names=['Cxcl13','Ccl19','Ccl21a'], groupby='group', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_TLS.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_TLS.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata, var_names=['Cxcl13','Ccl19','Ccl21a'], groupby=['cell_type','group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_TLS_sub.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_TLS_sub.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

genesets = {
    "GOBP_REGULATION_OF_FATTY_ACID_OXIDATION": ["Abcb11","Abcd1","Abcd2","Acacb","Acadl","Acadvl","Acsl5","Akt1","Akt2","Appl2","C1qtnf2","Cnr1","Cpt1a","Dbi","Dgat1","Dgat2","Etfbkmt","Fabp1","Fabp3","Fmo1","Fmo2","Fmo4","Gdf15","Irs1","Irs2","Klhl25","Lonp2","Mfsd2a","Mir199a-2","Mir214","Mlycd","Mtln","Mtor","Nucb2","Obp2a","Pdk4","Plin5","Ppara","Ppard","Pparg","Ppargc1a","Sirt4","Sox9","Twist1","Tysnd1"],
    "Glycolysis_Gluconeogenesis": ["Hk1","Hk2","Hk3","Hkdc1","Gck","Gpi1","Pfkl","Pfkm","Pfkp","Fbp2","Fbp1","Aldoa","Aldoc","Aldob","Aldoart1","Aldoart2","Tpi1","Gm29667","Gapdh","Gapdhs","Pgk1","Pgk2","Pgam1","Pgam2","Eno1","Eno2","Eno3","Eno1b","Eno4","Pkm","Pklr","Pdha1","Pdha2","Pdhb","Dlat","Dld","Pdhx","Ldhal6b","Ldha","Ldhb","Ldhc","Adh1","Adh7","Adh4","Adh5","Akr1a1","Aldh2","Aldh3a2","Aldh1b1","Aldh7a1","Aldh9a1","Aldh3a1","Aldh3b2","Aldh3b1","Aldh3b3","Acss2","Acss1","Galm","Pgm1","Pgm2","G6pc1","G6pc2","G6pc3","Adpgk","Bpgm","Minpp1","Pck1","Pck2"],
    "Sphingolipid_metabolism": ["Sptlc2","Sptlc3","Sptlc1","Kdsr","Cers1","Cers4","Cers2","Cers3","Cers6","Cers5","Asah1","Asah2","Acer1","Acer2","Acer3","Degs1","Degs1l","Degs2","Sgms1","Sgms2","Smpd1","Smpd2","Smpd3","Smpd4","Enpp7","Cerk","Plpp1","Plpp2","Plpp3","Sgpp1","Sgpp2","Samd8","Sphk1","Sphk2","Sgpl1","Ugcg","Gba1","Gba2","B4galt6","B4galt5","Glb1","Ugt8a","Galc","Gal3st1","Arsa","Neu1","Neu2","Neu4","Neu3","Gla","B4galnt1","Hexa","Hexb","Gm2a","Psap","Psapl1"]
}
for set_name, genes in genesets.items():
    sc.tl.score_genes(adata, gene_list=genes, score_name=set_name)

for set_name, genes in genesets.items():
    sc.pl.umap(adata, color=set_name, size=25, show=False,cmap='RdBu_r')
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_{set_name}_score', fig_size={'w': 7, 'h': 6})
    sc.pl.violin(adata, keys=set_name, groupby='cell_type', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_{set_name}_score_cell_type', fig_size={'w': 10, 'h': 6})
    sc.pl.violin(adata, keys=set_name, groupby='group', rotation=45, show=False, cmap='RdBu_r')          # 颜色同步倒转
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'violin_{set_name}_score_group', fig_size={'w': 10, 'h': 6})
    sc.pl.dotplot(adata, var_names=[set_name], groupby='group', standard_scale='var',  cmap='RdBu_r',  show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_score', fig_size={'w': 6, 'h': 6})


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2604012626-02/01.Mission/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2604012626-02/02.Mission/mapping_raw/"
adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2604012626-02/01.Mission/B_Plasma/B_Plasma.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2604012626-02/01.Mission/Endothelial/Endothelial.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

# adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['lowQ']), :]

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
    plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white')
    plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white')

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')