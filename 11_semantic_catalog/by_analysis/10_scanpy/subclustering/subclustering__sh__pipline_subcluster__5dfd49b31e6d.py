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

# 对Endothelial注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2604012626-02/01.Mission/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'B Cell', 'DC', 'Endothelial', 'Epithelial', 'Fibroblast', 'Macrophage', 'Mast Cell', 'Neutrophil', 'Plasma Cell', 'T Cell'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2604012626-02/01.Mission/Endothelial/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(["Endothelial"]), :]
adata_sub.obs.rename(columns={'Sample': 'sample'}, inplace=True)
adata_sub.obs.rename(columns={'Group': 'group'}, inplace=True)
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

sc.pl.umap(adata, color=['Sema3c','Kit','Gpihbp1','Ccl21a', 'Prox1','Hey1', 'Igfbp3','Cd36', 'Car4','Ackr1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['Sema3c','Kit','Gpihbp1','Ccl21a', 'Prox1','Hey1', 'Igfbp3','Cd36', 'Car4','Ackr1'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

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
    '0':'arteries',
    '1':'veins',
    '2':'arteries',
    '3':'arteries',
    '4':'LECs',
    '5':'arteries',
    '6':'capillaries',
    '7':'arteries',
    '8':'capillaries'
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
adata.write(os.path.join(out_dir, 'Endothelial.h5ad'))


#### 评分/基因表达
out_dir="/mnt/NexData_file/home/sh/Projects/SC2604012626-02/01.Mission/Endothelial_plot/"
sc.pl.umap(adata, color=['Cd36'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_Cd36.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_Cd36.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
with rc_context({'figure.figsize': (5, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.violin(adata, keys=['Cd36'], groupby='cell_type', show=False)
    plt.savefig(os.path.join(out_dir, 'violin_Cd36.png'))
    plt.savefig(os.path.join(out_dir, 'violin_Cd36.pdf'))

sc.pl.dotplot(adata, var_names=['Cxcl13','Ccl19','Ccl21a','Madcam1','Glycam1','Icam1','Vcam1'], groupby='group', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_HEV-TLS.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_HEV-TLS.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata, var_names=['Hif1a','Pparg','Srebf1'], groupby='group', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_keygene.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_keygene.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
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


# 对B cell,Plasma Cell注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2604012626-02/01.Mission/result.h5ad")
adata_raw.obs.cell_type.cat.categories
# 'B Cell', 'DC', 'Endothelial', 'Epithelial', 'Fibroblast', 'Macrophage', 'Mast Cell', 'Neutrophil', 'Plasma Cell', 'T Cell'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2604012626-02/01.Mission/B_Plasma/"
adata_sub = adata_raw[adata_raw.obs['cell_type'].isin(['B Cell','Plasma Cell']), :]
adata_sub.obs.rename(columns={'Sample': 'sample'}, inplace=True)
adata_sub.obs.rename(columns={'Group': 'group'}, inplace=True)
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

Bcell_marker = {
    "GCB": ["Mme", "Aicda","Sugct","Mki67"],
    "Activate B": ["Il2ra", "Tnfrsf8","Cd69"],
    "Naive B": ["Fcer2a", "Blk","Cxcr5"],
    "Memory B": ["Cd27", "Aim2","Tnfrsf13b","Cd80","Pdcd1lg2","Cxcr3"],
    "Plasma": ["Tnfrsf17","Jchain","Mzb1"]
}
marker=["Mme", "Aicda","Sugct","Mki67","Il2ra", "Tnfrsf8","Cd69","Fcer2a", "Blk", "Cxcr5","Cd27", "Aim2","Tnfrsf13b","Cd80","Pdcd1lg2","Cxcr3","Tnfrsf17","Jchain","Mzb1"]
sc.pl.umap(adata, color=marker, cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_marker.png'))
plt.savefig(os.path.join(out_dir, 'umap_marker.pdf'))
sc.pl.dotplot(adata, var_names=Bcell_marker, groupby='leiden', dendrogram=False, show=False,standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_marker.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_marker.pdf'), bbox_inches="tight")

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
    '1':'Plasma cell',
    '2':'Plasma cell',
    '3':'Activate B',
    '4':'Plasma cell',
    '5':'Memory B'
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
adata.write(os.path.join(out_dir, 'B_Plasma.h5ad'))

adata_sub=adata[adata.obs['cell_type'].isin(['Activate B', 'Memory B', 'Naive B']), :]
adata_sub.write(os.path.join(out_dir, 'Bcells.h5ad'))
#### 评分/基因表达
out_dir="/mnt/NexData_file/home/sh/Projects/SC2604012626-02/01.Mission/B_Plasma_plot/"

sc.pl.dotplot(adata, var_names=['Bcl6','Prdm1','Aicda','Cd80','Cd86','Ciita'], groupby='group', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_Differentiation-Activation.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_Differentiation-Activation.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata, var_names=['Bax','Ddit3','Atf4','Hspa1a','Hspa1b','Hspa2'], groupby='group', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_Apoptosis-Stress.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_Apoptosis-Stress.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata, var_names=['S1pr1','S1pr2','S1pr3'], groupby='group', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_S1P.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_S1P.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata, var_names=['Bcl6','Prdm1','Irf4','Xbp1'], groupby=['group'], dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_keygene.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dotplot_keygene.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
# mito_genes = adata.var_names[adata.var['mt']].tolist() if 'mt' in adata.var.columns else []
genesets = {
    "GOBP_REGULATION_OF_FATTY_ACID_OXIDATION": ["Abcb11","Abcd1","Abcd2","Acacb","Acadl","Acadvl","Acsl5","Akt1","Akt2","Appl2","C1qtnf2","Cnr1","Cpt1a","Dbi","Dgat1","Dgat2","Etfbkmt","Fabp1","Fabp3","Fmo1","Fmo2","Fmo4","Gdf15","Irs1","Irs2","Klhl25","Lonp2","Mfsd2a","Mir199a-2","Mir214","Mlycd","Mtln","Mtor","Nucb2","Obp2a","Pdk4","Plin5","Ppara","Ppard","Pparg","Ppargc1a","Sirt4","Sox9","Twist1","Tysnd1"],
    "Glycolysis_Gluconeogenesis": ["Hk1","Hk2","Hk3","Hkdc1","Gck","Gpi1","Pfkl","Pfkm","Pfkp","Fbp2","Fbp1","Aldoa","Aldoc","Aldob","Aldoart1","Aldoart2","Tpi1","Gm29667","Gapdh","Gapdhs","Pgk1","Pgk2","Pgam1","Pgam2","Eno1","Eno2","Eno3","Eno1b","Eno4","Pkm","Pklr","Pdha1","Pdha2","Pdhb","Dlat","Dld","Pdhx","Ldhal6b","Ldha","Ldhb","Ldhc","Adh1","Adh7","Adh4","Adh5","Akr1a1","Aldh2","Aldh3a2","Aldh1b1","Aldh7a1","Aldh9a1","Aldh3a1","Aldh3b2","Aldh3b1","Aldh3b3","Acss2","Acss1","Galm","Pgm1","Pgm2","G6pc1","G6pc2","G6pc3","Adpgk","Bpgm","Minpp1","Pck1","Pck2"],
    "mito_module": ['mt-Nd2', 'mt-Nd3', 'mt-Nd6', 'mt-Cytb']
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


# 映射所有细分 后面没映射
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/02.Mission/mapping_raw/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/mapping_raw/"
# adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/Monocytes/Macrophage.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD4T/CD4T.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/AE2511046228-02/06.Mission/CD8T/CD8T.h5ad")
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