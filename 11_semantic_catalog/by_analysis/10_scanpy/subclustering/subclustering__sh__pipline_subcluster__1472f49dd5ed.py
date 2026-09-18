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

# 调整Main注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/01.Mission/Main/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Proximal tubule cells', 'Macrophages', 'Endothelial cells', 'Fibroblasts', 'Distal tubule cells', 'Podocytes', 'T cells', 'Intercalated cells', 'Principal cells', 'Urothelial cells'
adata_raw.obs['sample'].cat.categories
adata_raw = adata_raw[adata_raw.obs['sample'].isin(["CF8","K8","JS3"]), :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/03.Mission/Main/"
adata_recluster = recluster(adata=adata_raw,species='mouse',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

marker_dict = {
    'Proximal tubule cells': ['Slc34a1', 'Lrp2', 'Cubn'],  # 近端肾小管
    'Macrophages': ['Cd74', 'Ctss', 'Mrc1'],  # 巨噬细胞
    'Endothelial cells': ['Flt1','Emcn', 'Ldb2'],  # 内皮细胞
    'Fibroblasts': ['Cald1', 'Mylk', 'Pdgfra'],  # 成纤维细胞
    'Loop of henle': ['Slc12a1', 'Umod', 'Egf'],  # 亨勒环
    'Distal tubule cells': ['Slc12a3', 'Trpm6', 'Klhl3'],  # 远端小管
    'Podocytes': ['Nphs1', 'Ptpro', 'Nebl'],  # 足细胞
    'Intercalated cells': ['Atp6v0d2', 'Atp6v1g3', 'Slc26a4'],  # 集合管插入细胞
    'Principal cells': ['Scnn1b', 'Fxyd4', 'Phactr1'],  # 主细胞
    'T cells': ['Cd3d', 'Cd3e', 'Ms4a4b',],  # T细胞
    'Urothelial cells': ['Upk1b', 'Tmprss2', 'Upk3a'],  # 尿路上皮细胞
    'Adipocytes': ['Pparg', 'Plin1', 'Adipor2']  # 脂肪细胞
}

all_markers = []
for cell_type, markers in marker_dict.items():
    all_markers.extend(markers)

sc.pl.umap(adata, color=all_markers, cmap='coolwarm', ncols=3)
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=marker_dict, groupby='leiden', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'Proximal tubule cells',
    '1':'Proximal tubule cells',
    '2':'Proximal tubule cells',
    '3':'Endothelial cells',
    '4':'Proximal tubule cells',
    '5':'LOH',
    '6':'Principal cells',
    '7':'Proximal tubule cells',
    '8':'Distal tubule cells',
    '9':'Proximal tubule cells',
    '10':'Fibroblasts',
    '11':'Macrophages',
    '12':'Urothelial cells',
    '13':'Proximal tubule cells',
    '14':'Proximal tubule cells'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sample', title='sample', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='group', title='group', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.pdf'))

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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':5, 'h':5})

# method = "cell_type"
# keys = adata_raw.obs[method].cat.categories.to_list()
# values = list(adata_raw.uns[method+'_colors'])
# color_dict = dict(zip(keys,values))
# df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# # df.columns.names = ["cell_type"]
# df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

# sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
# sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True)
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')
# find_group_marker(adata=adata_raw,species='mouse',output=out_dir)
adata.write(os.path.join(out_dir, 'Main.h5ad'))
marker_dict = {
    'Distal tubule cells': ['Slc12a3', 'Trpm6', 'Klhl3'],'Endothelial cells': ['Flt1','Emcn', 'Ldb2'],'Fibroblasts': ['Cald1', 'Mylk', 'Pdgfra'],'Loop of henle': ['Slc12a1', 'Umod', 'Egf'],'Macrophages': ['Cd74', 'Ctss', 'Mrc1'],'Principal cells': ['Scnn1b', 'Fxyd4', 'Phactr1'],
    'Proximal tubule cells': ['Slc34a1', 'Lrp2', 'Cubn'],'Urothelial cells': ['Upk1b', 'Tmprss2', 'Upk3a']
}
sc.pl.stacked_violin(
    adata, 
    var_names=marker_dict,
    groupby='cell_type',  # 替换为您的细胞类型列名
    # palette=palette,  # 指定颜色映射
    figsize=(18, 8),
    var_group_positions=None,
    var_group_labels=None,
    dendrogram=False,
    show=False
)
# 保存图片
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'marker_violin_stacked.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'marker_violin_stacked.png'), dpi=300, bbox_inches='tight')
plt.show()
sc.pl.dotplot(adata, var_names=marker_dict, groupby=['cell_type'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'marker_dotplot_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'marker_dotplot_cell_type.pdf'),bbox_inches='tight')


gene_sets = {
    "TLR4_recognize": ["Tlr4", "Cd14", "Ly96", "Myd88"],
    "NFKB_pathway": ["Rela", "p65", "Nfkb1", "Ikkb", "Nfkbia"],
    "inflammasome": ["Nlrp3", "Pycard", "Asc", "Casp1"],
    "inflammation_factors": ["Il1b", "Il18", "Tnf", "Il6", "Ccl2", "Cxcl10"],
    "pyroptosis": ["Gsdmd", "Casp4", "Casp11"]
}
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/03.Mission/Main_plot/"
adata=adata.raw.to_adata()
sc.tl.score_genes(adata, gene_list=gene_sets["TLR4_recognize"], score_name='TLR4_recognize')
sc.tl.score_genes(adata, gene_list=gene_sets["NFKB_pathway"], score_name='NFKB_pathway')
sc.tl.score_genes(adata, gene_list=gene_sets["inflammasome"], score_name='inflammasome')
sc.tl.score_genes(adata, gene_list=gene_sets["inflammation_factors"], score_name='inflammation_factors')
sc.tl.score_genes(adata, gene_list=gene_sets["pyroptosis"], score_name='pyroptosis')
sc.pl.umap(adata, color=['TLR4_recognize','NFKB_pathway','inflammasome','inflammation_factors','pyroptosis'], size=25, show=False,cmap='RdBu_r',vmax='p99',ncols=5)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_genesets', fig_size={'w': 30, 'h': 6})

sc.pl.dotplot(adata, var_names=['TLR4_recognize','NFKB_pathway','inflammasome','inflammation_factors','pyroptosis'], groupby=['cell_type'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['TLR4_recognize','NFKB_pathway','inflammasome','inflammation_factors','pyroptosis'], groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type_group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['TLR4_recognize','NFKB_pathway','inflammasome','inflammation_factors','pyroptosis'], groupby=['group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_group.pdf'),bbox_inches='tight')

all_genes = [
    "Tlr4", "Cd14", "Ly96", "Myd88",
    "Rela", "Nfkb1", "Nfkbia",
    "Nlrp3", "Pycard", "Casp1",
    "Il1b", "Il18", "Tnf", "Il6", "Ccl2", "Cxcl10",
    "Gsdmd", "Casp4"]
sc.pl.umap(adata, color=all_genes, size=25, show=False,cmap='RdBu_r',vmax='p99',ncols=4)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_genes', fig_size={'w': 24, 'h': 30})

sc.pl.dotplot(adata, var_names=all_genes, groupby=['cell_type'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=all_genes, groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type_group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=all_genes, groupby=['group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.pdf'),bbox_inches='tight')


# 调整Proximal tubule cells注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2512227671-01/03.Mission/Main/Main.h5ad")
adata_raw = adata_raw[adata_raw.obs['cell_type'].isin(["Proximal tubule cells"]), :]
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/03.Mission/PTC/"
adata_recluster = recluster(adata=adata_raw,species='mouse',neighbors=10,resolution=0.8,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

adata_sub = adata[~adata.obs['leiden'].isin(["14"]), :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
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

sc.pl.umap(adata, color=['Slc34a1','Slc5a2','Lrp2','Havcr1', 'Lcn2','Spp1','Ccl2','Cxcl10','Vcam1','Nfkbia','Hspa1a','Ddit3','Atf3','Fn1','Col1a1','Tgfb1'], cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['Slc34a1','Slc5a2','Lrp2','Havcr1', 'Lcn2','Spp1','Ccl2','Cxcl10','Vcam1','Nfkbia','Hspa1a','Ddit3','Atf3','Fn1','Col1a1','Tgfb1'], groupby='leiden', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'PT_Slc34a1',
    '1':'PT_Cmss1',
    '2':'PT_Kap',
    '3':'PT_Slc5a12',
    '4':'PT_Slc34a1',
    '5':'PT_Ghr',
    '6':'PT_Slc13a1',
    '7':'PT_Slc34a1',
    '8':'PT_Slc34a1',
    '9':'PT_Slc5a12'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'PTC.h5ad'))

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
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':5, 'h':5})

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

gene_sets = {
    "TLR4_recognize": ["Tlr4", "Cd14", "Ly96", "Myd88"],
    "NFKB_pathway": ["Rela", "p65", "Nfkb1", "Ikkb", "Nfkbia"],
    "inflammasome": ["Nlrp3", "Pycard", "Asc", "Casp1"]
}
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/03.Mission/PTC_plot/"
adata=adata.raw.to_adata()
sc.tl.score_genes(adata, gene_list=gene_sets["TLR4_recognize"], score_name='TLR4_recognize')
sc.tl.score_genes(adata, gene_list=gene_sets["NFKB_pathway"], score_name='NFKB_pathway')
sc.tl.score_genes(adata, gene_list=gene_sets["inflammasome"], score_name='inflammasome')
sc.pl.umap(adata, color=['TLR4_recognize','NFKB_pathway','inflammasome'], size=25, show=False,cmap='RdBu_r',vmax='p99',ncols=3)
save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'umap_genesets', fig_size={'w': 18, 'h': 6})

sc.pl.dotplot(adata, var_names=['TLR4_recognize','NFKB_pathway','inflammasome'], groupby=['cell_type'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['TLR4_recognize','NFKB_pathway','inflammasome'], groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_cell_type_group.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['TLR4_recognize','NFKB_pathway','inflammasome'], groupby=['group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genesets_group.pdf'),bbox_inches='tight')

gene_sets = {
    "TLR4_recognize": ["Tlr4", "Cd14", "Ly96", "Myd88"],
    "NFKB_pathway": ["Rela","Nfkb1", "Nfkbia"],
    "inflammasome": ["Nlrp3", "Pycard","Casp1"]
}
sc.pl.violin(adata, ['TLR4_recognize','NFKB_pathway','inflammasome'], jitter=0, groupby='group', rotation=45, multi_panel=True, show=False)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'gene_sets_violin.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'gene_sets_violin.png'), dpi=300, bbox_inches='tight')
sc.pl.stacked_violin(
    adata, 
    var_names=gene_sets,
    groupby='group',  # 替换为您的细胞类型列名
    # palette=palette,  # 指定颜色映射
    figsize=(18, 8),
    var_group_positions=None,
    var_group_labels=None,
    dendrogram=False,
    show=False
)
# 保存图片
plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'gene_sets_violin_stacked.pdf'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'gene_sets_violin_stacked.png'), dpi=300, bbox_inches='tight')
plt.show()