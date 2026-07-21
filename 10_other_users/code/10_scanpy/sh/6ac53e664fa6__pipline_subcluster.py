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
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/AE2512227671-01/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Proximal tubule cells', 'Macrophages', 'Endothelial cells', 'Fibroblasts', 'Distal tubule cells', 'Podocytes', 'T cells', 'Intercalated cells', 'Principal cells', 'Urothelial cells'
out_dir="/mnt/NexData_file/home/sh/Projects/AE2512227671-01/01.Mission/Main/"
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,legend_fontsize=5,legend_loc='on data')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype_label.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sample', title='sample', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='group', title='group', show=False,legend_fontsize=5)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.pdf'))

plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group',fig_size={'w':10, 'h':5})
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':15, 'h':5})

method = "cell_type"
keys = adata_raw.obs[method].cat.categories.to_list()
values = list(adata_raw.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

sc.tl.rank_genes_groups(adata_raw, groupby='cell_type', method='wilcoxon', pts=True)
sc.pl.rank_genes_groups_heatmap(adata_raw, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')
find_group_marker(adata=adata_raw,species='mouse',output=out_dir)
adata_raw.write(os.path.join(out_dir, 'result.h5ad'))
# marker_dict = {
#     'Proximal tubule cells': ['Slc34a1', 'Lrp2', 'Cubn'],  # 近端肾小管
#     'Macrophages': ['Cd74', 'Ctss', 'Mrc1'],  # 巨噬细胞
#     'Endothelial cells': ['Flt1','Emcn', 'Ldb2'],  # 内皮细胞
#     'Fibroblasts': ['Cald1', 'Mylk', 'Pdgfra'],  # 成纤维细胞
#     'Loop of henle': ['Slc12a1', 'Umod', 'Egf'],  # 亨勒环
#     'Distal tubule cells': ['Slc12a3', 'Trpm6', 'Klhl3'],  # 远端小管
#     'Podocytes': ['Nphs1', 'Ptpro', 'Nebl'],  # 足细胞
#     'Intercalated cells': ['Atp6v0d2', 'Atp6v1g3', 'Slc26a4'],  # 集合管插入细胞
#     'Principal cells': ['Scnn1b', 'Fxyd4', 'Phactr1'],  # 主细胞
#     'T cells': ['Cd3d', 'Cd3e', 'Ms4a4b',],  # T细胞
#     'Urothelial cells': ['Upk1b', 'Tmprss2', 'Upk3a'],  # 尿路上皮细胞
#     'Adipocytes': ['Pparg', 'Plin1', 'Adipor2']  # 脂肪细胞
# }

# all_markers = []
# for cell_type, markers in marker_dict.items():
#     all_markers.extend(markers)

# sc.pl.stacked_violin(
#     adata_raw, 
#     var_names=marker_dict,
#     groupby='cell_type',  # 替换为您的细胞类型列名
#     # palette=palette,  # 指定颜色映射
#     figsize=(18, 8),
#     var_group_positions=None,
#     var_group_labels=None,
#     dendrogram=False,
#     show=False
# )
# # 保存图片
# plt.tight_layout()
# plt.savefig(os.path.join(out_dir, 'marker_violin_stacked.pdf'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'marker_violin_stacked.png'), dpi=300, bbox_inches='tight')
# plt.show()


# 调整Proximal tubule cells注释    #### 后面都没做
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2512227671/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Proximal tubule cells', 'Macrophages', 'Endothelial cells', 'Fibroblasts', 'Loop of henle', 'Distal tubule cells', 'Podocytes', 'Intercalated cells', 'Principal cells', 'T cells', 'Urothelial cells', 'Adipocytes'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512227671-05/01.Mission/PTC/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Proximal tubule cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
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

# sc.pl.umap(adata, color=['Lyz2','Apoe','Cd14','Fcgr3','Vcan','Mki67','Top2a','Cd86','Mrc1','Clec9a','Cd209a','Siglech','Tcf4'], cmap='coolwarm', ncols=3)
# plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
# plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
# sc.pl.dotplot(adata, var_names=['Lyz2','Apoe','Cd14','Fcgr3','Vcan','Mki67','Top2a','Cd86','Mrc1','Clec9a','Cd209a','Siglech','Tcf4'], groupby='leiden', dendrogram=False, show=False, standard_scale='var')
# plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'PT1_Slc5a12',
    '1':'PT2_Cmss1',
    '2':'PT3_Gpc6',
    '3':'PT4_Kap',
    '4':'PT5_Anks1b',
    '5':'PT6_Gramd1b',
    '6':'PT4_Kap',
    '7':'PT5_Anks1b',
    '8':'PT7_Atp11a'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

# sc.pl.rank_genes_groups_dotplot(
#     adata,
#     groupby='cell_type',
#     key='rank_genes_groups',
#     n_genes=3,
#     show=False)
# plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
# plt.savefig(os.path.join(out_dir, 'top3_diff_genes_dotplot.pdf'),bbox_inches='tight')
# plt.close()

marker_dict = {
    'PT1_Slc5a12': ['Slc5a12','Itpr2','Gatm'],
    'PT2_Cmss1': ['Cmss1','Camk1d','Gm19951'], 
    'PT3_Gpc6': ['Gpc6','Nat8f6','Hsd3b3'], 
    'PT4_Kap': ['Kap','Gpx3','Fth1'],
    'PT5_Anks1b': ['Anks1b','Esrrg','Erc2'], 
    'PT6_Gramd1b': ['Gramd1b','Ghr','Cyp7b1'],
    'PT7_Atp11a': ['Atp11a','Col4a3','Gm42397']
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

plot_celltype_proportion(adata, out_dir=out_dir,type='group',fig_size={'w':15, 'h':5})
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':15, 'h':5})

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

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

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})

sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var', swap_axes=True,show_gene_labels=True)
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

sc.pl.umap(adata, color=['Srebf1','Acaca','Fasn','Scd1','Elovl6'], cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['Srebf1','Acaca','Fasn','Scd1','Elovl6'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Srebf1','Acaca','Fasn','Scd1','Elovl6'], groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_group.pdf'),bbox_inches='tight')


# 映射所有细分 后面没细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2511266943-01/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2511266943-05/06.Mission/Main_mapping/"

adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2511266943-05/06.Mission/Macrophages/Macrophages.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2511266943-05/06.Mission/NK/NK.h5ad")
bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2511266943-05/06.Mission/Tcells/Tcells.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

adata_raw = adata_raw[~adata_raw.obs['sub_cell_type'].isin(['B cells']), :]
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