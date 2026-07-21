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

# 调整Myofibroblasts注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2509094832_ERV/result.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'DC', 'Endothelial cells', 'Fibroblasts', 'Macrophages', 'Muscle cells', 'Pericytes', 'Plasma cells', 'Satellite cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Main_mod/"
sc.pl.umap(adata_raw, color='cell_type', title='Cell Type', show=False)
save_fig(plt, out_dir, file_name="dimplot_umap_celltype", fig_size={'w':6, 'h':5})
sc.pl.umap(adata_raw, color='cluster', title='cluster', show=False,legend_loc='on data')
save_fig(plt, out_dir, file_name="dimplot_umap_cluster", fig_size={'w':5, 'h':5})

adata_raw = adata_raw[~adata_raw.obs['cluster'].isin(["0"]), :]
sc.tl.umap(adata_raw)

anno_dict={
    '1':'Myofibroblasts'
}
mask = adata_raw.obs['cluster'].isin(anno_dict.keys())
# 获取新的类别列表
new_categories = list(adata_raw.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata_raw.obs['cell_type'].cat.categories]
# 更新 categories
adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].cat.set_categories(new_categories)
adata_raw.obs.loc[mask, 'cell_type'] = adata_raw.obs.loc[mask, 'cluster'].map(anno_dict)
sc.pl.umap(adata_raw, color='cell_type', show=False)
save_fig(plt, out_dir, file_name='sub_cluster_umap')
sc.pl.umap(adata_raw, color='cell_type',legend_loc ='on data', show=False)
save_fig(plt, out_dir, file_name='sub_cluster_umap_withlabel')

sc.pl.umap(adata_raw, color=['LTR/ERVK','LTR/ERV1','LTR/ERVL','LTR/ERVL-MaLR','CNN1','ACTA2','TAGLN','C1R','COL1A2','DCN','PDGFRA'], cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.dotplot(adata_raw, var_names=['LTR/ERVK','LTR/ERV1','LTR/ERVL','LTR/ERVL-MaLR','CNN1','ACTA2','TAGLN','C1R','COL1A2','DCN','PDGFRA'], groupby='cluster', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata_raw, var_names=['LTR/ERVK','LTR/ERV1','LTR/ERVL','LTR/ERVL-MaLR','CNN1','ACTA2','TAGLN','C1R','COL1A2','DCN','PDGFRA'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_diffgenes_cell_type.pdf'),bbox_inches='tight')
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))


#### 基因表达展示
dsRNA=['IFIH1','RIGI','EIF2AK2','TLR3']
dsDNA=['CGAS','STING1']
Signal_Transduction_Transcription=['MAVS','TBK1','IKBKE','IRF3','IRF7']
Product_IFN=['IFNB1']
ISG=['ISG15','MX1','MX2','IFIT1','IFIT3','OAS1','RSAD2','USP18']
Interferon_receptor=['IFNAR1','IFNAR2']
pathways = {'dsRNA': dsRNA,'dsDNA': dsDNA,'Signal_Transduction_Transcription': Signal_Transduction_Transcription,'Product_IFN': Product_IFN,'ISG': ISG,'Interferon_receptor': Interferon_receptor}

base_dir = "/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Main_geneplot/"
for name, genes in pathways.items():
    out_dir = os.path.join(base_dir, name)
    os.makedirs(out_dir, exist_ok=True)
    sc.pl.dotplot(adata_raw, var_names=genes, groupby='cell_type', dendrogram=False, show=False)
    plt.savefig(os.path.join(out_dir, 'genelist_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'genelist_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
    sc.pl.umap(adata_raw, color=genes, cmap='coolwarm', ncols=3)
    plt.savefig(os.path.join(out_dir, 'genelist_umap.png'))
    plt.savefig(os.path.join(out_dir, 'genelist_umap.pdf'))
    for g in genes:
        sc.pl.dotplot(adata_raw, var_names=g, groupby=['cell_type','group'], dendrogram=False, show=False)
        plt.savefig(os.path.join(out_dir, f'{g}_dotplot_group_celltype.png'), facecolor='white', bbox_inches='tight', dpi=300)
        plt.savefig(os.path.join(out_dir, f'{g}_dotplot_group_celltype.pdf'), facecolor='white', bbox_inches='tight', dpi=300)
        plt.close()



# 对Macrophages注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Main_mod/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'DC', 'Endothelial cells', 'Fibroblasts', 'Macrophages', 'Muscle cells', 'Pericytes', 'Plasma cells', 'Satellite cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Macrophages/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Macrophages", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.2,n_pcs = 15,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
adata.obs['cell_type']=adata.obs['leiden']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False,legend_loc ='on data',legend_fontsize=5)
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

sc.pl.umap(adata, color=['LYZ','APOE','CD14','FCGR3A','VCAN','MKI67','TOP2A','CD86','MRC1','CLEC9A','CLEC4C'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))

anno_dict={
    '0':'HLA-A Macrophages',
    '1':'F13A1 Macrophages',
    '2':'cDC'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Macrophages.h5ad'))

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

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 对Plasma cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Main_mod/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'DC', 'Endothelial cells', 'Fibroblasts', 'Macrophages', 'Muscle cells', 'Pericytes', 'Plasma cells', 'Satellite cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Plasma/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Plasma cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.2,n_pcs = 10,output=out_dir,batch_effect = 2)
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

sc.pl.umap(adata, color=['IGHA1','IGHG1','IGHM','IGHD','IRF4','PRDM1','CD38'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))

# anno_dict={
#     '0':'HLA-A Macrophages',
#     '1':'F13A1 Macrophages',
#     '2':'cDC'
# }
# adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
# adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

# sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
# plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
# plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

# sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
# find_group_marker(adata=adata,species='human',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Plasma.h5ad'))

# groups = adata.obs['group'].cat.categories.to_list()
# if len(groups) > 1:
#     axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
#     fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
#     for i in range(len(groups)):
#         if i == len(groups) - 1:
#             legend = 'right margin'
#         else:
#             legend = None
#         sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

# plot_celltype_proportion(adata, out_dir=out_dir,type='group')
# plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

# method = "cell_type"
# keys = adata.obs[method].cat.categories.to_list()
# values = list(adata.uns[method+'_colors'])
# color_dict = dict(zip(keys,values))
# df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# # df.columns.names = ["cell_type"]
# df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Main_mod/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Main_mapping/"

adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Macrophages/Macrophages.h5ad")
cell_type_mapping = bdata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

# adata_raw = adata_raw[~adata_raw.obs['cell_type'].isin(['lowQ']), :]
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')



out="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Macrophages/cluster2/"
species='human'
os.system(f'/mnt/NexData_file/share/scRNA-tools/envs/pip_diffgene_enrich/bin/Rscript /mnt/NexData_file/home/lst/scrna-code/diff/do_enrichment.R -f {os.path.join(out, "differentially_expressed_genes.xls")} -s {species} -o /mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Macrophages/cluster2/')