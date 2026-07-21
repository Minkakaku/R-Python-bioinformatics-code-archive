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

# adata_sub=adata[adata.obs['sample'].isin(['D4_1','D7_1','D14_1','D21_1','D4_2','D7_2','D14_2','D21_2'])]
# 对pericytes注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2512027072-01/01.Mission/Main_mapping/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'Keratinocytes', 'Macrophages', 'Neutrophils', 'Schwann cells', 'Smooth muscle cells', 'T cells', 'pericytes'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512027072-01/02.Mission/pericytes/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="pericytes", :]
# adata_sub.uns['log1p']['base'] = 2
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.5,n_pcs = 10,output=out_dir)
adata = adata_recluster.sub_cluster()

# （1）血管舒张相关周细胞（vdPCs）：GUCY1A3','GUCY1B3
# （2）脂肪细胞样周细胞（adiPCs）：APOE','APOD','PTGDS
# （3）肌样周细胞（myoPCs）：ACTA2','TAGLN
# （4）基质生成型周细胞（matPCs）：COL1A1','COL3A1','POSTN','BASP1','SPP1
sc.pl.umap(adata, color=['Gucy1a1','Gucy1b1','Apoe','Apod','Ptgds','Acta2','Tagln','Col1a1','Col3a1','Postn','Basp1','Spp1'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_markergenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_markergenes.pdf'))
sc.pl.dotplot(adata, var_names=['Gucy1a1','Gucy1b1','Apoe','Apod','Ptgds','Acta2','Tagln','Col1a1','Col3a1','Postn','Basp1','Spp1'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'), bbox_inches="tight")

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
    '0':'vdPCs',
    '1':'matPCs',
    '2':'adiPCs',
    '3':'myoPCs',
    '4':'myoPCs',
    '5':'vdPCs',
    '6':'matPCs',
    '7':'adiPCs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'pericytes.h5ad'))

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

desired_order=['D4', 'D7', 'D14', 'D21']
adata.obs['group'] = (
    pd.Categorical(adata.obs['group'],
                   categories=desired_order,
                   ordered=True)
)
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})

method = "cell_type"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='cell_type')

adata.obs['group']=adata.obs['sample']
anno_dict={
    'D4_1':'group_1',
    'D4_2':'group_2',
    'D4_3':'group_3',
    'D7_1':'group_1',
    'D7_2':'group_2',
    'D7_3':'group_3',
    'D14_1':'group_1',
    'D14_2':'group_2',
    'D14_3':'group_3',
    'D21_1':'group_1',
    'D21_2':'group_2',
    'D21_3':'group_3'
}
adata.obs['group'] = [anno_dict[clust] for clust in adata.obs['sample']]
adata.obs['group'] = adata.obs['group'].astype('category')
adata=adata[adata.obs['group'].isin(['group_1', 'group_2'])]
adata.obs.to_csv(os.path.join(out_dir,"all_cells_features_sub.xls"), sep='\t')
method = "group"
keys = adata.obs[method].cat.categories.to_list()
values = list(adata.uns[method+'_colors'])
color_dict = dict(zip(keys,values))
df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# df.columns.names = ["cell_type"]
df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')
adata.write(os.path.join(out_dir, 'pericytes_sub.h5ad'))
# adata_non_mito = adata[:, ~adata.var_names.str.startswith('MT-')]
# sc.tl.rank_genes_groups(adata_non_mito, groupby='cell_type')
# sc.pl.rank_genes_groups_heatmap(adata_non_mito, n_genes=10, dendrogram=False, standard_scale='var', swap_axes=True, show_gene_labels=True, cmap='viridis')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')
# sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, dendrogram=False,standard_scale='var',swap_axes=True,show_gene_labels=True,cmap='viridis')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.png'),bbox_inches='tight')
# plt.savefig(os.path.join(out_dir, 'heatmap_top10genes.pdf'),bbox_inches='tight')

# top50 = sc.get.rank_genes_groups_df(adata, group=None)   # group=None = 全部群
# top50 = top50.groupby('group').head(50).reset_index(drop=True)
# top50.to_excel(os.path.join(out_dir, 'top50_diff_genes.xlsx'), index=False)
# # method = "group"
# # keys = adata.obs[method].cat.categories.to_list()
# # values = list(adata.uns[method+'_colors'])
# # color_dict = dict(zip(keys,values))
# # df = pd.DataFrame.from_dict(color_dict,orient='index',columns=['color'])
# # # df.columns.names = ["cell_type"]
# # df.to_csv(os.path.join(out_dir,method+'color_dict.xls'),sep="\t",index=True,index_label='group')
# adata_sub=adata[adata.obs['cell_type'].isin(['HTC1','HTC2']), :]
# out_dir="/mnt/NexData_file/home/sh/Projects/SC2509265318-02/03.Mission/Chondrocyte/HTC1_HTC2/"
# sc.tl.rank_genes_groups(adata_sub, groupby='cell_type', method='wilcoxon', pts=True)
# find_group_marker(adata=adata_sub,species='mouse',output=out_dir)



# 映射所有细分 
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/SC2512027072/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512027072-01/01.Mission/Main_mapping/"

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)

adata_raw.obs['cell_type'] = adata_raw.obs['cell_type'].astype('category')
adata_raw.obs['cell_type'].cat.categories

# sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
# plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
# plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.umap(adata_raw, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')

plot_celltype_proportion(adata_raw, out_dir=out_dir,type='group')
plot_celltype_proportion(adata_raw, out_dir=out_dir,type='sample',fig_size={'w':10, 'h':5})

adata=adata_raw[adata_raw.obs['sample'].isin(['D21_1', 'D21_2'])]
out_dir="/mnt/NexData_file/home/sh/Projects/SC2512027072-01/01.Mission/Main_mapping_sub/"
sc.pl.umap(adata, color='cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata.write(os.path.join(out_dir, 'Main_mapping.h5ad'))
adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')

plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')

sc.pl.violin(adata, ['Rgs5','Cspg4','Abcc9','Mcam','Myh11','Mcam','Acta2'], jitter=0, groupby='cell_type', rotation=45, multi_panel=True, show=False)
plt.savefig(os.path.join(out_dir, 'violin_markergenes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'violin_markergenes.pdf'), bbox_inches="tight")
sc.pl.stacked_violin(adata, var_names=['Rgs5','Cspg4','Abcc9','Mcam','Myh11','Mcam','Acta2'], groupby='cell_type', show=False, swap_axes=True, dendrogram=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'stacked_violin_markergenes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'stacked_violin_markergenes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['Abcc9','Mcam','Myh11','Acta2'], groupby='cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.png'), bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_markergenes.pdf'), bbox_inches="tight")