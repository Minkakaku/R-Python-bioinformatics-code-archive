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

# Muscle cells注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Main_mapping/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'DC', 'Endothelial cells', 'Fibroblasts', 'Macrophages', 'Muscle cells', 'Pericytes', 'Plasma cells', 'Satellite cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/02.Mission/Muscle/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Muscle cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=10,resolution=0.3,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

# adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
# adata.obs['cell_type']=adata.obs['leiden']
# cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
# adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)
# adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
# adata_raw.obs['sub_cell_type'].cat.categories
# with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
#     sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False,legend_loc ='on data',legend_fontsize=5)
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
#     plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color='sample', title='Sample', show=False,legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_sample.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_sample.pdf'),bbox_inches='tight')
sc.pl.umap(adata_raw, color=['PAX7','MYF5','CD34','CALCR','VCAM1','SDC3','SDC4','ITGA7'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes_raw.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes_raw.pdf'))
sc.pl.dotplot(adata_raw, var_names=['PAX7','MYF5','CD34','CALCR','VCAM1','SDC3','SDC4','ITGA7'], groupby='cell_type', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_genes_raw.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_raw.pdf'),bbox_inches='tight')
sc.pl.umap(adata, color=['PAX7','MYF5','CD34','CALCR','VCAM1','SDC3','SDC4','ITGA7','MYOD1','MYF5','DES','CDH15','MYOG','MYF6','MEF2A','MEF2C','MEF2D','DES','MYH3','MYH1','MYH2','MYH4','ACTA1','TNNT1','TNNT3','CKM','DMD'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['PAX7','MYF5','CD34','CALCR','VCAM1','SDC3','SDC4','ITGA7','MYOD1','MYF5','DES','CDH15','MYOG','MYF6','MEF2A','MEF2C','MEF2D','DES','MYH3','MYH1','MYH2','MYH4','ACTA1','TNNT1','TNNT3','CKM','DMD'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'TNNT1 Myotubes',
    '1':'Myocytes',
    '2':'MYH1 Myotubes',
    '3':'MYH2 Myotubes',
    '4':'Myoblasts',
    '5':'Myoblasts'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Muscle.h5ad'))

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

out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/02.Mission/Muscle/ERV_LTR/"
sc.pl.umap(adata, color=['LTR/ERVK','LTR/ERV1','LTR/ERVL','LTR/ERVL-MaLR'], cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['LTR/ERVK','LTR/ERV1','LTR/ERVL','LTR/ERVL-MaLR'], groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['LTR/ERVK','LTR/ERV1','LTR/ERVL','LTR/ERVL-MaLR'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.pdf'),bbox_inches='tight')

# Endothelial cells注释
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Main_mapping/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B cells', 'DC', 'Endothelial cells', 'Fibroblasts', 'Macrophages', 'Muscle cells', 'Pericytes', 'Plasma cells', 'Satellite cells', 'T cells'
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/02.Mission/Endothelial/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Endothelial cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=15,resolution=0.5,n_pcs = 15,output=out_dir,batch_effect = 1)
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

sc.pl.umap(adata, color=['SEMA3C','KIT','GPIHBP1','CCL21', 'PROX1','HEY1', 'IGFBP3','CD36', 'CA4','ACKR1'], cmap='coolwarm',vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['SEMA3C','KIT','GPIHBP1','CCL21', 'PROX1','HEY1', 'IGFBP3','CD36', 'CA4','ACKR1'], groupby='leiden', dendrogram=False, show=False,standard_scale="var")
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')

anno_dict={
    '0':'PKD1L1 Capillaries',
    '1':'IFI44L Capillaries',
    '2':'vein ECs',
    '3':'RGS5 Capillaries',
    '4':'artery ECs',
    '5':'lymphatic ECs'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.png"),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir,"dimplot_umap_celltype.pdf"),bbox_inches='tight')

sc.tl.rank_genes_groups(adata, groupby="cell_type", method='wilcoxon', pts=True)
find_group_marker(adata=adata,species='human',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata.write(os.path.join(out_dir, 'Endothelial.h5ad'))

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

out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/02.Mission/Endothelial/ERV_LTR/"
sc.pl.umap(adata, color=['LTR/ERVK','LTR/ERV1','LTR/ERVL','LTR/ERVL-MaLR'], cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_genes.png'))
plt.savefig(os.path.join(out_dir, 'umap_genes.pdf'))
sc.pl.dotplot(adata, var_names=['LTR/ERVK','LTR/ERV1','LTR/ERVL','LTR/ERVL-MaLR'], groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=['LTR/ERVK','LTR/ERV1','LTR/ERVL','LTR/ERVL-MaLR'], groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_genes_cell_type.pdf'),bbox_inches='tight')

# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/01.Mission/Main_mapping/Main_mapping.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/02.Mission/Main_mapping/"
# adata_raw.obs['sub_cell_type']=adata_raw.obs['cell_type']
cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

bdata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2509094832-01/02.Mission/Muscle/Muscle.h5ad")
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



adata=sc.read("/mnt/NexData_file/scRNA/scVision_output/SC2509094832_ERV_all/result.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/SC2509094832-01/02.Mission/ERV_all/"
geneset1=['ERV_0002_LTR103_Mam','ERV_0005_MER21B','ERV_0006_MLT1E1A','ERV_0007_MLT2B3','ERV_0010_THE1C-int','ERV_0011_LTR1D1','ERV_0014_MER39','ERV_0015_MLT1B','ERV_0018_LTR12','ERV_0019_THE1B-int','ERV_0020_LTR8','ERV_0022_MER4E1','ERV_0024_MER50','ERV_0027_LTR8']
geneset2=['ERV_0035_MER21C','ERV_0039_LTR2','ERV_0042_MLT1H','ERV_0043_MLT1H','ERV_0045_MLT1B','ERV_0047_HERVL74-int','ERV_0054_HUERS-P3-int','ERV_0055_MER101-int','ERV_0057_MLT1B','ERV_0058_MLT1B-int','ERV_0060_LTR5_Hs','ERV_0061_MLT1H-int','ERV_0063_LTR49-int','ERV_0065_MLT1J']
geneset3=['ERV_0069_MLT1F1','ERV_0071_MLT1A0','ERV_0072_MLT1B','ERV_0073_MLT1H2','ERV_0076_MER4D','ERV_0077_MER51B','ERV_0079_THE1B-int','ERV_0080_THE1B','ERV_0089_LTR57','ERV_0090_MER65C','ERV_0092_MER57F','ERV_0093_MLT1A1','ERV_0094_MamGypLTR2c','ERV_0095_LTR5_Hs','ERV_0096_HUERS-P2-int','ERV_0099_THE1C']
geneset4=['ERV_0101_MER51A','ERV_0102_PABL_B','ERV_0104_LTR12C','ERV_0105_MLT1D','ERV_0106_MLT1A0','ERV_0107_LTR56','ERV_0108_MER21C','ERV_0113_LTR8B','ERV_0114_LTR64','ERV_0119_Harlequin-int','ERV_0120_LTR103b_Mam','ERV_0123_LTR8A','ERV_0129_LTR2','ERV_0130_MLT1L']
sc.pl.umap(adata, color=geneset1, cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_geneset1.png'))
plt.savefig(os.path.join(out_dir, 'umap_geneset1.pdf'))
sc.pl.dotplot(adata, var_names=geneset1, groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=geneset1, groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset1_cell_type.pdf'),bbox_inches='tight')

sc.pl.umap(adata, color=geneset2, cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_geneset2.png'))
plt.savefig(os.path.join(out_dir, 'umap_geneset2.pdf'))
sc.pl.dotplot(adata, var_names=geneset2, groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset2.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset2.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=geneset2, groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset2_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset2_cell_type.pdf'),bbox_inches='tight')

sc.pl.umap(adata, color=geneset3, cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_geneset3.png'))
plt.savefig(os.path.join(out_dir, 'umap_geneset3.pdf'))
sc.pl.dotplot(adata, var_names=geneset3, groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset3.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset3.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=geneset3, groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset3_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset3_cell_type.pdf'),bbox_inches='tight')

sc.pl.umap(adata, color=geneset4, cmap='coolwarm', ncols=3,vmax='p99')
plt.savefig(os.path.join(out_dir, 'umap_geneset4.png'))
plt.savefig(os.path.join(out_dir, 'umap_geneset4.pdf'))
sc.pl.dotplot(adata, var_names=geneset4, groupby=['cell_type','group'], dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset4.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset4.pdf'),bbox_inches='tight')
sc.pl.dotplot(adata, var_names=geneset4, groupby='cell_type', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset4_cell_type.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'dotplot_geneset4_cell_type.pdf'),bbox_inches='tight')