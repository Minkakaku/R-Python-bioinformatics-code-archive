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

# 对CD8T进行细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/Tcells.h5ad")
adata_raw = adata_raw[adata_raw.obs.cell_type=="CD8T cell", :]
out_dir = "/mnt/NexData_file/home/sh/Projects/FWSC20250579/03.Mission/CD8T/"
adata_sub = recluster(adata=adata_raw,species='human',neighbors=20,resolution=0.2,n_pcs = 10,output=out_dir)
adata = adata_sub.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['Ptprc','Ccr10','Cxcr5','Ccr4','Ccr7','Sell','Cd44','Fasl'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_diffgenes.pdf'))
sc.pl.umap(adata, color=['Ccr7','Lef1','Sell','Tcf7','Cd27','Cd28','S1pr1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_naive_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_naive_diffgenes.pdf'))
sc.pl.umap(adata, color=['Ccr7','Sell','Il7r','Cd27','Cd28','Prf1','Gzma','Ccl5','Gpr183','S1pr1'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_cm_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_cm_diffgenes.pdf'))
sc.pl.umap(adata, color=['Klrg1','Cx3cr1','Fcgr3','Prf1','Gzmc','Tbx21','Eomes','S1pr1','S1pr5'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_eff_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_eff_diffgenes.pdf'))
sc.pl.umap(adata, color=['Gzmk','Cxcr4','Cxcr3','Cd44'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_em_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_em_diffgenes.pdf'))
sc.pl.umap(adata, color=['Cd6','Xcl1','Xcl1','Myadm','Capg','Rora','Nr4a1','Cd69','Itgae'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_rm_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_rm_diffgenes.pdf'))
sc.pl.umap(adata, color=['Cd160','Klrc1','Nr4a1','Ikzf2','Entpd1','Cd69','Itgae'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_iel_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_iel_diffgenes.pdf'))
sc.pl.umap(adata, color=['Havcr2','Cxcl13','Pdcd1','Layn','Tox','Ifng','Gzmb','Mir155hg','Tnfrsf9','Itgae'], cmap='coolwarm')
plt.savefig(os.path.join(out_dir, 'umap_ex_diffgenes.png'))
plt.savefig(os.path.join(out_dir, 'umap_ex_diffgenes.pdf'))


anno_dict={
    '0':'CD8 Teff',
    '1':'CD8 Tex',
    '2':'CD8 Trm'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

adata.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'CD8T.h5ad'))

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


# 映射所有细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/FWSC20250579/01.Mission/result_subcluster_Tcells/Tcells.h5ad")
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20250579/03.Mission/mapping_T/"
adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['sub_cell_type']), axis=1)

adata_raw.obs['sub_cell_type'] = adata_raw.obs['sub_cell_type'].astype('category')
adata_raw.obs['sub_cell_type'].cat.categories

sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin')
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_raw.write(os.path.join(out_dir, 'T_mapping.h5ad'))

adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')