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

adata_raw = sc.read_h5ad("/04.Mission/mapping_raw/Main_mapping.h5ad")
adata_raw.obs['cell_type'].cat.categories
# 'B/Plasma', 'doublets', 'Mast cells', 'Fibroblasts', 'T cells', 'AT2', 'Endothelial cells', 'Macrophages', 'AT1', 'pDCs', 'Ciliated cells', 'Neutrophils'
output_path="/09.Mission/Main_human/"
adata=adata_raw[adata_raw.obs['group'].isin(['EP','S', 'NS']), :]
sc.pl.umap(adata, color='cell_type', show=False)
save_fig(plt,file_name= 'dimplot_umap_celltype',out_dir=output_path,fig_size={'w':7,'h':6})

mainmarkers=['CNN1', 'ACTA2','TAGLN','C1R', 'COL1A2','DCN', 'PDGFRA']
sc.pl.umap(adata, color=mainmarkers, cmap='coolwarm')
plt.savefig(os.path.join(output_path, 'umap_mainmarkers.png'))
plt.savefig(os.path.join(output_path, 'umap_mainmarkers.pdf'))

sc.pl.umap(adata, color='leiden', show=False)
save_fig(plt,file_name= 'dimplot_umap_leiden',out_dir=output_path,fig_size={'w':7,'h':6})

anno_dict={
    '23':'SMCs'
}
mask = adata.obs['leiden'].isin(anno_dict.keys())
# 获取新的类别列表
new_categories = list(adata.obs['cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata.obs['cell_type'].cat.categories]
adata.obs['cell_type'] = adata.obs['cell_type'].cat.set_categories(new_categories)
adata.obs.loc[mask, 'cell_type'] = adata.obs.loc[mask, 'leiden'].map(anno_dict)
# 更新 categories
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(output_path, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(output_path, 'dimplot_umap_celltype.pdf'))

new_categories = list(adata.obs['sub_cell_type'].cat.categories) + [v for v in anno_dict.values() if v not in adata.obs['sub_cell_type'].cat.categories]
adata.obs['sub_cell_type'] = adata.obs['sub_cell_type'].cat.set_categories(new_categories)
adata.obs.loc[mask, 'sub_cell_type'] = adata.obs.loc[mask, 'leiden'].map(anno_dict)
# 更新 categories
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(output_path, 'umap_mapping.png'))
    plt.savefig(os.path.join(output_path, 'umap_mapping.pdf'))

adata.write(os.path.join(output_path, 'Main_mapping.h5ad'))


output_path="/09.Mission/Markerplot/"

adata=adata.raw.to_adata()
marker_genes=['HGF','MET']
sc.pl.umap(adata, color=marker_genes, show=False, vmax='p99')
save_fig(plt,file_name= 'Marker_umap',out_dir=output_path,fig_size={'w':12,'h':6})

adata.obs['group'] = pd.Categorical(adata.obs['group'], categories=['EP','S', 'NS'], ordered=True)
sc.pl.dotplot(adata, var_names=marker_genes, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'Marker_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'Marker_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata, var_names='HGF', groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'HGF_dotplot_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'HGF_dotplot_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata, var_names='MET', groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'MET_dotplot_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'MET_dotplot_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

gene = 'MET'                       # 换成你的基因名
dot_df = sc.get.obs_df(
            adata,
            keys=[gene,'cell_type', 'group'],
            use_raw=False                     # 如需 raw 改为 True
         )
# 1) 平均表达量
mean_df = (dot_df
           .groupby(['cell_type', 'group'])[gene]
           .mean()
           .rename('mean'))

# 2) 检测比例（>0 视为检测到）
pct_df = (dot_df
          .assign(detected=dot_df[gene] > 0)
          .groupby(['cell_type', 'group'])['detected']
          .mean()
          .mul(100)
          .rename('pct'))

# 3) 合并成最终表格（行索引 cell_type，列索引 group）
plot_df = pd.concat([mean_df, pct_df], axis=1).reset_index()

# 4. 画气泡图
plt.figure(figsize=(5,5))
sns.scatterplot(
    data=plot_df,
    x='group', y='cell_type',
    size='pct', hue='mean',
    sizes=(0, 300), palette='Reds', edgecolor='black', linewidth=0.5
)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title=None, frameon=False)
plt.title(''); plt.xlabel(''); plt.ylabel('')
plt.tight_layout()
plt.margins(x=0.25)
# sc.pl.dotplot(adata, var_names=genelist, groupby='cell_type', dendrogram=False, mean_only_expressed=True, show=False)
plt.savefig(os.path.join(output_path, 'MET_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'MET_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
plot_df.to_csv(os.path.join(output_path, 'MET_plotdata.xls'), sep='\t', index=None)