import os, sys
import numpy as np
import scanpy as sc
import matplotlib
# from scanpy_fun_subcluster import recluster
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
# from scanpy_fun_subcluster import find_group_marker
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

adata_raw = sc.read_h5ad("/adata_cluster.h5ad")
df = pd.read_table("/06.Mission/cellchat_all.txt", header=0)
adata_raw.obs['barcode'] = (
    adata_raw.obs['sample'].astype(str) + '_' +
    adata_raw.obs_names.str[:17]
)
mapping = df.set_index('orig.ident')['cell_type'].to_dict()
adata_raw.obs['cell_type'] = (
    adata_raw.obs['barcode']
    .map(mapping)
    .fillna('unknown')
    .astype('category')
)
adata = adata_raw[~adata_raw.obs['cell_type'].isin(['unknown']), :]
# adata_raw.obs['cell_type'].cat.categories
output_path="/09.Mission/Markerplot_mouse/"

df_sub = pd.read_table("/08.Mission/Markerplot_mouse/Epi_cells_features.txt", header=0, index_col=0)   # 第 0 列就是 orig.ident
mapping = df_sub['cell_type'].to_dict()
adata.obs['cell_type'] = (
    adata.obs['barcode']
    .map(mapping)               # 映射
    .combine_first(adata.obs['cell_type'])  # 缺失的用原值
    .astype('category')          # 先转成 category
)
adata.obs['cell_type'].cat.categories

sc.pl.umap(adata, color='cell_type', show=False)
save_fig(plt,file_name= 'dimplot_umap_celltype',out_dir=output_path,fig_size={'w':9,'h':6})
# adata=adata.raw.to_adata()
marker_genes=['Hgf','Met']
sc.pl.umap(adata, color=marker_genes, show=False, vmax='p99')
save_fig(plt,file_name= 'Marker_umap',out_dir=output_path,fig_size={'w':12,'h':6})

# adata=adata[adata.obs['group'].isin(['ELA', 'Control']), :]

sc.pl.dotplot(adata, var_names=marker_genes, groupby='cell_type', dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'Marker_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'Marker_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata, var_names='Hgf', groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'Hgf_dotplot_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'Hgf_dotplot_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
sc.pl.dotplot(adata, var_names='Met', groupby=['cell_type','group'], dendrogram=False, show=False)
plt.savefig(os.path.join(output_path, 'Met_dotplot_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'Met_dotplot_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)


gene = 'Hgf'                       # 换成你的基因名
dot_df = sc.get.obs_df(
            adata,
            keys=[gene,'cell_type', 'group'],
            use_raw=True                     # 如需 raw 改为 True
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
plt.savefig(os.path.join(output_path, 'Hgf_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(output_path, 'Hgf_dotplot.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
plot_df.to_csv(os.path.join(output_path, 'Hgf_plotdata.xls'), sep='\t', index=None)