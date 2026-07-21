import pandas as pd
import anndata as ad
import numpy as np
import gc
import sys
gc.collect()
# 保存图片
import os
import matplotlib.pyplot as plt
plt.rcParams.update({
    'figure.figsize': (4,3),
    'font.size': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12
})
import scanpy as sc
from matplotlib import font_manager as fm
from matplotlib.pyplot import rc_context
sys.path.append('/mnt/NexData_file/share/scRNA-tools/scanpy_fun')
from subcluster import *
from draw_group import draw_group
sys.path.append('/mnt/NexData_file/home/zwl/Project/personal')
from commen_def import *
from find_marker import find_marker_gene

os.getcwd()

import itertools

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

def pairwise_find_marker_genes(adata,group_key: str,species: str,outdir_base: str,exclude_groups: list = None):
    """
    对 `group_key` 中的所有组进行两两差异表达分析
    
    参数:
        adata: AnnData 对象
        group_key: `.obs` 中的分组列名称
        species: 物种名称 (human/mouse)
        outdir_base: 输出结果的根目录
        exclude_groups: 可选，排除某些 group 不参与比较
    """
    groups = adata.obs[group_key].unique().tolist()
    if exclude_groups:
        groups = [g for g in groups if g not in exclude_groups]
    
    # 所有两两组合
    pair_combinations = list(itertools.combinations(groups, 2))
    for group1, group2 in pair_combinations:
        comp_name = f"{group1}_vs_{group2}"
        outdir = os.path.join(outdir_base, comp_name)
        os.makedirs(outdir, exist_ok=True)
        # 保留这两个组的数据
        adata_subset = adata[adata.obs[group_key].isin([group1, group2])].copy()
        print(f"▶ 正在处理: {group1} vs {group2}，细胞数: {adata_subset.n_obs}")
        find_marker_gene(
            adata=adata_subset,
            group=group_key,
            species=species,
            out_dir=outdir
        )

adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/Epithelial/Epithelial.h5ad")
for ci in adata.obs['cell_type'].unique():
    print(ci)
    out_dir=os.path.join('/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/Epithelial/groupDEG',ci)
    os.makedirs(out_dir,exist_ok=True)
    subadata = adata[adata.obs['cell_type'].isin([ci])]
    subadata.uns['log1p']['base'] = None
    pairwise_find_marker_genes(subadata,group_key='group',species='human',
                               outdir_base=out_dir)


adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/Myeloids_delete/cell_type/adata_anno.h5ad")
for ci in adata.obs['cell_type'].unique():
    print(ci)
    out_dir=os.path.join('/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/Myeloids_delete/groupDEG',ci)
    os.makedirs(out_dir,exist_ok=True)
    subadata = adata[adata.obs['cell_type'].isin([ci])]
    subadata.uns['log1p']['base'] = None
    pairwise_find_marker_genes(subadata,group_key='group',species='human',
                               outdir_base=out_dir)


adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/Neutrophils/cell_type/adata_anno.h5ad")
for ci in adata.obs['cell_type'].unique():
    print(ci)
    out_dir=os.path.join('/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/Neutrophils/groupDEG',ci)
    os.makedirs(out_dir,exist_ok=True)
    subadata = adata[adata.obs['cell_type'].isin([ci])]
    subadata.uns['log1p']['base'] = None
    pairwise_find_marker_genes(subadata,group_key='group',species='human',
                               outdir_base=out_dir)


adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/B/cell_type/adata_anno.h5ad")
for ci in adata.obs['cell_type'].unique():
    print(ci)
    out_dir=os.path.join('/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/B/groupDEG',ci)
    os.makedirs(out_dir,exist_ok=True)
    subadata = adata[adata.obs['cell_type'].isin([ci])]
    subadata.uns['log1p']['base'] = None
    pairwise_find_marker_genes(subadata,group_key='group',species='human',
                               outdir_base=out_dir)


adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD4T/CD4T.h5ad")
for ci in adata.obs['cell_type'].unique():
    print(ci)
    out_dir=os.path.join("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD4T/groupDEG/",ci)
    os.makedirs(out_dir,exist_ok=True)
    subadata = adata[adata.obs['cell_type'].isin([ci])]
    subadata.uns['log1p']['base'] = None
    pairwise_find_marker_genes(subadata,group_key='group',species='human',
                               outdir_base=out_dir)


adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD8T/CD8T.h5ad")
for ci in adata.obs['cell_type'].unique():
    print(ci)
    out_dir=os.path.join("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD8T/groupDEG/",ci)
    os.makedirs(out_dir,exist_ok=True)
    subadata = adata[adata.obs['cell_type'].isin([ci])]
    subadata.uns['log1p']['base'] = None
    pairwise_find_marker_genes(subadata,group_key='group',species='human',
                               outdir_base=out_dir)


#### CD4T 注释映射
mapping_df = pd.read_csv("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD4T/CD4T_cell_type_umap_mapping.csv", index_col=0)
adata = sc.read_h5ad("/mnt/NexData_file/home/wmx/Project/AE2507033447-08/20260413/subcluster/T_cells/result.h5ad")
out_dir = "/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD4T/"
adata_sub = adata[adata.obs_names.isin(mapping_df.index), :].copy()
# 按细胞名匹配更新 cell_type（mapping_df 的索引就是细胞名）
adata_sub.obs['cell_type'] = adata_sub.obs_names.map(mapping_df['cell_type'])

# 更新 UMAP 坐标
umap_dict = mapping_df[['UMAP_1', 'UMAP_2']].to_dict(orient='index')
umap_coords = np.array([
    [umap_dict[cell]['UMAP_1'], umap_dict[cell]['UMAP_2']] 
    for cell in adata_sub.obs_names
], dtype=np.float32)
# 赋值给 adata
adata_sub.obsm['X_umap'] = umap_coords
adata_sub.write("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD4T/updated_adata.h5ad")

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='group', title='Group', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.pdf'))

groups = ['non_pcr_B', 'pcr_B','non_pcr_A', 'pcr_A']
adata_sub.obs['group'] = (
    pd.Categorical(adata_sub.obs['group'],
                   categories=groups,
                   ordered=True)
)
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_sub[adata_sub.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_sub.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_sub, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata_sub, out_dir=out_dir,type='sample',fig_size={'w':15, 'h':5})
adata_sub.write(os.path.join(out_dir, 'CD4T.h5ad'))


#### CD8T 注释映射
mapping_df = pd.read_csv("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD8T/CD8T_cell_type_umap_mapping.csv", index_col=0)
adata = sc.read_h5ad("/mnt/NexData_file/home/wmx/Project/AE2507033447-08/20260413/subcluster/T_cells/result.h5ad")
out_dir = "/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD8T/"
adata_sub = adata[adata.obs_names.isin(mapping_df.index), :].copy()
# 按细胞名匹配更新 cell_type（mapping_df 的索引就是细胞名）
adata_sub.obs['cell_type'] = adata_sub.obs_names.map(mapping_df['cell_type'])
adata_sub.obs['RNA_snn_res.0.8'] = adata_sub.obs_names.map(mapping_df['RNA_snn_res.0.8'])
adata_sub.obs['cell_type'] = adata_sub.obs['cell_type'].astype('category')
adata_sub.obs['RNA_snn_res.0.8'] = adata_sub.obs['RNA_snn_res.0.8'].astype('category')
adata_sub.obs['cell_type'] = adata_sub.obs['cell_type'].cat.add_categories(['CD8_Tm_C1', 'CD8_Tm_C2', 'CD8_Tm_C3'])
mask = adata_sub.obs['RNA_snn_res.0.8'] == 0
adata_sub.obs.loc[mask, 'cell_type'] = 'CD8_Tm_C1'
mask = adata_sub.obs['RNA_snn_res.0.8'].isin([3,5])
adata_sub.obs.loc[mask, 'cell_type'] = 'CD8_Tm_C2'
mask = adata_sub.obs['RNA_snn_res.0.8'].isin([6,7])
adata_sub.obs.loc[mask, 'cell_type'] = 'CD8_Tm_C3'
adata_sub.obs['cell_type'] = adata_sub.obs['cell_type'].cat.remove_unused_categories()

# 更新 UMAP 坐标
umap_dict = mapping_df[['UMAP_1', 'UMAP_2']].to_dict(orient='index')
umap_coords = np.array([
    [umap_dict[cell]['UMAP_1'], umap_dict[cell]['UMAP_2']] 
    for cell in adata_sub.obs_names
], dtype=np.float32)
# 赋值给 adata
adata_sub.obsm['X_umap'] = umap_coords
adata_sub.write("/mnt/NexData_file/home/sh/Projects/AE2507033447-08/18.Mission/CD8T/updated_adata.h5ad")

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'),bbox_inches='tight')
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'),bbox_inches='tight')

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='sample', title='Sample', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_sample.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_sub, color='group', title='Group', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group.pdf'))

groups = ['non_pcr_B', 'pcr_B','non_pcr_A', 'pcr_A']
adata_sub.obs['group'] = (
    pd.Categorical(adata_sub.obs['group'],
                   categories=groups,
                   ordered=True)
)
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata_sub[adata_sub.obs["group"]==groups[i]], color="cell_type", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_celltype.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

adata_sub.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
plot_celltype_proportion(adata_sub, out_dir=out_dir,type='group',fig_size={'w':5, 'h':5})
plot_celltype_proportion(adata_sub, out_dir=out_dir,type='sample',fig_size={'w':15, 'h':5})
adata_sub.write(os.path.join(out_dir, 'CD8T.h5ad'))