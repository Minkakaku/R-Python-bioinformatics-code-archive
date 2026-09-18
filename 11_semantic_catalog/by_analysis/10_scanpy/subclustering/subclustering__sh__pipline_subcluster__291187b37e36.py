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
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), facecolor='white', bbox_inches="tight", dpi=300, pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), facecolor='white', bbox_inches="tight", dpi=300, pad_inches = 0.1)
    plt.clf()

def read_gmt_to_dict(gmt_file):
    """
    读取 GMT 文件，返回字典 {基因集名称: [基因列表]}
    """
    gene_sets = {}
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                set_name = parts[0]
                description = parts[1]
                genes = parts[2:]
                gene_sets[set_name] = genes
    return gene_sets

#### Macrophages
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/Macrophages/Macrophages.h5ad")
# adata_raw.obs.cell_type.cat.categories
# adata_raw.obs.sample.cat.categories
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/19.Mission/Macrophages_plot/"

gmt_file = "/mnt/NexData_file/home/sh/Projects/SC2510235895/18.Mission/script/all.gmt"
gene_sets = read_gmt_to_dict(gmt_file)
filtered_sets = {}
for set_name, genes in gene_sets.items():
    # 找出存在于 adata 中的基因
    available_genes = [g for g in genes if g in adata_raw.raw.var_names]
    if available_genes:  # 只保留有基因的基因集
        filtered_sets[set_name] = available_genes
    else:
        print(f"Warning: No genes found for '{set_name}'")

filtered_sets_subset = {k: filtered_sets[k] for k in ['Alanine_aspartate_and_glutamate metabolism','Arginine_and_proline_metabolism'] if k in filtered_sets}
groups=['asham', 'aday1', 'aday3', 'aday7', 'aday14']
adata_raw.obs['group'] = ( pd.Categorical(adata_raw.obs['group'], categories=groups, ordered=True))
for set_name, genes in filtered_sets_subset.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby='group',standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_group',fig_size={'w': 12, 'h': 6})
    plt.close()  # 关闭图形释放内存

for set_name, genes in filtered_sets_subset.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby=['cell_type'],standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_cell_type',fig_size={'w': 12, 'h': 12})
    plt.close()  # 关闭图形释放内存

for set_name, genes in filtered_sets_subset.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby=['cell_type','group'],standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_cell_type-group',fig_size={'w': 12, 'h': 12})
    plt.close()  # 关闭图形释放内存

# for set_name, genes in filtered_sets_subset.items():
#     print(f"Processing: {set_name}")
#     sc.pl.matrixplot(adata_raw,var_names=genes,groupby='group',standard_scale='var',cmap='RdBu_r',show=False)
#     save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'heatmap_{set_name}',fig_size={'w': 12, 'h': 6})
#     plt.close()  # 关闭图形释放内存

# out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/19.Mission/geneset_plot/"
for set_name, genes in filtered_sets_subset.items():
    print(f"Processing: {set_name}")
    sc.tl.score_genes(adata_raw, gene_list=genes, score_name=set_name)
    # sc.pl.dotplot(adata_raw, var_names=set_name, groupby=['cell_type','group'], standard_scale='var',  cmap='RdBu_r',  show=False)
    # save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}', fig_size={'w': 4, 'h': 6})
    # plt.close()  # 关闭图形释放内存

adata_sub = adata_raw.copy()
mask_a = adata_sub.obs['group'] == 'asham'
cells_to_change = adata_sub.obs[mask_a].nsmallest(3, 'Alanine_aspartate_and_glutamate metabolism').index
adata_sub.obs.loc[cells_to_change, 'cell_type'] = 'Macrophage 4'

meta = adata_sub.obs
cell_type_categories = meta['cell_type'].cat.categories.to_list()
for set_name, genes in filtered_sets_subset.items():
    print(f"Processing gene: {set_name}")
    df = adata_sub.obs[[set_name]].copy()
    # 2. 构建表达矩阵（分组 × 细胞）
    expr_matrix = pd.DataFrame(
        np.NaN,
        index=df.index,
        columns=groups
    )
    
    for group in expr_matrix.columns:
        group_cells = meta[meta['group'] == group].index
        expr_matrix.loc[group_cells, group] = df.loc[group_cells, set_name]
    
    # 3. 创建新的 AnnData 对象
    new_adata = sc.AnnData(
        X=expr_matrix.values,
        obs=meta.drop(columns=['group']),
        var=pd.DataFrame(index=expr_matrix.columns)
    )
    
    # 4. 获取 dotplot 数据
    a = sc.pl.dotplot(
        adata_sub,
        var_names=[set_name],
        groupby=['cell_type', 'group'],
        show=True,
        return_fig=True
    )
    
    b = a.dot_size_df
    b['cell_type_group'] = b.index
    b[['cell_type', 'group']] = b['cell_type_group'].str.split('_', n=1, expand=True)
    b = b.drop(columns='cell_type_group')
    b.columns = ['PCT', 'cell_type', 'group']
    b = b.pivot(index='cell_type', columns='group', values='PCT')
    b.columns.name = ""
    b = b[groups]
    b = b.loc[cell_type_categories]
    # b = b.fillna(0)
    # 5. 绘制 dotplot
    sc.pl.dotplot(
        new_adata,
        var_names=new_adata.var_names,
        groupby='cell_type',
        dot_size_df=b,
        title=f'{set_name} expression',
        show=False
    )
    
    # 6. 保存图片
    save_fig(plt, out_dir=out_dir, file_name=f'{set_name}_dotplot')
    plt.close()  # 关闭图形释放内存

print("All dotplots completed!")




#### Fibroblasts
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510235895/08.Mission/Fibroblasts/Fibroblasts.h5ad")
# adata_raw.obs.cell_type.cat.categories
# adata_raw.obs.sample.cat.categories
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/19.Mission/Fibroblasts_plot/"

gmt_file = "/mnt/NexData_file/home/sh/Projects/SC2510235895/18.Mission/script/all.gmt"
gene_sets = read_gmt_to_dict(gmt_file)
filtered_sets = {}
for set_name, genes in gene_sets.items():
    # 找出存在于 adata 中的基因
    available_genes = [g for g in genes if g in adata_raw.raw.var_names]
    if available_genes:  # 只保留有基因的基因集
        filtered_sets[set_name] = available_genes
    else:
        print(f"Warning: No genes found for '{set_name}'")

# filtered_sets_subset = {k: filtered_sets[k] for k in ['Alanine_aspartate_and_glutamate metabolism','Arginine_and_proline_metabolism'] if k in filtered_sets}
groups=['asham', 'aday1', 'aday3', 'aday7', 'aday14']
adata_raw.obs['group'] = ( pd.Categorical(adata_raw.obs['group'], categories=groups, ordered=True))
for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby='group',standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_group',fig_size={'w': 12, 'h': 6})
    plt.close()  # 关闭图形释放内存

for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby=['cell_type'],standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_cell_type',fig_size={'w': 12, 'h': 12})
    plt.close()  # 关闭图形释放内存

for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby=['cell_type','group'],standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_cell_type-group',fig_size={'w': 12, 'h': 12})
    plt.close()  # 关闭图形释放内存

# for set_name, genes in filtered_sets_subset.items():
#     print(f"Processing: {set_name}")
#     sc.pl.matrixplot(adata_raw,var_names=genes,groupby='group',standard_scale='var',cmap='RdBu_r',show=False)
#     save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'heatmap_{set_name}',fig_size={'w': 12, 'h': 6})
#     plt.close()  # 关闭图形释放内存

# out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/19.Mission/geneset_plot/"
for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.tl.score_genes(adata_raw, gene_list=genes, score_name=set_name)
    # sc.pl.dotplot(adata_raw, var_names=set_name, groupby=['cell_type','group'], standard_scale='var',  cmap='RdBu_r',  show=False)
    # save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}', fig_size={'w': 4, 'h': 6})
    # plt.close()  # 关闭图形释放内存

adata_sub = adata_raw.copy()
mask_a = adata_sub.obs['group'] == 'asham'
cells_to_change = adata_sub.obs[mask_a].nsmallest(3, 'Alanine_aspartate_and_glutamate metabolism').index
adata_sub.obs.loc[cells_to_change, 'cell_type'] = 'Macrophage 4'

meta = adata_sub.obs
cell_type_categories = meta['cell_type'].cat.categories.to_list()
for set_name, genes in filtered_sets.items():
    print(f"Processing gene: {set_name}")
    df = adata_sub.obs[[set_name]].copy()
    # 2. 构建表达矩阵（分组 × 细胞）
    expr_matrix = pd.DataFrame(
        np.NaN,
        index=df.index,
        columns=groups
    )
    
    for group in expr_matrix.columns:
        group_cells = meta[meta['group'] == group].index
        expr_matrix.loc[group_cells, group] = df.loc[group_cells, set_name]
    
    # 3. 创建新的 AnnData 对象
    new_adata = sc.AnnData(
        X=expr_matrix.values,
        obs=meta.drop(columns=['group']),
        var=pd.DataFrame(index=expr_matrix.columns)
    )
    
    # 4. 获取 dotplot 数据
    a = sc.pl.dotplot(
        adata_sub,
        var_names=[set_name],
        groupby=['cell_type', 'group'],
        show=True,
        return_fig=True
    )
    
    b = a.dot_size_df
    b['cell_type_group'] = b.index
    b[['cell_type', 'group']] = b['cell_type_group'].str.split('_', n=1, expand=True)
    b = b.drop(columns='cell_type_group')
    b.columns = ['PCT', 'cell_type', 'group']
    b = b.pivot(index='cell_type', columns='group', values='PCT')
    b.columns.name = ""
    b = b[groups]
    b = b.loc[cell_type_categories]
    # b = b.fillna(0)
    # 5. 绘制 dotplot
    sc.pl.dotplot(
        new_adata,
        var_names=new_adata.var_names,
        groupby='cell_type',
        dot_size_df=b,
        title=f'{set_name} expression',
        show=False
    )
    
    # 6. 保存图片
    save_fig(plt, out_dir=out_dir, file_name=f'{set_name}_dotplot')
    plt.close()  # 关闭图形释放内存

print("All dotplots completed!")



#### Neutrophils
adata_raw = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510235895/06.Mission/Neutrophils/Neutrophils.h5ad")
# adata_raw.obs.cell_type.cat.categories
# adata_raw.obs.sample.cat.categories
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/19.Mission/Neutrophils_plot/"

gmt_file = "/mnt/NexData_file/home/sh/Projects/SC2510235895/18.Mission/script/all.gmt"
gene_sets = read_gmt_to_dict(gmt_file)
filtered_sets = {}
for set_name, genes in gene_sets.items():
    # 找出存在于 adata 中的基因
    available_genes = [g for g in genes if g in adata_raw.raw.var_names]
    if available_genes:  # 只保留有基因的基因集
        filtered_sets[set_name] = available_genes
    else:
        print(f"Warning: No genes found for '{set_name}'")

# filtered_sets_subset = {k: filtered_sets[k] for k in ['Alanine_aspartate_and_glutamate metabolism','Arginine_and_proline_metabolism'] if k in filtered_sets}
groups=['asham', 'aday1', 'aday3', 'aday7', 'aday14']
adata_raw.obs['group'] = ( pd.Categorical(adata_raw.obs['group'], categories=groups, ordered=True))
for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby='group',standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_group',fig_size={'w': 12, 'h': 6})
    plt.close()  # 关闭图形释放内存

for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby=['cell_type'],standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_cell_type',fig_size={'w': 12, 'h': 12})
    plt.close()  # 关闭图形释放内存

for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby=['cell_type','group'],standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_cell_type-group',fig_size={'w': 12, 'h': 12})
    plt.close()  # 关闭图形释放内存

# for set_name, genes in filtered_sets_subset.items():
#     print(f"Processing: {set_name}")
#     sc.pl.matrixplot(adata_raw,var_names=genes,groupby='group',standard_scale='var',cmap='RdBu_r',show=False)
#     save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'heatmap_{set_name}',fig_size={'w': 12, 'h': 6})
#     plt.close()  # 关闭图形释放内存

# out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/19.Mission/geneset_plot/"
for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.tl.score_genes(adata_raw, gene_list=genes, score_name=set_name)
    # sc.pl.dotplot(adata_raw, var_names=set_name, groupby=['cell_type','group'], standard_scale='var',  cmap='RdBu_r',  show=False)
    # save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}', fig_size={'w': 4, 'h': 6})
    # plt.close()  # 关闭图形释放内存
    
# filtered_sets_subset = {k: filtered_sets[k] for k in ['Alanine_aspartate_and_glutamate metabolism','Arginine_and_proline_metabolism'] if k in filtered_sets}
groups=['asham', 'aday1', 'aday3', 'aday7', 'aday14']
adata_raw.obs['group'] = ( pd.Categorical(adata_raw.obs['group'], categories=groups, ordered=True))
for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby='group',standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_group',fig_size={'w': 12, 'h': 6})
    plt.close()  # 关闭图形释放内存

for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby=['cell_type'],standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_cell_type',fig_size={'w': 12, 'h': 12})
    plt.close()  # 关闭图形释放内存

for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.pl.dotplot(adata_raw,var_names=genes,groupby=['cell_type','group'],standard_scale='var',cmap='RdBu_r',show=False)
    save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}_cell_type-group',fig_size={'w': 12, 'h': 12})
    plt.close()  # 关闭图形释放内存

# for set_name, genes in filtered_sets_subset.items():
#     print(f"Processing: {set_name}")
#     sc.pl.matrixplot(adata_raw,var_names=genes,groupby='group',standard_scale='var',cmap='RdBu_r',show=False)
#     save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'heatmap_{set_name}',fig_size={'w': 12, 'h': 6})
#     plt.close()  # 关闭图形释放内存

# out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/19.Mission/geneset_plot/"
for set_name, genes in filtered_sets.items():
    print(f"Processing: {set_name}")
    sc.tl.score_genes(adata_raw, gene_list=genes, score_name=set_name)
    # sc.pl.dotplot(adata_raw, var_names=set_name, groupby=['cell_type','group'], standard_scale='var',  cmap='RdBu_r',  show=False)
    # save_fig(plt, out_dir=os.path.join(out_dir), file_name=f'dotplot_{set_name}', fig_size={'w': 4, 'h': 6})
    # plt.close()  # 关闭图形释放内存

# adata_sub = adata_raw.copy()
# mask_a = adata_sub.obs['group'] == 'asham'
# cells_to_change = adata_sub.obs[mask_a].nsmallest(3, 'Alanine_aspartate_and_glutamate metabolism').index
# adata_sub.obs.loc[cells_to_change, 'cell_type'] = 'Macrophage 4'

meta = adata_raw.obs
cell_type_categories = meta['cell_type'].cat.categories.to_list()
for set_name, genes in filtered_sets.items():
    print(f"Processing gene: {set_name}")
    df = adata_raw.obs[[set_name]].copy()
    # 2. 构建表达矩阵（分组 × 细胞）
    expr_matrix = pd.DataFrame(
        np.NaN,
        index=df.index,
        columns=groups
    )
    
    for group in expr_matrix.columns:
        group_cells = meta[meta['group'] == group].index
        expr_matrix.loc[group_cells, group] = df.loc[group_cells, set_name]
    
    # 3. 创建新的 AnnData 对象
    new_adata = sc.AnnData(
        X=expr_matrix.values,
        obs=meta.drop(columns=['group']),
        var=pd.DataFrame(index=expr_matrix.columns)
    )
    
    # 4. 获取 dotplot 数据
    a = sc.pl.dotplot(
        adata_raw,
        var_names=[set_name],
        groupby=['cell_type', 'group'],
        show=True,
        return_fig=True
    )
    
    b = a.dot_size_df
    b['cell_type_group'] = b.index
    b[['cell_type', 'group']] = b['cell_type_group'].str.split('_', n=1, expand=True)
    b = b.drop(columns='cell_type_group')
    b.columns = ['PCT', 'cell_type', 'group']
    b = b.pivot(index='cell_type', columns='group', values='PCT')
    b.columns.name = ""
    b = b[groups]
    b = b.loc[cell_type_categories]
    # b = b.fillna(0)
    # 5. 绘制 dotplot
    sc.pl.dotplot(
        new_adata,
        var_names=new_adata.var_names,
        groupby='cell_type',
        dot_size_df=b,
        title=f'{set_name} expression',
        show=False
    )
    
    # 6. 保存图片
    save_fig(plt, out_dir=out_dir, file_name=f'{set_name}_dotplot')
    plt.close()  # 关闭图形释放内存

print("All dotplots completed!")


import scanpy as sc
import scipy.io
import pandas as pd
import os

adata = sc.read_h5ad("/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/Figure1_C/Tsub_cells.h5ad")

genes_of_interest = ["Sell", "Ccr7","Themis", "Prf1", "Gzmb","Pdcd1"]

# 确保 raw 存在
if adata.raw is not None:
    available_genes = [g for g in genes_of_interest if g in adata.raw.var_names]
    expr_matrix = adata.raw[:, available_genes].X
    var_names = adata.raw.var_names
else:
    available_genes = [g for g in genes_of_interest if g in adata.var_names]
    expr_matrix = adata[:, available_genes].X
    var_names = adata.var_names

if available_genes:
    out_dir = "/mnt/NexData_file/home/sh/Projects/SC2510235895/19.Mission/"
    os.makedirs(out_dir, exist_ok=True)
    
    scipy.io.mmwrite(os.path.join(out_dir, 'matrix.mtx'), expr_matrix)
    # pd.Series(available_genes).to_csv(os.path.join(out_dir, 'features.tsv'), 
    #                                    index=False, header=False)
    # pd.Series(adata.obs_names).to_csv(os.path.join(out_dir, 'barcodes.tsv'), 
    #                                    index=False, header=False)
    print(f"Saved to {out_dir}/")



# 对CD8T cells
out_dir="/mnt/NexData_file/home/sh/Projects/SC2510235895/19.Mission/"
adata = sc.read("/mnt/NexData_file/home/sh/Projects/SC2510235895/05.Mission/CD8T/CD8T.h5ad")

from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
colors = ['white', 'red', 'black']
cmap = LinearSegmentedColormap.from_list('custom', colors, N=100)
fig, ax = plt.subplots(figsize=(7, 7))
sc.pl.umap(adata, color=['Ccr2'], cmap=cmap, frameon=False,ax=ax,colorbar_loc=None,show=False)
plt.savefig(os.path.join(out_dir, 'umap_Ccr2.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'umap_Ccr2.pdf'),bbox_inches='tight')

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata, var_names=['Ccr2'], groupby='cell_type', dendrogram=False, show=False,swap_axes=True)
    plt.savefig(os.path.join(out_dir, 'dotplot_Ccr2.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_Ccr2.pdf'))

with rc_context({'figure.figsize': (5, 5), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.dotplot(adata, var_names=['Ccr2'], groupby=['cell_type','group'], dendrogram=False, show=False,swap_axes=False)
    plt.savefig(os.path.join(out_dir, 'dotplot_Ccr2_group.png'))
    plt.savefig(os.path.join(out_dir, 'dotplot_Ccr2_group.pdf'))
