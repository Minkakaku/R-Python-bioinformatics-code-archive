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

# 对Proximal tubule cells注释细分
adata_raw = sc.read_h5ad("/mnt/NexData_file/scRNA/scVision_output/FWSC20253251/06-Cell_type_annotation/result.h5ad")
# 'Proximal tubule cells', 'Macrophages', 'Distal tubule cells','Proliferation', 'Endothelial cells', 'T cells', 'Neutrophils','Dendritic cells', 'Fibroblasts', 'Erythrocyte', 'Principal cells','B cells', 'Smooth muscle cells'
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/02.Mission/Proximal tubule cells/"
adata_sub = adata_raw[adata_raw.obs.cell_type=="Proximal tubule cells", :]
adata_recluster = recluster(adata=adata_sub,species='mouse',neighbors=20,resolution=0.5,n_pcs = 10,output=out_dir,batch_effect = 1)
adata = adata_recluster.sub_cluster()

sc.pl.rank_genes_groups_dotplot(
    adata,
    groupby='leiden',
    key='rank_genes_groups',
    n_genes=10,
    show=False)
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.png'),bbox_inches='tight')
plt.savefig(os.path.join(out_dir, 'top10_diff_genes_dotplot.pdf'),bbox_inches='tight')
plt.close()

sc.pl.umap(adata, color=['Slc5a2','Slc34a1','Slco4c1','Slc5a12','Slc22a8','Cyp2e1','Slc7a13'], cmap='coolwarm',show=False)
plt.savefig(os.path.join(out_dir, 'umap_healthy.png'))
plt.savefig(os.path.join(out_dir, 'umap_healthy.pdf'))
sc.pl.dotplot(adata, var_names=['Slc5a2','Slc34a1','Slco4c1','Slc5a12','Slc22a8','Cyp2e1','Slc7a13'], groupby='leiden', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_healthy.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_healthy.pdf'), facecolor='white', bbox_inches="tight")

sc.pl.umap(adata, color=['Fgg','Mt2','Krt20','Havcr1','Clcf1','Jag1','Gm20481','Vcam1','C3','Cd74'], cmap='coolwarm',show=False)
plt.savefig(os.path.join(out_dir, 'umap_injured.png'))
plt.savefig(os.path.join(out_dir, 'umap_injured.pdf'))
sc.pl.dotplot(adata, var_names=['Fgg','Mt2','Krt20','Havcr1','Clcf1','Jag1','Gm20481','Vcam1','C3','Cd74'], groupby='leiden', dendrogram=False, show=False, standard_scale='var')
plt.savefig(os.path.join(out_dir, 'dotplot_injured.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_injured.pdf'), facecolor='white', bbox_inches="tight")

anno_dict={
    '0':'Healthy_S3',
    '1':'Failed_repair_PT',
    '2':'Injured_S1',
    '3':'Injured_S2',
    '4':'Injured_S3',
    '5':'Injured_S2',
    '6':'Healthy_S2',
    '7':'Injured_S3',
    '8':'Injured_S2',
    '9':'Healthy_S1'
}
adata.obs['cell_type'] = [anno_dict[clust] for clust in adata.obs['leiden']]
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')

cell_type_mapping = adata.obs[['cell_type']].to_dict()['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs['cell_type']
adata_raw.obs['sub_cell_type'] = adata_raw.obs.apply(lambda row: cell_type_mapping.get(row.name, row['cell_type']), axis=1)
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata_raw, color='sub_cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))

with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cell_type', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_celltype.pdf'))


# groups = adata.obs['group'].drop_duplicates().tolist()
groups = ['H0', 'H2', 'H8', 'H24', 'H72']
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

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
adata.obs['group'] = adata.obs['group'].cat.reorder_categories(groups, ordered=True)
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')
adata.write(os.path.join(out_dir, 'PTC.h5ad'))

# 映射所有细分
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/02.Mission/mapping_raw/"
sc.pl.umap(adata_raw, color='sub_cell_type', legend_loc='right margin', show=False)
plt.savefig(os.path.join(out_dir, 'umap_mapping.png'), facecolor='white', bbox_inches="tight", dpi=300)
plt.savefig(os.path.join(out_dir, 'umap_mapping.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
adata_raw.obs.to_csv(os.path.join(out_dir,"all_cells_features.xls"), sep='\t')
adata_raw.write(os.path.join(out_dir, 'Main_mapping.h5ad'))



# 提取ARG2表达
out_dir="/mnt/NexData_file/home/sh/Projects/FWSC20253251/02.Mission/ARG2/"
adata_sub = adata_raw[adata_raw.obs.cell_type.isin(["Macrophages"]), :]
adata_sub.obs['group'] = adata_sub.obs['group'].cat.reorder_categories(groups, ordered=True)
adata = adata_sub
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='cluster', title='Cell Type', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_cluster.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_cluster.pdf'))

sc.pl.umap(adata, color=['Arg2'], cmap='coolwarm',show=False)
plt.savefig(os.path.join(out_dir, 'umap_Arg2.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_Arg2.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.umap(adata, color=['Arg2'], cmap='coolwarm',show=False,vmax=1)
plt.savefig(os.path.join(out_dir, 'umap_Arg2_fiter.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'umap_Arg2_fiter.pdf'), facecolor='white', bbox_inches="tight")
sc.pl.dotplot(adata, var_names=['Arg2'], groupby='cluster', dendrogram=False, show=False)
plt.savefig(os.path.join(out_dir, 'dotplot_Arg2-cluster.png'), facecolor='white', bbox_inches="tight")
plt.savefig(os.path.join(out_dir, 'dotplot_Arg2-cluster.pdf'), facecolor='white', bbox_inches="tight")


adata = adata_raw[adata_raw.obs.cell_type.isin(["Macrophages"]), :]
adata=adata.raw.to_adata()
expr = adata[:, 'Arg2'].X.toarray().ravel()
adata.obs['Arg2_gene'] = pd.Categorical(
    np.where(expr >= 1, 'high', 'low'),
    categories=['high', 'low']          # 固定顺序，防止后续错位
)
adata.uns['Arg2_gene_colors'] = ['#1f77b4', 'lightgray']
with rc_context({'figure.figsize': (4, 4), 'savefig.facecolor': 'white', 'savefig.bbox': 'tight', 'savefig.dpi': 300}):
    sc.pl.umap(adata, color='Arg2_gene', title='Arg2_gene', show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Arg2_gene.png'))
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_Arg2_gene.pdf'))

groups = ['H0', 'H2', 'H8', 'H24', 'H72']
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.umap(adata[adata.obs["group"]==groups[i]], color="Arg2_gene", legend_loc=legend, title=groups[i], ax=axs[i], show=False)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_Arg2_gene.png'), facecolor='white', bbox_inches="tight", dpi=300)
    plt.savefig(os.path.join(out_dir, 'dimplot_umap_group_Arg2_gene.pdf'), facecolor='white', bbox_inches="tight", dpi=300)

sc.tl.rank_genes_groups(adata, groupby='Arg2_gene',method='wilcoxon',pts=True)
find_group_marker(adata=adata,species='mouse',output=out_dir)

adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep='\t')
adata.obs['group'] = adata.obs['group'].cat.reorder_categories(groups, ordered=True)
adata.write(os.path.join(out_dir, 'Arg2_gene.h5ad'))
adata.obs['cell_type'] = adata.obs['Arg2_gene'] 
adata.uns['cell_type_colors'] = adata.uns['Arg2_gene_colors'] 
plot_celltype_proportion(adata, out_dir=out_dir,type='group')
plot_celltype_proportion(adata, out_dir=out_dir,type='sample')


adata = sc.read("/mnt/NexData_file/home/sh/Projects/FWSC20253251/02.Mission/Proximal tubule cells/PTC.h5ad")
out_dir = "/mnt/NexData_file/home/sh/Projects/FWSC20253251/02.Mission/"
from pyecharts import options as opts
from pyecharts.charts import Sankey

group_order = ['H0', 'H2', 'H8', 'H24', 'H72']
cell_type_order = ['Healthy_S1', 'Healthy_S2', 'Healthy_S3', 'Injured_S1','Injured_S2', 'Injured_S3','Failed_repair_PT']

# 2. 只保留指定顺序的数据
df = adata.obs[["group", "cell_type"]].astype(str)
df = df[df["group"].isin(group_order) & df["cell_type"].isin(cell_type_order)]
df["group"] = pd.Categorical(df["group"], categories=group_order, ordered=True)
df["cell_type"] = pd.Categorical(df["cell_type"], categories=cell_type_order, ordered=True)

# 3. 节点：按指定顺序生成
nodes = [{"name": f"{g}_{c}"} for g in group_order for c in cell_type_order
         if f"{g}_{c}" in {f"{row.group}_{row.cell_type}" for row in df.itertuples()}]

# 4. 边：相邻组间
links = []
for g1, g2 in zip(group_order[:-1], group_order[1:]):
    src = df[df["group"] == g1].value_counts("cell_type").reindex(cell_type_order, fill_value=0)
    tgt = df[df["group"] == g2].value_counts("cell_type").reindex(cell_type_order, fill_value=0)
    for ct in cell_type_order:
        if src[ct] > 0 and tgt[ct] > 0:
            links.append({"source": f"{g1}_{ct}", "target": f"{g2}_{ct}", "value": int(min(src[ct], tgt[ct]))})

# 5. 绘图
Sankey().add(
    "",
    nodes,
    links,
    linestyle_opt=opts.LineStyleOpts(opacity=0.5, curve=0.3,color = 'source'),
    label_opts=opts.LabelOpts(position="left", font_size=10),node_gap = 20,
).set_global_opts(
    title_opts=opts.TitleOpts(title="Group-celltype flow"),
).render("flow.html")

# 3. 计算每个 group 内 cell_type 的比例
prop_df = (df.groupby(["group", "cell_type"], observed=True).size()
             .unstack(fill_value=0)
             .reindex(columns=cell_type_order, fill_value=0)
             .apply(lambda x: x / x.sum(), axis=1))

# 4. 堆积折线图
prop_df.plot.area(figsize=(8, 4), alpha=0.8)
for x in range(1, len(group_order)):
    plt.axvline(x - 0.5, color='white', lw=2)
plt.xlabel("Group")
plt.ylabel("Proportion")
plt.title("Stacked Area of Cell-Type Proportions across Groups")
plt.legend(title="Cell type", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("stacked_area.png", dpi=300)
plt.close()

# 4. 堆积柱状图
ax = prop_df.plot.bar(stacked=True, figsize=(10, 6), alpha=0.8, width=0.8)
# 为每个cell_type在组间添加连线
for i, cell_type in enumerate(cell_type_order):
    # 只对非零的cell_type绘制连线
    if prop_df[cell_type].sum() > 0:
        # 计算每个组的中心y坐标
        y_centers = prop_df.cumsum(axis=1)[cell_type] - prop_df[cell_type] / 2
        x_positions = np.arange(len(group_order))
        # 绘制组间连线
        ax.plot(x_positions, y_centers, color='black', lw=1.5, marker='o', markersize=4)
        # 添加填充区域（可选）
        ax.fill_between(x_positions, 
                       y_centers - prop_df[cell_type] / 2, 
                       y_centers + prop_df[cell_type] / 2,
                       color=ax.patches[i].get_facecolor(), 
                       alpha=0.3)

# 美化图表
plt.xlabel("Group", fontsize=12)
plt.ylabel("Proportion", fontsize=12)
plt.title("Stacked Bar Chart with Group Connectors", fontsize=14)
plt.legend(title="Cell type", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=0)
plt.tight_layout()
plt.savefig("stacked_bar_conn.png", dpi=300, bbox_inches='tight')
plt.show()
prop_df.to_csv(os.path.join(out_dir, "prop_df.xls"), sep='\t')