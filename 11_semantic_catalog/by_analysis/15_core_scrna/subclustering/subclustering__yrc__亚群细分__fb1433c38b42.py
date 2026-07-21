#设置路径
outdir = '/home/yrc/SC2507103606/9.3'
os.makedirs(outdir,exist_ok=True)

#输入对象
adata = sc.read_h5ad('/home/yrc/SC2507103606/recluster/adata_anno.h5ad')
sc.pl.umap(adata, color=['cell_type'], legend_loc='on data', add_outline=False,
           frameon=False, legend_fontsize=10, legend_fontoutline=2)

#细分函数
import random
random.seed(123)
adata_sub = adata[adata.obs.cell_type.isin(["Microglia","Myeloid"]), :]
os.makedirs(f"{outdir}/MYE/",exist_ok=True)
adata_sub = recluster(adata_sub,species='human',resolution=0.5,neighbors=20,n_pcs = 10,batch_effect = 2,output=f"{outdir}/MYE/")
adata_sub = adata_sub.sub_cluster()

#微调分辨率
sc.tl.leiden(adata_sub, resolution=1)
sc.pl.umap(adata_sub, color=['leiden'], show=False,size=25)
save_fig(plt,file_name='leiden_umap',out_dir=outdir,fig_size={'w':7,'h':5})

#初看分布
sc.pl.umap(adata_sub, color=['cell_type','leiden','sample'], show=False,size=25)
split = adata_sub.obs['sample'].unique().tolist()
n_split = len(split)

nrows = 1
ncols = 4
fig, ax = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False, figsize=(7 * ncols, 5 * nrows))
axes = ax.flatten()

for i in range(n_split):
     tmp = adata_sub[adata_sub.obs['sample'] == split[i], :]
     sc.pl.umap(tmp, color='cell_type', ax=axes[i], show=False)
     axes[i].set_title(f"{split[i]}: cell_type")

plt.subplots_adjust(wspace=0.6)
plt.show

#初染marker
sc.pl.umap(adata_sub, color=['LYZ','APOE','C1QA','C1QB','C1QC','HLA-DRA','VCAN','ACE','P2RY12','SALL1','ITGA4','THBD'], show=False,size=25)
save_fig(plt,file_name='marker_umap',out_dir=outdir,fig_size={'w':20,'h':15})

#marker映射
marker_dict ={'Microglia':['P2RY12','SALL1'],
              'Macrophage':['ITGA4','THBD']}
sc.pl.dotplot(adata_sub,var_names=marker_dict,groupby = 'cell_type',show=False)
save_fig(plt,file_name='marker_list_bubble',out_dir=outdir,fig_size={'w':10,'h':10})
unique_markers = list(set(gene for sublist in marker_dict.values() for gene in sublist))
sc.pl.umap(adata_sub, color=unique_markers,show=False)
save_fig(plt,file_name='marker_list_umap',out_dir=outdir,fig_size={'w':24,'h':5})
sc.pl.stacked_violin(adata_sub, marker_dict, groupby="cell_type", swap_axes=False, dendrogram=False,show=False)
save_fig(plt,file_name='marker_list_stack',out_dir=outdir)

#注释
anno_list = {'0': 'Microglia',
             '1': 'Microglia',
             '2': 'Microglia',
             '3': 'Microglia',
             '4': 'Macrophage',
             '5': 'Microglia',
             '6': 'Macrophage',
             '7': 'Macrophage',
             '8': 'Microglia',
             '9': 'Microglia',
             '10': 'Microglia',
             '11': 'Microglia',
             '12': 'Microglia',
             '13': 'Monocyte',
             '14': 'Macrophage',
}
adata_sub.obs['cell_type'] = [anno_list[clust] for clust in adata_sub.obs['leiden']]
sc.pl.umap(adata_sub,color='cell_type',show=False,size=20)
save_fig(plt,file_name='anno_umap',out_dir=outdir,fig_size={'w':6,'h':5})
adata_sub.write_h5ad(os.path.join(outdir,'adata_sub_anno.h5ad'))

#分组图
draw_group(adata_sub,out_dir = outdir,method = 'cell_type',size=20)

#密度图
adata_sub.obs.to_csv(os.path.join(outdir, 'all_cells_features.xls'),sep='\t')
sc.tl.embedding_density(adata_sub, basis='umap', groupby='group')
sc.pl.embedding_density(adata_sub, basis='umap', key='umap_density_group',show=False)
save_fig(plt,out_dir=outdir,file_name='group_density',fig_size={'w':20,'h':10})
sc.tl.embedding_density(adata_sub, basis='umap', groupby='sample')
sc.pl.embedding_density(adata_sub, basis='umap', key='umap_density_sample',show=False)
save_fig(plt,out_dir=outdir,file_name='sample_density',fig_size={'w':40,'h':10})

#统计
result = adata_sub.obs.groupby(['cell_type', 'sample']).size().reset_index(name='number').pivot(index='cell_type', columns='sample', values='number').reset_index()
result.columns.name = None
# result['Total'] = result.sum(axis=1)
# result['Percent'] = round((result['Total'] / sum(result['Total'])) * 100, 1)
result.to_csv(os.path.join(outdir, 'celltype_statistics.xls'), sep='\t', index=None)

#映射大群
# 备份原始注释
adata.obs['cell_type_original'] = adata.obs['cell_type'].copy()
# 将分类数据转换为字符串，然后赋值
adata.obs['cell_type'] = adata.obs['cell_type'].astype(str)
adata.obs.loc[adata_sub.obs_names, 'cell_type'] = adata_sub.obs['cell_type'].astype(str)
# 如果需要，可以再转换回分类数据
# adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')
sc.pl.umap(adata, color=['cell_type'], legend_loc='on data', add_outline=False,show=False,
           frameon=False, legend_fontsize=10, legend_fontoutline=2)
save_fig(plt,file_name='all_anno_umap',out_dir=outdir,fig_size={'w':5,'h':5})