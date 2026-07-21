import gseapy as gp
import os
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt
import scanpy.external as sce

class recluster_new:
    def __init__(self,adata=None,species = None,neighbors=20,resolution=0.5,n_pcs = 10,method='leiden',batch_effect = 0,output=os.getcwd()):
        '''
        species  def is None,can use human or mouse
        neighbors 30
        resulution 0.5
        method leiden or louvain
        batch_effect 0
        output output diectory
        '''
        self.adata = adata
        self.species =species
        self.neighbors = neighbors
        self.resolution = resolution
        self.output = output
        self.method = method
        self.batch  = batch_effect
        self.n_pcs = n_pcs
    
    def sub_cluster(self):
        
        neighbors = self.neighbors
        resolution = self.resolution
        output = self.output
        method = self.method
        species = self.species
        batch_effect = self.batch
        n_pcs = self.n_pcs
        self.adata = self.adata.raw.to_adata()
        adata_tmp = self.adata.copy()
        # sc.pp.normalize_total(self.adata, target_sum=1e4)
        # sc.pp.log1p(self.adata)
        # self.adata.uns['log1p']['base'] = None
        self.adata.uns['log1p']={'base':None}
        sc.pp.highly_variable_genes(self.adata,n_top_genes=2000)
        self.adata.raw = adata_tmp
        self.adata = self.adata[:,self.adata.var.highly_variable]
        print("neighbors:"+str(neighbors))
        print("resolution:"+str(resolution))
        print("output"+output)
        # sc.pp.log1p(self.adata)
        # sc.pp.highly_variable_genes(self.adata, n_top_genes=2000)
        sc.pp.scale(self.adata,max_value=10)
        sc.tl.pca(self.adata, svd_solver='arpack')
        sc.pl.pca(self.adata, color='sample', show=False)
        save_fig(plt, output, file_name='pca_sample')
        sc.pl.pca_variance_ratio(self.adata, n_pcs=30, show=False)
        save_fig(plt, output, file_name='pca_variance')
        # Umap Tsne
        if batch_effect == 2:
            self.adata.obsm["X_pca_raw"] = self.adata.obsm["X_pca"].copy()
            sce.pp.harmony_integrate(self.adata, key ='sample',adjusted_basis="X_pca")
        sc.pp.neighbors(self.adata, n_neighbors=neighbors, n_pcs=n_pcs)
        if batch_effect == 1:
            self.adata.obsm["X_pca_raw"] = self.adata.obsm["X_pca"].copy()
            sce.pp.bbknn(self.adata, batch_key='sample', n_pcs=n_pcs)  # running bbknn 1.3.6
        sc.tl.umap(self.adata)
        if method == 'leiden':
            sc.tl.leiden(self.adata, resolution=resolution)
        else:
            sc.tl.louvain(self.adata, resolution=resolution)
        # plot cluster
        sc.pl.umap(self.adata, color=method, show=False)
        save_fig(plt, output, file_name='sub_cluster_umap')
        sc.pl.umap(self.adata, color=method,legend_loc ='on data', show=False)
        save_fig(plt, output, file_name='sub_cluster_umap_withlabel')
#        sc.tl.tsne(self.adata)
#        sc.pl.tsne(self.adata, color=method, size=25, show=False)
#        save_fig(plt, output, file_name='sub_cluster_tsne')
#        sc.pl.tsne(self.adata, color=method,size=25, legend_loc ='on data', show=False)
#        save_fig(plt, output, file_name='sub_cluster_tsne_withlabel')
        # cal marker gene
        # sc.pp.log1p(self.adata)
        print(self.adata.obs.columns.tolist())
        sc.tl.rank_genes_groups(self.adata,
                                groupby=f"{method}",
                                method='wilcoxon',use_raw=True,
                                pts=True)
        sc.pl.rank_genes_groups(self.adata, n_genes=10, sharey=False,show=False)
        save_fig(plt, output, file_name='rank_genes_groups', fig_size={'width':20, 'height':12})
        self.adata.write_h5ad(os.path.join(output,"result.h5ad"))
        # enrich
        # find_group_marker(adata=self.adata,species=species,output=output)
        from matplotlib.pyplot import rc_context
        sc.set_figure_params(dpi=80, color_map = 'viridis_r')

        # group = self.adata.obs['group'].drop_duplicates().tolist()
        # if len(group) > 1:
        #     axs = tuple(map(lambda x: f'ax{x+1}', range(len(group))))
        #     fig, axs = plt.subplots(nrows=1, ncols=len(group), figsize=(len(group)*5, 6))
        #     # keys = self.adata.obs['cell_type'].cat.categories.to_list()
        #     # values = list(self.adata.uns['cell_type_colors'])
        #     # color_dict = dict(zip(keys,values))
        #     for i in range(len(group)):
        #         if i == len(group) - 1:
        #             legend = 'right margin'
        #         else:
        #             legend = None
        #         sc.pl.umap(self.adata[self.adata.obs["group"]==group[i]], size = 20,color=method, legend_loc=legend,
        #                 #  cmap = color_dict,
        #                    title=group[i], ax=axs[i], show=False)
        #     plt.savefig(os.path.join(output, 'umap_groups.png'), facecolor='white', bbox_inches="tight", dpi=300)
        #     plt.savefig(os.path.join(output, 'umap_groups.pdf'), facecolor='white', bbox_inches="tight", dpi=300)
        group = self.adata.obs['group'].drop_duplicates().tolist()
        if len(group) > 1:
            axs = tuple(map(lambda x: f'ax{x+1}', range(len(group))))
            fig, axs = plt.subplots(nrows=1, ncols=len(group), figsize=(len(group)*6, 6))
            keys = self.adata.obs[method].cat.categories.to_list()
            values = list(self.adata.uns[method+'_colors'])
            color_dict = dict(zip(keys,values))
            for i in range(len(group)):
                tmp = self.adata[self.adata.obs["group"]==group[i]]
                sc.pl.umap(tmp, size = 20,color=method, 
                        legend_loc='on data',palette = color_dict, title=group[i], ax=axs[i], show=False)
        save_fig(plt,out_dir=output,file_name='umap_split_groups', fig_size={'w':len(group)*6, 'h':6})
        
        samples = self.adata.obs['sample'].drop_duplicates().tolist()
        print(samples)
        if len(samples) > 1:
            axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
            fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*6, 6))
            keys = self.adata.obs[method].cat.categories.to_list()
            values = list(self.adata.uns[method+'_colors'])
            color_dict = dict(zip(keys,values))
            for i in range(len(samples)):
                tmp = self.adata[self.adata.obs["sample"]==samples[i]]
                sc.pl.umap(tmp, size = 20,color=method, 
                        legend_loc='on data',palette = color_dict, title=samples[i], ax=axs[i], show=False)
        save_fig(plt,out_dir=output,file_name='umap_split_sample', fig_size={'w':len(samples)*6, 'h':6})

        sc.pl.umap(self.adata, color='sample', title='', show=False)
        save_fig(plt, output, file_name='umap_sample_dimplot')
        sc.pl.umap(self.adata, color='group', title='', show=False)
        save_fig(plt, output, file_name='umap_group_dimplot')

        return self.adata

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