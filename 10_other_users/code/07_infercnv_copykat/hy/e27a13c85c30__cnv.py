import pandas as pd
import infercnvpy as cnv
import os
import matplotlib.pyplot as plt
import scanpy as sc

adata = sc.read_h5ad("merged_data.h5ad")
cc = adata.obs["celltype"].value_counts()
print(cc)

adata = adata.raw.to_adata()

def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()


os.makedirs('./CNV')
out_dir='./CNV'

#human
gene_pos = pd.read_csv('/lvdata/wzb/pipline/InferCNV/inferCNV_R/GRch38_2020A.txt',sep="\t",index_col=0,header=None)
#mouse
#gene_pos = pd.read_csv('/lvdata/wzb/pipline/InferCNV/inferCNV_R/mm10_2020A.txt',sep="\t",index_col=0,header=None)

gene_pos.columns = ["chromosome", "start", "end"]
gene_pos.index.name = None
for i in ["chromosome", "start", "end"]:
    adata.var[i] = gene_pos[i]
    

cnv.tl.infercnv(
    adata,
    reference_key='celltype',
    reference_cat= ['TNK'],
    window_size=250
)

#cnv.pl.chromosome_heatmap(adata, groupby="celltype",show=False, save="heatmap_cluster.png")
#cnv.pl.chromosome_heatmap_summary(adata, groupby="celltype",show=False, save="heatmap_cluster_s.png")
cnv.pl.chromosome_heatmap(adata, groupby="celltype",show=False,save="heatmap_celltype.png")
cnv.pl.chromosome_heatmap_summary(adata, groupby="celltype",show=False,save="heatmap_celltype_s.png")

cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
#cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True)
#save_fig(plt,out_dir,file_name='chromosome_heatmap_cnvleiden',fig_size={'w':10,'h':7})

cnv.tl.umap(adata)
adata.obs['cellname'] = adata.obs.index
cnv.tl.cnv_score(adata,groupby='cellname')
#cnv.tl.cnv_score(adata)

sc.pl.umap(
    adata,
    color=["cnv_leiden",'celltype','cnv_score'],
    show=False,
    size= 20,
    wspace = 0.5
)
save_fig(plt,out_dir,file_name='cnv_score_celltype',fig_size={'w':30,'h':7})

sc.pl.umap(
    adata,
    color=["cnv_leiden",'celltype','cnv_score'],
    legend_loc="on data",
    show=False,
    legend_fontoutline=2,
    size= 20,
    wspace = 0.5
)
save_fig(plt,out_dir,file_name='cnv_score_celltype2',fig_size={'w':30,'h':7})

cnv.pl.umap(
    adata,
    color=["cnv_leiden",'celltype','cnv_score'],
    show=False,
    size= 20,
    wspace = 0.5
)
save_fig(plt,out_dir,file_name='cnv_score_cnvleiden',fig_size={'w':30,'h':7})

cnv.pl.umap(
    adata,
    color=["cnv_leiden",'celltype','cnv_score'],
    legend_loc="on data",
    show=False,
    legend_fontoutline=2,
    size= 20,
    wspace = 0.5
)
save_fig(plt,out_dir,file_name='cnv_score_cnvleiden2',fig_size={'w':30,'h':7})

#注释
adata.obs["cnv_status"] = "normal"
adata.obs.loc[adata.obs["cnv_leiden"].isin(["6", "8", "15", "19"]), "cnv_status"] = ("tumor")
sc.pl.umap(adata, color=["cnv_status"],title='leiden', add_outline=False,size= 20,
           frameon=False, legend_fontsize=10, legend_fontoutline=2)
save_fig(plt,out_dir,file_name='umap_cnv',fig_size={'w':5,'h':5})

#meta
meta=pd.DataFrame(data=adata.obs)
meta=meta[['cell_id','celltype','cnv_score']]
meta.to_csv('./CNV/metadata_a.csv')
