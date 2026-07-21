import pandas as pd
import infercnvpy as cnv
import os
import matplotlib.pyplot as plt
import scanpy as sc
def save_fig(plt, out_dir = os.getcwd(), file_name=None, fig_size={'w':5, 'h':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

#work_dir
out_path='/mnt/NexData_file/home/lst/project/AE2606294867/inferCNV'
#annodata_path
data_path="/mnt/NexData_file/home/lst/project/AE2606294867/infercnv_use.h5ad"
#sub par
celltype_col = "cell_type"
#species
species = "mouse"

if species == "human":
    gene_pos = "/mnt/NexData_file/home/lyh/Script/scRNA/inferCNV/GRch38_2020A.txt"
else:
    gene_pos = "/mnt/NexData_file/home/lyh/Script/scRNA/inferCNV/mm10_2020A.txt"

adata = sc.read_h5ad(data_path)
adata = adata.raw.to_adata()
gene_pos = pd.read_csv(gene_pos,sep="\t",index_col=0,header=None)
gene_pos.columns = ["chromosome", "start", "end"]
gene_pos.index.name = None
for i in ["chromosome", "start", "end"]:
    adata.var[i] = gene_pos[i]

reference_cats = ['T Cell']
cnv.tl.infercnv(
    adata,
    reference_key=celltype_col,
    reference_cat= reference_cats,
    window_size=250
)

outdir = '/mnt/NexData_file/home/lst/project/AE2606294867/inferCNV'
os.makedirs(outdir,exist_ok=True)
cnv.pl.chromosome_heatmap(adata, groupby=celltype_col,show=False)
save_fig(plt,out_dir=outdir,file_name='chromosome_heatmap',fig_size={'w':15,'h':10})

cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
cnv.tl.umap(adata)
adata.obs['cellname'] = adata.obs.index
cnv.tl.cnv_score(adata,groupby='cellname')

sc.pl.embedding(
    adata,
    color=["cnv_score",'cell_type','cnv_leiden'],
    show=False,
    size= 20,
    basis="X_cnv_umap",
    wspace = 0.5
)
save_fig(plt,out_dir=outdir,file_name='CNV_score',fig_size={'w':25,'h':6})

adata.obs[["cnv_score"]].to_csv("/mnt/NexData_file/home/lst/project/AE2606294867/inferCNV/cnv_score.csv")