import pandas as pd
import seaborn as sns
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

path = pd.read_csv("/data/NAS02/sh/Projects/FWSC2024-192/05.Mission/Enrich/enrich_select.csv",sep="\t")
adata = 

name_path = path.columns.to_list()
name_path
for i in name_path:
    sc.tl.score_genes(adata,path[i].dropna(how = 'all').to_list(),score_name=i)
    outdir = '/lvdata/wzb/scRNA/FW2024-318/res'
    max_value = adata.obs[i].max()
    sc.pl.umap(adata[adata.obs['group']=='Patient2-1'],color=i,show=False,cmap='RdBu_r',vmax=max_value)
    save_fig(plt,out_dir=outdir,file_name=f"Patient2-1_{i}_score_umap",fig_size={'w':5,'h':5})
    sc.pl.umap(adata[adata.obs['group']=='Patient2-2'],color=i,show=False,cmap='RdBu_r',vmax=max_value)
    save_fig(plt,out_dir=outdir,file_name=f"Patient2-2_{i}_score_umap",fig_size={'w':5,'h':5})
sc.pl.dotplot(adata,name_path,groupby=['group'],show=False,swap_axes=True)
save_fig(plt,out_dir=outdir,file_name="path_score_dotplot",fig_size={'w':10,'h':5})