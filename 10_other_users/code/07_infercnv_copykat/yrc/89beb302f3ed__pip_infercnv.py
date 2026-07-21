import os
import glob
import subprocess
import argparse
import pandas as pd
import infercnvpy as cnv
import matplotlib.pyplot as plt
import scanpy as sc


def save_fig(plt, fig_out, file_name=None, fig_size={'w':5, 'h':5}):

    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['w'], fig_size['h'])
    except:
        pass
    plt.savefig(os.path.join(fig_out, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(fig_out, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()



def run_cnv(h5ad, output, species, reference_celltype):

    # 设置基因组路径
    if species == 'human':
        gene_pos = pd.read_csv('/mnt/NexData_file/share/scRNA-tools/InferCNV/inferCNV_R/GRch38_2020A.txt',sep="\t",index_col=0,header=None)
    elif species == 'mouse':
        gene_pos = pd.read_csv('/mnt/NexData_file/share/scRNA-tools/InferCNV/inferCNV_R/mm10_2020A.txt',sep="\t",index_col=0,header=None)
    else:
        raise ValueError("Invalid species. Must be 'human' or 'mouse'.")
    
    # 确保输出目录存在
    os.makedirs(os.path.join(output, 'CNV'), exist_ok=True)
    out_dir = os.path.join(output, 'CNV')

    adata = sc.read_h5ad(h5ad)
    adata = adata.raw.to_adata()
    #adata.obs['cell_type'] = adata.obs['celltype']
    gene_pos.columns = ["chromosome", "start", "end"]
    gene_pos.index.name = None
    for i in ["chromosome", "start", "end"]:
        adata.var[i] = gene_pos[i]


    reference_celltype_list = reference_celltype.split(',')
    
    cnv.tl.infercnv(
        adata,
        reference_key='cell_type',
        reference_cat= reference_celltype_list,
        window_size=250
    )

    cnv.pl.chromosome_heatmap(adata, groupby="cell_type",show=False)
    save_fig(plt,out_dir,file_name='chromosome_heatmap',fig_size={'w':10,'h':7})

    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    cnv.tl.leiden(adata)
    cnv.tl.umap(adata)

    
    ####by cell####
    adata.obs['cellname'] = adata.obs.index
    cnv.tl.cnv_score(adata,groupby='cellname')
    adata.obs['cnv_score_cell'] = adata.obs.cnv_score

    #violinplot
    sc.pl.violin(adata, ['cnv_score'], groupby='cnv_leiden')
    save_fig(plt,out_dir,file_name='cnv_vlnplot',fig_size={'w':30,'h':7})

    sc.pl.umap(
        adata,
        color=["cnv_leiden",'cell_type','cnv_score'],
        show=False,
        size= 20,
        wspace = 0.5
    )
    save_fig(plt,out_dir,file_name='cnv_celltype_bycell',fig_size={'w':30,'h':7})

    sc.pl.umap(
        adata,
        color=["cnv_leiden",'cell_type','cnv_score'],
        legend_loc="on data",
        show=False,
        legend_fontoutline=2,
        size= 20,
        wspace = 0.5
    )
    save_fig(plt,out_dir,file_name='cnv_celltype_bycell2',fig_size={'w':30,'h':7})

    cnv.pl.umap(
        adata,
        color=["cnv_leiden",'cell_type','cnv_score'],
        show=False,
        size= 20,
        wspace = 0.5
    )
    save_fig(plt,out_dir,file_name='cnv_cnvleiden_bycell',fig_size={'w':30,'h':7})

    cnv.pl.umap(
        adata,
        color=["cnv_leiden",'cell_type','cnv_score'],
        legend_loc="on data",
        show=False,
        legend_fontoutline=2,
        size= 20,
        wspace = 0.5
    )
    save_fig(plt,out_dir,file_name='cnv_cnvleiden_bycell2',fig_size={'w':30,'h':7})

    ####by cluster####
    cnv.tl.cnv_score(adata)
    adata.obs['cnv_score_cluster'] = adata.obs.cnv_score

    sc.pl.umap(
        adata,
        color=["cnv_leiden",'cell_type','cnv_score'],
        show=False,
        size= 20,
        wspace = 0.5
    )
    save_fig(plt,out_dir,file_name='cnv_celltype_bycluster',fig_size={'w':30,'h':7})

    sc.pl.umap(
        adata,
        color=["cnv_leiden",'cell_type','cnv_score'],
        legend_loc="on data",
        show=False,
        legend_fontoutline=2,
        size= 20,
        wspace = 0.5
    )
    save_fig(plt,out_dir,file_name='cnv_celltype_bycluster2',fig_size={'w':30,'h':7})

    cnv.pl.umap(
        adata,
        color=["cnv_leiden",'cell_type','cnv_score'],
        show=False,
        size= 20,
        wspace = 0.5
    )
    save_fig(plt,out_dir,file_name='cnv_cnvleiden_bycluster',fig_size={'w':30,'h':7})

    cnv.pl.umap(
        adata,
        color=["cnv_leiden",'cell_type','cnv_score'],
        legend_loc="on data",
        show=False,
        legend_fontoutline=2,
        size= 20,
        wspace = 0.5
    )
    save_fig(plt,out_dir,file_name='cnv_cnvleiden_bycluster2',fig_size={'w':30,'h':7})

    #meta
    meta=pd.DataFrame(data=adata.obs)
    meta=meta[['cell_type','cnv_leiden','cnv_score_cell','cnv_score_cluster']]
    meta.to_csv(os.path.join(out_dir, f'metadata.csv'))





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process inferCNVpy for h5ad.")
    parser.add_argument("--h5ad", required=True, help="h5ad object")
    parser.add_argument("--output", required=True, help="output directory")
    parser.add_argument("--species", required=True, choices=["human", "mouse"], help="Species")
    parser.add_argument("--reference_celltype", required=True, help="reference cell types")

    args = parser.parse_args()
    
    # 确保输出路径存在
    os.makedirs(args.output, exist_ok=True)
    
    run_cnv(args.h5ad, args.output, args.species, args.reference_celltype)
