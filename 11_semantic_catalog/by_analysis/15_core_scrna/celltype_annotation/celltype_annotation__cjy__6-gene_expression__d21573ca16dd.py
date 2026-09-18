import os
import shutil
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
def save_fig(plt, out_dir, file_name, fig_size={'width':5, 'height':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['width'], fig_size['height'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

def visua_gene_express(adata, genes, out_dir):
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(adata, color=genes, s=15, ncols=3, cmap='coolwarm', use_raw=True, show=False)
        plt.savefig(os.path.join(out_dir, 'marker_list_umap.png'), facecolor='white', bbox_inches="tight", dpi=200)
        plt.savefig(os.path.join(out_dir, 'marker_list_umap.pdf'))

    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.dotplot(adata, genes,groupby='cluster',cmap='coolwarm', use_raw=True, show=False,swap_axes=True)
        plt.savefig(os.path.join(out_dir, 'marker_list_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=200)
        plt.savefig(os.path.join(out_dir, 'marker_list_dotplot.pdf'))

    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.dotplot(adata, genes,groupby='cell_type',cmap='coolwarm', use_raw=True, show=False)
        plt.savefig(os.path.join(out_dir, 'anno_dotplot.png'), facecolor='white', bbox_inches="tight", dpi=200)

    if len(genes) > 1:
        fig, axs = plt.subplots(len(genes), 1, figsize=(10, len(genes)*4), sharey=True)
        for i in range(len(genes)):
            sc.pl.violin(adata, genes[i], groupby="cluster", stripplot=False, xlabel='', inner='box', ax=axs[i], show=False)
    else:
        fig, ax = plt.subplots(figsize=(4, 4))
        sc.pl.violin(adata, genes[0], groupby="cluster", stripplot=False, xlabel='', inner='box', ax=ax, show=False)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="gene expression")
    parser.add_argument("--input_file", required=True, help="输入文件 h5ad")
    parser.add_argument("--out_dir", required=True, help="输出文件夹路径")
    parser.add_argument("--marker_file", required=True, help="输入marker文件")
    args = parser.parse_args()
    gene_express_out_dir = os.path.join(args.out_dir, '06-Cell_type_annotation')
    if not os.path.exists(marker_gene_out_dir):
        os.makedirs(marker_gene_out_dir)
    adata = sc.read_h5ad(args.input_file)
    df = pd.read_csv("args.marker_file", sep="\t")

    genes = (df["markers"]
    .dropna()
    .str.split(",")   # 按逗号拆分
    .explode()        # 拉平成一列
    .str.strip()      # 去空格
    .tolist()
    )
    genes_unique = list(dict.fromkeys(genes))
    pd.DataFrame({"markers": genes_unique}).to_csv(os.path.join(gene_express_out_dir,"markers.xls"),sep='\t', index=False)
    visua_gene_express(adata, genes_unique, gene_express_out_dir)