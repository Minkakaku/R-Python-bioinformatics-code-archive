#!/usr/bin/env python

import os
import argparse
import scanpy as sc
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("--input", required=True)
parser.add_argument("--outdir", required=True)

args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)
def save_fig(plt, out_dir, file_name, fig_size={'width':5, 'height':5}):
    plt.tight_layout()
    try:
        plt.gcf().set_size_inches(fig_size['width'], fig_size['height'])
    except:
        pass
    plt.savefig(os.path.join(out_dir, f'{file_name}.png'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.savefig(os.path.join(out_dir, f'{file_name}.pdf'), bbox_inches = 'tight', pad_inches = 0.1)
    plt.clf()

print("Loading:", args.input)

adata = sc.read_h5ad(args.input)

print(
    f"Cells: {adata.n_obs}, Genes: {adata.n_vars}"
)

print("Running TSNE...")

# 如果有Harmony结果优先用Harmony
if "X_pca_harmony" in adata.obsm.keys():

    sc.tl.tsne(
        adata,
        use_rep="X_pca_harmony",
        random_state=1234
    )

# 否则用PCA
elif "X_pca" in adata.obsm.keys():

    sc.tl.tsne(
        adata,
        use_rep="X_pca",
        random_state=1234
    )

else:

    sc.tl.tsne(
        adata,
        random_state=1234
    )

adata.write(
    os.path.join(
        args.outdir,
        "adata_tsne.h5ad"
    )
)

sc.pl.tsne(adata, 
           color='cluster',
           title='Cluster',
           legend_loc='on data',
           legend_fontsize=10,
           legend_fontoutline=2, 
           show=False)
save_fig(plt, args.outdir, file_name='tsne_dimplot')
sc.pl.tsne(adata, color='sample', title='', show=False)
save_fig(plt, args.outdir, file_name='TSNE_sample_dimplot')
sc.pl.tsne(adata, color='group', title='', show=False)
save_fig(plt, args.outdir, file_name='TSNE_group_dimplot')

groups = adata.obs['group'].drop_duplicates().tolist()
if len(groups) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(groups))))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups)*5, 5))
    for i in range(len(groups)):
        if i == len(groups) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.tsne(adata[adata.obs["group"]==groups[i]], color="cluster", title=groups[i], ax=axs[i], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
    save_fig(plt, args.outdir, file_name='TSNE_split_group_cluster', fig_size={'width':len(groups)*5, 'height':5})

samples = adata.obs['sample'].drop_duplicates().tolist()
if len(samples) > 1:
    axs = tuple(map(lambda x: f'ax{x+1}', range(len(samples))))
    fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples)*5, 5))
    for i in range(len(samples)):
        if i == len(samples) - 1:
            legend = 'right margin'
        else:
            legend = None
        sc.pl.tsne(adata[adata.obs["sample"]==samples[i]], color="cluster", title=samples[i], ax=axs[i], wspace=0.1, legend_loc='on data', legend_fontsize=12, legend_fontweight='normal', legend_fontoutline=2, show=False)
    save_fig(plt, args.outdir, file_name='TSNE_split_sample_cluster', fig_size={'width':len(samples)*5, 'height':5})

print("Done")