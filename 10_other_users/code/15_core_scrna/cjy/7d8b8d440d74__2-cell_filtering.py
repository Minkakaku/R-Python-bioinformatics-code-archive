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

def cell_filtering(adata, min_genes, max_genes, min_cells, total_counts,pct_counts_mt, pct_counts_ribo, expected_doublet_rate,doublet_score , out_dir):    
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[(adata.obs.n_genes_by_counts < max_genes)]
    adata = adata[adata.obs['pct_counts_mt'] < pct_counts_mt, :]
    adata = adata[adata.obs['pct_counts_ribo'] < pct_counts_ribo, :]
    adata = adata[adata.obs['total_counts'] < total_counts, :]
    
    sc.external.pp.scrublet(adata, expected_doublet_rate=expected_doublet_rate / 100, threshold=doublet_score, batch_key='sample', random_state=1234)
    samples = sorted(adata.obs['sample'].drop_duplicates().tolist())
    if len(samples) > 1:
        sc.external.pl.scrublet_score_distribution(adata,show=False)
        save_fig(plt, out_dir, file_name='scrublet_score_distribution', fig_size={'width':12, 'height':2.5*len(samples)})  
    adata = adata[adata.obs['predicted_doublet']==False, :]

    pd.DataFrame(adata.obs).to_csv(os.path.join(out_dir, 'all_cell_features_filtered.xls'), sep='\t')
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby = 'sample', rotation=90, multi_panel=True, show=False)
    save_fig(plt, out_dir, file_name='filtered_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})  
    
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
    save_fig(fig, out_dir, file_name='filtered_nCount_nFeature_scatter')      
    
    for s in adata.obs['sample'].drop_duplicates().tolist():
        os.makedirs(os.path.join(out_dir, s), exist_ok=True)
        sc.pl.violin(adata[adata.obs['sample']==s], ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo'], jitter=0, groupby = 'sample', rotation= 90, multi_panel=True, show=False)
        save_fig(plt, os.path.join(out_dir, s), file_name='filtered_nCount_nFeature_voilon', fig_size={'width':25, 'height':10})
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
        sc.pl.scatter(adata[adata.obs['sample']==s], x='total_counts', y='pct_counts_mt', color="sample", legend_loc='none', ax=ax1, show=False)
        sc.pl.scatter(adata[adata.obs['sample']==s], x='total_counts', y='n_genes_by_counts', color="sample", ax=ax2, show=False)
        save_fig(fig, os.path.join(out_dir, s), file_name='filtered_nCount_nFeature_scatter')       
        pd.DataFrame(adata[adata.obs['sample']==s].obs).to_csv(os.path.join(out_dir, s, 'cell_features_filtered.xls'), sep='\t')
    adata.write(os.path.join(out_dir, 'adata_counts.h5ad'))

    return adata

# if __name__ == "__main__":
#     input_file = snakemake.input[0]
#     out_dir = snakemake.params["out_dir"]
#     Cell_filter_out_dir = os.path.join(out_dir, '03-Cell_filter')
#     if not os.path.exists(Cell_filter_out_dir):
#         os.makedirs(Cell_filter_out_dir)
#     else:
#         os.system(f'rm -rf {Cell_filter_out_dir}')
#         os.makedirs(Cell_filter_out_dir)
#     params = snakemake.params

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Single-cell filtering script")
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--min_genes", type=int, required=True)
    parser.add_argument("--max_genes", type=int, required=True)
    parser.add_argument("--min_cells", type=int, required=True)
    parser.add_argument("--total_counts", type=int, required=True)
    parser.add_argument("--pct_counts_mt", type=float, required=True)
    parser.add_argument("--pct_counts_ribo", type=float, required=True)
    parser.add_argument("--expected_doublet_rate", type=float, required=True)
    parser.add_argument("--doublet_score", type=float, required=True)
    args = parser.parse_args()

    Cell_filter_out_dir = os.path.join(args.out_dir, '03-Cell_filter')
    if not os.path.exists(Cell_filter_out_dir):
        os.makedirs(Cell_filter_out_dir)
    else:
        os.system(f'rm -rf {Cell_filter_out_dir}')
        os.makedirs(Cell_filter_out_dir)
    adata = sc.read_h5ad(args.input_file)
    cell_filtering(
        adata,
        args.min_genes, args.max_genes, args.min_cells,
        args.total_counts, args.pct_counts_mt, args.pct_counts_ribo,
        args.expected_doublet_rate, args.doublet_score, Cell_filter_out_dir
    )

    csv_path = os.path.join(args.out_dir, "Parameter_dict.csv")
    if os.path.exists(csv_path):
        param_df = pd.read_csv(csv_path)
        Parameter_dict = dict(zip(param_df["Key"], param_df["Value"]))
    else:
        Parameter_dict = {}

    VALUE_KEYS = [
    "min_genes",
    "max_genes",
    "min_cells",
    "total_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "expected_doublet_rate",
    "doublet_score"
    ]

    run_info = {k: getattr(args, k) for k in VALUE_KEYS}
    Parameter_dict.update(run_info)

    pd.DataFrame(
        Parameter_dict.items(),
        columns=["Key", "Value"]
        ).to_csv(csv_path, index=False)

