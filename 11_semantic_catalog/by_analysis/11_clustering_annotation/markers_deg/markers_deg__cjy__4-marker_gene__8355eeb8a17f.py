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

def find_marker_gene(adata, method, out_dir, logfc=1, fdr=0.05, pct=0.25):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        os.makedirs(out_dir, exist_ok=True)
   
    sc.tl.rank_genes_groups(adata, 'cluster', method=method, pts=True)
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    save_fig(plt, out_dir, file_name='rank_genes_groups', fig_size={'width':20, 'height':12})
    
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    
    all_diff = []
    for c in groups:
        out = os.path.join(out_dir, f'cluster{c}')
        os.makedirs(out, exist_ok=True)
        tmp = sc.get.rank_genes_groups_df(adata, group=c)
        tmp = tmp[(tmp['logfoldchanges'].abs()>logfc) & (tmp['pvals_adj']<fdr) & (tmp['pct_nz_group'] >= pct)]
        tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
        tmp['cluster'] = c
        all_diff.append(tmp)


        sc.pl.rank_genes_groups_violin(adata, groups=c, gene_names=list(tmp['names'].head(5)), jitter=0, size=0, show=False)
        save_fig(plt, out, file_name='top_5_genes_violin')

    all_diff = pd.concat(all_diff, axis=0)
    all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
    all_diff.to_csv(os.path.join(out_dir, 'all_cluster_marker.xls'), sep='\t', index=None)
    all_diff['cluster'].value_counts().sort_index().to_csv(os.path.join(out_dir, 'cluster_diff_statistics.xls'), sep='\t')

    sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes=3,
        values_to_plot="logfoldchanges",
        vmin=-3,
        vmax=3,
        min_logfoldchange=1,
        cmap='bwr',
        colorbar_title='log fold change',
        show=False
    )
    save_fig(plt, out_dir, file_name='dotplot_top3_marker_gene', fig_size={'width':20, 'height':8})

    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=10,
        values_to_plot="logfoldchanges",
        min_logfoldchange=1,
        swap_axes=True,
        cmap='bwr',
        standard_scale='var',
        show_gene_labels=True,
        dendrogram=True,
        show=False
    )
    save_fig(plt, out_dir, file_name='heatmap_top10_marker_gene', fig_size={'width':26, 'height':20})

    adata.write(os.path.join(out_dir, 'adata_cluster.h5ad'))

    return adata

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Finde marker genes")
    parser.add_argument("--input_file", required=True, help="输入文件 h5ad")
    parser.add_argument("--out_dir", required=True, help="输出文件夹路径")
    parser.add_argument("--method_test", default="wilcoxon", help="选择计算Marker基因方法：wilcoxon 或 t-test")
    args = parser.parse_args()
    marker_gene_out_dir = os.path.join(args.out_dir, '05-Marker_gene_analysis')
    if not os.path.exists(marker_gene_out_dir):
        os.makedirs(marker_gene_out_dir)
    else:
        os.system(f'rm -rf {marker_gene_out_dir}')
        os.makedirs(marker_gene_out_dir)
    adata = sc.read_h5ad(args.input_file)
    adata_marker = find_marker_gene(adata, method='wilcoxon', out_dir=marker_gene_out_dir, logfc=0.5, fdr=0.05, pct=0.25)



