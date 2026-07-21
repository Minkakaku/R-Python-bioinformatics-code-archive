import argparse
from Velocity_module import run_velocity_pipeline

# 定义主函数
def main(
    input_h5ad,
    loom_dir,
    output_dir,
    rename_file=None,
    sample_col="sample",
    groupby="cell_type",
    add_sample_prefix=True,
    # mode="stochastic",
    basis="umap",
    # n_top_genes=2000,
    # n_pcs=30,
    # n_neighbors=30,
    # min_shared_counts=20,
    prefix="scvelo",
    rank_velocity_genes=True,
    top_gene_n=5,
):
    """Main function: read h5ad + loom files, merge, run scVelo, and save plots/results."""
    run_velocity_pipeline(
        input_h5ad=input_h5ad,
        loom_dir=loom_dir,
        output_dir=output_dir,
        rename_file=rename_file,
        sample_col=sample_col,
        groupby=groupby,
        add_sample_prefix=add_sample_prefix,
        # mode=mode,
        basis=basis,
        # n_top_genes=n_top_genes,
        # n_pcs=n_pcs,
        # n_neighbors=n_neighbors,
        # min_shared_counts=min_shared_counts,
        prefix=prefix,
        rank_velocity_genes=rank_velocity_genes,
        top_gene_n=top_gene_n,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run RNA velocity analysis with scVelo.")

    parser.add_argument("--input_h5ad", type=str, required=True, help="Input h5ad file for selected cells.")  # 这个是你所需要做分析的h5ad文件
    parser.add_argument("--loom_dir", type=str, required=True, help="Directory containing loom files.") # loom文件所在的文件夹
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory for h5ad, tables and plots.") # 输出路径

    parser.add_argument("--rename_file", type=str, default=None, required=False,
                        help="Optional sample rename csv. First column: old loom sample; second column: new sample.") # csv格式的重命名文件，方便匹配loom文件的barcode和h5ad的，这个如果不需要重命名直接不传就行
    parser.add_argument("--sample_col", type=str, default="sample", required=False,
                        help="Sample column in input h5ad obs. Default: sample.") # h5ad文件样品名所在列名称,默认是sample列
    parser.add_argument("--groupby", type=str, default="cell_type", required=False, 
                        help="Group column for plotting and summaries. Default: cell_type.") # 绘图指定颜色的列的名称，包括计算velocity相关基因的时候所指定的列，默认是cell_type

    parser.add_argument("--add_sample_prefix", action="store_true", default=False, # 是否对h5ad进行barcode调整，默认是需要的，目前测试的文件需要
                        help="Add sample prefix to h5ad obs_names as sample:barcode.")

    # parser.add_argument("--mode", type=str, default="stochastic", choices=["stochastic", "deterministic", "dynamical"], # 默认"stochastic"
    #                     help="Velocity mode. Default: stochastic.")
    parser.add_argument("--basis", type=str, default="umap", help="Embedding basis. Default: umap.") # 绘图坐标，默认umap

    # parser.add_argument("--n_top_genes", type=int, default=2000, help="Top genes for scVelo normalization.") # velocity分析参数，默认2000
    # parser.add_argument("--n_pcs", type=int, default=30, help="Number of PCs for moments.") # velocity分析参数，默认30
    # parser.add_argument("--n_neighbors", type=int, default=30, help="Number of neighbors for moments.") # velocity分析参数，默认30
    # parser.add_argument("--min_shared_counts", type=int, default=20, help="min_shared_counts for scv.pp.filter_and_normalize.") # velocity分析参数，默认20

    parser.add_argument("--prefix", type=str, default="scvelo", help="Output file prefix.") # 绘图文件前缀，可以修改
    parser.add_argument("--top_gene_n", type=int, default=5, help="Top velocity genes per group to plot.")  # 最后velocity基因绘图数量，默认5
    parser.add_argument("--no_rank_velocity_genes", dest="rank_velocity_genes", action="store_false", # 是否进行velocity基因绘图
                        help="Skip rank_velocity_genes and gene velocity plots.")
    parser.set_defaults(rank_velocity_genes=True) # 默认绘制elocity基因

    args = parser.parse_args()

    main(
        input_h5ad=args.input_h5ad, 
        loom_dir=args.loom_dir,
        output_dir=args.output_dir,
        rename_file=args.rename_file,
        sample_col=args.sample_col,
        groupby=args.groupby,
        add_sample_prefix=args.add_sample_prefix,
        # mode=args.mode,
        basis=args.basis,
        # n_top_genes=args.n_top_genes,
        # n_pcs=args.n_pcs,
        # n_neighbors=args.n_neighbors,
        # min_shared_counts=args.min_shared_counts,
        prefix=args.prefix,
        rank_velocity_genes=args.rank_velocity_genes,
        top_gene_n=args.top_gene_n,
    )
