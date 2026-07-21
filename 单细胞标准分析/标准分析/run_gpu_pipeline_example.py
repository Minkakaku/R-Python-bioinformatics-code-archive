import argparse
import os

import scanpy as sc

from gpu_utils import check_gpu_available
from scanpy_pipeline_gpu import dim_reduction_cluster_gpu, find_marker_gene_gpu


def parse_args():
    parser = argparse.ArgumentParser(description="Run the GPU/CUDA Scanpy standard analysis steps.")
    parser.add_argument("--adata", help="Input h5ad file, usually 03-Cell_filter/adata_counts.h5ad")
    parser.add_argument("--out-root", help="Standard analysis output root directory")
    parser.add_argument("--method", default="leiden", choices=["leiden", "louvain"])
    parser.add_argument("--n-pcs", type=int, default=20)
    parser.add_argument("--n-neighbors", type=int, default=20)
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--batch-effect", default="**harmony**", choices=["none", "**harmony**", "**bbknn**"])
    parser.add_argument("--marker-method", default="wilcoxon")
    parser.add_argument("--logfc", type=float, default=1)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--pct", type=float, default=0.25)
    parser.add_argument("--neighbors-algorithm", default="brute", choices=["brute", "cagra"])
    parser.add_argument("--marker-method-fast", action="store_true", help="Use wilcoxon_binned for very large datasets.")
    parser.add_argument("--check-env", action="store_true", help="Only check GPU availability and import status.")
    return parser.parse_args()


def main():
    args = parse_args()
    gpu_ok = check_gpu_available()
    if args.check_env:
        print(f"GPU available: {gpu_ok}")
        return

    if not args.adata or not args.out_root:
        raise SystemExit("--adata and --out-root are required unless --check-env is used")

    batch_effect = None if args.batch_effect == "none" else args.batch_effect
    marker_method = "wilcoxon_binned" if args.marker_method_fast else args.marker_method

    adata = sc.read_h5ad(args.adata)
    adata = dim_reduction_cluster_gpu(
        _adata=adata,
        method=args.method,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        resolution=args.resolution,
        batch_effect=batch_effect,
        out_dir=os.path.join(args.out_root, "04-Dimension_reduction_cluster"),
        neighbors_algorithm=args.neighbors_algorithm,
    )
    adata = find_marker_gene_gpu(
        _adata=adata,
        method=marker_method,
        out_dir=os.path.join(args.out_root, "05-Marker_gene_analysis"),
        logfc=args.logfc,
        fdr=args.fdr,
        pct=args.pct,
        use_raw=True,
        positive_only=True,
    )
    adata.write(os.path.join(args.out_root, "result_gpu.h5ad"))


if __name__ == "__main__":
    main()
