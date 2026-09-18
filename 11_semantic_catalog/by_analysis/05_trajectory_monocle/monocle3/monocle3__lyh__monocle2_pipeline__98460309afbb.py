from __future__ import annotations

import argparse
import os
import types
import warnings
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the OmicVerse Monocle2 trajectory workflow."
    )
    parser.add_argument("-i", "--input", required=True, help="Input .h5ad file.")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory.")
    parser.add_argument("--threads", type=int, default=24, help="CPU threads.")
    parser.add_argument("--mean-expr-thresh", type=float, default=0.1)
    parser.add_argument("--max-genes", type=int, default=1000)
    parser.add_argument("--max-components", type=int, default=2)
    parser.add_argument("--method", default="fast", choices=["fast", "exact"])
    parser.add_argument("--ncenter", type=int, default=None)
    parser.add_argument("--param-gamma", type=float, default=10)
    parser.add_argument("--sigma", type=float, default=None)
    parser.add_argument("--lambda-param", type=float, default=None)
    parser.add_argument("--root-state", default=None)
    parser.add_argument("--reverse", action="store_true")
    parser.add_argument("--state-key", default="State")
    parser.add_argument("--subtype-key", default="subtype")
    parser.add_argument("--group-key", default="group")
    parser.add_argument(
        "--trajectory-colors",
        default="State,subtype,group",
        help="Comma-separated obs columns for trajectory plots.",
    )
    parser.add_argument("--cell-size", type=float, default=1)
    parser.add_argument("--split-size", type=float, default=40)
    parser.add_argument("--de-qval", type=float, default=0.01)
    parser.add_argument("--heatmap-top-features", type=int, default=10000)
    parser.add_argument("--heatmap-splits", type=int, default=3)
    parser.add_argument("--skip-de", action="store_true")
    parser.add_argument("--skip-h5ad", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def setup_environment(threads: int) -> None:
    os.environ.setdefault("OPENBLAS_NUM_THREADS", str(threads))
    os.environ.setdefault("OMP_NUM_THREADS", str(threads))
    os.environ.setdefault("MKL_NUM_THREADS", str(threads))
    os.environ.setdefault("NUMEXPR_NUM_THREADS", str(threads))
    os.environ.setdefault("MPLBACKEND", "Agg")


def extract_dynamic_heatmap_gene_cluster(
    adata,
    pseudotime,
    var_names=None,
    *,
    dynamic_heatmap_func=None,
    lineage_key=None,
    lineages=None,
    max_lineages=2,
    top_features=None,
    cell_bins=100,
    use_cell_columns=True,
    use_fitted=True,
    aggregate="mean",
    smooth_window=15,
    score_smooth_window=21,
    fitted_window=31,
    layer=None,
    use_raw=False,
    gene_symbols=None,
    standard_scale="var",
    order_by="peak",
    reverse_ht=None,
    n_split=4,
    split_method="kmeans-peaktime",
    cluster_features_by=None,
    row_cluster=False,
    cell_annotation=None,
    verbose=True,
):
    import inspect
    import pandas as pd
    import omicverse as ov

    if n_split is None or int(n_split) <= 0:
        raise ValueError("`n_split` must be a positive integer.")

    if dynamic_heatmap_func is None:
        dynamic_heatmap_func = ov.pl.dynamic_heatmap

    dynamic_heatmap_func = inspect.unwrap(dynamic_heatmap_func)
    G = dynamic_heatmap_func.__globals__

    _dynamic_hvg_candidates = G["_dynamic_hvg_candidates"]
    _build_gene_groups = G["_build_gene_groups"]
    _sanitize_var_names = G["_sanitize_var_names"]
    _prepare_dynamic_matrix = G["_prepare_dynamic_matrix"]
    _fit_dynamic_matrix = G["_fit_dynamic_matrix"]
    _smooth_dynamic_matrix = G["_smooth_dynamic_matrix"]
    _compute_dynamic_feature_metadata = G["_compute_dynamic_feature_metadata"]
    _score_dynamic_features = G["_score_dynamic_features"]
    _filter_dynamic_candidates = G["_filter_dynamic_candidates"]
    _split_dynamic_features = G["_split_dynamic_features"]
    _cluster_dynamic_rows = G["_cluster_dynamic_rows"]

    if layer is not None:
        use_raw = False
    if use_raw and adata.raw is None:
        raise ValueError("use_raw=True was requested, but adata.raw is not available.")

    if var_names is None:
        selected_var_names = _dynamic_hvg_candidates(adata, use_raw)
        if not selected_var_names:
            selected_var_names = (
                list(adata.raw.var_names) if use_raw else list(adata.var_names)
            )
    else:
        selected_var_names, _ = _build_gene_groups(var_names)
        if selected_var_names is None:
            raise ValueError("Please provide `var_names`.")
        selected_var_names = _sanitize_var_names(adata, selected_var_names, use_raw)

    annotation_keys = None
    if cell_annotation is not None:
        annotation_keys = (
            [cell_annotation] if isinstance(cell_annotation, str) else list(cell_annotation)
        )

    matrix, metadata = _prepare_dynamic_matrix(
        adata,
        selected_var_names,
        pseudotime,
        lineage_key=lineage_key,
        lineages=lineages,
        use_raw=use_raw,
        layer=layer,
        gene_symbols=gene_symbols,
        n_bins=cell_bins,
        aggregate=aggregate,
        standard_scale=standard_scale,
        reverse_ht=reverse_ht,
        use_cell_columns=use_cell_columns,
        annotation_keys=annotation_keys,
    )

    if lineage_key is not None and "lineage" in metadata.columns and max_lineages is not None:
        lineage_order = list(dict.fromkeys(metadata["lineage"].astype(str)))
        keep_lineages = lineage_order[: max(int(max_lineages), 1)]
        keep_columns = metadata.index[metadata["lineage"].astype(str).isin(keep_lineages)]
        metadata = metadata.loc[keep_columns]
        matrix = matrix.loc[:, keep_columns]

    if use_fitted:
        matrix = _fit_dynamic_matrix(matrix, metadata, window=fitted_window)

    matrix = _smooth_dynamic_matrix(matrix, metadata, window=smooth_window)

    order_by = str(order_by).lower()
    if order_by not in {"peak", "peaktime", "valley", "valleytime"}:
        raise ValueError("order_by must be one of {'peak', 'peaktime', 'valley', 'valleytime'}.")
    order_mode = "peak" if order_by in {"peak", "peaktime"} else "valley"

    feature_metadata = _compute_dynamic_feature_metadata(
        matrix,
        metadata,
        order_by=order_mode,
        cluster_features_by=cluster_features_by,
    )

    if top_features is not None or var_names is None:
        dynamic_score = _score_dynamic_features(
            matrix,
            metadata,
            smooth_window=score_smooth_window,
        )
        if var_names is None:
            dynamic_keep = _filter_dynamic_candidates(dynamic_score)
            if dynamic_keep:
                matrix = matrix.loc[dynamic_keep]
                feature_metadata = feature_metadata.loc[dynamic_keep]
                dynamic_score = dynamic_score.loc[dynamic_keep]

        if top_features is not None and top_features > 0 and matrix.shape[0] > top_features:
            keep_features = dynamic_score.index[: int(top_features)]
            matrix = matrix.loc[keep_features]
            feature_metadata = feature_metadata.loc[keep_features]

    ordered_rows = feature_metadata.sort_values(
        ["cluster_time", "peak_lineage"],
        kind="stable",
    ).index
    matrix = matrix.loc[ordered_rows]
    feature_metadata = feature_metadata.loc[ordered_rows]

    row_split = _split_dynamic_features(
        matrix,
        feature_metadata,
        n_split=n_split,
        split_method=split_method,
    )
    if row_split is None:
        raise ValueError("No expression-pattern clusters were generated.")

    matrix = matrix.loc[row_split.index]
    feature_metadata = feature_metadata.loc[row_split.index]

    if row_cluster:
        clustered_rows = _cluster_dynamic_rows(matrix, feature_metadata)
        row_split = row_split.loc[clustered_rows]

    gene_cluster = pd.DataFrame(
        {
            "gene": row_split.index.astype(str),
            "cluster": row_split.astype(str).to_numpy(),
        }
    )

    if verbose:
        print(
            "Extracted gene-cluster table: "
            f"{gene_cluster.shape[0]} genes x {gene_cluster['cluster'].nunique()} clusters"
        )
    return gene_cluster


def save_current_figure(path: Path, width: float, height: float, dpi: int = 300) -> None:
    import matplotlib.pyplot as plt

    plt.gcf().set_size_inches(width, height)
    plt.savefig(path.with_suffix(".png"), format="png", bbox_inches="tight", dpi=dpi)
    plt.savefig(path.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
    plt.close()


def has_obs_key(adata, key: str) -> bool:
    return bool(key) and key in adata.obs.columns


def clean_uns_for_h5ad(adata) -> None:
    def clean_mapping(mapping):
        for key, value in list(mapping.items()):
            if isinstance(value, dict):
                clean_mapping(value)
            elif isinstance(value, types.FunctionType):
                mapping.pop(key, None)

    clean_mapping(adata.uns)
    if "monocle" in adata.uns:
        for key in ["mst", "pr_graph_cell_proj_tree", "disp_func"]:
            adata.uns["monocle"].pop(key, None)


def main() -> None:
    args = parse_args()
    setup_environment(args.threads)
    warnings.filterwarnings("ignore", category=FutureWarning)

    import math

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import omicverse as ov
    from omicverse.single import Monocle

    np.random.seed(42)
    ov.plot_set(font_path="Arial")

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[monocle2] Reading: {input_path}")
    adata = ov.read(str(input_path))
    adata.obs = adata.obs.astype(str)

    mono = Monocle(adata)
    mono.preprocess()
    mono.select_ordering_genes(
        mean_expr_thresh=args.mean_expr_thresh,
        max_genes=args.max_genes,
    )
    print(f"[monocle2] Ordering genes: {mono.adata.var['use_for_ordering'].sum()}")

    mono.plot_ordering_genes(figsize=(5, 4))
    save_current_figure(outdir / "1_Monocle_ordering_genes", 5, 4)

    reduce_kwargs = {
        "max_components": args.max_components,
        "verbose": args.verbose,
        "method": args.method,
        "param_gamma": args.param_gamma,
    }
    if args.ncenter is not None:
        reduce_kwargs["ncenter"] = args.ncenter
    if args.sigma is not None:
        reduce_kwargs["sigma"] = args.sigma
    if args.lambda_param is not None:
        reduce_kwargs["lambda_param"] = args.lambda_param

    mono.reduce_dimension(**reduce_kwargs)
    root_state = args.root_state
    if root_state is not None:
        try:
            root_state = int(root_state)
        except ValueError:
            pass
    mono.order_cells(root_state=root_state, reverse=args.reverse or None)
    print(mono)

    trajectory_colors = [
        key.strip() for key in args.trajectory_colors.split(",") if key.strip()
    ]
    for color in trajectory_colors:
        if not has_obs_key(mono.adata, color):
            print(f"[monocle2] Skip trajectory color '{color}': not found in adata.obs")
            continue
        ov.pl.trajectory(
            mono.adata,
            method="monocle",
            basis="X_DDRTree",
            color=color,
            show=False,
        )
        if color == args.state_key:
            name = "1_Monocle_State"
        elif color == args.subtype_key:
            name = "1_Monocle_celltype"
        elif color == args.group_key:
            name = "1_Monocle_group"
        else:
            name = f"1_Monocle_{color}"
        save_current_figure(outdir / name, 6, 4)

    if has_obs_key(mono.adata, args.subtype_key):
        ov.pl.trajectory_tree(
            mono.adata,
            method="monocle",
            color=args.subtype_key,
            cell_size=args.cell_size,
        )
        save_current_figure(outdir / "1_Monocle_celltype_tree", 5, 5)

    if has_obs_key(mono.adata, args.group_key) and has_obs_key(mono.adata, args.subtype_key):
        adata_plot = mono.adata
        group_key = args.group_key
        groups = adata_plot.obs[group_key].dropna().unique().tolist()
        n = len(groups)
        if n > 0:
            ncols = 3
            nrows = math.ceil(n / ncols)
            fig, axes = plt.subplots(
                nrows=nrows,
                ncols=ncols,
                figsize=(6 * ncols, 4 * nrows),
                squeeze=False,
            )
            for i, group in enumerate(groups):
                ax = axes.flat[i]
                adata_sub = adata_plot[adata_plot.obs[group_key] == group].copy()
                ov.pl.trajectory(
                    adata_sub,
                    method="monocle",
                    basis="X_DDRTree",
                    color=args.subtype_key,
                    show=False,
                    ax=ax,
                    frameon=False,
                    title=str(group),
                    show_tree=True,
                    show_branch_points=True,
                    size=args.split_size,
                )
                ax.set_title(str(group), fontsize=12)
                ax.set_xlabel("")
                ax.set_ylabel("")
            for j in range(i + 1, nrows * ncols):
                fig.delaxes(axes.flat[j])
            plt.tight_layout()
            plt.savefig(
                outdir / "1_Monocle_group_split.png",
                format="png",
                bbox_inches="tight",
                dpi=300,
            )
            plt.savefig(
                outdir / "1_Monocle_group_split.pdf",
                format="pdf",
                bbox_inches="tight",
            )
            plt.close(fig)

    if has_obs_key(mono.adata, args.subtype_key):
        state_int = mono.adata.obs[args.state_key].astype(int)
        state_subtype = pd.crosstab(state_int, mono.adata.obs[args.subtype_key])
        state_subtype_norm = state_subtype.div(state_subtype.sum(axis=1), axis=0)
        fig, ax = plt.subplots(figsize=(4, 8))
        im = ax.imshow(state_subtype_norm.values, aspect="auto", cmap="RdBu")
        ax.set_xticks(range(state_subtype_norm.shape[1]))
        ax.set_xticklabels(state_subtype_norm.columns, rotation=30, ha="right")
        ax.set_yticks(range(state_subtype_norm.shape[0]))
        ax.set_yticklabels(state_subtype_norm.index)
        ax.set_xlabel(args.subtype_key)
        ax.set_ylabel(args.state_key)
        ax.set_title(f"{args.state_key} / {args.subtype_key} correspondence")
        plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
        plt.savefig(
            outdir / "1_Monocle_correspondence_heatmap.png",
            format="png",
            bbox_inches="tight",
            dpi=300,
        )
        plt.savefig(
            outdir / "1_Monocle_correspondence_heatmap.pdf",
            format="pdf",
            bbox_inches="tight",
        )
        plt.close(fig)

        try:
            ov.pl.branch_streamplot(
                mono.adata,
                group_key=args.subtype_key,
                pseudotime_key="Pseudotime",
                show=False,
            )
            save_current_figure(outdir / "1_Monocle_subtype_Pseudotime", 6, 4)
        except Exception as exc:
            print(f"[monocle2] Skip branch_streamplot: {exc}")

    if not args.skip_de:
        ordering_genes = mono.adata.var_names[mono.adata.var["use_for_ordering"]].tolist()
        mono_ord = Monocle(mono.adata[:, ordering_genes].copy())
        de = mono_ord.differential_gene_test(cores=args.threads)
        de.to_csv(outdir / "2_Monocle_differential_gene_test.csv")

        sig = de[(de["qval"] < args.de_qval) & (de["status"] == "OK")]
        all_gene = sig.sort_values("pval").index.tolist()
        print(f"[monocle2] Significant dynamic genes: {len(all_gene)}")

        if all_gene:
            mono.adata.obs[args.state_key] = mono.adata.obs[args.state_key].astype(str)
            annotations = [
                key for key in [args.state_key, args.subtype_key] if has_obs_key(mono.adata, key)
            ]
            ov.pl.dynamic_heatmap(
                mono.adata,
                pseudotime="Pseudotime",
                var_names=all_gene,
                cell_annotation=annotations or None,
                use_cell_columns=False,
                use_fitted=True,
                cell_bins=100,
                figsize=(7, 7),
                show_row_names=False,
                standard_scale="var",
                cmap="RdBu_r",
                order_by="peak",
                split_method="kmeans-peaktime",
                show=False,
                n_split=args.heatmap_splits,
            )
            save_current_figure(outdir / "2_Monocle_Pesudotime_heatmap", 7, 7)

            gene_cluster = extract_dynamic_heatmap_gene_cluster(
                mono.adata,
                pseudotime="Pseudotime",
                var_names=all_gene,
                top_features=args.heatmap_top_features,
                n_split=args.heatmap_splits,
                split_method="kmeans-peaktime",
            )
            gene_cluster.to_csv(outdir / "Pesudotime_heatmap_gene.csv", index=False)

    if not args.skip_h5ad:
        adata_save = mono.adata.copy()
        clean_uns_for_h5ad(adata_save)
        adata_save.write_h5ad(outdir / "add_monocle.h5ad")

    print(f"[monocle2] Done. Outputs: {outdir}")


if __name__ == "__main__":
    main()
