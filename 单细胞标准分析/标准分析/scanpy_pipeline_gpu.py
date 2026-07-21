import os
import shutil
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

from gpu_utils import Timer, check_gpu_available, safe_import_rapids_singlecell, to_cpu, to_gpu


sc.settings.verbosity = 3
sc.settings.cachedir = "/mnt/NexData_file/share/scRNA_platform/scVision4/tmp/.scanpy_cached"
sc.settings.set_figure_params(figsize=(5, 5), facecolor="white")


def _cpu_func(name):
    from scanpy_pipeline import __dict__ as cpu_namespace

    return cpu_namespace[name]


def make_report_dir(out_dir):
    return _cpu_func("make_report_dir")(out_dir)


def get_10x_mtx(data_path, sample_dict, out_dir):
    return _cpu_func("get_10x_mtx")(data_path, sample_dict, out_dir)


def data_qc(_adata, out_dir):
    return _cpu_func("data_qc")(_adata, out_dir)


def cell_filtering(
    _adata,
    min_genes,
    min_cells,
    pct_counts_mt,
    pct_counts_ribo,
    upper_lim,
    total_counts,
    expected_doublet_rate,
    threshold,
    out_dir,
):
    return _cpu_func("cell_filtering")(
        _adata,
        min_genes,
        min_cells,
        pct_counts_mt,
        pct_counts_ribo,
        upper_lim,
        total_counts,
        expected_doublet_rate,
        threshold,
        out_dir,
    )


def celltype_annotation(_adata, database, organism, out_dir, cell_type=None):
    return _cpu_func("celltype_annotation")(_adata, database, organism, out_dir, cell_type=cell_type)


def save_fig(plt_obj, out_dir, file_name, fig_size={"width": 5, "height": 5}):
    os.makedirs(out_dir, exist_ok=True)
    plt_obj.tight_layout()
    try:
        plt_obj.gcf().set_size_inches(fig_size["width"], fig_size["height"])
    except Exception:
        pass
    plt_obj.savefig(os.path.join(out_dir, f"{file_name}.png"), bbox_inches="tight", pad_inches=0.1)
    plt_obj.savefig(os.path.join(out_dir, f"{file_name}.pdf"), bbox_inches="tight", pad_inches=0.1)
    plt_obj.clf()
    plt_obj.close("all")


def _cpu_dim_reduction_cluster(_adata, method, n_pcs, n_neighbors, resolution, batch_effect, out_dir):
    return _cpu_func("dim_reduction_cluster")(
        _adata=_adata,
        method=method,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        resolution=resolution,
        batch_effect=batch_effect,
        out_dir=out_dir,
    )


def _cpu_find_marker_gene(_adata, method, out_dir, logfc, fdr, pct):
    return _cpu_func("find_marker_gene")(
        _adata=_adata,
        method=method,
        out_dir=out_dir,
        logfc=logfc,
        fdr=fdr,
        pct=pct,
    )


def _run_hvg(rsc, adata, n_top_genes):
    try:
        rsc.pp.highly_variable_genes(
            adata,
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
            n_top_genes=n_top_genes,
        )
    except TypeError:
        rsc.pp.highly_variable_genes(
            adata,
            min_mean=0.0125,
            max_mean=3,
            min_disp=0.5,
        )


def _run_neighbors(rsc, adata, n_neighbors, n_pcs, neighbors_algorithm):
    try:
        rsc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            algorithm=neighbors_algorithm,
        )
    except TypeError:
        rsc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)


def _run_harmony(rsc, adata, key):
    if hasattr(rsc.pp, "harmony_integrate"):
        try:
            rsc.pp.harmony_integrate(adata, key=key, basis="X_pca", adjusted_basis="X_pca")
            return adata
        except TypeError:
            rsc.pp.harmony_integrate(adata, key, adjusted_basis="X_pca")
            return adata
    to_cpu(adata)
    sce.pp.harmony_integrate(adata, key, adjusted_basis="X_pca")
    to_gpu(adata)
    return adata


def _run_bbknn(rsc, adata, batch_key, n_pcs):
    if hasattr(rsc.pp, "bbknn"):
        rsc.pp.bbknn(adata, batch_key=batch_key, n_pcs=n_pcs)
        return adata
    to_cpu(adata)
    sce.pp.bbknn(adata, batch_key=batch_key, n_pcs=n_pcs)
    to_gpu(adata)
    return adata


def dim_reduction_cluster_gpu(
    _adata,
    method,
    n_pcs,
    n_neighbors,
    resolution,
    batch_effect,
    out_dir,
    target_sum=None,
    n_top_genes=2000,
    random_state=1234,
    neighbors_algorithm="brute",
    fallback_cpu=True,
):
    os.makedirs(out_dir, exist_ok=True)
    rsc = safe_import_rapids_singlecell()
    if rsc is None or not check_gpu_available():
        if not fallback_cpu:
            raise RuntimeError("GPU unavailable and fallback_cpu=False")
        return _cpu_dim_reduction_cluster(_adata, method, n_pcs, n_neighbors, resolution, batch_effect, out_dir)

    if method not in {"leiden", "louvain"}:
        raise ValueError("method must be 'leiden' or 'louvain'")

    with Timer("save counts layer"):
        _adata.layers["counts"] = _adata.X.copy()

    with Timer("transfer AnnData to GPU"):
        to_gpu(_adata)

    with Timer("GPU normalize_total + log1p"):
        if target_sum is None:
            rsc.pp.normalize_total(_adata)
        else:
            rsc.pp.normalize_total(_adata, target_sum=target_sum)
        rsc.pp.log1p(_adata)

    with Timer("save raw log-normalized full gene matrix"):
        to_cpu(_adata)
        _adata.raw = _adata.copy()
        to_gpu(_adata)

    with Timer("GPU highly_variable_genes"):
        _run_hvg(rsc, _adata, n_top_genes)

    with Timer("plot highly variable genes"):
        to_cpu(_adata)
        sc.pl.highly_variable_genes(_adata, show=False)
        save_fig(plt, out_dir, file_name="highly_variable_genes", fig_size={"width": 16, "height": 8})
        to_gpu(_adata)

    with Timer("subset highly variable genes"):
        _adata = _adata[:, _adata.var["highly_variable"]].copy()

    with Timer("GPU regress_out + scale"):
        rsc.pp.regress_out(_adata, keys=["total_counts", "pct_counts_mt"])
        rsc.pp.scale(_adata, max_value=10)

    with Timer("GPU PCA"):
        try:
            rsc.pp.pca(_adata, n_comps=max(50, n_pcs), random_state=random_state)
        except TypeError:
            rsc.pp.pca(_adata, n_comps=max(50, n_pcs))

    with Timer("plot PCA"):
        to_cpu(_adata)
        sc.pl.pca(_adata, color="sample", show=False)
        save_fig(plt, out_dir, file_name="pca_sample")
        sc.pl.pca(_adata, color="group", show=False)
        save_fig(plt, out_dir, file_name="pca_group")
        sc.pl.pca_variance_ratio(_adata, log=True, show=False)
        save_fig(plt, out_dir, file_name="pca_variance_ratio")
        to_gpu(_adata)

    if batch_effect == "**harmony**":
        with Timer("GPU Harmony integrate by sample"):
            _run_harmony(rsc, _adata, "sample")

    with Timer("GPU neighbors"):
        _run_neighbors(rsc, _adata, n_neighbors, n_pcs, neighbors_algorithm)

    if batch_effect == "**bbknn**":
        with Timer("GPU BBKNN by sample"):
            _run_bbknn(rsc, _adata, "sample", n_pcs)

    with Timer("GPU UMAP"):
        try:
            rsc.tl.umap(_adata, random_state=random_state)
        except TypeError:
            rsc.tl.umap(_adata)

    with Timer(f"GPU clustering: {method}"):
        if method == "leiden":
            rsc.tl.leiden(_adata, resolution=resolution, key_added="cluster")
        else:
            rsc.tl.louvain(_adata, resolution=resolution, key_added="cluster")

    with Timer("transfer AnnData to CPU for plotting and tables"):
        to_cpu(_adata)

    with Timer("plot UMAP"):
        sc.pl.umap(
            _adata,
            color="cluster",
            title="Cluster",
            legend_loc="on data",
            legend_fontsize=10,
            legend_fontoutline=2,
            show=False,
        )
        save_fig(plt, out_dir, file_name="umap_dimplot")
        sc.pl.umap(_adata, color="sample", title="", show=False)
        save_fig(plt, out_dir, file_name="umap_sample_dimplot")
        sc.pl.umap(_adata, color="group", title="", show=False)
        save_fig(plt, out_dir, file_name="umap_group_dimplot")

    mt_column = "total_counts_mt" if "total_counts_mt" in _adata.obs.columns else "pct_counts_mt"
    with Timer("plot QC UMAP"):
        with rc_context({"figure.figsize": (4, 4)}):
            sc.pl.umap(
                _adata,
                color=["cluster", "total_counts", "n_genes_by_counts", mt_column],
                title=["Cluster", "total_counts", "n_genes_by_counts", mt_column],
                legend_loc="on data",
                frameon=True,
                legend_fontsize=9,
                legend_fontoutline=3,
                wspace=0.2,
                ncols=2,
                cmap="coolwarm",
                use_raw=True,
                show=False,
            )
            plt.savefig(
                os.path.join(out_dir, "cluster_total_counts_n_genes_by_counts.png"),
                facecolor="white",
                bbox_inches="tight",
                dpi=300,
            )
            plt.close("all")

    with Timer("write metadata and cluster statistics"):
        _adata.obs.to_csv(os.path.join(out_dir, "all_cells_features.xls"), sep="\t")
        result = (
            _adata.obs.groupby(["cluster", "sample"])
            .size()
            .reset_index(name="number")
            .pivot(index="cluster", columns="sample", values="number")
            .reset_index()
        )
        result.columns.name = None
        result.to_csv(os.path.join(out_dir, "cluster_cells_statistics.xls"), sep="\t", index=None)

    with Timer("plot cluster correlation"):
        with rc_context({"figure.figsize": (4, 4)}):
            try:
                sc.pl.correlation_matrix(_adata, "cluster", show=False)
                plt.savefig(
                    os.path.join(out_dir, "cluster_correlation.png"),
                    facecolor="white",
                    bbox_inches="tight",
                    dpi=300,
                )
            finally:
                plt.close("all")

    groups = _adata.obs["group"].drop_duplicates().tolist()
    if len(groups) > 1:
        with Timer("plot UMAP split by group"):
            fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups) * 5, 5))
            axs = np.atleast_1d(axs)
            for i, group in enumerate(groups):
                sc.pl.umap(
                    _adata[_adata.obs["group"] == group],
                    color="cluster",
                    title=group,
                    ax=axs[i],
                    wspace=0.1,
                    legend_loc="on data",
                    legend_fontsize=12,
                    legend_fontweight="normal",
                    legend_fontoutline=2,
                    show=False,
                )
            save_fig(plt, out_dir, file_name="umap_split_group_cluster", fig_size={"width": len(groups) * 5, "height": 5})

    samples = _adata.obs["sample"].drop_duplicates().tolist()
    if len(samples) > 1:
        with Timer("plot UMAP split by sample"):
            fig, axs = plt.subplots(nrows=1, ncols=len(samples), figsize=(len(samples) * 5, 5))
            axs = np.atleast_1d(axs)
            for i, sample in enumerate(samples):
                sc.pl.umap(
                    _adata[_adata.obs["sample"] == sample],
                    color="cluster",
                    title=sample,
                    ax=axs[i],
                    wspace=0.1,
                    legend_loc="on data",
                    legend_fontsize=12,
                    legend_fontweight="normal",
                    legend_fontoutline=2,
                    show=False,
                )
            save_fig(plt, out_dir, file_name="umap_split_sample_cluster", fig_size={"width": len(samples) * 5, "height": 5})

    return _adata


def _rank_genes_groups_gpu(rsc, adata, method, use_raw, chunk_size):
    kwargs = {
        "groupby": "cluster",
        "method": method,
        "pts": True,
        "use_raw": use_raw,
        "corr_method": "benjamini-hochberg",
    }
    try:
        rsc.tl.rank_genes_groups(adata, chunk_size=chunk_size, **kwargs)
    except TypeError:
        rsc.tl.rank_genes_groups(adata, **kwargs)


def _filter_marker_table(tmp, logfc, fdr, pct, positive_only):
    if "pct_nz_group" in tmp.columns:
        pct_col = "pct_nz_group"
    elif "pts" in tmp.columns:
        pct_col = "pts"
    else:
        pct_col = None

    if positive_only:
        mask = tmp["logfoldchanges"] > logfc
    else:
        mask = tmp["logfoldchanges"].abs() > logfc
    mask = mask & (tmp["pvals_adj"] < fdr)
    if pct_col is not None:
        mask = mask & (tmp[pct_col] >= pct)
    return tmp[mask].copy()


def find_marker_gene_gpu(
    _adata,
    method,
    out_dir,
    logfc=1,
    fdr=0.05,
    pct=0.25,
    use_raw=True,
    positive_only=True,
    chunk_size=128,
    fallback_cpu=True,
):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    rsc = safe_import_rapids_singlecell()
    if rsc is None or not check_gpu_available():
        if not fallback_cpu:
            raise RuntimeError("GPU unavailable and fallback_cpu=False")
        return _cpu_find_marker_gene(_adata, method, out_dir, logfc, fdr, pct)

    with Timer("transfer AnnData to GPU"):
        to_gpu(_adata)

    with Timer(f"GPU rank_genes_groups: {method}"):
        _rank_genes_groups_gpu(rsc, _adata, method, use_raw, chunk_size)

    with Timer("transfer AnnData to CPU for marker table and plotting"):
        to_cpu(_adata)

    with Timer("plot rank_genes_groups"):
        sc.pl.rank_genes_groups(_adata, n_genes=25, sharey=False, show=False)
        save_fig(plt, out_dir, file_name="rank_genes_groups", fig_size={"width": 20, "height": 12})

    result = _adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    all_diff = []

    for cluster in groups:
        cluster_dir = os.path.join(out_dir, f"cluster{cluster}")
        os.makedirs(cluster_dir, exist_ok=True)
        tmp = sc.get.rank_genes_groups_df(_adata, group=cluster)
        tmp = _filter_marker_table(tmp, logfc=logfc, fdr=fdr, pct=pct, positive_only=positive_only)
        tmp.to_csv(os.path.join(cluster_dir, "differentially_expressed_genes.xls"), sep="\t", index=None)
        tmp["cluster"] = cluster
        all_diff.append(tmp)

        if len(tmp) > 0:
            try:
                sc.pl.rank_genes_groups_violin(
                    _adata,
                    groups=cluster,
                    gene_names=list(tmp["names"].head(5)),
                    jitter=0,
                    size=0,
                    show=False,
                )
                save_fig(plt, cluster_dir, file_name="top_5_genes_violin")
            except Exception as exc:
                warnings.warn(f"top_5_genes_violin failed for cluster {cluster}: {exc}")

    if all_diff:
        all_diff_df = pd.concat(all_diff, axis=0)
    else:
        all_diff_df = pd.DataFrame()
    if "names" in all_diff_df.columns:
        all_diff_df.rename(columns={"names": "Symbol"}, inplace=True)
    all_diff_df.to_csv(os.path.join(out_dir, "all_cluster_marker.xls"), sep="\t", index=None)
    if "cluster" in all_diff_df.columns:
        all_diff_df["cluster"].value_counts().sort_index().to_csv(
            os.path.join(out_dir, "cluster_diff_statistics.xls"),
            sep="\t",
        )
    else:
        pd.Series(dtype=int).to_csv(os.path.join(out_dir, "cluster_diff_statistics.xls"), sep="\t")

    try:
        sc.pl.rank_genes_groups_dotplot(
            _adata,
            n_genes=3,
            values_to_plot="logfoldchanges",
            vmin=-3,
            vmax=3,
            min_logfoldchange=1,
            cmap="bwr",
            colorbar_title="log fold change",
            show=False,
        )
        save_fig(plt, out_dir, file_name="dotplot_top3_marker_gene", fig_size={"width": 20, "height": 8})
    except Exception as exc:
        warnings.warn(f"rank_genes_groups_dotplot failed: {exc}")

    try:
        sc.pl.rank_genes_groups_heatmap(
            _adata,
            n_genes=10,
            values_to_plot="logfoldchanges",
            min_logfoldchange=1,
            swap_axes=True,
            cmap="bwr",
            standard_scale="var",
            show_gene_labels=True,
            dendrogram=True,
            show=False,
        )
        save_fig(plt, out_dir, file_name="heatmap_top10_marker_gene", fig_size={"width": 26, "height": 20})
    except Exception as exc:
        warnings.warn(f"rank_genes_groups_heatmap failed: {exc}")

    return _adata


def visua_gene_express(_adata, genes, out_dir):
    return _cpu_func("visua_gene_express")(_adata, genes, out_dir)
