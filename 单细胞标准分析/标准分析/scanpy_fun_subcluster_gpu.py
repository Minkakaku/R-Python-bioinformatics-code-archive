import os
import sys
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from matplotlib import pyplot as plt

from gpu_utils import Timer, check_gpu_available, safe_import_rapids_singlecell, to_cpu, to_gpu


_HERE = os.path.dirname(os.path.abspath(__file__))
_LEGACY_SUBCLUSTER_DIR = os.path.abspath(os.path.join(_HERE, "..", "未整理脚本"))
if _LEGACY_SUBCLUSTER_DIR not in sys.path:
    sys.path.append(_LEGACY_SUBCLUSTER_DIR)


def save_fig(plt_obj, out_dir=os.getcwd(), file_name=None, fig_size={"w": 5, "h": 5}):
    os.makedirs(out_dir, exist_ok=True)
    plt_obj.tight_layout()
    try:
        width = fig_size.get("w", fig_size.get("width", 5))
        height = fig_size.get("h", fig_size.get("height", 5))
        plt_obj.gcf().set_size_inches(width, height)
    except Exception:
        pass
    plt_obj.savefig(os.path.join(out_dir, f"{file_name}.png"), bbox_inches="tight", pad_inches=0.1)
    plt_obj.savefig(os.path.join(out_dir, f"{file_name}.pdf"), bbox_inches="tight", pad_inches=0.1)
    plt_obj.clf()
    plt_obj.close("all")


def _cpu_recluster(*args, **kwargs):
    from scanpy_fun_subcluster import recluster

    return recluster(*args, **kwargs)


def _run_neighbors(rsc, adata, n_neighbors, n_pcs, neighbors_algorithm):
    try:
        rsc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, algorithm=neighbors_algorithm)
    except TypeError:
        rsc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)


def _run_harmony(rsc, adata, key):
    if hasattr(rsc.pp, "harmony_integrate"):
        try:
            rsc.pp.harmony_integrate(adata, key=key, basis="X_pca", adjusted_basis="X_pca")
            return
        except TypeError:
            rsc.pp.harmony_integrate(adata, key, adjusted_basis="X_pca")
            return
    to_cpu(adata)
    sce.pp.harmony_integrate(adata, key=key, adjusted_basis="X_pca")
    to_gpu(adata)


def _run_bbknn(rsc, adata, batch_key, n_pcs):
    if hasattr(rsc.pp, "bbknn"):
        rsc.pp.bbknn(adata, batch_key=batch_key, n_pcs=n_pcs)
        return
    to_cpu(adata)
    sce.pp.bbknn(adata, batch_key=batch_key, n_pcs=n_pcs)
    to_gpu(adata)


class recluster_gpu:
    def __init__(
        self,
        adata=None,
        species=None,
        neighbors=30,
        resolution=0.5,
        n_pcs=10,
        method="leiden",
        batch_effect=0,
        output=os.getcwd(),
        n_top_genes=2000,
        random_state=1234,
        neighbors_algorithm="brute",
        fallback_cpu=True,
    ):
        self.adata = adata
        self.species = species
        self.neighbors = neighbors
        self.resolution = resolution
        self.output = output
        self.method = method
        self.batch = batch_effect
        self.n_pcs = n_pcs
        self.n_top_genes = n_top_genes
        self.random_state = random_state
        self.neighbors_algorithm = neighbors_algorithm
        self.fallback_cpu = fallback_cpu

    def sub_cluster(self):
        rsc = safe_import_rapids_singlecell()
        if rsc is None or not check_gpu_available():
            if not self.fallback_cpu:
                raise RuntimeError("GPU unavailable and fallback_cpu=False")
            cpu_runner = _cpu_recluster(
                adata=self.adata,
                species=self.species,
                neighbors=self.neighbors,
                resolution=self.resolution,
                n_pcs=self.n_pcs,
                method=self.method,
                batch_effect=self.batch,
                output=self.output,
            )
            return cpu_runner.sub_cluster()

        os.makedirs(self.output, exist_ok=True)

        with Timer("restore raw to adata"):
            self.adata = self.adata.raw.to_adata()
            adata_tmp = self.adata.copy()
            self.adata.uns.setdefault("log1p", {})
            self.adata.uns["log1p"]["base"] = None

        with Timer("transfer subcluster AnnData to GPU"):
            to_gpu(self.adata)

        with Timer("GPU highly_variable_genes for subcluster"):
            try:
                rsc.pp.highly_variable_genes(self.adata, n_top_genes=self.n_top_genes)
            except TypeError:
                rsc.pp.highly_variable_genes(self.adata)

        with Timer("save raw and subset highly variable genes"):
            to_cpu(self.adata)
            self.adata.raw = adata_tmp
            to_gpu(self.adata)
            self.adata = self.adata[:, self.adata.var.highly_variable].copy()

        with Timer("GPU PCA for subcluster"):
            try:
                rsc.pp.pca(self.adata, n_comps=max(50, self.n_pcs), random_state=self.random_state)
            except TypeError:
                rsc.pp.pca(self.adata, n_comps=max(50, self.n_pcs))

        with Timer("plot PCA for subcluster"):
            to_cpu(self.adata)
            sc.pl.pca(self.adata, color="sample", show=False)
            save_fig(plt, self.output, file_name="pca_sample")
            sc.pl.pca_variance_ratio(self.adata, n_pcs=30, show=False)
            save_fig(plt, self.output, file_name="pca_variance")
            to_gpu(self.adata)

        with Timer("GPU neighbors for subcluster"):
            _run_neighbors(rsc, self.adata, self.neighbors, self.n_pcs, self.neighbors_algorithm)

        if self.batch == 1:
            with Timer("GPU BBKNN subcluster by sample"):
                _run_bbknn(rsc, self.adata, "sample", self.n_pcs)
        elif self.batch == 2:
            with Timer("GPU Harmony subcluster by sample"):
                _run_harmony(rsc, self.adata, "sample")
            with Timer("GPU neighbors after Harmony"):
                _run_neighbors(rsc, self.adata, self.neighbors, self.n_pcs, self.neighbors_algorithm)
        elif self.batch == 3:
            with Timer("GPU Harmony subcluster by group"):
                _run_harmony(rsc, self.adata, "group")
            with Timer("GPU neighbors after Harmony"):
                _run_neighbors(rsc, self.adata, self.neighbors, self.n_pcs, self.neighbors_algorithm)

        with Timer("GPU UMAP for subcluster"):
            try:
                rsc.tl.umap(self.adata, random_state=self.random_state)
            except TypeError:
                rsc.tl.umap(self.adata)

        with Timer(f"GPU subcluster clustering: {self.method}"):
            if self.method == "leiden":
                rsc.tl.leiden(self.adata, resolution=self.resolution)
            else:
                rsc.tl.louvain(self.adata, resolution=self.resolution)

        with Timer("GPU rank_genes_groups for subcluster"):
            rsc.tl.rank_genes_groups(
                self.adata,
                groupby=self.method,
                method="wilcoxon",
                pts=True,
                use_raw=True,
                corr_method="benjamini-hochberg",
            )

        with Timer("transfer subcluster AnnData to CPU for output"):
            to_cpu(self.adata)

        sc.pl.umap(self.adata, color=self.method, show=False)
        save_fig(plt, self.output, file_name="sub_cluster_umap")
        sc.pl.umap(self.adata, color=self.method, legend_loc="on data", show=False)
        save_fig(plt, self.output, file_name="sub_cluster_umap_withlabel")

        sc.pl.rank_genes_groups(self.adata, n_genes=10, sharey=False, show=False)
        save_fig(plt, self.output, file_name="rank_genes_groups", fig_size={"width": 20, "height": 12})

        find_group_marker_gpu(adata=self.adata, species=self.species, output=self.output, groupby_key=self.method)

        groups = self.adata.obs["group"].drop_duplicates().tolist()
        if len(groups) > 1:
            fig, axs = plt.subplots(nrows=1, ncols=len(groups), figsize=(len(groups) * 6, 6))
            axs = np.atleast_1d(axs)
            keys = self.adata.obs[self.method].cat.categories.to_list()
            values = list(self.adata.uns[self.method + "_colors"])
            color_dict = dict(zip(keys, values))
            for i, group in enumerate(groups):
                tmp = self.adata[self.adata.obs["group"] == group]
                sc.pl.umap(
                    tmp,
                    size=20,
                    color=self.method,
                    legend_loc="on data",
                    palette=color_dict,
                    title=group,
                    ax=axs[i],
                    show=False,
                )
            save_fig(plt, out_dir=self.output, file_name="umap_split_groups", fig_size={"w": len(groups) * 6, "h": 6})

        self.adata.write_h5ad(os.path.join(self.output, "result.h5ad"))
        return self.adata


def find_group_marker_gpu(adata, species="human", output=os.getcwd(), groupby_key="leiden", logfc=1, fdr=0.05, pct=0.25):
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    all_diff = []
    for cluster in groups:
        out = os.path.join(output, f"cluster{cluster}")
        os.makedirs(out, exist_ok=True)
        tmp = sc.get.rank_genes_groups_df(adata, group=cluster)
        if "pct_nz_group" in tmp.columns:
            tmp = tmp[(tmp["logfoldchanges"] > logfc) & (tmp["pvals_adj"] < fdr) & (tmp["pct_nz_group"] >= pct)].copy()
        else:
            tmp = tmp[(tmp["logfoldchanges"] > logfc) & (tmp["pvals_adj"] < fdr)].copy()
        tmp.to_csv(os.path.join(out, "differentially_expressed_genes.xls"), sep="\t", index=None)
        tmp["cluster"] = cluster
        all_diff.append(tmp)

        try:
            import gseapy as gp

            if len(tmp) > 0 and species == "human":
                gp.enrichr(
                    gene_list=list(tmp["names"]),
                    gene_sets=["GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "KEGG_2021_Human"],
                    organism="Human",
                    cutoff=0.05,
                    no_plot=False,
                    outdir=out,
                    verbose=False,
                )
            elif len(tmp) > 0 and species == "mouse":
                gp.enrichr(
                    gene_list=list(tmp["names"]),
                    gene_sets=["GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023", "KEGG_2019_Mouse"],
                    organism="Mouse",
                    cutoff=0.05,
                    no_plot=False,
                    outdir=out,
                    verbose=False,
                )
        except Exception as exc:
            warnings.warn(f"gseapy enrichment skipped for cluster {cluster}: {exc}")

    all_diff_df = pd.concat(all_diff, axis=0) if all_diff else pd.DataFrame()
    if "names" in all_diff_df.columns:
        all_diff_df.rename(columns={"names": "Symbol"}, inplace=True)
    all_diff_df.to_csv(os.path.join(output, "all_cluster_marker.xls"), sep="\t", index=None)
    if "cluster" in all_diff_df.columns:
        all_diff_df["cluster"].value_counts().sort_index().to_csv(os.path.join(output, "cluster_diff_statistics.xls"), sep="\t")
    else:
        pd.Series(dtype=int).to_csv(os.path.join(output, "cluster_diff_statistics.xls"), sep="\t")
