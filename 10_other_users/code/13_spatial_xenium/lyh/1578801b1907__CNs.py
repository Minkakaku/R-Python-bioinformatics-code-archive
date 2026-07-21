import pathlib
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import MiniBatchKMeans
from sklearn.neighbors import NearestNeighbors
from yellowbrick.cluster import KElbowVisualizer


def _normalize_mode(mode: str) -> str:
    normalized = str(mode).lower()
    if normalized == "dadius":
        normalized = "radius"
    if normalized not in {"knn", "radius"}:
        raise ValueError(f"Unknown mode: {mode}. Must be 'knn' or 'radius'.")
    return normalized


def _validate_inputs(df, region_col, cluster_col, X, Y, mode, k, radius, n_neighborhoods):
    required_columns = [X, Y, cluster_col, region_col]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise KeyError(f"Missing required columns: {missing_columns}")

    if k is not None and k <= 0:
        raise ValueError("`k` must be a positive integer")
    if mode == "radius" and radius is None:
        raise ValueError("`radius` must be set when `mode` is 'radius'")
    if radius is not None and radius <= 0:
        raise ValueError("`radius` must be positive")
    if n_neighborhoods <= 0:
        raise ValueError("`n_neighborhoods` must be a positive integer")


def _build_neighborhood_matrix(
    df: pd.DataFrame,
    region_col: str,
    cluster_col: str,
    X: str,
    Y: str,
    mode: str,
    k: int,
    radius: float,
):
    cells = df[[X, Y, cluster_col, region_col]].copy()
    cells[region_col] = cells[region_col].astype(str)
    dummies = pd.get_dummies(cells[cluster_col])
    sum_cols = dummies.columns.tolist()
    values = dummies.to_numpy()

    windows = []
    for _, idx in cells.groupby(region_col).groups.items():
        coords = df.loc[idx, [X, Y]].to_numpy()
        if len(coords) == 0:
            continue

        if mode == "knn":
            n_neighbors = min(int(k), len(coords))
            nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm="auto").fit(coords)
            neigh_idx = nbrs.kneighbors(return_distance=False)
        else:
            nbrs = NearestNeighbors(algorithm="auto").fit(coords)
            neigh_idx = nbrs.radius_neighbors(coords, radius=radius, return_distance=False)

        region_pos_idx = df.index.get_indexer_for(idx)
        region_values = values[region_pos_idx]

        window_vals = []
        for local_idx, neighbors in enumerate(neigh_idx):
            neighbors = np.asarray(neighbors, dtype=int)
            if neighbors.size == 0:
                neighbors = np.array([local_idx], dtype=int)
            window_vals.append(region_values[neighbors].sum(axis=0))

        tmp = pd.DataFrame(np.asarray(window_vals, dtype=np.float32), index=idx, columns=sum_cols)
        tmp[[X, Y, region_col, cluster_col]] = df.loc[idx, [X, Y, region_col, cluster_col]]
        windows.append(tmp)

    if not windows:
        raise ValueError("No valid neighborhoods were generated from the input data")

    window_df = pd.concat(windows).loc[df.index]
    X_vals = window_df[sum_cols].to_numpy()
    if len(X_vals) == 0:
        raise ValueError("No rows available for clustering after neighborhood construction")

    return window_df, X_vals


def determine_optimal_cluster_count(
    df: pd.DataFrame,
    region_col: str,
    cluster_col: str,
    X: str = "x",
    Y: str = "y",
    mode: str = "knn",
    k: int = 35,
    radius: float = None,
    max_clusters: int = 30,
    metric: str = "distortion",
):
    mode = _normalize_mode(mode)
    _validate_inputs(df, region_col, cluster_col, X, Y, mode, k, radius, max_clusters)

    _, X_vals = _build_neighborhood_matrix(
        df=df,
        region_col=region_col,
        cluster_col=cluster_col,
        X=X,
        Y=Y,
        mode=mode,
        k=k,
        radius=radius,
    )

    effective_max_clusters = min(int(max_clusters), len(X_vals))
    if effective_max_clusters < 2:
        raise ValueError("Elbow mode requires at least 2 neighborhoods to evaluate cluster counts")

    km = MiniBatchKMeans(n_clusters=2, random_state=0)
    visualizer = KElbowVisualizer(
        km,
        k=list(range(2, effective_max_clusters + 1)),
        metric=metric,
        timings=False,
    )
    visualizer.fit(X_vals)
    visualizer.show()

    best_k = visualizer.elbow_value_
    if best_k is None:
        best_k = effective_max_clusters
        warnings.warn(
            f"Failed to infer an elbow point automatically; falling back to {best_k} clusters.",
            RuntimeWarning,
        )

    return {
        "selected_clusters": int(best_k),
        "effective_max_clusters": int(effective_max_clusters),
        "metric": metric,
    }


def neighborhood_analysis_df(
    df: pd.DataFrame,
    region_col: str,
    cluster_col: str,
    X: str = "x",
    Y: str = "y",
    mode: str = "knn",
    k: int = 35,
    radius: float = None,
    n_neighborhoods: int = 30,
    elbow: bool = False,
    metric: str = "distortion",
    return_metadata: bool = False,
):
    """
    Compute Cellular Neighborhoods (CNs) from a pandas DataFrame.
    """
    mode = _normalize_mode(mode)
    _validate_inputs(df, region_col, cluster_col, X, Y, mode, k, radius, n_neighborhoods)

    _, X_vals = _build_neighborhood_matrix(
        df=df,
        region_col=region_col,
        cluster_col=cluster_col,
        X=X,
        Y=Y,
        mode=mode,
        k=k,
        radius=radius,
    )

    effective_max_clusters = min(int(n_neighborhoods), len(X_vals))
    if effective_max_clusters < 1:
        raise ValueError("`n_neighborhoods` must be at least 1")

    if elbow:
        if effective_max_clusters < 2:
            raise ValueError("Elbow mode requires at least 2 neighborhoods to evaluate cluster counts")
        km = MiniBatchKMeans(n_clusters=2, random_state=0)
        visualizer = KElbowVisualizer(
            km,
            k=list(range(2, effective_max_clusters + 1)),
            metric=metric,
            timings=False,
        )
        visualizer.fit(X_vals)
        visualizer.show()
        selected_clusters = visualizer.elbow_value_
        if selected_clusters is None:
            selected_clusters = effective_max_clusters
            warnings.warn(
                f"Failed to infer an elbow point automatically; falling back to {selected_clusters} clusters.",
                RuntimeWarning,
            )
    else:
        selected_clusters = effective_max_clusters
        if selected_clusters != n_neighborhoods:
            warnings.warn(
                f"`n_neighborhoods` was reduced from {n_neighborhoods} to {selected_clusters} "
                f"because only {len(X_vals)} neighborhoods are available.",
                RuntimeWarning,
            )

    km = MiniBatchKMeans(n_clusters=int(selected_clusters), random_state=0)
    labels = km.fit_predict(X_vals)

    if mode == "knn":
        neighborhood_name = f"CN_k{k}_n{int(selected_clusters)}"
        mode_value = k
    else:
        neighborhood_name = f"CN_r{radius:g}_n{int(selected_clusters)}"
        mode_value = radius

    result_df = df.copy()
    result_df[neighborhood_name] = labels

    metadata = {
        "mode": mode,
        "mode_value": mode_value,
        "selected_clusters": int(selected_clusters),
        "column_name": neighborhood_name,
        "metric": metric,
    }

    if return_metadata:
        return result_df, metadata
    return result_df


def hf_generate_random_colors(n: int, rand_seed: int = 1):
    np.random.seed(rand_seed)
    return sns.color_palette("hls", n)


def cn_exp_heatmap_df(
    df: pd.DataFrame,
    cluster_col: str,
    cn_col: str,
    palette=None,
    savefig: bool = False,
    output_fname: str = "",
    output_dir: str = "./",
    row_clus: bool = True,
    col_clus: bool = True,
    rand_seed: int = 1,
    figsize=(10, 5),
    vmin=None,
    vmax=None,
):
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if cluster_col not in df.columns or cn_col not in df.columns:
        raise KeyError(f"`{cluster_col}` and `{cn_col}` must both exist in the input DataFrame")

    if palette is None:
        unique_cns = np.sort(df[cn_col].unique())
        colors = hf_generate_random_colors(len(unique_cns), rand_seed=rand_seed)
        palette = dict(zip(unique_cns, colors))

    neigh_data = pd.DataFrame({cn_col: list(palette.keys()), "color": list(palette.values())})
    neigh_data.set_index(cn_col, inplace=True)

    cluster_dummies = pd.get_dummies(df[cluster_col])
    df3 = pd.concat([df, cluster_dummies], axis=1)
    sum_cols_sorted = sorted(cluster_dummies.columns.tolist())

    subset = df3[sum_cols_sorted + [cn_col]]
    niche_sub = subset.groupby(cn_col, observed=False).sum()
    niche_df = niche_sub.apply(lambda x: x / x.sum() * 10, axis=1)
    neigh_clusters = niche_df.to_numpy()

    tissue_avgs = df3[sum_cols_sorted].to_numpy().mean(axis=0)
    fc_2 = np.log2(
        (
            (neigh_clusters + tissue_avgs)
            / (neigh_clusters + tissue_avgs).sum(axis=1, keepdims=True)
        ) / tissue_avgs
    )
    fc_2 = pd.DataFrame(fc_2, columns=sum_cols_sorted, index=niche_df.index)

    if vmin is None or vmax is None:
        max_abs_val = np.nanmax(np.abs(fc_2.to_numpy()))
        if max_abs_val == 0:
            max_abs_val = 1
        calc_vmin = -max_abs_val if vmin is None else vmin
        calc_vmax = max_abs_val if vmax is None else vmax
    else:
        calc_vmin = vmin
        calc_vmax = vmax

    s = sns.clustermap(
        fc_2,
        vmin=calc_vmin,
        vmax=calc_vmax,
        cmap="bwr",
        figsize=figsize,
        row_colors=[neigh_data.reindex(fc_2.index)["color"]],
        cbar_pos=(0.03, 0.06, 0.03, 0.1),
    )

    s.ax_row_dendrogram.set_visible(row_clus)
    s.ax_col_dendrogram.set_visible(col_clus)
    s.ax_heatmap.set_ylabel("", labelpad=25)
    s.ax_heatmap.tick_params(axis="y", pad=42)
    s.ax_heatmap.yaxis.set_ticks_position("right")
    s.ax_heatmap.text(
        -0.25,
        -0.85,
        "log2 fold change\nover tissue average",
        transform=s.ax_heatmap.transAxes,
        fontsize=12,
        verticalalignment="bottom",
        horizontalalignment="center",
    )

    if savefig:
        base_name = output_fname or str(cn_col)
        s.figure.savefig(output_dir / f"{base_name}.pdf", bbox_inches="tight")
        s.figure.savefig(output_dir / f"{base_name}.png", bbox_inches="tight", dpi=300)

    return s
