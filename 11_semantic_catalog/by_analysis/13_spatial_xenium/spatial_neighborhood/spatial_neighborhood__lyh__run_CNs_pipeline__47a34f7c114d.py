import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CN_MODULE_DIR = PROJECT_ROOT / "CNs_V2"
if str(CN_MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(CN_MODULE_DIR))

from CNs import determine_optimal_cluster_count, neighborhood_analysis_df  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(description="Run single-mode CN analysis and export metadata results.")
    parser.add_argument("--meta_path", required=True, help="Path to the metadata CSV used for CN analysis.")
    parser.add_argument("--output_dir", required=True, help="Directory used to store CN outputs.")
    parser.add_argument("--index_col", default="0", help="CSV index column passed to pandas. Use 'none' to disable.")
    parser.add_argument("--region_col", default="sample", help="Column identifying independent regions.")
    parser.add_argument("--cluster_col", default="predicted.celltype", help="Column containing cell type labels.")
    parser.add_argument("--x_col", default="X", help="X coordinate column.")
    parser.add_argument("--y_col", default="Y", help="Y coordinate column.")
    parser.add_argument("--mode", default="radius", help="CN mode: knn or radius. 'dadius' is also accepted as radius.")
    parser.add_argument("--radius", type=float, default=80, help="Radius used for radius mode.")
    parser.add_argument("--knn_k", type=int, default=40, help="k used for knn mode.")
    parser.add_argument("--n_clusters", default="auto", help="Cluster count. Use an integer or 'auto' for elbow-based selection.")
    parser.add_argument("--elbow_max", type=int, default=20, help="Max clusters tested when n_clusters=auto.")
    parser.add_argument("--metric", default="distortion", help="Metric passed to Yellowbrick elbow selection.")
    return parser.parse_args()


def resolve_index_col(index_col_value: str):
    if isinstance(index_col_value, str) and index_col_value.lower() == "none":
        return None
    try:
        return int(index_col_value)
    except (TypeError, ValueError):
        return index_col_value


def normalize_mode(mode: str) -> str:
    normalized = str(mode).lower()
    if normalized == "dadius":
        normalized = "radius"
    if normalized not in {"knn", "radius"}:
        raise ValueError(f"Unsupported mode: {mode}. Use 'knn' or 'radius'.")
    return normalized


def resolve_clusters(n_clusters_value: str):
    if str(n_clusters_value).lower() == "auto":
        return None
    try:
        resolved = int(n_clusters_value)
    except ValueError as exc:
        raise ValueError("--n_clusters must be an integer or 'auto'") from exc
    if resolved <= 0:
        raise ValueError("--n_clusters must be a positive integer")
    return resolved


def save_elbow_plot(output_dir: Path, mode: str, radius: float, knn_k: int):
    if mode == "radius":
        prefix = f"radius_r{radius:g}"
    else:
        prefix = f"knn_k{knn_k}"
    plt.savefig(output_dir / f"{prefix}_elbow.pdf", bbox_inches="tight")
    plt.savefig(output_dir / f"{prefix}_elbow.png", bbox_inches="tight", dpi=300)
    plt.close("all")


def main():
    args = parse_args()
    meta_path = Path(args.meta_path).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not meta_path.exists():
        raise FileNotFoundError(f"Metadata CSV not found: {meta_path}")

    mode = normalize_mode(args.mode)
    manual_clusters = resolve_clusters(args.n_clusters)
    index_col = resolve_index_col(args.index_col)

    print(f"Loading metadata from: {meta_path}")
    df = pd.read_csv(meta_path, index_col=index_col)
    print(f"Loaded {df.shape[0]} rows and {df.shape[1]} columns")

    required_columns = [args.region_col, args.cluster_col, args.x_col, args.y_col]
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise KeyError(f"Missing required columns in metadata: {missing_columns}")

    if manual_clusters is None:
        print(f"Running elbow-based cluster selection for {mode} mode...")
        elbow_info = determine_optimal_cluster_count(
            df,
            region_col=args.region_col,
            cluster_col=args.cluster_col,
            X=args.x_col,
            Y=args.y_col,
            mode=mode,
            k=args.knn_k,
            radius=args.radius,
            max_clusters=args.elbow_max,
            metric=args.metric,
        )
        selected_clusters = int(elbow_info["selected_clusters"])
        save_elbow_plot(output_dir, mode=mode, radius=args.radius, knn_k=args.knn_k)
        selection_method = "auto"
    else:
        selected_clusters = manual_clusters
        elbow_info = None
        selection_method = "manual"

    print(f"Running {mode} CN analysis with {selected_clusters} clusters...")
    df, run_info = neighborhood_analysis_df(
        df,
        region_col=args.region_col,
        cluster_col=args.cluster_col,
        mode=mode,
        elbow=False,
        k=args.knn_k,
        radius=args.radius,
        n_neighborhoods=selected_clusters,
        X=args.x_col,
        Y=args.y_col,
        metric=args.metric,
        return_metadata=True,
    )

    output_csv = output_dir / "add_CNs.csv"
    df.to_csv(output_csv)

    config = {
        "meta_path": str(meta_path),
        "output_csv": str(output_csv),
        "mode": mode,
        "region_col": args.region_col,
        "cluster_col": args.cluster_col,
        "x_col": args.x_col,
        "y_col": args.y_col,
        "radius": args.radius,
        "knn_k": args.knn_k,
        "cluster_selection": selection_method,
        "requested_n_clusters": args.n_clusters,
        "selected_n_clusters": run_info["selected_clusters"],
        "cn_column": run_info["column_name"],
        "metric": args.metric,
        "elbow": elbow_info,
    }
    (output_dir / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    print(f"Saved CN metadata to: {output_csv}")
    print(f"Selected CN column: {run_info['column_name']}")
    print(f"Finished. Outputs stored in: {output_dir}")


if __name__ == "__main__":
    main()
