from pathlib import Path
import argparse
import anndata as ad


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default=r"d:\develop\TRAE\SP_pipe\FWSP20255000\Macrophages\adata_cluster.h5ad",
    )
    parser.add_argument("--groupby", default="res_0.3")
    parser.add_argument(
        "--output",
        default=r"d:\develop\TRAE\SP_pipe\FWSP20255000\Macrophages\res_0_3_cluster_cell_counts.tsv",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    adata = ad.read_h5ad(input_path)

    if args.groupby not in adata.obs.columns:
        raise ValueError(f"{args.groupby} 不在 adata.obs 中，可用列: {list(adata.obs.columns)}")

    counts = adata.obs[args.groupby].astype(str).value_counts(sort=False)
    counts = counts.sort_index()
    result = counts.rename_axis("cluster").reset_index(name="cell_count")
    result["cell_ratio"] = result["cell_count"] / int(result["cell_count"].sum())

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)
    print(f"已输出: {output_path}")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()

