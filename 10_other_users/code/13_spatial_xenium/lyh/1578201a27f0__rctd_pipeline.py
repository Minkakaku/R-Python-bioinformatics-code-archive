#!/usr/bin/env python3
import argparse
import gc
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rctd
import scanpy as sc
import scipy.sparse as sp
import torch
from rctd import Reference, run_rctd

import rctd._likelihood
import rctd._rctd

Q_MATRICES_PATH = "/mnt/NexData_file/home/lyh/Software/rctdpy_need/q_matrices.npz"


def log(message):
    print(f"[{time.strftime('%F %T')}] {message}", flush=True)


def matrix_count_check(X, label, max_values=1_000_000, tol=1e-6):
    if sp.issparse(X):
        values = X.data
    else:
        arr = np.asarray(X)
        values = arr.ravel()

    if values.size == 0:
        log(f"{label}: matrix has no stored nonzero values; treating as count-compatible")
        return

    if values.size > max_values:
        rng = np.random.default_rng(0)
        values = values[rng.choice(values.size, size=max_values, replace=False)]

    values = np.asarray(values)
    finite = np.isfinite(values)
    if not finite.all():
        raise ValueError(f"{label}.X contains non-finite values; expected raw count matrix")

    min_value = float(values.min())
    max_fractional = float(np.max(np.abs(values - np.round(values))))
    if min_value < -tol:
        raise ValueError(f"{label}.X contains negative values; expected raw count matrix")
    if max_fractional > tol:
        raise ValueError(
            f"{label}.X does not look like integer counts; max fractional deviation={max_fractional:.3g}"
        )

    log(f"{label}.X count check passed on {values.size} sampled/stored values")


def patch_q_matrices(local_npz_path):
    def load_local_q_matrices(*args, **kwargs):
        log(f"Loading local q_matrices: {local_npz_path}")
        return dict(np.load(local_npz_path))

    rctd._rctd.load_cached_q_matrices = load_local_q_matrices
    rctd._likelihood.load_cached_q_matrices = load_local_q_matrices


def patch_run_full_mode_low_vram():
    from rctd._irwls import solve_irwls_batch_shared
    from rctd._types import FullResult, resolve_device
    import rctd._doublet as rctd_doublet
    import rctd._full as rctd_full
    import rctd._multi as rctd_multi

    def run_full_mode_low_vram(
        spatial_counts,
        spatial_numi,
        norm_profiles,
        cell_type_names,
        q_mat,
        sq_mat,
        x_vals,
        batch_size=1000,
        device="auto",
        pixel_mask=None,
    ):
        N, _ = spatial_counts.shape
        device_obj = resolve_device(device)
        log(f"low_vram full fit: N={N}, batch_size={batch_size}, device={device_obj}")

        P_gpu = torch.as_tensor(norm_profiles, device=device_obj)
        Q_gpu = torch.as_tensor(q_mat, device=device_obj)
        SQ_gpu = torch.as_tensor(sq_mat, device=device_obj)
        X_gpu = torch.as_tensor(x_vals, device=device_obj)

        all_weights = []
        all_converged = []
        for start in range(0, N, batch_size):
            end = min(start + batch_size, N)
            if start == 0 or end == N or (start // batch_size) % 10 == 0:
                log(f"full fit batch {start}:{end}/{N}")

            Y_batch = torch.as_tensor(spatial_counts[start:end], device=device_obj)
            nUMI_batch = torch.as_tensor(spatial_numi[start:end], device=device_obj)
            weights, converged = solve_irwls_batch_shared(
                P=P_gpu,
                Y_batch=Y_batch,
                nUMI_batch=nUMI_batch,
                Q_mat=Q_gpu,
                SQ_mat=SQ_gpu,
                x_vals=X_gpu,
                max_iter=50,
                min_change=0.001,
                constrain=False,
                bulk_mode=False,
            )
            all_weights.append(weights.cpu().numpy())
            all_converged.append(converged.cpu().numpy())

            del Y_batch, nUMI_batch, weights, converged
            if device_obj.type == "cuda":
                torch.cuda.empty_cache()

        return FullResult(
            weights=np.vstack(all_weights),
            cell_type_names=cell_type_names,
            converged=np.concatenate(all_converged),
            pixel_mask=pixel_mask,
        )

    rctd_full.run_full_mode = run_full_mode_low_vram
    rctd_multi.run_full_mode = run_full_mode_low_vram
    rctd_doublet.run_full_mode = run_full_mode_low_vram


def fitted_result_dataframe(adata, result, mode):
    fit_ids = adata.obs_names[result.pixel_mask].astype(str)
    cell_type_names = np.asarray(result.cell_type_names, dtype=object)
    dominant_idx = result.weights.argmax(axis=1)

    df = pd.DataFrame(index=fit_ids)
    df.index.name = "cell_id"
    df["rctd_mode"] = mode
    df["rctd_dominant_type"] = cell_type_names[dominant_idx]
    df["rctd_dominant_weight"] = result.weights.max(axis=1)

    if mode == "full":
        df["rctd_converged"] = result.converged

    elif mode == "doublet":
        spot_class_names = np.array(["reject", "singlet", "doublet_certain", "doublet_uncertain"], dtype=object)
        df["rctd_spot_class"] = spot_class_names[result.spot_class]
        df["rctd_first_type"] = cell_type_names[result.first_type]
        df["rctd_second_type"] = cell_type_names[result.second_type]
        df["rctd_first_weight"] = result.weights_doublet[:, 0]
        df["rctd_second_weight"] = result.weights_doublet[:, 1]
        df["rctd_min_score"] = result.min_score
        df["rctd_singlet_score"] = result.singlet_score
        df["rctd_first_class"] = result.first_class
        df["rctd_second_class"] = result.second_class
        if result.first_class_name is not None:
            df["rctd_first_class_name"] = result.first_class_name
        if result.second_class_name is not None:
            df["rctd_second_class_name"] = result.second_class_name

    elif mode == "multi":
        for i in range(result.sub_weights.shape[1]):
            idx = result.cell_type_indices[:, i].astype(int)
            df[f"rctd_multi_type_{i + 1}"] = np.where(idx >= 0, cell_type_names[idx], pd.NA)
            df[f"rctd_multi_weight_{i + 1}"] = result.sub_weights[:, i]
            df[f"rctd_multi_confident_{i + 1}"] = result.conf_list[:, i]
        df["rctd_n_types"] = result.n_types
        df["rctd_min_score"] = result.min_score
        if result.sub_weights.shape[1] >= 1:
            df["rctd_multi_type"] = df["rctd_multi_type_1"]
            df["rctd_multi_weight"] = df["rctd_multi_weight_1"]
            df["rctd_confident"] = df["rctd_multi_confident_1"]

    for i, cell_type in enumerate(result.cell_type_names):
        df[f"rctd_weight__{cell_type}"] = result.weights[:, i]
    return df




def sanitize_for_obs_write(df):
    clean = pd.DataFrame(index=df.index)
    for col in df.columns:
        series = df[col]
        non_na = series.dropna()
        unique_values = set(non_na.unique().tolist()) if len(non_na) else set()
        if pd.api.types.is_bool_dtype(series) or unique_values <= {True, False}:
            clean[col] = series.map({True: 1.0, False: 0.0})
        elif pd.api.types.is_numeric_dtype(series):
            clean[col] = series
        else:
            clean[col] = series.astype("string").fillna("NA").astype(str)
    return clean

def expand_to_obs(adata, fitted_df):
    obs_df = pd.DataFrame(index=adata.obs_names)
    for col in fitted_df.columns:
        obs_df[col] = fitted_df[col].reindex(adata.obs_names)
    return obs_df


def main():
    parser = argparse.ArgumentParser(description="Run low-VRAM rctd-py on processed spatial/ref h5ad files and save outputs under a working directory.")
    parser.add_argument("--spatial", required=True, help="Processed spatial h5ad with raw count matrix in .X")
    parser.add_argument("--reference", required=True, help="Reference h5ad with raw count matrix in .X")
    parser.add_argument("--cell-type-col", required=True, help="Column in reference.obs containing cell type labels")
    parser.add_argument("--workdir", required=True, help="Working directory for RCTD outputs. Created if missing.")
    parser.add_argument("--mode", default="multi", choices=["full", "doublet", "multi"])
    parser.add_argument("--max-multi-types", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=2000)
    parser.add_argument("--umi-min", type=int, default=100)
    parser.add_argument("--counts-min", type=int, default=10)
    parser.add_argument("--device", default="cuda", choices=["cuda", "cpu", "auto"])
    parser.add_argument("--skip-count-check", action="store_true")
    args = parser.parse_args()

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    output_csv = workdir / f"{Path(args.spatial).stem}.rctd_{args.mode}.csv"
    output_h5ad = workdir / f"{Path(args.spatial).stem}.rctd_{args.mode}.h5ad"

    log(f"Python={os.sys.executable}")
    log(f"rctd={getattr(rctd, '__version__', 'unknown')} file={Path(rctd.__file__).resolve()}")
    log(f"torch={torch.__version__}, cuda={torch.version.cuda}, cuda_available={torch.cuda.is_available()}")

    patch_q_matrices(Q_MATRICES_PATH)
    patch_run_full_mode_low_vram()

    log(f"Reading spatial: {args.spatial}")
    spatial = sc.read_h5ad(args.spatial)
    spatial.var_names_make_unique()
    spatial.obs_names_make_unique()
    log(f"Spatial shape={spatial.shape}")

    log(f"Reading reference: {args.reference}")
    ref = sc.read_h5ad(args.reference)
    ref.var_names_make_unique()
    ref.obs_names_make_unique()
    log(f"Reference shape={ref.shape}")

    if args.cell_type_col not in ref.obs.columns:
        raise ValueError(f"--cell-type-col '{args.cell_type_col}' not found in reference.obs")

    if not args.skip_count_check:
        matrix_count_check(spatial.X, "spatial")
        matrix_count_check(ref.X, "reference")

    reference = Reference(ref, cell_type_col=args.cell_type_col)
    log(f"Reference cell types={reference.n_types}: {reference.cell_type_names}")

    config = rctd.RCTDConfig(
        device=args.device,
        MAX_MULTI_TYPES=args.max_multi_types,
        dtype="float32",
        compile=False,
        UMI_min=args.umi_min,
        counts_MIN=args.counts_min,
        DOUBLET_THRESHOLD=20.0,
        CONFIDENCE_THRESHOLD=5.0,
    )

    start = time.time()
    result = run_rctd(spatial, reference, mode=args.mode, config=config, batch_size=args.batch_size)
    log(f"RCTD {args.mode} finished in {time.time() - start:.1f}s")

    fitted_df = fitted_result_dataframe(spatial, result, args.mode)
    obs_result = expand_to_obs(spatial, fitted_df)

    csv_df = obs_result.copy()
    csv_df.insert(0, "cell_id", csv_df.index.astype(str))
    csv_df.to_csv(output_csv, index=False)
    log(f"Wrote CSV: {output_csv} rows={csv_df.shape[0]} cols={csv_df.shape[1]}")

    obs_for_h5ad = sanitize_for_obs_write(obs_result)
    for col in obs_for_h5ad.columns:
        spatial.obs[col] = obs_for_h5ad[col]
    spatial.uns["rctd_mode"] = args.mode
    spatial.uns["rctd_cell_type_col"] = args.cell_type_col
    spatial.uns["rctd_cell_type_names"] = list(result.cell_type_names)

    spatial.write_h5ad(output_h5ad)
    log(f"Wrote h5ad with RCTD obs columns: {output_h5ad}")

    del result, spatial, ref, reference
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
