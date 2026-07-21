import os
import sys

# ---------- 第1步：必须在任何 torch 相关导入前屏蔽 GPU ----------
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["PYTORCH_NO_CUDA_MEMORY_CACHING"] = "1"

import torch
torch.set_default_device('cpu')
torch.set_default_tensor_type(torch.FloatTensor)

# ---------- 第2步：正常导入 ----------
import argparse
import logging
import time
import pickle
from pathlib import Path

import pandas as pd
import scanpy as sc
from rctd import RCTDConfig, Reference, run_rctd


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run RCTD deconvolution on spatial transcriptomics data."
    )
    parser.add_argument("--ref_path", type=str, required=True)
    parser.add_argument("--spatial_path", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--celltype_col", type=str, default="Celltype")
    parser.add_argument(
        "--mode",
        type=str,
        default="doublet",
        choices=["doublet", "full", "multi"],
    )
    parser.add_argument("--batch_size", type=int, default=10000,
                        help="Batch size for GPU processing.")
    parser.add_argument("--sigma_override", type=int, default=None,
                        help="Override platform-effect sigma with a known value.")
    # RCTDConfig 参数
    parser.add_argument("--gene_cutoff", type=float, default=0.000125)
    parser.add_argument("--fc_cutoff", type=float, default=0.5)
    parser.add_argument("--gene_cutoff_reg", type=float, default=0.0002)
    parser.add_argument("--fc_cutoff_reg", type=float, default=0.75)
    parser.add_argument("--UMI_min", type=int, default=100)
    parser.add_argument("--UMI_max", type=int, default=20_000_000)
    parser.add_argument("--counts_MIN", type=int, default=10)
    parser.add_argument("--UMI_min_sigma", type=int, default=300)
    parser.add_argument("--MIN_OBS", type=int, default=3)
    parser.add_argument("--CONFIDENCE_THRESHOLD", type=float, default=5.0)
    parser.add_argument("--DOUBLET_THRESHOLD", type=float, default=20.0)
    parser.add_argument("--max_iter", type=int, default=50)
    parser.add_argument("--step_size", type=float, default=0.3)
    return parser.parse_args()


def safe_array(arr, expected_length, default=0):
    """
    将某个可能为 None 或长度不足的数组/列表补全到 expected_length。
    若为 None，则返回全部默认值的列表。
    """
    if arr is None:
        return [default] * expected_length
    arr = list(arr)
    if len(arr) < expected_length:
        arr += [default] * (expected_length - len(arr))
    elif len(arr) > expected_length:
        arr = arr[:expected_length]   # 理论上不应发生，但做个截断
    return arr


def main() -> None:
    setup_logging()
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    os.chdir(args.output_dir)
    logging.info("Output directory: %s", os.getcwd())

    sample_id = Path(args.spatial_path).stem
    logging.info("Processing sample: %s", sample_id)

    # 加载参考数据
    logging.info("Loading reference data from %s", args.ref_path)
    ref_adata = sc.read_h5ad(args.ref_path)
    ref_adata.X = ref_adata.X.astype(int)
    logging.info("Reference cell types: %s",
                 ref_adata.obs[args.celltype_col].unique().tolist())

    # 加载空间数据
    logging.info("Loading spatial data from %s", args.spatial_path)
    spat_adata = sc.read_h5ad(args.spatial_path)
    spat_adata.X = spat_adata.X.astype(int)

    # 基因交集
    logging.info("Intersecting genes...")
    common_genes = ref_adata.var_names.intersection(spat_adata.var_names)
    ref_adata = ref_adata[:, common_genes].copy()
    spat_adata = spat_adata[:, common_genes].copy()
    logging.info("Common genes after intersection: %d", len(common_genes))

    # 构建 Reference
    reference = Reference(ref_adata, cell_type_col=args.celltype_col)
    logging.info("Reference profiles shape: %s", reference.profiles.shape)
    logging.info("Cell type names: %s", reference.cell_type_names)

    # 创建 RCTDConfig
    config = RCTDConfig(
        gene_cutoff=args.gene_cutoff,
        fc_cutoff=args.fc_cutoff,
        gene_cutoff_reg=args.gene_cutoff_reg,
        fc_cutoff_reg=args.fc_cutoff_reg,
        UMI_min=args.UMI_min,
        UMI_max=args.UMI_max,
        counts_MIN=args.counts_MIN,
        UMI_min_sigma=args.UMI_min_sigma,
        MIN_OBS=args.MIN_OBS,
        CONFIDENCE_THRESHOLD=args.CONFIDENCE_THRESHOLD,
        DOUBLET_THRESHOLD=args.DOUBLET_THRESHOLD,
        max_iter=args.max_iter,
        step_size=args.step_size,
    )

    # ========== 运行 RCTD ==========
    logging.info("Starting RCTD in %s mode (CPU only)...", args.mode)
    t0 = time.time()
    result = run_rctd(
        spat_adata,
        reference,
        mode=args.mode,
        config=config,
        batch_size=args.batch_size,
        sigma_override=args.sigma_override,
    )
    elapsed = time.time() - t0
    logging.info("RCTD finished in %.1f seconds", elapsed)

    # ========== 立即保存原始 result（防止后续构建 DataFrame 失败导致重跑）==========
    raw_result_path = Path(args.output_dir) / f"{sample_id}_rctd_result.pkl"
    with open(raw_result_path, 'wb') as f:
        pickle.dump({
            'result': result,
            'cell_type_names': result.cell_type_names,
            'obs_names': list(spat_adata.obs_names)
        }, f, protocol=pickle.HIGHEST_PROTOCOL)
    logging.info("Raw result saved to %s", raw_result_path)

    # ========== 构建最终 CSV ==========
    try:
        cell_types = result.cell_type_names
        n_spots = spat_adata.n_obs

        if args.mode == "doublet":
            # 安全获取各个数组
            spot_class = safe_array(result.spot_class, n_spots, default=0)
            first_type = safe_array(result.first_type, n_spots, default=-1)
            second_type = safe_array(result.second_type, n_spots, default=-1)
            first_class = safe_array(result.first_class, n_spots, default=float('nan'))
            second_class = safe_array(result.second_class, n_spots, default=float('nan'))
            min_score = safe_array(result.min_score, n_spots, default=0.0)
            singlet_score = safe_array(result.singlet_score, n_spots, default=0.0)

            class_labels = {0: "reject", 1: "singlet", 2: "doublet_certain", 3: "doublet_uncertain"}
            df = pd.DataFrame({
                "barcode": spat_adata.obs_names,
                "spot_class": [class_labels.get(c, "unknown") for c in spot_class],
                "first_type": [cell_types[t] if t != -1 else "NA" for t in first_type],
                "second_type": [cell_types[t] if t != -1 else "NA" for t in second_type],
                "first_class": first_class,
                "second_class": second_class,
                "min_score": min_score,
                "singlet_score": singlet_score,
            })

        elif args.mode == "full":
            # full 模式通常长度一致，但这里也加一层保护
            weights = result.weights
            if weights.shape[0] != n_spots:
                logging.warning("Weight matrix rows (%d) != number of spots (%d), using row names if possible.",
                                weights.shape[0], n_spots)
                # 如果权重矩阵行数与 spots 不符，尝试用共同索引
                if hasattr(result, 'barcodes'):
                    common_bcs = [b for b in spat_adata.obs_names if b in result.barcodes]
                    df = pd.DataFrame(weights, index=result.barcodes, columns=cell_types)
                    df = df.loc[common_bcs]
                    df.index.name = "barcode"
                else:
                    raise ValueError("Mismatch between weights rows and spot number.")
            else:
                df = pd.DataFrame(weights, index=spat_adata.obs_names, columns=cell_types)
                df.index.name = "barcode"

        elif args.mode == "multi":
            n_types = result.n_types
            max_types = result.weights.shape[1]
            # n_types 应与 spot 数量一致
            if len(n_types) != n_spots:
                n_types = safe_array(n_types, n_spots, default=0)
            records = []
            for i, barcode in enumerate(spat_adata.obs_names):
                cur_n = n_types[i]
                if cur_n == 0:
                    records.append({"barcode": barcode, "n_types": 0})
                    continue
                type_indices = result.cell_type_indices[i, :cur_n]
                weights_arr = result.weights[i, :cur_n]
                type_names = [cell_types[idx] for idx in type_indices]
                row = {"barcode": barcode, "n_types": cur_n}
                for j in range(max_types):
                    row[f"type_{j+1}"] = type_names[j] if j < len(type_names) else ""
                    row[f"weight_{j+1}"] = weights_arr[j] if j < len(weights_arr) else 0.0
                records.append(row)
            df = pd.DataFrame(records)

        # 保存 CSV
        output_file = f"{sample_id}_RCTD.csv"
        df.to_csv(output_file, index=False)
        logging.info("Saved result CSV to %s", output_file)

    except Exception as e:
        logging.error("Failed to build final CSV: %s", e)
        logging.info("Raw result is available at %s, you can process it manually.", raw_result_path)
        raise

if __name__ == "__main__":
    main()
