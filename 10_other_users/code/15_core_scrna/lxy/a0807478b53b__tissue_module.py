import argparse
import os
import pickle
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scanpy as sc
import tensorly as tl
from tensorly.decomposition import non_negative_tucker

def build_tensors_by_group(adata, group_col='group', sample_col='sample',
                           cn_col='CNs_10', celltype_col='celltype',
                           cn_sort_key=None):
    if cn_sort_key is None:
        cn_sort_key = lambda x: int(str(x).replace('CN', '')) if str(x).startswith('CN') else x

    # 1. 提取全局唯一值（所有 group 共用相同的 CN 和 celltype 类别）
    cns_cats = sorted(adata.obs[cn_col].dropna().unique(), key=cn_sort_key)
    ctype_cats = sorted(adata.obs[celltype_col].dropna().unique())

    # 2. 确定分组
    if group_col is None:
        groups = {'All': sorted(adata.obs[sample_col].dropna().unique())}
    else:
        sample_group_map = (
            adata.obs[[sample_col, group_col]]
            .drop_duplicates()
            .dropna(subset=[sample_col, group_col])
        )
        groups = {}
        for grp in sorted(sample_group_map[group_col].unique()):
            samples = sorted(sample_group_map.loc[sample_group_map[group_col] == grp, sample_col].unique())
            groups[grp] = samples

    # 3. 为每个 group 构建张量
    dats = []
    labels = []
    patients = []
    for grp_name, grp_samples in groups.items():
        adata_grp = adata[adata.obs[sample_col].isin(grp_samples)]

        multi_idx = pd.MultiIndex.from_product(
            [grp_samples, cns_cats, ctype_cats],
            names=[sample_col, cn_col, celltype_col]
        )

        counts = adata_grp.obs.groupby([sample_col, cn_col, celltype_col]).size()
        counts_filled = counts.reindex(multi_idx, fill_value=0)

        group_sums = counts_filled.groupby(level=[sample_col, cn_col]).transform('sum')
        proportions = (counts_filled / group_sums).fillna(0)

        tensor_np = proportions.values.reshape(len(grp_samples), len(cns_cats), len(ctype_cats))

        dats.append(tensor_np)
        labels.append(str(grp_name))
        patients.append(list(grp_samples))
        print(f"Group '{grp_name}': {len(grp_samples)} samples -> tensor shape {tensor_np.shape}")

    return dats, labels, cns_cats, ctype_cats, patients

def evaluate_ranks(dat, num_tissue_modules=8, group_label="All", out_dir='./'):
    os.makedirs(out_dir, exist_ok=True)
    num_tissue_modules = num_tissue_modules + 1
    # pal = sns.color_palette('bright', 10)
    # palg = sns.color_palette('Greys', 10)
    
    mat1 = np.zeros((num_tissue_modules, 15))
    for i in range(2, 15):
        for j in range(1, num_tissue_modules):
            facs_overall = non_negative_tucker(dat, rank=[j, i, i], random_state=0)
            mat1[j, i] = np.mean((dat - tl.tucker_to_tensor((facs_overall[0], facs_overall[1])))**2)
            
    plt.figure()
    for j in range(1, num_tissue_modules):
        plt.plot(2 + np.arange(13), mat1[j][2:], label=f'rank = ({j},x,x)')
        
    plt.xlabel('x (CN/CT module rank)')
    plt.ylabel('reconstruction error')
    plt.title(f'Evaluate Ranks - Group: {group_label}')
    plt.legend()
    
    safe_label = "".join(c if c.isalnum() or c in ('-', '_') else '_' for c in group_label)
    save_path = os.path.join(out_dir, f"evaluate_ranks_{safe_label}")
    plt.savefig(f"{save_path}.pdf", bbox_inches='tight')
    plt.savefig(f"{save_path}.png", bbox_inches='tight', dpi=300)
    print(f"Saved rank evaluation plots for {group_label} to {out_dir}")
    plt.close()

def tissue_module_plots(dat, person_rank, rank, patients, nbs, cells_of_interest,
                        labels=None, random_state=0, out_dir='./'):
    os.makedirs(out_dir, exist_ok=True)
    
    if isinstance(dat, np.ndarray):
        dats = [dat]
        patients_list = [patients]
    else:
        dats = list(dat)
        if patients and isinstance(patients[0], (list, np.ndarray)):
            patients_list = patients
        else:
            patients_list = [patients] * len(dats)

    n_dats = len(dats)
    if labels is None:
        labels = [f'Group {i + 1}' for i in range(n_dats)]

    person_rank_list = [int(person_rank)] * n_dats if isinstance(person_rank, (int, np.integer)) else [int(v) for v in person_rank]
    rank_list = [int(rank)] * n_dats if isinstance(rank, (int, np.integer)) else [int(v) for v in rank]

    results = []

    for g, d in enumerate(dats):
        cur_pr = person_rank_list[g]
        cur_r = rank_list[g]
        cur_patients = patients_list[g]
        lbl = labels[g]
        safe_label = "".join(c if c.isalnum() or c in ('-', '_') else '_' for c in lbl)
        prefix = f'[{lbl}] ' if n_dats > 1 else ''

        print(f"{prefix}person_rank={cur_pr}, rank={cur_r}")
        facs_overall = non_negative_tucker(d, rank=[cur_pr, cur_r, cur_r], random_state=random_state)
        print(f"{prefix}core tensor shape: {facs_overall[0].shape}")

        def save_and_close(name):
            save_path = os.path.join(out_dir, f"{safe_label}_{name}")
            plt.savefig(f"{save_path}.pdf", bbox_inches='tight')
            plt.savefig(f"{save_path}.png", bbox_inches='tight', dpi=300)
            plt.close()

        # Optimize heatmap figure sizes dynamically based on matrix shapes
        # Add baseline margins to accommodate axis labels, colorbar, and titles properly.
        w_pr = max(cur_pr * 0.8 + 3, 4)
        w_r = max(cur_r * 0.8 + 3, 4)

        plt.figure(figsize=(w_pr, max(len(cur_patients) * 0.4 + 2, 3)))
        sns.heatmap(pd.DataFrame(facs_overall[1][0], index=cur_patients))
        plt.yticks(rotation=0)
        plt.ylabel('Donor')
        plt.xlabel('Donor module')
        plt.title('Donor modules')
        save_and_close("donor_modules")

        plt.figure(figsize=(w_r, max(len(nbs) * 0.4 + 2, 3)))
        sns.heatmap(pd.DataFrame(facs_overall[1][1], index=nbs))
        plt.yticks(rotation=0)
        plt.ylabel('CN')
        plt.xlabel('CN module')
        plt.title('CN modules')
        save_and_close("cn_modules")

        plt.figure(figsize=(w_r, max(len(cells_of_interest) * 0.4 + 2, 3)))
        sns.heatmap(pd.DataFrame(facs_overall[1][2], index=cells_of_interest))
        plt.yticks(rotation=0)
        plt.ylabel('CT')
        plt.xlabel('CT module')
        plt.title('CT modules')
        save_and_close("ct_modules")

        for i in range(cur_pr):
            plt.figure(figsize=(w_r, max(cur_r * 0.6 + 2, 3)))
            sns.heatmap(pd.DataFrame(facs_overall[0][i]))
            plt.yticks(rotation=0)
            plt.ylabel('CN module')
            plt.xlabel('CT module')
            plt.title(f'Tissue module {i}')
            save_and_close(f"tissue_module_{i}")

        results.append(facs_overall)

    return results if n_dats > 1 else results[0]

def tensor_plots(dats, ind, labels, cells, person_rank=2, rank=6, cells_plot=None,
                 tissue_modules=[], scale=1, savename=None, out_dir='./',
                 edge_threshold=0.1, show=False):
    pal = sns.color_palette('bright', 10)
    n_dats = len(dats)
    
    person_rank_list = [int(person_rank)] * n_dats if isinstance(person_rank, (int, np.integer)) else [int(v) for v in person_rank]
    rank_list = [int(rank)] * n_dats if isinstance(rank, (int, np.integer)) else [int(v) for v in rank]
    edge_threshold_list = [float(edge_threshold)] * n_dats if isinstance(edge_threshold, (int, float, np.integer, np.floating)) else [float(v) for v in edge_threshold]

    if cells_plot is None:
        cells_plot = cells

    nbs = list(range(len(ind)))
    offset = len(nbs)

    if not tissue_modules:
        tissue_modules = [[''] * pr for pr in person_rank_list]

    all_axes = []

    for dat_i, dat in enumerate(dats):
        cur_person_rank = person_rank_list[dat_i]
        cur_rank = rank_list[dat_i]
        cur_edge_threshold = edge_threshold_list[dat_i]
        print(f"Group '{labels[dat_i]}': person_rank={cur_person_rank}, rank={cur_rank}, edge_threshold={cur_edge_threshold}")
        facs_overall = non_negative_tucker(dat, rank=[cur_person_rank, cur_rank, cur_rank], random_state=0)

        for p in range(cur_person_rank):
            nb_scatter_size = scale * scale * 45
            cel_scatter_size = scale * scale * 45
            offset_all = 0
            idx2 = 0
            idx_idx2 = {}

            visible_count = cur_rank
            fig_h = (2.5 * visible_count + 6) / 4
            fig_h = fig_h if fig_h > 2 else 2

            fig, ax = plt.subplots(figsize=(len(cells_plot) / 8 * scale * 5, fig_h * scale))

            for idx in reversed(range(cur_rank)):
                idx_idx2[idx] = idx2

                max_row = np.max(facs_overall[0][p][idx, :])
                max_col = np.max(facs_overall[0][p][:, idx])
                an = 1.0 if max_row > cur_edge_threshold else 0.05
                ac = 1.0 if max_col > cur_edge_threshold else 0.05

                nb_fac = facs_overall[1][1][:, idx]
                cel_fac = facs_overall[1][2][:, idx]
                cel_fac = np.array([cel_fac[cells.index(cell)] for cell in cells_plot if cell in cells])

                cols_alpha = [(*pal[nb % 10], an * np.minimum(nb_fac, 1.0)[i]) for i, nb in enumerate(nbs)]
                cell_cols_alpha = [(0, 0, 0, ac * np.minimum(cel_fac, 1.0)[i]) for i, _ in enumerate(cel_fac)]

                an_alpha = float(np.clip(2 * an, 0.1, 1.0))
                ac_alpha = float(np.clip(2 * ac, 0.1, 1.0))

                ax.scatter(1 * np.arange(len(nb_fac)) + 0.5,
                           -2.5 * idx2 + np.zeros(len(nb_fac)) - offset_all,
                           c=cols_alpha, s=nb_scatter_size * 2,
                           edgecolors=len(cols_alpha) * [(0, 0, 0, an_alpha)], linewidths=0.5)
                ax.scatter(1 * np.arange(len(cel_fac)) + offset + 1.5,
                           -2.5 * idx2 + np.zeros(len(cel_fac)) - offset_all,
                           c=cell_cols_alpha, s=cel_scatter_size * 2, vmax=0.5,
                           edgecolors=len(cell_cols_alpha) * [(0, 0, 0, ac_alpha)], linewidths=0.5)

                rect = plt.Rectangle((0, -2.5 * idx2 - 1 - offset_all), len(nbs), 2,
                                     linewidth=scale * scale * 1, edgecolor='black', facecolor='none',
                                     zorder=0, alpha=an, linestyle='--', clip_on=False)
                ax.add_artist(rect)
                ax.scatter([offset], [-2.5 * idx2 - offset_all], c='black', marker='D',
                           s=scale * scale * 5, zorder=5, alpha=an)
                ax.scatter([offset + 1], [-2.5 * idx2 - offset_all], c='black', marker='D',
                           s=scale * scale * 5, zorder=5, alpha=ac)

                rect = plt.Rectangle((offset + 1, -2.5 * idx2 - 1 - offset_all), len(cells_plot), 2,
                                     linewidth=scale * scale * 1, edgecolor='black', facecolor='none',
                                     zorder=0, alpha=ac, linestyle='-.', clip_on=False)
                ax.add_artist(rect)

                ax.text(-0.25, -2.5 * idx2 - offset_all, f'{cur_rank - idx}',
                        ha='center', va='center', color='black', fontsize=scale * 12, alpha=1)
                idx2 += 1

            current_label = labels[dat_i] if dat_i < len(labels) else f"Dataset {dat_i}"
            safe_label = "".join(c if c.isalnum() or c in ('-', '_') else '_' for c in current_label)
            
            p_label = tissue_modules[dat_i][p] if p < len(tissue_modules[dat_i]) else str(p+1)
            ax.text((len(nbs) + len(cells_plot) + 1) / 2, 4 - offset_all,
                    f'{current_label} tissue module {p + 1}: {p_label}',
                    ha='center', va='center', color='black', fontsize=scale * 12, alpha=1)
            ax.text(len(nbs) / 2, 2 - offset_all, 'CN modules',
                    ha='center', va='center', color='black', fontsize=scale * 12, alpha=1)
            ax.text(offset + len(cells_plot) / 2 + 1, 2 - offset_all, 'CT modules',
                    ha='center', va='center', color='black', fontsize=scale * 12, alpha=1)
            rect = plt.Rectangle((-0.5, -0.5 - 2.5 * idx2 - offset_all + 1.5),
                                 len(nbs) + len(cells_plot) + 1.75, 2 + 2.5 * idx2,
                                 linewidth=scale * scale * 1, edgecolor='black', facecolor='none',
                                 zorder=0, alpha=1, clip_on=False)
            ax.add_artist(rect)

            for nb_i in range(cur_rank):
                for cel_i in range(cur_rank):
                    if nb_i in idx_idx2 and cel_i in idx_idx2 and facs_overall[0][p][nb_i, cel_i] > cur_edge_threshold:
                        val = facs_overall[0][p][nb_i, cel_i]
                        line_w = 2 * scale * scale * 1 * float(np.clip(val, 0.0, 1.0))
                        line_alpha = float(np.clip(10 * val, 0.0, 1.0))
                        ax.plot([offset, offset + 1],
                                [-2.5 * idx_idx2[nb_i] - offset_all, -2.5 * idx_idx2[cel_i] - offset_all],
                                color='black', linewidth=line_w, alpha=line_alpha)

            # Fix dynamic y-coordinate to correctly scale with rank, preventing overlap with boxes
            bottom_y = -2.5 * idx2 - offset_all + 0.5
            ax.text(-1, bottom_y, 'CN', fontsize=scale * 12, ha='center', va='top')
            for i, k in enumerate(nbs):
                ax.annotate(f'{k + 1}', (i + 0.5, bottom_y),
                            bbox={'boxstyle': 'circle', 'color': pal[k % 10], 'pad': 0.05, 'alpha': 0.8},
                            annotation_clip=False, size=scale * 12, ha='center', va='top')
            for i, chk in enumerate(cells_plot):
                ax.text(i + offset + 1.5, bottom_y, chk, rotation=90, color='black',
                        ha='center', va='top', fontsize=scale * 12, alpha=1)

            ax.set_ylim(bottom_y - 3, 6)
            ax.axis('off')

            if savename:
                os.makedirs(out_dir, exist_ok=True)
                base_name = os.path.splitext(savename)[0]
                save_base = f"{base_name}_{safe_label}_module{p + 1}" if n_dats > 1 else f"{base_name}_module{p + 1}"
                pdf_path = os.path.join(out_dir, f"{save_base}.pdf")
                png_path = os.path.join(out_dir, f"{save_base}.png")
                plt.savefig(pdf_path, bbox_inches='tight')
                plt.savefig(png_path, bbox_inches='tight', dpi=300)
                print(f"Saved tensor plot: {pdf_path}")

            all_axes.append(ax)
            if show:
                plt.show()
            else:
                plt.close(fig)

    return all_axes

def main():
    parser = argparse.ArgumentParser(description="Tissue Module Analysis Pipeline")
    parser.add_argument("--mode", type=str, required=True, choices=["evaluate", "decompose"],
                        help="Execution mode: 'evaluate' for rank evaluation without decomposition, 'decompose' for tensor decomposition and plots.")
    parser.add_argument("--data_input", type=str, default=None,
                        help="Input Scanpy .h5ad file. Required for evaluate mode.")
    parser.add_argument("--tensors_input", type=str, default=None,
                        help="Path to pre-computed tensors.pkl file for decompose mode if skipping h5ad loading.")
    parser.add_argument("--outdir", type=str, default="./results", help="Output directory to save results.")
    
    # Optional parameters for tensor construction
    parser.add_argument("--group_col", type=str, default="group", help="Column in adata.obs for grouping (e.g. 'condition').")
    parser.add_argument("--sample_col", type=str, default="sample", help="Column specifying samples.")
    parser.add_argument("--cn_col", type=str, default="CNs_10", help="Column specifying Cellular Neighborhoods.")
    parser.add_argument("--celltype_col", type=str, default="celltype", help="Column specifying cell types.")
    
    # Parameters for 'evaluate' mode
    parser.add_argument("--num_tissue_modules", type=int, default=8, help="Number of max tissue modules for rank evaluation.")
    
    # Parameters for 'decompose' mode
    parser.add_argument("--person_rank", type=str, default="2", help="Tucker decomposition rank 1 (tissue modules). Int or comma-separated.")
    parser.add_argument("--rank", type=str, default="6", help="Tucker decomposition rank 2 & 3 (CN/CT modules). Int or comma-separated.")
    parser.add_argument("--edge_threshold", type=str, default="0.1", help="Edge threshold masking. Float or comma-separated floats.")
    parser.add_argument("--tissue_modules_labels", type=str, default=None, help="Labels for tissue modules (e.g. 'Inflamed,Resting'), separated by comma.")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    
    cache_path = args.tensors_input if args.tensors_input else os.path.join(args.outdir, "tensors_data.pkl")

    if args.mode == "evaluate":
        if not args.data_input:
            raise ValueError("--data_input is required for 'evaluate' mode.")
        
        print(f"Loading {args.data_input} ...")
        adata = sc.read_h5ad(args.data_input)
        
        print("Building tensors by group...")
        dats, labels, cns_cats, ctype_cats, patients = build_tensors_by_group(
            adata, group_col=args.group_col, sample_col=args.sample_col,
            cn_col=args.cn_col, celltype_col=args.celltype_col
        )
        
        # Cache for next run
        print(f"Saving tensors cache to {cache_path} ...")
        with open(cache_path, "wb") as f:
            pickle.dump({
                "dats": dats,
                "labels": labels,
                "cns_cats": cns_cats,
                "ctype_cats": ctype_cats,
                "patients": patients
            }, f)
        
        eval_dir = os.path.join(args.outdir, "evaluate_ranks")
        print("Evaluating ranks...")
        for i, dat in enumerate(dats):
            evaluate_ranks(dat, num_tissue_modules=args.num_tissue_modules, group_label=labels[i], out_dir=eval_dir)
            
        print("Evaluate mode finished successfully! Next, run '--mode decompose' with chosen ranks.")

    elif args.mode == "decompose":
        if not os.path.exists(cache_path):
            raise FileNotFoundError(f"Tensors cache not found at {cache_path}. Please run evaluate mode first or specify --tensors_input.")
        
        print(f"Loading tensors from {cache_path} ...")
        with open(cache_path, "rb") as f:
            saved_data = pickle.load(f)
            dats = saved_data["dats"]
            labels = saved_data["labels"]
            cns_cats = saved_data["cns_cats"]
            ctype_cats = saved_data["ctype_cats"]
            patients = saved_data["patients"]
            
        # Parse list arguments
        def parse_arg(arg_str, typ):
            vals = [typ(x) for x in arg_str.split(",")]
            return vals if len(vals) > 1 else vals[0]
            
        pr = parse_arg(args.person_rank, int)
        rk = parse_arg(args.rank, int)
        th = parse_arg(args.edge_threshold, float)
        
        # Compute decomposition heatmaps
        heatmap_dir = os.path.join(args.outdir, "heatmaps")
        print("Generating Tissue Module Heatmaps...")
        results = tissue_module_plots(
            dat=dats, person_rank=pr, rank=rk, patients=patients,
            nbs=cns_cats, cells_of_interest=ctype_cats, labels=labels, out_dir=heatmap_dir
        )
        if not isinstance(results, list):
            results = [results]
            
        # Export CSVs
        csv_dir = os.path.join(args.outdir, "csv")
        print("Exporting decomposition matrices to CSV...")
        for g, res in enumerate(results):
            grp_csv_dir = os.path.join(csv_dir, labels[g])
            os.makedirs(grp_csv_dir, exist_ok=True)
            
            # extract and save
            core_tensor = res[0]
            person_rank, cn_rank, ct_rank = core_tensor.shape
            
            # core tensor
            for core_i in range(person_rank):
                cn_mod_names = [f"CN_module_{i+1}" for i in range(cn_rank)]
                ct_mod_names = [f"CT_module_{i+1}" for i in range(ct_rank)]
                pd.DataFrame(core_tensor[core_i], index=cn_mod_names, columns=ct_mod_names).to_csv(
                    os.path.join(grp_csv_dir, f"core_tensor_module_{core_i+1}.csv")
                )
            # factors
            donor_cols = [f"Donor_module_{i+1}" for i in range(person_rank)]
            cn_cols = [f"CN_module_{i+1}" for i in range(cn_rank)]
            ct_cols = [f"CT_module_{i+1}" for i in range(ct_rank)]
            
            pd.DataFrame(res[1][0], index=patients[g], columns=donor_cols).to_csv(os.path.join(grp_csv_dir, "donor_modules.csv"))
            pd.DataFrame(res[1][1], index=cns_cats, columns=cn_cols).to_csv(os.path.join(grp_csv_dir, "cn_modules.csv"))
            pd.DataFrame(res[1][2], index=ctype_cats, columns=ct_cols).to_csv(os.path.join(grp_csv_dir, "ct_modules.csv"))

        # Tissue Module Plots (Network representation)
        tensor_dir = os.path.join(args.outdir, "tensor_plots")
        print("Generating Full Tensor Plots...")
        tm_labels = [args.tissue_modules_labels.split(",")] if args.tissue_modules_labels else []
        
        tensor_plots(
            dats=dats, ind=cns_cats, labels=labels, cells=ctype_cats,
            person_rank=pr, rank=rk, cells_plot=ctype_cats,
            tissue_modules=tm_labels, scale=1, savename="tensor", out_dir=tensor_dir,
            edge_threshold=th, show=False
        )
        print("Decompose mode finished successfully!")

if __name__ == "__main__":
    main()