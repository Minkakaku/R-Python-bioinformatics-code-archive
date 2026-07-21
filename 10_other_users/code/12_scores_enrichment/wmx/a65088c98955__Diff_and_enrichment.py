import os
import scanpy as sc
import subprocess
import pandas as pd

rscript = "/mnt/NexData_file/share/scRNA-tools/envs/pip_diffgene_enrich/bin/Rscript"
r_code = "/mnt/NexData_file/home/lst/scrna-code/diff/do_enrichment.R"


# --- 1. Cell Type 之间的 Marker 分析 ---
def find_marker_gene(adata, group_col='cell_type', out_dir=os.getcwd()):
    print(f"[INFO] Starting Marker gene analysis for {group_col}...")

    counts = adata.obs[group_col].value_counts()
    keep = counts[counts > 1].index
    adata_sub = adata[adata.obs[group_col].isin(keep)].copy()

    if adata_sub.n_obs == 0:
        print(f"[SKIP] No cells left for marker analysis.")
        return

    sc.tl.rank_genes_groups(
        adata_sub,
        group_col,
        method='wilcoxon',
        pts=True
    )

    groups = adata_sub.uns['rank_genes_groups']['names'].dtype.names

    for c in groups:
        target_dir = os.path.join(out_dir, 'cell_type', str(c).replace('/', '_').replace(' ', '_'))
        os.makedirs(target_dir, exist_ok=True)

        df = sc.get.rank_genes_groups_df(adata_sub, group=c)

        df_filter = df[
            (df['logfoldchanges'] > 0.25) &
            (df['pvals_adj'] < 0.05) &
            (df['pct_nz_group'] >= 0.1)
        ]

        df_filter.to_csv(
            os.path.join(target_dir, 'differentially_expressed_genes.xls'),
            sep='\t',
            index=None
        )

        print(f"[SAVE] cell_type marker: {c}, genes = {df_filter.shape[0]}")


# --- 2. 全部细胞 + 每个 Cell Type 内：各组 vs K ---
def find_group_degs(adata, cell_type_col='cell_type', cond_col='group', out_dir=os.getcwd()):
    print(f"[INFO] Starting Group DE analysis: each group vs K...")

    control_group = 'K'

    group_base_dir = os.path.join(out_dir, 'group')
    os.makedirs(group_base_dir, exist_ok=True)

    adata.obs[cond_col] = adata.obs[cond_col].astype(str)
    adata.obs[cell_type_col] = adata.obs[cell_type_col].astype(str)

    # 增加 All：全部细胞的 6 组比较
    cell_types = ['All'] + list(adata.obs[cell_type_col].dropna().unique())

    for ct in cell_types:

        if ct == 'All':
            adata_ct = adata.copy()
            ct_name = 'All'
        else:
            adata_ct = adata[adata.obs[cell_type_col] == ct].copy()
            ct_name = str(ct).replace('/', '_').replace(' ', '_')

        groups = [g for g in adata_ct.obs[cond_col].unique() if pd.notna(g)]

        if control_group not in groups:
            print(f"[SKIP] {ct}: control group {control_group} not found.")
            continue

        n_control = (adata_ct.obs[cond_col] == control_group).sum()

        if n_control < 3:
            print(f"[SKIP] {ct}: control group {control_group} cells too small: {n_control}")
            continue

        for g in groups:

            if g == control_group:
                continue

            n_g = (adata_ct.obs[cond_col] == g).sum()

            if n_g < 3:
                print(
                    f"[SKIP] {ct} ({g} vs {control_group}): "
                    f"Sample size too small ({n_g} vs {n_control})"
                )
                continue

            g_name = str(g).replace('/', '_').replace(' ', '_')
            control_name = str(control_group).replace('/', '_').replace(' ', '_')

            comparison_name = f"{ct_name}_{g_name}_vs_{control_name}"
            target_dir = os.path.join(group_base_dir, comparison_name)
            os.makedirs(target_dir, exist_ok=True)

            print(f"[RUN] Comparing {comparison_name}")

            adata_tmp = adata_ct[
                adata_ct.obs[cond_col].isin([g, control_group])
            ].copy()

            adata_tmp.obs[cond_col] = adata_tmp.obs[cond_col].astype('category')

            sc.tl.rank_genes_groups(
                adata_tmp,
                cond_col,
                groups=[g],
                reference=control_group,
                method='wilcoxon',
                pts=True
            )

            df = sc.get.rank_genes_groups_df(adata_tmp, group=g)

            # 保持原代码逻辑：只保存一个 DEG 文件，然后用于富集
            df_filter = df[
                (df['logfoldchanges'].abs() > 0.25) &
                (df['pvals_adj'] < 0.05)
            ]

            df_filter.to_csv(
                os.path.join(target_dir, 'differentially_expressed_genes.xls'),
                sep='\t',
                index=None
            )

            print(f"[SAVE] {comparison_name}: {df_filter.shape[0]} genes")

    print(f"[INFO] DE analysis finished.")


# --- 3. 递归运行富集分析 ---
def run_enrichment_recursive(base_dir, species):
    for root, dirs, files in os.walk(base_dir):
        if 'differentially_expressed_genes.xls' in files:
            deg_file = os.path.join(root, 'differentially_expressed_genes.xls')

            try:
                df = pd.read_csv(deg_file, sep='\t')
                if df.shape[0] == 0:
                    print(f"[SKIP ENRICH] Empty DEG file: {deg_file}")
                    continue
            except Exception as e:
                print(f"[SKIP ENRICH] Cannot read {deg_file}: {e}")
                continue

            print(f"[ENRICH] Processing: {root}")

            cmd = [
                rscript,
                r_code,
                "-f", deg_file,
                "-s", species,
                "-o", root
            ]

            subprocess.run(cmd)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--group_col", default='group', help="组别所在的列名")

    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    adata = sc.read_h5ad(args.input_file)

    # 确保 log1p 基础设置正确，避免 scanpy 报错
    if 'log1p' not in adata.uns:
        sc.pp.log1p(adata)

    adata.uns['log1p']['base'] = None

    # 1. 按细胞类型找 Marker
    find_marker_gene(
        adata,
        group_col='cell_type',
        out_dir=args.out_dir
    )

    # 2. 每个 cell_type 内，各组 vs K；同时增加 All 全部细胞
    find_group_degs(
        adata,
        cell_type_col='cell_type',
        cond_col=args.group_col,
        out_dir=args.out_dir
    )

    # 3. 对所有生成的结果运行富集分析
    run_enrichment_recursive(
        os.path.join(args.out_dir, 'cell_type'),
        args.species
    )

    run_enrichment_recursive(
        os.path.join(args.out_dir, 'group'),
        args.species
    )