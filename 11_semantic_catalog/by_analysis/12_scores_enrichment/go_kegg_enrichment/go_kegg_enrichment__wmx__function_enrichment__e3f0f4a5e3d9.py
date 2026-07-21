import os
import scanpy as sc
import subprocess
import pandas as pd
from itertools import combinations

rscript = "/mnt/NexData_file/share/scRNA-tools/envs/pip_diffgene_enrich/bin/Rscript"
r_code = "/mnt/NexData_file/home/lst/scrna-code/diff/do_enrichment.R"

# --- 1. Cell Type 之间的 Marker 分析 ---
def find_marker_gene(adata, group_col='cell_type', out_dir=os.getcwd()):
    print(f"[INFO] Starting Marker gene analysis for {group_col}...")
    # 剔除细胞数过少的群
    counts = adata.obs[group_col].value_counts()
    keep = counts[counts > 1].index
    adata_sub = adata[adata.obs[group_col].isin(keep)].copy()
    
    sc.tl.rank_genes_groups(adata_sub, group_col, method='wilcoxon', pts=True)
    
    groups = adata_sub.uns['rank_genes_groups']['names'].dtype.names
    for c in groups:
        # 结果存放在 cell_type 目录下
        target_dir = os.path.join(out_dir, 'cell_type', str(c))
        os.makedirs(target_dir, exist_ok=True)
        
        df = sc.get.rank_genes_groups_df(adata_sub, group=c)
        # 过滤差异基因
        df_filter = df[(df['logfoldchanges'] > 0.25) & (df['pvals_adj'] < 0.05) & (df['pct_nz_group'] >= 0.1)]
        df_filter.to_csv(os.path.join(target_dir, 'differentially_expressed_genes.xls'), sep='\t', index=None)

# --- 2. 全部细胞 或 同一 Cell Type 内的多组别两两对比 ---
def find_group_degs(adata, cell_type_col='cell_type', cond_col='group', out_dir=os.getcwd()):
    print(f"[INFO] Starting Group DE analysis within each cell type...")
    
    # 建立输出的根目录
    group_base_dir = os.path.join(out_dir, 'group')
    os.makedirs(group_base_dir, exist_ok=True)
    groups = [g for g in adata.obs[cond_col].unique() if pd.notna(g)]

    cell_types = adata.obs[cell_type_col].unique()
    
    for ct in cell_types:
        # 提取当前细胞类型的子集
        adata_ct = adata[adata.obs[cell_type_col] == ct]
        groups = [g for g in adata_ct.obs[cond_col].unique() if pd.notna(g)]
        
        if len(groups) < 2:
            print(f"[SKIP] {ct}: Less than 2 groups found.")
            continue
            
        # 生成两两对比组合
        for g1, g2 in combinations(groups, 2):
            n_g1 = (adata_ct.obs[cond_col] == g1).sum()
            n_g2 = (adata_ct.obs[cond_col] == g2).sum()
            
            # 统计学样本量检查
            if n_g1 < 3 or n_g2 < 3:
                print(f"[SKIP] {ct} ({g1} vs {g2}): Sample size too small ({n_g1} vs {n_g2})")
                continue
            
            # --- 1. 执行 g1 vs g2 (以 g2 为对照) ---
            comparison_name_1 = f"{ct}_{g1}_vs_{g2}"
            target_dir_1 = os.path.join(group_base_dir, comparison_name_1)
            os.makedirs(target_dir_1, exist_ok=True)

            print(f"[RUN] Comparing {comparison_name_1}")
            sc.tl.rank_genes_groups(adata_ct, cond_col, groups=[g1], reference=g2, method='wilcoxon', pts=True)
            
            # 获取 g1 的差异结果
            df1 = sc.get.rank_genes_groups_df(adata_ct, group=g1)
            df_filter1 = df1[(df1['logfoldchanges'].abs() > 0.25) & (df1['pvals_adj'] < 0.05)]
            df_filter1.to_csv(os.path.join(target_dir_1, 'differentially_expressed_genes.xls'), sep='\t', index=None)

            # --- 2. 执行 g2 vs g1 (以 g1 为对照) ---
            comparison_name_2 = f"{ct}_{g2}_vs_{g1}"
            target_dir_2 = os.path.join(group_base_dir, comparison_name_2)
            os.makedirs(target_dir_2, exist_ok=True)

            print(f"[RUN] Comparing {comparison_name_2}")
            # 这里的 groups 设为 [g2]
            sc.tl.rank_genes_groups(adata_ct, cond_col, groups=[g2], reference=g1, method='wilcoxon', pts=True)
            
            # 提取时必须指定 group=g2，否则会取错数据
            df2 = sc.get.rank_genes_groups_df(adata_ct, group=g2)
            df_filter2 = df2[(df2['logfoldchanges'].abs() > 0.25) & (df2['pvals_adj'] < 0.05)]
            df_filter2.to_csv(os.path.join(target_dir_2, 'differentially_expressed_genes.xls'), sep='\t', index=None)

    print(f"[INFO] DE analysis finished.")

# --- 3. 递归运行富集分析 ---
def run_enrichment_recursive(base_dir, species):
    # 遍历 base_dir 下的所有子文件夹，寻找包含 xls 的目录进行富集
    for root, dirs, files in os.walk(base_dir):
        if 'differentially_expressed_genes.xls' in files:
            print(f"[ENRICH] Processing: {root}")
            deg_file = os.path.join(root, 'differentially_expressed_genes.xls')
            cmd = [
                rscript, r_code,
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

    adata = sc.read_h5ad(args.input_file)
    # 确保 log1p 基础设置正确，避免 scanpy 报错
    if 'log1p' not in adata.uns:
        sc.pp.log1p(adata)
    adata.uns['log1p']['base'] = None

    # 1. 按细胞类型找 Marker
    find_marker_gene(adata, group_col='cell_type', out_dir=args.out_dir)
    
    # 2. 按细胞类型内部的组别找 DEGs
    find_group_degs(adata, cell_type_col='cell_type', cond_col=args.group_col, out_dir=args.out_dir)
    
    # 3. 对所有生成的结果运行富集分析
    # 分别对 cell_type 和 group 文件夹进行递归处理
    run_enrichment_recursive(os.path.join(args.out_dir, 'cell_type'), args.species)
    run_enrichment_recursive(os.path.join(args.out_dir, 'group'), args.species)
    