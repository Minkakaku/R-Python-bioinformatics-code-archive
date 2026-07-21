import os
import scanpy as sc
import subprocess

rscript = "/mnt/NexData_file/share/scRNA-tools/envs/pip_diffgene_enrich/bin/Rscript"
r_code = "/mnt/NexData_file/home/lst/scrna-code/diff/do_enrichment.R"


def find_marker_gene(adata, group='cell_type',out_dir=os.getcwd(),filter = True, logfoldchanges=0.25,pvals_adj=0.05,pct_nz_group=0.1, method='wilcoxon', reference = 'rest',gene_annot=None):
    
    sample_counts = adata.obs[group].value_counts()
    problematic_groups = sample_counts[sample_counts == 1].index.to_list()
    adata = adata[~adata.obs[group].isin(problematic_groups)]
    adata.uns['log1p']['base'] = None
    sc.tl.rank_genes_groups(adata, group, method=method,reference= reference,pts=True)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    all_diff = []
    for c in groups:
        out = os.path.join(out_dir, f'{c}')
        os.makedirs(out, exist_ok=True)
        tmp = sc.get.rank_genes_groups_df(adata, group=c)
        #tmp[(tmp['pct_nz_group'] >=0.10)].to_csv(os.path.join(out, 'All_expressed_genes.xls'), sep='\t', index=None)
        if filter :
            tmp = tmp[(tmp['logfoldchanges']>logfoldchanges) & (tmp['pvals_adj']<pvals_adj) & (tmp['pct_nz_group'] >=pct_nz_group)]
        tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)


def function_enrichment(base_dir,species):
    for level1 in os.listdir(base_dir):
        out = os.path.join(base_dir, level1)
        deg_file = os.path.join(out, "differentially_expressed_genes.xls")
        if not os.path.exists(deg_file):
            print(f"[SKIP] no DEG file: {out}")
            continue
        cmd = [
            rscript,
            r_code,
            "-f", deg_file,
            "-s", species,
            "-o", out
        ]
        subprocess.run(cmd)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="annotation")
    parser.add_argument("--input_file", required=True, help="输入文件 h5ad")
    parser.add_argument("--out_dir", required=True, help="输出文件夹路径")
    parser.add_argument("--species", required=True, help="物种")
    args = parser.parse_args()
    anno_out_dir = os.path.join(args.out_dir, '06-Cell_type_annotation')
    cluster_dir = os.path.join(args.out_dir, '05-Marker_gene_analysis')
    if not os.path.exists(anno_out_dir):
        os.makedirs(anno_out_dir)
    adata = sc.read_h5ad(args.input_file)
    find_marker_gene(adata,group='cell_type', method='wilcoxon', out_dir=anno_out_dir, logfoldchanges=0.25, pvals_adj=0.05, pct_nz_group=0.1)
    function_enrichment(anno_out_dir,species)
    function_enrichment(cluster_dir,species)