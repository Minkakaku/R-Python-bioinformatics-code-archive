import os
import subprocess

base_dir = "/mnt/NexData_file/scRNA/scVision_output/SC2509225099/05-Marker_gene_analysis"
species = "mouse"

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

function_enrichment(base_dir,species)