import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import cosg as cosg
import numpy as np
import gseapy as gp
import scanpy as sc
import os

def find_marker_gene(data, group,species=None,out_dir=os.getcwd(),  method='wilcoxon', reference = 'rest',gene_annot=None):
    """
    -data: input h5ad data
    -group: cell_type group leiden
    -species: human or mouse
    -out_dir: output dir
    -method: wilcoxon
    -reference: rest
    """
    adata = ad.read_h5ad(data)
    sample_counts = adata.obs[group].value_counts()
    problematic_groups = sample_counts[sample_counts == 1].index.to_list()
    adata = adata[~adata.obs[group].isin(problematic_groups)]
    adata.uns['log1p']['base'] = 2
    sc.tl.rank_genes_groups(adata, group, method=method,reference= reference,pts=True)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names

    all_diff = []
    for c in groups:
        out = os.path.join(out_dir, f'cluster{c}')
        os.makedirs(out, exist_ok=True)
        tmp = sc.get.rank_genes_groups_df(adata, group=c)
        tmp = tmp[(tmp['logfoldchanges']>0.25) & (tmp['pvals_adj']<0.05) & (tmp['pct_nz_group'] >= 0.25)]
        tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
        tmp['cluster'] = c
        all_diff.append(tmp)
        try:
            if list(tmp['names']) != []:
                if species == 'human':
                    enr_go = gp.enrichr(gene_list=list(tmp['names']),
                        gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human'],
                        # background='hsapiens_gene_ensembl',
                        organism='Human', 
                        cutoff=0.05,
                        no_plot=False,
                        outdir=out,
                        verbose=False
                    )
                if species == 'mouse':
                    enr_go = gp.enrichr(gene_list=list(tmp['names']),
                        gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2019_Mouse'],
                        # background='hsapiens_gene_ensembl',
                        organism='Mouse', 
                        cutoff=0.05,
                        no_plot=False,
                        outdir=out,
                        verbose=False
                        )
        except:
            try:
                if list(tmp['names']) != []:
                    if species == 'human':
                        enr_go = gp.enrichr(gene_list=list(tmp['names']),
                            gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human'],
                            # background='hsapiens_gene_ensembl',
                            organism='Human', 
                            cutoff=0.05,
                            no_plot=True,
                            outdir=out,
                            verbose=False
                        )
                    if species == 'mouse':
                        enr_go = gp.enrichr(gene_list=list(tmp['names']),
                            gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2019_Mouse'],
                            # background='hsapiens_gene_ensembl',
                            organism='Mouse', 
                            cutoff=0.05,
                            no_plot=True,
                            outdir=out,
                            verbose=False
                            )
            except:
                continue
    all_diff = pd.concat(all_diff, axis=0)
    all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
    all_diff.to_csv(os.path.join(out_dir, 'all_cluster_marker.xls'), sep='\t', index=None)
    all_diff['cluster'].value_counts().sort_index().to_csv(os.path.join(out_dir, 'cluster_diff_statistics.xls'), sep='\t')