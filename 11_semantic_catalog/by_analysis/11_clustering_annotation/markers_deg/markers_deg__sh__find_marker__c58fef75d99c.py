import gseapy as gp
import os
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt
import scanpy.external as sce

result = pd.read_csv("/mnt/NexData_file/home/sh/Projects/FWAE20254246/02.Mission/script/all_markers_-_RNA_snn_res.0.6.csv", sep=',', header='infer')
groups = result['cluster'].drop_duplicates().tolist()
groups = [9,12,13]
logfc=1
fdr = 0.05
pct = 0.1
species = 'human'
output = "/mnt/NexData_file/home/sh/Projects/FWAE20254246/02.Mission/enrich/"
all_diff = []
for c in groups:
    out = os.path.join(output, f'cluster{c}')
    os.makedirs(out, exist_ok=True)
    tmp = result[(result['cluster']==c)]
    # tmp = tmp[(tmp['logfoldchanges'].abs()>logfc) & (tmp['pvals_adj']<fdr) & (tmp['pct_nz_group'] >= pct)]
    tmp = tmp[(tmp['avg_log2FC']>logfc) & (tmp['p_val_adj']<fdr) & (tmp['pct.1'] >= pct)]
    tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
    # tmp['cluster'] = c
    # all_diff.append(tmp)
    try:
        if list(tmp['cluster']) != []:
            if species == 'human':
                enr_go = gp.enrichr(gene_list=list(tmp['gene']),
                    gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human'],
                    # background='hsapiens_gene_ensembl',
                    organism='Human', 
                    cutoff=0.05,
                    no_plot=False,
                    outdir=out,
                    verbose=False
                )
            if species == 'mouse':
                enr_go = gp.enrichr(gene_list=list(tmp['gene']),
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
            if list(tmp['cluster']) != []:
                if species == 'human':
                    enr_go = gp.enrichr(gene_list=list(tmp['gene']),
                        gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2021_Human'],
                        # background='hsapiens_gene_ensembl',
                        organism='Human', 
                        cutoff=0.05,
                        no_plot=True,
                        outdir=out,
                        verbose=False
                    )
                if species == 'mouse':
                    enr_go = gp.enrichr(gene_list=list(tmp['gene']),
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

# all_diff = pd.concat(all_diff, axis=0)
# all_diff.rename(columns={'cluster': 'Symbol'}, inplace=True)
# all_diff.to_csv(os.path.join(output, 'all_cluster_marker.xls'), sep='\t', index=None)
# all_diff['cluster'].value_counts().sort_index().to_csv(os.path.join(output, 'cluster_diff_statistics.xls'), sep='\t')