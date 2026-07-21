import gseapy as gp
import os
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt
import scanpy.external as sce

def find_group_marker(adata, species = 'human' , output =os.getcwd() , logfc=1, fdr = 0.05, pct = 0.25):
        '''
        '''
        # result = adata.uns['rank_genes_groups']
        # groups = result['names'].dtype.names
        # groups = adata['cluster'].dtype.names
        groups = adata['cluster'].cat.categories.tolist()
        # all_diff = []
        for c in groups:
            out = os.path.join(output, f'cluster{c}')
            os.makedirs(out, exist_ok=True)
            # tmp = sc.get.rank_genes_groups_df(adata, group=c)
            tmp = adata[adata['cluster'] == c]
            # tmp = tmp[(tmp['logfoldchanges'].abs()>logfc) & (tmp['pvals_adj']<fdr) & (tmp['pct_nz_group'] >= pct)]
            # tmp = tmp[(tmp['logfoldchanges']>logfc) & (tmp['pvals_adj']<fdr) & (tmp['pct_nz_group'] >= pct)]
            # tmp.to_csv(os.path.join(out, 'differentially_expressed_genes.xls'), sep='\t', index=None)
            tmp['cluster'] = c
            # all_diff.append(tmp)
            try:
                if list(tmp['gene']) != []:
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
                    if list(tmp['gene']) != []:
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
        # all_diff.rename(columns={'names': 'Symbol'}, inplace=True)
        # all_diff.to_csv(os.path.join(output, 'all_cluster_marker.xls'), sep='\t', index=None)
        # all_diff['cluster'].value_counts().sort_index().to_csv(os.path.join(output, 'cluster_diff_statistics.xls'), sep='\t')

adata = pd.read_csv("/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/script/diffgene_anno.xls", sep='\t')
adata['cluster'] = adata['cluster'].astype(str).astype(object).astype("category")
find_group_marker(adata=adata,species='human',output="/mnt/NexData_file/home/sh/Projects/FWSC2024-192/11.Mission/mono_enrich/")